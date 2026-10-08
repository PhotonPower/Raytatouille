#include "rtt/analysis/ghosts.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "common.hpp"
#include "rtt/paraxial/paraxial.hpp"
#include "rtt/paraxial/seidel.hpp"
#include "rtt/trace/sequential.hpp"

namespace rtt::analysis {
namespace {

using compile::CompiledSystem;
using compile::PathId;
using trace::RayBatch;

/// A ghost whose paraxial focus lies farther than this factor times |y| from its last surface
/// (|u| <= kCollimatedSlopeRatio |y|) leaves collimated and has no focus offset; this keeps
/// rounding residues of u from giving foci at 1e17 mm.
constexpr double kCollimatedSlopeRatio = 1e-12;

/// Transmitted power (as path_transmission, #122: mean final weight over the launched rays),
/// weighted RMS radius about the centroid and arrived rays of `start` traced on `path`.
struct PathSpot {
  double power = 0.0;
  double rms = 0.0;
  std::size_t arrived = 0;
  RayLosses losses;
  std::vector<model::Diagnostic> warnings;
};

PathSpot path_spot(const CompiledSystem& system,
                   PathId path,
                   const RayBatch& start,
                   double lost_warning_fraction,
                   const trace::RunControl* control) {
  detail::LossCounter losses(system, path, lost_warning_fraction);
  RayBatch rays = start;
  const trace::SequentialTracer tracer;
  [[maybe_unused]] const auto stats = control != nullptr
                                          ? tracer.trace(system, path, rays, *control)
                                          : tracer.trace(system, path, rays);
  const std::uint32_t image = detail::image_surface(system, path);
  // Serial sums in batch order: independent of the number of threads.
  double sum_w = 0.0;
  double sum_x = 0.0;
  double sum_y = 0.0;
  PathSpot out;
  for (std::size_t i = 0; i < rays.size(); ++i) {
    if (!detail::arrived(rays, i, image)) continue;
    const double w = rays.weight()[i];
    const Point2 p = detail::local_point(system, rays, i, image);
    sum_w += w;
    sum_x += w * p.x;
    sum_y += w * p.y;
    ++out.arrived;
  }
  out.power = sum_w / static_cast<double>(rays.size());
  if (sum_w > 0.0) {
    const double cx = sum_x / sum_w;
    const double cy = sum_y / sum_w;
    double sum_r2 = 0.0;
    for (std::size_t i = 0; i < rays.size(); ++i) {
      if (!detail::arrived(rays, i, image)) continue;
      const Point2 p = detail::local_point(system, rays, i, image);
      sum_r2 += rays.weight()[i] * ((p.x - cx) * (p.x - cx) + (p.y - cy) * (p.y - cy));
    }
    out.rms = std::sqrt(sum_r2 / sum_w);
  }
  losses.add(rays, image);
  out.losses = losses.result();
  out.warnings = losses.warnings();
  return out;
}

/// Paraxial marginal ray of the base path (axial object point to the rim of the entrance
/// pupil), as in rtt-paraxial (marginal_start); none without a finite entrance pupil (e.g.
/// object-side telecentric), then the ranking has no paraxial diagnostics.
std::optional<paraxial::RayStart> marginal_start(const CompiledSystem& system,
                                                 PathId base,
                                                 std::uint16_t wl) {
  const paraxial::FirstOrder fo = paraxial::first_order(system, base, wl);
  if (!fo.entrance_pupil || !fo.entrance_pupil->z || !fo.entrance_pupil->diameter) {
    return std::nullopt;
  }
  const double z_ep = *fo.entrance_pupil->z;
  const double r_ep = 0.5 * *fo.entrance_pupil->diameter;
  if (system.object().at_infinity) return paraxial::RayStart{z_ep, r_ep, 0.0};
  const double z_obj = -system.object().distance.value;
  return paraxial::RayStart{z_obj, 0.0, r_ep / (z_ep - z_obj)};
}

GhostRanking rank(const compile::GhostSystem& ghosts,
                  std::uint16_t field,
                  std::uint16_t wavelength,
                  const GhostRankingOptions& options,
                  const trace::RunControl* control) {
  const CompiledSystem& system = ghosts.system;
  if (!(options.resolution_radius > 0.0) || !std::isfinite(options.resolution_radius)) {
    throw std::invalid_argument("ghost ranking: resolution_radius must be finite and > 0");
  }
  if (!(options.lost_warning_fraction >= 0.0 && options.lost_warning_fraction <= 1.0)) {
    throw std::invalid_argument("analysis: lost_warning_fraction must lie in [0, 1]");
  }
  if (ghosts.ghosts.empty()) {
    throw std::invalid_argument("ghost ranking: the system has no ghosts");
  }
  const PathId base = ghosts.ghosts.front().base;
  detail::check_path(system, base);
  detail::check_wavelength(system, wavelength);
  if (field >= system.fields().points.size()) {
    throw std::invalid_argument("analysis: field index " + std::to_string(field) +
                                " does not exist");
  }
  const std::uint16_t fields[] = {field};
  const RayBatch start =
      control != nullptr
          ? trace::make_rays(system, base, fields, wavelength, options.sampling, options.aiming,
                             *control)
          : trace::make_rays(system, base, fields, wavelength, options.sampling, options.aiming);
  if (start.size() == 0) throw std::invalid_argument("analysis: no start rays");
  const std::optional<paraxial::RayStart> marginal = marginal_start(system, base, wavelength);

  GhostRanking out;
  out.base = base;
  out.field = field;
  out.wavelength = wavelength;
  out.resolution_radius = options.resolution_radius;
  const PathSpot useful = path_spot(system, base, start, options.lost_warning_fraction, control);
  if (!(useful.power > 0.0)) {
    throw AnalysisError("ghost ranking: no ray of the base path arrives at its image surface");
  }
  out.base_power = useful.power;
  out.base_rms_radius = useful.rms;
  out.warnings = useful.warnings;
  const double r0 = options.resolution_radius;
  const double base_area = useful.rms * useful.rms + r0 * r0;

  out.entries.reserve(ghosts.ghosts.size());
  for (const compile::GhostPath& ghost : ghosts.ghosts) {
    // Ghosts warn about nothing: losses are normal for them (ADR 0027, addendum #124).
    const PathSpot spot = path_spot(system, ghost.path, start, 1.0, control);
    GhostEntry e;
    e.path = ghost.path;
    e.surface_j = ghost.surface_j;
    e.surface_i = ghost.surface_i;
    e.power = spot.power;
    e.relative_power = spot.power / useful.power;
    e.rms_radius = spot.rms;
    e.rays_arrived = spot.arrived;
    // r0 > 0: the area is positive (ADR 0027, addendum #124).
    const double area = spot.rms * spot.rms + r0 * r0;
    e.relative_irradiance = e.relative_power * base_area / area;
    e.losses = spot.losses;
    if (marginal) {
      // Paraxial diagnostics: the base marginal ray on the ghost path. The focus is where the
      // ray after the last event before the image meets the axis, by the transfer
      // y' = y + u t with y' = 0 (Greivenkamp, OPTI-201/202, Sec. 9, p. 9-2; docs/quellen.md).
      const auto ray = paraxial::trace_ray(system, ghost.path, wavelength, marginal->z, marginal->y,
                                           marginal->u);
      if (ray.size() >= 2) {
        const paraxial::RayAtEvent& before = ray[ray.size() - 2];
        if (std::abs(before.u) > kCollimatedSlopeRatio * std::abs(before.y)) {
          e.focus_offset = before.z - before.y / before.u - ray.back().z;
        }
      }
      e.paraxial_blur_radius = std::abs(ray.back().y);
    }
    out.entries.push_back(e);
  }
  // Rank value descending; stable: ties keep the ghost order.
  std::stable_sort(out.entries.begin(), out.entries.end(),
                   [](const GhostEntry& a, const GhostEntry& b) {
                     return a.relative_irradiance > b.relative_irradiance;
                   });
  return out;
}

}  // namespace

GhostRanking ghost_ranking(const compile::GhostSystem& ghosts,
                           std::uint16_t field,
                           std::uint16_t wavelength,
                           const GhostRankingOptions& options) {
  return rank(ghosts, field, wavelength, options, nullptr);
}

GhostRanking ghost_ranking(const compile::GhostSystem& ghosts,
                           std::uint16_t field,
                           std::uint16_t wavelength,
                           const GhostRankingOptions& options,
                           const trace::RunControl& control) {
  return rank(ghosts, field, wavelength, options, &control);
}

}  // namespace rtt::analysis
