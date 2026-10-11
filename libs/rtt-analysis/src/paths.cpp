#include "rtt/analysis/paths.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

#include "common.hpp"
#include "rtt/math/types.hpp"

namespace rtt::analysis {
namespace {

using compile::CompiledSystem;
using compile::PathId;
using trace::RayBatch;
using trace::RayStatus;

/// Status of ray i on a path: Alive if it arrived, its own status if it was lost on the way,
/// Vignetted if it stopped Alive on another surface (as in the other analyses).
RayStatus final_status(const RayBatch& rays, std::size_t i, std::uint32_t image) {
  if (detail::arrived(rays, i, image)) return RayStatus::Alive;
  return rays.status()[i] == RayStatus::Alive ? RayStatus::Vignetted : rays.status()[i];
}

/// A copy of `start` traced on `path` (rays that are not Alive at the start are left as they
/// are by the tracer), with or without a control.
RayBatch traced(const CompiledSystem& system,
                PathId path,
                const RayBatch& start,
                const trace::RunControl* control) {
  RayBatch rays = start;
  const trace::SequentialTracer tracer;
  [[maybe_unused]] const auto stats = control != nullptr
                                          ? tracer.trace(system, path, rays, *control)
                                          : tracer.trace(system, path, rays);
  return rays;
}

void check_start(const RayBatch& start) {
  if (start.size() == 0) throw std::invalid_argument("analysis: no start rays");
}

void check_options(const PathOptions& options) {
  if (!(options.lost_warning_fraction >= 0.0 && options.lost_warning_fraction <= 1.0)) {
    throw std::invalid_argument("analysis: lost_warning_fraction must lie in [0, 1]");
  }
}

/// Checks both paths of an OPL difference and returns their common image surface.
std::uint32_t common_image(const CompiledSystem& system, PathId path_a, PathId path_b) {
  detail::check_path(system, path_a);
  detail::check_path(system, path_b);
  const std::uint32_t image = detail::image_surface(system, path_a);
  if (detail::image_surface(system, path_b) != image) {
    throw std::invalid_argument(
        "analysis: an OPL difference needs two paths that end on the same surface");
  }
  return image;
}

PathTransmission transmission(const CompiledSystem& system,
                              PathId path,
                              const RayBatch& start,
                              const PathOptions& options,
                              const trace::RunControl* control) {
  detail::check_path(system, path);
  check_start(start);
  detail::LossCounter losses(system, path, options.lost_warning_fraction);
  const RayBatch rays = traced(system, path, start, control);
  PathTransmission out;
  out.path = path;
  out.image_surface = detail::image_surface(system, path);
  out.rays_launched = rays.size();
  out.rays.reserve(rays.size());
  // Serial sums in batch order: independent of the number of threads.
  double sum = 0.0;
  for (std::size_t i = 0; i < rays.size(); ++i) {
    PathRay r;
    r.px = start.pupil_x()[i];
    r.py = start.pupil_y()[i];
    r.status = final_status(rays, i, out.image_surface);
    if (r.status == RayStatus::Alive) {
      r.weight = rays.weight()[i];
      out.min = out.rays_arrived == 0 ? r.weight : std::min(out.min, r.weight);
      out.max = out.rays_arrived == 0 ? r.weight : std::max(out.max, r.weight);
      ++out.rays_arrived;
      sum += r.weight;
    }
    out.rays.push_back(r);
  }
  out.mean = sum / static_cast<double>(rays.size());
  losses.add(rays, out.image_surface);
  out.losses = losses.result();
  out.warnings = losses.warnings();
  return out;
}

PathOplDifference difference(const CompiledSystem& system,
                             PathId path_a,
                             PathId path_b,
                             const RayBatch& start,
                             const PathOptions& options,
                             const trace::RunControl* control) {
  const std::uint32_t image = common_image(system, path_a, path_b);
  check_start(start);
  detail::LossCounter losses_a(system, path_a, options.lost_warning_fraction);
  detail::LossCounter losses_b(system, path_b, options.lost_warning_fraction);
  const RayBatch a = traced(system, path_a, start, control);
  const RayBatch b = traced(system, path_b, start, control);
  PathOplDifference out;
  out.path_a = path_a;
  out.path_b = path_b;
  out.image_surface = image;
  out.points.reserve(start.size());
  for (std::size_t i = 0; i < start.size(); ++i) {
    OplDifferencePoint p;
    p.px = start.pupil_x()[i];
    p.py = start.pupil_y()[i];
    const RayStatus status_a = final_status(a, i, image);
    const RayStatus status_b = final_status(b, i, image);
    p.status = status_a != RayStatus::Alive ? status_a : status_b;
    if (p.status == RayStatus::Alive) {
      p.delta = b.opl()[i] - a.opl()[i];
      if (!out.chief && p.px == 0.0 && p.py == 0.0) out.chief = p.delta;
    }
    out.points.push_back(p);
  }
  losses_a.add(a, image);
  losses_b.add(b, image);
  out.losses_a = losses_a.result();
  out.losses_b = losses_b.result();
  out.warnings = losses_a.warnings();
  for (auto& w : losses_b.warnings()) out.warnings.push_back(std::move(w));
  return out;
}

/// Rays of trace::make_rays for the convenience forms; the cheap checks of the options, the
/// path and the wavelength come before the aiming.
RayBatch start_rays(const CompiledSystem& system,
                    PathId path,
                    std::uint16_t field,
                    std::uint16_t wavelength,
                    const PathOptions& options,
                    const trace::RunControl* control) {
  check_options(options);
  detail::reject_gauss(options.sampling, "path analysis");
  detail::check_path(system, path);
  detail::check_wavelength(system, wavelength);
  const std::uint16_t fields[] = {field};
  return control != nullptr
             ? trace::make_rays(system, path, fields, wavelength, options.sampling, options.aiming,
                                *control)
             : trace::make_rays(system, path, fields, wavelength, options.sampling, options.aiming);
}

/// #210: the transfer of wave b needs the index along the ray direction, so the image medium
/// of each path must be isotropic (a crystal mode has its own wave normal and index; H6, #35).
void check_isotropic_image(const CompiledSystem& system, PathId path) {
  const auto& last = system.path(path).events.back();
  if (system.media()[last.medium_after].optic_axis.has_value()) {
    throw std::invalid_argument(
        "analysis: an interferogram needs an isotropic image medium, not a crystal (#35)");
  }
}

/// Checks of the interferogram that come before any ray is traced or aimed.
std::uint32_t interferogram_image(const CompiledSystem& system, PathId path_a, PathId path_b) {
  const std::uint32_t image = common_image(system, path_a, path_b);
  check_isotropic_image(system, path_a);
  check_isotropic_image(system, path_b);
  return image;
}

math::Vec3 position(const RayBatch& r, std::size_t i) {
  return {r.pos_x()[i], r.pos_y()[i], r.pos_z()[i]};
}

Interferogram interfere(const CompiledSystem& system,
                        PathId path_a,
                        PathId path_b,
                        const RayBatch& start,
                        const PathOptions& options,
                        const trace::RunControl* control) {
  const std::uint32_t image = interferogram_image(system, path_a, path_b);
  check_start(start);
  check_options(options);
  detail::LossCounter losses_a(system, path_a, options.lost_warning_fraction);
  detail::LossCounter losses_b(system, path_b, options.lost_warning_fraction);
  const RayBatch a = traced(system, path_a, start, control);
  const RayBatch b = traced(system, path_b, start, control);
  const auto& medium_b = system.media()[system.path(path_b).events.back().medium_after];
  Interferogram out;
  out.path_a = path_a;
  out.path_b = path_b;
  out.image_surface = image;
  out.points.reserve(start.size());
  // Serial in batch order: independent of the number of threads.
  for (std::size_t i = 0; i < start.size(); ++i) {
    InterferogramPoint p;
    p.px = start.pupil_x()[i];
    p.py = start.pupil_y()[i];
    p.status_a = final_status(a, i, image);
    p.status_b = final_status(b, i, image);
    const bool on_a = p.status_a == RayStatus::Alive;
    const bool on_b = p.status_b == RayStatus::Alive;
    if (on_a && on_b) {
      // Wave b carried from its hit x_b to x_a along d_b: OPL_b(x_a) = OPL_b(x_b) +
      // n'_b d_b . (x_a - x_b), exact for a plane wave (BSc Optics, Sec. 1.5.1, Eqs.
      // (1.37)-(1.39)), first order for a curved one (Sec. 2.1).
      const math::Vec3 xa = position(a, i);
      const math::Vec3 xb = position(b, i);
      const math::Vec3 db(b.dir_x()[i], b.dir_y()[i], b.dir_z()[i]);
      const std::uint16_t wl = a.wl()[i];
      const double n_b = std::abs(medium_b.index[wl].real());
      p.delta = b.opl()[i] + n_b * db.dot(xa - xb) - a.opl()[i];
      const double lambda_mm = system.wavelengths_um()[wl] * 1e-3;
      p.waves = p.delta / lambda_mm;
      p.phase = 2.0 * std::numbers::pi * p.waves;
      // Two-beam interference, BSc Optics, Sec. 5.2, Eq. (5.12); equal amplitudes Eqs.
      // (5.10), (5.11).
      const double c = std::cos(p.phase);
      p.fringe = 1.0 + c;
      const double w_a = a.weight()[i];
      const double w_b = b.weight()[i];
      p.intensity = w_a + w_b + 2.0 * std::sqrt(w_a * w_b) * c;
      p.shear = (xa - xb).norm();
      out.max_shear = std::max(out.max_shear, p.shear);
      p.arrival = InterferogramArrival::Both;
      const Point2 local = detail::local_point(system, a, i, image);
      p.x = local.x;
      p.y = local.y;
    } else if (on_a) {
      p.arrival = InterferogramArrival::OnlyA;
      p.intensity = a.weight()[i];
      const Point2 local = detail::local_point(system, a, i, image);
      p.x = local.x;
      p.y = local.y;
    } else if (on_b) {
      p.arrival = InterferogramArrival::OnlyB;
      p.intensity = b.weight()[i];
      const Point2 local = detail::local_point(system, b, i, image);
      p.x = local.x;
      p.y = local.y;
    } else {
      p.arrival = InterferogramArrival::None;
    }
    out.points.push_back(p);
  }
  losses_a.add(a, image);
  losses_b.add(b, image);
  out.losses_a = losses_a.result();
  out.losses_b = losses_b.result();
  out.warnings = losses_a.warnings();
  for (auto& w : losses_b.warnings()) out.warnings.push_back(std::move(w));
  return out;
}

}  // namespace

PathTransmission path_transmission(const compile::CompiledSystem& system,
                                   compile::PathId path,
                                   const trace::RayBatch& start,
                                   const PathOptions& options) {
  return transmission(system, path, start, options, nullptr);
}

PathTransmission path_transmission(const compile::CompiledSystem& system,
                                   compile::PathId path,
                                   std::uint16_t field,
                                   std::uint16_t wavelength,
                                   const PathOptions& options) {
  return transmission(system, path, start_rays(system, path, field, wavelength, options, nullptr),
                      options, nullptr);
}

PathOplDifference opl_difference(const compile::CompiledSystem& system,
                                 compile::PathId path_a,
                                 compile::PathId path_b,
                                 const trace::RayBatch& start,
                                 const PathOptions& options) {
  return difference(system, path_a, path_b, start, options, nullptr);
}

PathOplDifference opl_difference(const compile::CompiledSystem& system,
                                 compile::PathId path_a,
                                 compile::PathId path_b,
                                 std::uint16_t field,
                                 std::uint16_t wavelength,
                                 const PathOptions& options) {
  [[maybe_unused]] const std::uint32_t image = common_image(system, path_a, path_b);
  return difference(system, path_a, path_b,
                    start_rays(system, path_a, field, wavelength, options, nullptr), options,
                    nullptr);
}

PathTransmission path_transmission(const compile::CompiledSystem& system,
                                   compile::PathId path,
                                   const trace::RayBatch& start,
                                   const PathOptions& options,
                                   const trace::RunControl& control) {
  return transmission(system, path, start, options, &control);
}

PathTransmission path_transmission(const compile::CompiledSystem& system,
                                   compile::PathId path,
                                   std::uint16_t field,
                                   std::uint16_t wavelength,
                                   const PathOptions& options,
                                   const trace::RunControl& control) {
  return transmission(system, path, start_rays(system, path, field, wavelength, options, &control),
                      options, &control);
}

PathOplDifference opl_difference(const compile::CompiledSystem& system,
                                 compile::PathId path_a,
                                 compile::PathId path_b,
                                 const trace::RayBatch& start,
                                 const PathOptions& options,
                                 const trace::RunControl& control) {
  return difference(system, path_a, path_b, start, options, &control);
}

PathOplDifference opl_difference(const compile::CompiledSystem& system,
                                 compile::PathId path_a,
                                 compile::PathId path_b,
                                 std::uint16_t field,
                                 std::uint16_t wavelength,
                                 const PathOptions& options,
                                 const trace::RunControl& control) {
  [[maybe_unused]] const std::uint32_t image = common_image(system, path_a, path_b);
  return difference(system, path_a, path_b,
                    start_rays(system, path_a, field, wavelength, options, &control), options,
                    &control);
}

Interferogram interferogram(const compile::CompiledSystem& system,
                            compile::PathId path_a,
                            compile::PathId path_b,
                            const trace::RayBatch& start,
                            const PathOptions& options) {
  return interfere(system, path_a, path_b, start, options, nullptr);
}

Interferogram interferogram(const compile::CompiledSystem& system,
                            compile::PathId path_a,
                            compile::PathId path_b,
                            std::uint16_t field,
                            std::uint16_t wavelength,
                            const PathOptions& options) {
  [[maybe_unused]] const std::uint32_t image = interferogram_image(system, path_a, path_b);
  return interfere(system, path_a, path_b,
                   start_rays(system, path_a, field, wavelength, options, nullptr), options,
                   nullptr);
}

Interferogram interferogram(const compile::CompiledSystem& system,
                            compile::PathId path_a,
                            compile::PathId path_b,
                            const trace::RayBatch& start,
                            const PathOptions& options,
                            const trace::RunControl& control) {
  return interfere(system, path_a, path_b, start, options, &control);
}

Interferogram interferogram(const compile::CompiledSystem& system,
                            compile::PathId path_a,
                            compile::PathId path_b,
                            std::uint16_t field,
                            std::uint16_t wavelength,
                            const PathOptions& options,
                            const trace::RunControl& control) {
  [[maybe_unused]] const std::uint32_t image = interferogram_image(system, path_a, path_b);
  return interfere(system, path_a, path_b,
                   start_rays(system, path_a, field, wavelength, options, &control), options,
                   &control);
}

}  // namespace rtt::analysis
