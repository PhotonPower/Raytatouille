#include "rtt/analysis/spot.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <string>

#include "rtt/trace/sequential.hpp"

namespace rtt::analysis {

namespace {

using compile::CompiledSystem;
using compile::PathId;

void check_path(const CompiledSystem& system, PathId path) {
  if (path.index >= system.paths().size() || system.path(path).events.empty()) {
    throw std::invalid_argument("analysis: path index " + std::to_string(path.index) +
                                " does not exist or has no events");
  }
}

void check_wavelength(const CompiledSystem& system, std::uint16_t wavelength) {
  if (wavelength >= system.wavelengths_um().size()) {
    throw std::invalid_argument("analysis: wavelength index " + std::to_string(wavelength) +
                                " does not exist");
  }
}

/// Image surface: the surface of the last event of the path (decided for #28).
std::uint32_t image_surface(const CompiledSystem& system, PathId path) {
  return system.path(path).events.back().surface;
}

/// Traces rays of one field and wavelength through the whole path.
trace::RayBatch trace_rays(const CompiledSystem& system,
                           PathId path,
                           std::uint16_t field,
                           std::uint16_t wavelength,
                           const trace::PupilSampling& sampling,
                           trace::Aiming aiming) {
  const std::uint16_t fields[] = {field};
  trace::RayBatch rays = trace::make_rays(system, path, fields, wavelength, sampling, aiming);
  [[maybe_unused]] const auto stats = trace::SequentialTracer().trace(system, path, rays);
  return rays;
}

/// True if ray i ended Alive on the image surface.
bool arrived(const trace::RayBatch& rays, std::size_t i, std::uint32_t image) {
  return rays.status()[i] == trace::RayStatus::Alive && rays.last_surface()[i] == image;
}

/// Local x, y of ray i on the image surface, mm.
Point2 local_point(const CompiledSystem& system,
                   const trace::RayBatch& rays,
                   std::size_t i,
                   std::uint32_t image) {
  const math::Vec3 p = system.surfaces()[image].to_local.apply_point(
      math::Vec3(rays.pos_x()[i], rays.pos_y()[i], rays.pos_z()[i]));
  return {p.x(), p.y()};
}

/// Chief ray (pupil centre, reference wavelength) on the image surface.
Point2 chief_point(const CompiledSystem& system,
                   PathId path,
                   std::uint16_t field,
                   trace::Aiming aiming) {
  const std::uint32_t image = image_surface(system, path);
  const trace::RayBatch rays = trace_rays(system, path, field, system.reference_wavelength(),
                                          trace::SinglePupilPoint{0.0, 0.0}, aiming);
  if (!arrived(rays, 0, image)) {
    throw AnalysisError("analysis: the chief ray of field " + std::to_string(field) +
                        " does not reach the image surface");
  }
  return local_point(system, rays, 0, image);
}

}  // namespace

SpotStatistics spot_statistics(std::span<const SpotPoint> points, Point2 chief) {
  // Definitions decided for #28; sums in the order of `points`.
  double sum_w = 0.0;
  double cx = 0.0;
  double cy = 0.0;
  for (const SpotPoint& p : points) {
    sum_w += p.weight;
    cx += p.weight * p.x;
    cy += p.weight * p.y;
  }
  if (points.empty() || !(sum_w > 0.0)) {
    throw AnalysisError("analysis: no spot points with positive weight");
  }
  SpotStatistics s;
  s.centroid = {cx / sum_w, cy / sum_w};
  double sum_c = 0.0;
  double sum_chief = 0.0;
  for (const SpotPoint& p : points) {
    const double dcx = p.x - s.centroid.x;
    const double dcy = p.y - s.centroid.y;
    const double dhx = p.x - chief.x;
    const double dhy = p.y - chief.y;
    sum_c += p.weight * (dcx * dcx + dcy * dcy);
    sum_chief += p.weight * (dhx * dhx + dhy * dhy);
    if (p.weight > 0.0) {  // points without weight do not belong to the spot
      s.geo_centroid = std::max(s.geo_centroid, std::hypot(dcx, dcy));
      s.geo_chief = std::max(s.geo_chief, std::hypot(dhx, dhy));
    }
  }
  s.rms_centroid = std::sqrt(sum_c / sum_w);
  s.rms_chief = std::sqrt(sum_chief / sum_w);
  return s;
}

SpotDiagram spot(const compile::CompiledSystem& system,
                 compile::PathId path,
                 std::uint16_t field,
                 std::optional<std::uint16_t> wavelength,
                 const SpotOptions& options) {
  check_path(system, path);
  if (field >= system.fields().points.size()) {
    throw std::invalid_argument("analysis: field index " + std::to_string(field) +
                                " does not exist");
  }
  std::vector<std::uint16_t> wavelengths;
  std::vector<double> weights;
  if (wavelength) {
    check_wavelength(system, *wavelength);
    wavelengths = {*wavelength};
    weights = {1.0};
  } else {
    const auto& w = system.wavelength_weights();
    const double total = std::accumulate(w.begin(), w.end(), 0.0);
    if (!(total > 0.0)) {
      throw std::invalid_argument(
          "analysis: the wavelength weights do not sum to a positive "
          "value");
    }
    for (std::size_t i = 0; i < w.size(); ++i) {
      wavelengths.push_back(static_cast<std::uint16_t>(i));
      weights.push_back(w[i] / total);  // normalised to sum 1 (decided for #28)
    }
  }

  SpotDiagram d;
  d.field = field;
  d.wavelength = wavelength;
  d.image_surface = image_surface(system, path);
  d.chief = chief_point(system, path, field, options.aiming);
  for (std::size_t k = 0; k < wavelengths.size(); ++k) {
    const trace::RayBatch rays =
        trace_rays(system, path, field, wavelengths[k], options.sampling, options.aiming);
    d.rays_launched += rays.size();
    for (std::size_t i = 0; i < rays.size(); ++i) {
      if (!arrived(rays, i, d.image_surface)) continue;
      const Point2 p = local_point(system, rays, i, d.image_surface);
      d.points.push_back({p.x, p.y, wavelengths[k], weights[k] * rays.weight()[i]});
    }
  }
  d.rays_arrived = d.points.size();
  d.vignetted_fraction = d.rays_launched == 0
                             ? 0.0
                             : static_cast<double>(d.rays_launched - d.rays_arrived) /
                                   static_cast<double>(d.rays_launched);
  d.stats = spot_statistics(d.points, d.chief);
  return d;
}

RayFan ray_fan(const compile::CompiledSystem& system,
               compile::PathId path,
               std::uint16_t field,
               std::uint16_t wavelength,
               const FanOptions& options) {
  check_path(system, path);
  check_wavelength(system, wavelength);
  if (options.points < 1) throw std::invalid_argument("analysis: fan points must be >= 1");
  RayFan fan;
  fan.field = field;
  fan.wavelength = wavelength;
  fan.image_surface = image_surface(system, path);
  fan.chief = chief_point(system, path, field, options.aiming);

  const auto make_fan = [&](const trace::PupilSampling& sampling, bool tangential) {
    const trace::RayBatch rays =
        trace_rays(system, path, field, wavelength, sampling, options.aiming);
    std::vector<FanPoint> points;
    points.reserve(rays.size());
    for (std::size_t i = 0; i < rays.size(); ++i) {
      FanPoint p;
      p.p = tangential ? rays.pupil_y()[i] : rays.pupil_x()[i];
      if (arrived(rays, i, fan.image_surface)) {
        const Point2 q = local_point(system, rays, i, fan.image_surface);
        p.ex = q.x - fan.chief.x;
        p.ey = q.y - fan.chief.y;
      } else {
        // Alive but stopped before the image surface cannot happen with a complete trace;
        // count it as vignetted all the same.
        p.status = rays.status()[i] == trace::RayStatus::Alive ? trace::RayStatus::Vignetted
                                                               : rays.status()[i];
      }
      points.push_back(p);
    }
    return points;
  };
  fan.tangential = make_fan(trace::FanYPupil{options.points}, true);
  fan.sagittal = make_fan(trace::FanXPupil{options.points}, false);
  return fan;
}

}  // namespace rtt::analysis
