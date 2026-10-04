#include "rtt/analysis/chromatic.hpp"

#include <cmath>
#include <string>

#include "common.hpp"
#include "rtt/paraxial/paraxial.hpp"

namespace rtt::analysis {

namespace {

using detail::arrived;
using detail::check_path;
using detail::image_surface;
using detail::local_point;
using detail::trace_aimed;

}  // namespace

LongitudinalColour longitudinal_colour(const compile::CompiledSystem& system,
                                       compile::PathId path,
                                       const ChromaticOptions& options) {
  check_path(system, path);
  const auto count = static_cast<std::uint16_t>(system.wavelengths_um().size());
  LongitudinalColour lc;
  lc.pair =
      options.pair.value_or(paraxial::ChromaticPair{0, static_cast<std::uint16_t>(count - 1)});
  if (lc.pair.first >= count || lc.pair.second >= count) {
    throw std::invalid_argument("analysis: wavelength pair index does not exist");
  }
  if (!(options.zone > 0.0 && options.zone <= 1.0)) {
    throw std::invalid_argument("analysis: zone must lie in (0, 1]");
  }
  const std::uint32_t image = image_surface(system, path);
  const model::Field axis{0.0, 0.0, 1.0};
  int direction = 1;
  for (std::uint16_t wl = 0; wl < count; ++wl) {
    FocusPosition f;
    f.wavelength = wl;
    // Paraxial focus: rear focal point (object at infinity) or paraxial image of the axial
    // object point.
    const paraxial::FirstOrder fo = paraxial::first_order(system, path, wl);
    const std::optional<double> z = system.object().at_infinity ? fo.rear_focal_z : fo.image_z;
    if (!z) throw AnalysisError("analysis: no paraxial focus (afocal system)");
    f.paraxial_z = *z;
    if (wl == system.reference_wavelength()) direction = fo.image_direction;
    // Real focus: axis crossing of the meridional ray at (0, zone) of the on-axis field.
    const trace::AimedRay aimed =
        trace::aim_ray(system, path, axis, wl, 0.0, options.zone, options.aiming);
    const trace::RayBatch rays = trace_aimed(system, path, aimed, wl);
    if (!arrived(rays, 0, image)) {
      throw AnalysisError("analysis: the zone ray at wavelength " + std::to_string(wl) +
                          " does not reach the image surface");
    }
    const double dy = rays.dir_y()[0];
    if (dy == 0.0) throw AnalysisError("analysis: the zone ray does not cross the axis");
    f.real_z = rays.pos_z()[0] - rays.pos_y()[0] * rays.dir_z()[0] / dy;
    lc.foci.push_back(f);
  }
  const FocusPosition& a = lc.foci[lc.pair.first];
  const FocusPosition& b = lc.foci[lc.pair.second];
  lc.paraxial = (a.paraxial_z - b.paraxial_z) * direction;
  lc.real = (a.real_z - b.real_z) * direction;
  return lc;
}

LateralColour lateral_colour(const compile::CompiledSystem& system,
                             compile::PathId path,
                             std::uint16_t field,
                             trace::Aiming aiming) {
  check_path(system, path);
  LateralColour lat;
  lat.field = field;
  const std::uint32_t image = image_surface(system, path);
  const auto count = static_cast<std::uint16_t>(system.wavelengths_um().size());
  for (std::uint16_t wl = 0; wl < count; ++wl) {
    const trace::AimedRay aimed = trace::aim_ray(system, path, field, wl, 0.0, 0.0, aiming);
    const trace::RayBatch rays = trace_aimed(system, path, aimed, wl);
    if (!arrived(rays, 0, image)) {
      throw AnalysisError("analysis: the chief ray at wavelength " + std::to_string(wl) +
                          " does not reach the image surface");
    }
    lat.chief.push_back(local_point(system, rays, 0, image));
  }
  const Point2 ref = lat.chief[system.reference_wavelength()];
  for (const Point2& c : lat.chief) lat.offset.push_back({c.x - ref.x, c.y - ref.y});
  return lat;
}

}  // namespace rtt::analysis
