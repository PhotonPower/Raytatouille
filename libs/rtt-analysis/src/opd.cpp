#include "rtt/analysis/opd.hpp"

#include <algorithm>
#include <cmath>
#include <optional>
#include <stdexcept>
#include <string>

#include "common.hpp"
#include "rtt/compile/errors.hpp"
#include "rtt/paraxial/paraxial.hpp"

namespace rtt::analysis {

namespace {

using compile::CompiledSystem;
using compile::PathId;
using detail::arrived;
using detail::check_path;
using detail::check_wavelength;
using detail::image_surface;
using detail::trace_rays;

math::Vec3 global_pos(const trace::RayBatch& rays, std::size_t i) {
  return {rays.pos_x()[i], rays.pos_y()[i], rays.pos_z()[i]};
}

math::Vec3 global_dir(const trace::RayBatch& rays, std::size_t i) {
  return {rays.dir_x()[i], rays.dir_y()[i], rays.dir_z()[i]};
}

/// Reference sphere and the data needed to measure OPL up to it (decided for #29).
struct Reference {
  ReferenceSphere sphere;
  math::Vec3 exit_pupil = math::Vec3::Zero();  ///< centre of the paraxial exit pupil, mm
  double lambda_mm = 0.0;                      ///< reference wavelength, mm
};

Reference make_reference(const CompiledSystem& system,
                         PathId path,
                         std::uint16_t field,
                         trace::Aiming aiming) {
  compile::require_stop(system, path);  // before any ray is traced (ADR 0022)
  const std::uint16_t ref = system.reference_wavelength();
  const std::uint32_t image = image_surface(system, path);
  const trace::RayBatch chief =
      trace_rays(system, path, field, ref, trace::SinglePupilPoint{0.0, 0.0}, aiming);
  if (!arrived(chief, 0, image)) {
    detail::throw_lost(system, chief, 0, field,
                       "analysis: the chief ray of field " + std::to_string(field) +
                           " does not reach the image surface");
  }
  const paraxial::FirstOrder fo = paraxial::first_order(system, path, ref);
  // With a stop on the path (require_stop above) first_order() always gives the exit pupil.
  if (!fo.exit_pupil) throw std::logic_error("analysis: no exit pupil despite a stop");
  const paraxial::Pupil& xp = *fo.exit_pupil;
  if (!xp.z) {
    throw AnalysisError(
        "analysis: the exit pupil is at infinity (image-space telecentric); OPD against a "
        "reference sphere is not supported for it yet");
  }
  Reference r;
  r.sphere.centre = global_pos(chief, 0);
  r.exit_pupil = math::Vec3(0.0, 0.0, *xp.z);
  r.sphere.radius = (r.sphere.centre - r.exit_pupil).norm();
  if (!(r.sphere.radius > 0.0)) {
    throw AnalysisError("analysis: the exit pupil lies on the image surface");
  }
  r.lambda_mm = system.wavelengths_um()[ref] * 1e-3;
  return r;
}

/// OPL of ray i up to the reference sphere: OPL at the image surface minus |n'| s, where s is
/// the signed path from the sphere to the image surface along the ray. Of the two crossings of
/// the ray line with the sphere the one closer to the exit-pupil centre is taken (the sphere
/// passes through it). None if the line misses the sphere.
std::optional<double> opl_to_sphere(const trace::RayBatch& rays,
                                    std::size_t i,
                                    const Reference& ref,
                                    double n_image) {
  const math::Vec3 p = global_pos(rays, i);
  const math::Vec3 d = global_dir(rays, i);
  const math::Vec3 pc = p - ref.sphere.centre;
  // |p - s d - C|^2 = R^2  ->  s^2 - 2 b s + |p - C|^2 - R^2 = 0 with b = d . (p - C).
  const double b = d.dot(pc);
  const double disc = b * b - (pc.squaredNorm() - ref.sphere.radius * ref.sphere.radius);
  if (!(disc >= 0.0)) return std::nullopt;
  const double root = std::sqrt(disc);
  const double s1 = b + root;
  const double s2 = b - root;
  const double e1 = (p - s1 * d - ref.exit_pupil).squaredNorm();
  const double e2 = (p - s2 * d - ref.exit_pupil).squaredNorm();
  const double s = e1 <= e2 ? s1 : s2;
  return rays.opl()[i] - n_image * s;
}

/// OPD of every ray of `sampling` at `wavelength`, in waves at the reference wavelength.
std::vector<OpdPoint> opd_points(const CompiledSystem& system,
                                 PathId path,
                                 std::uint16_t field,
                                 std::uint16_t wavelength,
                                 const trace::PupilSampling& sampling,
                                 trace::Aiming aiming,
                                 const Reference& ref) {
  const std::uint32_t image = image_surface(system, path);
  const auto& last = system.path(path).events.back();
  const double n_image = std::abs(system.media()[last.medium_after].index[wavelength].real());

  // Reference ray: chief ray of the same wavelength (W(0, 0) = 0, decided for #29).
  const trace::RayBatch chief =
      trace_rays(system, path, field, wavelength, trace::SinglePupilPoint{0.0, 0.0}, aiming);
  if (!arrived(chief, 0, image)) {
    detail::throw_lost(system, chief, 0, field,
                       "analysis: the chief ray of field " + std::to_string(field) +
                           " does not reach the image surface at wavelength " +
                           std::to_string(wavelength));
  }
  const std::optional<double> opl_chief = opl_to_sphere(chief, 0, ref, n_image);
  if (!opl_chief) throw AnalysisError("analysis: the chief ray misses the reference sphere");

  const trace::RayBatch rays = trace_rays(system, path, field, wavelength, sampling, aiming);
  std::vector<OpdPoint> points;
  points.reserve(rays.size());
  for (std::size_t i = 0; i < rays.size(); ++i) {
    OpdPoint q;
    q.px = rays.pupil_x()[i];
    q.py = rays.pupil_y()[i];
    std::optional<double> opl;
    if (arrived(rays, i, image)) opl = opl_to_sphere(rays, i, ref, n_image);
    if (opl) {
      // Wyant & Creath, Sec. I: W > 0 if the wavefront leads, W = OPL_ref - OPL_ray.
      q.w = (*opl_chief - *opl) / ref.lambda_mm;
    } else {
      q.status = rays.status()[i] == trace::RayStatus::Alive ? trace::RayStatus::Vignetted
                                                             : rays.status()[i];
    }
    points.push_back(q);
  }
  return points;
}

}  // namespace

OpdMap opd_map(const compile::CompiledSystem& system,
               compile::PathId path,
               std::uint16_t field,
               std::uint16_t wavelength,
               const OpdOptions& options) {
  check_path(system, path);
  check_wavelength(system, wavelength);
  if (options.grid < 1) throw std::invalid_argument("analysis: OPD grid must be >= 1");
  OpdMap map;
  map.field = field;
  map.wavelength = wavelength;
  const Reference ref = make_reference(system, path, field, options.aiming);
  map.sphere = ref.sphere;
  map.points = opd_points(system, path, field, wavelength, trace::GridPupil{options.grid},
                          options.aiming, ref);

  // Statistics over the arrived points, serial in grid order: mean, standard deviation
  // (piston removed, tilt kept), peak to valley.
  double sum = 0.0;
  double lo = 0.0;
  double hi = 0.0;
  for (const OpdPoint& q : map.points) {
    if (q.status != trace::RayStatus::Alive) {
      ++map.vignetted;
      continue;
    }
    lo = map.arrived == 0 ? q.w : std::min(lo, q.w);
    hi = map.arrived == 0 ? q.w : std::max(hi, q.w);
    sum += q.w;
    ++map.arrived;
  }
  if (map.arrived == 0) throw AnalysisError("analysis: no ray reaches the image surface");
  const double mean = sum / static_cast<double>(map.arrived);
  double var = 0.0;
  for (const OpdPoint& q : map.points) {
    if (q.status == trace::RayStatus::Alive) var += (q.w - mean) * (q.w - mean);
  }
  map.rms = std::sqrt(var / static_cast<double>(map.arrived));
  map.pv = hi - lo;
  return map;
}

OpdFan opd_fan(const compile::CompiledSystem& system,
               compile::PathId path,
               std::uint16_t field,
               std::uint16_t wavelength,
               const OpdOptions& options) {
  check_path(system, path);
  check_wavelength(system, wavelength);
  if (options.fan_points < 1) throw std::invalid_argument("analysis: fan points must be >= 1");
  OpdFan fan;
  fan.field = field;
  fan.wavelength = wavelength;
  const Reference ref = make_reference(system, path, field, options.aiming);
  fan.sphere = ref.sphere;
  fan.tangential = opd_points(system, path, field, wavelength, trace::FanYPupil{options.fan_points},
                              options.aiming, ref);
  fan.sagittal = opd_points(system, path, field, wavelength, trace::FanXPupil{options.fan_points},
                            options.aiming, ref);
  return fan;
}

}  // namespace rtt::analysis
