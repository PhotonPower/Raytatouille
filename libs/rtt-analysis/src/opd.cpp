#include "rtt/analysis/opd.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

#include "common.hpp"
#include "reference.hpp"

namespace rtt::analysis {

namespace {

using compile::CompiledSystem;
using compile::PathId;
using detail::check_path;
using detail::check_wavelength;
using detail::make_reference;
using detail::points_against;
using detail::Reference;

/// Through-focus shift of the options (ADR 0032, point 9: not finite is an input error).
void check_focus_shift(const OpdOptions& options) {
  if (!std::isfinite(options.focus_shift)) {
    throw std::invalid_argument("analysis: focus_shift must be finite");
  }
}

}  // namespace

OpdPupilPoints opd_points(const compile::CompiledSystem& system,
                          compile::PathId path,
                          std::uint16_t field,
                          std::uint16_t wavelength,
                          const trace::PupilSampling& sampling,
                          const OpdOptions& options) {
  return opd_points(system, path, field, wavelength, sampling, options, trace::RunControl{});
}

OpdPupilPoints opd_points(const compile::CompiledSystem& system,
                          compile::PathId path,
                          std::uint16_t field,
                          std::uint16_t wavelength,
                          const trace::PupilSampling& sampling,
                          const OpdOptions& options,
                          const trace::RunControl& control) {
  check_path(system, path);
  check_wavelength(system, wavelength);
  OpdPupilPoints result;
  result.field = field;
  result.wavelength = wavelength;
  check_focus_shift(options);
  const Reference ref = make_reference(system, path, field, options.aiming, options.focus_shift);
  result.sphere = ref.sphere;
  detail::LossCounter losses(system, path, options.lost_warning_fraction);
  result.points = points_against(system, path, field, wavelength, sampling, options.aiming, ref,
                                 losses, control);
  result.losses = losses.result();
  result.warnings = losses.warnings();
  return result;
}

OpdMap opd_map(const compile::CompiledSystem& system,
               compile::PathId path,
               std::uint16_t field,
               std::uint16_t wavelength,
               const OpdOptions& options) {
  return opd_map(system, path, field, wavelength, options, trace::RunControl{});
}

OpdMap opd_map(const compile::CompiledSystem& system,
               compile::PathId path,
               std::uint16_t field,
               std::uint16_t wavelength,
               const OpdOptions& options,
               const trace::RunControl& control) {
  check_path(system, path);
  check_wavelength(system, wavelength);
  if (options.grid < 1) throw std::invalid_argument("analysis: OPD grid must be >= 1");
  OpdMap map;
  map.field = field;
  map.wavelength = wavelength;
  check_focus_shift(options);
  const Reference ref = make_reference(system, path, field, options.aiming, options.focus_shift);
  map.sphere = ref.sphere;
  detail::LossCounter losses(system, path, options.lost_warning_fraction);
  map.points = points_against(system, path, field, wavelength, trace::GridPupil{options.grid},
                              options.aiming, ref, losses, control);
  map.losses = losses.result();
  map.warnings = losses.warnings();

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
  return opd_fan(system, path, field, wavelength, options, trace::RunControl{});
}

OpdFan opd_fan(const compile::CompiledSystem& system,
               compile::PathId path,
               std::uint16_t field,
               std::uint16_t wavelength,
               const OpdOptions& options,
               const trace::RunControl& control) {
  check_path(system, path);
  check_wavelength(system, wavelength);
  if (options.fan_points < 1) throw std::invalid_argument("analysis: fan points must be >= 1");
  OpdFan fan;
  fan.field = field;
  fan.wavelength = wavelength;
  check_focus_shift(options);
  const Reference ref = make_reference(system, path, field, options.aiming, options.focus_shift);
  fan.sphere = ref.sphere;
  detail::LossCounter losses(system, path, options.lost_warning_fraction);
  fan.tangential =
      points_against(system, path, field, wavelength, trace::FanYPupil{options.fan_points},
                     options.aiming, ref, losses, control);
  fan.sagittal =
      points_against(system, path, field, wavelength, trace::FanXPupil{options.fan_points},
                     options.aiming, ref, losses, control);
  fan.losses = losses.result();
  fan.warnings = losses.warnings();
  return fan;
}

}  // namespace rtt::analysis
