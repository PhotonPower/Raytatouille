// STUB for the red run of the tests (#208): replaced by the implementation in the next commit.
#include "rtt/analysis/jones.hpp"

namespace rtt::analysis {

JonesPupil jones_pupil(const compile::CompiledSystem& system,
                       compile::PathId path,
                       std::uint16_t field,
                       std::uint16_t wavelength,
                       const JonesPupilOptions& options) {
  return jones_pupil(system, path, field, wavelength, options, trace::RunControl{});
}

JonesPupil jones_pupil(const compile::CompiledSystem& /*system*/,
                       compile::PathId /*path*/,
                       std::uint16_t /*field*/,
                       std::uint16_t /*wavelength*/,
                       const JonesPupilOptions& /*options*/,
                       const trace::RunControl& /*control*/) {
  return {};
}

std::array<math::Vec3, 2> input_chief_frame(const math::Vec3& /*k0_c*/) noexcept {
  return {math::Vec3::UnitX(), math::Vec3::UnitY()};
}

std::optional<std::array<math::Vec3, 2>> double_pole_basis(const math::Vec3& x_c,
                                                           const math::Vec3& y_c,
                                                           const math::Vec3& /*k_c*/,
                                                           const math::Vec3& /*k*/) {
  return std::array<math::Vec3, 2>{x_c, y_c};
}

std::optional<std::array<math::Vec3, 2>> dipole_basis(const math::Vec3& /*y_c*/,
                                                      const math::Vec3& /*k*/) {
  return std::nullopt;
}

PolarizationEllipse polarization_ellipse(std::complex<double> /*a*/,
                                         std::complex<double> /*b*/) noexcept {
  return {};
}

}  // namespace rtt::analysis
