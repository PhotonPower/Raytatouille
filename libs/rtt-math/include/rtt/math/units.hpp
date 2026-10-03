#pragma once

/// @file units.hpp
/// Unit conventions of the whole code base (see docs/architecture.md, "Konventionen").
///
/// - Lengths: millimetre (mm)
/// - Wavelengths: micrometre (um)
/// - Angles: radian internally; API parameters in degrees carry the suffix `_deg`
/// - Temperature: degree Celsius; pressure: atm
/// - OPL: mm; OPD and wavefront: waves at the reference wavelength

#include <numbers>

namespace rtt::math {

inline constexpr double kPi = std::numbers::pi;

/// Converts degrees to radians.
[[nodiscard]] constexpr double deg_to_rad(double deg) noexcept {
  return deg * (kPi / 180.0);
}

/// Converts radians to degrees.
[[nodiscard]] constexpr double rad_to_deg(double rad) noexcept {
  return rad * (180.0 / kPi);
}

/// Converts a wavelength in micrometre to millimetre (the internal length unit).
[[nodiscard]] constexpr double um_to_mm(double um) noexcept {
  return um * 1.0e-3;
}

/// Converts a length in millimetre to micrometre.
[[nodiscard]] constexpr double mm_to_um(double mm) noexcept {
  return mm * 1.0e3;
}

}  // namespace rtt::math
