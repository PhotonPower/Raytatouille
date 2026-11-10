#pragma once

/// @file frames.hpp
/// Frames of the ZMX import and export: the decomposition of a rotation into the tilts of a
/// coordinate break (ADR 0044, addendum PR C, rule 3; ADR 0045, addendum #366).

#include <array>
#include <cmath>

#include "rtt/math/types.hpp"
#include "rtt/math/units.hpp"

namespace rtt::io::detail {

/// Angles (a, b, c) in degrees with R = Rx(a) Ry(b) Rz(c) (math::rotation_xyz_intrinsic).
/// Derived in #317 (ADR 0044, addendum PR C) by multiplying out: row 0 of R is
/// (cos b cos c, -cos b sin c, sin b) and column 2 is (sin b, -sin a cos b, cos a cos b), so
/// b = atan2(R02, hypot(R00, R01)) and a = atan2(-R12, R22); row 1 of Rx(a)^T R = Ry(b) Rz(c)
/// is (sin c, cos c, 0), so c = atan2(cos a R10 + sin a R20, cos a R11 + sin a R21). That form
/// needs no case at b = +-90 deg, where only a +- c is defined: c takes up the noise in a.
inline std::array<double, 3> xyz_angles_deg(const math::Mat3& r) {
  const double b = std::atan2(r(0, 2), std::hypot(r(0, 0), r(0, 1)));
  const double a = std::atan2(-r(1, 2), r(2, 2));
  const double ca = std::cos(a);
  const double sa = std::sin(a);
  const double c = std::atan2(ca * r(1, 0) + sa * r(2, 0), ca * r(1, 1) + sa * r(2, 1));
  // + 0.0 turns a -0.0 into 0.0 for the file
  return {math::rad_to_deg(a) + 0.0, math::rad_to_deg(b) + 0.0, math::rad_to_deg(c) + 0.0};
}

}  // namespace rtt::io::detail
