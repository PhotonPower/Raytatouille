#pragma once

/// @file thickness.hpp
/// Internal to rtt-analysis: the centre thickness of the dimension report (reports.cpp), split
/// out so that the tests reach it with vertices that are not at the origin (ADR 0041, point 6).
/// Not installed, not public API.

#include "rtt/math/isometry.hpp"
#include "rtt/math/types.hpp"

namespace rtt::analysis::detail {

/// Distance from the vertex of the first surface to the vertex of the second along the z axis of
/// the first, mm. A vertex is the local point (0, 0, sag(0, 0)) of its surface (ADR 0041,
/// point 6). For sag(0, 0) == 0 on both surfaces it is the distance of the origins bit for bit,
/// also the sign of a zero.
/// @param first_to_local   global -> local coordinates of the first surface
/// @param first_sag0       sag(0, 0) of the first surface, mm
/// @param second_to_global local coordinates of the second surface -> global
/// @param second_sag0      sag(0, 0) of the second surface, mm
[[nodiscard]] inline double centre_thickness(const math::Isometry3& first_to_local,
                                             double first_sag0,
                                             const math::Isometry3& second_to_global,
                                             double second_sag0) {
  const math::Vec3 vertex =
      second_sag0 == 0.0 ? math::Vec3::Zero() : math::Vec3(0.0, 0.0, second_sag0);
  const double z = first_to_local.apply_point(second_to_global.apply_point(vertex)).z();
  return first_sag0 == 0.0 ? z : z - first_sag0;
}

}  // namespace rtt::analysis::detail
