#pragma once

/// @file intersect.hpp
/// Analytic ray-surface intersection in the local coordinate system of a surface.

#include <cmath>
#include <cstdint>

#include "rtt/geom/conic.hpp"
#include "rtt/geom/plane.hpp"
#include "rtt/math/types.hpp"

namespace rtt::geom {

/// Outcome of an intersection. The intersection never throws and reports problems as a status
/// instead of NaN (ADR 0009: status flags in tracing loops).
enum class HitStatus : std::uint8_t {
  Hit,     ///< t, point and normal are valid
  Missed,  ///< no intersection with t > t_min on the surface; t, point and normal are zero
};

/// Result of a ray-surface intersection in local coordinates.
template <rtt::math::Real T>
struct Intersection {
  HitStatus status = HitStatus::Missed;
  T t = T(0);  ///< ray parameter in mm: point = origin + t * direction
  rtt::math::Vec3T<T> point = rtt::math::Vec3T<T>::Zero();   ///< hit point, mm
  rtt::math::Vec3T<T> normal = rtt::math::Vec3T<T>::Zero();  ///< unit normal, +z at the vertex
};

/// Default lower bound for the ray parameter in mm: solutions with t <= kDefaultTMin are
/// ignored so that a ray starting on a surface does not hit it again at t = 0.
inline constexpr double kDefaultTMin = 1e-9;

/// Intersection of a ray with the conic z = c r^2 / (1 + sqrt(1 - (1 + k) c^2 r^2)).
/// @param c         curvature in 1/mm, positive if the centre of curvature is on the +z side
/// @param k         conic constant, dimensionless
/// @param origin    ray origin in local coordinates, mm
/// @param direction direction in local coordinates; must have length 1 so that t is in mm
/// @param t_min     smallest accepted ray parameter in mm
/// @return nearest hit with t > t_min on the sheet through the vertex, or Missed.
template <rtt::math::Real T>
[[nodiscard]] Intersection<T> intersect_conic(T c,
                                              T k,
                                              const rtt::math::Vec3T<T>& origin,
                                              const rtt::math::Vec3T<T>& direction,
                                              T t_min = T(kDefaultTMin)) noexcept {
  using std::sqrt;
  // Rearranging the sag z = c r^2 / (1 + s), s = sqrt(1 - (1 + k) c^2 r^2) (Welford,
  // Aberrations of Optical Systems, Ch. 2) gives the quadric
  //   F(x, y, z) = c (x^2 + y^2) + (1 + k) c z^2 - 2 z = 0,
  // whose solutions are z = (1 -+ s) / ((1 + k) c). The sag is the "-" root, i.e. the sheet
  // through the vertex, on which w = 1 - (1 + k) c z = s >= 0.
  // Inserting p = o + t d yields a t^2 + 2 b t + g = 0 with the coefficients below.
  const T q1 = T(1) + k;
  const T a = c * (direction.x() * direction.x() + direction.y() * direction.y() +
                   q1 * direction.z() * direction.z());
  const T b = c * (origin.x() * direction.x() + origin.y() * direction.y() +
                   q1 * origin.z() * direction.z()) -
              direction.z();
  const T g =
      c * (origin.x() * origin.x() + origin.y() * origin.y() + q1 * origin.z() * origin.z()) -
      T(2) * origin.z();
  const T disc = b * b - a * g;
  if (!(disc >= T(0))) {  // also catches NaN input
    return {};
  }
  // Roots without cancellation: q = -(b + sign(b) sqrt(disc)), t = q / a and t = g / q
  // (Press et al., Numerical Recipes, 3rd ed., Sec. 5.6, Eqs. 5.6.4-5.6.5). For a = 0 (plane,
  // paraboloid or hyperboloid along an asymptote) only g / q remains, the linear solution.
  const T q = -(b + (b < T(0) ? -sqrt(disc) : sqrt(disc)));
  Intersection<T> best;
  const auto consider = [&](T t) {
    if (!(t > t_min) || (best.status == HitStatus::Hit && t >= best.t)) {
      return;
    }
    const rtt::math::Vec3T<T> p = origin + t * direction;
    const T w = T(1) - q1 * c * p.z();
    if (!(w >= T(0)) || !p.allFinite()) {
      return;  // other sheet of the quadric, not part of the surface
    }
    // grad F = 2 (c x, c y, -w); with w = s this is -2 s (-dz/dx, -dz/dy, 1), so the normal
    // equals (-dz/dx, -dz/dy, 1) normalised (docs/architecture.md) but stays finite at the rim.
    const rtt::math::Vec3T<T> n(-c * p.x(), -c * p.y(), w);
    best.status = HitStatus::Hit;
    best.t = t;
    best.point = p;
    best.normal = n.normalized();
  };
  if (a != T(0)) {
    consider(q / a);
  }
  if (q != T(0)) {
    consider(g / q);
  }
  return best;
}

/// Intersection of a ray with a conic, see intersect_conic() for units and conventions.
template <rtt::math::Real T>
[[nodiscard]] Intersection<T> intersect(const Conic<T>& conic,
                                        const rtt::math::Vec3T<T>& origin,
                                        const rtt::math::Vec3T<T>& direction,
                                        T t_min = T(kDefaultTMin)) noexcept {
  return intersect_conic(conic.curvature(), conic.conic_constant(), origin, direction, t_min);
}

/// Intersection of a ray with the plane z = 0. Normal is always +z.
/// @param origin    ray origin in local coordinates, mm
/// @param direction direction in local coordinates; must have length 1 so that t is in mm
/// @param t_min     smallest accepted ray parameter in mm
template <rtt::math::Real T>
[[nodiscard]] Intersection<T> intersect(const Plane<T>& /*plane*/,
                                        const rtt::math::Vec3T<T>& origin,
                                        const rtt::math::Vec3T<T>& direction,
                                        T t_min = T(kDefaultTMin)) noexcept {
  // o_z + t d_z = 0. A ray parallel to the plane (d_z = 0) never reaches it.
  if (direction.z() == T(0)) {
    return {};
  }
  const T t = -origin.z() / direction.z();
  if (!(t > t_min)) {
    return {};
  }
  Intersection<T> hit;
  hit.point = origin + t * direction;
  hit.point.z() = T(0);
  if (!hit.point.allFinite()) {
    return {};
  }
  hit.status = HitStatus::Hit;
  hit.t = t;
  hit.normal = rtt::math::Vec3T<T>::UnitZ();
  return hit;
}

}  // namespace rtt::geom
