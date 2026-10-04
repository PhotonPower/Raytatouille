#pragma once

/// @file intersect.hpp
/// Ray-surface intersection in the local coordinate system of a surface: analytic for plane and
/// conic, Newton iteration for general shapes.

#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>

#include "rtt/geom/conic.hpp"
#include "rtt/geom/plane.hpp"
#include "rtt/geom/shape.hpp"
#include "rtt/math/types.hpp"

namespace rtt::geom {

/// Outcome of an intersection. The intersection never throws and reports problems as a status
/// instead of NaN (ADR 0009: status flags in tracing loops).
enum class HitStatus : std::uint8_t {
  Hit,            ///< t, point and normal are valid
  Missed,         ///< no intersection with t > t_min on the surface; t, point, normal are zero
  NoConvergence,  ///< Newton iteration did not converge; t, point and normal are zero
};

/// Result of a ray-surface intersection in local coordinates.
template <rtt::math::Real T>
struct Intersection {
  HitStatus status = HitStatus::Missed;
  T t = T(0);  ///< ray parameter in mm: point = origin + t * direction
  rtt::math::Vec3T<T> point = rtt::math::Vec3T<T>::Zero();   ///< hit point, mm
  rtt::math::Vec3T<T> normal = rtt::math::Vec3T<T>::Zero();  ///< unit normal, +z at the vertex
  int iterations = 0;  ///< Newton steps taken (0 for analytic intersections), for diagnosis
};

/// Default lower bound for the ray parameter in mm: solutions with t <= kDefaultTMin are
/// ignored so that a ray starting on a surface does not hit it again at t = 0.
inline constexpr double kDefaultTMin = 1e-9;

/// Newton iteration has converged when |F(t)| = |z(t) - sag(x(t), y(t))| <= kNewtonTolerance,
/// in mm (docs/architecture.md, Physik-Module, Geometrie), or when the last step reached the
/// rounding limit |dt| <= kNewtonStepFactor * eps * (1 + |t| + |origin|) with eps the machine
/// epsilon. The second criterion covers far ray origins (|origin| ~ 1e5 mm), for which the
/// rounding error of origin + t * direction alone exceeds kNewtonTolerance. Stopping on the
/// step size: Press et al., Numerical Recipes, 3rd ed., Sec. 9.4.
inline constexpr double kNewtonTolerance = 1e-12;

/// Factor of the rounding-limit step criterion, see kNewtonTolerance. Estimate: evaluating F
/// at o + t d has a rounding error of about eps (|o| + |t|) (1 + |grad sag|), so the Newton step
/// carries noise of about that error / |F'| (|direction| = 1). The criterion assumes
/// (1 + |grad sag|) / |F'| = O(1) and the factor 8 gives margin for that; at grazing incidence
/// or on steep flanks the step may stay above the limit and the result is NoConvergence.
inline constexpr double kNewtonStepFactor = 8.0;

/// Maximum number of Newton steps before the status is NoConvergence (docs/architecture.md).
inline constexpr int kMaxNewtonIterations = 30;

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
  // Rearranging the sag z = c r^2 / (1 + s), s = sqrt(1 - (1 + k) c^2 r^2) (Forbes, Opt.
  // Express 19(10), 9923-9942 (2011), Eq. (2.1), s = phi there; see also Welford, Aberrations of
  // Optical Systems, Ch. 2) gives the quadric
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

/// General intersection with any shape: start at the analytic hit of the base conic (or of the
/// vertex plane if the base conic is missed), then Newton on F(t) = z(t) - sag(x(t), y(t)).
/// @param shape     surface shape in local coordinates
/// @param origin    ray origin in local coordinates, mm
/// @param direction direction in local coordinates; must have length 1 so that t is in mm
/// @param t_min     smallest accepted ray parameter in mm
/// @return Hit (converged as described at kNewtonTolerance) with the number of Newton steps;
///         Missed if no start point exists or the start lies outside the shape's domain;
///         NoConvergence after
///         kMaxNewtonIterations steps, if an iterate leaves the domain, if F or the gradient is
///         not finite, if F'(t) = 0, or if Newton converges onto a crossing with t <= t_min.
template <rtt::math::Real T>
[[nodiscard]] Intersection<T> intersect(const Shape<T>& shape,
                                        const rtt::math::Vec3T<T>& origin,
                                        const rtt::math::Vec3T<T>& direction,
                                        T t_min = T(kDefaultTMin)) {
  using std::abs;
  using std::isfinite;
  using Vec = rtt::math::Vec3T<T>;
  // Start value: analytic hit of the base conic (docs/architecture.md); if that is missed,
  // the hit of the vertex plane z = 0.
  const auto [c, k] = shape.base_conic();
  T t = T(0);
  if (const Intersection<T> start = intersect_conic(c, k, origin, direction, t_min);
      start.status == HitStatus::Hit) {
    t = start.t;
  } else if (const Intersection<T> plane = intersect(Plane<T>(), origin, direction, t_min);
             plane.status == HitStatus::Hit) {
    t = plane.t;
  } else {
    return {};
  }
  // sag() and grad() are only evaluated inside the domain of the shape (NaN outside).
  const std::optional<T> max_radius = shape.max_radius();
  const auto in_domain = [&max_radius](const Vec& p) {
    return p.allFinite() &&
           (!max_radius || p.x() * p.x() + p.y() * p.y() <= *max_radius * *max_radius);
  };
  Vec p = origin + t * direction;
  if (!in_domain(p)) {
    return {};  // no surface at the start point: the ray passes outside the shape
  }
  Intersection<T> failed;
  failed.status = HitStatus::NoConvergence;
  const T origin_norm = origin.norm();
  T last_step = std::numeric_limits<T>::infinity();
  // Newton's method on F(t) = z(t) - sag(x(t), y(t)) with
  // F'(t) = d_z - dz/dx d_x - dz/dy d_y (Press et al., Numerical Recipes, 3rd ed., Sec. 9.4).
  for (int iteration = 0;; ++iteration) {
    failed.iterations = iteration;
    const T f = p.z() - shape.sag(p.x(), p.y());
    const auto [gx, gy] = shape.grad(p.x(), p.y());
    if (!isfinite(f) || !isfinite(gx) || !isfinite(gy)) {
      return failed;
    }
    const T step_limit =
        T(kNewtonStepFactor) * std::numeric_limits<T>::epsilon() * (T(1) + abs(t) + origin_norm);
    if (abs(f) <= T(kNewtonTolerance) || abs(last_step) <= step_limit) {
      if (!(t > t_min)) {
        return failed;  // converged behind the ray; a crossing ahead is not ruled out
      }
      // Normal from the gradient: (-dz/dx, -dz/dy, 1) normalised (docs/architecture.md).
      Intersection<T> hit;
      hit.status = HitStatus::Hit;
      hit.t = t;
      hit.point = p;
      hit.normal = Vec(-gx, -gy, T(1)).normalized();
      hit.iterations = iteration;
      return hit;
    }
    if (iteration == kMaxNewtonIterations) {
      return failed;
    }
    const T df = direction.z() - gx * direction.x() - gy * direction.y();
    if (!(df != T(0))) {
      return failed;  // ray tangent to the surface
    }
    last_step = f / df;
    t -= last_step;
    p = origin + t * direction;
    if (!in_domain(p)) {
      failed.iterations = iteration + 1;
      return failed;  // iteration left the domain of the shape
    }
  }
}

}  // namespace rtt::geom
