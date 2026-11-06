#pragma once

/// @file axicon.hpp
/// Axicon: a right circular cone about the local z axis with its apex at the origin (lengths in
/// mm; ADR 0043, point 1). Its analytic intersection is in intersect.hpp.

#include <cmath>
#include <optional>
#include <utility>

#include "rtt/geom/shape.hpp"

namespace rtt::geom {

/// Axicon z(r) = r tan(alpha), r = sqrt(x^2 + y^2), with the base angle alpha between the mantle
/// and the plane z = 0 (derived in ADR 0043, point 1; docs/quellen.md). The sign follows that of
/// a curvature: alpha > 0 lets the mantle rise towards +z; alpha = 0 is the plane. It equals the
/// odd asphere with c = 0 and A1 = tan(alpha) (OddAsphere), which also covers a cone with a conic.
///
/// The normal off the apex is (-s x / r, -s y / r, 1) normalized with s = tan(alpha); at the apex
/// the gradient is (0, 0) and the normal +z (ADR 0043, F8).
template <rtt::math::Real T>
class Axicon final {
 public:
  /// @param angle base angle alpha in rad, |alpha| < pi/2 (validate() checks angle_deg)
  explicit Axicon(T angle) : angle_(angle), slope_(std::tan(angle)) {}

  /// Base angle alpha in rad.
  [[nodiscard]] T angle() const noexcept { return angle_; }

  /// Slope s = tan(alpha) of the mantle along r, dimensionless.
  [[nodiscard]] T slope() const noexcept { return slope_; }

  /// Sag z = s r in mm at (x, y) in mm; defined everywhere.
  [[nodiscard]] T sag(T x, T y) const noexcept {
    using std::hypot;
    return slope_ * hypot(x, y);
  }

  /// Slopes (s x / r, s y / r); (0, 0) at the apex r = 0.
  [[nodiscard]] std::pair<T, T> grad(T x, T y) const noexcept {
    using std::hypot;
    const T r = hypot(x, y);
    if (!(r > T(0))) return {T(0), T(0)};
    return {slope_ * (x / r), slope_ * (y / r)};
  }

  /// (0, 0): only for the concept IntersectableShape; the intersection is analytic (ADR 0043,
  /// point 1).
  [[nodiscard]] std::pair<T, T> base_conic() const { return {T(0), T(0)}; }

  /// The domain is unbounded.
  [[nodiscard]] std::optional<T> max_radius() const { return std::nullopt; }

  /// Vertex expansion (ADR 0043, point 3): order 0 for alpha != 0 (no paraxial model at the apex,
  /// F6); for alpha = 0 the flat expansion of the plane (order 4, Hessian 0, cap (0, 0)). value 0,
  /// slope 0 and rotational in both cases.
  [[nodiscard]] VertexExpansion<T> vertex_expansion() const noexcept {
    VertexExpansion<T> e;
    e.order = slope_ == T(0) ? 4 : 0;
    e.rotational4 = true;
    return e;
  }

 private:
  T angle_;
  T slope_;
};

}  // namespace rtt::geom
