#pragma once

/// @file conic.hpp
/// Conic surface of revolution (sphere for k = 0) in local coordinates (lengths in mm).

#include <cmath>
#include <limits>
#include <optional>
#include <utility>

#include "rtt/geom/shape.hpp"

namespace rtt::geom {

/// Conic of revolution about the local z axis with vertex at the origin.
///
/// Sag (W. T. Welford, Aberrations of Optical Systems, Ch. 2):
///   z(r) = c r^2 / (1 + sqrt(1 - (1 + k) c^2 r^2)),  r^2 = x^2 + y^2.
/// k = 0 sphere, k = -1 paraboloid, k < -1 hyperboloid, -1 < k < 0 prolate and k > 0 oblate
/// ellipsoid. c = 1/R in 1/mm; c > 0 places the centre of curvature on the +z side.
template <rtt::math::Real T>
class Conic final : public Shape<T> {
 public:
  /// @param curvature      c = 1/R in 1/mm (sign convention see class comment)
  /// @param conic_constant k, dimensionless
  Conic(T curvature, T conic_constant) : c_(curvature), k_(conic_constant) {}

  /// Curvature c in 1/mm.
  [[nodiscard]] T curvature() const noexcept { return c_; }

  /// Conic constant k, dimensionless.
  [[nodiscard]] T conic_constant() const noexcept { return k_; }

  /// Sag z in mm at (x, y) in mm; NaN for r > max_radius().
  [[nodiscard]] T sag(T x, T y) const override;

  /// Slopes (dz/dx, dz/dy); NaN for r > max_radius(), +-infinity at r = max_radius().
  [[nodiscard]] std::pair<T, T> grad(T x, T y) const override;

  /// (c, k) of this conic.
  [[nodiscard]] std::pair<T, T> base_conic() const override { return {c_, k_}; }

  /// 1 / (|c| sqrt(1 + k)) in mm for c != 0 and k > -1, otherwise unbounded.
  [[nodiscard]] std::optional<T> max_radius() const override;

 private:
  T c_;
  T k_;
};

template <rtt::math::Real T>
T Conic<T>::sag(T x, T y) const {
  using std::sqrt;
  // Welford, Aberrations of Optical Systems, Ch. 2. The form with the square root in the
  // denominator stays finite for all k and avoids cancellation near the vertex.
  const T r2 = x * x + y * y;
  const T s2 = T(1) - (T(1) + k_) * c_ * c_ * r2;
  if (s2 < T(0)) {
    return std::numeric_limits<T>::quiet_NaN();
  }
  return c_ * r2 / (T(1) + sqrt(s2));
}

template <rtt::math::Real T>
std::pair<T, T> Conic<T>::grad(T x, T y) const {
  using std::sqrt;
  // Differentiating the sag above gives dz/dr = c r / sqrt(1 - (1 + k) c^2 r^2), hence
  // dz/dx = c x / s and dz/dy = c y / s with s = sqrt(1 - (1 + k) c^2 r^2).
  const T s2 = T(1) - (T(1) + k_) * c_ * c_ * (x * x + y * y);
  if (s2 < T(0)) {
    const T nan = std::numeric_limits<T>::quiet_NaN();
    return {nan, nan};
  }
  const T s = sqrt(s2);
  const T cx = c_ * x;
  const T cy = c_ * y;
  if (s == T(0)) {
    // Rim of a sphere or ellipsoid: the surface is vertical. Avoid 0/0 for a zero component.
    const T inf = std::numeric_limits<T>::infinity();
    const auto vertical = [inf](T v) { return v == T(0) ? T(0) : (v > T(0) ? inf : -inf); };
    return {vertical(cx), vertical(cy)};
  }
  return {cx / s, cy / s};
}

template <rtt::math::Real T>
std::optional<T> Conic<T>::max_radius() const {
  using std::abs;
  using std::sqrt;
  // Domain of the sag: (1 + k) c^2 r^2 <= 1.
  if (c_ == T(0) || T(1) + k_ <= T(0)) {
    return std::nullopt;
  }
  return T(1) / (abs(c_) * sqrt(T(1) + k_));
}

}  // namespace rtt::geom
