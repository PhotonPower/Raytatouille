#pragma once

/// @file odd_asphere.hpp
/// Odd asphere: conic of revolution plus a polynomial in r with every power from r^1 (lengths in
/// mm; ADR 0043, point 1).

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <optional>
#include <utility>
#include <vector>

#include "rtt/geom/conic.hpp"
#include "rtt/geom/shape.hpp"

namespace rtt::geom {

/// Odd asphere z(r) = c r^2 / (1 + phi) + A1 r + A2 r^2 + A3 r^3 + ...,
/// phi = sqrt(1 - (1 + k) c^2 r^2), r = sqrt(x^2 + y^2) (OpticStudio User Guide 2025 R1, "Odd
/// Asphere", there with beta_1 ... beta_8; conic part: G. W. Forbes, Opt. Express 19(10),
/// 9923-9942 (2011), Eqs. (2.1), (2.2); docs/quellen.md). Sign conventions of c and k as for
/// Conic; c = 0 is a plane base. The coefficients follow the model (rtt/model/surface.hpp):
/// coefficients[i] = A_(i + 1) in mm^(-i).
///
/// The gradient is (dz/dx, dz/dy) = (x, y) g(r) with g = c / phi + sum_{n >= 2} n A_n r^(n - 2)
/// + A1 / r; the term with A1 is formed as A1 (x / r) with r = hypot(x, y), which stays finite for
/// subnormal r, and is 0 at r = 0 (ADR 0043, point 1, F8): a cone has no gradient at its apex,
/// and 0 gives an axial ray the normal +z.
template <rtt::math::Real T>
class OddAsphere final {
 public:
  /// @param curvature      c of the base conic in 1/mm (0: plane base)
  /// @param conic_constant k of the base conic, dimensionless
  /// @param coefficients   A1, A2, A3, ... (coefficients[0] = A1), A_n in mm^(1 - n)
  OddAsphere(T curvature, T conic_constant, std::vector<T> coefficients)
      : conic_(curvature, conic_constant),
        coefficients_(std::move(coefficients)),
        has_polynomial_(std::any_of(
            coefficients_.begin(), coefficients_.end(), [](T a) { return a != T(0); })) {}

  /// Polynomial coefficients A1, A2, ... in mm^(1 - n).
  [[nodiscard]] const std::vector<T>& coefficients() const noexcept { return coefficients_; }

  /// Sag z in mm at (x, y) in mm; NaN for r > max_radius().
  [[nodiscard]] T sag(T x, T y) const noexcept;

  /// Slopes (dz/dx, dz/dy); NaN for r > max_radius(), +-infinity at r = max_radius().
  [[nodiscard]] std::pair<T, T> grad(T x, T y) const noexcept;

  /// (c, k) of the base conic: the start value of the intersection, not the vertex curvature
  /// (ADR 0043, point 3).
  [[nodiscard]] std::pair<T, T> base_conic() const { return conic_.base_conic(); }

  /// Domain of the base conic; the polynomial is defined everywhere.
  [[nodiscard]] std::optional<T> max_radius() const { return conic_.max_radius(); }

  /// Vertex expansion (ADR 0043, point 3; ADR 0041, point 1): value 0, slope 0, order 0 with
  /// A1 != 0 (|r| is not differentiable at the vertex), else 2 with A3 != 0 (|r|^3 has no third
  /// derivative), else 4; rotational. Hessian (c_v, 0, c_v) with the vertex curvature
  /// c_v = c + 2 A2 (z = (c / 2 + A2) r^2 + ...). Cap with the native K = k and
  /// A4_eff = A4 + (1 + k)(c^3 - c_v^3) / 8, so that K c_v^3 + 8 A4_eff = 8 b4 - c_v^3 with the r^4
  /// coefficient b4 = (1 + k) c^3 / 8 + A4. With A2 = 0 it is (k, A4) directly, bit for bit as
  /// EvenAsphere.
  [[nodiscard]] VertexExpansion<T> vertex_expansion() const noexcept;

 private:
  /// A_n, 0 beyond the list.
  [[nodiscard]] T coefficient(std::size_t n) const noexcept {
    return n >= 1 && n <= coefficients_.size() ? coefficients_[n - 1] : T(0);
  }

  Conic<T> conic_;
  std::vector<T> coefficients_;
  /// A coefficient is not zero: only then is the polynomial part added (ADR 0043, point 5.1),
  /// so that without one the shape is the conic bit for bit, the sign of a zero sag included.
  bool has_polynomial_;
};

template <rtt::math::Real T>
T OddAsphere<T>::sag(T x, T y) const noexcept {
  using std::hypot;
  const T base = conic_.sag(x, y);
  if (!has_polynomial_) return base;
  // sum_n A_n r^n = r (A1 + r (A2 + r (A3 + ...))), Horner's scheme in r.
  const T r = hypot(x, y);
  T acc = T(0);
  for (auto it = coefficients_.rbegin(); it != coefficients_.rend(); ++it) {
    acc = acc * r + *it;
  }
  return base + acc * r;
}

template <rtt::math::Real T>
std::pair<T, T> OddAsphere<T>::grad(T x, T y) const noexcept {
  using std::hypot;
  const auto [cx, cy] = conic_.grad(x, y);
  if (!has_polynomial_) return {cx, cy};
  // d/dx of A_n r^n is n A_n r^(n - 2) x: for n >= 2 the sum sum_n n A_n r^(n - 2) by Horner's
  // scheme in r, times x; the term n = 1 is A1 x / r.
  const T r = hypot(x, y);
  T acc = T(0);
  for (std::size_t n = coefficients_.size(); n >= 2; --n) {
    acc = acc * r + static_cast<T>(n) * coefficients_[n - 1];
  }
  T gx = cx + x * acc;
  T gy = cy + y * acc;
  if (const T a1 = coefficient(1); a1 != T(0) && r > T(0)) {
    gx += a1 * (x / r);
    gy += a1 * (y / r);
  }
  return {gx, gy};
}

template <rtt::math::Real T>
VertexExpansion<T> OddAsphere<T>::vertex_expansion() const noexcept {
  VertexExpansion<T> e = conic_.vertex_expansion();
  const auto [c, k] = conic_.base_conic();
  const T a2 = coefficient(2);
  const T a4 = coefficient(4);
  e.order = coefficient(1) != T(0) ? 0 : coefficient(3) != T(0) ? 2 : 4;
  if (a2 == T(0)) {
    // As EvenAsphere: the curvature of the conic and the cap (k, A4).
    if (coefficients_.size() >= 4) e.a4 = a4;
    return e;
  }
  const T cv = c + T(2) * a2;
  e.hessian = {cv, T(0), cv};
  e.conic = k;
  e.a4 = a4 + (T(1) + k) * (c * c * c - cv * cv * cv) / T(8);
  return e;
}

}  // namespace rtt::geom
