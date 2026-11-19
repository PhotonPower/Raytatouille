#pragma once

/// @file terms.hpp
/// Additive sag terms of a shape stack (ADR 0041, point 1; ADR 0051): Zernike terms after Noll
/// (orthonormal) and in the Fringe set (without normalisation), and the XY polynomial. A term is
/// evaluated in the local x, y of its surface, in mm, also beyond its normalisation radius
/// (ADR 0051, point 1).

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <optional>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

#include "rtt/geom/shape.hpp"
#include "rtt/math/zernike.hpp"

namespace rtt::geom {

/// Sum of Zernike polynomials z = sum_k a_k Z_k(rho, theta), rho = r / R, theta = atan2(y, x)
/// (from +x counter-clockwise), coefficients in mm (ADR 0051, points 2 and 3). Noll terms carry
/// their normalisation factor sqrt(n + 1) or sqrt(2 (n + 1)) (Niu/Tian Eq. (3)), Fringe terms none
/// (Eq. (22)); both are stored as (n, m, angular factor, a times factor).
template <rtt::math::Real T>
class ZernikeTerm final {
 public:
  /// One polynomial: R_n^m(rho) times cos(m theta) (or 1 for m = 0) or sin(m theta), with the
  /// coefficient in mm including the normalisation factor.
  struct Entry {
    int n = 0;
    int m = 0;
    bool sine = false;
    T coefficient = T(0);
  };

  /// Noll terms: `coefficients[k]` belongs to Z_{k+1}, in mm.
  /// @throws std::invalid_argument unless R > 0 is finite and there are at most 231 coefficients
  [[nodiscard]] static ZernikeTerm noll(T normalization_radius, std::span<const T> coefficients) {
    if (coefficients.size() > static_cast<std::size_t>(rtt::math::kMaxZernikeTerms)) {
      throw std::invalid_argument("Zernike: at most 231 Noll terms");
    }
    std::vector<Entry> entries;
    for (std::size_t k = 0; k < coefficients.size(); ++k) {
      const int j = static_cast<int>(k) + 1;
      const rtt::math::NollIndex nm = rtt::math::noll_to_nm(j);
      const T factor = T(rtt::math::zernike_noll_factor(nm.n, nm.m));
      entries.push_back({nm.n, nm.m, nm.m != 0 && j % 2 == 1, factor * coefficients[k]});
    }
    return ZernikeTerm(normalization_radius, std::move(entries), false);
  }

  /// Fringe terms: `coefficients[k]` belongs to Fringe term k + 1 (Niu/Tian Table 4, j = k), mm.
  /// @throws std::invalid_argument unless R > 0 is finite and there are at most 37 coefficients
  [[nodiscard]] static ZernikeTerm fringe(T normalization_radius, std::span<const T> coefficients) {
    if (coefficients.size() > static_cast<std::size_t>(rtt::math::kFringeTerms)) {
      throw std::invalid_argument("Zernike: at most 37 Fringe terms");
    }
    std::vector<Entry> entries;
    for (std::size_t k = 0; k < coefficients.size(); ++k) {
      const rtt::math::FringeIndex f = rtt::math::fringe_to_nm(static_cast<int>(k) + 1);
      entries.push_back({f.n, f.m, f.sine, coefficients[k]});
    }
    return ZernikeTerm(normalization_radius, std::move(entries), true);
  }

  /// z in mm at the local point (x, y) in mm.
  [[nodiscard]] T sag(T x, T y) const noexcept {
    const Polar p = polar(x, y);
    T z = T(0);
    for (const Entry& e : entries_) {
      const T r = T(rtt::math::zernike_radial_unchecked(e.n, e.m, p.rho));
      z += e.coefficient * r * angular(e, p).first;
    }
    return z;
  }

  /// (dz/dx, dz/dy), dimensionless. With rho = r/R and theta: dz/dx = (1/R) (R' Phi cos theta -
  /// (R/rho) Phi' sin theta), dz/dy = (1/R) (R' Phi sin theta + (R/rho) Phi' cos theta), Phi the
  /// angular factor (derived from the chain rule, ADR 0051, point 2). R/rho is evaluated as a
  /// polynomial, so the vertex needs no special case: there theta = 0 by convention.
  [[nodiscard]] std::pair<T, T> grad(T x, T y) const noexcept {
    const Polar p = polar(x, y);
    T gx = T(0);
    T gy = T(0);
    for (const Entry& e : entries_) {
      const auto [over, d] = rtt::math::zernike_radial_over_rho_and_derivative(e.n, e.m, p.rho);
      const auto [phi, dphi] = angular(e, p);
      const T a = e.coefficient * T(d) * phi;
      const T b = e.coefficient * T(over) * dphi;
      gx += a * p.cos - b * p.sin;
      gy += a * p.sin + b * p.cos;
    }
    return {gx / radius_, gy / radius_};
  }

  /// Unbounded: polynomials are defined everywhere.
  [[nodiscard]] std::optional<T> max_radius() const noexcept { return std::nullopt; }

  /// Taylor expansion at the vertex from the coefficients of rho^0 ... rho^4 of the radial
  /// polynomials (ADR 0051, point 2): value from m = 0 (rho^0), slope from m = 1 (rho^1),
  /// Hessian from m = 0 (rho^2: (x^2 + y^2)/R^2) and m = 2 (rho^2 cos 2 theta = (x^2 - y^2)/R^2,
  /// rho^2 sin 2 theta = 2 x y / R^2), b4 from m = 0 (rho^4). Rotational to fourth order exactly
  /// when every coefficient with 1 <= m <= 4 is 0 (terms with m >= 5 start at rho^5). Order 4;
  /// as a term K = 0 and A4 = b4 - c^3 / 8 (ADR 0041, point 1).
  [[nodiscard]] VertexExpansion<T> vertex_expansion() const {
    VertexExpansion<T> e;
    e.order = 4;
    e.rotational4 = true;
    const T r2 = radius_ * radius_;
    T b4 = T(0);
    T hxx = T(0);
    T hyy = T(0);
    for (const Entry& t : entries_) {
      if (t.m >= 1 && t.m <= 4 && t.coefficient != T(0)) e.rotational4 = false;
      const std::span<const double> c = rtt::math::zernike_radial_coefficients(t.n, t.m);
      const std::size_t k = c.size() - 1;  // K: c[K] belongs to rho^m, c[K - 1] to rho^(m + 2)
      if (t.m == 0) {
        e.value += t.coefficient * T(c[k]);
        if (k >= 1) {
          const T q = T(2) * t.coefficient * T(c[k - 1]) / r2;
          hxx += q;
          hyy += q;
        }
        if (k >= 2) b4 += t.coefficient * T(c[k - 2]) / (r2 * r2);
      } else if (t.m == 1) {
        const T s = t.coefficient * T(c[k]) / radius_;
        (t.sine ? e.slope.second : e.slope.first) += s;
      } else if (t.m == 2) {
        const T q = T(2) * t.coefficient * T(c[k]) / r2;
        if (t.sine) {
          e.hessian[1] += q;
        } else {
          hxx += q;
          hyy -= q;
        }
      }
    }
    e.hessian[0] = hxx;
    e.hessian[2] = e.rotational4 ? hxx : hyy;
    const T cc = e.hessian[0];
    e.a4 = b4 - cc * cc * cc / T(8);
    return e;
  }

  [[nodiscard]] T normalization_radius() const noexcept { return radius_; }
  [[nodiscard]] std::span<const Entry> entries() const noexcept { return entries_; }
  /// True for the Fringe set, false for Noll (for the layout data, ADR 0051, point 6).
  [[nodiscard]] bool fringe() const noexcept { return fringe_; }

 private:
  ZernikeTerm(T normalization_radius, std::vector<Entry> entries, bool fringe)
      : radius_(normalization_radius), entries_(std::move(entries)), fringe_(fringe) {
    if (!(radius_ > T(0)) || !std::isfinite(static_cast<double>(radius_))) {
      throw std::invalid_argument("Zernike: normalization radius must be finite and > 0 mm");
    }
    for (const Entry& e : entries_) max_m_ = std::max(max_m_, e.m);
  }

  struct Polar {
    double rho = 0.0;
    T cos = T(1);  ///< cos theta, 1 at the vertex by convention
    T sin = T(0);
    std::array<T, rtt::math::kMaxZernikeDegree + 1> cos_m{};  ///< cos(m theta)
    std::array<T, rtt::math::kMaxZernikeDegree + 1> sin_m{};  ///< sin(m theta)
  };

  /// rho and the multiples of theta by the recurrence of (cos + i sin)^m, without atan2.
  [[nodiscard]] Polar polar(T x, T y) const noexcept {
    using std::hypot;
    Polar p;
    const T r = hypot(x, y);
    p.rho = static_cast<double>(r / radius_);
    if (r > T(0)) {
      p.cos = x / r;
      p.sin = y / r;
    }
    p.cos_m[0] = T(1);
    p.sin_m[0] = T(0);
    for (int m = 1; m <= max_m_; ++m) {
      const auto i = static_cast<std::size_t>(m);
      p.cos_m[i] = p.cos_m[i - 1] * p.cos - p.sin_m[i - 1] * p.sin;
      p.sin_m[i] = p.sin_m[i - 1] * p.cos + p.cos_m[i - 1] * p.sin;
    }
    return p;
  }

  /// (Phi, dPhi/dtheta) of an entry: 1 and 0 for m = 0, cos(m theta) and -m sin(m theta), or
  /// sin(m theta) and m cos(m theta).
  [[nodiscard]] static std::pair<T, T> angular(const Entry& e, const Polar& p) noexcept {
    const auto i = static_cast<std::size_t>(e.m);
    if (e.m == 0) return {T(1), T(0)};
    if (e.sine) return {p.sin_m[i], T(e.m) * p.cos_m[i]};
    return {p.cos_m[i], -T(e.m) * p.sin_m[i]};
  }

  T radius_;
  std::vector<Entry> entries_;
  bool fringe_ = false;
  int max_m_ = 0;
};

/// Largest degree of the XY polynomial (ADR 0051, point 4).
inline constexpr int kMaxXyDegree = 20;
/// Number of XY polynomial terms up to degree 20: sum_{d=1}^{20} (d + 1) = 230 (ADR 0051, point 4).
inline constexpr int kMaxXyTerms = 230;

/// Exponents (i, j) of x^i y^j of XY polynomial term k = 1 ... 230 (ADR 0051, point 4; [UG]
/// "Extended Polynomial"): degree d = 1, 2, ... ascending, within a degree x^d, x^(d-1) y, ...,
/// y^d. So k = 1 is x, 2 is y, 3 is x^2, 4 is xy, 5 is y^2, 12 is x^2 y^2, 13 is x y^3.
/// @throws std::invalid_argument for k outside [1, kMaxXyTerms]
[[nodiscard]] inline std::pair<int, int> xy_polynomial_powers(int k) {
  if (k < 1 || k > kMaxXyTerms) {
    throw std::invalid_argument("XY polynomial: term must lie in [1, 230]");
  }
  // Terms before degree d: sum_{e=1}^{d-1} (e + 1) = (d - 1)(d + 2)/2.
  int d = 1;
  while (d * (d + 3) / 2 < k) ++d;              // d (d + 3)/2 terms up to degree d
  const int p = k - 1 - (d - 1) * (d + 2) / 2;  // position within degree d, 0 ... d
  return {d - p, p};
}

/// XY polynomial z = sum_k c_k (x/R)^i_k (y/R)^j_k without a constant term, coefficients in mm,
/// R the normalisation radius in mm; the order of xy_polynomial_powers() (ADR 0051, point 4; [UG]
/// "Extended Polynomial", the polynomial part without the base conic).
template <rtt::math::Real T>
class XyPolynomialTerm final {
 public:
  /// `coefficients[k]` belongs to term k + 1, in mm.
  /// @throws std::invalid_argument unless R > 0 is finite and there are at most 230 coefficients
  XyPolynomialTerm(T normalization_radius, std::vector<T> coefficients)
      : radius_(normalization_radius), c_(std::move(coefficients)) {
    if (!(radius_ > T(0)) || !std::isfinite(static_cast<double>(radius_))) {
      throw std::invalid_argument("XY polynomial: normalization radius must be finite and > 0 mm");
    }
    if (c_.size() > static_cast<std::size_t>(kMaxXyTerms)) {
      throw std::invalid_argument("XY polynomial: at most 230 coefficients");
    }
    for (std::size_t k = 0; k < c_.size(); ++k) {
      const auto [i, j] = xy_polynomial_powers(static_cast<int>(k) + 1);
      powers_.push_back({i, j});
      degree_ = std::max(degree_, i + j);
    }
  }

  /// z in mm at the local point (x, y) in mm.
  [[nodiscard]] T sag(T x, T y) const noexcept {
    const Powers p = powers(x, y);
    T z = T(0);
    for (std::size_t k = 0; k < c_.size(); ++k) {
      const auto [i, j] = powers_[k];
      z += c_[k] * p.u[static_cast<std::size_t>(i)] * p.v[static_cast<std::size_t>(j)];
    }
    return z;
  }

  /// (dz/dx, dz/dy), dimensionless: dz/dx = (1/R) sum_k c_k i u^(i-1) v^j and dz/dy = (1/R)
  /// sum_k c_k j u^i v^(j-1), u = x/R, v = y/R.
  [[nodiscard]] std::pair<T, T> grad(T x, T y) const noexcept {
    const Powers p = powers(x, y);
    T gx = T(0);
    T gy = T(0);
    for (std::size_t k = 0; k < c_.size(); ++k) {
      const auto [i, j] = powers_[k];
      const auto ui = static_cast<std::size_t>(i);
      const auto vj = static_cast<std::size_t>(j);
      if (i >= 1) gx += c_[k] * T(i) * p.u[ui - 1] * p.v[vj];
      if (j >= 1) gy += c_[k] * T(j) * p.u[ui] * p.v[vj - 1];
    }
    return {gx / radius_, gy / radius_};
  }

  /// Unbounded: a polynomial is defined everywhere.
  [[nodiscard]] std::optional<T> max_radius() const noexcept { return std::nullopt; }

  /// Taylor expansion at the vertex, exact (ADR 0051, point 4): value 0, slope (c_x, c_y)/R,
  /// Hessian (2 c_x2, c_xy, 2 c_y2)/R^2; b4 = c_x4 / R^4. Rotational to fourth order exactly when
  /// c_x = c_y = c_xy = 0, c_x2 = c_y2, every third-degree term 0, c_x4 = c_y4 = c_x2y2 / 2 and
  /// c_x3y = c_xy3 = 0, since (x^2 + y^2)^2 = x^4 + 2 x^2 y^2 + y^4. Order 4; as a term K = 0 and
  /// A4 = b4 - c^3 / 8 (ADR 0041, point 1).
  [[nodiscard]] VertexExpansion<T> vertex_expansion() const {
    // c[k]: coefficient of term k = 1 ... 14 (degree <= 4), 0 beyond the list.
    std::array<T, 15> c{};
    for (std::size_t k = 0; k < c_.size() && k < 14; ++k) c[k + 1] = c_[k];
    VertexExpansion<T> e;
    e.order = 4;
    const T r2 = radius_ * radius_;
    e.slope = {c[1] / radius_, c[2] / radius_};
    e.hessian[0] = T(2) * c[3] / r2;
    e.hessian[1] = c[4] / r2;
    e.hessian[2] = T(2) * c[5] / r2;
    e.rotational4 = c[1] == T(0) && c[2] == T(0) && c[4] == T(0) && c[3] == c[5] && c[6] == T(0) &&
                    c[7] == T(0) && c[8] == T(0) && c[9] == T(0) && c[10] == c[14] &&
                    c[12] == T(2) * c[10] && c[11] == T(0) && c[13] == T(0);
    const T b4 = c[10] / (r2 * r2);
    const T cc = e.hessian[0];
    e.a4 = b4 - cc * cc * cc / T(8);
    return e;
  }

  [[nodiscard]] T normalization_radius() const noexcept { return radius_; }
  /// c_1, c_2, ... in mm.
  [[nodiscard]] std::span<const T> coefficients() const noexcept { return c_; }

 private:
  struct Powers {
    std::array<T, kMaxXyDegree + 1> u{};  ///< (x/R)^i
    std::array<T, kMaxXyDegree + 1> v{};  ///< (y/R)^j
  };

  [[nodiscard]] Powers powers(T x, T y) const noexcept {
    Powers p;
    const T u = x / radius_;
    const T v = y / radius_;
    p.u[0] = T(1);
    p.v[0] = T(1);
    for (std::size_t i = 1; i <= static_cast<std::size_t>(degree_); ++i) {
      p.u[i] = p.u[i - 1] * u;
      p.v[i] = p.v[i - 1] * v;
    }
    return p;
  }

  T radius_;
  std::vector<T> c_;
  std::vector<std::pair<int, int>> powers_;
  int degree_ = 0;
};

}  // namespace rtt::geom
