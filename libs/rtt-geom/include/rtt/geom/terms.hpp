#pragma once

/// @file terms.hpp
/// Additive sag terms of a shape stack (ADR 0041, point 1; ADR 0051): Zernike terms after Noll
/// (orthonormal) and in the Fringe set (without normalisation), the XY polynomial and the grid
/// sag. A term is evaluated in the local x, y of its surface, in mm; the polynomial terms also
/// beyond their normalisation radius (ADR 0051, point 1), the grid sag only on its rectangle.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
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

/// Largest number of grid points per axis (ADR 0051, point 5).
inline constexpr int kMaxGridPoints = 4096;

/// Grid sag with bicubic Hermite interpolation (ADR 0051, point 5): nx by ny points with spacing
/// (dx, dy) in mm, centred on (x0, y0) in mm, x_i = x0 + (i - (nx - 1)/2) dx, y_j likewise, value k
/// = j nx + i (i along x, j along y, both ascending). In the cell [x_i, x_{i+1}] x [y_j, y_{j+1}]
/// with t = (x - x_i)/dx, u = (y - y_j)/dy:
///   p = sum_{a, b in {0, 1}} [H_a(t) H_b(u) z + dx G_a(t) H_b(u) z_x + dy H_a(t) G_b(u) z_y
///                             + dx dy G_a(t) G_b(u) z_xy] at (i + a, j + b),
/// H_0 = 2t^3 - 3t^2 + 1, H_1 = 3t^2 - 2t^3, G_0 = t^3 - 2t^2 + t, G_1 = t^3 - t^2 (the cubic
/// Hermite polynomial, Burden/Faires section 3.4 as in docs/quellen.md, as a tensor product). p is
/// C^1 and meets values and derivatives at the points. A point on an inner grid line belongs to
/// the cell to the right (above), the last line to the last cell. Outside the rectangle sag and
/// grad are NaN and in_domain is false (ADR 0042, point 1).
template <rtt::math::Real T>
class GridSagTerm final {
 public:
  /// @param nx, ny       points per axis, 3 ... 4096
  /// @param dx, dy       spacing in mm, finite and > 0
  /// @param x0, y0       centre of the grid in mm, finite
  /// @param sag          nx ny values in mm
  /// @param dz_dx, dz_dy nx ny slopes, dimensionless; d2z_dxdy nx ny values in 1/mm. The three
  ///                     lists together or none: without them the derivatives are estimated by
  ///                     finite differences of second order (estimate_derivatives()).
  /// @throws std::invalid_argument for sizes, spacings or list lengths out of range
  GridSagTerm(int nx,
              int ny,
              T dx,
              T dy,
              T x0,
              T y0,
              std::vector<T> sag,
              std::vector<T> dz_dx,
              std::vector<T> dz_dy,
              std::vector<T> d2z_dxdy)
      : nx_(nx),
        ny_(ny),
        dx_(dx),
        dy_(dy),
        x0_(x0),
        y0_(y0),
        z_(std::move(sag)),
        zx_(std::move(dz_dx)),
        zy_(std::move(dz_dy)),
        zxy_(std::move(d2z_dxdy)) {
    if (nx_ < 3 || ny_ < 3 || nx_ > kMaxGridPoints || ny_ > kMaxGridPoints) {
      throw std::invalid_argument("grid sag: nx and ny must lie in [3, 4096]");
    }
    if (!(dx_ > T(0)) || !(dy_ > T(0)) || !std::isfinite(static_cast<double>(dx_)) ||
        !std::isfinite(static_cast<double>(dy_))) {
      throw std::invalid_argument("grid sag: spacing must be finite and > 0 mm");
    }
    const std::size_t n = static_cast<std::size_t>(nx_) * static_cast<std::size_t>(ny_);
    if (z_.size() != n) throw std::invalid_argument("grid sag: sag needs nx ny values");
    derivatives_given_ = !zx_.empty() || !zy_.empty() || !zxy_.empty();
    if (derivatives_given_) {
      if (zx_.size() != n || zy_.size() != n || zxy_.size() != n) {
        throw std::invalid_argument("grid sag: the three derivative lists need nx ny values each");
      }
    } else {
      estimate_derivatives();
    }
    x_min_ = node(x0_, dx_, nx_, 0);
    x_max_ = node(x0_, dx_, nx_, nx_ - 1);
    y_min_ = node(y0_, dy_, ny_, 0);
    y_max_ = node(y0_, dy_, ny_, ny_ - 1);
  }

  /// z in mm at the local point (x, y) in mm; NaN outside the grid.
  [[nodiscard]] T sag(T x, T y) const noexcept {
    if (!in_domain(x, y)) return std::numeric_limits<T>::quiet_NaN();
    const Cell c = cell(x, y);
    const std::array<T, 2> ht = h(c.t);
    const std::array<T, 2> gt = g(c.t);
    const std::array<T, 2> hu = h(c.u);
    const std::array<T, 2> gu = g(c.u);
    T z = T(0);
    for (std::size_t a = 0; a < 2; ++a) {
      for (std::size_t b = 0; b < 2; ++b) {
        const std::size_t k = c.k + b * static_cast<std::size_t>(nx_) + a;
        z += ht[a] * hu[b] * z_[k] + dx_ * gt[a] * hu[b] * zx_[k] + dy_ * ht[a] * gu[b] * zy_[k] +
             dx_ * dy_ * gt[a] * gu[b] * zxy_[k];
      }
    }
    return z;
  }

  /// (dz/dx, dz/dy), dimensionless, the derivative of the interpolant; NaN outside the grid.
  [[nodiscard]] std::pair<T, T> grad(T x, T y) const noexcept {
    if (!in_domain(x, y)) {
      const T nan = std::numeric_limits<T>::quiet_NaN();
      return {nan, nan};
    }
    const Cell c = cell(x, y);
    const std::array<T, 2> ht = h(c.t);
    const std::array<T, 2> gt = g(c.t);
    const std::array<T, 2> hu = h(c.u);
    const std::array<T, 2> gu = g(c.u);
    const std::array<T, 2> dht = dh(c.t);
    const std::array<T, 2> dgt = dg(c.t);
    const std::array<T, 2> dhu = dh(c.u);
    const std::array<T, 2> dgu = dg(c.u);
    T px = T(0);  // dp/dt
    T py = T(0);  // dp/du
    for (std::size_t a = 0; a < 2; ++a) {
      for (std::size_t b = 0; b < 2; ++b) {
        const std::size_t k = c.k + b * static_cast<std::size_t>(nx_) + a;
        px += dht[a] * hu[b] * z_[k] + dx_ * dgt[a] * hu[b] * zx_[k] +
              dy_ * dht[a] * gu[b] * zy_[k] + dx_ * dy_ * dgt[a] * gu[b] * zxy_[k];
        py += ht[a] * dhu[b] * z_[k] + dx_ * gt[a] * dhu[b] * zx_[k] +
              dy_ * ht[a] * dgu[b] * zy_[k] + dx_ * dy_ * gt[a] * dgu[b] * zxy_[k];
      }
    }
    return {px / dx_, py / dy_};
  }

  /// The largest distance of a grid corner from the origin, mm (ADR 0051, point 5).
  [[nodiscard]] std::optional<T> max_radius() const noexcept {
    using std::hypot;
    using std::max;
    return max(max(hypot(x_min_, y_min_), hypot(x_min_, y_max_)),
               max(hypot(x_max_, y_min_), hypot(x_max_, y_max_)));
  }

  /// True on the closed rectangle [x_0, x_{nx-1}] x [y_0, y_{ny-1}] (false for NaN).
  [[nodiscard]] bool in_domain(T x, T y) const noexcept {
    return x >= x_min_ && x <= x_max_ && y >= y_min_ && y <= y_max_;
  }

  /// Order 1 (the interpolant is only C^1, ADR 0051, point 5): value and slope at the vertex,
  /// no Hessian; paraxially "no model at the vertex" (ADR 0041, point 2). NaN if the origin lies
  /// outside the grid (compile reports shape.grid_vertex_outside before).
  [[nodiscard]] VertexExpansion<T> vertex_expansion() const noexcept {
    VertexExpansion<T> e;
    e.order = 1;
    e.value = sag(T(0), T(0));
    e.slope = grad(T(0), T(0));
    return e;
  }

  [[nodiscard]] int nx() const noexcept { return nx_; }
  [[nodiscard]] int ny() const noexcept { return ny_; }
  [[nodiscard]] T spacing_x() const noexcept { return dx_; }
  [[nodiscard]] T spacing_y() const noexcept { return dy_; }
  [[nodiscard]] T decenter_x() const noexcept { return x0_; }
  [[nodiscard]] T decenter_y() const noexcept { return y0_; }
  /// True if the file gave the derivatives, false if they were estimated.
  [[nodiscard]] bool derivatives_given() const noexcept { return derivatives_given_; }
  [[nodiscard]] std::span<const T> values() const noexcept { return z_; }
  [[nodiscard]] std::span<const T> dz_dx() const noexcept { return zx_; }
  [[nodiscard]] std::span<const T> dz_dy() const noexcept { return zy_; }
  [[nodiscard]] std::span<const T> d2z_dxdy() const noexcept { return zxy_; }

 private:
  struct Cell {
    std::size_t k = 0;  ///< index of the lower left point (i, j)
    T t = T(0);         ///< (x - x_i) / dx in [0, 1]
    T u = T(0);         ///< (y - y_j) / dy in [0, 1]
  };

  /// x_i = x0 + (i - (n - 1)/2) d.
  [[nodiscard]] static T node(T centre, T d, int n, int i) noexcept {
    return centre + (T(i) - T(n - 1) / T(2)) * d;
  }

  /// The cell of the point (inside the rectangle) and t, u.
  [[nodiscard]] Cell cell(T x, T y) const noexcept {
    using std::floor;
    const auto locate = [](T v, T lo, T d, int n, std::size_t& index, T& frac) {
      const T s = (v - lo) / d;
      int i = static_cast<int>(floor(s));
      i = std::clamp(i, 0, n - 2);  // the last grid line belongs to the last cell
      index = static_cast<std::size_t>(i);
      frac = s - T(i);
    };
    Cell c;
    std::size_t i = 0;
    std::size_t j = 0;
    locate(x, x_min_, dx_, nx_, i, c.t);
    locate(y, y_min_, dy_, ny_, j, c.u);
    c.k = j * static_cast<std::size_t>(nx_) + i;
    return c;
  }

  /// (H_0, H_1) and (G_0, G_1) at t, and their derivatives in t.
  [[nodiscard]] static std::array<T, 2> h(T t) noexcept {
    const T h1 = t * t * (T(3) - T(2) * t);
    return {T(1) - h1, h1};
  }
  [[nodiscard]] static std::array<T, 2> g(T t) noexcept {
    const T s = T(1) - t;
    return {t * s * s, -t * t * s};
  }
  [[nodiscard]] static std::array<T, 2> dh(T t) noexcept {
    const T d = T(6) * t * (T(1) - t);  // H_1' = 6t - 6t^2, H_0' = -H_1'
    return {-d, d};
  }
  [[nodiscard]] static std::array<T, 2> dg(T t) noexcept {
    return {(T(1) - t) * (T(1) - T(3) * t), t * (T(3) * t - T(2))};
  }

  /// Finite differences of second order (ADR 0051, point 5): central inside,
  /// (-3 z_0 + 4 z_1 - z_2) / (2 d) and (3 z_{n-1} - 4 z_{n-2} + z_{n-3}) / (2 d) at the edges;
  /// z_x along x, z_y along y, z_xy as the y difference of z_x.
  void estimate_derivatives() {
    const auto nx = static_cast<std::size_t>(nx_);
    const auto ny = static_cast<std::size_t>(ny_);
    // The difference along a line of n values with stride `stride` from `first`.
    const auto diff = [](const std::vector<T>& v, std::size_t first, std::size_t stride,
                         std::size_t n, std::size_t i, T d) {
      const auto at = [&](std::size_t m) { return v[first + m * stride]; };
      if (i == 0) return (T(-3) * at(0) + T(4) * at(1) - at(2)) / (T(2) * d);
      if (i == n - 1) return (T(3) * at(n - 1) - T(4) * at(n - 2) + at(n - 3)) / (T(2) * d);
      return (at(i + 1) - at(i - 1)) / (T(2) * d);
    };
    zx_.assign(nx * ny, T(0));
    zy_.assign(nx * ny, T(0));
    zxy_.assign(nx * ny, T(0));
    for (std::size_t j = 0; j < ny; ++j) {
      for (std::size_t i = 0; i < nx; ++i) zx_[j * nx + i] = diff(z_, j * nx, 1, nx, i, dx_);
    }
    for (std::size_t i = 0; i < nx; ++i) {
      for (std::size_t j = 0; j < ny; ++j) {
        zy_[j * nx + i] = diff(z_, i, nx, ny, j, dy_);
        zxy_[j * nx + i] = diff(zx_, i, nx, ny, j, dy_);
      }
    }
  }

  int nx_;
  int ny_;
  T dx_;
  T dy_;
  T x0_;
  T y0_;
  std::vector<T> z_;
  std::vector<T> zx_;
  std::vector<T> zy_;
  std::vector<T> zxy_;
  bool derivatives_given_ = false;
  T x_min_ = T(0);
  T x_max_ = T(0);
  T y_min_ = T(0);
  T y_max_ = T(0);
};

}  // namespace rtt::geom
