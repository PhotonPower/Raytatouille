#pragma once

/// @file forbes.hpp
/// Forbes Q-con and Q-bfs aspheres (lengths in mm; ADR 0043, point 1) and the sums of their
/// polynomials. Sources (docs/quellen.md):
/// - [F07] G. W. Forbes, "Shape specification for axially symmetric optical surfaces", Opt.
///   Express 15(8), 5218-5226 (2007);
/// - [F10a] G. W. Forbes, "Robust and fast computation for the polynomials of optics", Opt.
///   Express 18(13), 13851-13862 (2010);
/// - [F10b] G. W. Forbes, "Robust, efficient computational methods for axially symmetric optical
///   aspheres", Opt. Express 18(19), 19700-19712 (2010);
/// - [F11] G. W. Forbes, "Manufacturability estimates for optical aspheres", Opt. Express
///   19(10), 9923-9942 (2011).
/// The normalized radial coordinate is u = r / rho with the normalization radius rho ([F07]
/// rho_max, half the clear aperture); beyond rho (u > 1) the polynomials are extrapolated (ADR
/// 0043, point 2, F3).

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <span>
#include <utility>
#include <vector>

#include "rtt/geom/conic.hpp"
#include "rtt/geom/shape.hpp"

namespace rtt::geom {

/// Coefficients of the three-term recurrence of Q-con, Q_(n+1) = (a_n + b_n x) Q_n - c_n Q_(n-1)
/// ([F10a], Eq. (5.2a,b,c)): a_n = -(2n + 5)(n^2 + 5n + 10) / ((n + 1)(n + 2)(n + 5)),
/// b_n = 2(n + 3)(2n + 5) / ((n + 1)(n + 5)), c_n = (n + 3)(n + 4) n / ((n + 1)(n + 2)(n + 5)).
/// For n = 0 they give Q_1 = 6x - 5 from Q_0 = 1 ([F07], Eq. (9)). The products are formed in T
/// (exact while below 2^53 for double, n up to about 1.6e5), each coefficient is then one rounded
/// division; no integer overflow for any n.
template <rtt::math::Real T>
struct QConRecurrence {
  T a;
  T b;
  T c;
};

/// The coefficients of [F10a], Eq. (5.2a,b,c) for n >= 0.
template <rtt::math::Real T>
[[nodiscard]] QConRecurrence<T> qcon_recurrence(std::size_t n) noexcept {
  const T m = static_cast<T>(n);
  const T den = (m + T(1)) * (m + T(2)) * (m + T(5));
  return {-(T(2) * m + T(5)) * (m * m + T(5) * m + T(10)) / den,
          T(2) * (m + T(3)) * (T(2) * m + T(5)) / ((m + T(1)) * (m + T(5))),
          (m + T(3)) * (m + T(4)) * m / den};
}

namespace detail {

/// Clenshaw's sum with the recurrence `rec` (rec[n] for n = 0 ... size - 1, one more than
/// coefficients): S(x) = sum s_n Q_n(x) and S'(x) = dS/dx. [F10a], Eq. (2.2) (alpha_n =
/// s_n + (a_n + b_n x) alpha_(n+1) - c_(n+1) alpha_(n+2), S = alpha_0 by Eq. (A.10), since
/// Q_0 = 1 and Q_1 = (a_0 + b_0 x) Q_0) and Eq. (2.3) with j = 1 (alpha'_n = b_n alpha_(n+1) +
/// (a_n + b_n x) alpha'_(n+1) - c_(n+1) alpha'_(n+2), S' = alpha'_0). Working down from n = M
/// with alpha_(M+1) = alpha_(M+2) = 0 gives the start values of the paper.
template <rtt::math::Real T>
[[nodiscard]] std::pair<T, T> qcon_clenshaw(std::span<const T> s,
                                            std::span<const QConRecurrence<T>> rec,
                                            T x) noexcept {
  T a1 = T(0);  // alpha_(n+1)
  T a2 = T(0);  // alpha_(n+2)
  T d1 = T(0);
  T d2 = T(0);
  for (std::size_t k = s.size(); k-- > 0;) {
    const T p = rec[k].a + rec[k].b * x;
    const T c1 = rec[k + 1].c;
    const T a0 = s[k] + p * a1 - c1 * a2;
    const T d0 = rec[k].b * a1 + p * d1 - c1 * d2;
    a2 = a1;
    a1 = a0;
    d2 = d1;
    d1 = d0;
  }
  return {a1, d1};
}

}  // namespace detail

/// S(x) = sum_m s_m Q^con_m(x) and S'(x) = dS/dx (x = u^2) by Clenshaw's method ([F10a],
/// Eqs. (2.2), (2.3), (5.2)), without evaluating a single polynomial. Q^con_m(x) =
/// P_m^(0,4)(2x - 1) ([F07], Eq. (8)).
/// @param s the coefficients s_0, s_1, ... (any unit; S has the same)
/// @param x the argument u^2, dimensionless; any value (u > 1 extrapolates)
template <rtt::math::Real T>
[[nodiscard]] std::pair<T, T> qcon_sum(std::span<const T> s, T x) {
  std::vector<QConRecurrence<T>> rec;
  rec.reserve(s.size() + 1);
  for (std::size_t n = 0; n <= s.size(); ++n) rec.push_back(qcon_recurrence<T>(n));
  return detail::qcon_clenshaw<T>(s, rec, x);
}

/// The band of the change-of-basis matrix L between Q-bfs and the auxiliary basis P_m
/// ([F10b], Eqs. (2.5), (A.12)): P_m = f_m Q_m + g_(m-1) Q_(m-1) + h_(m-2) Q_(m-2), with f_m :=
/// L_mm, g_m := L_(m+1)m, h_m := L_(m+2)m. Start values f_0 = 2, f_1 = sqrt(19)/2,
/// g_0 = -1/2 and, from m = 2 on and in this order ([F10b], text before Eq. (A.14) and
/// Eqs. (A.14)-(A.16)): h_(m-2) = -m(m - 1)/(2 f_(m-2)), g_(m-1) = -(1 + g_(m-2) h_(m-2)) /
/// f_(m-1), f_m = sqrt(m(m + 1) + 3 - g_(m-1)^2 - h_(m-2)^2). f has m_max + 1 entries, g and h
/// as many (the last one or two are 0 and unused).
template <rtt::math::Real T>
struct QBfsMatrix {
  std::vector<T> f;
  std::vector<T> g;
  std::vector<T> h;
};

/// The matrix band of QBfsMatrix up to m_max.
template <rtt::math::Real T>
[[nodiscard]] QBfsMatrix<T> qbfs_matrix(std::size_t m_max) {
  using std::sqrt;
  // Built with push_back: GCC at -O2 reports a false -Wnull-dereference for element stores
  // into the vector members of a freshly constructed aggregate.
  std::vector<T> f;
  std::vector<T> g;
  std::vector<T> h;
  f.reserve(m_max + 1);
  g.reserve(m_max + 1);
  h.reserve(m_max + 1);
  f.push_back(T(2));
  g.push_back(m_max >= 1 ? T(-1) / T(2) : T(0));
  if (m_max >= 1) f.push_back(sqrt(T(19)) / T(2));
  for (std::size_t m = 2; m <= m_max; ++m) {
    const T mm = static_cast<T>(m);                         // products in T: no integer overflow
    h.push_back(-(mm * (mm - T(1))) / (T(2) * f[m - 2]));   // h_(m-2)
    g.push_back(-(T(1) + g[m - 2] * h[m - 2]) / f[m - 1]);  // g_(m-1)
    f.push_back(sqrt(mm * (mm + T(1)) + T(3) - g[m - 1] * g[m - 1] - h[m - 2] * h[m - 2]));  // f_m
  }
  g.resize(m_max + 1, T(0));
  h.resize(m_max + 1, T(0));
  return {std::move(f), std::move(g), std::move(h)};
}

/// The coefficients b_m of S = sum b_m P_m in the auxiliary basis from a_m of S = sum a_m Q_m
/// ([F10b], Eqs. (3.4a,b), (3.5)): b_M = a_M / f_M, b_(M-1) = (a_(M-1) - g_(M-1) b_M) / f_(M-1),
/// b_m = (a_m - g_m b_(m+1) - h_m b_(m+2)) / f_m for m = M - 2 down to 0.
template <rtt::math::Real T>
[[nodiscard]] std::vector<T> qbfs_auxiliary(std::span<const T> a) {
  if (a.empty()) return {};
  const std::size_t m_max = a.size() - 1;
  const QBfsMatrix<T> l = qbfs_matrix<T>(m_max);
  std::vector<T> b(a.size(), T(0));
  b[m_max] = a[m_max] / l.f[m_max];
  if (m_max >= 1) b[m_max - 1] = (a[m_max - 1] - l.g[m_max - 1] * b[m_max]) / l.f[m_max - 1];
  if (m_max >= 2) {
    for (std::size_t m = m_max - 1; m-- > 0;) {
      b[m] = (a[m] - l.g[m] * b[m + 1] - l.h[m] * b[m + 2]) / l.f[m];
    }
  }
  return b;
}

namespace detail {

/// S(x) = sum b_m P_m(x) and S'(x) in the auxiliary basis P_(m+1) = (2 - 4x) P_m - P_(m-1),
/// P_0 = 2, P_1 = 6 - 8x ([F10b], Eq. (2.6)), by Clenshaw's method: alpha_M = b_M,
/// alpha_(M-1) = b_(M-1) + (2 - 4x) alpha_M, alpha_m = b_m + (2 - 4x) alpha_(m+1) - alpha_(m+2),
/// S = 2 (alpha_0 + alpha_1) ([F10b], Eqs. (3.7a,b), (3.8), (3.9)); the derivative with
/// alpha'_M = 0, alpha'_(M-1) = -4 alpha_M, alpha'_m = (2 - 4x) alpha'_(m+1) - alpha'_(m+2) -
/// 4 alpha_(m+1), S' = 2 (alpha'_0 + alpha'_1) ([F10b], Eqs. (3.10a,b), (3.11), (3.12), j = 1).
/// Working down from m = M with zeros beyond M gives these start values.
template <rtt::math::Real T>
[[nodiscard]] std::pair<T, T> qbfs_clenshaw(std::span<const T> b, T x) noexcept {
  const T p = T(2) - T(4) * x;
  T a1 = T(0);
  T a2 = T(0);
  T d1 = T(0);
  T d2 = T(0);
  for (std::size_t m = b.size(); m-- > 0;) {
    const T a0 = b[m] + p * a1 - a2;
    const T d0 = p * d1 - d2 - T(4) * a1;
    a2 = a1;
    a1 = a0;
    d2 = d1;
    d1 = d0;
  }
  return {T(2) * (a1 + a2), T(2) * (d1 + d2)};
}

}  // namespace detail

/// S(x) = sum_m a_m Q^bfs_m(x) and S'(x) = dS/dx (x = u^2) through the auxiliary basis
/// ([F10b], Section 3): a -> b by qbfs_auxiliary(), then Clenshaw. Q^bfs_m as [F07], Eq. (14)
/// (Q_0 = 1, Q_1 = (13 - 16x)/sqrt(19), ...).
/// @param a the coefficients a_0, a_1, ... (any unit; S has the same)
/// @param x the argument u^2, dimensionless; any value (u > 1 extrapolates)
template <rtt::math::Real T>
[[nodiscard]] std::pair<T, T> qbfs_sum(std::span<const T> a, T x) {
  const std::vector<T> b = qbfs_auxiliary<T>(a);
  return detail::qbfs_clenshaw<T>(b, x);
}

/// Q-con asphere ([F07], Eqs. (2), (3), (8); [F10a], Eq. (5.1)):
///   z(r) = c r^2 / (1 + phi) + u^4 sum_m a_m Q^con_m(u^2), phi = sqrt(1 - (1 + k) c^2 r^2),
/// u = r / rho. [F07] writes the conic with epsilon = 1 + k. c in 1/mm (0: plane base), the a_m
/// in mm. Gradient (dz/dx, dz/dy) = (x, y) g with g = c / phi + u^2 (4 S + 2 u^2 S') / rho^2,
/// [F10a], Eq. (5.4) divided by r (S and S' at x = u^2). The polynomial part is added only if a
/// coefficient is not zero (ADR 0043, point 5.1), so that without one the shape is the conic bit
/// for bit. Domain, max_radius() and start value are those of the base conic (c, k).
template <rtt::math::Real T>
class QCon final {
 public:
  /// @param curvature            c of the base conic in 1/mm
  /// @param conic_constant       k of the base conic, dimensionless
  /// @param normalization_radius rho in mm, > 0
  /// @param coefficients         a_0, a_1, ... in mm
  QCon(T curvature, T conic_constant, T normalization_radius, std::vector<T> coefficients)
      : conic_(curvature, conic_constant),
        rho_(normalization_radius),
        coefficients_(std::move(coefficients)),
        has_polynomial_(std::any_of(
            coefficients_.begin(), coefficients_.end(), [](T a) { return a != T(0); })) {
    recurrence_.reserve(coefficients_.size() + 1);
    for (std::size_t n = 0; n <= coefficients_.size(); ++n) {
      recurrence_.push_back(qcon_recurrence<T>(n));
    }
  }

  /// The coefficients a_0, a_1, ... in mm.
  [[nodiscard]] const std::vector<T>& coefficients() const noexcept { return coefficients_; }

  /// The normalization radius rho in mm.
  [[nodiscard]] T normalization_radius() const noexcept { return rho_; }

  /// Sag z in mm at (x, y) in mm; NaN outside the domain of the base conic.
  [[nodiscard]] T sag(T x, T y) const noexcept {
    const T base = conic_.sag(x, y);
    if (!has_polynomial_) return base;
    const T u2 = (x * x + y * y) / (rho_ * rho_);
    const T s = detail::qcon_clenshaw<T>(coefficients_, recurrence_, u2).first;
    return base + u2 * u2 * s;
  }

  /// Slopes (dz/dx, dz/dy); NaN outside the domain, +-infinity on its rim.
  [[nodiscard]] std::pair<T, T> grad(T x, T y) const noexcept {
    const auto [cx, cy] = conic_.grad(x, y);
    if (!has_polynomial_) return {cx, cy};
    const T rho2 = rho_ * rho_;
    const T u2 = (x * x + y * y) / rho2;
    const auto [s, ds] = detail::qcon_clenshaw<T>(coefficients_, recurrence_, u2);
    const T g = u2 * (T(4) * s + T(2) * u2 * ds) / rho2;
    return {cx + x * g, cy + y * g};
  }

  /// (c, k) of the base conic.
  [[nodiscard]] std::pair<T, T> base_conic() const { return conic_.base_conic(); }

  /// Domain of the base conic; the polynomial is defined everywhere.
  [[nodiscard]] std::optional<T> max_radius() const { return conic_.max_radius(); }

  /// Vertex expansion (ADR 0043, point 3): that of the base conic (curvature c, order 4,
  /// rotational) with the native K = k and A4 = sum a_m (-1)^m C(m + 4, 4) / rho^4, the r^4
  /// coefficient of the polynomial part: u^4 S(0) with Q^con_m(0) = P_m^(0,4)(-1) = (-1)^m
  /// C(m + 4, 4). Without a non-zero coefficient it is the conic's bit for bit.
  [[nodiscard]] VertexExpansion<T> vertex_expansion() const noexcept {
    VertexExpansion<T> e = conic_.vertex_expansion();
    if (!has_polynomial_) return e;
    T s0 = T(0);
    for (std::size_t m = 0; m < coefficients_.size(); ++m) {
      // C(m + 4, 4) in T: exact while the product of the four factors stays below 2^53 (m up to
      // about 9000 for double; the product is divisible by 4! = 24), no integer overflow beyond.
      const T mm = static_cast<T>(m);
      const T binomial = (mm + T(1)) * (mm + T(2)) * (mm + T(3)) * (mm + T(4)) / T(24);
      const T q0 = m % 2 == 0 ? binomial : -binomial;
      s0 += q0 * coefficients_[m];
    }
    const T rho2 = rho_ * rho_;
    e.a4 = s0 / (rho2 * rho2);
    return e;
  }

 private:
  Conic<T> conic_;
  T rho_;
  std::vector<T> coefficients_;
  bool has_polynomial_;
  std::vector<QConRecurrence<T>> recurrence_;  ///< [F10a], Eq. (5.2) for n = 0 ... M + 1
};

/// Q-bfs asphere with a conic base ([F11], Eq. (2.4); for k = 0 [F07], Eqs. (10), (11) and
/// [F10b], Eq. (2.1)):
///   z(r) = c r^2 / (1 + phi) + u^2 (1 - u^2) / sigma * sum_m a_m Q^bfs_m(u^2),
/// phi = sqrt(1 - (1 + k) c^2 r^2) ([F11], Eq. (2.2)), 1/sigma = sqrt(1 - k c^2 r^2) / phi
/// ([F11], Eq. (2.3): sigma is the cosine between the normal of the base and the axis), u =
/// r / rho. c in 1/mm (0: plane base), k dimensionless, the a_m in mm. The polynomials do not
/// depend on k ([F11], p. 9927).
///
/// Gradient (dz/dx, dz/dy) = (x, y) g with
///   g = c / phi + (S (2 - 4u^2) + 2 u^2 (1 - u^2) S') / (sigma rho^2)
///       + u^2 (1 - u^2) S c^2 sigma / phi^4,
/// [F11], Eq. (2.11) divided by r; for k = 0 equal to [F10b], Eq. (3.14). The domain is
/// phi > 0, strictly, because 1/sigma diverges at phi = 0 (ADR 0043, point 1): in_domain(x, y) is
/// (1 + k) c^2 r^2 < 1; outside it sag and grad are NaN. max_radius() and the start value are
/// those of the conic (c, k). The polynomial part is added only if a coefficient is not zero.
template <rtt::math::Real T>
class QBfs final {
 public:
  /// @param curvature            c of the base conic in 1/mm
  /// @param conic_constant       k of the base conic, dimensionless
  /// @param normalization_radius rho in mm, > 0
  /// @param coefficients         a_0, a_1, ... in mm
  QBfs(T curvature, T conic_constant, T normalization_radius, std::vector<T> coefficients)
      : conic_(curvature, conic_constant),
        rho_(normalization_radius),
        coefficients_(std::move(coefficients)),
        auxiliary_(qbfs_auxiliary<T>(coefficients_)),
        has_polynomial_(std::any_of(
            coefficients_.begin(), coefficients_.end(), [](T a) { return a != T(0); })) {}

  /// The coefficients a_0, a_1, ... in mm.
  [[nodiscard]] const std::vector<T>& coefficients() const noexcept { return coefficients_; }

  /// The normalization radius rho in mm.
  [[nodiscard]] T normalization_radius() const noexcept { return rho_; }

  /// The domain phi > 0: (1 + k) c^2 r^2 < 1.
  [[nodiscard]] bool in_domain(T x, T y) const noexcept {
    const auto [c, k] = conic_.base_conic();
    return (T(1) + k) * c * c * (x * x + y * y) < T(1);
  }

  /// Sag z in mm at (x, y) in mm; NaN outside in_domain().
  [[nodiscard]] T sag(T x, T y) const noexcept {
    if (!in_domain(x, y)) return std::numeric_limits<T>::quiet_NaN();
    const T base = conic_.sag(x, y);
    if (!has_polynomial_) return base;
    const T r2 = x * x + y * y;
    const T u2 = r2 / (rho_ * rho_);
    const T s = detail::qbfs_clenshaw<T>(auxiliary_, u2).first;
    return base + u2 * (T(1) - u2) * s * inverse_sigma(r2);
  }

  /// Slopes (dz/dx, dz/dy); NaN outside in_domain().
  [[nodiscard]] std::pair<T, T> grad(T x, T y) const noexcept {
    using std::sqrt;
    if (!in_domain(x, y)) {
      const T nan = std::numeric_limits<T>::quiet_NaN();
      return {nan, nan};
    }
    const auto [cx, cy] = conic_.grad(x, y);
    if (!has_polynomial_) return {cx, cy};
    const auto [c, k] = conic_.base_conic();
    const T r2 = x * x + y * y;
    const T rho2 = rho_ * rho_;
    const T u2 = r2 / rho2;
    const auto [s, ds] = detail::qbfs_clenshaw<T>(auxiliary_, u2);
    const T phi = sqrt(T(1) - (T(1) + k) * c * c * r2);
    const T inv_sigma = sqrt(T(1) - k * c * c * r2) / phi;
    const T phi2 = phi * phi;
    const T g = (s * (T(2) - T(4) * u2) + T(2) * u2 * (T(1) - u2) * ds) * inv_sigma / rho2 +
                u2 * (T(1) - u2) * s * c * c / (inv_sigma * phi2 * phi2);
    return {cx + x * g, cy + y * g};
  }

  /// (c, k) of the base conic.
  [[nodiscard]] std::pair<T, T> base_conic() const { return conic_.base_conic(); }

  /// The outer extent of the domain: that of the conic (in_domain() excludes its rim).
  [[nodiscard]] std::optional<T> max_radius() const { return conic_.max_radius(); }

  /// Vertex expansion (ADR 0043, point 3): value 0, slope 0, order 4, rotational; the vertex
  /// curvature c_v = c + 2 S_0 / rho^2 ([F10b], Eq. (3.16); for k != 0 the same, since
  /// 1/sigma(0) = 1) and the cap with the native K = k and A4 = (1 + k)(c^3 - c_v^3) / 8 +
  /// (S_1 - S_0) / rho^4 + c^2 S_0 / (2 rho^2), with S_0 = S(0) and S_1 = S'(0): the r^4
  /// coefficient of the polynomial part is (S_1 - S_0) / rho^4 + c^2 S_0 / (2 rho^2) (1/sigma =
  /// 1 + c^2 r^2 / 2 + ...), and K c_v^3 + 8 A4 = 8 b4 - c_v^3 (ADR 0041, point 1). Without a
  /// non-zero coefficient it is the conic's bit for bit.
  [[nodiscard]] VertexExpansion<T> vertex_expansion() const noexcept {
    VertexExpansion<T> e = conic_.vertex_expansion();
    if (!has_polynomial_) return e;
    const auto [c, k] = conic_.base_conic();
    const auto [s0, s1] = detail::qbfs_clenshaw<T>(auxiliary_, T(0));
    const T rho2 = rho_ * rho_;
    const T cv = c + T(2) * s0 / rho2;
    e.hessian = {cv, T(0), cv};
    e.conic = k;
    e.a4 = (T(1) + k) * (c * c * c - cv * cv * cv) / T(8) + (s1 - s0) / (rho2 * rho2) +
           c * c * s0 / (T(2) * rho2);
    return e;
  }

 private:
  /// 1/sigma = sqrt(1 - k c^2 r^2) / phi ([F11], Eq. (2.3)) at r^2 inside the domain.
  [[nodiscard]] T inverse_sigma(T r2) const noexcept {
    using std::sqrt;
    const auto [c, k] = conic_.base_conic();
    return sqrt(T(1) - k * c * c * r2) / sqrt(T(1) - (T(1) + k) * c * c * r2);
  }

  Conic<T> conic_;
  T rho_;
  std::vector<T> coefficients_;
  std::vector<T> auxiliary_;  ///< b_m of [F10b], Eqs. (3.4), (3.5)
  bool has_polynomial_;
};

}  // namespace rtt::geom
