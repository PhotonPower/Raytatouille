#include <algorithm>
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "rtt/geom/conic.hpp"
#include "rtt/geom/forbes.hpp"
#include "rtt/geom/intersect.hpp"
#include "rtt/geom/shape.hpp"

// Forbes Q-con and Q-bfs of ADR 0043, point 1, and the acceptance 1, 2, 4, 5 and 11 of ADR 0043
// for the shapes themselves. Sources: [F07] Opt. Express 15(8), 5218 (2007); [F10a] Opt. Express
// 18(13), 13851 (2010); [F10b] Opt. Express 18(19), 19700 (2010); [F11] Opt. Express 19(10),
// 9923 (2011) (docs/quellen.md).

#include "forbes_tables.inc"  // kQConValue, kQConDerivative, kQBfsValue, kQBfsDerivative

using rtt::geom::Conic;
using rtt::geom::DomainShape;
using rtt::geom::HitStatus;
using rtt::geom::IntersectableShape;
using rtt::geom::Intersection;
using rtt::geom::ParaxialShape;
using rtt::geom::QBfs;
using rtt::geom::QCon;
using rtt::geom::VertexExpansion;
using Vec3 = rtt::math::Vec3T<double>;

static_assert(IntersectableShape<QCon<double>, double>);
static_assert(ParaxialShape<QCon<double>, double>);
static_assert(!DomainShape<QCon<double>, double>);
static_assert(IntersectableShape<QBfs<double>, double>);
static_assert(ParaxialShape<QBfs<double>, double>);
// Q-bfs is defined for phi > 0 only (ADR 0043, point 1).
static_assert(DomainShape<QBfs<double>, double>);

namespace {

constexpr double kEps = std::numeric_limits<double>::epsilon();

void require_same(const VertexExpansion<double>& a, const VertexExpansion<double>& b) {
  CHECK(a.value == b.value);
  CHECK(a.slope == b.slope);
  CHECK(a.hessian == b.hessian);
  CHECK(a.order == b.order);
  CHECK(a.rotational4 == b.rotational4);
  CHECK(a.conic == b.conic);
  CHECK(a.a4 == b.a4);
}

/// max over the table of |row|.
double max_abs(std::span<const double> row) {
  double m = 0.0;
  for (const double v : row) m = std::max(m, std::abs(v));
  return m;
}

/// The unit coefficient vector e_m.
std::vector<double> unit(std::size_t m) {
  std::vector<double> e(m + 1, 0.0);
  e[m] = 1.0;
  return e;
}

/// The exact Cartesian oval of [F07], Section 4: an object point 10 mm in front of the vertex in
/// glass n = 1.5 is imaged without error to a point 10 mm behind it in air, n |P - O| + |P - I| =
/// n d1 + d2 with O = (0, -d1), I = (0, d2). The surface z(r) through the vertex: the first sign
/// change of the optical path difference on a scan from z = 0 to -d1, refined by bisection.
double oval(double r) {
  const double n = 1.5;
  const double d1 = 10.0;
  const double d2 = 10.0;
  if (r == 0.0) return 0.0;
  const auto f = [&](double z) {
    return n * std::hypot(r, z + d1) + std::hypot(r, z - d2) - (n * d1 + d2);
  };
  double hi = 0.0;
  double f_hi = f(hi);
  for (int i = 1; i <= 2000; ++i) {
    const double lo = -d1 * i / 2000.0;
    const double f_lo = f(lo);
    if ((f_lo < 0.0) != (f_hi < 0.0)) {
      double a = lo;
      double b = hi;
      const double fa = f_lo;
      for (int k = 0; k < 200; ++k) {
        const double mid = a + (b - a) / 2.0;
        if (mid == a || mid == b) break;
        if ((f(mid) < 0.0) == (fa < 0.0)) {
          a = mid;
        } else {
          b = mid;
        }
      }
      return a + (b - a) / 2.0;
    }
    hi = lo;
    f_hi = f_lo;
  }
  // A throw, not FAIL plus return: MSVC flags the return as unreachable (C4702).
  throw std::logic_error("no root of the oval at r = " + std::to_string(r));
}

}  // namespace

TEST_CASE("Q-con polynomials and their derivatives against exact values", "[geom][forbes]") {
  // ADR 0043, acceptance 2: Clenshaw with [F10a], Eqs. (2.2), (2.3), (5.2) against the exact
  // values of forbes_tables.py (Gram-Schmidt with the weight x^4, [F07]) at x = k/64, rounded
  // once. Empirical bound of the form n eps max_x |Q_m| with n = 4 (m + 1), for the derivative
  // with max(1, max_x |Q'_m|); measured 0.29 (value) and 0.23 (derivative) of it.
  double worst = 0.0;
  for (std::size_t m = 0; m <= 20; ++m) {
    INFO("m = " << m);
    const std::vector<double> e = unit(m);
    const double bound_v = 4.0 * static_cast<double>(m + 1) * kEps * max_abs(kQConValue[m]);
    const double bound_d =
        4.0 * static_cast<double>(m + 1) * kEps * std::max(1.0, max_abs(kQConDerivative[m]));
    for (std::size_t k = 0; k <= 64; ++k) {
      const auto [s, ds] = rtt::geom::qcon_sum<double>(e, static_cast<double>(k) / 64.0);
      CHECK(std::abs(s - kQConValue[m][k]) <= bound_v);
      CHECK(std::abs(ds - kQConDerivative[m][k]) <= bound_d);
      worst = std::max(worst, std::abs(s - kQConValue[m][k]) / bound_v);
    }
  }
  // Content guard: [F07], Eq. (9), Q_1 = 6x - 5 and Q_2 = 15 - 14x (3 - 2x) at x = 0.3.
  CHECK(std::abs(rtt::geom::qcon_sum<double>(unit(1), 0.3).first - (6.0 * 0.3 - 5.0)) <= 1e-15);
  CHECK(std::abs(rtt::geom::qcon_sum<double>(unit(2), 0.3).first -
                 (15.0 - 14.0 * 0.3 * (3.0 - 2.0 * 0.3))) <= 1e-14);
  CHECK(worst > 0.0);  // the comparison is not trivially exact
}

TEST_CASE("Q-bfs polynomials and their derivatives against exact values", "[geom][forbes]") {
  // ADR 0043, acceptance 2: [F10b], Eqs. (3.4), (3.5), (3.7)-(3.12) against the exact values of
  // forbes_tables.py (Gram-Schmidt in the slope, [F10b], Eqs. (2.2)-(2.4)), m <= 16. Bound
  // empirical, n eps max_x |Q_m| with n = 8 (m + 1), for the derivative with max(1, max_x |Q'_m|);
  // measured 0.67 (value) and 0.29 (derivative) of it.
  for (std::size_t m = 0; m <= 16; ++m) {
    INFO("m = " << m);
    const std::vector<double> e = unit(m);
    const double bound_v = 8.0 * static_cast<double>(m + 1) * kEps * max_abs(kQBfsValue[m]);
    const double bound_d =
        8.0 * static_cast<double>(m + 1) * kEps * std::max(1.0, max_abs(kQBfsDerivative[m]));
    for (std::size_t k = 0; k <= 64; ++k) {
      const auto [s, ds] = rtt::geom::qbfs_sum<double>(e, static_cast<double>(k) / 64.0);
      CHECK(std::abs(s - kQBfsValue[m][k]) <= bound_v);
      CHECK(std::abs(ds - kQBfsDerivative[m][k]) <= bound_d);
    }
  }
  // Content guard: [F07], Eq. (14), Q_1 = (13 - 16x)/sqrt(19) at x = 0.3.
  CHECK(std::abs(rtt::geom::qbfs_sum<double>(unit(1), 0.3).first -
                 (13.0 - 16.0 * 0.3) / std::sqrt(19.0)) <= 1e-15);
}

TEST_CASE("Q-bfs change of basis: the printed block of [F10b], Eq. (A.17)", "[geom][forbes]") {
  // Content guard of Eqs. (A.14)-(A.16): f_0 ... f_5, g_0 ... g_4 and h_0 ... h_3.
  const auto l = rtt::geom::qbfs_matrix<double>(5);
  const auto s = [](double v) { return std::sqrt(v); };
  const std::array<double, 6> f = {2.0,
                                   s(19.0 / 4.0),
                                   4.0 * s(10.0 / 19.0),
                                   0.5 * s(509.0 / 10.0),
                                   6.0 * s(259.0 / 509.0),
                                   0.5 * s(25607.0 / 259.0)};
  const std::array<double, 5> g = {-0.5, -5.0 / (2.0 * s(19.0)), -17.0 / (2.0 * s(190.0)),
                                   -91.0 / (2.0 * s(5090.0)), -473.0 / (2.0 * s(131831.0))};
  const std::array<double, 4> h = {-0.5, -6.0 / s(19.0), -1.5 * s(19.0 / 10.0),
                                   -20.0 * s(10.0 / 509.0)};
  for (std::size_t i = 0; i < f.size(); ++i) CHECK(std::abs(l.f[i] - f[i]) <= 1e-14);
  for (std::size_t i = 0; i < g.size(); ++i) CHECK(std::abs(l.g[i] - g[i]) <= 1e-14);
  for (std::size_t i = 0; i < h.size(); ++i) CHECK(std::abs(l.h[i] - h[i]) <= 1e-14);
}

TEST_CASE("Q-con and Q-bfs without non-zero coefficients are the conic bit for bit",
          "[geom][forbes]") {
  // ADR 0043, acceptance 1: the polynomial part is only added with a non-zero coefficient.
  for (const auto& [c, k] :
       {std::pair{-1.0 / 32.0, -0.4}, std::pair{0.025, 0.0}, std::pair{0.03, -1.5}}) {
    const Conic<double> conic(c, k);
    for (const std::vector<double>& a :
         {std::vector<double>{}, std::vector<double>{0.0, -0.0, 0.0}}) {
      const QCon<double> qcon(c, k, 10.0, a);
      const QBfs<double> qbfs(c, k, 10.0, a);
      for (const auto& [x, y] : {std::pair{0.0, 0.0}, std::pair{3.0, -4.0}, std::pair{-9.5, 2.25},
                                 std::pair{1e-200, 0.0}}) {
        for (const auto& [zs, gs] : {std::pair{qcon.sag(x, y), qcon.grad(x, y)},
                                     std::pair{qbfs.sag(x, y), qbfs.grad(x, y)}}) {
          const double zc = conic.sag(x, y);
          const auto [hx, hy] = conic.grad(x, y);
          CHECK(zs == zc);
          CHECK(std::signbit(zs) == std::signbit(zc));
          CHECK(gs.first == hx);
          CHECK(gs.second == hy);
        }
      }
      require_same(qcon.vertex_expansion(), conic.vertex_expansion());
      require_same(qbfs.vertex_expansion(), conic.vertex_expansion());
      CHECK(qcon.base_conic() == conic.base_conic());
      CHECK(qbfs.max_radius() == conic.max_radius());
    }
  }
}

TEST_CASE("Q-con and Q-bfs with zero coefficients intersect as the analytic conic",
          "[geom][forbes][newton]") {
  // ADR 0043, acceptance 1: Newton against the closed form within kNewtonTolerance / |F'| plus
  // the rounding 8 eps (|o| + t).
  const double c = 1.0 / 40.0;
  const double k = -0.6;
  const Conic<double> conic(c, k);
  const QCon<double> qcon(c, k, 10.0, {0.0, 0.0});
  const QBfs<double> qbfs(c, k, 10.0, {0.0});
  int compared = 0;
  for (int i = -4; i <= 4; ++i) {
    for (int j = -4; j <= 4; ++j) {
      const Vec3 o(2.0 * i, 2.0 * j, -5.0);
      const Vec3 d = Vec3(0.01 * j, -0.02 * i, 1.0).normalized();
      const Intersection<double> ref = rtt::geom::intersect(conic, o, d);
      REQUIRE(ref.status == HitStatus::Hit);
      const auto [gx, gy] = conic.grad(ref.point.x(), ref.point.y());
      const double df = d.z() - gx * d.x() - gy * d.y();
      const double bound =
          rtt::geom::kNewtonTolerance / std::abs(df) + 8.0 * kEps * (o.norm() + ref.t);
      for (const Intersection<double>& hit :
           {rtt::geom::intersect(qcon, o, d), rtt::geom::intersect(qbfs, o, d)}) {
        REQUIRE(hit.status == HitStatus::Hit);
        CHECK(std::abs(hit.t - ref.t) <= bound);
      }
      ++compared;
    }
  }
  CHECK(compared == 81);
}

TEST_CASE("Q-con vertex expansion: A4 from S(0) with Q_m(0) = (-1)^m C(m + 4, 4)",
          "[geom][forbes][vertex]") {
  // ADR 0043, acceptance 4: the cap against S(0) from the table (x = 0 is k = 0); c_v = c.
  const double c = 1.0 / 25.0;
  const double k = -0.7;
  const double rho = 10.0;
  const std::vector<double> a = {0.012, -0.004, 0.0021, 7e-4, -3e-4};
  const QCon<double> qcon(c, k, rho, a);
  const VertexExpansion<double> e = qcon.vertex_expansion();
  double s0 = 0.0;
  double scale = 0.0;
  for (std::size_t m = 0; m < a.size(); ++m) {
    s0 += a[m] * kQConValue[m][0];
    scale += std::abs(a[m] * kQConValue[m][0]);
  }
  const double rho4 = rho * rho * rho * rho;
  CHECK(e.hessian == Conic<double>(c, k).vertex_expansion().hessian);
  CHECK(e.conic == k);
  CHECK(e.order == 4);
  CHECK(e.rotational4);
  CHECK(std::abs(e.a4 - s0 / rho4) <= (static_cast<double>(a.size()) + 2.0) * kEps * scale / rho4);
  CHECK(e.a4 != 0.0);
}

TEST_CASE("Q-bfs vertex expansion: c + 2 S_0 / rho^2 and the cap with native K",
          "[geom][forbes][vertex]") {
  // ADR 0043, point 3 and acceptance 4: S_0 and S_1 from the tables ([F10b], Eq. (3.16) for
  // c_v); then a fit of the sag r^2, r^4, r^6 at r = 1/8, 1/4, 1/2 as independent check of the
  // cap K c_v^3 + 8 A4 = 8 b4 - c_v^3 (ADR 0041, point 1). The neglected r^8 terms bias b4 by
  // about 1e-7 relative here (1.5e-6 with r up to 1; a 60-digit fit agrees with the closed
  // form to 1e-37), the rounding of the sag contributes 1e-10.
  const double c = 1.0 / 30.0;
  const double k = -0.6;
  const double rho = 12.0;
  const std::vector<double> a = {3e-3, -1.2e-3, 4e-4, -1e-4};
  const QBfs<double> qbfs(c, k, rho, a);
  const VertexExpansion<double> e = qbfs.vertex_expansion();
  double s0 = 0.0;
  double s1 = 0.0;
  double scale = 0.0;
  for (std::size_t m = 0; m < a.size(); ++m) {
    s0 += a[m] * kQBfsValue[m][0];
    s1 += a[m] * kQBfsDerivative[m][0];
    scale += std::abs(a[m]) * (max_abs(kQBfsValue[m]) + max_abs(kQBfsDerivative[m]));
  }
  // |S - S_table| <= 8 (M + 1) eps sum |a_m| max |Q_m| (acceptance 2) for S and S'.
  const double ds = 8.0 * static_cast<double>(a.size()) * kEps * scale;
  const double rho2 = rho * rho;
  const double cv = c + 2.0 * s0 / rho2;
  CHECK(std::abs(e.hessian[0] - cv) <= 2.0 * ds / rho2 + 4.0 * kEps * std::abs(cv));
  CHECK(e.hessian[0] == e.hessian[2]);
  CHECK(e.hessian[1] == 0.0);
  CHECK(e.conic == k);
  CHECK(e.order == 4);
  const double a4 = (1.0 + k) * (c * c * c - cv * cv * cv) / 8.0 + (s1 - s0) / (rho2 * rho2) +
                    c * c * s0 / (2.0 * rho2);
  CHECK(std::abs(e.a4 - a4) <= 1e-12 * std::abs(a4));
  // Independent: the r^4 coefficient of the sag by a three-term fit.
  const std::array<double, 3> h = {0.125, 0.25, 0.5};
  std::array<std::array<double, 3>, 3> m{};
  std::array<double, 3> z{};
  for (std::size_t i = 0; i < 3; ++i) {
    m[i] = {std::pow(h[i], 2.0), std::pow(h[i], 4.0), std::pow(h[i], 6.0)};
    z[i] = qbfs.sag(h[i], 0.0);
  }
  const auto det3 = [](const std::array<std::array<double, 3>, 3>& q) {
    return q[0][0] * (q[1][1] * q[2][2] - q[1][2] * q[2][1]) -
           q[0][1] * (q[1][0] * q[2][2] - q[1][2] * q[2][0]) +
           q[0][2] * (q[1][0] * q[2][1] - q[1][1] * q[2][0]);
  };
  auto ma = m;
  auto mb = m;
  for (std::size_t i = 0; i < 3; ++i) {
    ma[i][0] = z[i];
    mb[i][1] = z[i];
  }
  const double fit_c = 2.0 * det3(ma) / det3(m);
  const double b4 = det3(mb) / det3(m);
  CHECK(std::abs(fit_c - e.hessian[0]) <= 1e-9 * std::abs(e.hessian[0]));
  const double cap = e.conic * cv * cv * cv + 8.0 * e.a4;
  CHECK(std::abs(cap - (8.0 * b4 - cv * cv * cv)) <= 1e-6 * std::abs(cap));
  CHECK(std::abs(cv - c) > 1e-6);  // content guard: S_0 moves the vertex curvature
}

TEST_CASE("Q-con and Q-bfs gradients against central differences, also beyond rho",
          "[geom][forbes]") {
  // ADR 0043, acceptance 5: up to u = 1.3 (extrapolation, F3). Step 1e-6 mm.
  const std::vector<double> a = {0.012, -0.004, 0.0021, 7e-4, -3e-4};
  const QCon<double> qcon(1.0 / 25.0, -0.7, 10.0, a);
  const QBfs<double> qbfs(1.0 / 40.0, -0.6, 10.0, a);
  const QBfs<double> plane_base(0.0, 0.0, 10.0, a);
  int checked = 0;
  for (const auto& [x, y] : {std::pair{0.5, 0.0}, std::pair{3.0, -4.0}, std::pair{-6.0, 6.5},
                             std::pair{0.0, -10.0}, std::pair{9.0, 8.0}, std::pair{-12.0, 5.0}}) {
    const double hh = 1e-6;
    const auto check = [&](const auto& shape) {
      const auto [gx, gy] = shape.grad(x, y);
      const double nx = (shape.sag(x + hh, y) - shape.sag(x - hh, y)) / (2.0 * hh);
      const double ny = (shape.sag(x, y + hh) - shape.sag(x, y - hh)) / (2.0 * hh);
      CHECK(std::abs(gx - nx) <= 1e-8);
      CHECK(std::abs(gy - ny) <= 1e-8);
    };
    INFO("(" << x << ", " << y << ")");
    check(qcon);
    check(qbfs);
    check(plane_base);
    ++checked;
  }
  CHECK(checked == 6);
}

TEST_CASE("Q-bfs domain: phi > 0 strictly, finite up to the rim, NaN beyond", "[geom][forbes]") {
  // ADR 0043, point 1 and acceptance 11.
  const double c = 1.0 / 20.0;
  const double k = 0.5;  // rim at r = 1 / (c sqrt(1 + k))
  const QBfs<double> qbfs(c, k, 10.0, {1e-3, -2e-4});
  const double rim = 1.0 / (c * std::sqrt(1.0 + k));
  REQUIRE(qbfs.max_radius().has_value());
  CHECK(*qbfs.max_radius() == Conic<double>(c, k).max_radius());
  const double inside = rim * (1.0 - 1e-9);
  CHECK(qbfs.in_domain(inside, 0.0));
  CHECK(std::isfinite(qbfs.sag(inside, 0.0)));
  CHECK(std::isfinite(qbfs.grad(inside, 0.0).first));
  CHECK_FALSE(qbfs.in_domain(0.0, rim));
  CHECK(std::isnan(qbfs.sag(0.0, rim)));
  CHECK(std::isnan(qbfs.sag(rim * 1.1, 0.0)));
  CHECK(std::isnan(qbfs.grad(rim * 1.1, 0.0).second));
}

TEST_CASE("[F07] Table 1 as Q-con data against the exact Cartesian oval", "[geom][forbes]") {
  // ADR 0043, acceptance 2c: epsilon = 1 + k = 0.156284, c = -0.5 1/mm, rho = 3.5 mm, a_m in nm
  // ([F07], Table 1, column "a_m for Eq. (3)"). Bound 12 nm: fit error of [F07], Fig. 5(b)
  // (6.5 nm) + rounding of the eight printed coefficients (8 x 0.5 nm, |u^4 Q_m| <= 1) + the
  // rounding of epsilon (< 0.1 nm). Measured 5.9 nm.
  const std::vector<double> a = {70858e-6, -57065e-6, -9992e-6, -2825e-6,
                                 -705e-6,  -195e-6,   -54e-6,   -16e-6};
  const QCon<double> qcon(-0.5, 0.156284 - 1.0, 3.5, a);
  double worst = 0.0;
  for (int i = 0; i <= 350; ++i) {
    const double r = 3.5 * i / 350.0;
    worst = std::max(worst, std::abs(qcon.sag(0.0, r) - oval(r)));
  }
  CHECK(worst <= 12e-6);
  // Content guard: the polynomial part matters; the conic alone (it meets the oval at the rim
  // r = 3.5 mm) departs from it by up to tens of um inside ([F07], Section 4: 60 um).
  double conic_only = 0.0;
  for (int i = 0; i <= 35; ++i) {
    const double r = 0.1 * i;
    conic_only =
        std::max(conic_only, std::abs(Conic<double>(-0.5, 0.156284 - 1.0).sag(0.0, r) - oval(r)));
  }
  CHECK(conic_only > 30e-3);
}

TEST_CASE("[F10b] Eq. (4.8) as Q-bfs data against the parabola", "[geom][forbes]") {
  // ADR 0043, acceptance 2d: c = 1/25 1/mm (best-fit sphere), rho = 20 mm, a_m in nm. Bound 5 nm:
  // fit error about 2 nm ([F10b], Fig. 4) + rounding 7 x 0.5 nm x max|u^2 (1 - u^2) Q_m| / phi_min
  // (0.42 / 0.6). Measured 2.17 nm. The vertex curvature is 1/20 1/mm within 1e-5 relative: the
  // truncated b_7 = 2.6 nm contributes (4 / rho^2) 15 b_7 = 3.9e-7 1/mm.
  const std::vector<double> a = {2019004e-6, 7143e-6, -13944e-6, 4190e-6, -1095e-6, 283e-6, -68e-6};
  const QBfs<double> qbfs(1.0 / 25.0, 0.0, 20.0, a);
  double worst = 0.0;
  for (int i = 0; i <= 400; ++i) {
    const double r = 20.0 * i / 400.0;
    worst = std::max(worst, std::abs(qbfs.sag(r, 0.0) - r * r / 40.0));
  }
  CHECK(worst <= 5e-6);
  const double cv = qbfs.vertex_expansion().hessian[0];
  CHECK(std::abs(cv - 0.05) <= 1e-5 * 0.05);
  CHECK(std::abs(cv - 1.0 / 25.0) > 1e-3);  // content guard: not the best-fit sphere
}
