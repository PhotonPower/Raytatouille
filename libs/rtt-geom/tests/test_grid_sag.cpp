// The grid sag term in rtt-geom (#316, ADR 0051, points 5 and 9, PR 4): bicubic Hermite
// interpolation of a sphere against the derived bound, the order of convergence with finite
// differences, exact values at the grid points, continuity of the gradient, the domain and the
// intersection at the edge of the grid.

#include <algorithm>
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

#include "rtt/geom/intersect.hpp"
#include "rtt/geom/plane.hpp"
#include "rtt/geom/stacked.hpp"
#include "rtt/geom/terms.hpp"

using rtt::geom::GridSagTerm;
using rtt::geom::Plane;
using rtt::geom::Stacked;
using rtt::geom::TermVariant;

namespace {

constexpr double kU = 0x1.0p-53;  // unit roundoff

double unit(std::mt19937_64& g) {
  return static_cast<double>(g() >> 11U) * 0x1.0p-53;
}

/// The sphere z = R - sqrt(R^2 - x^2 - y^2) (vertex at the origin, R > 0) and its derivatives.
struct Sphere {
  double r = 100.0;
  [[nodiscard]] double s(double x, double y) const { return std::sqrt(r * r - x * x - y * y); }
  [[nodiscard]] double z(double x, double y) const { return r - s(x, y); }
  [[nodiscard]] double zx(double x, double y) const { return x / s(x, y); }
  [[nodiscard]] double zy(double x, double y) const { return y / s(x, y); }
  [[nodiscard]] double zxy(double x, double y) const {
    const double q = s(x, y);
    return x * y / (q * q * q);
  }
};

/// An n by n grid of `f` with spacing d, centred on (x0, y0); with or without derivatives.
GridSagTerm<double> sampled(
    const Sphere& f, int n, double d, bool derivatives, double x0 = 0.0, double y0 = 0.0) {
  std::vector<double> z;
  std::vector<double> zx;
  std::vector<double> zy;
  std::vector<double> zxy;
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < n; ++i) {
      const double x = x0 + (i - (n - 1) / 2.0) * d;
      const double y = y0 + (j - (n - 1) / 2.0) * d;
      z.push_back(f.z(x, y));
      if (derivatives) {
        zx.push_back(f.zx(x, y));
        zy.push_back(f.zy(x, y));
        zxy.push_back(f.zxy(x, y));
      }
    }
  }
  return {n, n, d, d, x0, y0, std::move(z), std::move(zx), std::move(zy), std::move(zxy)};
}

}  // namespace

TEST_CASE("grid sag of a sphere within the derived Hermite bound (ADR 0051, point 9)", "[grid]") {
  // ADR 0051, point 9: f - H_x H_y f = (f - H_x f) + H_x (f - H_y f), the one-dimensional
  // remainder f''''(xi) / 4! (x - a)^2 (x - b)^2 <= h^4 / 384 max |f''''| (Burden/Faires section
  // 3.4, docs/quellen.md) and |H_x g| <= max |g| + dx / 4 max |g_x| (|H_0| + |H_1| = 1,
  // |G_0| + |G_1| = t (1 - t) <= 1/4); with exact derivatives d/dx commutes with H_y. So
  //   eps = dx^4 / 384 M_x4 + dy^4 / 384 (M_y4 + dx / 4 M_xy4) + 16 u max |z|.
  // The Taylor series of the sphere at the vertex has only non-negative coefficients
  // (R - sqrt(R^2 - rho^2) = R sum c_k (rho^2 / R^2)^k, c_k > 0), so every derivative takes its
  // largest magnitude over the grid at the corner (X, Y) = (max |x|, max |y|). Closed forms
  // (sympy, s = sqrt(R^2 - x^2 - y^2)):
  //   z_xxxx = 3 (R^2 - y^2)(R^2 + 4 x^2 - y^2) / s^7,  z_yyyy = 3 (R^2 - x^2)(R^2 - x^2 + 4 y^2) /
  //   s^7, z_xyyyy = 3 x (3 R^4 - 6 R^2 x^2 + 24 R^2 y^2 + 3 x^4 - 24 x^2 y^2 + 8 y^4) / s^9.
  const Sphere f;
  constexpr int kN = 41;
  constexpr double kD = 0.5;
  const double big_x = (kN - 1) / 2.0 * kD;  // 10 mm
  const double r2 = f.r * f.r;
  const double x2 = big_x * big_x;
  const double s = f.s(big_x, big_x);
  const double s7 = std::pow(s, 7.0);
  const double m_x4 = 3.0 * (r2 - x2) * (r2 + 4.0 * x2 - x2) / s7;
  const double m_y4 = 3.0 * (r2 - x2) * (r2 - x2 + 4.0 * x2) / s7;
  const double m_xy4 = 3.0 * big_x *
                       (3.0 * r2 * r2 - 6.0 * r2 * x2 + 24.0 * r2 * x2 + 3.0 * x2 * x2 -
                        24.0 * x2 * x2 + 8.0 * x2 * x2) /
                       std::pow(s, 9.0);
  const double h4 = kD * kD * kD * kD / 384.0;
  const double eps = h4 * m_x4 + h4 * (m_y4 + kD / 4.0 * m_xy4) + 16.0 * kU * f.z(big_x, big_x);
  const GridSagTerm<double> grid = sampled(f, kN, kD, true);
  std::mt19937_64 g(316U);
  double worst = 0.0;
  for (int k = 0; k < 10000; ++k) {
    const double x = (2.0 * unit(g) - 1.0) * big_x;
    const double y = (2.0 * unit(g) - 1.0) * big_x;
    worst = std::max(worst, std::abs(grid.sag(x, y) - f.z(x, y)));
  }
  CHECK(worst <= eps);
  CHECK(worst > 0.0);
  // Content guard: slopes 1 % off break the bound by far.
  std::vector<double> zx(grid.dz_dx().begin(), grid.dz_dx().end());
  for (double& v : zx) v *= 1.01;
  const GridSagTerm<double> off(
      kN, kN, kD, kD, 0.0, 0.0, std::vector<double>(grid.values().begin(), grid.values().end()), zx,
      std::vector<double>(grid.dz_dy().begin(), grid.dz_dy().end()),
      std::vector<double>(grid.d2z_dxdy().begin(), grid.d2z_dxdy().end()));
  CHECK(std::abs(off.sag(7.3, -2.1) - f.z(7.3, -2.1)) > 100.0 * eps);
}

TEST_CASE("grid sag without derivatives converges with the third order (ADR 0051, point 9)",
          "[grid]") {
  // Finite differences err by O(d^2) in z_x and z_y and enter with the weight d G (|G| <= 1/4),
  // so the interpolant errs by O(d^3); halving the spacing on the same rectangle must reduce the
  // largest error at the same points by at least 6 (asymptotically 8).
  const Sphere f{30.0};
  std::mt19937_64 g(2026U);
  std::vector<std::pair<double, double>> points;
  for (int k = 0; k < 4000; ++k)
    points.emplace_back((2.0 * unit(g) - 1.0) * 10.0, (2.0 * unit(g) - 1.0) * 10.0);
  std::array<double, 3> err{};
  for (std::size_t level = 0; level < 3; ++level) {
    const int n = 20 * (1 << level) + 1;  // 21, 41, 81 points on [-10, 10]
    const double d = 20.0 / (n - 1);
    const GridSagTerm<double> grid = sampled(f, n, d, false);
    REQUIRE_FALSE(grid.derivatives_given());
    for (const auto& [x, y] : points)
      err[level] = std::max(err[level], std::abs(grid.sag(x, y) - f.z(x, y)));
  }
  INFO("errors " << err[0] << ", " << err[1] << ", " << err[2]);
  CHECK(err[0] / err[1] >= 6.0);
  CHECK(err[1] / err[2] >= 6.0);
}

TEST_CASE("finite differences are exact for a quadratic, edges included", "[grid]") {
  // Central and one-sided differences of second order are exact for polynomials of degree 2, so
  // the estimated derivatives equal the analytic ones up to rounding.
  const auto q = [](double x, double y) {
    return 0.01 * x * x - 0.02 * x * y + 0.005 * y * y + 0.1 * x;
  };
  constexpr int kN = 5;
  constexpr double kD = 2.0;
  std::vector<double> z;
  for (int j = 0; j < kN; ++j) {
    for (int i = 0; i < kN; ++i) z.push_back(q((i - 2) * kD + 1.0, (j - 2) * kD - 0.5));
  }
  const GridSagTerm<double> grid(kN, kN, kD, kD, 1.0, -0.5, z, {}, {}, {});
  for (int j = 0; j < kN; ++j) {
    for (int i = 0; i < kN; ++i) {
      const double x = (i - 2) * kD + 1.0;
      const double y = (j - 2) * kD - 0.5;
      const auto k = static_cast<std::size_t>(j * kN + i);
      INFO("i " << i << ", j " << j);
      CHECK(std::abs(grid.dz_dx()[k] - (0.02 * x - 0.02 * y + 0.1)) <= 1e-15);
      CHECK(std::abs(grid.dz_dy()[k] - (-0.02 * x + 0.01 * y)) <= 1e-15);
      CHECK(std::abs(grid.d2z_dxdy()[k] - -0.02) <= 1e-15);
    }
  }
}

TEST_CASE("grid values and derivatives are met, the gradient is continuous", "[grid]") {
  const Sphere f;
  constexpr int kN = 9;
  constexpr double kD = 1.5;
  const GridSagTerm<double> grid = sampled(f, kN, kD, true, 0.7, -0.4);
  for (int j = 0; j < kN; ++j) {
    for (int i = 0; i < kN; ++i) {
      const double x = 0.7 + (i - 4) * kD;
      const double y = -0.4 + (j - 4) * kD;
      INFO("i " << i << ", j " << j);
      // A few roundings of the cell position and of the basis at t = 0 or 1.
      CHECK(std::abs(grid.sag(x, y) - f.z(x, y)) <= 8.0 * kU * (1.0 + f.z(x, y)));
      const auto [gx, gy] = grid.grad(x, y);
      CHECK(std::abs(gx - f.zx(x, y)) <= 1e-13);
      CHECK(std::abs(gy - f.zy(x, y)) <= 1e-13);
    }
  }
  // Across an inner grid line the gradient changes by O(delta) only (C^1), while the second
  // derivative of the interpolant may jump.
  const double xl = 0.7 + kD;  // grid line i = 5
  for (const double y : {-3.1, 0.2, 2.9}) {
    const auto [l1, l2] = grid.grad(xl - 1e-9, y);
    const auto [r1, r2] = grid.grad(xl + 1e-9, y);
    CHECK(std::abs(l1 - r1) <= 1e-10);
    CHECK(std::abs(l2 - r2) <= 1e-10);
  }
}

TEST_CASE("grid sag domain: NaN outside, the last line belongs to the last cell", "[grid]") {
  const Sphere f;
  const GridSagTerm<double> grid =
      sampled(f, 5, 1.0, true, 1.0, 0.0);  // x in [-1, 3], y in [-2, 2]
  CHECK(grid.in_domain(-1.0, -2.0));
  CHECK(grid.in_domain(3.0, 2.0));
  CHECK_FALSE(grid.in_domain(3.0 + 1e-12, 0.0));
  CHECK_FALSE(grid.in_domain(0.0, -2.0 - 1e-12));
  CHECK_FALSE(grid.in_domain(std::numeric_limits<double>::quiet_NaN(), 0.0));
  CHECK(std::isnan(grid.sag(3.5, 0.0)));
  CHECK(std::isnan(grid.grad(0.0, 2.5).first));
  CHECK(std::abs(grid.sag(3.0, 2.0) - f.z(3.0, 2.0)) <= 1e-15);
  const double corner = std::hypot(3.0, 2.0);
  CHECK(grid.max_radius() == std::optional<double>(corner));
  const auto e = grid.vertex_expansion();
  CHECK(e.order == 1);
  CHECK_FALSE(e.rotational4);
  CHECK(std::abs(e.value) <= 1e-15);  // the sphere's vertex is a grid point here
  CHECK_THROWS_AS(
      GridSagTerm<double>(2, 5, 1.0, 1.0, 0.0, 0.0, std::vector<double>(10), {}, {}, {}),
      std::invalid_argument);
  CHECK_THROWS_AS(GridSagTerm<double>(3, 3, 1.0, 1.0, 0.0, 0.0, std::vector<double>(9),
                                      std::vector<double>(9), {}, {}),
                  std::invalid_argument);
}

TEST_CASE("intersection at the edge of the grid (ADR 0051, point 1)", "[grid][stacked]") {
  const Sphere f;
  const GridSagTerm<double> grid = sampled(f, 21, 1.0, true);  // [-10, 10]^2
  const Stacked<double> stack(Plane<double>{}, std::vector<TermVariant<double>>{grid});
  // Start and hit outside the grid: missed.
  const auto miss = rtt::geom::intersect_newton<double>(stack, rtt::math::Vec3(12.0, 0.0, -5.0),
                                                        rtt::math::Vec3(0.0, 0.0, 1.0));
  CHECK(miss.status == rtt::geom::HitStatus::Missed);
  // An oblique ray close to the edge whose hit lies inside: it crosses z = 0 at x = 9.95 and
  // meets the sag (about 0.5 mm) further inside.
  const rtt::math::Vec3 d = rtt::math::Vec3(-0.05, 0.0, 1.0).normalized();
  const rtt::math::Vec3 o = rtt::math::Vec3(9.95, 0.0, 0.0) - 5.0 * d;
  const auto hit = rtt::geom::intersect_newton<double>(stack, o, d);
  REQUIRE(hit.status == rtt::geom::HitStatus::Hit);
  const rtt::math::Vec3 p = o + hit.t * d;
  CHECK(p.x() < 9.95);
  CHECK(std::abs(p.z() - stack.sag(p.x(), p.y())) <= 1e-12);
}
