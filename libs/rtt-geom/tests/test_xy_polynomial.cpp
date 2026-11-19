// The XY polynomial term in rtt-geom (#316, ADR 0051, points 4 and 9, PR 3): the order of the
// terms after [UG] "Extended Polynomial", sag and gradient, the vertex expansion against the
// biconic and the stack.

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

#include "rtt/geom/anamorphic.hpp"
#include "rtt/geom/intersect.hpp"
#include "rtt/geom/plane.hpp"
#include "rtt/geom/stacked.hpp"
#include "rtt/geom/terms.hpp"

using rtt::geom::Biconic;
using rtt::geom::Plane;
using rtt::geom::Stacked;
using rtt::geom::TermVariant;
using rtt::geom::xy_polynomial_powers;
using rtt::geom::XyPolynomialTerm;

namespace {

double unit(std::mt19937_64& g) {
  return static_cast<double>(g() >> 11U) * 0x1.0p-53;
}

/// Coefficients 0 ... 0 with c_k = value at the 1-based positions given.
std::vector<double> only(std::size_t size,
                         std::initializer_list<std::pair<std::size_t, double>> c) {
  std::vector<double> a(size, 0.0);
  for (const auto& [k, v] : c) a[k - 1] = v;
  return a;
}

}  // namespace

TEST_CASE("XY polynomial order after [UG] Extended Polynomial (ADR 0051, point 4)", "[xy]") {
  // [UG]: x, y, x^2, xy, y^2, ...; term 12 is x^2 y^2, term 13 is x y^3; 230 terms to degree 20.
  CHECK(xy_polynomial_powers(1) == std::pair{1, 0});
  CHECK(xy_polynomial_powers(2) == std::pair{0, 1});
  CHECK(xy_polynomial_powers(3) == std::pair{2, 0});
  CHECK(xy_polynomial_powers(4) == std::pair{1, 1});
  CHECK(xy_polynomial_powers(5) == std::pair{0, 2});
  CHECK(xy_polynomial_powers(6) == std::pair{3, 0});
  CHECK(xy_polynomial_powers(12) == std::pair{2, 2});
  CHECK(xy_polynomial_powers(13) == std::pair{1, 3});
  CHECK(xy_polynomial_powers(14) == std::pair{0, 4});
  CHECK(xy_polynomial_powers(15) == std::pair{5, 0});
  CHECK(xy_polynomial_powers(230) == std::pair{0, 20});
  // Every (i, j) with 1 <= i + j <= 20 exactly once, degree ascending, x power descending.
  int k = 0;
  for (int d = 1; d <= 20; ++d) {
    for (int i = d; i >= 0; --i) CHECK(xy_polynomial_powers(++k) == std::pair{i, d - i});
  }
  CHECK(k == rtt::geom::kMaxXyTerms);
  CHECK_THROWS_AS((void)xy_polynomial_powers(0), std::invalid_argument);
  CHECK_THROWS_AS((void)xy_polynomial_powers(231), std::invalid_argument);
  CHECK_THROWS_AS(XyPolynomialTerm<double>(10.0, std::vector<double>(231, 0.0)),
                  std::invalid_argument);
  CHECK_THROWS_AS(XyPolynomialTerm<double>(0.0, {1.0}), std::invalid_argument);
}

TEST_CASE("XY polynomial sag is the sum of the monomials", "[xy]") {
  constexpr double kR = 8.0;
  const XyPolynomialTerm<double> t(kR, only(13, {{1, 1e-3}, {4, -2e-3}, {12, 5e-4}, {13, 3e-4}}));
  for (const auto& [x, y] : std::array<std::pair<double, double>, 4>{
           {{0.0, 0.0}, {2.0, -3.0}, {-7.5, 1.25}, {12.0, 9.0}}}) {
    INFO("(" << x << ", " << y << ")");
    const double u = x / kR;
    const double v = y / kR;
    const double expected = 1e-3 * u - 2e-3 * u * v + 5e-4 * u * u * v * v + 3e-4 * u * v * v * v;
    CHECK(std::abs(t.sag(x, y) - expected) <= 1e-15 * (1.0 + std::abs(expected)));
  }
  CHECK(t.sag(0.0, 0.0) == 0.0);  // no constant term
}

TEST_CASE("XY polynomial gradient against central differences", "[xy]") {
  constexpr double kR = 5.0;
  std::mt19937_64 g(316U);
  std::vector<double> a(65);  // up to degree 10
  for (double& v : a) v = (2.0 * unit(g) - 1.0) * 1e-3;
  const XyPolynomialTerm<double> t(kR, a);
  for (const auto& [x, y] : std::array<std::pair<double, double>, 5>{
           {{0.0, 0.0}, {1.0, 2.0}, {-3.0, 0.5}, {4.9, -0.1}, {6.0, 3.0}}}) {
    INFO("(" << x << ", " << y << ")");
    const double h = 1e-5;
    const double dx = (t.sag(x + h, y) - t.sag(x - h, y)) / (2.0 * h);
    const double dy = (t.sag(x, y + h) - t.sag(x, y - h)) / (2.0 * h);
    const auto [gx, gy] = t.grad(x, y);
    CHECK(std::abs(gx - dx) <= 1e-7);
    CHECK(std::abs(gy - dy) <= 1e-7);
  }
}

TEST_CASE("XY polynomial vertex expansion (ADR 0051, points 4 and 9)", "[xy]") {
  constexpr double kR = 10.0;
  SECTION("x^2 and y^2 have the Hessian of the biconic with c = 2 c_k / R^2") {
    // ADR 0051, point 9: the test takes R and c so that 1 / radius == 2 c / (R R), then the
    // Hessians agree bit for bit.
    const double cx2 = 0.25;
    const double cy2 = -0.125;
    const double rx = (kR * kR) / (2.0 * cx2);  // 200
    const double ry = (kR * kR) / (2.0 * cy2);  // -400
    REQUIRE(1.0 / rx == 2.0 * cx2 / (kR * kR));
    REQUIRE(1.0 / ry == 2.0 * cy2 / (kR * kR));
    const auto e = XyPolynomialTerm<double>(kR, only(5, {{3, cx2}, {5, cy2}})).vertex_expansion();
    const auto b = Biconic<double>(1.0 / rx, 1.0 / ry, 0.0, 0.0).vertex_expansion();
    CHECK(e.hessian == b.hessian);
    CHECK(e.value == 0.0);
    CHECK(e.slope == std::pair{0.0, 0.0});
    CHECK_FALSE(e.rotational4);
    CHECK(e.order == 4);
  }
  SECTION("tilt and the mixed term") {
    const auto e = XyPolynomialTerm<double>(kR, only(4, {{1, 2e-3}, {2, -1e-3}, {4, 3e-3}}))
                       .vertex_expansion();
    CHECK(e.slope.first == 2e-3 / kR);
    CHECK(e.slope.second == -1e-3 / kR);
    CHECK(e.hessian[1] == 3e-3 / (kR * kR));
    CHECK_FALSE(e.rotational4);
  }
  SECTION("(x^2 + y^2) a + (x^2 + y^2)^2 b is rotational with b4 = b / R^4") {
    const double a = 2e-3;
    const double b = -5e-4;
    const auto e =
        XyPolynomialTerm<double>(kR, only(14, {{3, a}, {5, a}, {10, b}, {12, 2.0 * b}, {14, b}}))
            .vertex_expansion();
    CHECK(e.rotational4);
    const double c = 2.0 * a / (kR * kR);
    CHECK(e.hessian[0] == c);
    CHECK(e.hessian[2] == c);
    CHECK(e.hessian[1] == 0.0);
    CHECK(e.a4 == b / (kR * kR * kR * kR) - c * c * c / 8.0);
  }
  SECTION("each break of the symmetry to fourth order clears rotational4") {
    // The rotational set of the section above, then one coefficient changed.
    for (std::size_t k = 1; k <= 14; ++k) {
      INFO("term " << k);
      std::vector<double> c = only(14, {{3, 2e-3}, {5, 2e-3}, {10, 1e-3}, {12, 2e-3}, {14, 1e-3}});
      c[k - 1] += 1e-4;
      CHECK_FALSE(XyPolynomialTerm<double>(kR, c).vertex_expansion().rotational4);
    }
    // Degree 5 and higher do not enter the fourth order.
    const std::vector<double> c = only(20, {{3, 2e-3}, {5, 2e-3}, {15, 1e-3}, {20, -1e-3}});
    CHECK(XyPolynomialTerm<double>(kR, c).vertex_expansion().rotational4);
  }
}

TEST_CASE("a stack with an XY polynomial is intersected by Newton (ADR 0051, point 1)",
          "[xy][stacked]") {
  const XyPolynomialTerm<double> t(10.0, only(5, {{1, 1e-2}, {3, 0.2}, {4, -0.05}, {5, 0.1}}));
  const Stacked<double> stack(Plane<double>{}, std::vector<TermVariant<double>>{t});
  CHECK(stack.sag(1.0, 2.0) == t.sag(1.0, 2.0));
  CHECK(stack.grad(1.0, 2.0) == t.grad(1.0, 2.0));
  const rtt::math::Vec3 o(1.0, 2.0, -1.0);
  const rtt::math::Vec3 d(0.0, 0.0, 1.0);
  const auto hit = rtt::geom::intersect_newton<double>(stack, o, d);
  REQUIRE(hit.status == rtt::geom::HitStatus::Hit);
  CHECK(std::abs(hit.t - (1.0 + stack.sag(1.0, 2.0))) <= 1e-12);
  CHECK(stack.sag(1.0, 2.0) != 0.0);  // content guard
}
