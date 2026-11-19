// Zernike terms and the shape stack in rtt-geom (#316, ADR 0051, points 1 to 3 and 9, PR 2).

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numbers>
#include <optional>
#include <random>
#include <vector>

#include "rtt/geom/intersect.hpp"
#include "rtt/geom/stacked.hpp"
#include "rtt/geom/terms.hpp"
#include "rtt/math/zernike.hpp"

using rtt::geom::Conic;
using rtt::geom::Plane;
using rtt::geom::Stacked;
using rtt::geom::TermVariant;
using rtt::geom::ZernikeTerm;

namespace {

constexpr double kU = std::numeric_limits<double>::epsilon() / 2.0;  // unit roundoff

/// gamma_k = k u / (1 - k u) (Higham, as in ADR 0032, point 6).
double gamma(int k) {
  return k * kU / (1.0 - k * kU);
}

double unit(std::mt19937_64& g) {
  return static_cast<double>(g() >> 11U) * 0x1.0p-53;
}

ZernikeTerm<double> noll(double radius, const std::vector<double>& a) {
  return ZernikeTerm<double>::noll(radius, a);
}

}  // namespace

TEST_CASE("Fringe and Noll give the same sag (ADR 0051, point 9)", "[zernike]") {
  // Random Fringe coefficients b_k; the Noll coefficient of the same (n, m, cos/sin) is
  // b_k / N_nm (Z_Noll = N_nm Z_Fringe, Niu/Tian Eq. (3) against (22)). Tolerance: the rounding
  // of both sums, (gamma_37 + gamma_79 + 4u) sum_k |b_k Z_k|, with the same radial polynomial.
  constexpr double kR = 7.5;
  std::mt19937_64 g(316U);
  std::vector<double> fringe(37);
  for (double& b : fringe) b = (2.0 * unit(g) - 1.0) * 1e-3;
  std::vector<double> noll_a(79, 0.0);  // (12, 0) is Noll j = 79
  for (int k = 1; k <= 37; ++k) {
    const rtt::math::FringeIndex f = rtt::math::fringe_to_nm(k);
    const int j = rtt::math::noll_from_nm(f.n, f.m, f.sine);
    noll_a[static_cast<std::size_t>(j - 1)] =
        fringe[static_cast<std::size_t>(k - 1)] / rtt::math::zernike_noll_factor(f.n, f.m);
  }
  const ZernikeTerm<double> zf = ZernikeTerm<double>::fringe(kR, fringe);
  const ZernikeTerm<double> zn = noll(kR, noll_a);
  // The unconverted coefficients (Noll with the Fringe values) as content guard.
  std::vector<double> wrong(79, 0.0);
  for (int k = 1; k <= 37; ++k) {
    const rtt::math::FringeIndex f = rtt::math::fringe_to_nm(k);
    wrong[static_cast<std::size_t>(rtt::math::noll_from_nm(f.n, f.m, f.sine) - 1)] =
        fringe[static_cast<std::size_t>(k - 1)];
  }
  const ZernikeTerm<double> zw = noll(kR, wrong);
  double worst_ratio = 0.0;
  double guard_ratio = std::numeric_limits<double>::infinity();
  for (int i = 0; i < 10000; ++i) {
    const double r = kR * std::sqrt(unit(g));
    const double phi = 2.0 * std::numbers::pi * unit(g);
    const double x = r * std::cos(phi);
    const double y = r * std::sin(phi);
    double scale = 0.0;  // sum_k |b_k Z_k|
    for (int k = 1; k <= 37; ++k) {
      const rtt::math::FringeIndex f = rtt::math::fringe_to_nm(k);
      scale += std::abs(fringe[static_cast<std::size_t>(k - 1)] *
                        rtt::math::zernike_radial(f.n, f.m, r / kR));
    }
    const double tol = (gamma(37) + gamma(79) + 4.0 * kU) * scale;
    const double d = std::abs(zf.sag(x, y) - zn.sag(x, y));
    worst_ratio = std::max(worst_ratio, d / tol);
    guard_ratio = std::min(guard_ratio, std::abs(zw.sag(x, y) - zf.sag(x, y)) / tol);
  }
  CHECK(worst_ratio <= 1.0);
  CHECK(guard_ratio > 1e6);  // content guard: the normalisation factors matter
}

TEST_CASE("vertex expansion of Zernike terms (ADR 0051, point 9)", "[zernike]") {
  constexpr double kR = 10.0;
  constexpr double kA = 2e-3;
  SECTION("Z4 defocus: value -sqrt(3) a, curvature 4 sqrt(3) a / R^2, rotational") {
    const auto e = noll(kR, {0.0, 0.0, 0.0, kA}).vertex_expansion();
    // The same expressions as the code: coefficient sqrt(3) a times c_K = -1, and
    // 2 coefficient c_{K-1} / R^2 with c_{K-1} = 2.
    const double coef = std::sqrt(3.0) * kA;
    CHECK(e.value == coef * -1.0);
    CHECK(e.hessian[0] == 2.0 * coef * 2.0 / (kR * kR));
    CHECK(e.hessian[2] == e.hessian[0]);
    CHECK(e.hessian[1] == 0.0);
    CHECK(e.rotational4);
    CHECK(e.order == 4);
    CHECK(std::abs(e.hessian[0] - 4.0 * std::sqrt(3.0) * kA / (kR * kR)) <= 1e-18);
  }
  SECTION("Z6 astigmatism (cos 2 theta): Hessian not isotropic, not rotational") {
    const auto e = noll(kR, {0.0, 0.0, 0.0, 0.0, 0.0, kA}).vertex_expansion();
    CHECK(e.hessian[0] == -e.hessian[2]);
    CHECK(e.hessian[0] > 0.0);
    CHECK_FALSE(e.rotational4);
  }
  SECTION("Z2 tilt x: slope a * 2 / R") {
    const auto e = noll(kR, {0.0, kA}).vertex_expansion();
    CHECK(e.slope.first == 2.0 * kA / kR);  // Noll factor sqrt(4) = 2, R_1^1 = rho
    CHECK(e.slope.second == 0.0);
    CHECK_FALSE(e.rotational4);
  }
  SECTION("Z28 (n = 6, m = 6) alone stays rotational to fourth order") {
    std::vector<double> a(28, 0.0);
    a[27] = kA;
    REQUIRE(rtt::math::noll_to_nm(28).m == 6);
    const auto e = noll(kR, a).vertex_expansion();
    CHECK(e.rotational4);
    CHECK(e.value == 0.0);
    CHECK(e.hessian[0] == 0.0);
  }
  SECTION("Z11 spherical: b4 = sqrt(5) a 6 / R^4 enters A4") {
    std::vector<double> a(11, 0.0);
    a[10] = kA;
    const auto e = noll(kR, a).vertex_expansion();
    const double coef = std::sqrt(5.0) * kA;  // R_4^0 = 6 rho^4 - 6 rho^2 + 1
    const double c = 2.0 * coef * -6.0 / (kR * kR);
    CHECK(e.hessian[0] == c);
    CHECK(e.a4 == coef * 6.0 / (kR * kR * kR * kR) - c * c * c / 8.0);
    CHECK(e.value == coef * 1.0);
  }
}

TEST_CASE("Zernike gradient against central differences (ADR 0051, point 2)", "[zernike]") {
  constexpr double kR = 5.0;
  std::mt19937_64 g(2026U);
  std::vector<double> a(45);
  for (double& v : a) v = (2.0 * unit(g) - 1.0) * 1e-3;
  std::vector<double> f(37);
  for (double& v : f) v = (2.0 * unit(g) - 1.0) * 1e-3;
  for (const ZernikeTerm<double>& z : {noll(kR, a), ZernikeTerm<double>::fringe(kR, f)}) {
    for (const auto& [x, y] : std::array<std::pair<double, double>, 6>{
             {{0.0, 0.0}, {1e-7, -2e-7}, {1.0, 2.0}, {-3.0, 0.5}, {4.9, -0.1}, {6.0, 3.0}}}) {
      INFO("(" << x << ", " << y << ")");
      const double h = 1e-5;
      const double dx = (z.sag(x + h, y) - z.sag(x - h, y)) / (2.0 * h);
      const double dy = (z.sag(x, y + h) - z.sag(x, y - h)) / (2.0 * h);
      const auto [gx, gy] = z.grad(x, y);
      CHECK(std::abs(gx - dx) <= 1e-7);
      CHECK(std::abs(gy - dy) <= 1e-7);
    }
  }
}

TEST_CASE("a stack sums base and terms and is intersected by Newton (ADR 0051, point 1)",
          "[zernike][stacked]") {
  constexpr double kR = 10.0;
  const Conic<double> conic(1.0 / 50.0, -0.5);
  const ZernikeTerm<double> defocus = noll(kR, {0.0, 0.0, 0.0, 1e-3});
  const Stacked<double> stack(conic, std::vector<TermVariant<double>>{defocus});
  CHECK(stack.sag(1.0, 2.0) == conic.sag(1.0, 2.0) + defocus.sag(1.0, 2.0));
  const auto [gx, gy] = stack.grad(1.0, 2.0);
  CHECK(gx == conic.grad(1.0, 2.0).first + defocus.grad(1.0, 2.0).first);
  CHECK(gy == conic.grad(1.0, 2.0).second + defocus.grad(1.0, 2.0).second);
  CHECK(stack.base_conic() == conic.base_conic());
  CHECK(stack.max_radius() == conic.max_radius());
  // in_domain: the disk of the conic's max_radius (no own domain), terms unbounded.
  CHECK(stack.in_domain(1.0, 2.0));
  const std::optional<double> r_max = conic.max_radius();
  REQUIRE(r_max.has_value());
  CHECK_FALSE(stack.in_domain(r_max.value_or(0.0) * 1.01, 0.0));
  // Vertex expansion: the curvatures add, the cap is K = 0, A4 = b4 - c^3 / 8.
  const auto e = stack.vertex_expansion();
  CHECK(e.hessian[0] ==
        conic.vertex_expansion().hessian[0] + defocus.vertex_expansion().hessian[0]);
  CHECK(e.rotational4);
  CHECK(e.conic == 0.0);
  // A ray along +z from z = -1 at (1, 2) hits at z = sag(1, 2).
  const rtt::math::Vec3 o(1.0, 2.0, -1.0);
  const rtt::math::Vec3 d(0.0, 0.0, 1.0);
  const auto hit = rtt::geom::intersect_newton<double>(stack, o, d);
  REQUIRE(hit.status == rtt::geom::HitStatus::Hit);
  CHECK(std::abs(hit.t - (1.0 + stack.sag(1.0, 2.0))) <= 1e-12);
  // Over a plane the stack is the term alone.
  const Stacked<double> flat(Plane<double>{}, std::vector<TermVariant<double>>{defocus});
  CHECK(flat.sag(3.0, -1.0) == defocus.sag(3.0, -1.0));
  CHECK_FALSE(flat.max_radius().has_value());
}
