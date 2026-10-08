#include <Eigen/Geometry>
#include <algorithm>
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <utility>
#include <vector>

#include "rtt/geom/conic.hpp"
#include "rtt/geom/phase.hpp"
#include "rtt/math/types.hpp"

// Phase functions of the phase layers and their gradient (#126, ADR 0025, point 1):
// LinearGrating phi = 2 pi G (x cos psi + y sin psi), RadialPhase phi = sum_k c_k rho^(2k).
// Reference values are hand values or closed forms independent of the implementation; the
// curved-surface case is checked against M. Mansuripur, Proc. SPIE 6620, 66200N (2007),
// Eqs. (9)-(11) (docs/quellen.md).

using rtt::geom::Conic;
using rtt::geom::LinearGratingPhase;
using rtt::geom::PhaseFunction;
using rtt::geom::RadialPhasePolynomial;
using rtt::math::Vec3;

namespace {

constexpr double kPi = std::numbers::pi;

bool near_rel(double a, double b, double rel) {
  return std::abs(a - b) <= rel * std::max(std::abs(a), std::abs(b));
}

bool near_abs(double a, double b, double tol) {
  return std::abs(a - b) <= tol;
}

/// R = 2 mm, c = (1, 0.5, -0.25) rad: the hand example of the plan for #126.
RadialPhasePolynomial<double> example_radial() {
  return {2.0, {1.0, 0.5, -0.25}};
}

}  // namespace

TEST_CASE("linear grating: |grad phi| = 2 pi G along the grating vector (#126)", "[phase]") {
  // ADR 0025, point 1: phi = 2 pi G (x cos psi + y sin psi), grating vector (cos psi, sin psi).
  // G = 300/mm, psi = 0: grad phi = (600 pi, 0) rad/mm exactly up to the rounding of 2 pi G.
  const LinearGratingPhase<double> g0(300.0, 0.0);
  REQUIRE(near_rel(g0.grad(0.0, 0.0).first, 600.0 * kPi, 1e-15));
  REQUIRE(g0.grad(0.0, 0.0).second == 0.0);
  REQUIRE(g0.phase(0.0, 0.0) == 0.0);
  REQUIRE(near_rel(g0.phase(2.0, 7.0), 1200.0 * kPi, 1e-15));  // independent of y
  // Rotated: |grad| = 2 pi G, direction (cos psi, sin psi). 1e-12 relative to |grad| covers
  // the rounding of cos and sin (not correctly rounded, a few ulp).
  for (const double psi_deg : {30.0, 90.0, -45.0, 135.0}) {
    INFO("psi = " << psi_deg << " deg");
    const double psi = psi_deg * kPi / 180.0;
    const LinearGratingPhase<double> g(300.0, psi);
    const double k = 600.0 * kPi;
    const auto [gx, gy] = g.grad(1.5, -2.5);
    REQUIRE(near_rel(std::hypot(gx, gy), k, 1e-12));
    REQUIRE(near_abs(gx, k * std::cos(psi), 1e-12 * k));
    REQUIRE(near_abs(gy, k * std::sin(psi), 1e-12 * k));
    REQUIRE(g.phase(0.0, 0.0) == 0.0);
    // phi is constant along the grooves, direction (-sin psi, cos psi). Points up to 10 mm:
    // |phi| <= 2 pi G 10 ~ 2e4 rad, rounding a few ulp of that, ~1e-11; tolerance 1e-10 rad.
    for (const double t : {-10.0, 3.0, 7.5}) {
      REQUIRE(near_abs(g.phase(1.0 - t * std::sin(psi), -2.0 + t * std::cos(psi)),
                       g.phase(1.0, -2.0), 1e-10));
    }
  }
  // psi = 90 deg at (2, 3): phi = 2 pi 300 * 3 = 1800 pi (cos psi is 6e-17, not 0).
  const LinearGratingPhase<double> g90(300.0, kPi / 2.0);
  REQUIRE(near_rel(g90.phase(2.0, 3.0), 1800.0 * kPi, 1e-12));
}

TEST_CASE("radial phase: hand values of phi and its gradient (#126)", "[phase]") {
  // phi = c1 u + c2 u^2 + c3 u^3, u = (x^2 + y^2) / R^2, R = 2, c = (1, 0.5, -0.25);
  // dphi/dx = dphi/du * 2 x / R^2, dphi/du = c1 + 2 c2 u + 3 c3 u^2.
  // (1, 1): u = 0.5, phi = 0.5 + 0.125 - 0.03125 = 0.59375,
  //         dphi/du = 1 + 0.5 - 0.1875 = 1.3125, grad = 1.3125 * 2 / 4 * (1, 1) = (0.65625,
  //         0.65625).
  // (2, 0): u = 1, phi = 1 + 0.5 - 0.25 = 1.25, dphi/du = 1 + 1 - 0.75 = 1.25,
  //         grad = (1.25 * 2 * 2 / 4, 0) = (1.25, 0).
  // All binary fractions: relative 1e-12 is far above the rounding of the few operations.
  const auto f = example_radial();
  REQUIRE(near_rel(f.phase(1.0, 1.0), 0.59375, 1e-12));
  REQUIRE(near_rel(f.grad(1.0, 1.0).first, 0.65625, 1e-12));
  REQUIRE(near_rel(f.grad(1.0, 1.0).second, 0.65625, 1e-12));
  REQUIRE(near_rel(f.phase(2.0, 0.0), 1.25, 1e-12));
  REQUIRE(near_rel(f.grad(2.0, 0.0).first, 1.25, 1e-12));
  REQUIRE(f.grad(2.0, 0.0).second == 0.0);
  REQUIRE(near_rel(f.grad(0.0, -2.0).second, -1.25, 1e-12));  // rotational symmetry
  // Vertex: phi = 0, gradient 0 (regular at r = 0).
  REQUIRE(f.phase(0.0, 0.0) == 0.0);
  REQUIRE(f.grad(0.0, 0.0) == std::pair{0.0, 0.0});
  // No coefficients: phi = 0 everywhere.
  const RadialPhasePolynomial<double> empty(5.0, {});
  REQUIRE(empty.phase(3.0, 4.0) == 0.0);
  REQUIRE(empty.grad(3.0, 4.0) == std::pair{0.0, 0.0});
}

TEST_CASE("radial phase: gradient equals the central difference of phi (#126)", "[phase]") {
  // Central difference with h = 1e-4 mm: truncation h^2 / 6 |phi'''| and rounding about
  // eps |phi| / h. For the example (R = 2, c = (1, 0.5, -0.25)) at r <= 1.6 mm,
  // |phi'''| < 3.2 rad/mm^3 (c2 term 0.5 r^4 / 16: 24 x 0.5 / 16 ~ 1.2; c3 term 0.25 r^6 / 64:
  // 120 x^3 0.25 / 64 ~ 1.9; with opposite signs even less), so truncation < 5.4e-9 and
  // rounding ~ 2e-12 rad/mm.
  // Tolerance 1e-7 rad/mm.
  const auto f = example_radial();
  constexpr double h = 1e-4;
  for (const auto& [x, y] : std::array{std::pair{0.7, -1.3}, std::pair{1.5, 0.4},
                                       std::pair{-0.2, 0.1}, std::pair{1.1, 1.1}}) {
    INFO("(" << x << ", " << y << ")");
    const auto [gx, gy] = f.grad(x, y);
    // Not trivially 0 = 0: |grad| = 2 r / R^2 dphi/du >= 0.05 rad/mm at these points
    // (r >= 0.22 mm, dphi/du = 1 + u - 0.75 u^2 >= 1 for u <= 0.6).
    REQUIRE(std::hypot(gx, gy) > 0.05);
    REQUIRE(near_abs(gx, (f.phase(x + h, y) - f.phase(x - h, y)) / (2.0 * h), 1e-7));
    REQUIRE(near_abs(gy, (f.phase(x, y + h) - f.phase(x, y - h)) / (2.0 * h), 1e-7));
  }
}

TEST_CASE("radial phase against Mansuripur's polynomial in r, Eq. (9) (#126)", "[phase]") {
  // Mansuripur (2007), Eq. (9): F(r) = sum_n a_n r^n in periods; phi = 2 pi F gives
  // a_n = c_(n/2) / (2 pi R^n) for even n. Example R = 2, c = (1, 0.5, -0.25):
  // a2 = 1 / (8 pi), a4 = 0.5 / (32 pi), a6 = -0.25 / (128 pi). Hand value at r = sqrt(2):
  // F = (0.25 + 0.0625 - 0.015625) / pi = 0.296875 / pi, 2 pi F = 0.59375 = phi(1, 1), and
  // 2 pi dF/dr = sqrt(2) (0.5 + 0.25 - 0.09375) = 0.65625 sqrt(2) = |grad phi(1, 1)|.
  const auto f = example_radial();
  const double a2 = 1.0 / (2.0 * kPi * 4.0);
  const double a4 = 0.5 / (2.0 * kPi * 16.0);
  const double a6 = -0.25 / (2.0 * kPi * 64.0);
  REQUIRE(near_rel(a2, 1.0 / (8.0 * kPi), 1e-15));
  REQUIRE(near_rel(2.0 * kPi * (a2 * 2.0 + a4 * 4.0 + a6 * 8.0), 0.59375, 1e-14));
  for (const double r : {std::sqrt(2.0), 0.3, 1.0, 2.5}) {
    INFO("r = " << r);
    const double big_f = a2 * std::pow(r, 2) + a4 * std::pow(r, 4) + a6 * std::pow(r, 6);
    const double df_dr = 2.0 * a2 * r + 4.0 * a4 * std::pow(r, 3) + 6.0 * a6 * std::pow(r, 5);
    const double alpha = 0.6;  // any azimuth
    const double x = r * std::cos(alpha);
    const double y = r * std::sin(alpha);
    REQUIRE(near_rel(f.phase(x, y), 2.0 * kPi * big_f, 1e-12));
    const auto [gx, gy] = f.grad(x, y);
    // Radial component of the gradient: dphi/dr = 2 pi dF/dr (signed).
    REQUIRE(near_rel(gx * std::cos(alpha) + gy * std::sin(alpha), 2.0 * kPi * df_dr, 1e-12));
    REQUIRE(near_abs(-gx * std::sin(alpha) + gy * std::cos(alpha), 0.0,
                     1e-12 * std::abs(2.0 * kPi * df_dr)));
  }
}

TEST_CASE("phase of a surface is the sum of its layers (#126)", "[phase]") {
  // ADR 0025, point 1. The sum starts at 0 and adds the layers in order, so it is bit-equal
  // to the same sum written out.
  const LinearGratingPhase<double> grating(150.0, 0.3);
  const auto radial = example_radial();
  const std::vector<PhaseFunction<double>> layers{grating, radial};
  const double x = 0.8;
  const double y = -0.6;
  // Closed form: 2 pi 150 (0.8 cos 0.3 - 0.6 sin 0.3) for the grating; u = 0.25 for the
  // radial layer, 0.25 + 0.5 * 0.0625 - 0.25 * 0.015625 = 0.27734375. Relative 1e-12.
  const double expected =
      2.0 * kPi * 150.0 * (0.8 * std::cos(0.3) - 0.6 * std::sin(0.3)) + 0.27734375;
  REQUIRE(near_rel(rtt::geom::phase<double>(layers, x, y), expected, 1e-12));
  REQUIRE(rtt::geom::phase<double>(layers, x, y) == 0.0 + grating.phase(x, y) + radial.phase(x, y));
  const auto [gx, gy] = rtt::geom::phase_grad<double>(layers, x, y);
  REQUIRE(gx == 0.0 + grating.grad(x, y).first + radial.grad(x, y).first);
  REQUIRE(gy == 0.0 + grating.grad(x, y).second + radial.grad(x, y).second);
  // No layers: phase 0, gradient 0.
  const std::vector<PhaseFunction<double>> none;
  REQUIRE(rtt::geom::phase<double>(none, x, y) == 0.0);
  REQUIRE(rtt::geom::phase_grad<double>(none, x, y) == std::pair{0.0, 0.0});
}

TEST_CASE("tangential gradient: projection into the tangent plane (#126)", "[phase]") {
  // ADR 0025, point 2: g_par = (I - N N^T) g with g = (gx, gy, 0).
  const std::pair g{3.0, -4.0};
  // Plane surface, N = +z or -z: g_par = g exactly.
  for (const double s : {1.0, -1.0}) {
    const Vec3 p = rtt::geom::tangential_gradient(g, Vec3(0.0, 0.0, s));
    REQUIRE(p == Vec3(3.0, -4.0, 0.0));
  }
  // Tilted normal: g_par . N = 0 and g_par . dr = g . dr for tangential dr (the surface
  // gradient). Relative 1e-12 of |g| = 5 rad/mm.
  const Vec3 n = Vec3(0.3, -0.4, 0.8).normalized();
  const Vec3 p = rtt::geom::tangential_gradient(g, n);
  REQUIRE(near_abs(p.dot(n), 0.0, 5e-12));
  const Vec3 lateral(3.0, -4.0, 0.0);
  for (const Vec3& e : {Vec3(1.0, 0.0, 0.0), Vec3(0.0, 1.0, 0.0), Vec3(0.2, 0.7, -0.1)}) {
    const Vec3 dr = n.cross(e).normalized();  // tangential
    REQUIRE(near_abs(p.dot(dr), lateral.dot(dr), 5e-12));
  }
}

TEST_CASE("tangential gradient on a surface of revolution equals Mansuripur Eq. (11) (#126)",
          "[phase]") {
  // Mansuripur (2007), Sec. 4, Eqs. (10), (11): on a surface of revolution with sag h(r),
  // dF/ds = (dF/dr) / sqrt(1 + (dh/dr)^2) along the meridional arc length s, i.e. along the
  // tangent t = (r_hat + h' z_hat) / sqrt(1 + h'^2). With phi = 2 pi F the projected gradient
  // must be (dphi/dr) / sqrt(1 + h'^2) * t. h' from the conic slope (rtt-geom), dphi/dr from
  // the radial phase; relative 1e-12.
  const RadialPhasePolynomial<double> f(10.0, {3.0, -1.5, 0.4});
  for (const auto& [c, k] : std::array{std::pair{1.0 / 50.0, 0.0}, std::pair{-1.0 / 30.0, -0.5},
                                       std::pair{1.0 / 25.0, 1.2}}) {
    const Conic<double> shape(c, k);
    for (const double r : {2.0, 7.5, 12.0}) {
      INFO("c = " << c << ", k = " << k << ", r = " << r);
      const double alpha = -1.1;
      const double x = r * std::cos(alpha);
      const double y = r * std::sin(alpha);
      const auto [zx, zy] = shape.grad(x, y);
      const Vec3 n = Vec3(-zx, -zy, 1.0).normalized();
      const Vec3 p = rtt::geom::tangential_gradient(f.grad(x, y), n);
      const Vec3 r_hat(std::cos(alpha), std::sin(alpha), 0.0);
      const double h_prime = zx * std::cos(alpha) + zy * std::sin(alpha);  // dh/dr
      // dphi/dr = sum_k c_k 2 k r^(2k-1) / R^(2k), independent of the implementation.
      const double dphi_dr = 3.0 * 2.0 * r / 100.0 - 1.5 * 4.0 * std::pow(r, 3) / 1e4 +
                             0.4 * 6.0 * std::pow(r, 5) / 1e6;
      REQUIRE(std::abs(dphi_dr) > 0.1);  // not trivially 0 = 0
      const double norm = std::sqrt(1.0 + h_prime * h_prime);
      const Vec3 t = (r_hat + Vec3(0.0, 0.0, h_prime)) / norm;
      const Vec3 expected = (dphi_dr / norm) * t;  // Eq. (11) along the tangent
      REQUIRE(near_rel(p.norm(), std::abs(dphi_dr) / norm, 1e-12));
      REQUIRE((p - expected).norm() <= 1e-12 * expected.norm());
    }
  }
}

TEST_CASE("phase functions reject invalid parameters (#126)", "[phase]") {
  const double inf = std::numeric_limits<double>::infinity();
  const double nan = std::numeric_limits<double>::quiet_NaN();
  REQUIRE_THROWS_AS(LinearGratingPhase<double>(inf, 0.0), std::invalid_argument);
  REQUIRE_THROWS_AS(LinearGratingPhase<double>(300.0, nan), std::invalid_argument);
  REQUIRE_NOTHROW(LinearGratingPhase<double>(0.0, 0.0));     // no grating
  REQUIRE_NOTHROW(LinearGratingPhase<double>(-300.0, 0.0));  // reversed grating vector
  REQUIRE_THROWS_AS(RadialPhasePolynomial<double>(0.0, {1.0}), std::invalid_argument);
  REQUIRE_THROWS_AS(RadialPhasePolynomial<double>(-2.0, {1.0}), std::invalid_argument);
  REQUIRE_THROWS_AS(RadialPhasePolynomial<double>(nan, {1.0}), std::invalid_argument);
  REQUIRE_THROWS_AS(RadialPhasePolynomial<double>(2.0, {1.0, inf}), std::invalid_argument);
  REQUIRE_NOTHROW(RadialPhasePolynomial<double>(2.0, {}));
}
