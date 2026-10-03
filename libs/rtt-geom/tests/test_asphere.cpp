#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <utility>
#include <vector>

#include "rtt/geom/asphere.hpp"
#include "rtt/geom/conic.hpp"

using rtt::geom::Conic;
using rtt::geom::EvenAsphere;

// Even asphere z = conic(r) + A4 r^4 + A6 r^6 + ... with the conic sag of
// W. T. Welford, Aberrations of Optical Systems, Ch. 2; coefficient order as in the model
// (rtt/model/surface.hpp, coefficients[0] = A4).

namespace {
bool near(double a, double b, double tol) {
  return std::abs(a - b) <= tol;
}

bool near_rel(double a, double b, double rel) {
  return std::abs(a - b) <= rel * std::max(std::abs(a), std::abs(b));
}
}  // namespace

TEST_CASE("asphere without coefficients equals the conic", "[asphere]") {
  // Issue #4: identical sag and gradient, absolute 1e-15 mm.
  for (const auto& [c, k] : {std::pair{0.02, 0.0}, std::pair{-0.03, -1.0}, std::pair{0.01, -2.5},
                             std::pair{0.025, 0.6}}) {
    const Conic<double> conic(c, k);
    for (const auto& coefficients : {std::vector<double>{}, std::vector<double>{0.0, 0.0, 0.0}}) {
      const EvenAsphere<double> asphere(c, k, coefficients);
      for (const auto& [x, y] :
           {std::pair{0.0, 0.0}, std::pair{3.0, -4.0}, std::pair{-12.0, 7.5}}) {
        REQUIRE(near(asphere.sag(x, y), conic.sag(x, y), 1e-15));
        const auto [ax, ay] = asphere.grad(x, y);
        const auto [cx, cy] = conic.grad(x, y);
        REQUIRE(near(ax, cx, 1e-15));
        REQUIRE(near(ay, cy, 1e-15));
      }
      REQUIRE(asphere.base_conic() == conic.base_conic());
      REQUIRE(asphere.max_radius() == conic.max_radius());
    }
  }
}

TEST_CASE("asphere sag and gradient equal the hand-evaluated polynomial", "[asphere]") {
  // Paraboloid base (k = -1, conic part c r^2 / 2) with c = 0.01 /mm, A4 = 1e-6 /mm^3,
  // A6 = -2e-9 /mm^5. At (x, y) = (3, 4), r = 5:
  //   z     = 0.01 * 25 / 2 + 1e-6 * 625 - 2e-9 * 15625 = 0.12559375 mm
  //   dz/dr = c r + 4 A4 r^3 + 6 A6 r^5 = 0.05 + 0.0005 - 0.0000375 = 0.0504625
  //   dz/dx = dz/dr * x / r = 0.0302775, dz/dy = dz/dr * y / r = 0.04037
  const EvenAsphere<double> asphere(0.01, -1.0, {1e-6, -2e-9});
  REQUIRE(near(asphere.sag(3.0, 4.0), 0.12559375, 1e-15));
  const auto [gx, gy] = asphere.grad(3.0, 4.0);
  REQUIRE(near(gx, 0.0302775, 1e-15));
  REQUIRE(near(gy, 0.04037, 1e-15));
  // At r = 10 (on the y axis): z = 0.5 + 0.01 - 0.002 = 0.508 mm,
  // dz/dy = 0.1 + 0.004 - 0.0012 = 0.1028, dz/dx = 0.
  REQUIRE(near(asphere.sag(0.0, 10.0), 0.508, 1e-15));
  const auto [hx, hy] = asphere.grad(0.0, 10.0);
  REQUIRE(hx == 0.0);
  REQUIRE(near(hy, 0.1028, 1e-15));
}

TEST_CASE("asphere coefficients are A4, A6, A8 in this order", "[asphere]") {
  // Flat base, single coefficient each: z = A r^(2 i + 4) for coefficients[i] = A.
  const double r = 2.0;
  REQUIRE(near(EvenAsphere<double>(0.0, 0.0, {1.0}).sag(r, 0.0), 16.0, 1e-15));
  REQUIRE(near(EvenAsphere<double>(0.0, 0.0, {0.0, 1.0}).sag(r, 0.0), 64.0, 1e-15));
  REQUIRE(near(EvenAsphere<double>(0.0, 0.0, {0.0, 0.0, 1.0}).sag(r, 0.0), 256.0, 1e-15));
}

TEST_CASE("asphere gradient matches central finite difference", "[asphere]") {
  const double h = 1e-5;
  const EvenAsphere<double> asphere(1.0 / 40.0, -0.5, {-2e-6, 1e-9, -3e-12, 1e-15});
  for (const auto& [x, y] :
       {std::pair{3.0, -4.0}, std::pair{10.0, 7.0}, std::pair{-12.0, 1.5}, std::pair{0.0, 18.0}}) {
    const auto [gx, gy] = asphere.grad(x, y);
    const double fx = (asphere.sag(x + h, y) - asphere.sag(x - h, y)) / (2.0 * h);
    const double fy = (asphere.sag(x, y + h) - asphere.sag(x, y - h)) / (2.0 * h);
    if (x != 0.0) {
      REQUIRE(near_rel(gx, fx, 1e-7));
    } else {
      REQUIRE(near(gx, fx, 1e-12));
    }
    REQUIRE(near_rel(gy, fy, 1e-7));
  }
}

TEST_CASE("asphere is NaN outside the domain of its base conic", "[asphere]") {
  const EvenAsphere<double> asphere(0.02, 0.0, {1e-6});
  REQUIRE(asphere.max_radius().has_value());
  REQUIRE(std::isnan(asphere.sag(0.0, 50.5)));
  const auto [gx, gy] = asphere.grad(0.0, 50.5);
  REQUIRE(std::isnan(gx));
  REQUIRE(std::isnan(gy));
}
