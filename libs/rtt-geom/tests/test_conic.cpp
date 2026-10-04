#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <utility>

#include "rtt/geom/conic.hpp"

using rtt::geom::Conic;

// Reference sag values are computed from closed forms that are independent of the
// implementation's formula. Source of the conic sag
//   z = c r^2 / (1 + sqrt(1 - (1 + k) c^2 r^2)):
// G. W. Forbes, Manufacturability estimates for optical aspheres, Opt. Express 19(10),
// 9923-9942 (2011), Eq. (2.1) with the definition of phi that follows it (see also
// W. T. Welford, Aberrations of Optical Systems, Ch. 2).
// Rearranged, every point of the surface satisfies c r^2 - 2 z + (1 + k) c z^2 = 0.

namespace {
constexpr double kTol = 1e-12;

bool near(double a, double b, double tol = kTol) {
  return std::abs(a - b) <= tol;
}

bool near_rel(double a, double b, double rel) {
  return std::abs(a - b) <= rel * std::max(std::abs(a), std::abs(b));
}
}  // namespace

TEST_CASE("sphere sag equals R - sqrt(R^2 - r^2)", "[conic][sag]") {
  // Sphere of radius R with centre on the axis at z = R (Forbes 2011, Eq. (2.1), k = 0).
  for (const double radius : {50.0, -50.0, 12.5}) {
    const Conic<double> sphere(1.0 / radius, 0.0);
    for (const auto& [x, y] :
         {std::pair{0.0, 0.0}, std::pair{3.0, -4.0}, std::pair{0.0, 10.0}, std::pair{-7.0, 2.0}}) {
      const double r2 = x * x + y * y;
      const double expected = radius - std::copysign(std::sqrt(radius * radius - r2), radius);
      REQUIRE(near(sphere.sag(x, y), expected));
    }
  }
}

TEST_CASE("paraboloid sag equals r^2 / (2 R)", "[conic][sag]") {
  // k = -1 removes the square root: z = c r^2 / 2 (Forbes 2011, Eq. (2.1)).
  const double radius = 100.0;
  const Conic<double> parabola(1.0 / radius, -1.0);
  // 1e-12 mm scaled with (1 + z): floating-point resolution at large z (ulp(5000) ~ 9e-13 mm).
  for (const double r : {0.0, 1.0, 20.0, 150.0, 1000.0}) {
    const double z = r * r / (2.0 * radius);
    REQUIRE(near(parabola.sag(0.0, r), z, 1e-12 * (1.0 + z)));
  }
}

TEST_CASE("hyperboloid sag equals (sqrt(1 + c^2 r^2) - 1) / c for k = -2", "[conic][sag]") {
  // For k = -2 the implicit form c r^2 - 2 z - c z^2 = 0 has the root below.
  const double c = 0.02;
  const Conic<double> hyperbola(c, -2.0);
  for (const double r : {0.0, 5.0, 30.0, 400.0}) {
    const double expected = (std::sqrt(1.0 + c * c * r * r) - 1.0) / c;
    REQUIRE(near(hyperbola.sag(r, 0.0), expected));
  }
}

TEST_CASE("ellipsoid sag satisfies the implicit conic equation", "[conic][sag]") {
  for (const double k : {0.5, 3.0, -0.4}) {
    const double c = 0.02;
    const Conic<double> conic(c, k);
    for (const double r : {1.0, 10.0, 0.9 / (c * std::sqrt(1.0 + k))}) {
      const double z = conic.sag(r, 0.0);
      REQUIRE(near(c * r * r - 2.0 * z + (1.0 + k) * c * z * z, 0.0));
    }
  }
}

TEST_CASE("conic gradient matches central finite difference", "[conic][grad]") {
  // Tolerance from the issue: relative 1e-7 against a central difference.
  const double h = 1e-5;
  for (const auto& [c, k] : {std::pair{0.02, 0.0}, std::pair{-0.03, 0.0}, std::pair{0.01, -1.0},
                             std::pair{0.02, -2.5}, std::pair{0.02, 0.7}}) {
    const Conic<double> conic(c, k);
    for (const auto& [x, y] : {std::pair{3.0, -4.0}, std::pair{10.0, 7.0}, std::pair{-12.0, 1.5}}) {
      const auto [gx, gy] = conic.grad(x, y);
      const double fx = (conic.sag(x + h, y) - conic.sag(x - h, y)) / (2.0 * h);
      const double fy = (conic.sag(x, y + h) - conic.sag(x, y - h)) / (2.0 * h);
      REQUIRE(near_rel(gx, fx, 1e-7));
      REQUIRE(near_rel(gy, fy, 1e-7));
    }
  }
}

TEST_CASE("conic gradient vanishes at the vertex", "[conic][grad]") {
  const Conic<double> sphere(0.02, 0.0);
  const auto [gx, gy] = sphere.grad(0.0, 0.0);
  REQUIRE(gx == 0.0);
  REQUIRE(gy == 0.0);
}

TEST_CASE("conic reports its base conic and domain", "[conic]") {
  const Conic<double> sphere(0.02, 0.0);
  REQUIRE(sphere.base_conic() == std::pair{0.02, 0.0});
  REQUIRE(sphere.curvature() == 0.02);
  REQUIRE(sphere.conic_constant() == 0.0);
  // Sphere: the surface ends at r = R.
  REQUIRE(sphere.max_radius().has_value());
  REQUIRE(near(*sphere.max_radius(), 50.0));
  // Ellipsoid: r_max = 1 / (|c| sqrt(1 + k)).
  const Conic<double> oblate(-0.02, 0.5);
  REQUIRE(near(*oblate.max_radius(), 1.0 / (0.02 * std::sqrt(1.5))));
  // Paraboloid, hyperboloid and the flat limit are unbounded.
  REQUIRE_FALSE(Conic<double>(0.02, -1.0).max_radius().has_value());
  REQUIRE_FALSE(Conic<double>(0.02, -2.0).max_radius().has_value());
  REQUIRE_FALSE(Conic<double>(0.0, 0.0).max_radius().has_value());
}

TEST_CASE("conic sag outside the domain is NaN, not a value of the other sheet", "[conic][sag]") {
  const Conic<double> sphere(0.02, 0.0);
  REQUIRE(std::isnan(sphere.sag(0.0, 50.5)));
  const auto [gx, gy] = sphere.grad(0.0, 50.5);
  REQUIRE(std::isnan(gx));
  REQUIRE(std::isnan(gy));
}
