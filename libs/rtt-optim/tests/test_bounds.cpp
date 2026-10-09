// Bound transformations of the optimizer (ADR 0030, point 5). Source: F. James, M. Winkler,
// "MINUIT User's Guide", CERN 2004, section 1.3.1, eqs. (1.1)-(1.6) (docs/quellen.md).

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <limits>
#include <numbers>

#include "rtt/optim/bounds.hpp"

using rtt::optim::at_bound;
using rtt::optim::Bounds;
using rtt::optim::to_external;
using rtt::optim::to_internal;

namespace {

// Relative difference in units of eps_M.
double ulps(double a, double b) {
  return std::abs(a - b) / (std::numeric_limits<double>::epsilon() * std::max(std::abs(b), 1.0));
}

}  // namespace

TEST_CASE("bounds: without bounds the transformation is the identity", "[optim][bounds]") {
  const Bounds none;
  for (const double p : {-3.5, 0.0, 1e-300, 2.0, 7e12}) {
    CHECK(to_internal(p, none) == p);
    CHECK(to_external(p, none) == p);
    CHECK_FALSE(at_bound(p, none));
  }
}

TEST_CASE("bounds: two bounds follow MINUIT eqs. (1.1) and (1.2)", "[optim][bounds]") {
  const Bounds b{2.0, 6.0};
  // (1.1): theta = arcsin(2 (p - a)/(b - a) - 1): -pi/2 at a, 0 in the middle, pi/2 at b.
  CHECK(ulps(to_internal(2.0, b), -std::numbers::pi / 2.0) <= 1.0);
  CHECK(to_internal(4.0, b) == 0.0);
  CHECK(ulps(to_internal(6.0, b), std::numbers::pi / 2.0) <= 1.0);
  // (1.2): p = a + (b - a)/2 (sin theta + 1).
  CHECK(to_external(0.0, b) == 4.0);
  CHECK(ulps(to_external(std::numbers::pi / 6.0, b), 2.0 + 2.0 * 1.5) < 4.0);  // sin = 1/2
  // Periodic in theta and never outside [a, b].
  for (const double t : {-10.0, -1.0, 0.3, 1.5707963267948966, 2.0, 25.0}) {
    const double p = to_external(t, b);
    CHECK(p >= 2.0);
    CHECK(p <= 6.0);
  }
}

TEST_CASE("bounds: a lower bound follows MINUIT eqs. (1.3) and (1.4)", "[optim][bounds]") {
  const Bounds b{1.0, std::nullopt};
  // (1.3): theta = sqrt((p - a + 1)^2 - 1); p - a + 1 = sqrt(2) gives theta = 1 (branch >= 0).
  const double p1 = 1.0 + std::numbers::sqrt2 - 1.0;
  CHECK(ulps(to_internal(p1, b), 1.0) < 4.0);
  CHECK(to_internal(1.0, b) == 0.0);
  // (1.4): p = a - 1 + sqrt(theta^2 + 1); theta and -theta give the same p.
  CHECK(ulps(to_external(1.0, b), p1) < 4.0);
  CHECK(to_external(-1.0, b) == to_external(1.0, b));
  CHECK(to_external(0.0, b) == 1.0);
  for (const double t : {-30.0, -1e-3, 0.0, 1e-9, 0.5, 40.0}) CHECK(to_external(t, b) >= 1.0);
}

TEST_CASE("bounds: an upper bound follows MINUIT eqs. (1.5) and (1.6)", "[optim][bounds]") {
  const Bounds b{std::nullopt, -2.0};
  // (1.5): theta = sqrt((b - p + 1)^2 - 1); b - p + 1 = sqrt(2) gives theta = 1.
  const double p1 = -2.0 + 1.0 - std::numbers::sqrt2;
  CHECK(ulps(to_internal(p1, b), 1.0) < 4.0);
  CHECK(to_internal(-2.0, b) == 0.0);
  // (1.6): p = b + 1 - sqrt(theta^2 + 1).
  CHECK(ulps(to_external(1.0, b), p1) < 4.0);
  CHECK(to_external(-1.0, b) == to_external(1.0, b));
  for (const double t : {-30.0, 0.0, 1e-9, 40.0}) CHECK(to_external(t, b) <= -2.0);
}

TEST_CASE("bounds: the round trip p -> theta -> p is exact to a few ulp", "[optim][bounds]") {
  // Small values near a one-sided bound (natural size 1e-4, ADR 0030 point 14) included: the
  // quadratic part of the transformation must not lose the value to cancellation.
  const Bounds lower{0.0, std::nullopt};
  const Bounds upper{std::nullopt, 0.0};
  const Bounds both{-1.0, 3.0};
  for (const double p : {1e-12, 1e-4, 5e-4, 0.3, 1.0, 2.5, 1e3}) {
    INFO("p = " << p);
    CHECK(ulps(to_external(to_internal(p, lower), lower), p) * std::max(1.0, 1.0 / p) < 16.0);
    CHECK(ulps(to_external(to_internal(-p, upper), upper), -p) * std::max(1.0, 1.0 / p) < 16.0);
  }
  for (const double p : {-0.999, -0.5, 0.0, 1.0, 2.75, 2.999}) {
    INFO("p = " << p);
    CHECK(ulps(to_external(to_internal(p, both), both), p) < 16.0);
  }
}

TEST_CASE("bounds: at_bound uses the thresholds of ADR 0030", "[optim][bounds]") {
  const Bounds both{0.0, 10.0};  // 1e-6 (b - a) = 1e-5
  CHECK(at_bound(0.0, both));
  CHECK(at_bound(0.9e-5, both));
  CHECK_FALSE(at_bound(1.1e-5, both));
  CHECK(at_bound(10.0 - 0.9e-5, both));
  CHECK_FALSE(at_bound(5.0, both));
  const Bounds lower{1000.0, std::nullopt};  // 1e-6 max(1, |a|) = 1e-3
  CHECK(at_bound(1000.0009, lower));
  CHECK_FALSE(at_bound(1000.0011, lower));
  const Bounds small{0.0, std::nullopt};  // 1e-6 max(1, 0) = 1e-6
  CHECK(at_bound(0.9e-6, small));
  CHECK_FALSE(at_bound(1.1e-6, small));
  const Bounds upper{std::nullopt, -5.0};  // 1e-6 max(1, 5) = 5e-6
  CHECK(at_bound(-5.0 - 4e-6, upper));
  CHECK_FALSE(at_bound(-5.0 - 6e-6, upper));
}
