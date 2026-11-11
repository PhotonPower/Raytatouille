// Bessel functions J0 and J1 (rtt/math/bessel.hpp, #373) against the power series of NIST DLMF,
// Version 1.2.8, Eq. (10.2.2): J_nu(x) = (x/2)^nu sum_k (-1)^k (x^2/4)^k / (k! Gamma(nu+k+1)).

#include <catch2/catch_test_macros.hpp>
#include <cmath>

#include "rtt/math/bessel.hpp"

namespace {

// Partial sum of (10.2.2) for integer nu = 0, 1; 60 terms leave a remainder far below 1e-16 for
// x <= 6. Its rounding error is below 1e-15 absolute here (checked against exact rational
// arithmetic: at most 5.3e-16, J1 at x = 6). The bound 1e-13 leaves room for the library
// (libstdc++ evaluates small x as exp(nu log(x/2) - lgamma(nu + 1)) times a series, a few ulp).
double series_j(int nu, double x) {
  const double q = 0.25 * x * x;
  double term = nu == 0 ? 1.0 : 0.5 * x;  // k = 0: (x/2)^nu / Gamma(nu + 1)
  double sum = term;
  for (int k = 1; k < 60; ++k) {
    term *= -q / (static_cast<double>(k) * static_cast<double>(k + nu));
    sum += term;
  }
  return sum;
}

}  // namespace

TEST_CASE("bessel_j0 and bessel_j1 at the origin", "[bessel]") {
  REQUIRE(rtt::math::bessel_j0(0.0) == 1.0);
  REQUIRE(rtt::math::bessel_j1(0.0) == 0.0);
}

TEST_CASE("bessel_j0 and bessel_j1 follow the power series (DLMF 10.2.2)", "[bessel]") {
  for (const double x :
       {1e-3, 0.1, 0.5, 1.0, 2.0, 2.404825557695773, 3.0, 3.831705970207512, 5.0, 6.0}) {
    CAPTURE(x);
    REQUIRE(std::abs(rtt::math::bessel_j0(x) - series_j(0, x)) <= 1e-13);
    REQUIRE(std::abs(rtt::math::bessel_j1(x) - series_j(1, x)) <= 1e-13);
  }
}
