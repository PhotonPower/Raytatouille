// Principal square root of a complex number (rtt/math/complex.hpp, ADR 0046).
// Reference: M. Abramowitz, I. A. Stegun, Handbook of Mathematical Functions (1964), 3.7.3 and
// 3.7.27: z^(1/2) = [(r + x)/2]^(1/2) +- i [(r - x)/2]^(1/2) = u +- i v with r = |z|, 2 u v = y,
// the sign that of y; the principal root for -pi < arg z <= pi.

#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <complex>
#include <limits>

#include "rtt/math/complex.hpp"

using rtt::math::principal_sqrt;
using Complex = std::complex<double>;

TEST_CASE("principal_sqrt of a negative real is purely imaginary, bit for bit", "[complex]") {
  // r = |x|, so u = |y| / (2 v) = 0 and v = sqrt((r - x)/2) = sqrt(|x|): no rounding in u. A
  // square root through the polar form gives u = sqrt(|x|) cos(pi/2) ~ 6e-17 sqrt(|x|) instead
  // (libc++), which hides total internal reflection (#373).
  for (const double x : {-0.5, -2.25, -1e-300, -3.0e200}) {
    const Complex q = principal_sqrt(Complex(x, 0.0));
    REQUIRE(q.real() == 0.0);
    REQUIRE_FALSE(std::signbit(q.real()));
    REQUIRE(q.imag() == std::sqrt(-x));
  }
}

TEST_CASE("principal_sqrt of a non-negative real is real, bit for bit", "[complex]") {
  for (const double x : {0.25, 2.25, 1.0, 1e-300, 3.0e200}) {
    const Complex q = principal_sqrt(Complex(x, 0.0));
    REQUIRE(q.real() == std::sqrt(x));
    REQUIRE(q.imag() == 0.0);
    REQUIRE_FALSE(std::signbit(q.imag()));
  }
  const Complex zero = principal_sqrt(Complex(0.0, 0.0));
  REQUIRE(zero.real() == 0.0);
  REQUIRE(zero.imag() == 0.0);
}

TEST_CASE("principal_sqrt of Pythagorean squares is exact", "[complex]") {
  // (3 + 4i)^2 = -7 + 24i, (4 + 3i)^2 = 7 + 24i, (1 + i)^2 = 2i: every step of 3.7.27 is exact,
  // provided that hypot is exact here. C11 Annex F does not require a correctly rounded hypot,
  // so the precondition is checked at run time (volatile: not folded at compile time).
  volatile double seven = 7.0;
  volatile double twenty_four = 24.0;
  volatile double two = 2.0;
  REQUIRE(std::hypot(-seven, twenty_four) == 25.0);
  REQUIRE(std::hypot(0.0, two) == 2.0);
  REQUIRE(principal_sqrt(Complex(-7.0, 24.0)) == Complex(3.0, 4.0));
  REQUIRE(principal_sqrt(Complex(7.0, 24.0)) == Complex(4.0, 3.0));
  REQUIRE(principal_sqrt(Complex(0.0, 2.0)) == Complex(1.0, 1.0));
  // The sign of v follows y (3.7.27): lower half plane.
  REQUIRE(principal_sqrt(Complex(-7.0, -24.0)) == Complex(3.0, -4.0));
  REQUIRE(principal_sqrt(Complex(7.0, -24.0)) == Complex(4.0, -3.0));
  REQUIRE(principal_sqrt(Complex(0.0, -2.0)) == Complex(1.0, -1.0));
}

TEST_CASE("principal_sqrt takes the sign of a zero imaginary part like C99 csqrt", "[complex]") {
  // On the cut, y = -0 is the lower side: sqrt(x - i0) = -i sqrt(|x|) for x < 0 (C11 G.6.4.2:
  // csqrt(conj(z)) = conj(csqrt(z))).
  const Complex lower = principal_sqrt(Complex(-4.0, -0.0));
  REQUIRE(lower.real() == 0.0);
  REQUIRE(lower.imag() == -2.0);
  const Complex upper = principal_sqrt(Complex(-4.0, 0.0));
  REQUIRE(upper.imag() == 2.0);
  const Complex positive = principal_sqrt(Complex(4.0, -0.0));
  REQUIRE(positive.real() == 2.0);
  REQUIRE(positive.imag() == 0.0);
  REQUIRE(std::signbit(positive.imag()));
}

TEST_CASE("principal_sqrt squares back and lies in the right half plane", "[complex]") {
  constexpr double kEps = std::numeric_limits<double>::epsilon();
  for (const double x : {-3.7, -1.0, -1e-3, 0.0, 1e-3, 0.8, 5.0}) {
    for (const double y : {-2.5, -1e-8, 1e-8, 0.3, 4.0}) {
      const Complex z(x, y);
      const Complex q = principal_sqrt(z);
      REQUIRE(q.real() >= 0.0);
      REQUIRE(std::signbit(q.imag()) == std::signbit(y));
      // u and v carry a few ulp each (one sqrt, one division); q^2 then deviates by a few ulp
      // of |z|.
      REQUIRE(std::abs(q * q - z) <= 8.0 * kEps * std::abs(z));
    }
  }
}
