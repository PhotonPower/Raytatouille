#pragma once

/// @file complex.hpp
/// Complex functions whose result must not depend on the standard library (ADR 0046).

#include <cmath>
#include <complex>

#include "rtt/math/real.hpp"

namespace rtt::math {

/// Principal square root of z = x + i y, computed in Cartesian form (M. Abramowitz,
/// I. A. Stegun, Handbook of Mathematical Functions, 1964, 3.7.3 and 3.7.27):
/// sqrt(z) = u + i v with r = |z| = (x^2 + y^2)^(1/2), u = [(r + x)/2]^(1/2) >= 0,
/// |v| = [(r - x)/2]^(1/2), 2 u v = y, and v carrying the sign of y; the principal root for
/// -pi < arg z <= pi.
///
/// Only the component without cancellation comes from a square root; the other one from
/// 2 u v = y. A negative real x (y = +0) therefore gives exactly u = +0 and v = sqrt(|x|);
/// `std::sqrt(std::complex)` of libc++ goes through the polar form and returns
/// u = sqrt(|x|) cos(pi/2) != 0 there, which hides total internal reflection (#373).
///
/// Signed zeros follow C99/C11 `csqrt` (C11 Annex G.6.4.2): sqrt(conj z) = conj sqrt(z), so
/// y = -0 is the lower side of the cut (sqrt(x - i0) = -i sqrt(|x|) for x < 0), and
/// sqrt(+-0 + i y) for y = +-0 is +0 + i y.
///
/// Accuracy: a few ulp in each component (one square root, one division), within the
/// precondition below. Outside it the result may be wrong: subnormal x or y can give NaN or
/// infinity (r / 2 underflows to 0), and |z| > the largest finite T overflows in `hypot`. The
/// physics callers stay far inside (|n^2 - xi^2| is 0 or of order 1).
/// @param z complex argument, any unit U
/// @return principal square root u + i v with u >= 0, unit sqrt(U)
/// @pre x and y each zero or normal (|.| >= std::numeric_limits<T>::min()), and
///      |z| <= std::numeric_limits<T>::max()
template <Real T>
[[nodiscard]] std::complex<T> principal_sqrt(std::complex<T> z) noexcept {
  using std::abs;
  using std::copysign;
  using std::hypot;
  using std::sqrt;
  const T x = z.real();
  const T y = z.imag();
  if (x == T(0) && y == T(0)) return {T(0), y};
  const T r = hypot(x, y);  // 3.7.3, without overflow of x^2 + y^2
  // Halves taken first, so that r + |x| does not overflow for |z| up to the largest finite T.
  if (x >= T(0)) {
    const T u = sqrt(r / T(2) + x / T(2));
    return {u, y / (T(2) * u)};
  }
  const T v = sqrt(r / T(2) - x / T(2));
  return {abs(y) / (T(2) * v), copysign(v, y)};
}

}  // namespace rtt::math
