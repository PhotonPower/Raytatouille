#pragma once

/// @file bessel.hpp
/// Bessel functions of the first kind J0 and J1 on every supported standard library (#373).
///
/// J_nu(x) = (x/2)^nu sum_k (-1)^k (x^2/4)^k / (k! Gamma(nu + k + 1)) (NIST DLMF, Version 1.2.8,
/// Eq. (10.2.2)). No own approximation: the values come from the platform.
/// - With the C++17 mathematical special functions (`__cpp_lib_math_special_functions`, set
///   by libstdc++ and MSVC): `std::cyl_bessel_j`.
/// - Otherwise (libc++, e.g. Apple clang on macOS, which does not implement them): `j0` and
///   `j1` of the C library (POSIX, in the libm of macOS and glibc).
///
/// Like `sin` or `cos`, the last bits may differ between platforms; bit equality is guaranteed
/// per platform only (ADR 0011).

#include <cmath>
#include <version>

#if defined(__cpp_lib_math_special_functions)
// std::cyl_bessel_j from <cmath>
#elif defined(_LIBCPP_VERSION)
#include <math.h>  // NOLINT(modernize-deprecated-headers): POSIX j0, j1 are not in <cmath>
#else
#error "rtt/math/bessel.hpp: neither std::cyl_bessel_j nor libc++ with POSIX j0/j1"
#endif

namespace rtt::math {

/// Bessel function of the first kind of order 0, J0(x) (DLMF (10.2.2) with nu = 0).
/// @param x argument, dimensionless, x >= 0 (the domain of `std::cyl_bessel_j`)
/// @return J0(x), dimensionless
[[nodiscard]] inline double bessel_j0(double x) {
#if defined(__cpp_lib_math_special_functions)
  return std::cyl_bessel_j(0.0, x);
#else
  return ::j0(x);
#endif
}

/// Bessel function of the first kind of order 1, J1(x) (DLMF (10.2.2) with nu = 1).
/// @param x argument, dimensionless, x >= 0 (the domain of `std::cyl_bessel_j`)
/// @return J1(x), dimensionless
[[nodiscard]] inline double bessel_j1(double x) {
#if defined(__cpp_lib_math_special_functions)
  return std::cyl_bessel_j(1.0, x);
#else
  return ::j1(x);
#endif
}

}  // namespace rtt::math
