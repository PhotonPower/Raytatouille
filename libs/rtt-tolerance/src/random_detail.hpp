#pragma once

/// @file random_detail.hpp
/// Kernels of the normal distribution of random.hpp (ADR 0036, point 5), separate so that the
/// tests can feed them chosen numbers. Not part of the API.

#include <cmath>
#include <numbers>
#include <utility>

namespace rtt::tolerance::detail {

/// Box-Muller pair (E. W. Weisstein, "Box-Muller Transformation", MathWorld, Eqs. (1) and (2))
/// from two uniform numbers in [0, 1): x1 = 1 - u_a in (0, 1] (the logarithm needs x1 > 0) and
/// x2 = u_b; z1 = sqrt(-2 ln x1) cos(2 pi x2), z2 = sqrt(-2 ln x1) sin(2 pi x2).
[[nodiscard]] inline std::pair<double, double> box_muller(double u_a, double u_b) {
  const double r = std::sqrt(-2.0 * std::log(1.0 - u_a));
  const double angle = 2.0 * std::numbers::pi * u_b;
  return {r * std::cos(angle), r * std::sin(angle)};
}

/// Truncated normal value with sigma = limit / 2 by rejection (ADR 0036, point 5): `normal(k)`
/// gives draw z_k of a standard normal distribution; the first sigma z_k with |z_k| <= 2,
/// k < max_draws, is returned. Otherwise limit (2 u - 1) with u = fallback() in [0, 1).
/// `normal` is called with k = 0, 1, 2, ... in ascending order, each k at most once, and stops
/// at the first accepted draw; random.cpp relies on this (it computes the pair of a call at its
/// even k and keeps it for the odd k).
template <class Normal, class Fallback>
[[nodiscard]] double truncated_normal(double limit,
                                      int max_draws,
                                      Normal&& normal,
                                      Fallback&& fallback) {
  const double sigma = limit / 2.0;
  for (int k = 0; k < max_draws; ++k) {
    const double z = normal(k);
    if (std::abs(z) <= 2.0) return sigma * z;
  }
  return limit * (2.0 * fallback() - 1.0);
}

}  // namespace rtt::tolerance::detail
