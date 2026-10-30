#pragma once

/// @file sensitivity_detail.hpp
/// Internal reductions of sensitivity.cpp, tested directly (tests/test_sensitivity.cpp).

#include <cmath>
#include <span>

namespace rtt::tolerance::detail {

/// Phi_RSS of a quadratic criterion (ADR 0036, point 6; Burge, slides 5 and 9).
struct QuadraticRss {
  double value = 0.0;     ///< sqrt(max(0, Phi_0^2 + sum q_i)), unit of the criterion
  bool negative = false;  ///< Phi_0^2 + sum q_i < 0: reported as 0 with a warning
};

/// Phi_RSS = sqrt(max(0, Phi_0^2 + sum q_i)), the sum in the order of `q` (ADR 0036, point 8).
[[nodiscard]] inline QuadraticRss quadratic_rss(double nominal, std::span<const double> q) {
  double sum = nominal * nominal;
  for (const double x : q) sum += x;
  if (sum < 0.0) return {0.0, true};
  return {std::sqrt(sum), false};
}

}  // namespace rtt::tolerance::detail
