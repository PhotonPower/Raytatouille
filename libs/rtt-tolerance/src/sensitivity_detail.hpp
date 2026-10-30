#pragma once

/// @file sensitivity_detail.hpp
/// Internal parts of sensitivity.cpp, tested directly (tests/test_sensitivity.cpp): the quadratic
/// reductions, the warning tolerances.rss_negative and the distribution of the allowance.

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <string_view>

#include "rtt/model/tolerance.hpp"
#include "rtt/model/validate.hpp"
#include "rtt/tolerance/sensitivity.hpp"

namespace rtt::tolerance::detail {

/// A sum of squares that is reported as 0 when it is negative (ADR 0036, point 6).
struct QuadraticRss {
  double value = 0.0;     ///< sqrt(max(0, sum)) or max(0, sum), see the function
  bool negative = false;  ///< the sum is negative: reported as 0 with a warning
};

/// Phi_RSS = sqrt(max(0, Phi_0^2 + sum q_i)), the sum in the order of `q` (ADR 0036, point 8).
[[nodiscard]] inline QuadraticRss quadratic_rss(double nominal, std::span<const double> q) {
  double sum = nominal * nominal;
  for (const double x : q) sum += x;
  if (sum < 0.0) return {0.0, true};
  return {std::sqrt(sum), false};
}

/// E[Phi^2]_RSS = max(0, Phi_0^2 + sum c_i^2 q_i) (ADR 0036, point 6), the sum in the order of
/// `q`; `ratios` holds c_i of the same items. Unit of the criterion squared. It can be negative
/// although Phi_0^2 + sum q_i is not, since c_i <= 1 weights the items differently.
[[nodiscard]] inline QuadraticRss predicted_mean_square(double nominal,
                                                        std::span<const double> q,
                                                        std::span<const double> ratios) {
  double sum = nominal * nominal;
  for (std::size_t i = 0; i < q.size(); ++i) sum += ratios[i] * ratios[i] * q[i];
  if (sum < 0.0) return {0.0, true};
  return {sum, false};
}

/// Which quadratic sum of tolerances.rss_negative is negative.
enum class NegativeSum : std::uint8_t {
  Rss,         ///< Phi_0^2 + sum q_i: rss reported as 0
  Prediction,  ///< Phi_0^2 + sum c_i^2 q_i: predicted_mean_square reported as 0
};

/// The warning tolerances.rss_negative of criterion k at /tolerances/criteria/<k>; the message
/// names the criterion, the part of the run (`stage`: "sensitivity" or "confirmation") and the
/// sum.
[[nodiscard]] model::Diagnostic rss_warning(const model::Criterion& criterion,
                                            std::size_t k,
                                            std::string_view stage,
                                            NegativeSum sum);

/// A distributed allowance: A in the unit of the criterion and the share of one item (unit of
/// the criterion, or its square for a quadratic criterion).
struct Distribution {
  double allowance = 0.0;
  double share = 0.0;
};

/// The allowance A of a criterion and the share of one of `n` items (ADR 0036, point 6). Sets
/// allowance, note and share of `inverse`; returns them if a positive A is distributed. Notes:
/// "no items", "no bound" (linear), "no upper bound" (quadratic), "nominal value on or outside
/// a bound" (A = 0).
[[nodiscard]] std::optional<Distribution> distribute(const model::Criterion& criterion,
                                                     CriterionKind kind,
                                                     double phi0,
                                                     std::size_t n,
                                                     InverseSensitivity& inverse);

}  // namespace rtt::tolerance::detail
