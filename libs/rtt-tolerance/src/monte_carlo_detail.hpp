#pragma once

/// @file monte_carlo_detail.hpp
/// Internal parts of monte_carlo.cpp, tested directly (tests/test_monte_carlo.cpp): the nearest
/// rank, the Wilson interval and the two-pass moments (ADR 0036, point 7).

#include <cstddef>
#include <optional>
#include <span>

#include "rtt/tolerance/monte_carlo.hpp"

namespace rtt::tolerance::detail {

/// r = ceil(p n) with the integer catch of ADR 0036, point 7: x = p n in double; if x lies at most
/// 4 eps x from the nearest integer k (eps = 2^-52), r = k, else r = ceil(x); then clamped to
/// [1, n]. 1-based. Requires p in (0, 1] and n >= 1.
[[nodiscard]] std::size_t nearest_rank(double p, std::size_t n) noexcept;

/// The yield good / total with the Wilson interval for z > 0 (NIST e-Handbook, Sec. 7.2.4.1;
/// ADR 0036, point 7): lower exactly 0 for good = 0, upper exactly 1 for good = total, both
/// clamped to [0, 1]. Without values for total = 0.
[[nodiscard]] Yield wilson(std::size_t good, std::size_t total, double z) noexcept;

/// Mean and standard deviation (n - 1) in two passes in the order of `x`: first the mean, then
/// the sum of the squared deviations (ADR 0036, point 7; not sum x^2 - n mean^2).
struct Moments {
  std::optional<double> mean;                ///< from n = 1
  std::optional<double> standard_deviation;  ///< from n = 2
};
[[nodiscard]] Moments moments(std::span<const double> x) noexcept;

}  // namespace rtt::tolerance::detail
