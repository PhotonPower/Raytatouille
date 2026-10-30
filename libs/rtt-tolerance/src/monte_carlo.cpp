#include "rtt/tolerance/monte_carlo.hpp"

#include <cstddef>
#include <span>

#include "monte_carlo_detail.hpp"

// Stub for the red run of the tests (#259); the implementation follows in the next commit.

namespace rtt::tolerance {

namespace detail {

std::size_t nearest_rank(double /*p*/, std::size_t /*n*/) noexcept {
  return 0;
}

Yield wilson(std::size_t /*good*/, std::size_t /*total*/, double /*z*/) noexcept {
  return {};
}

Moments moments(std::span<const double> /*x*/) noexcept {
  return {};
}

}  // namespace detail

MonteCarloResult monte_carlo(const TrialContext& /*context*/,
                             const MonteCarloOptions& /*options*/,
                             const trace::RunControl* /*control*/) {
  return {};
}

MonteCarloResult monte_carlo(const model::System& /*system*/,
                             const material::MaterialLibrary& /*materials*/,
                             const coating::CoatingLibrary* /*coatings*/,
                             const MonteCarloOptions& /*options*/,
                             const trace::RunControl* /*control*/) {
  return {};
}

}  // namespace rtt::tolerance
