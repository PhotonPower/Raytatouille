#include "rtt/tolerance/trial.hpp"

#include <memory>
#include <optional>
#include <span>
#include <utility>
#include <vector>

// Stub (#257): the interface only, for the red run of the tests.

namespace rtt::tolerance {

namespace detail {

struct ContextData {
  model::System nominal;
  std::vector<model::ToleranceItem> items;
  TrialResult nominal_result;
};

struct ContextAccess {
  static const ContextData& data(const TrialContext& context) { return *context.data_; }
};

}  // namespace detail

std::span<const model::ToleranceItem> TrialContext::items() const noexcept {
  return data_->items;
}

std::span<const model::Criterion> TrialContext::criteria() const noexcept {
  return data_->nominal.tolerances.criteria;
}

std::span<const model::Compensator> TrialContext::compensators() const noexcept {
  return data_->nominal.tolerances.compensators;
}

const TrialResult& TrialContext::nominal() const noexcept {
  return data_->nominal_result;
}

TrialContext prepare(const model::System& nominal,
                     const material::MaterialLibrary& /*materials*/,
                     const coating::CoatingLibrary* /*coatings*/,
                     std::span<const model::ToleranceItem> items,
                     const TrialOptions& /*options*/,
                     const trace::RunControl* /*control*/) {
  auto d = std::make_shared<detail::ContextData>();
  d->nominal = nominal;
  d->items.assign(items.begin(), items.end());
  return TrialContext(std::move(d));
}

TrialResult evaluate(const TrialContext& /*context*/,
                     std::span<const double> /*values*/,
                     const trace::RunControl* /*control*/) {
  return {};
}

std::optional<bool> satisfies(const model::Criterion& /*criterion*/,
                              std::optional<double> /*value*/,
                              std::optional<double> /*nominal_value*/) {
  return std::nullopt;
}

Applied apply(const TrialContext& context, std::span<const double> /*values*/) {
  return {detail::ContextAccess::data(context).nominal, {}};
}

}  // namespace rtt::tolerance
