// Stub for the red run (#258): only the interface; the implementation follows in the next commit.
#include "rtt/tolerance/sensitivity.hpp"

namespace rtt::tolerance {

CriterionKind criterion_kind(const model::Operand& /*operand*/) noexcept {
  return CriterionKind::Linear;
}

double sigma_ratio(model::ToleranceDistribution /*distribution*/) {
  return 0.0;
}

ToleranceSensitivity sensitivity(const TrialContext& /*context*/,
                                 const SensitivityOptions& /*options*/,
                                 const trace::RunControl* /*control*/) {
  return {};
}

ToleranceSensitivity sensitivity(const model::System& /*system*/,
                                 const material::MaterialLibrary& /*materials*/,
                                 const coating::CoatingLibrary* /*coatings*/,
                                 const SensitivityOptions& /*options*/,
                                 const trace::RunControl* /*control*/) {
  return {};
}

}  // namespace rtt::tolerance
