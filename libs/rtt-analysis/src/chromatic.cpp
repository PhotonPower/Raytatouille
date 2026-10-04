#include "rtt/analysis/chromatic.hpp"

namespace rtt::analysis {

LongitudinalColour longitudinal_colour(const compile::CompiledSystem& /*system*/,
                                       compile::PathId /*path*/,
                                       const ChromaticOptions& /*options*/) {
  return {};  // TODO(issue-31): implement
}

LateralColour lateral_colour(const compile::CompiledSystem& /*system*/,
                             compile::PathId /*path*/,
                             std::uint16_t /*field*/,
                             trace::Aiming /*aiming*/) {
  return {};  // TODO(issue-31): implement
}

}  // namespace rtt::analysis
