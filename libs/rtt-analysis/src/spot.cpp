#include "rtt/analysis/spot.hpp"

namespace rtt::analysis {

SpotStatistics spot_statistics(std::span<const SpotPoint> /*points*/, Point2 /*chief*/) {
  return {};  // TODO(issue-28): implement
}

SpotDiagram spot(const compile::CompiledSystem& /*system*/,
                 compile::PathId /*path*/,
                 std::uint16_t /*field*/,
                 std::optional<std::uint16_t> /*wavelength*/,
                 const SpotOptions& /*options*/) {
  return {};  // TODO(issue-28): implement
}

RayFan ray_fan(const compile::CompiledSystem& /*system*/,
               compile::PathId /*path*/,
               std::uint16_t /*field*/,
               std::uint16_t /*wavelength*/,
               const FanOptions& /*options*/) {
  return {};  // TODO(issue-28): implement
}

}  // namespace rtt::analysis
