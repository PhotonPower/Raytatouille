#include "rtt/analysis/ghosts.hpp"

// STUB for the red run of the #124 tests: an empty ranking.

namespace rtt::analysis {

GhostRanking ghost_ranking(const compile::GhostSystem& ghosts,
                           std::uint16_t field,
                           std::uint16_t wavelength,
                           const GhostRankingOptions& options) {
  GhostRanking out;
  out.base = ghosts.ghosts.empty() ? compile::PathId{0} : ghosts.ghosts.front().base;
  out.field = field;
  out.wavelength = wavelength;
  out.resolution_radius = options.resolution_radius;
  return out;
}

GhostRanking ghost_ranking(const compile::GhostSystem& ghosts,
                           std::uint16_t field,
                           std::uint16_t wavelength,
                           const GhostRankingOptions& options,
                           const trace::RunControl& /*control*/) {
  return ghost_ranking(ghosts, field, wavelength, options);
}

}  // namespace rtt::analysis
