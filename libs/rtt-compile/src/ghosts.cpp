#include "rtt/compile/ghosts.hpp"

// STUB for the red run of the #123 tests: no ghosts.

namespace rtt::compile {

std::vector<model::Path> ghost_paths(const CompiledSystem& /*system*/,
                                     PathId /*base*/,
                                     const GhostOptions& /*options*/) {
  return {};
}

GhostSystem compile_with_ghosts(const model::System& system,
                                std::string_view /*base*/,
                                const material::MaterialLibrary& materials,
                                const coating::CoatingLibrary& coatings,
                                const GhostOptions& /*options*/) {
  return {compile(system, materials, coatings), {}};
}

}  // namespace rtt::compile
