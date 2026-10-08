#include "rtt/compile/ghosts.hpp"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace rtt::compile {
namespace {

/// One ghost before naming: the indices of the reflecting events in the base path.
struct GhostPair {
  std::size_t j = 0;
  std::size_t i = 0;
};

/// Name part of event k: the surface id, with "#k" if the surface occurs more than once.
std::string surface_label(const CompiledSystem& system, const CompiledPath& base, std::size_t k) {
  const std::uint32_t surface = base.events[k].surface;
  std::size_t count = 0;
  for (const auto& e : base.events) count += e.surface == surface ? 1 : 0;
  std::string label = system.surfaces()[surface].id.str();
  if (count > 1) label += "#" + std::to_string(k);
  return label;
}

model::Event model_event(const CompiledSystem& system,
                         const CompiledEvent& e,
                         model::EventKind kind) {
  return model::Event{system.surfaces()[e.surface].id, kind, e.order};
}

/// Pairs (j, i) of Refract events, j ascending, then i ascending (ADR 0027).
std::vector<GhostPair> ghost_pairs(const CompiledPath& base, std::size_t max_paths) {
  std::vector<std::size_t> refract;
  for (std::size_t k = 0; k < base.events.size(); ++k) {
    if (base.events[k].kind == model::EventKind::Refract) refract.push_back(k);
  }
  const std::size_t n = refract.size();
  const std::size_t count = n < 2 ? 0 : n * (n - 1) / 2;
  if (count > max_paths) {
    throw std::invalid_argument("ghosts: path '" + base.name + "' has " + std::to_string(count) +
                                " ghosts, more than max_paths = " + std::to_string(max_paths));
  }
  std::vector<GhostPair> pairs;
  pairs.reserve(count);
  for (std::size_t b = 1; b < n; ++b) {
    for (std::size_t a = 0; a < b; ++a) pairs.push_back({refract[b], refract[a]});
  }
  return pairs;
}

}  // namespace

std::vector<model::Path> ghost_paths(const CompiledSystem& system,
                                     PathId base,
                                     const GhostOptions& options) {
  if (base.index >= system.paths().size()) {
    throw std::invalid_argument("ghosts: path index " + std::to_string(base.index) +
                                " does not exist");
  }
  const CompiledPath& path = system.path(base);
  const std::vector<GhostPair> pairs = ghost_pairs(path, options.max_paths);
  std::vector<model::Path> out;
  out.reserve(pairs.size());
  for (const GhostPair& p : pairs) {
    model::Path ghost;
    ghost.name = path.name + " ghost " + surface_label(system, path, p.j) + "/" +
                 surface_label(system, path, p.i);
    if (system.find_path(ghost.name)) {
      throw std::invalid_argument("ghosts: the system already has a path named '" + ghost.name +
                                  "'");
    }
    // base[0..j-1], Reflect at j, base[i+1..j-1] reversed, Reflect at i, base[i+1..end].
    for (std::size_t k = 0; k < p.j; ++k) {
      ghost.events.push_back(model_event(system, path.events[k], path.events[k].kind));
    }
    ghost.events.push_back(model_event(system, path.events[p.j], model::EventKind::Reflect));
    for (std::size_t k = p.j - 1; k > p.i; --k) {
      ghost.events.push_back(model_event(system, path.events[k], path.events[k].kind));
    }
    ghost.events.push_back(model_event(system, path.events[p.i], model::EventKind::Reflect));
    for (std::size_t k = p.i + 1; k < path.events.size(); ++k) {
      ghost.events.push_back(model_event(system, path.events[k], path.events[k].kind));
    }
    out.push_back(std::move(ghost));
  }
  return out;
}

GhostSystem compile_with_ghosts(const model::System& system,
                                std::string_view base,
                                const material::MaterialLibrary& materials,
                                const coating::CoatingLibrary& coatings,
                                const GhostOptions& options) {
  const CompiledSystem plain = compile(system, materials, coatings);
  const std::optional<PathId> base_id = plain.find_path(base);
  if (!base_id) {
    throw std::invalid_argument("ghosts: the system has no path named '" + std::string(base) + "'");
  }
  const std::vector<model::Path> ghosts = ghost_paths(plain, *base_id, options);
  const std::vector<GhostPair> pairs = ghost_pairs(plain.path(*base_id), options.max_paths);
  model::System copy = system;
  for (const model::Path& g : ghosts) copy.paths.push_back(g);
  GhostSystem out{compile(copy, materials, coatings), {}};
  out.ghosts.reserve(ghosts.size());
  const CompiledPath& base_path = out.system.path(*base_id);
  for (std::size_t k = 0; k < ghosts.size(); ++k) {
    // `copy` holds every ghost name, so compile kept it; anything else is a broken invariant.
    const std::optional<PathId> id = out.system.find_path(ghosts[k].name);
    if (!id) throw std::logic_error("ghosts: compiled ghost path not found: " + ghosts[k].name);
    GhostPath info;
    info.path = *id;
    info.base = *base_id;
    info.surface_j = base_path.events[pairs[k].j].surface;
    info.surface_i = base_path.events[pairs[k].i].surface;
    info.event_j = pairs[k].j;
    info.event_i = pairs[k].i;
    out.ghosts.push_back(info);
  }
  return out;
}

}  // namespace rtt::compile
