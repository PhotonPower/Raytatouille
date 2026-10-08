#pragma once

/// @file ghosts.hpp
/// Ghost generator (#123, ADR 0027): the two-reflection ghost paths of a compiled path, as
/// explicit model paths. A ghost reflects at refracting surface j back and at refracting
/// surface i < j forward again, then continues to the image surface. The ghosts are derived
/// data: the generator never writes them into a system file; compile determines their media
/// with the rules for explicit paths (ADR 0027).

#include <cstddef>
#include <cstdint>
#include <string_view>
#include <vector>

#include "rtt/coating/catalog.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/path.hpp"
#include "rtt/model/system.hpp"

namespace rtt::compile {

/// Options of the ghost generator.
struct GhostOptions {
  /// Largest number of ghost paths; more give std::invalid_argument (never a silent cut).
  /// N refracting events give N (N - 1) / 2 ghosts.
  std::size_t max_paths = 10000;
};

/// Two-reflection ghost paths of the path `base` (ADR 0027):
/// - Ghost surfaces are the Refract events of `base`; Reflect, Transmit and Diffract events are
///   not. In M4 `base` must have no diffraction orders and no crystal modes (Ordinary,
///   Extraordinary): ghosts at gratings and crystals come later (ADR 0027).
/// - The two ghost reflections are specular (order 0).
/// - For each pair of Refract events i < j (event indices in `base`) the ghost is
///   base[0..j-1], Reflect at j, base[i+1..j-1] in reverse order, Reflect at i,
///   base[i+1..end]. Events on the way back keep their kind (and diffraction order).
/// - Order: j ascending, then i ascending. Names "<base> ghost <surface j>/<surface i>"; a
///   surface that occurs more than once in `base` gets the event index appended ("P.S1#3").
/// @param system  compiled system that contains `base`
/// @param base    path to derive the ghosts from, e.g. an automatic path
/// @param options limit of the number of ghosts
/// @return the ghosts as explicit model paths (no automatic flag), ready to append to a copy
///         of the model and to compile
/// @throws std::invalid_argument for an invalid `base`, a `base` with a diffraction order or a
///         crystal mode, more than options.max_paths ghosts, a ghost name that is already the
///         name of a path of `system`, or two ghosts with the same name (surface ids that
///         contain '/' or '#')
[[nodiscard]] std::vector<model::Path> ghost_paths(const CompiledSystem& system,
                                                   PathId base,
                                                   const GhostOptions& options = {});

/// One ghost path in a system compiled with compile_with_ghosts().
struct GhostPath {
  PathId path;                  ///< the ghost path
  PathId base;                  ///< the path it was derived from
  std::uint32_t surface_j = 0;  ///< surface of the first ghost reflection (back), surfaces() index
  std::uint32_t surface_i = 0;  ///< surface of the second reflection (forward), surfaces() index
  std::size_t event_j = 0;      ///< index of the reflecting event j in the base path
  std::size_t event_i = 0;      ///< index of the reflecting event i in the base path
};

/// A system compiled with its ghost paths, and which path is which ghost.
struct GhostSystem {
  CompiledSystem system;          ///< all paths of the model, then the ghosts of `base`
  std::vector<GhostPath> ghosts;  ///< in the order of ghost_paths()
};

/// Compiles `system`, derives the ghosts of the path named `base` with ghost_paths(), appends
/// them to a copy of the model and compiles that copy (ADR 0027). The model itself is not
/// changed.
/// @param system    model to compile; not changed
/// @param base      name of the path to derive the ghosts from
/// @param materials material library, as for compile()
/// @param coatings  coating library, as for compile()
/// @param options   limit of the number of ghosts
/// @return the compiled copy (all paths of the model, then the ghosts) and one GhostPath per
///         ghost
/// @throws CompileError as compile(), std::invalid_argument for an unknown `base` and as
///         ghost_paths()
[[nodiscard]] GhostSystem compile_with_ghosts(const model::System& system,
                                              std::string_view base,
                                              const material::MaterialLibrary& materials,
                                              const coating::CoatingLibrary& coatings,
                                              const GhostOptions& options = {});

}  // namespace rtt::compile
