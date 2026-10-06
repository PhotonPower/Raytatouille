#pragma once

/// @file compiled_element.hpp
/// Elements of a compiled system (#81): which surfaces belong to a body and which media fill
/// it, for drawing lens sections (rtt/compile/layout.hpp).

#include <cstdint>
#include <string>
#include <vector>

#include "rtt/model/element.hpp"

namespace rtt::compile {

/// One element of the model in tree order (CompiledSystem::elements()).
struct CompiledElement {
  std::string name;  ///< element name from the model
  model::ElementKind kind = model::ElementKind::Lens;
  std::uint32_t first_surface = 0;  ///< surfaces [first_surface, + surface_count)
  std::uint32_t surface_count = 0;
  /// Media of the body, indices into CompiledSystem::media(): Lens and Plate one per segment
  /// (segment i between surfaces i and i + 1, ADR 0017), Mirror the substrate as one body;
  /// empty if the element has no material.
  std::vector<std::uint32_t> media;
  /// True if the media follow the segments (every Lens, and a Plate with different segment
  /// materials); otherwise a refraction toggles between inside and environment (ADR 0017).
  bool segmented = false;
};

}  // namespace rtt::compile
