#pragma once

/// @file element.hpp
/// Hierarchy: Assembly -> Element -> Surface (docs/architecture.md, "Hierarchie").

#include <cstdint>
#include <optional>
#include <string>
#include <variant>
#include <vector>

#include "rtt/model/pose.hpp"
#include "rtt/model/surface.hpp"

namespace rtt::model {

enum class ElementKind : std::uint8_t {
  Lens,         ///< >= 2 surfaces, material required (cemented groups allowed)
  Mirror,       ///< >= 1 surface, substrate material optional
  Plate,        ///< >= 2 plane surfaces, material required (plates, cubes, prisms)
  ThinElement,  ///< exactly 1 surface without thickness (ideal polarizer, retarder, splitter)
  Stop,         ///< exactly 1 surface, aperture stop
  Detector,     ///< exactly 1 surface
};

/// A physical body. Inside it the medium is `material`, outside the surrounding medium.
struct Element {
  std::string name;
  ElementKind kind = ElementKind::Lens;
  Pose pose;
  std::optional<std::string> material;  ///< catalog reference, e.g. "SCHOTT:N-BK7"
  std::vector<Surface> surfaces;
  bool operator==(const Element&) const = default;
};

struct Node;

/// Group of nodes moved and toleranced as one rigid body.
struct Assembly {
  std::string name;
  Pose pose;
  std::vector<Node> children;
  bool operator==(const Assembly&) const;
};

/// Tree node: either an assembly or an element.
struct Node {
  std::variant<Assembly, Element> value;
  bool operator==(const Node&) const;
};

// Both comparisons are defined out of line in element.cpp: defaulting them on the
// recursive Assembly <-> Node pair crashes Clang 18.

}  // namespace rtt::model
