#pragma once

/// @file element.hpp
/// Hierarchy: Assembly -> Element -> Surface (docs/architecture.md, "Hierarchie").

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <variant>
#include <vector>

#include "rtt/model/pose.hpp"
#include "rtt/model/surface.hpp"

namespace rtt::model {

enum class ElementKind : std::uint8_t {
  Lens,         ///< >= 2 surfaces, one material per segment (cemented groups allowed, ADR 0017)
  Mirror,       ///< >= 1 surface, substrate material optional (shorthand only)
  Plate,        ///< >= 2 plane surfaces, one material per segment (plates, cubes, prisms)
  ThinElement,  ///< exactly 1 surface without thickness (ideal polarizer, retarder, splitter)
  Stop,         ///< exactly 1 surface, aperture stop
  Detector,     ///< exactly 1 surface
};

/// Uniaxial crystal (ADR 0026): two ordinary catalog references for the principal indices n_O
/// and n_E (e.g. the Zemax convention X / X-E). File form: "material": {"ordinary": ...,
/// "extraordinary": ...}; Lens and Plate only, one crystal for all segments.
struct CrystalMaterial {
  std::string ordinary;       ///< catalog reference for n_O, e.g. "BIREFRINGENT:CALCITE"
  std::string extraordinary;  ///< catalog reference for n_E, e.g. "BIREFRINGENT:CALCITE-E"
  bool operator==(const CrystalMaterial&) const = default;
};

/// A physical body. Outside it the medium is the surrounding medium.
///
/// Inside, an element with N surfaces has N - 1 segments; segment i is the glass between
/// surface i and surface i + 1 (both zero-based, in the order of `surfaces`). The material of
/// the segments is given in exactly one of two forms (ADR 0017):
/// - `material`: shorthand, one catalog reference for all segments (file: a string);
/// - `segment_materials`: one catalog reference per segment, N - 1 entries (file: an array).
/// Only Lens and Plate accept `segment_materials`. Both forms are kept as read, so that writing
/// reproduces the file byte for byte. A crystal (ADR 0026) is given in `crystal` instead; then
/// `material` and `segment_materials` are empty.
struct Element {
  std::string name;
  ElementKind kind = ElementKind::Lens;
  Pose pose;
  std::optional<std::string> material;  ///< catalog reference for all segments, "SCHOTT:N-BK7"
  std::vector<Surface> surfaces;
  /// Catalog reference per segment (segment i between surfaces i and i + 1); empty if the
  /// shorthand `material` is used or the element has no material. Last member with an explicit
  /// initializer, so that aggregate initialisations {name, kind, pose, material, surfaces} stay
  /// valid and free of -Wmissing-field-initializers.
  std::vector<std::string> segment_materials{};  // NOLINT(readability-redundant-member-init)
  /// Uniaxial crystal for all segments (ADR 0026); Lens and Plate only.
  std::optional<CrystalMaterial> crystal{};  // NOLINT(readability-redundant-member-init)
  /// Direction of the optic axis of the crystal in element coordinates, dimensionless, not
  /// normalised; the sign has no meaning (ADR 0026, point 2). Required exactly for a crystal.
  std::optional<std::array<double, 3>> optic_axis{};  // NOLINT(readability-redundant-member-init)
  bool operator==(const Element&) const = default;

  /// Material of segment `segment` (between surfaces `segment` and `segment + 1`): the list
  /// entry if `segment_materials` is used, otherwise the shorthand `material`. Empty if the
  /// element has no isotropic material (none, or a crystal) or `segment + 1` is not a valid
  /// surface index.
  [[nodiscard]] std::optional<std::string> segment_material(std::size_t segment) const;
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
