#pragma once

/// @file tolerance_targets.hpp
/// Private: the anchors of a system (ADR 0035, points 2 and 3) with their JSON pointers into the
/// edit form, the Params by pointer, and the perturbed quantity of a tolerance. Shared by the
/// class enumeration and validate.

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include "rtt/model/system.hpp"
#include "rtt/model/tolerance.hpp"

namespace rtt::model::detail {

/// An element or an assembly below the root.
struct NodeEntry {
  std::string pointer;                 ///< e.g. "/root/children/1"
  const Assembly* assembly = nullptr;  ///< set for an assembly
  const Element* element = nullptr;    ///< set for an element
  /// Names of the enclosing assemblies below the root, outermost first.
  std::vector<std::string> ancestors;
};

struct SurfaceEntry {
  std::string pointer;  ///< e.g. "/root/children/1/surfaces/0"
  const Surface* surface = nullptr;
  const Element* element = nullptr;
  std::size_t index = 0;  ///< position in the element's surface list
};

/// The anchors of a system. With duplicate names (node.name_duplicate, surface.id_duplicate,
/// parameters.name_duplicate) the first in tree or table order counts.
class TargetIndex {
 public:
  explicit TargetIndex(const System& system);

  [[nodiscard]] const NodeEntry* node(const std::string& name) const;
  [[nodiscard]] const SurfaceEntry* surface(const std::string& id) const;
  [[nodiscard]] const ParameterRow* row(const std::string& name) const;
  /// The Param at an absolute pointer into the edit form, or null.
  [[nodiscard]] const Param* param(const std::string& pointer) const;
  /// True if `name` is the root assembly.
  [[nodiscard]] bool is_root(const std::string& name) const { return name == root_name_; }
  /// The pointer of the anchor's object (surface or node), or none for a row or an unknown or
  /// root anchor.
  [[nodiscard]] std::optional<std::string> anchor_pointer(const ToleranceTarget& target) const;

 private:
  void walk(const Assembly& a, const std::string& pointer, std::vector<std::string>& ancestors);

  std::string root_name_;
  std::unordered_map<std::string, NodeEntry> nodes_;
  std::unordered_map<std::string, SurfaceEntry> surfaces_;
  std::unordered_map<std::string, const ParameterRow*> rows_;
  std::unordered_map<std::string, const Param*> params_;
};

/// The quantity a tolerance perturbs (ADR 0035, point 2) as a key: the absolute pointer of the
/// Param for the dimensional kinds (radius and curvature_relative share the radius), "row:" plus
/// the name for a row, the anchor pointer plus kind for the rigid kinds, and the element pointer
/// plus kind and segment for the material kinds. None if the target does not resolve or does not
/// fit the kind (validate reports that).
[[nodiscard]] std::optional<std::string> quantity(const TargetIndex& index,
                                                  const ToleranceTarget& target,
                                                  ToleranceKind kind);

/// The class skip of a material kind for one material reference (ADR 0035, point 5): Metal for
/// CONST with kappa > 0, ConstAbbe for abbe at a real CONST material or VACUUM, otherwise
/// none. Catalog references cannot be judged without a catalog and give none.
enum class MaterialClass : std::uint8_t { Other, RealConst, Metal };
[[nodiscard]] MaterialClass classify_material(std::string_view reference);

/// True for the rigid kinds (decenter, shift, tilts).
[[nodiscard]] bool is_rigid(ToleranceKind kind) noexcept;

}  // namespace rtt::model::detail
