#pragma once

/// @file tolerance_classes.hpp
/// The items a tolerance class generates (ADR 0035, point 5): which kinds at which targets,
/// with stable ids, the effective class and the overrides that apply, and why a kind is skipped.
/// This depends only on the tree, the classes and the overrides, not on the limits of a class;
/// the limits are data in rtt-tolerance (#260). validate uses it for the hash collisions of
/// explicit and generated ids.

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "rtt/model/system.hpp"
#include "rtt/model/tolerance.hpp"

namespace rtt::model {

/// Level of a generated item; a class has one limit per class, kind and level (ADR 0035,
/// point 5).
enum class ToleranceLevel : std::uint8_t { Surface, Element, Assembly };

/// Why a class does not generate a kind at a target (ADR 0035, point 5). Not an error: the
/// reason is part of the result of the generation.
enum class ClassSkip : std::uint8_t {
  None,              ///< generated
  Plane,             ///< curvature_relative or conic at a surface without a radius (plane)
  Sphere,            ///< conic at a sphere (conic constant 0)
  Bound,             ///< the Param is bound to a row of the parameter table (ADR 0029)
  NoMaterialInBeam,  ///< index or abbe at an element other than Lens and Plate
  NoMaterial,        ///< index or abbe at a Lens or Plate without isotropic material
  Crystal,           ///< index or abbe at a crystal (ADR 0026)
  Metal,             ///< index or abbe at a CONST material with kappa > 0
  ConstAbbe,         ///< abbe at a real CONST material or VACUUM (n_F - n_C = 0)
};

/// One kind at one target that a class would generate.
struct ClassTolerance {
  /// "class:node:<name>:<kind>" or "class:surface:<SurfaceId>:<kind>", with ":<segment>" for
  /// the material kinds of an element with segment_materials (ADR 0035, point 5).
  std::string id;
  ToleranceTarget target;  ///< anchor Node or Surface; segment for material kinds
  ToleranceKind kind = ToleranceKind::Thickness;
  ToleranceLevel level = ToleranceLevel::Element;
  /// The effective class: the class override of the target or its nearest ancestor, else the
  /// class of the system; never None (targets under None are not listed).
  ToleranceClass tolerance_class = ToleranceClass::Normal;
  /// Index into ToleranceSet::overrides of the single value for this target, kind and segment.
  std::optional<std::size_t> value_override;
  /// Index into ToleranceSet::items of the explicit item that perturbs the same quantity
  /// (ADR 0035, point 2); it takes precedence.
  std::optional<std::size_t> replaced_by;
  ClassSkip skip = ClassSkip::None;
  /// True if the item takes part in an analysis: not skipped and not replaced.
  [[nodiscard]] bool generated() const noexcept { return skip == ClassSkip::None && !replaced_by; }
  bool operator==(const ClassTolerance&) const = default;
};

/// What became of a single-value override (ADR 0035, point 5).
enum class OverrideEffect : std::uint8_t {
  Applied,       ///< sets limit and distribution of a generated item
  ClassNone,     ///< no class at the target (none, or no class at all)
  Skipped,       ///< the class skips the kind there; see `skip`
  Replaced,      ///< an explicit item perturbs the same quantity
  NotGenerated,  ///< no class generates the kind there, or the target is unknown (an error
                 ///< of validate: tolerances.target_invalid, tolerances.unknown_target)
};

struct OverrideStatus {
  std::size_t index = 0;  ///< into ToleranceSet::overrides
  OverrideEffect effect = OverrideEffect::Applied;
  ClassSkip skip = ClassSkip::None;  ///< the reason for Skipped
  bool operator==(const OverrideStatus&) const = default;
};

/// The result of the generation.
struct ClassTolerances {
  /// In tree order (pre-order: assembly before its children, element before its surfaces,
  /// children and surfaces in list order), per target the kinds in the order of ADR 0035,
  /// point 5, material kinds per segment ascending. The root and targets under the class None
  /// are not listed.
  std::vector<ClassTolerance> items;
  /// One entry per single-value override, in the order of ToleranceSet::overrides.
  std::vector<OverrideStatus> value_overrides;
};

/// The items the classes of `system.tolerances` generate (ADR 0035, point 5). Unknown targets
/// and kinds that no class generates are left to validate. A catalog metal is not recognised
/// (no catalogs here); compile reports it (tolerances.material_unsupported).
[[nodiscard]] ClassTolerances class_tolerances(const System& system);

}  // namespace rtt::model
