#pragma once

/// @file tolerance.hpp
/// Tolerance set of a system (ADR 0035, file section "tolerances", docs/dateiformat.md): class,
/// overrides, explicit items, compensators, criteria and Monte Carlo defaults. Pure data: the
/// enumeration of the class items and the semantic checks are in validate (#253), the class
/// limits and the analysis in rtt-tolerance (ADR 0036).

#include <array>
#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <variant>
#include <vector>

#include "rtt/model/optimization.hpp"

namespace rtt::model {

/// What a tolerance changes (ADR 0035, point 2). Units: mm for lengths, degree for the kinds
/// ending in Deg, dimensionless otherwise; `Param` and `Parameter` in the unit of their target.
enum class ToleranceKind : std::uint8_t {
  Radius,             ///< R' = R + dR, mm (file: "radius")
  CurvatureRelative,  ///< C' = C (1 + d), dimensionless (file: "curvature_relative")
  Conic,              ///< conic constant k' = k + dk (file: "conic")
  Thickness,          ///< offset on pose.position[2] of the target, mm (file: "thickness")
  Param,              ///< offset on the value of the attached Param (file: "param")
  Parameter,          ///< offset on a row of the parameter table (file: "parameter")
  DecenterX,          ///< rigid shift along x of the nominal frame, mm (file: "decenter_x")
  DecenterY,          ///< rigid shift along y of the nominal frame, mm (file: "decenter_y")
  ShiftZ,             ///< rigid shift along z of the nominal frame, mm (file: "shift_z")
  TiltXDeg,           ///< rigid rotation about x through the pivot, degree (file: "tilt_x_deg")
  TiltYDeg,           ///< rigid rotation about y through the pivot, degree (file: "tilt_y_deg")
  TiltZDeg,           ///< rigid rotation about z through the pivot, degree (file: "tilt_z_deg")
  Index,              ///< refractive index change dn_d, dimensionless (file: "index")
  Abbe,               ///< Abbe number change dnu_d, dimensionless (file: "abbe")
};

/// Distribution of a tolerance within [-T, T] (ADR 0035, point 2; ADR 0036).
enum class ToleranceDistribution : std::uint8_t {
  Normal,    ///< sigma = T/2, truncated at +-T (default, file: "normal")
  Uniform,   ///< uniform on [-T, T] (file: "uniform")
  Endpoint,  ///< -T or +T with probability 1/2 each (file: "endpoint")
};

/// Tolerance class (ADR 0035, point 5): a fixed set; the limits per class are data in
/// rtt-tolerance (#260).
enum class ToleranceClass : std::uint8_t {
  None,     ///< no generated items (file: "none", or the key missing)
  Coarse,   ///< file: "coarse"
  Normal,   ///< file: "normal"
  Precise,  ///< file: "precise"
};

/// The kind of anchor a tolerance target uses (ADR 0035, points 2 and 3).
enum class ToleranceAnchor : std::uint8_t {
  Surface,    ///< a SurfaceId (file key "surface")
  Node,       ///< the name of an element or an assembly, not the root (file key "node")
  Parameter,  ///< the name of a row of the parameter table (file key "parameter")
};

/// Target of an item, an override or a compensator: a stable anchor, never an index (ADR 0035,
/// point 3).
struct ToleranceTarget {
  ToleranceAnchor anchor = ToleranceAnchor::Node;
  std::string name;  ///< SurfaceId, node name or row name, depending on `anchor`
  /// JSON pointer relative to the object of the anchor in the edit form (ADR 0024), pointing at
  /// a Param object, e.g. "/shape/base/coefficients/0"; only with the kind Param and at
  /// compensators, never with the anchor Parameter.
  std::optional<std::string> param;
  /// 0-based segment of an element with several materials (ADR 0017); only with the kinds
  /// Index and Abbe.
  std::optional<std::uint32_t> segment;
  bool operator==(const ToleranceTarget&) const = default;
};

/// An explicit tolerance (ADR 0035, point 2).
struct ToleranceItem {
  /// Not empty, unique in the set, not starting with "class:" (reserved for generated items);
  /// determines the random sequence of the item (ADR 0036, point 5).
  std::string id;
  ToleranceTarget target;
  ToleranceKind kind = ToleranceKind::Thickness;
  double limit = 0.0;  ///< T > 0 in the unit of the kind; the perturbation lies in [-T, T]
  ToleranceDistribution distribution = ToleranceDistribution::Normal;
  bool operator==(const ToleranceItem&) const = default;
};

/// Override of the class for a node and its subtree (file: {"node": ..., "class": ...}).
struct ToleranceClassOverride {
  std::string node;                                       ///< name of an element or an assembly
  ToleranceClass tolerance_class = ToleranceClass::None;  ///< None switches generation off
  bool operator==(const ToleranceClassOverride&) const = default;
};

/// A single value for one generated kind at one target (file: {"node" | "surface": ...,
/// "kind": ..., "limit": ..., "distribution": ..., "segment": ...}).
struct ToleranceValueOverride {
  ToleranceTarget target;  ///< anchor Node or Surface; `param` is not used
  ToleranceKind kind = ToleranceKind::Thickness;
  double limit = 0.0;  ///< as ToleranceItem::limit
  ToleranceDistribution distribution = ToleranceDistribution::Normal;
  bool operator==(const ToleranceValueOverride&) const = default;
};

/// The pivot of the tilts of one target, not inherited (file: {"node" | "surface": ...,
/// "pivot": [x, y, z]}).
struct TolerancePivotOverride {
  ToleranceTarget target;         ///< anchor Node or Surface; `param` and `segment` not used
  std::array<double, 3> pivot{};  ///< mm, in the coordinates of the target
  bool operator==(const TolerancePivotOverride&) const = default;
};

/// One entry of "overrides": exactly one of the three forms (ADR 0035, point 5).
using ToleranceOverride =
    std::variant<ToleranceClassOverride, ToleranceValueOverride, TolerancePivotOverride>;

/// A compensator (ADR 0035, point 6): an unbound Param (anchor Node or Surface with `param`) or a
/// row of the form value (anchor Parameter), with a travel relative to the value of the trial
/// before compensation, in the unit of the target.
struct Compensator {
  ToleranceTarget target;
  double travel_min = 0.0;  ///< < 0
  double travel_max = 0.0;  ///< > 0
  bool operator==(const Compensator&) const = default;
};

/// A criterion (ADR 0035, point 7): an operand of ADR 0030 with bounds on its value. The
/// operand's target and weight are not used (the file has neither key; validate rejects values
/// other than the defaults).
struct Criterion {
  std::string id;  ///< not empty, unique among the criteria
  Operand operand;
  std::optional<double> min;         ///< absolute lower bound, unit of the operand
  std::optional<double> max;         ///< absolute upper bound, unit of the operand
  std::optional<double> max_change;  ///< bound on |value - nominal value|, > 0
  bool operator==(const Criterion&) const = default;
};

/// Monte Carlo defaults (ADR 0035, point 8); the API may override both.
struct MonteCarloDefaults {
  std::optional<std::uint64_t> seed;    ///< in [0, 2^53 - 1] for the file (validate)
  std::optional<std::uint64_t> trials;  ///< in [1, 2^32] (validate)
  bool operator==(const MonteCarloDefaults&) const = default;
};

/// The tolerance set (ADR 0035, point 1). Empty: no tolerances.
struct ToleranceSet {
  ToleranceClass tolerance_class = ToleranceClass::None;  ///< file key "class"
  std::vector<ToleranceOverride> overrides;
  std::vector<ToleranceItem> items;
  std::vector<Compensator> compensators;
  /// Merit function of the compensators; none: the system's `optimization` is used.
  std::optional<Optimization> compensator_merit;
  std::vector<Criterion> criteria;
  std::optional<MonteCarloDefaults> monte_carlo;
  [[nodiscard]] bool empty() const noexcept {
    return tolerance_class == ToleranceClass::None && overrides.empty() && items.empty() &&
           compensators.empty() && !compensator_merit && criteria.empty() && !monte_carlo;
  }
  bool operator==(const ToleranceSet&) const = default;
};

/// The file name of a kind ("radius", ..., "abbe"; ADR 0035, point 2), as used in the ids of
/// generated items.
[[nodiscard]] std::string_view tolerance_kind_name(ToleranceKind kind) noexcept;

/// 64-bit FNV-1a hash of the UTF-8 bytes of a tolerance id (ADR 0036, point 5): hash =
/// offset_basis; per octet hash = (hash xor octet) * FNV_Prime mod 2^64, with offset_basis
/// 0xCBF29CE484222325 and FNV_Prime 0x00000100000001B3 (RFC 9923, "The FNV Non-Cryptographic
/// Hash Algorithm", 2026, Sec. 2 and Sec. 5; draft-eastlake-fnv-35, same sections). Two ids of
/// one set with equal hashes are an error (validate: tolerances.id_hash_collision).
[[nodiscard]] std::uint64_t tolerance_id_hash(std::string_view id) noexcept;

}  // namespace rtt::model
