#pragma once

/// @file classes.hpp
/// The limits of the tolerance classes (ADR 0035, point 5; ADR 0036, point 9; #260): one limit
/// per class, kind and level. Which kinds a class generates at which target is decided in
/// rtt-model (model::class_tolerances); this table only gives the numbers. The values and their
/// origin are in the addendum #260 of ADR 0036.

#include "rtt/model/tolerance.hpp"
#include "rtt/model/tolerance_classes.hpp"

namespace rtt::tolerance {

/// The limit T of `kind` at `level` in class `tolerance_class`, in the unit of the kind (mm for
/// lengths, degree for the tilts, dimensionless for curvature_relative, conic and index; ADR
/// 0035, point 2). The perturbation lies in [-T, T]. For abbe it is the share p of nu_d (SCHOTT
/// TIE-29, table 5-1); expand() gives the item the limit p * nu_d of its glass (ADR 0035,
/// addendum #260).
/// @throws std::invalid_argument for ToleranceClass::None and for a kind that no class generates
///         at that level (ADR 0035, point 5: e.g. radius, tilt_z_deg, param, index at an
///         assembly)
[[nodiscard]] double class_limit(model::ToleranceClass tolerance_class,
                                 model::ToleranceKind kind,
                                 model::ToleranceLevel level);

/// The share p of nu_d that limits abbe in class `tolerance_class` (SCHOTT TIE-29, section 5,
/// table 5-1: coarse step 4, normal step 3, precise step 1, i.e. 0.008, 0.005, 0.002;
/// dimensionless). The limit of a generated abbe item is p * nu_d of its glass (expand(); ADR
/// 0035, addendum #260). Equal to class_limit(tolerance_class, Abbe, Element).
/// @throws std::invalid_argument for ToleranceClass::None
[[nodiscard]] double abbe_share(model::ToleranceClass tolerance_class);

}  // namespace rtt::tolerance
