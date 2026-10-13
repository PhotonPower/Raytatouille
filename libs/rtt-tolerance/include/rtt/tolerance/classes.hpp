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
/// lengths, degree for the tilts, dimensionless for curvature_relative, conic, index and abbe;
/// ADR 0035, point 2). The perturbation lies in [-T, T].
/// @throws std::invalid_argument for ToleranceClass::None and for a kind that no class generates
///         at that level (ADR 0035, point 5: e.g. radius, tilt_z_deg, param, index at an
///         assembly)
[[nodiscard]] double class_limit(model::ToleranceClass tolerance_class,
                                 model::ToleranceKind kind,
                                 model::ToleranceLevel level);

}  // namespace rtt::tolerance
