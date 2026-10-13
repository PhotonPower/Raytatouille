#pragma once

/// @file parabasal.hpp
/// The paraxial operands of a perturbed trial from parabasal rays (ADR 0030, addendum #256;
/// ADR 0036, point 2), private to rtt-optim.

#include <cstdint>
#include <limits>
#include <string>

#include "overlay_basis.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/model/optimization.hpp"

namespace rtt::optim::detail {

/// Value of a paraxial operand, NaN with the reason if it is not defined.
struct ParabasalValue {
  double value = std::numeric_limits<double>::quiet_NaN();
  std::string undefined;
};

/// efl, bfl, image_fnumber or magnification of `op` in the trial `trial` from five rays with
/// real aiming and the nominal aiming data: the chief ray and the rays at (+-delta, 0) and
/// (0, +-delta) (addendum #256).
/// @param trial          the compiled trial (with its overlay)
/// @param trial_infinite the trial's variant "object at infinity" (variant_of()), compiled with
///                       the same overlay; required for efl and bfl with a finite object, else
///                       unused
/// @param nominal        the nominal data of the operand's path and wavelength
/// @param reference      the nominal aiming data of the path (analysis::NominalReference)
/// @param object_at_infinity whether the system's object is at infinity
/// @throws as trace::aim_ray() with a reference (input errors, a reference that does not fit)
[[nodiscard]] ParabasalValue parabasal_value(const model::FirstOrderOperand& op,
                                             const compile::CompiledSystem& trial,
                                             const compile::CompiledSystem* trial_infinite,
                                             const WavelengthData& nominal,
                                             const analysis::NominalReference& reference,
                                             bool object_at_infinity,
                                             compile::PathId path,
                                             std::uint16_t wavelength,
                                             double delta);

}  // namespace rtt::optim::detail
