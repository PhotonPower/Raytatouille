#pragma once

// Internal to rtt-paraxial: construction of the paraxial marginal and chief rays shared by
// seidel() (seidel.cpp) and prescription() (prescription.cpp), so that both use the same rays
// (#84). Not installed, not part of the public API.

#include <cstdint>
#include <string_view>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/paraxial/paraxial.hpp"
#include "rtt/paraxial/seidel.hpp"

namespace rtt::paraxial::detail {

/// Curvature, conic constant and A4 of a surface in global orientation.
struct ShapeData {
  double c = 0.0;
  double conic = 0.0;
  double a4 = 0.0;
};

/// The paraxial trace (paraxial.cpp) uses c in global orientation: 1/R times the sign of the
/// local z axis along global z. The fourth-order sag c^3 r^4 / 8 (K + 1) + A4 r^4 changes sign
/// with the axis as well, so A4 flips with c while K, a ratio, stays.
[[nodiscard]] ShapeData shape_data(const compile::CompiledSurface& s);

/// Marginal ray in object space: from the axial object point to the rim of the entrance pupil
/// at +y (object at infinity: parallel at height r_ep through the pupil plane).
[[nodiscard]] RayStart marginal_start(bool object_at_infinity,
                                      double z_obj,
                                      double z_ep,
                                      double r_ep);

/// Chief ray of the maximum field (see seidel.hpp), in object space, through the centre of the
/// entrance pupil at z_ep; paraxial image heights and field angles with a finite object are
/// converted at the reference wavelength (#35, as rtt-trace since #31/#50).
/// @param fo     first-order data at `wavelength`; z_ep its entrance pupil
/// @param caller name for error messages ("paraxial: <caller>: ...")
/// @throws ParaxialError for field angles outside (-90, 90) degree, an object height with an
///         object at infinity, a paraxial image height without a finite paraxial image, or an
///         entrance pupil at the reference wavelength at infinity or in the object plane where
///         the field conversion needs it
[[nodiscard]] RayStart chief_start(const compile::CompiledSystem& system,
                                   compile::PathId path,
                                   std::uint16_t wavelength,
                                   const FirstOrder& fo,
                                   double z_ep,
                                   double z_obj,
                                   std::string_view caller);

}  // namespace rtt::paraxial::detail
