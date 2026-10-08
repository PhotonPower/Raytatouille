#pragma once

/// @file seidel.hpp
/// Third-order (Seidel) aberration sums S_I ... S_V and the first-order chromatic terms C_L, C_T
/// per surface and for the system, from the paraxial marginal and chief rays
/// (docs/architecture.md, "Engine 1, paraxial"; sources in docs/quellen.md).
///
/// Conventions (paraxial conventions of #7, see paraxial.hpp):
/// - y is the ray height in mm, u = dy/dz the global slope, n the signed index (negative while
///   light travels towards -z, n' = -n at every reflection), c the vertex curvature in global
///   orientation in 1/mm. The source formulas (J. Sasian, OPTI 517 lecture notes) use the same
///   symbols with u = -y/s (a slope, positive for a ray rising towards +z) and n' = -n at
///   mirrors; they apply unchanged.
/// - Refraction invariants A = n (u + y c) (marginal ray) and A_bar = n (u_bar + y_bar c) (chief
///   ray), Lagrange invariant H = n (u_bar y - u y_bar) = A_bar y - A y_bar. A, A_bar and H do
///   not change at a refraction or reflection, H does not change along the path either.
/// - Marginal ray: from the axial object point to the rim of the paraxial entrance pupil at +y
///   (object at infinity: parallel to the axis at height +r_EP). Chief ray: from the maximum
///   field through the centre of the entrance pupil, in the meridional (y, z) plane, with the
///   field value counted positive as defined by the field type (docs/architecture.md,
///   "Feldwinkel und Pupille"): a field angle theta > 0 gives a chief ray rising towards +y
///   (slope tan theta), an object height h > 0 gives an object point at +y, a paraxial image
///   height h' > 0 an image point at +y. S_II and S_V change sign with the chief ray, S_I, S_III
///   and S_IV do not. For a finite object with an angle field the object point lies at -y (#8);
///   S_II and S_V then have the opposite sign of an equivalent object-height field.
/// - Normalisation: the wavefront aberration at normalised pupil vector rho (|rho| = 1 at the
///   rim of the pupil) and normalised field vector eta (|eta| = 1 at the maximum field) is
///     W = 1/8 S_I (rho.rho)^2 + 1/2 S_II (eta.rho)(rho.rho) + 1/2 S_III (eta.rho)^2
///       + 1/4 (S_III + S_IV)(eta.eta)(rho.rho) + 1/2 S_V (eta.eta)(eta.rho)
///       + 1/2 C_L (rho.rho) + C_T (eta.rho)
///   in mm (Sasian, OPTI 517 L4 p. 23: W040 = S_I/8, W131 = S_II/2, W222 = S_III/2,
///   W220 = (S_IV + S_III)/4, W311 = S_V/2, d_lambda W020 = C_L/2, d_lambda W111 = C_T).
///   All S and C values are in mm. Checked against real rays in our convention (#30): the
///   transverse aberration in the paraxial image plane is eps = (dW/drho) / (n'_K u'_K) with
///   n'_K u'_K of the marginal ray in image space. For a converging image (n'_K u'_K < 0)
///   S_I > 0 therefore means under-corrected spherical aberration (marginal focus before the
///   paraxial focus along the light).
/// - Chromatic terms use dn = n(first) - n(second) of the chosen wavelength pair, with the sign of
///   n (so dn'/n' = dn/n at a mirror and mirrors contribute nothing), and the rays and n of the
///   wavelength passed to seidel().

#include <cstdint>
#include <optional>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/paraxial/paraxial.hpp"

namespace rtt::paraxial {

/// Seidel sums and chromatic terms of one surface or of a system, all in mm (normalisation in
/// the file comment).
struct SeidelTerms {
  double s1 = 0.0;   ///< S_I, spherical aberration
  double s2 = 0.0;   ///< S_II, coma
  double s3 = 0.0;   ///< S_III, astigmatism
  double s4 = 0.0;   ///< S_IV, Petzval field curvature
  double s5 = 0.0;   ///< S_V, distortion
  double c_l = 0.0;  ///< C_L, longitudinal (axial) colour; 0 without a wavelength pair
  double c_t = 0.0;  ///< C_T, transverse (lateral) colour; 0 without a wavelength pair
};

/// Paraxial data of one surface as needed by surface_seidel(). Heights and slopes are those of
/// the marginal and chief ray at the vertex plane, in the medium before the surface.
struct SeidelSurfaceInput {
  double c = 0.0;         ///< vertex curvature in global orientation, 1/mm
  double conic = 0.0;     ///< conic constant K (0 = sphere, -1 = paraboloid)
  double a4 = 0.0;        ///< fourth-order asphere coefficient in global orientation, 1/mm^3
  double n = 1.0;         ///< signed index before the surface
  double n_after = 1.0;   ///< signed index after the surface (-n for a reflection)
  double y = 0.0;         ///< marginal ray height, mm
  double u = 0.0;         ///< marginal ray slope dy/dz before the surface
  double y_bar = 0.0;     ///< chief ray height, mm
  double u_bar = 0.0;     ///< chief ray slope dy/dz before the surface
  double dn = 0.0;        ///< signed index difference n(first) - n(second) before the surface
  double dn_after = 0.0;  ///< the same after the surface
};

/// Seidel sums and chromatic terms of one surface (formulas and sources in seidel.cpp).
///
/// Spherical part (Sasian, OPTI 517 L4 p. 23/24) with Delta(x) = x' - x across the surface and
/// P = c Delta(1/n):
///   S_I = -A^2 y Delta(u/n), S_II = -A A_bar y Delta(u/n), S_III = -A_bar^2 y Delta(u/n),
///   S_IV = -H^2 P, S_V = -A_bar [A_bar^2 Delta(1/n^2) y - (H + A_bar y) y_bar P],
///   C_L = A y Delta(dn/n), C_T = A_bar y Delta(dn/n).
/// Aspheric cap (Sasian, OPTI 518 L14 p. 16/17): with a = (K c^3 + 8 A4) Delta(n) y^4 the sums
/// grow by a, (y_bar/y) a, (y_bar/y)^2 a, 0 and (y_bar/y)^3 a; computed without dividing by y.
/// The u after the surface follows from n'u' = n u - y c (n' - n). A surface with n' = n (no
/// reflection) gives all terms 0.
[[nodiscard]] SeidelTerms surface_seidel(const SeidelSurfaceInput& in);

/// Pair of wavelength indices for the chromatic terms: dn = n(first) - n(second), e.g. F and C.
struct ChromaticPair {
  std::uint16_t first = 0;
  std::uint16_t second = 0;
};

/// Paraxial start of a ray in object space.
struct RayStart {
  double z = 0.0;  ///< global z of the start plane, mm
  double y = 0.0;  ///< height there, mm
  double u = 0.0;  ///< slope dy/dz
};

/// Contribution of one event of the path.
struct SeidelSurface {
  std::uint32_t surface = 0;  ///< index into CompiledSystem::surfaces()
  double y = 0.0;             ///< marginal ray height at the vertex plane, mm
  double y_bar = 0.0;         ///< chief ray height at the vertex plane, mm
  double a = 0.0;             ///< refraction invariant A = n (u + y c) of the marginal ray
  double a_bar = 0.0;         ///< refraction invariant A_bar of the chief ray
  /// Lagrange invariant n' (u_bar' y - u' y_bar) evaluated after the event, mm; equal to
  /// Seidel::lagrange up to rounding.
  double lagrange = 0.0;
  SeidelTerms terms;  ///< contribution of this event, mm
};

/// Seidel sums of a path.
struct Seidel {
  /// One entry per event of the path, in path order. Events that neither reflect nor change
  /// the index (stop, detector, ...) contribute 0.
  std::vector<SeidelSurface> surfaces;
  SeidelTerms sum;                         ///< sum over all events, mm
  double lagrange = 0.0;                   ///< H = n (u_bar y - u y_bar) in object space, mm
  RayStart marginal;                       ///< marginal ray used, in object space
  RayStart chief;                          ///< chief ray of the maximum field used, in object space
  std::optional<ChromaticPair> chromatic;  ///< pair used for C_L and C_T, if any
};

/// Seidel sums S_I ... S_V per event and for the path, and C_L, C_T if a wavelength pair is
/// given. Conventions, normalisation and sign of the chief ray in the file comment.
///
/// - Pupil: the paraxial entrance pupil of first_order() (position and diameter).
/// - Maximum field: the field point with the largest radial value, tan theta =
///   hypot(tan theta_x, tan theta_y) for field angles and hypot(x, y) for object or paraxial
///   image heights, placed in the meridional plane at +y. A paraxial image height is converted
///   with the paraxial chief ray of unit field value, linear in the field (as in rtt-trace,
///   #8). As in rtt-trace (#31, #50), a paraxial image height and a field angle with a finite
///   object are converted at the reference wavelength, into a slope or an object height; the
///   chief ray at `wavelength` passes through the centre of its own entrance pupil (#35). If
///   all field points lie on the axis, the chief ray is 0 and only S_I and C_L can be non-zero.
/// - Surface shape: vertex curvature, conic constant and A4 of the base shape (Conic,
///   EvenAsphere); higher asphere terms do not enter the third order.
/// @param system     compiled system
/// @param path       path to evaluate
/// @param wavelength index into system.wavelengths_um() for the paraxial rays and n
/// @param chromatic  wavelength pair for C_L and C_T; without it they are 0
/// @throws rtt::compile::NoStopError (a std::invalid_argument) if the path has no stop (checked
///         by rtt::compile::require_stop after the checks of first_order; ADR 0022)
/// @throws ParaxialError if the path is not rotationally symmetric, the stop aperture is not
///         circular, the entrance pupil is at infinity or has no diameter, an object height is
///         given for an object at infinity, a field angle is not in (-90, 90) degree, a paraxial
///         image height is given without a finite paraxial image, the field conversion needs
///         the entrance pupil at the reference wavelength and it lies at infinity or in the
///         object plane, or a path id or wavelength index does not exist
[[nodiscard]] Seidel seidel(const compile::CompiledSystem& system,
                            compile::PathId path,
                            std::uint16_t wavelength,
                            std::optional<ChromaticPair> chromatic = std::nullopt);

}  // namespace rtt::paraxial
