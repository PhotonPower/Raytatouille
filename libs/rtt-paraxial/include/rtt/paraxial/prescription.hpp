#pragma once

/// @file prescription.hpp
/// Paraxial system data as in a prescription report (#84): marginal and chief ray at every
/// event of a path, total track, paraxial working F-number, paraxial image-space NA and the
/// Lagrange invariant.
///
/// Conventions as in paraxial.hpp (decided for #7, docs/architecture.md, "Engine 1,
/// paraxial"): lengths in mm, global z along the optical axis, y the ray height, u = dy/dz the
/// slope in the global frame, n the signed index (negative while light travels towards -z), c
/// the vertex curvature in global orientation.
///
/// Rays: the same marginal and chief rays as seidel() (seidel.hpp, both built by one internal
/// function): the marginal ray from the axial object point to the rim of the paraxial entrance
/// pupil at +y, the chief ray of the maximum field through the centre of the entrance pupil.
///
/// Sources (docs/quellen.md), read in full text:
/// - J. Sasian, OPTI 517 lecture notes L4 "Seidel aberration coefficients", p. 24: refraction
///   invariant A = n i = n u + n y c, Lagrange invariant H (Zhe in the source) = n u_bar y -
///   n u y_bar.
/// - J. E. Greivenkamp, OPTI-502 lecture notes, Sec. 9 "Stops and Pupils", p. 9-41: optical
///   invariant n u_bar y - n u y_bar, invariant on refraction and transfer; with the marginal
///   and chief rays the Lagrange invariant. p. 9-34: NA = n sin U ~ n u; p. 9-35:
///   F/# = f / D_EP = 1 / (2 n' |u'|) for an object at infinity; p. 9-36: working
///   F/#_W = 1 / (2 NA) ~ 1 / (2 n |u|) for any cone of light; p. 9-37 for the image space.

#include <cstdint>
#include <optional>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/paraxial/paraxial.hpp"
#include "rtt/paraxial/seidel.hpp"

namespace rtt::paraxial {

/// One paraxial ray at one event.
struct PrescriptionRay {
  double y = 0.0;  ///< height at the vertex plane, mm
  double u = 0.0;  ///< slope dy/dz after the event (global frame)
  /// Paraxial angle of incidence i = u + y c in rad, with u the slope before the event and c
  /// the vertex curvature in global orientation (Sasian L4, p. 24: A = n i = n u + n y c). It is
  /// the ray slope minus the slope -y c of the surface normal at height y, so i > 0 when the
  /// ray is turned counter-clockwise from the normal in the y-z plane (y up, z to the right).
  /// Equal to A / n with A, n before the event (SeidelSurface::a). For plane surfaces
  /// (stop, detector) i = u.
  double i = 0.0;
};

/// Data at one event of the path.
struct PrescriptionSurface {
  std::uint32_t surface = 0;                ///< index into CompiledSystem::surfaces()
  double z = 0.0;                           ///< global z of the surface vertex, mm
  double n = 1.0;                           ///< signed index after the event (see paraxial.hpp)
  std::optional<PrescriptionRay> marginal;  ///< none if the marginal ray is not defined
  std::optional<PrescriptionRay> chief;     ///< none if the chief ray is not defined
  /// Lagrange invariant n' (u_bar' y - u' y_bar) after the event, mm (Greivenkamp p. 9-41,
  /// Sasian L4 p. 24); equal to Prescription::lagrange_invariant up to rounding. None without
  /// both rays.
  std::optional<double> lagrange;
};

/// Prescription data of a path at one wavelength.
struct Prescription {
  /// One entry per event of the path, in path order (also transmit events: stop, detector).
  std::vector<PrescriptionSurface> surfaces;
  /// Sum of |z_(k+1) - z_k| over consecutive events of the path, mm: the axial length from the
  /// first to the last event vertex. With mirrors this is the unfolded length along the path.
  double total_track = 0.0;
  /// Global z of the first event vertex minus global z of the object, mm; none for an object
  /// at infinity.
  std::optional<double> object_distance;
  /// Paraxial working F-number 1 / (2 |n'_K u'_K|) of the marginal ray after the last event
  /// (Greivenkamp p. 9-36, 9-37). For an object at infinity it equals EFL / EPD (p. 9-35).
  /// None for afocal systems, without a marginal ray, or if u'_K = 0. Paraxial only: a working
  /// F-number from a real marginal ray may follow later.
  std::optional<double> paraxial_working_f_number;
  /// Paraxial image-space NA |n'_K u'_K| of the marginal ray after the last event (Greivenkamp
  /// p. 9-34, sin U ~ u); also for afocal systems (about 0 for an object at infinity). None
  /// without a marginal ray.
  std::optional<double> paraxial_image_na;
  /// Lagrange invariant H = n (u_bar y - u y_bar) in object space, evaluated at the first
  /// event vertex with the slopes before it, mm (equal to Seidel::lagrange). None without both
  /// rays. Its sign follows the chief ray of seidel.hpp (field convention of #8): positive for
  /// a positive field angle, negative for a positive object height (the chief ray then falls
  /// towards the entrance pupil centre).
  std::optional<double> lagrange_invariant;
  std::optional<RayStart> marginal_start;  ///< marginal ray in object space, if defined
  std::optional<RayStart> chief_start;     ///< chief ray in object space, if defined
  /// First-order data of the path; lateral and angular magnification are taken from here.
  FirstOrder first_order;
};

/// Paraxial prescription data of a path (marginal and chief ray per event, system data).
///
/// Unlike seidel(), a missing stop or pupil is not an error (as in first_order()): the
/// affected values are none.
/// - Marginal ray: needs the entrance pupil diameter, with a stop from first_order(). For an
///   object at infinity it is the ray parallel to the axis at half that diameter, so it also
///   exists for an entrance pupil at infinity and without a stop; then the diameter comes from
///   the system aperture (EntrancePupilDiameter; ImageSpaceFNumber of a focal system,
///   EPD = EFL / F#). For a finite object it needs a stop and an entrance pupil at a finite
///   position not in the object plane. Otherwise none.
/// - Chief ray: needs a stop and an entrance pupil at a finite position not in the object
///   plane; the maximum field as in seidel(). If all field points lie on the axis it is the
///   zero ray (H = 0).
/// @param system     compiled system
/// @param path       path to evaluate
/// @param wavelength index into system.wavelengths_um() for the rays and n
/// @throws ParaxialError if the path is not rotationally symmetric, the stop aperture is not
///         circular, a path id or wavelength index does not exist, or, only when the chief
///         ray is constructed (stop and usable entrance pupil), the field definition is invalid
///         for it (field angle not in (-90, 90) degree, object height for an object at
///         infinity, paraxial image height without a finite paraxial image)
[[nodiscard]] Prescription prescription(const compile::CompiledSystem& system,
                                        compile::PathId path,
                                        std::uint16_t wavelength);

}  // namespace rtt::paraxial
