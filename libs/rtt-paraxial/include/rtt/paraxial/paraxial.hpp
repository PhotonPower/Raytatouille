#pragma once

/// @file paraxial.hpp
/// Paraxial (first-order) y-nu trace and first-order system data for rotationally symmetric
/// systems (docs/architecture.md, "Engine 1, paraxial").
///
/// Conventions (decided for #7, recorded in docs/architecture.md, "Engine 1, paraxial"; units
/// and axes as in "Konventionen"):
/// - Lengths in mm, positions as global z coordinates; the optical axis is the global z axis.
/// - y is the ray height in mm, u = dy/dz the slope in the global frame (paraxial: tan u ~ u).
/// - The index n carries a sign: positive while light travels towards +z, negative while it
///   travels towards -z. It changes sign at every reflection (n' = -n). Then n * u equals |n|
///   times the direction cosine along the propagation direction (the optical direction cosine),
///   so refraction n'u' = nu - y phi with phi = c (n' - n) and the transfer
///   y' = y + (z' - z) u with global z hold unchanged after mirrors.
/// - c is the vertex curvature in global orientation: 1/R times the sign of the surface's local
///   z axis along global z.
/// - n is the real part of the complex index at the chosen wavelength; kappa is ignored.
/// - Focal lengths are positive for converging systems, for lenses and mirrors alike (a concave
///   mirror has R < 0 for light from -z and EFL = -R/2 = |R|/2).

#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "rtt/compile/compiled_system.hpp"

namespace rtt::paraxial {

/// Thrown when paraxial data are requested for a path that is not rotationally symmetric about
/// the global z axis (decentred or tilted surface, a diffraction order other than 0 (ADR 0025;
/// order 0 at a phase surface is the surface without phase layer), birefringent events) or for an
/// invalid path id or wavelength index. The message names the surface; where the place in the
/// system file is known, surface() and location() give it (ADR 0022).
class ParaxialError : public std::runtime_error {
 public:
  using std::runtime_error::runtime_error;

  /// Error at `surface`: surface() is its id, location() its JSON pointer.
  ParaxialError(const std::string& message, const compile::CompiledSurface& surface)
      : std::runtime_error(message), surface_(surface.id), location_(surface.location) {}

  /// Error at a place in the system file that is not a surface, e.g. "/fields/points/3".
  ParaxialError(const std::string& message, std::string location)
      : std::runtime_error(message), location_(std::move(location)) {}

  /// Id of the surface the error is about, if any.
  [[nodiscard]] const std::optional<model::SurfaceId>& surface() const noexcept { return surface_; }
  /// JSON pointer into the system file, e.g. "/root/children/1/surfaces/0", if known.
  [[nodiscard]] const std::optional<std::string>& location() const noexcept { return location_; }

 private:
  std::optional<model::SurfaceId> surface_;
  std::optional<std::string> location_;
};

/// Paraxial ray just after one event of a path.
struct RayAtEvent {
  double z = 0.0;  ///< global z of the surface vertex, mm
  double y = 0.0;  ///< ray height at the vertex plane, mm
  double u = 0.0;  ///< slope dy/dz after the event (global frame)
  double n = 1.0;  ///< signed index after the event (see file comment)
};

/// Traces one paraxial ray along a path.
/// @param system     compiled system
/// @param path       path to follow
/// @param wavelength index into system.wavelengths_um()
/// @param z0         global z where the ray starts, mm; the ray starts in the medium before the
///                   first event and travels towards +z
/// @param y0         height at z0, mm
/// @param u0         slope dy/dz at z0
/// @return the ray after each event of the path, in path order
/// @throws ParaxialError if the path is not rotationally symmetric, or if the path id or the
///         wavelength index does not exist (see ParaxialError)
[[nodiscard]] std::vector<RayAtEvent> trace_ray(const compile::CompiledSystem& system,
                                                compile::PathId path,
                                                std::uint16_t wavelength,
                                                double z0,
                                                double y0,
                                                double u0);

/// Paraxial image of the aperture stop: position and diameter.
struct Pupil {
  std::optional<double> z;         ///< global z, mm; none if the pupil is at infinity
  std::optional<double> diameter;  ///< mm; none if the system aperture does not define it
};

/// First-order data of a path at one wavelength. Optional values are absent when they are not
/// defined (afocal system, object at infinity, no stop on the path, ...).
struct FirstOrder {
  double object_index = 1.0;  ///< |n| in object space
  double image_index = 1.0;   ///< |n| in image space
  /// +1 if light leaves the system towards +z, -1 towards -z (odd number of reflections).
  int image_direction = 1;

  /// Power Phi = -n'_K u'_K / y_1 for a ray from infinity, 1/mm, as computed; for an afocal
  /// system 0 up to rounding (below the afocal threshold of first_order()).
  double power = 0.0;
  std::optional<double> efl;                 ///< 1/Phi, mm; none if afocal
  std::optional<double> front_focal_length;  ///< |n| / Phi in object space (H -> F), mm
  std::optional<double> rear_focal_length;   ///< |n'| / Phi in image space (H' -> F'), mm
  std::optional<double> ffl;                 ///< first vertex -> F, positive if F lies upstream, mm
  std::optional<double> bfl;            ///< last vertex -> F', positive if F' lies downstream, mm
  std::optional<double> front_focal_z;  ///< global z of F, mm
  std::optional<double> rear_focal_z;   ///< global z of F', mm
  std::optional<double> front_principal_z;  ///< global z of H, mm
  std::optional<double> rear_principal_z;   ///< global z of H', mm

  /// Global z of the paraxial image of the axial object point (= rear_focal_z for an object at
  /// infinity); none if the image is at infinity.
  std::optional<double> image_z;
  /// Lateral magnification m = n u / (n' u') of the axial marginal ray; finite object only.
  std::optional<double> lateral_magnification;
  /// Ratio of the chief-ray slopes measured along the propagation direction, u'/u; needs a stop.
  /// For afocal systems this is the classic angular magnification.
  std::optional<double> angular_magnification;

  std::optional<Pupil> entrance_pupil;  ///< none if the path has no stop
  std::optional<Pupil> exit_pupil;      ///< none if the path has no stop
};

/// First-order data of a path (EFL, focal and principal points, pupils, magnification).
///
/// - Vertices: the first and last events that change the ray (reflection or change of index);
///   transmit events (stop, detector, thin elements) do not count.
/// - Afocal: |Phi| <= 16 N u (S + E) with u = 2^-53, N the number of events, S the power of the
///   same y-nu trace with every term by its magnitude and E the bound of the error from the
///   coordinates, traced alongside with the same magnitudes (#35, #280): Phi is zero up to the
///   rounding of the trace and of the coordinates. E grows linearly with the coordinates.
///   The threshold depends on the absolute position of the path, because a coordinate z carries
///   a rounding of about u |z|: a lens with Phi = 5e-13 / mm (R = +-1 mm) is focal at
///   z = 10 mm and afocal at z = 1000 mm. It is cautious for very small systems: a lens of
///   R = +-3 um at z = 10 mm counts as afocal for |Phi| up to about 1.3e-9 / mm.
/// - Stop: the first event at a surface of a Stop element. Its aperture must be circular (an
///   annulus counts with its outer radius). Entrance pupil = paraxial image of the stop through
///   the surfaces before it, exit pupil = through the surfaces after it.
/// - Pupil diameters from the system aperture: EntrancePupilDiameter directly; StopSize from
///   the stop radius; ImageSpaceFNumber as F-number at infinite conjugates, EPD = EFL / F#, also
///   for finite objects (the working F-number follows later); ObjectSpaceNA from u = NA / n in
///   object space, finite objects only (otherwise no diameter).
/// @param system     compiled system
/// @param path       path to evaluate
/// @param wavelength index into system.wavelengths_um()
/// @throws ParaxialError if the path is not rotationally symmetric, the stop aperture is not
///         circular, or the path id or the wavelength index does not exist
[[nodiscard]] FirstOrder first_order(const compile::CompiledSystem& system,
                                     compile::PathId path,
                                     std::uint16_t wavelength);

/// Section of a path for the paraxial computation of anamorphic systems (ADR 0042, point 2):
/// the plane x-z or y-z of the global frame through the z axis.
enum class Section : std::uint8_t {
  x,  ///< the x-z plane; the ray height is the global x, the slope dx/dz
  y,  ///< the y-z plane; the ray height is the global y, the slope dy/dz
};

/// Traces one paraxial ray along a path in one section (ADR 0042, point 2), with the same y-nu
/// trace as trace_ray(). Conventions as in the file comment, with x in place of y for
/// Section::x.
///
/// Per surface (axis +-z and position on the z axis as in first_order()):
/// - rotationally symmetric to second order (d2z/dx2 = d2z/dy2, d2z/dxdy = 0 in its local
///   frame): c exactly as in trace_ray(), the same in both sections;
/// - anamorphic: the global Hessian H_g = sigma Q H_l Q^T, with Q the 2x2 transverse block of
///   the rotation of the surface and sigma = +-1 the sign of its local z axis along global z;
///   the principal axes must be aligned with x and y (|H_g,xy| <= 1e-12 max(|H_g,xx|,
///   |H_g,yy|)), and c = H_g,xx in section x, H_g,yy in section y (1/mm);
/// - ideal cylinder lens (ADR 0031, point 4): power phi b_s^2 with phi = |n|/f and b = Q b_l
///   the global direction of power (b_l = (-sin psi, cos psi) in the local frame); b must lie
///   along x or y (|b_x b_y| <= 1e-12).
/// @param system     compiled system
/// @param path       path to follow
/// @param wavelength index into system.wavelengths_um()
/// @param z0         global z where the ray starts, mm; the ray starts in the medium before the
///                   first event and travels towards +z
/// @param h0         height at z0 in the section (global x or y), mm
/// @param u0         slope dh/dz at z0 in the global frame
/// @param section    the section
/// @return the ray after each event of the path, in path order; RayAtEvent::y is the height h
/// @throws ParaxialError as trace_ray(), except for anamorphic surfaces and the ideal cylinder
///         lens, and at a surface whose principal axes (or direction of power) are not aligned
///         with x and y
[[nodiscard]] std::vector<RayAtEvent> trace_ray_section(const compile::CompiledSystem& system,
                                                        compile::PathId path,
                                                        std::uint16_t wavelength,
                                                        double z0,
                                                        double h0,
                                                        double u0,
                                                        Section section);

/// First-order data of a path in one section (ADR 0042, point 2): the same computation as
/// first_order() with the curvatures and powers of trace_ray_section(). Pupils, magnification
/// and the afocal threshold are defined as there; a section without power is afocal (`efl`
/// empty). Positions are global z coordinates in mm, heights and pupil diameters are measured
/// along x (Section::x) or y (Section::y).
///
/// For a path without anamorphic surfaces (is_anamorphic() false) both sections are bit for
/// bit first_order(). For an anamorphic surface the curvature of a section is rounded relative
/// to the whole Hessian, not to the section alone (a cylinder rotated by 90 degree about z
/// gets the curvature cos^2(90 degree) c ~ 3.7e-33 c in the section without power); the
/// afocal threshold counts its magnitude as |H_xx| + 2 |H_xy| + |H_yy| of the local Hessian
/// (resp. |n|/|f| for the ideal cylinder lens) instead of |c|, and a factor b_s^2 <= 1e-12 of
/// the ideal cylinder lens is 0 (ADR 0042, point 2 and addendum #315 PR C1).
/// @param system     compiled system
/// @param path       path to evaluate
/// @param wavelength index into system.wavelengths_um()
/// @param section    the section
/// @throws ParaxialError as trace_ray_section(), or if the stop aperture is not circular
[[nodiscard]] FirstOrder first_order_section(const compile::CompiledSystem& system,
                                             compile::PathId path,
                                             std::uint16_t wavelength,
                                             Section section);

/// True if an event of the path is at a surface with d2z/dx2 != d2z/dy2 or d2z/dxdy != 0 at
/// its vertex in its local frame (biconic, toroid, cylinder, ...), or at an ideal cylinder lens
/// (ADR 0042, point 2, "Erkennung"). Such a path has no first_order(); first_order_section()
/// may compute it. The comparisons are exact, as in the vertex rule of ADR 0041, point 2.
/// Surfaces with order < 2 (a cone at the vertex: axicon, odd asphere with A1 != 0) do not
/// count, because their Hessian has no meaning; first_order() and first_order_section() reject
/// them alike ("no paraxial model at the vertex").
/// @param system compiled system
/// @param path   path to inspect
/// @throws ParaxialError if the path id does not exist
[[nodiscard]] bool is_anamorphic(const compile::CompiledSystem& system, compile::PathId path);

}  // namespace rtt::paraxial
