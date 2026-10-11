#pragma once

/// @file jones.hpp
/// Jones pupil and its polarization maps (ADR 0034, points 4 to 6, with the addenda #207 and
/// #208; #208).
///
/// Conventions:
/// - Pupil points p = (px, py) as in rtt-trace (sources.hpp); W in waves at the reference
///   wavelength as opd_points() (opd.hpp), bit for bit the same rays and reference sphere.
/// - All vectors global. k0 is the start direction of a ray, k its direction on the image
///   surface; index c the chief ray of the REFERENCE wavelength (addendum #208 to ADR 0034: both
///   chief frames come from it, like the frame of the PSF plane of ADR 0032, so that Jones pupils
///   of all wavelengths share their bases).
/// - J(p) = exp(i phi) B_out^T P B_in with phi = -2 pi W lambda_ref / lambda (Lam, Eq. (1.2),
///   fields ~ exp(i (k . r - omega t)), more OPL is more phase): a complex 2 x 2 matrix, rows in
///   the output basis (x', y'), columns in the input basis (x, y), row-major j[2 row + col].
///   P is power-normalised (ADR 0021), so J is too (the default of ADR 0034, point 5).
/// - Polarization states of a basis (a, b) as e = a e_1 + b e_2 have the orientation
///   psi = arg(S1 + i S2) / 2 in [0, pi) from e_1 towards e_2 and the ellipticity
///   chi = asin(S3 / S0) / 2 in [-pi/4, pi/4], with the Stokes parameters of the project
///   (docs/architecture.md, "Händigkeit und Stokes"; rtt/polar/ideal.hpp): chi > 0 is
///   right-handed.

#include <array>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

#include "rtt/analysis/opd.hpp"
#include "rtt/analysis/spot.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/math/types.hpp"
#include "rtt/model/validate.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/run_control.hpp"
#include "rtt/trace/sources.hpp"

namespace rtt::analysis {

/// Output basis of the Jones pupil (ADR 0034, point 5). The input basis is always the
/// double-pole basis about the chief ray.
enum class JonesBasis : std::uint8_t {
  /// Double-pole basis about the chief ray (Lam, Appendix B, pp. 279-280): B(k) =
  /// R(k_c -> k) (x'_c, y'_c); singular only at k = -k_c.
  DoublePole,
  /// Dipole basis with the poles at k parallel to +-y'_c (Lam, Appendix B, pp. 278-279, after the
  /// text, not the printed formulas): x' = y'_c x k / |y'_c x k|, y' = k x x'.
  Dipole,
};

/// Options of the Jones pupil.
struct JonesPupilOptions {
  trace::PupilSampling sampling = trace::GridPupil{33};  ///< pupil points, as opd_points()
  JonesBasis basis = JonesBasis::DoublePole;             ///< output basis
  /// Multiply every J by exp(-i psi), psi = arg(det J(0, 0)) / 2, the mean of the eigenphases of
  /// the chief ray at the evaluated wavelength (ADR 0034, point 5, decision F2). For
  /// |det J(0, 0)| <= 1e-12 ||J(0, 0)||_F^2 (e.g. an ideal polarizer at the chief ray) psi is the
  /// argument of the element of J(0, 0) with the largest modulus.
  bool remove_global_phase = false;
  trace::Aiming aiming = trace::Aiming::Real;  ///< aiming of all rays, chief rays included
  /// Warning "rays.lost" if more than this fraction of the launched rays is lost; in [0, 1].
  double lost_warning_fraction = 0.5;
};

/// Orientation and ellipticity of a polarization state in a basis (file comment).
struct PolarizationEllipse {
  double orientation = 0.0;  ///< psi in [0, pi), rad, from e_1 towards e_2
  double ellipticity = 0.0;  ///< chi in [-pi/4, pi/4], rad, > 0 right-handed
};

/// Jones pupil of one field at one wavelength (ADR 0034, points 5 and 6). Per pupil point; a
/// point whose ray did not arrive on the image surface or misses the reference sphere carries
/// the status of opd_points() and NaN in j and in every map, a point near a pole of a basis
/// (closer than 1e-6 rad) NaN in j and the maps with its status Alive, with the warning
/// "pupil.basis_pole". In the convenience form the rays stay close to the chief ray and no pole
/// is reached; with given start rays it can be.
struct JonesPupil {
  compile::PathId path{};
  std::uint16_t field = 0;               ///< index into CompiledSystem::fields().points
  std::uint16_t wavelength = 0;          ///< index into CompiledSystem::wavelengths_um()
  double wavelength_um = 0.0;            ///< evaluated wavelength in vacuum, um
  double reference_wavelength_um = 0.0;  ///< reference wavelength in vacuum, um (unit of W)
  JonesBasis basis = JonesBasis::DoublePole;
  bool remove_global_phase = false;
  double global_phase = 0.0;  ///< psi removed, rad; 0 without remove_global_phase
  /// Chief frame at the input, global: x_c, y_c = R(z -> k0_c) (x, y) (R = diag(1, -1, -1) for
  /// k0_c = -z), k0_c the start direction of the chief ray of the reference wavelength.
  math::Vec3 x_in = math::Vec3::UnitX();
  math::Vec3 y_in = math::Vec3::UnitY();
  math::Vec3 k_in = math::Vec3::UnitZ();  ///< k0_c
  /// Chief frame at the output, global, the frame of the PSF plane (ADR 0032, point 4):
  /// x'_c = (x_B - (x_B . k_c) k_c) / |...|, y'_c = k_c x x'_c, x_B the local x axis of the image
  /// surface, k_c the direction of the chief ray of the reference wavelength on it.
  math::Vec3 x_out = math::Vec3::UnitX();
  math::Vec3 y_out = math::Vec3::UnitY();
  math::Vec3 k_out = math::Vec3::UnitZ();  ///< k_c
  /// Of W, as opd_points(); in the form for given start rays the chief-ray point and radius NaN
  /// (no sphere there).
  ReferenceSphere sphere;
  // Per pupil point, in the order of the sampling:
  std::vector<double> px, py;
  std::vector<trace::RayStatus> status;  ///< as OpdPoint::status
  std::vector<double> w;                 ///< W in waves at lambda_ref; 0 if not arrived
  /// J, 4 complex values per point, row-major (2 row + col); NaN where not defined (above).
  std::vector<std::complex<double>> j;
  /// Transmission T = weight: the power for an unpolarized source with all factors (ADR 0021,
  /// point 2); NaN where not defined.
  std::vector<double> transmission;
  /// Diattenuation D = (s1^2 - s2^2) / (s1^2 + s2^2) from the singular values s1 >= s2 of J
  /// (Lam, Eq. (4.3)), and the input state of largest transmission in the input basis.
  std::vector<double> diattenuation;
  std::vector<double> diattenuation_orientation;  ///< psi of that state, rad
  std::vector<double> diattenuation_ellipticity;  ///< chi of that state, rad
  /// Physical retardance delta in [0, pi] from M = Q^T P (ADR 0034, point 4;
  /// rtt::polar::physical_retardance) and its fast eigenstate in the input basis.
  std::vector<double> retardance;
  std::vector<double> retardance_orientation;  ///< psi of the fast state, rad
  std::vector<double> retardance_ellipticity;  ///< chi of the fast state, rad
  RayLosses losses;  ///< every launched ray by the status of its point (ADR 0023)
  /// Warnings with stable codes: "rays.lost", "stop.clips_beam" as opd_points(),
  /// "pupil.basis_pole" (location /paths/i, the number of points in the message).
  std::vector<model::Diagnostic> warnings;
};

/// Jones pupil of `field` at `wavelength` on options.sampling (ADR 0034, points 5 and 6). The
/// rays are those of opd_points() with the same sampling and aiming, traced with Q requested
/// (RayBatch::enable_q), so W and status are bit for bit those of opd_points(); bit for bit at
/// every thread count. A point closer than 1e-6 rad to a pole of its input or output basis gets
/// J = NaN and NaN maps (not reached in this form, see JonesPupil).
/// @throws std::invalid_argument for an invalid path or wavelength, lost_warning_fraction outside
///         [0, 1], and as trace::make_rays (e.g. an invalid field or sampling)
/// @throws AnalysisError with location /paths/i for a path with an ordinary or extraordinary
///         event (crystal modes, ADR 0034, H6), before any ray is traced; with the location of
///         the image surface if the chief ray grazes it (|x_B x k_c| <= 1e-6); and as
///         opd_points() (chief ray lost or missing the reference sphere)
/// @throws rtt::paraxial::ParaxialError, rtt::compile::NoStopError as opd_points() (e.g. for a
///         path that is not rotationally symmetric)
[[nodiscard]] JonesPupil jones_pupil(const compile::CompiledSystem& system,
                                     compile::PathId path,
                                     std::uint16_t field,
                                     std::uint16_t wavelength,
                                     const JonesPupilOptions& options = {});

/// jones_pupil() with cancellation and progress (#83), stages "aim" and "trace"; with an empty
/// control the same result bit for bit.
[[nodiscard]] JonesPupil jones_pupil(const compile::CompiledSystem& system,
                                     compile::PathId path,
                                     std::uint16_t field,
                                     std::uint16_t wavelength,
                                     const JonesPupilOptions& options,
                                     const trace::RunControl& control);

/// Jones pupil of given start rays (ADR 0034, addendum #208, PR C): the form for paths that are
/// not rotationally symmetric (folds, tilts), without aiming and without paraxial data. The rays
/// of `start` (a copy, traced with Q requested) need one wavelength; `chief` is the index of the
/// chief ray among them, which gives both chief frames (its start direction and its direction on
/// the image surface, the surface of the last event) and the reference of the phase. px, py are
/// taken from start.pupil_x() and start.pupil_y(), `field` from the chief ray's field index. The
/// rays must come from a source, with P = I (and Q = I if enabled), as a fresh batch: P of a
/// batch traced before maps the transverse plane of its original start direction, while B_in is
/// formed at the current one, so J would mix two planes (not checked). Of the options, basis,
/// remove_global_phase and lost_warning_fraction apply; sampling and aiming do not.
///
/// The phase has no reference sphere: phi = 2 pi (OPL - OPL_c) / lambda with the OPL up to the
/// image surface, so w = (OPL_c - OPL) / lambda_ref in waves at lambda_ref (the sign of W) and
/// phi = -2 pi w lambda_ref / lambda as in the convenience form. This is exact for a collimated
/// beam on an image surface perpendicular to it; for a converging beam it is the OPD against the
/// image plane, not against a sphere. `sphere.radius` is NaN, `sphere.centre` the chief-ray
/// point on the image surface. A point arrives (status Alive) if its ray ends Alive on the image
/// surface (Alive elsewhere counts as Vignetted); a point that does not has w = 0 (as opd_points())
/// and J = NaN. A point within 1e-6 rad of a pole of a basis gets J = NaN and NaN maps, with the
/// warning "pupil.basis_pole" (location /paths/i, the number of points in the message).
/// @throws std::invalid_argument for an invalid path, an empty batch, chief >= start.size(), rays
///         of different wavelengths or an invalid wavelength index, lost_warning_fraction
///         outside [0, 1]
/// @throws AnalysisError with location /paths/i for a path with crystal modes (before tracing);
///         with the lost chief ray if it does not end Alive on the image surface; with the
///         location of the image surface if the chief ray grazes it (|x_B x k_c| <= 1e-6)
[[nodiscard]] JonesPupil jones_pupil(const compile::CompiledSystem& system,
                                     compile::PathId path,
                                     const trace::RayBatch& start,
                                     std::size_t chief,
                                     const JonesPupilOptions& options = {});

/// The form for given start rays with cancellation and progress (#83), stage "trace".
[[nodiscard]] JonesPupil jones_pupil(const compile::CompiledSystem& system,
                                     compile::PathId path,
                                     const trace::RayBatch& start,
                                     std::size_t chief,
                                     const JonesPupilOptions& options,
                                     const trace::RunControl& control);

/// Chief frame at the input: (x_c, y_c) = R(z -> k0_c) (x, y), global unit vectors; for
/// k0_c = -z the rotation by pi about x, R = diag(1, -1, -1) (ADR 0034, point 5). Near -z the
/// frame is not continuous: R(z -> k0_c) turns by almost pi about an axis that follows k0_c, so a
/// chief ray tilted towards y ends near (x, -y), one tilted towards x near (-x, y), which turns J
/// into -J. D, delta and psi (modulo pi) do not change; only systems whose rays start along -z
/// are affected (review of #208).
/// @param k0_c unit start direction of the chief ray, global
[[nodiscard]] std::array<math::Vec3, 2> input_chief_frame(const math::Vec3& k0_c) noexcept;

/// Double-pole basis at the unit direction k about the chief frame (x_c, y_c, k_c): the columns
/// R(k_c -> k) x_c and R(k_c -> k) y_c (Lam, Appendix B, pp. 279-280; R the smallest rotation,
/// rtt::trace::rotation_between). None within 1e-6 rad of the pole k = -k_c.
[[nodiscard]] std::optional<std::array<math::Vec3, 2>> double_pole_basis(const math::Vec3& x_c,
                                                                         const math::Vec3& y_c,
                                                                         const math::Vec3& k_c,
                                                                         const math::Vec3& k);

/// Dipole basis at the unit direction k with the poles at +-y_c (Lam, Appendix B, pp. 278-279,
/// after the text): x = y_c x k / |y_c x k|, y = k x x. None within 1e-6 rad of a pole.
[[nodiscard]] std::optional<std::array<math::Vec3, 2>> dipole_basis(const math::Vec3& y_c,
                                                                    const math::Vec3& k);

/// Orientation and ellipticity of the state e = a e_1 + b e_2 (file comment); for a = b = 0 both
/// are NaN.
[[nodiscard]] PolarizationEllipse polarization_ellipse(std::complex<double> a,
                                                       std::complex<double> b) noexcept;

}  // namespace rtt::analysis
