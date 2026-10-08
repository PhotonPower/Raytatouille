#pragma once

/// @file opd.hpp
/// Wavefront error (optical path difference, OPD) against a reference sphere, as data objects
/// (docs/architecture.md, Analyse; decided for #29).
///
/// Conventions:
/// - Sign (J. C. Wyant, K. Creath, "Basic Wavefront Aberration Theory for Optical Metrology",
///   Applied Optics and Optical Engineering XI (1992), Sec. I): W = OPL_ref - OPL_ray, both
///   measured up to the reference sphere. W > 0 means the wavefront leads the reference, i.e.
///   curves in more than the reference sphere.
/// - Reference sphere: centred on the real hit point of the chief ray (pupil centre) of the
///   reference wavelength on the image surface (the surface of the last path event); radius =
///   distance from that point to the centre of the paraxial exit pupil (rtt-paraxial). A ray's
///   OPL to the sphere is its OPL at the image surface minus |n'| s, with s the signed path from
///   the sphere to the image surface along the ray and n' the image-space index (real part).
///   It is evaluated as OPL - |n'| (s - R) on the crossing towards the exit pupil, free of
///   cancellation; the constant |n'| R is common to all rays and cancels in W (#102).
/// - Exit pupil at infinity (image-space telecentric, #102): the limit R -> infinity, radius
///   +infinity; W = OPL_chief - OPL_ray + |n'| d . (p - C) with p the hit, d the direction of
///   the ray and C the centre, i.e. the OPL to the foot of the perpendicular from C on the ray.
///   Continuous with a finite but distant exit pupil from both sides.
/// - Reference ray: the chief ray of the same wavelength, so W(0, 0) = 0 at every wavelength.
///   Piston and tilt are not removed; lateral colour shows up as tilt.
/// - Unit: waves at the reference wavelength, W = OPD[mm] / lambda_ref[mm], for every
///   wavelength (docs/architecture.md, Größen).
/// - Rays that do not arrive Alive on the image surface keep W = 0 and their status; they do
///   not count for RMS and PV.

#include <cstddef>
#include <cstdint>
#include <vector>

#include "rtt/analysis/spot.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/math/types.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/sources.hpp"

namespace rtt::analysis {

/// Reference sphere in global coordinates.
struct ReferenceSphere {
  math::Vec3 centre = math::Vec3::Zero();  ///< chief-ray point on the image surface, mm
  /// Distance to the exit-pupil centre, mm; +infinity if the exit pupil is at infinity.
  double radius = 0.0;
};

/// OPD at one pupil point.
struct OpdPoint {
  double px = 0.0;  ///< normalised pupil coordinate x
  double py = 0.0;  ///< normalised pupil coordinate y
  double w = 0.0;   ///< OPD in waves at the reference wavelength; 0 if the ray did not arrive
  /// Alive if the ray arrived on the image surface, otherwise its status; Vignetted also if it
  /// was Alive but stopped elsewhere or its line misses the reference sphere. Consumers must
  /// check it before using w.
  trace::RayStatus status = trace::RayStatus::Alive;
};

/// Options of OPD analyses.
struct OpdOptions {
  int grid = 33;        ///< map: GridPupil{grid}, points inside the unit circle
  int fan_points = 21;  ///< fans: points evenly spaced on [-1, 1]
  trace::Aiming aiming = trace::Aiming::Real;  ///< aiming of all rays, chief rays included
  /// Warning "rays.lost" if more than this fraction of the launched rays is lost (ADR 0023);
  /// in [0, 1]. Vignetting at the field edge is intended, hence the default of one half.
  double lost_warning_fraction = 0.5;
};

/// OPD map over the pupil.
struct OpdMap {
  std::uint16_t field = 0;       ///< index into CompiledSystem::fields().points
  std::uint16_t wavelength = 0;  ///< index into CompiledSystem::wavelengths_um()
  ReferenceSphere sphere;
  std::vector<OpdPoint> points;  ///< GridPupil order (rows from py = -1, px fastest)
  /// Standard deviation of W over the arrived points (piston removed, tilt NOT removed,
  /// unweighted over the uniform grid), waves.
  double rms = 0.0;
  double pv = 0.0;            ///< max W - min W over the arrived points, waves
  std::size_t arrived = 0;    ///< points with status Alive (they define rms and pv)
  std::size_t vignetted = 0;  ///< points that did not arrive (arrived + vignetted = size)
  /// Rays of the sampling by the final status of their points (a ray that misses the
  /// reference sphere counts as Vignetted), worst loss surface (ADR 0023).
  RayLosses losses;
  /// Warnings with stable codes (ADR 0022, 0023): "rays.lost" above lost_warning_fraction,
  /// "stop.clips_beam" if rays end Vignetted at the stop surface.
  std::vector<model::Diagnostic> warnings;
};

/// Tangential (px = 0) and sagittal (py = 0) OPD fans.
struct OpdFan {
  std::uint16_t field = 0;       ///< index into CompiledSystem::fields().points
  std::uint16_t wavelength = 0;  ///< index into CompiledSystem::wavelengths_um()
  ReferenceSphere sphere;
  std::vector<OpdPoint> tangential;
  std::vector<OpdPoint> sagittal;
  /// Rays of the sampling by the final status of their points (a ray that misses the
  /// reference sphere counts as Vignetted), worst loss surface (ADR 0023).
  RayLosses losses;
  /// Warnings with stable codes (ADR 0022, 0023): "rays.lost" above lost_warning_fraction,
  /// "stop.clips_beam" if rays end Vignetted at the stop surface.
  std::vector<model::Diagnostic> warnings;
};

/// OPD map of `field` at `wavelength`.
/// @throws std::invalid_argument for an invalid path, field, wavelength, grid < 1 or
///         options.lost_warning_fraction outside [0, 1] (and as rtt::trace::make_rays)
/// @throws rtt::paraxial::ParaxialError if the path is not rotationally symmetric or the stop
///         aperture is not circular (aiming and exit pupil)
/// @throws rtt::compile::NoStopError (a std::invalid_argument) if the path has no stop; checked
///         by rtt::compile::require_stop before any ray is traced (ADR 0022)
/// @throws AnalysisError if a chief ray does not reach the image surface or misses the
///         reference sphere, or no ray arrives (an exit pupil at infinity is supported, #102)
[[nodiscard]] OpdMap opd_map(const compile::CompiledSystem& system,
                             compile::PathId path,
                             std::uint16_t field,
                             std::uint16_t wavelength,
                             const OpdOptions& options = {});

/// OPD fans of `field` at `wavelength`; conventions as opd_map().
/// @throws as opd_map(), std::invalid_argument for fan_points < 1
[[nodiscard]] OpdFan opd_fan(const compile::CompiledSystem& system,
                             compile::PathId path,
                             std::uint16_t field,
                             std::uint16_t wavelength,
                             const OpdOptions& options = {});

}  // namespace rtt::analysis
