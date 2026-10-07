#pragma once

/// @file field.hpp
/// Field-dependent analyses: distortion and sagittal/tangential field curvature, as data
/// objects (docs/architecture.md, Analyse; decided for #31).
///
/// Conventions:
/// - Field sweep: relative field fraction f in [0, 1] along +y of the largest field point of the
///   system (largest hypot(x, y) of its field values); the field value is (0, f * max) in the
///   units of the system's field type (degree for angles, linear in the angle; mm for heights).
/// - Image surface as in spot.hpp: the surface of the last path event.
/// - Heights are signed: the offset from the image-surface vertex projected on the field
///   direction (x, y) / |(x, y)| in the global x-y plane (+y on axis), so an inverted image has
///   negative heights. The path is rotationally symmetric about z (required by the aiming).
/// - Distortion D = (h_real - h_par) / h_par * 100 %. h_real is the real chief ray (pupil centre)
///   on the image surface; h_par is the paraxial chief ray (rtt-paraxial) evaluated in the plane
///   of the image-surface vertex, NOT in the paraxial image plane (some programs use the
///   latter). D = 0 on axis. Real and paraxial chief ray start at the same field point: the field
///   value is converted into a direction or an object point at the reference wavelength (as in
///   rtt-trace, #50), and both rays are traced at `wavelength`.
/// - Field curvature, numerically with neighbour rays: for a field, two rays at normalised pupil
///   coordinates +-delta along the tangential direction (towards the field point; +y on axis)
///   and two along the sagittal direction are aimed and traced to the image surface. The focus
///   is the midpoint of the shortest connection of the two lines in image space. The result is
///   the z component of its offset from the image-surface vertex (parallel to the axis), with
///   the sign of the image-space propagation (positive = further along the light), mm.
///   The symmetric pair cancels terms odd in delta (e.g. coma); the error is O(delta^2).

#include <cstdint>
#include <vector>

#include "rtt/analysis/spot.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/model/model.hpp"
#include "rtt/trace/sources.hpp"

namespace rtt::analysis {

/// Options of a field sweep.
struct FieldSweepOptions {
  int samples = 11;  ///< field fractions i / (samples - 1), i = 0 .. samples - 1; >= 2
  trace::Aiming aiming = trace::Aiming::Real;  ///< aiming of all rays
};

/// Distortion at one field value.
struct DistortionPoint {
  double fraction = 0.0;         ///< relative field (sweep) or hypot(value) / max (field_at)
  model::Field field;            ///< field value in the units of the system's field type
  double real_height = 0.0;      ///< signed real chief-ray height on the image surface, mm
  double paraxial_height = 0.0;  ///< signed paraxial chief-ray height at the vertex plane, mm
  double percent = 0.0;          ///< D in percent; 0 on axis
};

/// Distortion over the field sweep at `wavelength`.
/// @throws std::invalid_argument for an invalid path or wavelength, samples < 2, or a system
///         without off-axis field point (and as rtt::trace::aim_ray)
/// @throws rtt::compile::NoStopError (a std::invalid_argument) if the path has no stop; checked
///         by rtt::compile::require_stop before any ray is traced (ADR 0022)
/// @throws rtt::paraxial::ParaxialError if the path is not rotationally symmetric
/// @throws AnalysisError if a chief ray does not reach the image surface
[[nodiscard]] std::vector<DistortionPoint> distortion(const compile::CompiledSystem& system,
                                                      compile::PathId path,
                                                      std::uint16_t wavelength,
                                                      const FieldSweepOptions& options = {});

/// Distortion at one field value; throws as distortion().
[[nodiscard]] DistortionPoint distortion_at(const compile::CompiledSystem& system,
                                            compile::PathId path,
                                            const model::Field& field,
                                            std::uint16_t wavelength,
                                            trace::Aiming aiming = trace::Aiming::Real);

/// Options of the field-curvature analysis.
struct FieldCurvatureOptions {
  int samples = 11;     ///< sweep as FieldSweepOptions::samples
  double delta = 1e-3;  ///< half separation of the neighbour rays, normalised pupil units
  trace::Aiming aiming = trace::Aiming::Real;  ///< aiming of all rays
};

/// Sagittal and tangential focus at one field value.
struct FieldCurvaturePoint {
  double fraction = 0.0;     ///< relative field as in DistortionPoint
  model::Field field;        ///< field value in the units of the system's field type
  double tangential = 0.0;   ///< tangential focus from the image-surface vertex, mm
  double sagittal = 0.0;     ///< sagittal focus from the image-surface vertex, mm
  double astigmatism = 0.0;  ///< tangential - sagittal, mm
};

/// Field curvature over the field sweep at `wavelength`.
/// @throws std::invalid_argument for an invalid path or wavelength, samples < 2, delta not in
///         (0, 1), or a system without off-axis field point (and as rtt::trace::aim_ray)
/// @throws rtt::compile::NoStopError (a std::invalid_argument) if the path has no stop; checked
///         by rtt::compile::require_stop before any ray is traced (ADR 0022)
/// @throws rtt::paraxial::ParaxialError if the path is not rotationally symmetric
/// @throws AnalysisError if a ray does not reach the image surface or the neighbour rays are
///         parallel in image space (afocal)
[[nodiscard]] std::vector<FieldCurvaturePoint> field_curvature(
    const compile::CompiledSystem& system,
    compile::PathId path,
    std::uint16_t wavelength,
    const FieldCurvatureOptions& options = {});

/// Field curvature at one field value; throws as field_curvature().
[[nodiscard]] FieldCurvaturePoint field_curvature_at(const compile::CompiledSystem& system,
                                                     compile::PathId path,
                                                     const model::Field& field,
                                                     std::uint16_t wavelength,
                                                     const FieldCurvatureOptions& options = {});

}  // namespace rtt::analysis
