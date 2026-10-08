#pragma once

/// @file paths.hpp
/// Evaluation of explicit ray paths (#122, M4): transmission per path and the optical path
/// difference of two paths for the same start rays. Basis of the Michelson acceptance; no
/// interferogram.
///
/// Conventions (decided for #122):
/// - Main form: the caller gives the start rays (trace::RayBatch). They are copied and traced
///   on each path with the SequentialTracer. This works for folded and tilted paths (beam
///   splitters, interferometers), which the paraxial aiming of trace::make_rays does not
///   accept. Convenience form: field, wavelength, pupil sampling and aiming through make_rays,
///   only for rotationally symmetric paths.
/// - Start rays: every ray of the batch counts as launched. A ray that is not Alive at the start
///   (e.g. NoConvergence from make_rays) is not traced and counts as lost with its status
///   (weight 0). Any system wavelength index is allowed, also mixed in one batch. The start
///   weight (e.g. a source apodization) and the start OPL are carried along.
/// - A ray arrived on a path if it is Alive after the trace and its last surface is the surface
///   of the path's last event (the image surface); every other ray is lost (ADR 0023).
/// - Transmission: weight is the power for an unpolarized source (ADR 0021), so the mean of the
///   final weights over all launched rays (lost rays 0) is the transmitted power fraction of a
///   uniformly illuminated pupil. No wavelength weights: rays count equally.
/// - Optical path difference in mm, not in waves: per start ray, the OPL on path b minus the
///   OPL on path a at the image surface, where the two rays may land at different points (in a
///   Michelson they coincide). A path length, not a wavefront; waves and the superposition at
///   one point belong to the interferogram (not part of #122). It is only defined where both
///   paths end on the same surface, so two paths with different image surfaces are an error.

#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

#include "rtt/analysis/spot.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/model/validate.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/run_control.hpp"
#include "rtt/trace/sources.hpp"

namespace rtt::analysis {

/// Options of the path evaluations.
struct PathOptions {
  /// Convenience form only: pupil sampling of the rays from trace::make_rays.
  trace::PupilSampling sampling = trace::HexapolarPupil{6};
  trace::Aiming aiming = trace::Aiming::Real;  ///< convenience form only: aiming of the rays
  /// Warning "rays.lost" if more than this fraction of the launched rays is lost on a path
  /// (ADR 0023); in [0, 1].
  double lost_warning_fraction = 0.5;
};

/// One launched ray of a path transmission.
struct PathRay {
  double px = 0.0;      ///< normalised pupil coordinate of the start ray (RayBatch::pupil_x)
  double py = 0.0;      ///< normalised pupil coordinate of the start ray (RayBatch::pupil_y)
  double weight = 0.0;  ///< final weight if the ray arrived, otherwise 0; dimensionless
  /// Alive if the ray arrived; otherwise its status (Vignetted if it stopped on another
  /// surface while Alive).
  trace::RayStatus status = trace::RayStatus::Alive;
};

/// Transmission of one path for a set of start rays.
struct PathTransmission {
  compile::PathId path{0};          ///< the evaluated path
  std::uint32_t image_surface = 0;  ///< surface of the path's last event, index into surfaces()
  std::vector<PathRay> rays;        ///< one entry per launched ray, in batch order
  std::size_t rays_launched = 0;    ///< rays of the start batch (= rays.size())
  std::size_t rays_arrived = 0;     ///< rays that reached the image surface Alive
  /// Mean final weight over all launched rays, lost rays counted as 0: the transmitted power
  /// fraction for a uniformly illuminated pupil and an unpolarized source. With start weights
  /// other than 1 it is the mean of the final weights, i.e. the apodized transmitted power, not
  /// a ratio to the start weights. Dimensionless.
  double mean = 0.0;
  double min = 0.0;  ///< smallest final weight of the arrived rays (0 if none), dimensionless
  double max = 0.0;  ///< largest final weight of the arrived rays (0 if none), dimensionless
  RayLosses losses;  ///< launched rays by final status, worst loss surface (ADR 0023)
  /// Warnings with stable codes (ADR 0022, 0023): "rays.lost" above lost_warning_fraction,
  /// "stop.clips_beam" if rays end Vignetted at the stop surface of the path.
  std::vector<model::Diagnostic> warnings;
};

/// Optical path difference of one start ray on two paths.
struct OplDifferencePoint {
  double px = 0.0;     ///< normalised pupil coordinate of the start ray
  double py = 0.0;     ///< normalised pupil coordinate of the start ray
  double delta = 0.0;  ///< OPL on path b minus OPL on path a at the image surface, mm; 0 if lost
  /// Alive if the ray arrived on both paths; otherwise its status on path a if it was lost
  /// there, else its status on path b.
  trace::RayStatus status = trace::RayStatus::Alive;
};

/// Optical path difference OPL_b - OPL_a of two paths with the same image surface.
struct PathOplDifference {
  compile::PathId path_a{0};               ///< reference path (subtracted)
  compile::PathId path_b{0};               ///< path whose OPL is taken positive
  std::uint32_t image_surface = 0;         ///< common surface of the last events
  std::vector<OplDifferencePoint> points;  ///< one entry per launched ray, in batch order
  /// delta of the first launched ray at pupil (0, 0) that arrived on both paths, mm; none if
  /// there is no such ray.
  std::optional<double> chief;
  RayLosses losses_a;                       ///< launched rays on path a by final status (ADR 0023)
  RayLosses losses_b;                       ///< launched rays on path b by final status (ADR 0023)
  std::vector<model::Diagnostic> warnings;  ///< of both paths, path a first (ADR 0022, 0023)
};

/// Transmission of `path` for the start rays `start` (main form, see file comment).
/// @param start start rays in global coordinates; wl must be system wavelength indices
/// @throws std::invalid_argument for an invalid path, an empty `start`, a wavelength index
///         that is not a system wavelength or an invalid status in `start` (as
///         trace::SequentialTracer::trace), or options.lost_warning_fraction outside [0, 1]
[[nodiscard]] PathTransmission path_transmission(const compile::CompiledSystem& system,
                                                 compile::PathId path,
                                                 const trace::RayBatch& start,
                                                 const PathOptions& options = {});

/// Transmission of `path` for the rays of trace::make_rays (convenience form): `field` at
/// `wavelength` with options.sampling and options.aiming.
/// @throws as the main form, and std::invalid_argument for an invalid field or wavelength
///         index (and as rtt::trace::make_rays); rtt::compile::NoStopError (a
///         std::invalid_argument) if the path has no stop; rtt::paraxial::ParaxialError if the
///         path is not rotationally symmetric (the paraxial aiming of make_rays; use the main
///         form for folded or tilted paths). Options and paths are checked before the aiming.
[[nodiscard]] PathTransmission path_transmission(const compile::CompiledSystem& system,
                                                 compile::PathId path,
                                                 std::uint16_t field,
                                                 std::uint16_t wavelength,
                                                 const PathOptions& options = {});

/// Optical path difference OPL_b - OPL_a for the start rays `start` (main form).
/// @throws std::invalid_argument if the last events of the two paths are on different surfaces,
///         and as path_transmission()
[[nodiscard]] PathOplDifference opl_difference(const compile::CompiledSystem& system,
                                               compile::PathId path_a,
                                               compile::PathId path_b,
                                               const trace::RayBatch& start,
                                               const PathOptions& options = {});

/// opl_difference() for the rays of trace::make_rays on path a (convenience form); the same
/// start rays are traced on path b.
/// @throws as opl_difference() and the convenience form of path_transmission()
[[nodiscard]] PathOplDifference opl_difference(const compile::CompiledSystem& system,
                                               compile::PathId path_a,
                                               compile::PathId path_b,
                                               std::uint16_t field,
                                               std::uint16_t wavelength,
                                               const PathOptions& options = {});

/// Overloads with cancellation and progress (#83, rtt/trace/run_control.hpp): run-time control,
/// never part of the options. Stages "aim" (convenience form) and "trace" per path. With an
/// empty `control` they are the functions above exactly; the results never depend on it.
/// @throws as the functions above, trace::Cancelled after a cancellation, or the exception of
///         the progress callback
[[nodiscard]] PathTransmission path_transmission(const compile::CompiledSystem& system,
                                                 compile::PathId path,
                                                 const trace::RayBatch& start,
                                                 const PathOptions& options,
                                                 const trace::RunControl& control);
[[nodiscard]] PathTransmission path_transmission(const compile::CompiledSystem& system,
                                                 compile::PathId path,
                                                 std::uint16_t field,
                                                 std::uint16_t wavelength,
                                                 const PathOptions& options,
                                                 const trace::RunControl& control);
[[nodiscard]] PathOplDifference opl_difference(const compile::CompiledSystem& system,
                                               compile::PathId path_a,
                                               compile::PathId path_b,
                                               const trace::RayBatch& start,
                                               const PathOptions& options,
                                               const trace::RunControl& control);
[[nodiscard]] PathOplDifference opl_difference(const compile::CompiledSystem& system,
                                               compile::PathId path_a,
                                               compile::PathId path_b,
                                               std::uint16_t field,
                                               std::uint16_t wavelength,
                                               const PathOptions& options,
                                               const trace::RunControl& control);

}  // namespace rtt::analysis
