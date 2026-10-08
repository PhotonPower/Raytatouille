#pragma once

/// @file ghosts.hpp
/// Ghost ranking (#124, ADR 0027 addendum): every ghost of a base path, traced with the start
/// rays of the base aiming, ranked by its irradiance at the image relative to the useful image.
///
/// Metric (ADR 0027, addendum #124): E = P / (pi (r^2 + r0^2)) with P the transmitted power
/// (mean final weight over the launched rays, #122), r the weighted RMS radius of the arrived
/// rays about their centroid on the image surface and r0 the resolution radius of the detector;
/// the rank value is rho = E_ghost / E_base = (P_g / P_b) (r_b^2 + r0^2) / (r_g^2 + r0^2).
/// - rho depends on r0: a model choice for a detector, not a physical constant.
/// - The RMS radius is geometric; diffraction is not taken into account (M6).
/// - Ghosts do not warn about lost rays (losses are normal for them); the warnings come from
///   the base path only.
/// - The paraxial focus position and blur radius of a ghost are diagnostics, not rank values.

#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

#include "rtt/analysis/spot.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/compile/ghosts.hpp"
#include "rtt/model/validate.hpp"
#include "rtt/trace/run_control.hpp"
#include "rtt/trace/sources.hpp"

namespace rtt::analysis {

/// Options of the ghost ranking.
struct GhostRankingOptions {
  trace::PupilSampling sampling = trace::HexapolarPupil{6};  ///< start rays of the base path
  trace::Aiming aiming = trace::Aiming::Real;                ///< aiming of the start rays
  /// Resolution radius r0 of the detector in mm (>= 0): spots smaller than this do not get
  /// brighter. A model choice, see the file comment.
  double resolution_radius = 0.005;
  /// Warning "rays.lost" if more than this fraction of the base path's rays is lost (ADR
  /// 0023); in [0, 1].
  double lost_warning_fraction = 0.5;
};

/// One ghost of a ranking.
struct GhostEntry {
  compile::PathId path;          ///< the ghost path
  std::uint32_t surface_j = 0;   ///< first ghost reflection (back), index into surfaces()
  std::uint32_t surface_i = 0;   ///< second ghost reflection (forward), index into surfaces()
  double power = 0.0;            ///< P_g: mean final weight over the launched rays, dimensionless
  double relative_power = 0.0;   ///< P_g / P_b, dimensionless
  double rms_radius = 0.0;       ///< r_g: weighted RMS radius on the image surface, mm (0: none)
  std::size_t rays_arrived = 0;  ///< rays of the ghost that reached the image surface
  /// rho = (P_g / P_b) (r_b^2 + r0^2) / (r_g^2 + r0^2): the rank value, dimensionless.
  double relative_irradiance = 0.0;
  /// Paraxial focus of the ghost minus the global z of the image surface vertex, mm; none if
  /// the ghost leaves collimated. Diagnostic (ADR 0027, addendum #124).
  std::optional<double> focus_offset;
  /// |height| of the paraxial marginal ray of the ghost at the image surface vertex plane, mm.
  /// Diagnostic.
  double paraxial_blur_radius = 0.0;
  RayLosses losses;  ///< launched rays of the ghost by final status (ADR 0023)
};

/// Ghosts of a base path, ranked by relative irradiance.
struct GhostRanking {
  compile::PathId base;             ///< the base path (useful image)
  std::uint16_t field = 0;          ///< field index of the start rays
  std::uint16_t wavelength = 0;     ///< wavelength index of the start rays
  double base_power = 0.0;          ///< P_b, dimensionless
  double base_rms_radius = 0.0;     ///< r_b, mm
  double resolution_radius = 0.0;   ///< r0 used, mm
  std::vector<GhostEntry> entries;  ///< by relative_irradiance descending, ties in ghost order
  std::vector<model::Diagnostic> warnings;  ///< of the base path only (ADR 0023)
};

/// Ranks the ghosts of `ghosts` (compile::compile_with_ghosts) for `field` at `wavelength`.
/// The start rays come from trace::make_rays on the base path with options.sampling and
/// options.aiming; every ghost is traced with the same rays (path_transmission, #122).
/// @param ghosts     system compiled with its ghosts
/// @param field      index into CompiledSystem::fields().points
/// @param wavelength index into CompiledSystem::wavelengths_um()
/// @param options    sampling, aiming, resolution radius r0 in mm, lost-ray warning threshold
/// @throws std::invalid_argument for a system without ghosts, an invalid field or wavelength,
///         resolution_radius < 0 or lost_warning_fraction outside [0, 1] (and as
///         trace::make_rays)
/// @throws rtt::paraxial::ParaxialError if the base path is not rotationally symmetric (aiming
///         and the paraxial diagnostics)
/// @throws AnalysisError if no ray of the base path arrives (no useful image to compare with)
[[nodiscard]] GhostRanking ghost_ranking(const compile::GhostSystem& ghosts,
                                         std::uint16_t field,
                                         std::uint16_t wavelength,
                                         const GhostRankingOptions& options = {});

/// ghost_ranking() with cancellation and progress (#83): stage "aim" for the start rays, then
/// "trace" for the base path and every ghost. With an empty control it is the function above
/// exactly.
/// @throws as above, trace::Cancelled after a cancellation, or the exception of the callback
[[nodiscard]] GhostRanking ghost_ranking(const compile::GhostSystem& ghosts,
                                         std::uint16_t field,
                                         std::uint16_t wavelength,
                                         const GhostRankingOptions& options,
                                         const trace::RunControl& control);

}  // namespace rtt::analysis
