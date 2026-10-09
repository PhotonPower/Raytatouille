#pragma once

/// @file generators.hpp
/// Merit generators rms_spot and rms_wavefront (ADR 0030, point 4 and addendum of #168): one
/// residual (wavefront) or two (spot) per ray of a Gaussian pupil quadrature, so that the sum of
/// their squares is the weighted mean square of the spot or the wavefront and stays smooth.
///
/// Expansion (deterministic): fields, then wavelengths, then rings, then arms (the order of
/// trace::pupil_points for trace::GaussPupil{rings, arms}); fields and wavelengths in the order
/// of the generator's lists, none: all in system order. The spot gives (x, y) per ray.
///
/// Residual of ray k of field f at wavelength l: sqrt(w W_f W_l q_k) times its deviation, with w
/// the weight of the generator, W_f and W_l the model weights of the chosen fields and
/// wavelengths normalised to sum 1 over the choice, q_k the quadrature weight
/// (trace::gauss_pupil_weights, sum 1). The sum of squares is w sum_f W_f RMS_f^2 (spot) and
/// w sum_f sum_l W_f W_l sigma_fl^2 (wavefront).
///
/// Deviations:
/// - rms_spot: (x_k - x_ref, y_k - y_ref) in the local coordinates of the image surface (the
///   surface of the last event), mm, as analysis::spot. Reference "centroid": per field over all
///   chosen wavelengths, weighted with W_l q_k over the arrived rays; "chief": the chief ray
///   (pupil centre) of the field at the reference wavelength.
/// - rms_wavefront: W_k - W_mean in waves at the reference wavelength, W from analysis::opd_points
///   (reference sphere about the chief ray of the reference wavelength, OPL reference the chief
///   ray of the same wavelength; tilt not removed, lateral colour shows as tilt). W_mean is
///   sum q_k W_k / sum q_k over the arrived rays of the field and wavelength (piston only).
///
/// Lost rays (F7): a ray that does not arrive Alive on the image surface (wavefront: or misses
/// the reference sphere) has the residuals 0, does not count for the reference and is counted in
/// GeneratorStats. Lost rays lower the sum of squares (no renormalisation): a known limit of
/// ADR 0030. If no ray of a field (spot) or of a field and wavelength (wavefront) arrives, or the
/// chief ray of a spot with reference "chief" does not, its residuals are NaN: the evaluation is
/// invalid (ADR 0030, point 10). All rays use real aiming.

#include <cstddef>
#include <span>
#include <string>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/model/optimization.hpp"

namespace rtt::optim {

/// Rays of one generator in one evaluation.
struct GeneratorStats {
  std::size_t rays_launched = 0;  ///< |fields| |wavelengths| rings arms
  std::size_t rays_lost = 0;      ///< rays with residuals 0 because they did not arrive
  /// Weighted mean square of the deviations, sum_f sum_l W_f W_l sum_k q_k |d_k|^2 (mm^2 for the
  /// spot, waves^2 for the wavefront), independent of the generator weight; NaN if a residual
  /// is NaN. The sum of squares of the residuals is weight times this.
  double mean_square = 0.0;
  /// Why the residuals are NaN (empty string if they are finite).
  std::string undefined;
};

/// Number of residuals of `generator` in `system`: |fields| |wavelengths| rings arms, times 2
/// for rms_spot. Fixed for a run: the field and wavelength counts do not depend on variables.
/// @throws std::invalid_argument for rings or arms < 1
[[nodiscard]] std::size_t generator_size(const model::Generator& generator,
                                         const compile::CompiledSystem& system);

/// Writes the residuals of `generator` in `system` to `out` in the order above.
/// @param out   generator_size(generator, system) values
/// @param stats rays and mean square of this evaluation (overwritten)
/// @throws std::invalid_argument for an `out` of another size, an unknown path, a field or
///         wavelength index out of range, or chosen weights that sum to 0 (the start check of
///         MeritFunction rejects these before a run)
/// @throws as analysis::opd_points (rms_wavefront: chief ray lost or missing the sphere) and
///         trace::make_rays
void generator_residuals(const model::Generator& generator,
                         const compile::CompiledSystem& system,
                         std::span<double> out,
                         GeneratorStats& stats);

}  // namespace rtt::optim
