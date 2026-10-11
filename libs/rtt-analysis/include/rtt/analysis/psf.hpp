#pragma once

/// @file psf.hpp
/// Diffraction point spread function by FFT (ADR 0033; #205): scalar, incoherent Fraunhofer
/// diffraction in the direction sines of the image space, from the pupil function of ADR 0032.
///
/// Conventions:
/// - PSF plane: perpendicular to the chief ray of the reference wavelength through C' (the
///   centre of the reference sphere after focus_shift), axes e_x, e_y from the local x axis of
///   the image surface (ADR 0032, point 4). Coordinates x in mm in (e_x, e_y).
/// - Cell value of the pupil u_j = f_j A_j exp(-2 pi i W_j lambda_ref / lambda) (ADR 0032),
///   amplitude U(x) ~ sum_j u_j exp(-2 pi i x . M p_j / lambda) (ADR 0033, point 2).
/// - Lattice: the FFT of length L gives U at x(k) = lambda / (L a) M^-T k, k in [-L/2, L/2)^2,
///   a = 2 / N (ADR 0033, point 3).
/// - Normalisation by the exact pupil energy E = sum f A^2 (ADR 0033, point 4): the intensity of a
///   lattice cell is its fraction of the energy, and the sum over the lattice is
///   1 - edge_energy_deficit.
/// - Limits: Fraunhofer/Debye approximation, scalar, incoherent, no obliquity factors, linear
///   mapping p -> s (ADR 0032, point 4); no vector PSF.

#include <array>
#include <cstdint>
#include <optional>
#include <span>
#include <vector>

#include "rtt/analysis/pupil.hpp"
#include "rtt/analysis/spot.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/math/types.hpp"
#include "rtt/model/validate.hpp"
#include "rtt/trace/run_control.hpp"

namespace rtt::analysis {

/// Options of psf().
struct PsfOptions {
  /// Pupil function (ADR 0032): grid N, apodization, focus_shift, aiming, warning thresholds.
  PupilFunctionOptions pupil;
  int padding = 4;  ///< Q >= 2: FFT length L = Q N (ADR 0033, point 3)
  /// Wavelength indices; none: all system wavelengths with their model weights.
  std::optional<std::vector<std::uint16_t>> wavelengths;
  /// Warning psf.aliasing when the energy in the outer L/16 rows and columns exceeds this factor
  /// times that of the same pupil without aberration (W = 0); >= 1.
  double aliasing_factor = 3.0;
  /// Warning psf.pupil_undersampled above this phase step between neighbouring arrived cells,
  /// waves at the evaluated wavelength; > 0.
  double undersampled_waves = 0.5;
};

/// One wavelength of a PSF.
struct PsfWavelength {
  std::uint16_t wavelength = 0;  ///< index into CompiledSystem::wavelengths_um()
  double weight = 0.0;           ///< model weight, normalised to sum 1 over the chosen wavelengths
  double energy = 0.0;           ///< E = sum f A^2 of its pupil function (cells)
  int length = 0;                ///< FFT length L_i (ADR 0033, point 7)
  /// Fraction of E on its own lattice that lies outside the reference field (|k| >= L/2) and
  /// is therefore not in Psf::intensity: > 0 only for a longer wavelength (L_i > L). Like
  /// edge_energy_deficit it is reported, not lost silently.
  double outside_field = 0.0;
};

/// Peak of the lattice, not interpolated.
struct PsfPeak {
  double value = 0.0;  ///< intensity of the cell (fraction of the energy)
  int k1 = 0;          ///< lattice index along the first lattice vector
  int k2 = 0;          ///< lattice index along the second lattice vector
  double x = 0.0;      ///< position along e_x, mm
  double y = 0.0;      ///< position along e_y, mm
};

/// Diffraction PSF of one field.
struct Psf {
  std::uint16_t field = 0;
  int grid = 0;  ///< N
  /// L of the reference lattice: Q N, larger (in steps of 2) only if a shorter wavelength needs
  /// L_i >= 2 N - 1 (ADR 0033, point 7). The reference is the system reference wavelength if
  /// chosen, otherwise the first chosen wavelength.
  int length = 0;
  math::Vec3 centre = math::Vec3::Zero();  ///< C', global, mm
  math::Vec3 e_x = math::Vec3::UnitX();    ///< first axis of the PSF plane, global
  math::Vec3 e_y = math::Vec3::UnitY();    ///< second axis of the PSF plane, global
  /// Lattice matrix B = lambda / (L a) M^-T of the reference wavelength, row-major, mm in
  /// (e_x, e_y): x(k) = (B[0] k1 + B[1] k2, B[2] k1 + B[3] k2).
  std::array<double, 4> lattice{};
  /// L x L intensities, row k2 + L/2, column k1 + L/2 (k from -L/2, integer division); fraction
  /// of the energy per cell. Polychromatic: wavelength i lies on its own lattice of length
  /// L_i = round(L lambda_i / lambda_ref sqrt(|det M_ref| / |det M_i|)) (nearest integer), which
  /// matches the reference within polychromatic_scale_error; a shorter wavelength with L_i < L
  /// adds nothing outside its own field (|k| beyond L_i / 2).
  std::vector<double> intensity;
  /// Strehl ratio: the intensity at C' (k = 0) relative to that of the same pupils without
  /// aberration (Wyant & Creath (64) for the discrete pupil). Monochromatic |sum u|^2 /
  /// |sum f A|^2. Polychromatic the peak ratio of the summed images with the weights of the
  /// intensity, omega_i = w_i E_i / sum w E:
  ///   sum_i omega_i |sum u|_i^2 / (L_i^2 E_i)  /  sum_i omega_i (sum f A)_i^2 / (L_i^2 E_i),
  /// so strehl times the denominator is intensity at k = 0 (review of #233). It is not the mean
  /// of the monochromatic S_i: the ideal peaks scale as 1 / lambda_i^2.
  double strehl = 0.0;
  PsfPeak peak;
  double edge_energy_deficit = 0.0;  ///< 1 - sum f^2 A^2 / sum f A^2, weighted
  /// max ||B_i - B_ref||_F / ||B_ref||_F over the wavelengths (ADR 0033, point 7): at most
  /// 1 / (2 L_i) from the rounding of L_i, plus the difference of the mappings M_i.
  double polychromatic_scale_error = 0.0;
  std::vector<PsfWavelength> wavelengths;
  RayLosses losses;
  /// rays.lost, stop.clips_beam, pupil.mapping_nonlinear (from the pupil functions),
  /// psf.aliasing, psf.pupil_undersampled.
  std::vector<model::Diagnostic> warnings;
};

/// Diffraction PSF of `field` on `path` (ADR 0033).
/// @throws std::invalid_argument for padding < 2, L > 8192, an unknown or repeated wavelength,
///         an empty wavelength list, aliasing_factor < 1 or undersampled_waves <= 0 (and as
///         pupil_function())
/// @throws AnalysisError (location /paths/i) for an afocal image space (no mapping), a path with
///         crystal modes or no arrived cell (and as pupil_function())
[[nodiscard]] Psf psf(const compile::CompiledSystem& system,
                      compile::PathId path,
                      std::uint16_t field,
                      const PsfOptions& options = {});

/// psf() with cancellation and progress (run_control.hpp).
[[nodiscard]] Psf psf(const compile::CompiledSystem& system,
                      compile::PathId path,
                      std::uint16_t field,
                      const PsfOptions& options,
                      const trace::RunControl& control);

/// Intensity at points `x` of the PSF plane (mm in e_x, e_y of `pupil.mapping`), monochromatic,
/// by the direct sum over the cells (O(N^2) per point), **in Strehl units**:
/// |sum_j u_j exp(-2 pi i x . M p_j / lambda)|^2 / |sum_j f_j A_j|^2, so 1 is the peak of the
/// same pupil without aberration and the value at x = 0 is the Strehl ratio. p_j is the centre of
/// cell j, the lattice of the FFT (ADR 0033, point 3), also for a boundary cell whose W and A were
/// sampled at its centroid (ADR 0032, point 1). The lattice values
/// of psf() are these values times |sum f A|^2 / (L^2 E). lambda and lambda_ref are
/// pupil.wavelength_um and pupil.reference_wavelength_um (vacuum, um), so no system is needed.
/// @throws AnalysisError if the pupil function has no mapping (afocal) or no arrived cell
[[nodiscard]] std::vector<double> psf_at(const PupilFunction& pupil, std::span<const Point2> x);

/// Strehl ratio through focus at one wavelength (ADR 0033, points 8 and 11).
struct ThroughFocus {
  std::uint16_t field = 0;
  std::uint16_t wavelength = 0;  ///< index into CompiledSystem::wavelengths_um()
  std::vector<double> shifts;    ///< focus shifts, mm along the chief ray (ADR 0032, point 5)
  /// |sum u|^2 / |sum f A|^2 per shift (monochromatic Strehl ratio, as Psf::strehl).
  std::vector<double> strehl;
  /// sum f over the cells that cross the shifted reference sphere / sum f over all cells, per
  /// shift: cells far from the focus can miss a strongly shifted sphere (refocus()).
  std::vector<double> arrived_fraction;
  /// Losses and warnings of the traced pupil function (focus shift 0): rays.lost,
  /// stop.clips_beam, pupil.mapping_nonlinear.
  RayLosses losses;
  std::vector<model::Diagnostic> warnings;
};

/// The Strehl ratio |sum u|^2 / |sum f A|^2 for every focus shift in `shifts` (mm along the chief
/// ray, ADR 0032, point 5), at `wavelength`: one trace, W recomputed per shift with refocus(),
/// without tracing again. options.focus_shift is ignored. The PSF per shift (ADR 0033, point 8,
/// optional there) is left to an addendum; psf() with focus_shift gives it with one trace each.
/// @throws std::invalid_argument for a non-finite shift (before any ray), and as
///         pupil_function()
/// @throws AnalysisError (location /paths/i) for an afocal image space, as psf(), a path with
///         crystal modes, or no arrived cell at a shift
[[nodiscard]] ThroughFocus through_focus(const compile::CompiledSystem& system,
                                         compile::PathId path,
                                         std::uint16_t field,
                                         std::uint16_t wavelength,
                                         std::span<const double> shifts,
                                         const PupilFunctionOptions& options = {});

/// through_focus() with cancellation and progress (run_control.hpp).
[[nodiscard]] ThroughFocus through_focus(const compile::CompiledSystem& system,
                                         compile::PathId path,
                                         std::uint16_t field,
                                         std::uint16_t wavelength,
                                         std::span<const double> shifts,
                                         const PupilFunctionOptions& options,
                                         const trace::RunControl& control);

}  // namespace rtt::analysis
