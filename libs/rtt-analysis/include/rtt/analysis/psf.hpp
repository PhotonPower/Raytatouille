#pragma once

/// @file psf.hpp
/// Diffraction point spread function by FFT (ADR 0033; #205), MTF and encircled energy (#206):
/// scalar, incoherent Fraunhofer diffraction in the direction sines of the image space, from the
/// pupil function of ADR 0032.
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
#include <complex>
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

// ---------------------------------------------------------------------------------------------
// MTF and encircled energy (ADR 0033, points 5-7, with the erratum on the signs; #206).
//
// Both use the autocorrelation of the cell values C(l) = sum_j u_(j+l) conj(u_j) for shifts l in
// cells, computed as the inverse FFT of |FFT(u)|^2 with L = 2 N >= 2 N - 1, so exact for the
// discrete model and independent of PsfOptions::padding. With the kernel of point 2 the
// intensity is I(x) ~ sum_l C(l) exp(-2 pi i nu(l) . x), nu(l) = M a l / lambda (cycles per mm in
// the PSF plane, components in e_x, e_y).
// ---------------------------------------------------------------------------------------------

/// Optional outputs of mtf().
struct MtfOutputs {
  /// 2D MTF over all shifts l in [-(N-1), N-1]^2 (H8); only for a single wavelength.
  bool map = false;
  /// Complex OTF next to the modulation (phase transfer).
  bool complex_otf = false;
};

/// Options of mtf().
struct MtfOptions {
  /// Pupil function (ADR 0032): grid N, apodization, focus_shift, aiming, warning thresholds.
  PupilFunctionOptions pupil;
  /// Wavelength indices; none: all system wavelengths with their model weights.
  std::optional<std::vector<std::uint16_t>> wavelengths;
  MtfOutputs outputs;
  /// Warning psf.pupil_undersampled above this phase step between neighbouring arrived cells,
  /// waves at the evaluated wavelength; > 0.
  double undersampled_waves = 0.5;
};

/// One MTF curve along a pupil direction, samples l = 0 ... N.
struct MtfCurve {
  /// Unit direction of nu(l) for l > 0 in (e_x, e_y) of the PSF plane (reference wavelength).
  std::array<double, 2> direction{};
  /// |nu(l)|, cycles per mm, increasing; the last sample (l = N) is the cut-off, where the
  /// shifted pupil no longer overlaps and the value is exactly 0.
  std::vector<double> frequency;
  /// MTF = |OTF|, 1 at l = 0 (C(0) is replaced by the exact energy E, ADR 0033, point 5).
  std::vector<double> modulation;
  /// If MtfOutputs::complex_otf: OTF(nu) = integral I exp(-2 pi i nu . x) / integral I, which is
  /// conj(C(l)) / E (erratum to point 5); a PSF moved by t multiplies it by exp(-2 pi i nu . t).
  /// Empty otherwise.
  std::vector<std::complex<double>> otf;
};

/// Diffraction MTF of one field.
struct Mtf {
  std::uint16_t field = 0;
  int grid = 0;                            ///< N
  math::Vec3 centre = math::Vec3::Zero();  ///< C', global, mm
  math::Vec3 e_x = math::Vec3::UnitX();    ///< first axis of the PSF plane, global
  math::Vec3 e_y = math::Vec3::UnitY();    ///< second axis of the PSF plane, global
  /// Shifts l = (0, l) in the pupil (py; meridional for a field along +-y, the convention of the
  /// project).
  MtfCurve tangential;
  /// Shifts l = (l, 0) in the pupil (px).
  MtfCurve sagittal;
  /// F = M a / lambda of the reference wavelength, row-major, cycles per mm per cell:
  /// nu(l) = (F[0] l1 + F[1] l2, F[2] l1 + F[3] l2) in (e_x, e_y).
  std::array<double, 4> frequency_matrix{};
  /// If MtfOutputs::map: (2N - 1)^2 values of |OTF|, row l2 + N - 1, column l1 + N - 1.
  std::vector<double> map;
  /// If MtfOutputs::map and complex_otf: the OTF in the same order.
  std::vector<std::complex<double>> otf_map;
  double edge_energy_deficit = 0.0;  ///< 1 - sum f^2 A^2 / sum f A^2, weighted
  /// Per wavelength: index, weight, energy E, length 2 N. Polychromatic: the OTF of each
  /// wavelength is interpolated linearly (real and imaginary part) in |nu| onto the frequencies of
  /// the reference wavelength (0 above its own cut-off) and summed with w_i E_i (ADR 0033, 7).
  std::vector<PsfWavelength> wavelengths;
  RayLosses losses;
  /// rays.lost, stop.clips_beam, pupil.mapping_nonlinear (from the pupil functions),
  /// psf.pupil_undersampled.
  std::vector<model::Diagnostic> warnings;
};

/// Diffraction MTF of `field` on `path` (ADR 0033, point 5).
/// @throws std::invalid_argument for an unknown wavelength, an empty wavelength list,
///         undersampled_waves <= 0, outputs.map with more than one wavelength (and as
///         pupil_function())
/// @throws AnalysisError (location /paths/i) for an afocal image space (no mapping), a path with
///         crystal modes or no arrived cell (and as pupil_function())
[[nodiscard]] Mtf mtf(const compile::CompiledSystem& system,
                      compile::PathId path,
                      std::uint16_t field,
                      const MtfOptions& options = {});

/// mtf() with cancellation and progress (run_control.hpp).
[[nodiscard]] Mtf mtf(const compile::CompiledSystem& system,
                      compile::PathId path,
                      std::uint16_t field,
                      const MtfOptions& options,
                      const trace::RunControl& control);

/// MTF of one pupil function (monochromatic), e.g. after refocus(). Losses and warnings are
/// copied from the pupil function; psf.pupil_undersampled is not checked.
/// @throws AnalysisError if the pupil function has no mapping (afocal) or no arrived cell
[[nodiscard]] Mtf mtf(const PupilFunction& pupil, const MtfOutputs& outputs = {});

/// Centre of the encircled-energy circle.
enum class EeCentre : std::uint8_t {
  /// C', the centre of the reference sphere (H11).
  Reference,
  /// The point x0 in the PSF plane given by the phase of the first harmonic of the periodic PSF:
  /// (a / lambda) M^T x0 = (arg C(e1), arg C(e2)) / (2 pi). Exact for a tilted wavefront (a moved
  /// PSF). It is **not** the first moment (centroid) of the PSF: for an asymmetric PSF such as
  /// coma it differs from it by about (2 pi / P)^2 <x^3> / 6, P = N lambda F#_w the period, i.e.
  /// about 1e-2 lambda F#_w at N = 64.
  LinearPhase,
};

/// Options of encircled_energy().
struct EncircledEnergyOptions {
  /// Pupil function (ADR 0032): grid N, apodization, focus_shift, aiming, warning thresholds.
  PupilFunctionOptions pupil;
  /// Wavelength indices; none: all system wavelengths with their model weights.
  std::optional<std::vector<std::uint16_t>> wavelengths;
  EeCentre centre = EeCentre::Reference;
  /// Warning psf.pupil_undersampled, as in MtfOptions; > 0.
  double undersampled_waves = 0.5;
};

/// Encircled energy of one field.
struct EncircledEnergy {
  std::uint16_t field = 0;
  int grid = 0;  ///< N
  /// Centre of the circles, global, mm: C' + offset.x e_x + offset.y e_y.
  math::Vec3 centre = math::Vec3::Zero();
  math::Vec3 e_x = math::Vec3::UnitX();  ///< first axis of the PSF plane, global
  math::Vec3 e_y = math::Vec3::UnitY();  ///< second axis of the PSF plane, global
  /// Centre relative to C' in the PSF plane, mm; 0 for EeCentre::Reference. Polychromatic:
  /// sum w_i E_i x0_i / sum w_i E_i.
  Point2 offset;
  std::vector<double> radii;  ///< mm, as given
  /// EE(r) = (|det M| a^2 / lambda^2) sum_l Re[C(l) exp(-2 pi i nu(l) . x0)] K_r(nu(l)) / E with
  /// K_r(nu) = r J1(2 pi r |nu|) / |nu|, K_r(0) = pi r^2 and C(0) = sum |u|^2 (ADR 0033, point 6
  /// with the erratum): fraction of the pupil energy E; polychromatic sum w_i E_i EE_i / sum w_i
  /// E_i at the same radius.
  std::vector<double> energy;
  /// Largest allowed radius, mm: a quarter of the shortest period lambda / a |M^-T e_i| of the
  /// periodic PSF (N lambda F#_w / 4 on the axis), the minimum over the wavelengths.
  double radius_max = 0.0;
  double edge_energy_deficit = 0.0;        ///< 1 - sum f^2 A^2 / sum f A^2, weighted
  std::vector<PsfWavelength> wavelengths;  ///< index, weight, energy E, length 2 N
  RayLosses losses;
  /// rays.lost, stop.clips_beam, pupil.mapping_nonlinear (from the pupil functions),
  /// psf.pupil_undersampled.
  std::vector<model::Diagnostic> warnings;
};

/// Encircled energy of `field` on `path` at `radii` (mm in the PSF plane; ADR 0033, point 6).
/// @throws std::invalid_argument for a radius < 0, not finite or > radius_max, an unknown
///         wavelength, an empty wavelength list, undersampled_waves <= 0 (and as
///         pupil_function()); an empty `radii` is allowed and gives empty vectors
/// @throws AnalysisError as mtf()
[[nodiscard]] EncircledEnergy encircled_energy(const compile::CompiledSystem& system,
                                               compile::PathId path,
                                               std::uint16_t field,
                                               std::span<const double> radii,
                                               const EncircledEnergyOptions& options = {});

/// encircled_energy() with cancellation and progress (run_control.hpp).
[[nodiscard]] EncircledEnergy encircled_energy(const compile::CompiledSystem& system,
                                               compile::PathId path,
                                               std::uint16_t field,
                                               std::span<const double> radii,
                                               const EncircledEnergyOptions& options,
                                               const trace::RunControl& control);

/// Encircled energy of one pupil function (monochromatic), e.g. after refocus(). Losses and
/// warnings are copied from the pupil function; psf.pupil_undersampled is not checked.
/// @throws std::invalid_argument for a radius < 0, not finite or > radius_max
/// @throws AnalysisError if the pupil function has no mapping (afocal) or no arrived cell
[[nodiscard]] EncircledEnergy encircled_energy(const PupilFunction& pupil,
                                               std::span<const double> radii,
                                               EeCentre centre = EeCentre::Reference);

}  // namespace rtt::analysis
