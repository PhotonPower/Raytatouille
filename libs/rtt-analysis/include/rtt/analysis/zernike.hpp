#pragma once

/// @file zernike.hpp
/// Zernike polynomials after Noll, the wavefront fit on the pupil function and the Strehl
/// approximation after Marechal (ADR 0032, points 6 and 7; #203).
///
/// Conventions:
/// - Noll indices and orthonormal polynomials after K. Niu, C. Tian, "Zernike polynomials and
///   their applications", J. Opt. 24 (2022) 123001: Z_j = sqrt(2 (n + 1)) R_n^m(rho) cos(m theta)
///   for even j and m != 0, sqrt(2 (n + 1)) R_n^m(rho) sin(m theta) for odd j and m != 0,
///   sqrt(n + 1) R_n^0(rho) for m = 0 (Eq. (3)); R_n^m after Eq. (4); j <-> (n, m) after Eqs. (8)
///   and (9); the mean of Z_j Z_k over the unit disk is delta_jk (Eq. (7)). j starts at 1
///   (piston), the order of model::ZernikeSag.
/// - rho = |p| <= 1 with the normalised pupil coordinates p = (px, py) (rim of the paraxial
///   entrance pupil); theta = atan2(py, px), from +x counter-clockwise (Niu/Tian, Sec. 2.1). Note:
///   trace::GaussPupil counts its arm angle from +y towards +x.
/// - Coefficients and RMS values in waves at the reference wavelength, like W.

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "rtt/analysis/pupil.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/model/validate.hpp"

namespace rtt::analysis {

/// Largest Noll index: n <= 20 completely (231 terms). Up to n = 20 the integer coefficients of
/// R_n^m are below 2^53 (the largest 2 333 760), so exact in double (ADR 0032, point 6).
inline constexpr int kMaxZernikeTerms = 231;

/// Radial degree n and azimuthal frequency m >= 0 of a Noll index; n - m is even. Whether the
/// angular factor is the cosine or the sine follows from the parity of j (Niu/Tian Eq. (3)).
struct NollIndex {
  int n = 0;
  int m = 0;
};

/// (n, m) of the Noll index j (Niu/Tian Eq. (8), evaluated in integers).
/// @throws std::invalid_argument for j outside [1, kMaxZernikeTerms]
[[nodiscard]] NollIndex noll_to_nm(int j);

/// The Noll index of (n, m) with the sine (`sine`, odd j) or the cosine (even j) for m != 0
/// (Niu/Tian Eq. (9): of the pair n(n+1)/2 + m, n(n+1)/2 + m + 1 the one with the parity);
/// `sine` is ignored for m = 0.
/// @throws std::invalid_argument unless 0 <= m <= n <= 20 and n - m is even
[[nodiscard]] int noll_from_nm(int n, int m, bool sine);

/// R_n^m(rho) after Niu/Tian Eq. (4), for 0 <= m <= n <= 20, n - m even. Evaluated as
/// rho^m sum_s c_s x^(K - s) with x = rho^2, K = (n - m) / 2, Horner in x, the integer c_s exact;
/// rounding error at most gamma_k sum_s |c_s| rho^(n - 2s), k = 3K + m + 1 (ADR 0032, point 6).
/// @throws std::invalid_argument for (n, m) outside that range
[[nodiscard]] double zernike_radial(int n, int m, double rho);

/// Z_j(rho, theta), orthonormal (Niu/Tian Eq. (3)); theta in rad from +x counter-clockwise.
/// @throws std::invalid_argument for j outside [1, kMaxZernikeTerms]
[[nodiscard]] double zernike(int j, double rho, double theta);

/// Options of the Zernike fit (ADR 0032, point 6).
struct ZernikeOptions {
  int terms = 37;  ///< J, the terms j = 1 ... J; in [1, kMaxZernikeTerms]
  /// Warning "zernike.pupil_incomplete" if 1 - arrived_fraction exceeds this; in [0, 1].
  double incomplete_warning_fraction = 0.01;
};

/// Zernike fit of a pupil function and the Strehl approximation (ADR 0032, points 6 and 7).
struct ZernikeFit {
  compile::PathId path{};
  std::uint16_t field = 0;
  std::uint16_t wavelength = 0;
  double wavelength_um = 0.0;            ///< evaluated wavelength, um
  double reference_wavelength_um = 0.0;  ///< unit of W and the coefficients, um
  /// c_j at index j - 1, waves at lambda_ref: weighted least squares over the arrived cells with
  /// f > 0, weight f (A^T F A c = A^T F W, Niu/Tian Eq. (54) with weights). On an incomplete
  /// pupil the Z_j are not orthogonal: the coefficients are a fit, not a decomposition.
  std::vector<double> coefficients;
  double residual_rms = 0.0;  ///< f-weighted RMS of W - sum c_j Z_j, waves
  double sigma = 0.0;         ///< f-weighted standard deviation of W (piston removed), waves
  /// f-weighted RMS of W after a fit of its own with Z_1 ... Z_4 (piston, tilt, defocus): the
  /// deviation about the diffraction focus (Wyant/Creath, Sec. IX), waves at lambda_ref.
  double sigma_bf = 0.0;
  /// Strehl ratio after Marechal, exp(-(2 pi sigma_BF lambda_ref / lambda)^2) at the evaluated
  /// wavelength lambda (Wyant/Creath Eq. (67), good down to about 0.1); uniform illumination
  /// assumed, the amplitude does not enter. Below 0.1 still given (documented, no warning).
  double strehl = 0.0;
  std::size_t cells = 0;          ///< arrived cells with f > 0 in the fit
  double arrived_fraction = 0.0;  ///< of the pupil function
  /// "zernike.pupil_incomplete" (location /paths/i) if 1 - arrived_fraction exceeds the option.
  std::vector<model::Diagnostic> warnings;
};

/// Zernike fit of `pupil` with options.terms terms (ADR 0032, point 6), with sigma, sigma_BF and
/// the Strehl approximation (point 7). Serial in cell order: the same result on every run.
/// @throws std::invalid_argument for options out of range
/// @throws AnalysisError with location /paths/i if fewer arrived cells with f > 0 than
///         max(terms, 4) exist, or a pivot of the Cholesky factorisation of A^T F A is at most
///         1e-8 times its largest diagonal element (the coefficients would lose about eight digits)
[[nodiscard]] ZernikeFit zernike_fit(const PupilFunction& pupil,
                                     const ZernikeOptions& options = {});

}  // namespace rtt::analysis
