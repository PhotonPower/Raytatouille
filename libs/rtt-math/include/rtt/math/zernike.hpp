#pragma once

/// @file zernike.hpp
/// Zernike circle polynomials: Noll indices with orthonormal polynomials, the Fringe set without
/// normalisation and the radial polynomials with their derivative (ADR 0032, point 6; ADR 0051,
/// points 2 and 3). Shared by the surface terms (rtt-geom) and the wavefront fit (rtt-analysis),
/// so that surface and wavefront coefficients agree.
///
/// Conventions:
/// - K. Niu, C. Tian, "Zernike polynomials and their applications", J. Opt. 24 (2022) 123001:
///   - Noll: Z_j = sqrt(2 (n + 1)) R_n^m(rho) cos(m theta) for even j and m != 0,
///     sqrt(2 (n + 1)) R_n^m(rho) sin(m theta) for odd j and m != 0, sqrt(n + 1) R_n^0(rho) for
///     m = 0 (Eq. (3)); j <-> (n, m) after Eqs. (8) and (9); j starts at 1 (piston).
///   - R_n^m after Eq. (4) (the same polynomial as Eq. (23) of the Fringe set).
///   - Fringe: Z = R_n^|m|(rho) cos(m theta) for m >= 0, R_n^|m|(rho) sin(|m| theta) for m < 0,
///     without normalisation factors (Eq. (22)); the order of Table 4 (section 2.1.3). The file
///     counts from 1 like Zemax and CODE V (section 2.1.3 and footnote a of Table 8), so Fringe
///     term k is j = k - 1 of Table 4.
/// - rho = r / R <= 1 for the normalisation radius R (any rho >= 0 is evaluated); theta from +x
///   counter-clockwise (Niu/Tian, Sec. 2.1).
/// - Functions throw std::invalid_argument for indices out of range (input errors at the API
///   boundary); the `_unchecked` evaluation for the hot path does not.

#include <cstddef>
#include <span>
#include <utility>

namespace rtt::math {

/// Largest radial degree n of the tables.
inline constexpr int kMaxZernikeDegree = 20;
/// Largest Noll index: n <= 20 completely (231 terms). Up to n = 20 the integer coefficients of
/// R_n^m are below 2^53 (the largest 2 333 760), so exact in double (ADR 0032, point 6).
inline constexpr int kMaxZernikeTerms = 231;
/// Number of Fringe terms (Niu/Tian, Table 4: j = 0 ... 36).
inline constexpr int kFringeTerms = 37;

/// Radial degree n and azimuthal frequency m >= 0 of a Noll index; n - m is even. Whether the
/// angular factor is the cosine or the sine follows from the parity of j (Niu/Tian Eq. (3)).
struct NollIndex {
  int n = 0;
  int m = 0;
};

/// Radial degree n, azimuthal frequency m >= 0 and the angular factor of a Fringe term.
struct FringeIndex {
  int n = 0;
  int m = 0;          ///< |m| of Niu/Tian Table 4
  bool sine = false;  ///< sin(m theta) for a negative m of Table 4, else cos(m theta)
};

/// (n, m) of the Noll index j (Niu/Tian Eq. (8), evaluated in integers).
/// @throws std::invalid_argument for j outside [1, kMaxZernikeTerms]
[[nodiscard]] NollIndex noll_to_nm(int j);

/// The Noll index of (n, m) with the sine (`sine`, odd j) or the cosine (even j) for m != 0
/// (Niu/Tian Eq. (9): of the pair n(n+1)/2 + m, n(n+1)/2 + m + 1 the one with the parity);
/// `sine` is ignored for m = 0.
/// @throws std::invalid_argument unless 0 <= m <= n <= 20 and n - m is even
[[nodiscard]] int noll_from_nm(int n, int m, bool sine);

/// (n, |m|, sine) of the Fringe term k = 1 ... 37 (file numbering, Niu/Tian Table 4 at j = k - 1):
/// j = 0 ... 35 are all (n, m) with N = (n + |m|)/2 <= 5, ordered by N, then n ascending, then m
/// descending (+|m| before -|m|); j = 36 is (12, 0).
/// @throws std::invalid_argument for k outside [1, kFringeTerms]
[[nodiscard]] FringeIndex fringe_to_nm(int k);

/// The integer coefficients c_0 ... c_K of R_n^m = rho^m sum_s c_s x^(K - s), x = rho^2,
/// K = (n - m)/2 (Niu/Tian Eq. (4): c_s = (-1)^s (n - s)! / (s! ((n + m)/2 - s)!
/// ((n - m)/2 - s)!)), exact as doubles; c_0 belongs to the highest power.
/// @throws std::invalid_argument unless 0 <= m <= n <= 20 and n - m is even
[[nodiscard]] std::span<const double> zernike_radial_coefficients(int n, int m);

/// R_n^m(rho) (Niu/Tian Eq. (4)), for 0 <= m <= n <= 20, n - m even. Evaluated as
/// rho^m sum_s c_s x^(K - s) with x = rho^2, Horner in x, the integer c_s exact; rounding error
/// at most gamma_k sum_s |c_s| rho^(n - 2s), k = 3K + m + 1 (ADR 0032, point 6).
/// @throws std::invalid_argument for (n, m) outside that range
[[nodiscard]] double zernike_radial(int n, int m, double rho);

/// R_n^m(rho) without argument checks (@pre as zernike_radial()).
[[nodiscard]] double zernike_radial_unchecked(int n, int m, double rho) noexcept;

/// (R_n^m(rho) / rho, dR_n^m/drho) without argument checks, for the gradient (ADR 0051, point 2):
/// both as polynomials, R / rho = rho^(m - 1) P(x) for m >= 1 (no cancellation near rho = 0), and
/// for m = 0 the first value is 0 (it is only used with m >= 1).
/// @pre 0 <= m <= n <= 20, n - m even
[[nodiscard]] std::pair<double, double> zernike_radial_over_rho_and_derivative(int n,
                                                                               int m,
                                                                               double rho) noexcept;

/// Z_j(rho, theta), orthonormal Noll polynomial (Niu/Tian Eq. (3)); theta in rad from +x
/// counter-clockwise.
/// @throws std::invalid_argument for j outside [1, kMaxZernikeTerms]
[[nodiscard]] double zernike(int j, double rho, double theta);

/// The normalisation factor of the Noll polynomial with (n, m): sqrt(n + 1) for m = 0, else
/// sqrt(2 (n + 1)) (Niu/Tian Eq. (3)). Z_Noll(n, m) = factor * Z_Fringe(n, m).
[[nodiscard]] double zernike_noll_factor(int n, int m) noexcept;

}  // namespace rtt::math
