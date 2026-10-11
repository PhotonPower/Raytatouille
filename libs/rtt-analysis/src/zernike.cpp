#include "rtt/analysis/zernike.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

#include "rtt/diagnostics/codes.hpp"

namespace rtt::analysis {

namespace {

constexpr int kMaxDegree = 20;  ///< n of kMaxZernikeTerms
/// Pivot of the Cholesky factorisation at most this times the largest diagonal element: the
/// normal equations are too badly conditioned (ADR 0032, point 6).
constexpr double kPivotFloor = 1e-8;

void check_noll(int j) {
  if (j < 1 || j > kMaxZernikeTerms) {
    throw std::invalid_argument("analysis: Noll index must lie in [1, " +
                                std::to_string(kMaxZernikeTerms) + "], got " + std::to_string(j));
  }
}

void check_nm(int n, int m) {
  if (n < 0 || n > kMaxDegree || m < 0 || m > n || (n - m) % 2 != 0) {
    throw std::invalid_argument("analysis: Zernike (n, m) = (" + std::to_string(n) + ", " +
                                std::to_string(m) + ") needs 0 <= m <= n <= 20, n - m even");
  }
}

/// c_s of R_n^m (Niu/Tian Eq. (4)): (-1)^s (n - s)! / (s! ((n + m)/2 - s)! ((n - m)/2 - s)!).
/// Exact in uint64 (20! < 2^64) and, for n <= 20, below 2^53, so exact as a double.
double radial_coefficient(int n, int m, int s) {
  std::uint64_t num = 1;
  for (int i = 2; i <= n - s; ++i) num *= static_cast<std::uint64_t>(i);
  std::uint64_t den = 1;
  for (int i = 2; i <= s; ++i) den *= static_cast<std::uint64_t>(i);
  for (int i = 2; i <= (n + m) / 2 - s; ++i) den *= static_cast<std::uint64_t>(i);
  for (int i = 2; i <= (n - m) / 2 - s; ++i) den *= static_cast<std::uint64_t>(i);
  const std::uint64_t quotient = num / den;  // exact: the coefficient is an integer
  const auto c = static_cast<double>(quotient);
  return s % 2 == 0 ? c : -c;
}

/// Coefficients c_0 ... c_K of every (n, m) up to n = 20, computed once.
struct RadialTable {
  std::array<std::array<std::vector<double>, kMaxDegree + 1>, kMaxDegree + 1> c;
  RadialTable() {
    for (int n = 0; n <= kMaxDegree; ++n) {
      for (int m = n % 2; m <= n; m += 2) {
        for (int s = 0; s <= (n - m) / 2; ++s) c[n][m].push_back(radial_coefficient(n, m, s));
      }
    }
  }
};

const RadialTable& radial_table() {
  static const RadialTable table;  // immutable after construction (thread-safe initialisation)
  return table;
}

/// R_n^m without argument checks.
double radial(int n, int m, double rho) {
  const std::vector<double>& c = radial_table().c[n][m];
  const double x = rho * rho;
  double acc = c[0];
  for (std::size_t s = 1; s < c.size(); ++s) acc = acc * x + c[s];
  for (int i = 0; i < m; ++i) acc *= rho;
  return acc;
}

/// Z_j without argument checks; (n, m) given.
double zernike_nm(int j, NollIndex nm, double rho, double theta) {
  const double r = radial(nm.n, nm.m, rho);
  if (nm.m == 0) return std::sqrt(nm.n + 1.0) * r;
  const double norm = std::sqrt(2.0 * (nm.n + 1.0));
  return j % 2 == 0 ? norm * r * std::cos(nm.m * theta) : norm * r * std::sin(nm.m * theta);
}

std::string path_pointer(compile::PathId path) {
  return "/paths/" + std::to_string(path.index);
}

/// Weighted least squares sum_k f_k (W_k - sum_j c_j A_kj)^2 over the rows `a` (cells x terms,
/// row-major): the normal equations A^T F A c = A^T F W (Niu/Tian Eq. (54) with the weights f;
/// derived in ADR 0032, point 6 as the Riemann sum of Eq. (7)), solved by the Cholesky
/// factorisation A^T F A = C^T C after K. Madsen, H. B. Nielsen, O. Tingleff, "Methods for
/// Non-Linear Least Squares Problems" (2004), Appendix A, Algorithm A.4 (as in rtt-optim), with
/// the pivot floor of ADR 0032 instead of d > 0. Sums in index order. False if a pivot is too
/// small.
bool weighted_least_squares(const std::vector<double>& a,
                            std::size_t terms,
                            const std::vector<double>& f,
                            const std::vector<double>& w,
                            std::vector<double>& c) {
  const std::size_t rows = f.size();
  std::vector<double> ata(terms * terms, 0.0);
  std::vector<double> atw(terms, 0.0);
  for (std::size_t k = 0; k < rows; ++k) {
    const double* row = a.data() + k * terms;
    for (std::size_t i = 0; i < terms; ++i) {
      const double fi = f[k] * row[i];
      atw[i] += fi * w[k];
      for (std::size_t j = i; j < terms; ++j) ata[i * terms + j] += fi * row[j];
    }
  }
  double max_diag = 0.0;
  for (std::size_t i = 0; i < terms; ++i) max_diag = std::max(max_diag, ata[i * terms + i]);
  // Factorisation into the upper triangle (MNT Alg. A.4).
  for (std::size_t k = 0; k < terms; ++k) {
    double d = ata[k * terms + k];
    for (std::size_t i = 0; i < k; ++i) d -= ata[i * terms + k] * ata[i * terms + k];
    if (!(d > kPivotFloor * max_diag)) return false;
    const double ckk = std::sqrt(d);
    ata[k * terms + k] = ckk;
    for (std::size_t j = k + 1; j < terms; ++j) {
      double s = ata[k * terms + j];
      for (std::size_t i = 0; i < k; ++i) s -= ata[i * terms + j] * ata[i * terms + k];
      ata[k * terms + j] = s / ckk;
    }
  }
  c.assign(terms, 0.0);
  for (std::size_t k = 0; k < terms; ++k) {  // C^T z = A^T F W, z in c
    double s = atw[k];
    for (std::size_t i = 0; i < k; ++i) s -= ata[i * terms + k] * c[i];
    c[k] = s / ata[k * terms + k];
  }
  for (std::size_t k = terms; k-- > 0;) {  // C c = z
    double s = c[k];
    for (std::size_t j = k + 1; j < terms; ++j) s -= ata[k * terms + j] * c[j];
    c[k] = s / ata[k * terms + k];
  }
  return true;
}

/// f-weighted RMS of W - A c over the rows.
double residual_rms(const std::vector<double>& a,
                    std::size_t stride,
                    std::size_t terms,
                    const std::vector<double>& f,
                    const std::vector<double>& w,
                    const std::vector<double>& c) {
  double sum_f = 0.0;
  double sum = 0.0;
  for (std::size_t k = 0; k < f.size(); ++k) {
    double r = w[k];
    for (std::size_t j = 0; j < terms; ++j) r -= c[j] * a[k * stride + j];
    sum_f += f[k];
    sum += f[k] * r * r;
  }
  return std::sqrt(sum / sum_f);
}

}  // namespace

NollIndex noll_to_nm(int j) {
  check_noll(j);
  // Niu/Tian Eq. (8): n = floor(sqrt(2j - 1) + 1/2) - 1, i.e. the n with n(n+1)/2 < j <=
  // (n+1)(n+2)/2, here in integers; m from the floors of Eq. (8).
  int n = 0;
  while ((n + 1) * (n + 2) / 2 < j) ++n;
  const int m =
      n % 2 == 0 ? 2 * ((2 * j + 1 - n * (n + 1)) / 4) : 2 * ((2 * (j + 1) - n * (n + 1)) / 4) - 1;
  return {n, m};
}

int noll_from_nm(int n, int m, bool sine) {
  check_nm(n, m);
  const int base = n * (n + 1) / 2;
  if (m == 0) return base + 1;
  // Niu/Tian Eq. (9): j in {base + m, base + m + 1}; the odd one carries the sine (Eq. (3)).
  const int first = base + m;
  const bool first_odd = first % 2 == 1;
  return first_odd == sine ? first : first + 1;
}

double zernike_radial(int n, int m, double rho) {
  check_nm(n, m);
  return radial(n, m, rho);
}

double zernike(int j, double rho, double theta) {
  check_noll(j);
  return zernike_nm(j, noll_to_nm(j), rho, theta);
}

ZernikeFit zernike_fit(const PupilFunction& pupil, const ZernikeOptions& options) {
  if (options.terms < 1 || options.terms > kMaxZernikeTerms) {
    throw std::invalid_argument("analysis: Zernike terms must lie in [1, " +
                                std::to_string(kMaxZernikeTerms) + "]");
  }
  if (!(options.incomplete_warning_fraction >= 0.0 && options.incomplete_warning_fraction <= 1.0)) {
    throw std::invalid_argument("analysis: incomplete_warning_fraction must lie in [0, 1]");
  }
  ZernikeFit fit;
  fit.path = pupil.path;
  fit.field = pupil.field;
  fit.wavelength = pupil.wavelength;
  fit.wavelength_um = pupil.wavelength_um;
  fit.reference_wavelength_um = pupil.reference_wavelength_um;
  fit.arrived_fraction = pupil.arrived_fraction;

  // Rows: the arrived cells with f > 0 in cell order; the design matrix with max(J, 4) columns
  // (the first four also serve the fit about the diffraction focus).
  const auto terms = static_cast<std::size_t>(options.terms);
  const std::size_t stride = std::max<std::size_t>(terms, 4);
  std::vector<NollIndex> nm(stride);
  for (std::size_t j = 0; j < stride; ++j) nm[j] = noll_to_nm(static_cast<int>(j + 1));
  std::vector<double> f;
  std::vector<double> w;
  std::vector<double> a;
  for (std::size_t k = 0; k < pupil.fraction.size(); ++k) {
    if (pupil.fraction[k] == 0.0 || pupil.status[k] != trace::RayStatus::Alive) continue;
    const double rho = std::hypot(pupil.px[k], pupil.py[k]);
    const double theta = std::atan2(pupil.py[k], pupil.px[k]);
    f.push_back(pupil.fraction[k]);
    w.push_back(pupil.w[k]);
    for (std::size_t j = 0; j < stride; ++j) {
      a.push_back(zernike_nm(static_cast<int>(j + 1), nm[j], rho, theta));
    }
  }
  fit.cells = f.size();
  if (fit.cells < stride) {
    throw AnalysisError("analysis: Zernike fit with " + std::to_string(options.terms) +
                            " terms needs at least " + std::to_string(stride) +
                            " arrived cells, got " + std::to_string(fit.cells),
                        path_pointer(pupil.path));
  }
  // The full fit on the columns 0 ... J - 1 of the rows (packed copy with stride J).
  std::vector<double> a_terms;
  a_terms.reserve(fit.cells * terms);
  for (std::size_t k = 0; k < fit.cells; ++k) {
    for (std::size_t j = 0; j < terms; ++j) a_terms.push_back(a[k * stride + j]);
  }
  if (!weighted_least_squares(a_terms, terms, f, w, fit.coefficients)) {
    throw AnalysisError(
        "analysis: the Zernike normal equations are singular (pivot below 1e-8 of "
        "the largest diagonal element); too few or badly placed cells for " +
            std::to_string(options.terms) + " terms",
        path_pointer(pupil.path));
  }
  fit.residual_rms = residual_rms(a, stride, terms, f, w, fit.coefficients);

  // sigma: f-weighted standard deviation of W (piston removed, Niu/Tian Table 10 on the disk).
  double sum_f = 0.0;
  double sum_fw = 0.0;
  for (std::size_t k = 0; k < fit.cells; ++k) {
    sum_f += f[k];
    sum_fw += f[k] * w[k];
  }
  const double mean = sum_fw / sum_f;
  double var = 0.0;
  for (std::size_t k = 0; k < fit.cells; ++k) var += f[k] * (w[k] - mean) * (w[k] - mean);
  fit.sigma = std::sqrt(var / sum_f);

  // sigma_BF: own fit with Z1 ... Z4 (ADR 0032, point 7): on the grid the Z_j are not exactly
  // orthogonal, so the coefficients of the large fit do not give the minimal variance.
  std::vector<double> a4;
  a4.reserve(fit.cells * 4);
  for (std::size_t k = 0; k < fit.cells; ++k) {
    for (std::size_t j = 0; j < 4; ++j) a4.push_back(a[k * stride + j]);
  }
  std::vector<double> c4;
  if (!weighted_least_squares(a4, 4, f, w, c4)) {
    throw AnalysisError(
        "analysis: the Zernike normal equations of piston, tilt and defocus are "
        "singular",
        path_pointer(pupil.path));
  }
  fit.sigma_bf = residual_rms(a, stride, 4, f, w, c4);
  // Wyant/Creath Eq. (67): S = exp(-(2 pi sigma)^2), sigma in waves at the evaluated wavelength.
  const double sigma_lambda = fit.sigma_bf * pupil.reference_wavelength_um / pupil.wavelength_um;
  fit.strehl =
      std::exp(-(2.0 * std::numbers::pi * sigma_lambda) * (2.0 * std::numbers::pi * sigma_lambda));

  if (1.0 - pupil.arrived_fraction > options.incomplete_warning_fraction) {
    constexpr diagnostics::DiagnosticCode kCode = "zernike.pupil_incomplete";
    fit.warnings.push_back(
        {kCode.severity(), path_pointer(pupil.path),
         "Zernike fit on an incomplete pupil: " +
             std::to_string(100.0 * (1.0 - pupil.arrived_fraction)) +
             " % of its area did not arrive; the Z_j are not orthogonal there, the coefficients "
             "are a fit, not a decomposition",
         std::string(kCode.str())});
  } else if (pupil.obscuration > 0.0) {
    // An annular pupil (ADR 0032, addendum #229): the Noll polynomials are not orthogonal on the
    // annulus either, although every cell of it arrived (annular polynomials would be a follow-up).
    constexpr diagnostics::DiagnosticCode kCode = "zernike.pupil_incomplete";
    fit.warnings.push_back(
        {kCode.severity(), path_pointer(pupil.path),
         "Zernike fit on an annular pupil (central obscuration " +
             std::to_string(pupil.obscuration) +
             " of the pupil radius): the Z_j are not orthogonal there, the coefficients are a "
             "fit, not a decomposition",
         std::string(kCode.str())});
  }
  return fit;
}

}  // namespace rtt::analysis
