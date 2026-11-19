// Zernike circle polynomials (ADR 0032, point 6; ADR 0051, points 2 and 3).

#include "rtt/math/zernike.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace rtt::math {

namespace {

void check_noll(int j) {
  if (j < 1 || j > kMaxZernikeTerms) {
    throw std::invalid_argument("Zernike: Noll index must lie in [1, " +
                                std::to_string(kMaxZernikeTerms) + "], got " + std::to_string(j));
  }
}

void check_nm(int n, int m) {
  if (n < 0 || n > kMaxZernikeDegree || m < 0 || m > n || (n - m) % 2 != 0) {
    throw std::invalid_argument("Zernike: (n, m) = (" + std::to_string(n) + ", " +
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
  std::array<std::array<std::vector<double>, kMaxZernikeDegree + 1>, kMaxZernikeDegree + 1> c;
  RadialTable() {
    for (int n = 0; n <= kMaxZernikeDegree; ++n) {
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

/// The Fringe set in the order of Niu/Tian Table 4 (section 2.1.3): N = (n + |m|)/2 ascending,
/// then n ascending, then m descending, for N <= 5 (36 terms), then (12, 0).
std::array<FringeIndex, kFringeTerms> make_fringe_table() {
  std::array<FringeIndex, kFringeTerms> t{};
  std::size_t k = 0;
  for (int big_n = 0; big_n <= 5; ++big_n) {
    for (int n = big_n; n <= 2 * big_n; ++n) {
      const int m = 2 * big_n - n;  // |m|
      t[k++] = {n, m, false};
      if (m != 0) t[k++] = {n, m, true};
    }
  }
  t[k] = {12, 0, false};
  return t;
}

const std::array<FringeIndex, kFringeTerms>& fringe_table() {
  static const std::array<FringeIndex, kFringeTerms> table = make_fringe_table();
  return table;
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

FringeIndex fringe_to_nm(int k) {
  if (k < 1 || k > kFringeTerms) {
    throw std::invalid_argument("Zernike: Fringe term must lie in [1, " +
                                std::to_string(kFringeTerms) + "], got " + std::to_string(k));
  }
  return fringe_table()[static_cast<std::size_t>(k - 1)];
}

std::span<const double> zernike_radial_coefficients(int n, int m) {
  check_nm(n, m);
  return radial_table().c[n][m];
}

double zernike_radial_unchecked(int n, int m, double rho) noexcept {
  const std::vector<double>& c = radial_table().c[n][m];
  const double x = rho * rho;
  double acc = c[0];
  for (std::size_t s = 1; s < c.size(); ++s) acc = acc * x + c[s];
  for (int i = 0; i < m; ++i) acc *= rho;
  return acc;
}

double zernike_radial(int n, int m, double rho) {
  check_nm(n, m);
  return zernike_radial_unchecked(n, m, rho);
}

std::pair<double, double> zernike_radial_over_rho_and_derivative(int n,
                                                                 int m,
                                                                 double rho) noexcept {
  // R = rho^m Q(x), x = rho^2, Q = sum_s c_s x^(K - s). Horner for Q and Q' together.
  // dR/drho = m rho^(m - 1) Q + rho^m 2 rho Q' = rho^(m - 1) (m Q + 2 x Q') for m >= 1, and
  // 2 rho Q' for m = 0 (derived from Niu/Tian Eq. (4)).
  const std::vector<double>& c = radial_table().c[n][m];
  const double x = rho * rho;
  double q = c[0];
  double dq = 0.0;
  for (std::size_t s = 1; s < c.size(); ++s) {
    dq = dq * x + q;
    q = q * x + c[s];
  }
  if (m == 0) return {0.0, 2.0 * rho * dq};
  double pm1 = 1.0;  // rho^(m - 1)
  for (int i = 1; i < m; ++i) pm1 *= rho;
  return {pm1 * q, pm1 * (m * q + 2.0 * x * dq)};
}

double zernike(int j, double rho, double theta) {
  check_noll(j);
  const NollIndex nm = noll_to_nm(j);
  const double r = zernike_radial_unchecked(nm.n, nm.m, rho);
  if (nm.m == 0) return std::sqrt(nm.n + 1.0) * r;
  const double norm = std::sqrt(2.0 * (nm.n + 1.0));
  return j % 2 == 0 ? norm * r * std::cos(nm.m * theta) : norm * r * std::sin(nm.m * theta);
}

double zernike_noll_factor(int n, int m) noexcept {
  return m == 0 ? std::sqrt(n + 1.0) : std::sqrt(2.0 * (n + 1.0));
}

}  // namespace rtt::math
