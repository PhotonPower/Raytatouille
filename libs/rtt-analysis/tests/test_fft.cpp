// Internal DFT of rtt-analysis (#204; ADR 0033, point 1): pocketfft behind detail::fft1/fft2.
//
// Error bound (derived for #204, docs/quellen.md): a DFT of length L computed in at most
// ceil(log2 L') butterfly stages with coefficients of modulus <= 1 and a bounded number of
// roundings per output, times 3 for Bluestein (two transforms of length L' >= 2L - 1 and the chirp
// products), stays within c_fft(L) u sum_j |x_j| per element, with c_fft(L) = 30 ceil(log2(2L))
// and u = 2^-53. The reference (naive DFT in double, Neumaier summation, twiddles from
// std::cos/std::sin of 2 pi ((j k) mod L)/L) adds at most 12 u sum_j |x_j|: twiddle (2 pi + 1) u,
// product and summation 4 u. The bound is coarse on purpose: it separates the errors that matter
// here (sign convention, index, normalisation, wrong axis), which are O(1). The measured error is
// reported with INFO (the probe for ADR 0033 measured <= 2.3e-16 max|X|).

#include <algorithm>
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstring>
#include <numbers>
#include <stdexcept>
#include <thread>
#include <vector>

#include "../src/fft.hpp"

using rtt::analysis::detail::fft1;
using rtt::analysis::detail::fft2;
using rtt::analysis::detail::FftDirection;
using C = std::complex<double>;

namespace {

constexpr double kU = 0x1p-53;

double c_fft(std::size_t n) {
  return 30.0 * std::ceil(std::log2(2.0 * static_cast<double>(n)));
}

/// Element bound of a 1D transform of length n against the naive reference.
double bound1(std::size_t n, double sum_abs) {
  return (c_fft(n) + 12.0) * kU * sum_abs;
}

/// Element bound of a 2D transform rows x cols (two passes) against the naive reference.
double bound2(std::size_t rows, std::size_t cols, double sum_abs) {
  return (c_fft(rows) + c_fft(cols) + 24.0) * kU * sum_abs;
}

/// Deterministic input in [-0.5, 0.5)^2 (xorshift64).
std::vector<C> input(std::size_t n, std::uint64_t seed) {
  std::vector<C> v(n);
  std::uint64_t s = 0x9E3779B97F4A7C15ull * (seed + 1);
  const auto next = [&s] {
    s ^= s << 13;
    s ^= s >> 7;
    s ^= s << 17;
    return static_cast<double>(s >> 11) * 0x1p-53 - 0.5;
  };
  for (C& z : v) {
    const double re = next();
    z = C(re, next());
  }
  return v;
}

double sum_abs(const std::vector<C>& x) {
  double s = 0.0;
  for (const C& z : x) s += std::abs(z);
  return s;
}

/// Neumaier summation of complex values.
class Sum {
 public:
  void add(C z) {
    add_part(re_, cre_, z.real());
    add_part(im_, cim_, z.imag());
  }
  [[nodiscard]] C value() const { return {re_ + cre_, im_ + cim_}; }

 private:
  static void add_part(double& s, double& c, double v) {
    const double t = s + v;
    c += std::abs(s) >= std::abs(v) ? (s - t) + v : (v - t) + s;
    s = t;
  }
  double re_ = 0.0, cre_ = 0.0, im_ = 0.0, cim_ = 0.0;
};

C twiddle(std::size_t jk, std::size_t n, double sign) {
  const double a =
      sign * 2.0 * std::numbers::pi * static_cast<double>(jk % n) / static_cast<double>(n);
  return {std::cos(a), std::sin(a)};
}

/// Naive forward DFT, row-major rows x cols.
std::vector<C> naive(const std::vector<C>& x, std::size_t rows, std::size_t cols) {
  std::vector<C> y(rows * cols);
  for (std::size_t kr = 0; kr < rows; ++kr) {
    for (std::size_t kc = 0; kc < cols; ++kc) {
      Sum s;
      for (std::size_t jr = 0; jr < rows; ++jr) {
        for (std::size_t jc = 0; jc < cols; ++jc) {
          s.add(x[jr * cols + jc] * twiddle(jr * kr, rows, -1.0) * twiddle(jc * kc, cols, -1.0));
        }
      }
      y[kr * cols + kc] = s.value();
    }
  }
  return y;
}

/// max_k |a_k - b_k|.
double max_diff(const std::vector<C>& a, const std::vector<C>& b) {
  double m = 0.0;
  for (std::size_t k = 0; k < a.size(); ++k) m = std::max(m, std::abs(a[k] - b[k]));
  return m;
}

double max_abs(const std::vector<C>& a) {
  double m = 0.0;
  for (const C& z : a) m = std::max(m, std::abs(z));
  return m;
}

}  // namespace

TEST_CASE("fft: forward convention exp(-2 pi i j k / L), a single frequency lands in its bin",
          "[fft]") {
  // x_j = exp(+2 pi i m j / L) gives y_k = L delta_{k m} for the forward transform; a wrong sign
  // would put it into bin L - m.
  const std::size_t n = 64;
  const std::size_t m = 5;
  std::vector<C> x(n);
  for (std::size_t j = 0; j < n; ++j) x[j] = twiddle(m * j, n, +1.0);
  std::vector<C> y = x;
  fft1(y, FftDirection::Forward);
  const double tol = bound1(n, sum_abs(x));
  INFO("tolerance " << tol);
  for (std::size_t k = 0; k < n; ++k) {
    INFO("bin " << k << ": " << y[k]);
    CHECK(std::abs(y[k] - (k == m ? C(static_cast<double>(n), 0.0) : C(0.0, 0.0))) <= tol);
  }
  // Backward(Forward(x)) = L x.
  fft1(y, FftDirection::Backward);
  for (std::size_t j = 0; j < n; ++j) {
    CHECK(std::abs(y[j] / static_cast<double>(n) - x[j]) <= 3.0 * c_fft(n) * kU * sum_abs(x));
  }
}

TEST_CASE("fft: 1D against a naive DFT, also a length with a large prime factor", "[fft]") {
  // 1018 = 2 * 509: Bluestein or mixed radix, as the polychromatic lengths need (ADR 0033, 7).
  for (const std::size_t n :
       {std::size_t{8}, std::size_t{64}, std::size_t{100}, std::size_t{128}, std::size_t{1018}}) {
    const std::vector<C> x = input(n, n);
    std::vector<C> y = x;
    fft1(y, FftDirection::Forward);
    const std::vector<C> ref = naive(x, 1, n);
    const double err = max_diff(y, ref);
    INFO("L = " << n << ": max error " << err << " = " << err / max_abs(ref) << " max|X|");
    CHECK(err <= bound1(n, sum_abs(x)));
    // Round trip.
    fft1(y, FftDirection::Backward);
    for (C& z : y) z /= static_cast<double>(n);
    CHECK(max_diff(y, x) <= 3.0 * c_fft(n) * kU * sum_abs(x));
  }
}

TEST_CASE("fft: 2D against a naive DFT, square and non-square", "[fft]") {
  for (const auto [rows, cols] :
       {std::array<std::size_t, 2>{8, 8}, std::array<std::size_t, 2>{16, 12}}) {
    const std::vector<C> x = input(rows * cols, rows * 100 + cols);
    std::vector<C> y = x;
    fft2(y, rows, cols, FftDirection::Forward);
    const std::vector<C> ref = naive(x, rows, cols);
    const double err = max_diff(y, ref);
    INFO(rows << " x " << cols << ": max error " << err);
    CHECK(err <= bound2(rows, cols, sum_abs(x)));
  }
}

TEST_CASE("fft: Parseval", "[fft]") {
  // sum |y|^2 = L sum |x|^2. With |y_hat - y| <= e per element (e = bound1): ||y_hat|^2 - |y|^2|
  // <= (2|y| + e) e, summed and with sum |y| <= sqrt(L) ||y||_2 = L ||x||_2 (Cauchy-Schwarz).
  for (const std::size_t n : {std::size_t{100}, std::size_t{1018}}) {
    const std::vector<C> x = input(n, 7 * n);
    std::vector<C> y = x;
    fft1(y, FftDirection::Forward);
    double sx = 0.0;
    double sy = 0.0;
    for (std::size_t k = 0; k < n; ++k) {
      sx += std::norm(x[k]);
      sy += std::norm(y[k]);
    }
    const double nl = static_cast<double>(n);
    const double e = bound1(n, sum_abs(x));
    const double tol = 2.0 * e * nl * std::sqrt(sx) + nl * e * e + 4.0 * nl * kU * nl * sx;
    INFO("L = " << n << ": |sum|y|^2 - L sum|x|^2| = " << std::abs(sy - nl * sx));
    CHECK(std::abs(sy - nl * sx) <= tol);
  }
}

TEST_CASE("fft: delta and constant", "[fft]") {
  for (const std::size_t n : {std::size_t{64}, std::size_t{100}, std::size_t{1018}}) {
    std::vector<C> d(n, C(0.0, 0.0));
    d[0] = C(1.0, 0.0);
    std::vector<C> y = d;
    fft1(y, FftDirection::Forward);
    CHECK(max_diff(y, std::vector<C>(n, C(1.0, 0.0))) <= bound1(n, 1.0));
    std::vector<C> c(n, C(1.0, 0.0));
    fft1(c, FftDirection::Forward);
    std::vector<C> expected(n, C(0.0, 0.0));
    expected[0] = C(static_cast<double>(n), 0.0);
    CHECK(max_diff(c, expected) <= bound1(n, static_cast<double>(n)));
  }
}

TEST_CASE("fft: the periodised Gaussian exp(-pi j^2 / L) is self-dual with factor sqrt(L)",
          "[fft]") {
  // Poisson summation (derived for #204): with g_j = sum_n exp(-pi (j - n L)^2 / L), the DFT is
  // sum_j g_j e^{-2 pi i j k / L} = sum_m f(m) e^{-2 pi i m k / L} = sqrt(L) sum_n
  // exp(-pi L (k/L - n)^2) = sqrt(L) g_k for f(x) = exp(-pi x^2 / L), whose Fourier transform is
  // sqrt(L) exp(-pi L xi^2). For j in [0, L) the terms with n outside {-1, 0, 1} are below
  // exp(-pi L), negligible for L >= 64.
  for (const std::size_t n : {std::size_t{64}, std::size_t{100}, std::size_t{128}}) {
    const double nl = static_cast<double>(n);
    std::vector<C> g(n);
    for (std::size_t j = 0; j < n; ++j) {
      double s = 0.0;
      for (const double shift : {-1.0, 0.0, 1.0}) {
        const double d = static_cast<double>(j) - shift * nl;
        s += std::exp(-std::numbers::pi * d * d / nl);
      }
      g[j] = C(s, 0.0);
    }
    std::vector<C> y = g;
    fft1(y, FftDirection::Forward);
    // Bound: the FFT, the rounding of g (exp and sum, <= 4 u per value) through the transform
    // and in sqrt(L) g_k.
    const double tol = bound1(n, sum_abs(g)) + 4.0 * kU * (sum_abs(g) + std::sqrt(nl) * max_abs(g));
    double err = 0.0;
    for (std::size_t k = 0; k < n; ++k) err = std::max(err, std::abs(y[k] - std::sqrt(nl) * g[k]));
    INFO("L = " << n << ": max error " << err << ", tolerance " << tol);
    CHECK(err <= tol);
  }
}

TEST_CASE("fft: 2D of an outer product is the outer product of the 1D transforms", "[fft]") {
  const std::size_t rows = 12;
  const std::size_t cols = 20;
  const std::vector<C> a = input(rows, 1);
  const std::vector<C> b = input(cols, 2);
  std::vector<C> ab(rows * cols);
  for (std::size_t r = 0; r < rows; ++r) {
    for (std::size_t c = 0; c < cols; ++c) ab[r * cols + c] = a[r] * b[c];
  }
  fft2(ab, rows, cols, FftDirection::Forward);
  std::vector<C> fa = a;
  std::vector<C> fb = b;
  fft1(fa, FftDirection::Forward);
  fft1(fb, FftDirection::Forward);
  // Content guard: the identity would pass the comparison below; fa must be the DFT of a.
  CHECK(max_diff(fa, naive(a, 1, rows)) <= bound1(rows, sum_abs(a)));
  const double scale = sum_abs(a) * sum_abs(b);
  const double tol = (c_fft(rows) + c_fft(cols) + c_fft(rows) + c_fft(cols) + 24.0) * kU * scale;
  double err = 0.0;
  for (std::size_t r = 0; r < rows; ++r) {
    for (std::size_t c = 0; c < cols; ++c) {
      err = std::max(err, std::abs(ab[r * cols + c] - fa[r] * fb[c]));
    }
  }
  INFO("max error " << err << ", tolerance " << tol);
  CHECK(err <= tol);
}

TEST_CASE("fft: bitwise the same on repetition and from parallel threads", "[fft]") {
  // No global state (ADR 0033, point 1): two threads transform their own copies at the same time
  // and get bitwise the serial result.
  const std::size_t rows = 128;
  const std::size_t cols = 100;
  const std::vector<C> x = input(rows * cols, 42);
  std::vector<C> serial = x;
  fft2(serial, rows, cols, FftDirection::Forward);
  // Content guard: the identity would be bitwise reproducible as well; bin (0, 0) is sum_j x_j.
  Sum dc;
  for (const C& z : x) dc.add(z);
  CHECK(std::abs(serial[0] - dc.value()) <= bound2(rows, cols, sum_abs(x)));
  std::vector<C> again = x;
  fft2(again, rows, cols, FftDirection::Forward);
  CHECK(std::memcmp(serial.data(), again.data(), serial.size() * sizeof(C)) == 0);
  std::vector<C> t1 = x;
  std::vector<C> t2 = x;
  std::thread a([&t1] { fft2(t1, 128, 100, FftDirection::Forward); });
  std::thread b([&t2] { fft2(t2, 128, 100, FftDirection::Forward); });
  a.join();
  b.join();
  CHECK(std::memcmp(serial.data(), t1.data(), serial.size() * sizeof(C)) == 0);
  CHECK(std::memcmp(serial.data(), t2.data(), serial.size() * sizeof(C)) == 0);
}

TEST_CASE("fft: wrong sizes are invalid_argument", "[fft]") {
  const std::vector<C> x = input(12, 3);
  std::vector<C> v = x;
  CHECK_THROWS_AS(fft2(v, 0, 12, FftDirection::Forward), std::invalid_argument);
  CHECK_THROWS_AS(fft2(v, 12, 0, FftDirection::Forward), std::invalid_argument);
  CHECK_THROWS_AS(fft2(v, 5, 3, FftDirection::Forward), std::invalid_argument);
  std::vector<C> empty;
  CHECK_THROWS_AS(fft1(empty, FftDirection::Forward), std::invalid_argument);
  // The rejected calls leave the data unchanged; a valid call transforms it.
  CHECK(std::memcmp(v.data(), x.data(), v.size() * sizeof(C)) == 0);
  CHECK_NOTHROW(fft2(v, 3, 4, FftDirection::Backward));
  // Content guard: bin (0, 0) of the backward transform is sum_j x_j as well.
  Sum dc;
  for (const C& z : x) dc.add(z);
  CHECK(std::abs(v[0] - dc.value()) <= bound2(3, 4, sum_abs(x)));
}
