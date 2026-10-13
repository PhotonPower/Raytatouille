// Random perturbations of tolerance items (ADR 0036, point 5; #259): Philox-4x32-10, FNV-1a 64,
// the 53-bit uniform numbers, Box-Muller and the three distributions; independence of order,
// other items and threads; statistics of the values against their analytic moments.
//
// Statistics: one item, N = 100000 trials with a fixed seed; the bands are derived here before
// the run, not chosen after it (the test is deterministic, the band says how far a correct
// implementation may lie from the expected value):
// - mean: |mean| <= 4 sigma / sqrt(N) (NIST/SEMATECH e-Handbook, Sec. 1.3.5.2, confidence limits
//   for the mean, Y +- t s/sqrt(N)).
// - variance: |s^2 - sigma^2| <= 4 sqrt(Var(s^2)) with Var(s^2) = (mu4 - (N-3)/(N-1) sigma^4)/N
//   for the estimator with N - 1. This is E. W. Weisstein, "Sample Variance Distribution",
//   MathWorld, Eq. (22), Var(m2) = (N-1)[(N-1) mu4 - (N-3) mu2^2]/N^3 for m2 with 1/N, times
//   (N/(N-1))^2. The band assumes that s^2 is close to normal (central limit theorem, N large);
//   it is not a Chebyshev bound.
// - mu4 per distribution, in units of the same sigma as sigma^2:
//   - uniform on [-T, T]: sigma^2 = T^2/3, kurtosis 9/5, so mu4 = 9/5 sigma^4 (NIST e-Handbook,
//     Sec. 1.3.6.6.2);
//   - normal truncated at +-a sigma_p (a = 2, sigma_p = T/2), in units of sigma_p: with the
//     density phi(x) = exp(-x^2/2)/sqrt(2 pi), phi' = -x phi and Z = 2 Phi(a) - 1, integration by
//     parts gives int_{-a}^{a} x^2 phi = Z - 2 a phi(a) and int_{-a}^{a} x^4 phi =
//     -2 a^3 phi(a) + 3 int x^2 phi, so Var = (1 - 2 a phi(a)/Z) sigma_p^2 = 0.773741 sigma_p^2
//     (ADR 0036, point 5) and mu4 = [3 Z - 2 a^3 phi(a) - 6 a phi(a)]/Z sigma_p^4
//     = 1.416189 sigma_p^4 (kurtosis 2.3655). Z = erf(a/sqrt(2)) from Phi(x) =
//     erfc(-x/sqrt(2))/2, which follows from the definitions of erf and erfc (NIST DLMF, Eqs.
//     (7.2.1), (7.2.2)) with t = x/sqrt(2). The numbers are computed at run time.
// - a share p: the expected share lies in the Wilson interval (NIST e-Handbook, Sec. 7.2.4.1)
//   with z = z_{1 - alpha/2}, alpha = 1e-4: P(|X| > z) = erfc(z/sqrt(2)) for a standard normal X,
//   z solved at run time by bisection (z = 3.8906). With 95 % a correct implementation would
//   fail for one seed in 20; this level rejects a correct one for about one seed in 10^4.

#include <oneapi/tbb/parallel_for.h>
#include <oneapi/tbb/task_arena.h>

#include <algorithm>
#include <array>
#include <bit>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "random_detail.hpp"
#include "rtt/model/tolerance.hpp"
#include "rtt/tolerance/random.hpp"

using rtt::model::ToleranceDistribution;
using rtt::model::ToleranceItem;
using rtt::tolerance::fnv1a64;
using rtt::tolerance::perturbation;
using rtt::tolerance::perturbations;
using rtt::tolerance::philox4x32_10;
using rtt::tolerance::uniform53;

namespace {

bool same(double a, double b) {
  return std::bit_cast<std::uint64_t>(a) == std::bit_cast<std::uint64_t>(b);
}

ToleranceItem item(const std::string& id, double limit, ToleranceDistribution distribution) {
  ToleranceItem i;
  i.id = id;
  i.limit = limit;
  i.distribution = distribution;
  return i;
}

/// z with P(|X| > z) = alpha for a standard normal X: erfc(z/sqrt(2)) = alpha by bisection on
/// [0, 40] (erfc decreases strictly).
double two_sided_z(double alpha) {
  double lo = 0.0;
  double hi = 40.0;
  for (int i = 0; i < 200; ++i) {
    const double mid = 0.5 * (lo + hi);
    if (std::erfc(mid / std::numbers::sqrt2) > alpha) {
      lo = mid;
    } else {
      hi = mid;
    }
  }
  return 0.5 * (lo + hi);
}

/// Wilson interval of the share k/n (NIST e-Handbook, Sec. 7.2.4.1) with z > 0.
std::pair<double, double> wilson(std::size_t k, std::size_t n, double z) {
  const double p = static_cast<double>(k) / static_cast<double>(n);
  const double nn = static_cast<double>(n);
  const double centre = p + z * z / (2.0 * nn);
  const double half = z * std::sqrt(p * (1.0 - p) / nn + z * z / (4.0 * nn * nn));
  const double denominator = 1.0 + z * z / nn;
  return {(centre - half) / denominator, (centre + half) / denominator};
}

/// Moments of a sample: mean and the variance with N - 1, both in two passes.
struct Moments {
  double mean = 0.0;
  double variance = 0.0;
};

Moments moments(const std::vector<double>& x) {
  double sum = 0.0;
  for (const double v : x) sum += v;
  const double n = static_cast<double>(x.size());
  const double mean = sum / n;
  double squares = 0.0;
  for (const double v : x) squares += (v - mean) * (v - mean);
  return {mean, squares / (n - 1.0)};
}

/// Checks mean 0 and variance `variance` of the sample within the bands above; mu4 is the fourth
/// central moment of the distribution.
void check_moments(const std::vector<double>& x, double variance, double mu4) {
  const double n = static_cast<double>(x.size());
  const Moments m = moments(x);
  INFO("mean " << m.mean << ", variance " << m.variance << " against " << variance);
  CHECK(std::abs(m.mean) <= 4.0 * std::sqrt(variance / n));
  const double var_s2 = (mu4 - (n - 3.0) / (n - 1.0) * variance * variance) / n;
  REQUIRE(var_s2 > 0.0);
  CHECK(std::abs(m.variance - variance) <= 4.0 * std::sqrt(var_s2));
}

constexpr std::uint64_t kSeed = 20261013;
constexpr std::uint32_t kTrials = 100000;

std::vector<double> sample(const ToleranceItem& i) {
  std::vector<double> x(kTrials);
  for (std::uint32_t t = 0; t < kTrials; ++t) x[t] = perturbation(kSeed, t, i);
  return x;
}

}  // namespace

TEST_CASE("random: Philox-4x32-10 gives the known-answer vectors", "[tolerance][random]") {
  // Known-answer vectors "philox4x32 10" from Random123, tests/kat_vectors
  // (github.com/DEShawResearch/random123): counter, key, expected output (docs/quellen.md,
  // ADR 0036, point 5). License of Random123, file LICENSE, verbatim:
  //
  // Copyright 2010-2012, D. E. Shaw Research.
  // All rights reserved.
  //
  // Redistribution and use in source and binary forms, with or without
  // modification, are permitted provided that the following conditions are
  // met:
  //
  // * Redistributions of source code must retain the above copyright
  //   notice, this list of conditions, and the following disclaimer.
  //
  // * Redistributions in binary form must reproduce the above copyright
  //   notice, this list of conditions, and the following disclaimer in the
  //   documentation and/or other materials provided with the distribution.
  //
  // * Neither the name of D. E. Shaw Research nor the names of its
  //   contributors may be used to endorse or promote products derived from
  //   this software without specific prior written permission.
  //
  // THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  // "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  // LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  // A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
  // OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  // SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
  // LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  // DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  // THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  // (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  // OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
  struct Kat {
    rtt::tolerance::PhiloxCounter counter;
    rtt::tolerance::PhiloxKey key;
    std::array<std::uint32_t, 4> expected;
  };
  const std::array<Kat, 3> kats{{
      {{0x00000000, 0x00000000, 0x00000000, 0x00000000},
       {0x00000000, 0x00000000},
       {0x6627e8d5, 0xe169c58d, 0xbc57ac4c, 0x9b00dbd8}},
      {{0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff},
       {0xffffffff, 0xffffffff},
       {0x408f276d, 0x41c83b0e, 0xa20bc7c6, 0x6d5451fd}},
      {{0x243f6a88, 0x85a308d3, 0x13198a2e, 0x03707344},
       {0xa4093822, 0x299f31d0},
       {0xd16cfe09, 0x94fdcceb, 0x5001e420, 0x24126ea1}},
  }};
  for (const Kat& k : kats) CHECK(philox4x32_10(k.counter, k.key) == k.expected);
}

TEST_CASE("random: FNV-1a 64 gives the values of the specification", "[tolerance][random]") {
  // draft-eastlake-fnv-35, Sec. 8.3 (FNV test code), FNV64svalues: the four test strings without
  // terminating NUL; "initial teststring is null so initial result is offset_basis" (Sec. 8.3).
  CHECK(fnv1a64("") == 0xcbf29ce484222325ULL);
  CHECK(fnv1a64("a") == 0xaf63dc4c8601ec8cULL);
  CHECK(fnv1a64("foobar") == 0x85944171f73967e8ULL);
  // Bytes >= 0x80 (ids are UTF-8): with a sign-extended char the hash would be
  // 0xbdb17e70953fffa1 instead.
  constexpr std::string_view kHigh{"Hello!\x01\xFF\xED", 9};
  CHECK(fnv1a64(kHigh) == 0xbd51ea7094ee6fa1ULL);
}

TEST_CASE("random: 53-bit uniform numbers from two words", "[tolerance][random]") {
  CHECK(same(uniform53(0, 0), 0.0));
  CHECK(same(uniform53(0xffffffffU, 0xffffffffU), 1.0 - std::ldexp(1.0, -53)));
  // w0 is widened before the shift: its top bit becomes 2^52, i.e. 1/2.
  CHECK(same(uniform53(0x80000000U, 0), 0.5));
  CHECK(same(uniform53(1, 0), std::ldexp(1.0, -32)));
  // Only the upper 21 bits of w1 count.
  CHECK(same(uniform53(0, 0x800U), std::ldexp(1.0, -53)));
  CHECK(same(uniform53(0, 0x7ffU), 0.0));
}

TEST_CASE("random: the words of a call and Box-Muller", "[tolerance][random]") {
  // Counter (trial, lower and upper half of the id hash, call), key (lower, upper half of the
  // seed); the values are recomputed here with the same library functions (no folded constants,
  // AGENTS.md, #218).
  const std::uint64_t seed = 0x0012345687654321ULL;
  const std::uint32_t trial = 7;
  const ToleranceItem u = item("L1.decenter", 0.05, ToleranceDistribution::Uniform);
  const std::uint64_t h = fnv1a64(u.id);
  const rtt::tolerance::PhiloxKey key{static_cast<std::uint32_t>(seed),
                                      static_cast<std::uint32_t>(seed >> 32)};
  const auto call_of = [&](std::uint32_t t, std::uint32_t j) {
    return philox4x32_10({t, static_cast<std::uint32_t>(h), static_cast<std::uint32_t>(h >> 32), j},
                         key);
  };
  const auto w0 = call_of(trial, 0);
  const double u_a = uniform53(w0[0], w0[1]);
  const double u_b = uniform53(w0[2], w0[3]);
  CHECK(same(perturbation(seed, trial, u), 0.05 * (2.0 * u_a - 1.0)));

  ToleranceItem e = u;
  e.distribution = ToleranceDistribution::Endpoint;
  CHECK(same(perturbation(seed, trial, e), u_a < 0.5 ? -0.05 : 0.05));

  // Normal: draw k from call floor(k/2), z1 for even and z2 for odd k; the first |z| <= 2.
  ToleranceItem n = u;
  n.distribution = ToleranceDistribution::Normal;
  // Reference of trial t: the value and the accepted k (-1 if all draws are rejected).
  const auto reference = [&](std::uint32_t t) {
    for (int k = 0; k < rtt::tolerance::kMaxNormalDraws; ++k) {
      const auto w = call_of(t, static_cast<std::uint32_t>(k / 2));
      const auto [z1, z2] =
          rtt::tolerance::detail::box_muller(uniform53(w[0], w[1]), uniform53(w[2], w[3]));
      const double z = k % 2 == 0 ? z1 : z2;
      if (std::abs(z) <= 2.0) return std::pair<double, int>{0.025 * z, k};
    }
    return std::pair<double, int>{std::numeric_limits<double>::quiet_NaN(), -1};
  };
  const auto [expected, k_trial] = reference(trial);
  REQUIRE(k_trial >= 0);
  CHECK(same(perturbation(seed, trial, n), expected));
  // Most trials accept at k = 0 (P = 0.954), which covers neither z2 nor call floor(k/2) for
  // k >= 2. Search trials that accept at an odd k and at a k >= 2 (P = 0.0436 and 0.0021 per
  // trial; 20000 trials miss one of them with P of about 1e-18).
  std::uint32_t trial_odd = 0;
  std::uint32_t trial_late = 0;
  int k_odd = -1;
  int k_late = -1;
  for (std::uint32_t t = 0; t < 20000 && (k_odd < 0 || k_late < 0); ++t) {
    const int k = reference(t).second;
    if (k_odd < 0 && k % 2 == 1) {
      k_odd = k;
      trial_odd = t;
    }
    if (k_late < 0 && k >= 2) {
      k_late = k;
      trial_late = t;
    }
  }
  INFO("odd k = " << k_odd << " at trial " << trial_odd << ", k >= 2: " << k_late << " at trial "
                  << trial_late);
  REQUIRE(k_odd >= 1);
  REQUIRE(k_late >= 2);
  CHECK(same(perturbation(seed, trial_odd, n), reference(trial_odd).first));
  CHECK(same(perturbation(seed, trial_late, n), reference(trial_late).first));

  // Box-Muller (MathWorld, Eqs. (1), (2)) with x1 = 1 - u_a and x2 = u_b.
  const auto [z1, z2] = rtt::tolerance::detail::box_muller(u_a, u_b);
  const double r = std::sqrt(-2.0 * std::log(1.0 - u_a));
  CHECK(same(z1, r * std::cos(2.0 * std::numbers::pi * u_b)));
  CHECK(same(z2, r * std::sin(2.0 * std::numbers::pi * u_b)));
  // u_a = 0 is allowed (x1 = 1, r = 0).
  const auto [zero1, zero2] = rtt::tolerance::detail::box_muller(0.0, 0.3);
  CHECK(zero1 == 0.0);
  CHECK(zero2 == 0.0);
}

TEST_CASE("random: rejection of the truncated normal distribution", "[tolerance][random]") {
  using rtt::tolerance::detail::truncated_normal;
  constexpr int kMax = rtt::tolerance::kMaxNormalDraws;
  int calls = 0;
  // The first |z| <= 2 is taken, |z| = 2 included; sigma = T/2.
  const std::array<double, 4> z{3.0, -2.5, -2.0, 1.0};
  CHECK(same(truncated_normal(
                 0.1, kMax,
                 [&](int k) {
                   ++calls;
                   return z[static_cast<std::size_t>(k)];
                 },
                 [] { return 0.0; }),
             0.05 * -2.0));
  CHECK(calls == 3);
  // All 64 draws rejected: T (2u - 1) with the fallback u.
  calls = 0;
  CHECK(same(truncated_normal(
                 0.1, kMax,
                 [&](int) {
                   ++calls;
                   return 2.5;
                 },
                 [] { return 0.75; }),
             0.1 * (2.0 * 0.75 - 1.0)));
  CHECK(calls == kMax);
  // The fallback is tested here, in the template, only: in perturbation() all 64 draws are
  // rejected with P = 0.0455^64, about 1e-86, so no seed reaches u_a of call 32 there.
}

TEST_CASE("random: values lie in [-T, T] and invalid limits are rejected", "[tolerance][random]") {
  for (const auto d : {ToleranceDistribution::Normal, ToleranceDistribution::Uniform,
                       ToleranceDistribution::Endpoint}) {
    const ToleranceItem i = item("x", 0.2, d);
    for (std::uint32_t t = 0; t < 2000; ++t) {
      const double v = perturbation(kSeed, t, i);
      REQUIRE(std::abs(v) <= 0.2);
      if (d == ToleranceDistribution::Endpoint) REQUIRE(std::abs(v) == 0.2);
    }
  }
  for (const double limit : {0.0, -1.0, std::numeric_limits<double>::infinity(),
                             std::numeric_limits<double>::quiet_NaN()}) {
    CHECK_THROWS_AS((void)perturbation(1, 0, item("x", limit, ToleranceDistribution::Normal)),
                    std::invalid_argument);
  }
  ToleranceItem bad = item("x", 1.0, ToleranceDistribution::Normal);
  bad.distribution = static_cast<ToleranceDistribution>(7);
  CHECK_THROWS_AS((void)perturbation(1, 0, bad), std::invalid_argument);
}

TEST_CASE("random: a value depends only on seed, trial and item", "[tolerance][random]") {
  const std::vector<ToleranceItem> items{item("A", 0.1, ToleranceDistribution::Normal),
                                         item("B", 0.2, ToleranceDistribution::Uniform),
                                         item("C", 1.0, ToleranceDistribution::Endpoint),
                                         item("D", 0.05, ToleranceDistribution::Normal)};
  const std::vector<ToleranceItem> fewer{items[3], items[0]};  // other order, B and C missing
  for (std::uint32_t t = 0; t < 50; ++t) {
    const std::vector<double> all = perturbations(kSeed, t, items);
    REQUIRE(all.size() == items.size());
    for (std::size_t i = 0; i < items.size(); ++i) {
      CHECK(same(all[i], perturbation(kSeed, t, items[i])));
    }
    const std::vector<double> some = perturbations(kSeed, t, fewer);
    CHECK(same(some[0], all[3]));
    CHECK(same(some[1], all[0]));
  }
  // Another seed, another trial or another id gives other values.
  const ToleranceItem& a = items[0];
  std::size_t differ_seed = 0;
  std::size_t differ_id = 0;
  ToleranceItem renamed = a;
  renamed.id = "A2";
  for (std::uint32_t t = 0; t < 50; ++t) {
    differ_seed += perturbation(kSeed, t, a) != perturbation(kSeed + 1, t, a) ? 1U : 0U;
    differ_id += perturbation(kSeed, t, a) != perturbation(kSeed, t, renamed) ? 1U : 0U;
  }
  CHECK(differ_seed == 50U);
  CHECK(differ_id == 50U);
  CHECK(perturbation(kSeed, 0, a) != perturbation(kSeed, 1, a));
  CHECK(perturbations(kSeed, 0, {}).empty());
}

TEST_CASE("random: the draws are bitwise the same for every number of threads",
          "[tolerance][random]") {
  const std::vector<ToleranceItem> items{item("A", 0.1, ToleranceDistribution::Normal),
                                         item("B", 0.2, ToleranceDistribution::Uniform),
                                         item("C", 1.0, ToleranceDistribution::Endpoint)};
  constexpr std::uint32_t kN = 2000;
  std::vector<double> serial;
  for (std::uint32_t t = 0; t < kN; ++t) {
    const std::vector<double> v = perturbations(kSeed, t, items);
    serial.insert(serial.end(), v.begin(), v.end());
  }
  const auto parallel = [&](int threads) {
    std::vector<double> out(serial.size());
    oneapi::tbb::task_arena arena(threads);
    arena.execute([&] {
      oneapi::tbb::parallel_for(std::uint32_t{0}, kN, [&](std::uint32_t t) {
        const std::vector<double> v = perturbations(kSeed, t, items);
        std::copy(v.begin(), v.end(), out.begin() + static_cast<std::ptrdiff_t>(t * items.size()));
      });
    });
    return out;
  };
  for (const int threads : {1, 4, oneapi::tbb::task_arena::automatic}) {
    INFO("threads " << threads);
    const std::vector<double> p = parallel(threads);
    bool all_same = true;
    for (std::size_t i = 0; i < serial.size(); ++i) all_same = all_same && same(p[i], serial[i]);
    CHECK(all_same);
  }
}

TEST_CASE("random: the statistics of the three distributions", "[tolerance][random]") {
  constexpr double kT = 0.3;
  SECTION("uniform: sigma = T/sqrt(3), mu4 = 9/5 sigma^4") {
    const double variance = kT * kT / 3.0;
    check_moments(sample(item("u", kT, ToleranceDistribution::Uniform)), variance,
                  9.0 / 5.0 * variance * variance);
  }
  SECTION("normal truncated at +-T: variance 0.773741 sigma^2, mu4 1.416189 sigma^4") {
    const double a = 2.0;
    const double sigma = kT / 2.0;
    const double phi = std::exp(-a * a / 2.0) / std::sqrt(2.0 * std::numbers::pi);
    const double z_mass = std::erf(a / std::numbers::sqrt2);  // 2 Phi(a) - 1
    const double var_ratio = 1.0 - 2.0 * a * phi / z_mass;
    const double mu4_ratio = (3.0 * z_mass - 2.0 * a * a * a * phi - 6.0 * a * phi) / z_mass;
    // Cross-check of the run-time values against ADR 0036, point 5 (6 digits there).
    REQUIRE(std::abs(var_ratio - 0.773741) <= 5e-7);
    REQUIRE(std::abs(mu4_ratio - 1.416189) <= 5e-7);
    const double s2 = sigma * sigma;
    check_moments(sample(item("n", kT, ToleranceDistribution::Normal)), var_ratio * s2,
                  mu4_ratio * s2 * s2);
  }
  SECTION("endpoint: only -T and +T, each with share 1/2") {
    const std::vector<double> x = sample(item("e", kT, ToleranceDistribution::Endpoint));
    const auto plus = static_cast<std::size_t>(std::count(x.begin(), x.end(), kT));
    const auto minus = static_cast<std::size_t>(std::count(x.begin(), x.end(), -kT));
    REQUIRE(plus + minus == x.size());
    const double z = two_sided_z(1e-4);
    REQUIRE(std::abs(z - 3.8906) <= 1e-4);  // cross-check of the bisection
    const auto [lower, upper] = wilson(plus, x.size(), z);
    CHECK(lower <= 0.5);
    CHECK(0.5 <= upper);
  }
}
