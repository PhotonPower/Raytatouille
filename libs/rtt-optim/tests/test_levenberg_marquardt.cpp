// Levenberg-Marquardt core (ADR 0030; #166). Reference problems from K. Madsen, H. B. Nielsen,
// O. Tingleff, "Methods for Non-Linear Least Squares Problems", 2nd ed., IMM DTU 2004 (MNT),
// the damping update of H. B. Nielsen, IMM-REP-1999-05, eq. (2.5), and the bound
// transformations of the MINUIT User's Guide (2004), eqs. (1.1)-(1.6); see docs/quellen.md.
//
// Tolerances are derived beforehand (solution_tolerance below), not measured.

#include <oneapi/tbb/task_arena.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <catch2/catch_test_macros.hpp>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

#include "rtt/optim/levenberg_marquardt.hpp"

using rtt::optim::Bounds;
using rtt::optim::levenberg_marquardt;
using rtt::optim::LmOptions;
using rtt::optim::LmResult;
using rtt::optim::LmStatus;
using rtt::optim::Residuals;

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

std::vector<Bounds> free_bounds(std::size_t n) {
  return std::vector<Bounds>(n);
}

// Rosenbrock as least squares (MNT examples 3.13 with lambda = 0 and 3.16):
// f = (10 (x2 - x1^2), 1 - x1), minimum (1, 1) with F = 0.
bool rosenbrock(std::span<const double> x, std::span<double> f) {
  f[0] = 10.0 * (x[1] - x[0] * x[0]);
  f[1] = 1.0 - x[0];
  return true;
}

// Chained Rosenbrock with n variables, m = 2 (n - 1) residuals, minimum at all ones (F = 0):
// f_2i = 10 (x_i+1 - x_i^2), f_2i+1 = 1 - x_i.
Residuals chained_rosenbrock() {
  return [](std::span<const double> x, std::span<double> f) {
    for (std::size_t i = 0; i + 1 < x.size(); ++i) {
      f[2 * i] = 10.0 * (x[i + 1] - x[i] * x[i]);
      f[2 * i + 1] = 1.0 - x[i];
    }
    return true;
  };
}

// The linear problem f = b - A x (MNT example 3.1); inconsistent, so F(x*) > 0.
constexpr std::array<std::array<double, 2>, 4> kLinearA{
    {{1.0, 2.0}, {3.0, 4.0}, {5.0, 7.0}, {1.0, -1.0}}};
constexpr std::array<double, 4> kLinearB{1.0, 2.0, 3.0, 4.0};
bool linear(std::span<const double> x, std::span<double> f) {
  for (std::size_t i = 0; i < 4; ++i) {
    f[i] = kLinearB[i] - (kLinearA[i][0] * x[0] + kLinearA[i][1] * x[1]);
  }
  return true;
}

// Least-squares solution of the linear problem from the normal equations A^T A x = A^T b
// (MNT eq. (3.5)), solved by Cramer's rule, independently of the solver.
std::array<double, 2> linear_solution() {
  double n00 = 0.0, n01 = 0.0, n11 = 0.0, r0 = 0.0, r1 = 0.0;
  for (std::size_t i = 0; i < 4; ++i) {
    n00 += kLinearA[i][0] * kLinearA[i][0];
    n01 += kLinearA[i][0] * kLinearA[i][1];
    n11 += kLinearA[i][1] * kLinearA[i][1];
    r0 += kLinearA[i][0] * kLinearB[i];
    r1 += kLinearA[i][1] * kLinearB[i];
  }
  const double det = n00 * n11 - n01 * n01;
  return {(r0 * n11 - n01 * r1) / det, (n00 * r1 - n01 * r0) / det};
}

// Dense matrix with rows x cols entries, row-major (a Jacobian with 1 or 2 columns).
struct Matrix {
  std::size_t rows = 0;
  std::size_t cols = 0;
  std::vector<double> values;
  [[nodiscard]] double at(std::size_t i, std::size_t j) const { return values[i * cols + j]; }
};

std::vector<double> column_norms(const Matrix& j) {
  std::vector<double> d(j.cols, 0.0);
  for (std::size_t c = 0; c < j.cols; ++c) {
    for (std::size_t i = 0; i < j.rows; ++i) d[c] += j.at(i, c) * j.at(i, c);
    d[c] = std::sqrt(d[c]);
  }
  return d;
}

// A priori bound on the error of the solution when the run stops with the step test
// (ADR 0030 point 8, MNT eq. (3.15b)): ||D h|| <= xtol (||D theta|| + xtol), h not applied.
// Near a solution the damped step in z = D theta is h_z = -(A_z + mu I)^-1 A_z e_z for the
// remaining error e_z (linearized, A_z = J_z^T J_z), so ||e_z|| <= ||h_z|| (1 + mu/lambda_min)
// with lambda_min the smallest eigenvalue of A_z. The test requires mu <= lambda_min/2 at the
// last solve, so ||e_z|| <= 1.5 ||h_z||; the bound uses the factor 2, the rest covers the
// nonlinear terms O(||e||^2) (relative 1e-8 here). In theta: ||e|| <= ||e_z|| / min_j d_j.
// d: column norms of the Jacobian (in theta) at the start; j_star: Jacobian at the solution
// (1 or 2 columns, so lambda_min has a closed form); theta_star: the solution in theta.
// Smallest eigenvalue of A_z = (J D^-1)^T (J D^-1) for J with 1 or 2 columns (closed form).
double scaled_lambda_min(const std::vector<double>& d, const Matrix& j_star) {
  const std::size_t n = j_star.cols;
  REQUIRE((n == 1 || n == 2));
  std::array<std::array<double, 2>, 2> a{};
  for (std::size_t r = 0; r < n; ++r) {
    for (std::size_t c = 0; c < n; ++c) {
      for (std::size_t i = 0; i < j_star.rows; ++i) {
        a[r][c] += j_star.at(i, r) / d[r] * (j_star.at(i, c) / d[c]);
      }
    }
  }
  return n == 1 ? a[0][0]
                : 0.5 * (a[0][0] + a[1][1]) - std::hypot(0.5 * (a[0][0] - a[1][1]), a[0][1]);
}

double solution_tolerance(const std::vector<double>& d,
                          const Matrix& j_star,
                          const std::vector<double>& theta_star,
                          double mu_last,
                          double xtol) {
  const std::size_t n = j_star.cols;
  const double lambda_min = scaled_lambda_min(d, j_star);
  INFO("lambda_min = " << lambda_min << ", mu at the last solve = " << mu_last);
  REQUIRE(mu_last <= lambda_min / 2.0);
  double scaled = 0.0;
  for (std::size_t j = 0; j < n; ++j) scaled += d[j] * theta_star[j] * d[j] * theta_star[j];
  return 2.0 * xtol * (std::sqrt(scaled) + xtol) / *std::min_element(d.begin(), d.end());
}

bool converged(LmStatus s) {
  return s == LmStatus::ConvergedStep || s == LmStatus::ConvergedGradient ||
         s == LmStatus::ConvergedMerit;
}

// Bitwise equality of two results (determinism).
void require_bitwise_equal(const LmResult& a, const LmResult& b) {
  REQUIRE(a.status == b.status);
  REQUIRE(a.p == b.p);
  REQUIRE(a.f == b.f);
  REQUIRE(std::bit_cast<std::uint64_t>(a.F) == std::bit_cast<std::uint64_t>(b.F));
  REQUIRE(a.iterations == b.iterations);
  REQUIRE(a.evaluations == b.evaluations);
  REQUIRE(a.failed_evaluations == b.failed_evaluations);
  REQUIRE(a.history.size() == b.history.size());
  for (std::size_t i = 0; i < a.history.size(); ++i) {
    const auto& x = a.history[i];
    const auto& y = b.history[i];
    REQUIRE(x.k == y.k);
    REQUIRE(std::bit_cast<std::uint64_t>(x.F) == std::bit_cast<std::uint64_t>(y.F));
    REQUIRE(std::bit_cast<std::uint64_t>(x.mu) == std::bit_cast<std::uint64_t>(y.mu));
    REQUIRE(std::bit_cast<std::uint64_t>(x.rho) == std::bit_cast<std::uint64_t>(y.rho));
    REQUIRE(std::bit_cast<std::uint64_t>(x.step_norm) == std::bit_cast<std::uint64_t>(y.step_norm));
    REQUIRE(x.accepted == y.accepted);
    REQUIRE(x.evaluations == y.evaluations);
  }
}

}  // namespace

TEST_CASE("LM: Rosenbrock converges to (1, 1) within the derived tolerance (T1)", "[optim][lm]") {
  const LmOptions options;
  const LmResult r = levenberg_marquardt(rosenbrock, 2, {-1.2, 1.0}, free_bounds(2), options);
  INFO("iterations = " << r.iterations);
  REQUIRE((r.status == LmStatus::ConvergedStep || r.status == LmStatus::ConvergedGradient));
  REQUIRE(r.iterations < options.max_iterations);
  // Analytic Jacobian [[-20 x1, 10], [-1, 0]]: at the start (-1.2, 1) and at the solution.
  const Matrix j0{2, 2, {24.0, 10.0, -1.0, 0.0}};
  const Matrix js{2, 2, {-20.0, 10.0, -1.0, 0.0}};
  const double tol =
      solution_tolerance(column_norms(j0), js, {1.0, 1.0}, r.history.back().mu, options.xtol);
  INFO("tolerance = " << tol);
  CHECK(std::hypot(r.p[0] - 1.0, r.p[1] - 1.0) <= tol);
  CHECK(r.changed == std::vector<std::uint8_t>{1, 1});
  CHECK(r.at_bound == std::vector<std::uint8_t>{0, 0});
  CHECK_FALSE(r.failed_variable.has_value());
}

TEST_CASE("LM: a linear problem gives rho = 1 and mu/3 per step, and the LS solution (T2)",
          "[optim][lm]") {
  // ftol = 0: only the step test stops the run, so the solution bound below applies.
  // tau = 1: strong damping at the start, so the first steps decrease F strongly.
  // Bound 1e-9 for rho - 1: f is linear, so the model L is exact and rho = 1 up to (a) the
  // rounding of F and F_new, about 2 m eps_M F = 1.8e-15 F for m = 4, divided by the predicted
  // decrease, and (b) the error of the central-difference Jacobian, eps_M |f| / h with
  // h = eps_M^(1/3) ~ 6e-6, i.e. a relative error of about 4e-11 |f| / |J| (|f| / |J| = O(1)
  // here). With a decrease of at least 1e-3 F per solve (checked below) (a) is below 2e-12;
  // the sum stays below 1e-10, so 1e-9 holds with a factor 10.
  LmOptions options;
  options.ftol = 0.0;
  options.tau = 1.0;
  const LmResult r = levenberg_marquardt(linear, 4, {0.0, 0.0}, free_bounds(2), options);
  REQUIRE(r.status == LmStatus::ConvergedStep);
  REQUIRE(r.history.size() >= 4);
  // Nielsen 1999 eq. (2.5) with beta = 2, gamma = 3, p = 3: rho = 1 gives
  // mu_new = mu max{1/3, 1 - (2 rho - 1)^3} = mu/3. Checked on the first three solves, where
  // the predicted decrease is far above the rounding of F.
  std::vector<double> f0(4);
  REQUIRE(linear(std::vector<double>{0.0, 0.0}, f0));
  double previous = 0.5 * (f0[0] * f0[0] + f0[1] * f0[1] + f0[2] * f0[2] + f0[3] * f0[3]);
  for (std::size_t i = 0; i < 3; ++i) {
    INFO("solve " << i + 1);
    REQUIRE(r.history[i].accepted);
    REQUIRE(r.history[i].F <= (1.0 - 1e-3) * previous);  // the decrease assumed above
    previous = r.history[i].F;
    CHECK(std::abs(r.history[i].rho - 1.0) <= 1e-9);
    CHECK(std::abs(r.history[i + 1].mu / r.history[i].mu - 1.0 / 3.0) <= 1e-9);
  }
  // Least-squares solution, independently of the solver.
  const std::array<double, 2> x_star = linear_solution();
  Matrix j{4, 2, {}};
  for (std::size_t i = 0; i < 4; ++i) {
    j.values.push_back(-kLinearA[i][0]);
    j.values.push_back(-kLinearA[i][1]);
  }
  const double tol = solution_tolerance(column_norms(j), j, {x_star[0], x_star[1]},
                                        r.history.back().mu, options.xtol);
  CHECK(std::hypot(r.p[0] - x_star[0], r.p[1] - x_star[1]) <= tol);
}

TEST_CASE("LM: rejected trials double nu, mu grows by 2 and then by 4 (T3)", "[optim][lm]") {
  // Calls in a fixed order: 1 start evaluation, 2n = 4 differences, then the trials one by
  // one. Calls 6 and 7, the first two trials, are invalid.
  std::atomic<int> calls{0};
  const Residuals f = [&calls](std::span<const double> x, std::span<double> out) {
    const int call = ++calls;
    if (call == 6 || call == 7) return false;
    return linear(x, out);
  };
  LmOptions options;
  options.tau = 1.0;
  const LmResult r = levenberg_marquardt(f, 4, {0.0, 0.0}, free_bounds(2), options);
  REQUIRE(r.history.size() >= 3);
  CHECK_FALSE(r.history[0].accepted);
  CHECK_FALSE(r.history[1].accepted);
  CHECK(r.history[2].accepted);
  CHECK(std::isnan(r.history[0].rho));
  CHECK(std::isnan(r.history[1].rho));
  // mu := mu nu; nu := 2 nu with nu starting at beta = 2 (exact in floating point).
  CHECK(r.history[1].mu == 2.0 * r.history[0].mu);
  CHECK(r.history[2].mu == 4.0 * r.history[1].mu);
  CHECK(r.failed_evaluations == 2);
  CHECK(converged(r.status));
}

TEST_CASE("LM: a minimum at a two-sided bound is reached from inside (T4)", "[optim][lm]") {
  // f = p - 2 with p in [0, 1]: the minimum of F is at the bound p = 1 (theta = pi/2).
  std::atomic<bool> outside{false};
  const Residuals f = [&outside](std::span<const double> p, std::span<double> out) {
    if (p[0] < 0.0 || p[0] > 1.0) outside = true;
    out[0] = p[0] - 2.0;
    return true;
  };
  const LmResult r = levenberg_marquardt(f, 1, {0.5}, {Bounds{0.0, 1.0}});
  INFO("status " << static_cast<int>(r.status) << ", p = " << r.p[0]);
  CHECK(converged(r.status));
  CHECK(1.0 - r.p[0] <= 1e-6);
  CHECK(r.at_bound[0] == 1);
  CHECK_FALSE(outside.load());
}

TEST_CASE("LM: minima at one-sided bounds are reached from inside (T5)", "[optim][lm]") {
  std::atomic<bool> outside{false};
  SECTION("lower bound: f = p + 1, p >= 0") {
    const Residuals f = [&outside](std::span<const double> p, std::span<double> out) {
      if (p[0] < 0.0) outside = true;
      out[0] = p[0] + 1.0;
      return true;
    };
    const LmResult r = levenberg_marquardt(f, 1, {2.0}, {Bounds{0.0, std::nullopt}});
    INFO("p = " << r.p[0]);
    CHECK(converged(r.status));
    CHECK(r.p[0] <= 1e-6);
    CHECK(r.at_bound[0] == 1);
  }
  SECTION("upper bound: f = p - 3, p <= 2") {
    const Residuals f = [&outside](std::span<const double> p, std::span<double> out) {
      if (p[0] > 2.0) outside = true;
      out[0] = p[0] - 3.0;
      return true;
    };
    const LmResult r = levenberg_marquardt(f, 1, {0.0}, {Bounds{std::nullopt, 2.0}});
    INFO("p = " << r.p[0]);
    CHECK(converged(r.status));
    CHECK(2.0 - r.p[0] <= 2e-6);
    CHECK(r.at_bound[0] == 1);
  }
  CHECK_FALSE(outside.load());
}

TEST_CASE("LM: an interior minimum next to a one-sided bound (T6)", "[optim][lm]") {
  // f = p - 0.5, p >= 0, start 3: theta = sqrt(p (p + 2)) (MINUIT eq. (1.3)), so
  // theta_0 = sqrt(15), theta* = sqrt(1.25); dp/dtheta = theta / sqrt(theta^2 + 1) (eq. (1.4)).
  const Residuals f = [](std::span<const double> p, std::span<double> out) {
    out[0] = p[0] - 0.5;
    return true;
  };
  const LmOptions options;
  const LmResult r = levenberg_marquardt(f, 1, {3.0}, {Bounds{0.0, std::nullopt}}, options);
  REQUIRE((r.status == LmStatus::ConvergedStep || r.status == LmStatus::ConvergedGradient));
  const auto dp = [](double t) { return t / std::sqrt(t * t + 1.0); };
  const double t0 = std::sqrt(15.0);
  const double ts = std::sqrt(1.25);
  const double tol_theta =
      solution_tolerance({dp(t0)}, Matrix{1, 1, {dp(ts)}}, {ts}, r.history.back().mu, options.xtol);
  CHECK(std::abs(r.p[0] - 0.5) <= dp(ts) * tol_theta * 1.01);  // 1 % for dp/dtheta over the step
  CHECK(r.at_bound[0] == 0);
}

TEST_CASE("LM: values of natural size 1e-4 next to a one-sided bound (T7)", "[optim][lm]") {
  // Everything lies in the quadratic part of MINUIT eq. (1.4): theta ~ sqrt(2 p) ~ 0.014.
  SECTION("interior minimum: f = 1e4 (p - 1e-4), p >= 0, start 5e-4") {
    const Residuals f = [](std::span<const double> p, std::span<double> out) {
      out[0] = 1e4 * (p[0] - 1e-4);
      return true;
    };
    const LmOptions options;
    const LmResult r = levenberg_marquardt(f, 1, {5e-4}, {Bounds{0.0, std::nullopt}}, options);
    REQUIRE((r.status == LmStatus::ConvergedStep || r.status == LmStatus::ConvergedGradient));
    const auto theta = [](double p) { return std::sqrt(p * (p + 2.0)); };
    const auto dp = [](double t) { return t / std::sqrt(t * t + 1.0); };
    const double t0 = theta(5e-4);
    const double ts = theta(1e-4);
    const double tol_theta = solution_tolerance({1e4 * dp(t0)}, Matrix{1, 1, {1e4 * dp(ts)}}, {ts},
                                                r.history.back().mu, options.xtol);
    CHECK(std::abs(r.p[0] - 1e-4) <= dp(ts) * tol_theta * 1.01);
    CHECK(r.at_bound[0] == 0);
  }
  SECTION("minimum at the bound: f = 1e4 (p + 1e-4), p >= 0, start 5e-4") {
    const Residuals f = [](std::span<const double> p, std::span<double> out) {
      out[0] = 1e4 * (p[0] + 1e-4);
      return true;
    };
    const LmResult r = levenberg_marquardt(f, 1, {5e-4}, {Bounds{0.0, std::nullopt}});
    INFO("p = " << r.p[0]);
    CHECK(converged(r.status));
    CHECK(r.p[0] <= 1e-6);
    CHECK(r.at_bound[0] == 1);
  }
}

TEST_CASE("LM: results are bitwise independent of threads and RunControl (T8)", "[optim][lm]") {
  // n = 7 variables: the 14 Jacobian evaluations do not divide evenly over 4 threads.
  const std::size_t n = 7;
  std::vector<double> x0(n);
  for (std::size_t i = 0; i < n; ++i) x0[i] = i % 2 == 0 ? -1.2 : 1.0;
  const Residuals chained = chained_rosenbrock();
  const auto run = [&](int threads, bool controlled) {
    rtt::trace::RunControl control;
    std::atomic<std::size_t> reports{0};
    if (controlled) {
      control.progress = [&reports](const rtt::trace::Progress&) { ++reports; };
      control.min_interval = std::chrono::milliseconds(0);
    }
    oneapi::tbb::task_arena arena(threads);
    LmResult chain;
    LmResult plain;
    arena.execute([&] {
      chain = levenberg_marquardt(chained, 2 * (n - 1), x0, free_bounds(n), {}, control);
      plain = levenberg_marquardt(rosenbrock, 2, {-1.2, 1.0}, free_bounds(2), {}, control);
    });
    if (controlled) REQUIRE(reports.load() > 0);
    return std::pair{chain, plain};
  };
  const auto reference = run(1, false);
  REQUIRE(converged(reference.first.status));
  for (const int threads : {1, 4, oneapi::tbb::task_arena::automatic}) {
    for (const bool controlled : {false, true}) {
      INFO("threads " << threads << ", controlled " << controlled);
      const auto other = run(threads, controlled);
      require_bitwise_equal(reference.first, other.first);
      require_bitwise_equal(reference.second, other.second);
    }
  }
}

TEST_CASE("LM: cancellation ends with the last accepted state (T9)", "[optim][lm]") {
  SECTION("cancelled before the start") {
    rtt::trace::RunControl control;
    control.cancel = rtt::trace::CancelToken{};
    control.cancel->request_cancel();
    const std::vector<double> p0{-1.2, 1.0};
    const LmResult r = levenberg_marquardt(rosenbrock, 2, p0, free_bounds(2), {}, control);
    CHECK(r.status == LmStatus::Cancelled);
    CHECK(r.p == p0);
    CHECK(r.evaluations == 0);
    CHECK(r.f.empty());
    CHECK(std::isnan(r.F));
    CHECK(r.changed == std::vector<std::uint8_t>{0, 0});
  }
  SECTION("cancelled after the first accepted solve equals a run with that k_max") {
    // K = the first accepted solve of an uncontrolled run (the first Rosenbrock solves from
    // (-1.2, 1) with mu_0 = 1e-3 overshoot and are rejected). The callback cancels when
    // "optimize" reports done >= K; the run then stops before the next trial evaluation and
    // must hold the state after solve K, the same as a run with k_max = K.
    const LmResult full = levenberg_marquardt(rosenbrock, 2, {-1.2, 1.0}, free_bounds(2));
    std::size_t first = 0;
    while (first < full.history.size() && !full.history[first].accepted) ++first;
    REQUIRE(first < full.history.size());
    const auto k_accepted = static_cast<int>(first + 1);
    rtt::trace::RunControl control;
    control.cancel = rtt::trace::CancelToken{};
    control.min_interval = std::chrono::milliseconds(0);
    rtt::trace::CancelToken token = *control.cancel;
    control.progress = [token, k_accepted](const rtt::trace::Progress& progress) mutable {
      if (progress.stage == std::string_view("optimize") &&
          progress.done >= static_cast<std::size_t>(k_accepted)) {
        token.request_cancel();
      }
    };
    const LmResult r = levenberg_marquardt(rosenbrock, 2, {-1.2, 1.0}, free_bounds(2), {}, control);
    LmOptions limited;
    limited.max_iterations = k_accepted;
    const LmResult k = levenberg_marquardt(rosenbrock, 2, {-1.2, 1.0}, free_bounds(2), limited);
    REQUIRE(k.status == LmStatus::MaxIterations);
    REQUIRE(k.history.size() == first + 1);
    REQUIRE(k.history.back().accepted);
    CHECK(r.status == LmStatus::Cancelled);
    CHECK(r.p == k.p);
    CHECK(r.f == k.f);
    CHECK(r.evaluations == k.evaluations);
    CHECK(r.p != std::vector<double>{-1.2, 1.0});
  }
  SECTION("an exception of the progress callback propagates") {
    rtt::trace::RunControl control;
    control.progress = [](const rtt::trace::Progress&) { throw std::runtime_error("stop"); };
    CHECK_THROWS_AS(levenberg_marquardt(rosenbrock, 2, {-1.2, 1.0}, free_bounds(2), {}, control),
                    std::runtime_error);
  }
}

TEST_CASE("LM: invalid evaluations (T10)", "[optim][lm]") {
  SECTION("invalid at the start point is an input error") {
    const Residuals f = [](std::span<const double>, std::span<double>) { return false; };
    CHECK_THROWS_AS(levenberg_marquardt(f, 1, {1.0}, free_bounds(1)), std::invalid_argument);
    const Residuals g = [](std::span<const double>, std::span<double> out) {
      out[0] = kNaN;
      return true;
    };
    CHECK_THROWS_AS(levenberg_marquardt(g, 1, {1.0}, free_bounds(1)), std::invalid_argument);
  }
  const std::vector<double> p0{0.5, -0.25, 2.0};
  const auto quadratic = [](std::span<const double> p, std::span<double> out) {
    out[0] = p[0] - 1.0;
    out[1] = p[1] + 1.0;
    out[2] = p[2] * p[2] - 2.0;
  };
  SECTION("an invalid Jacobian column ends with Failed and the smallest index") {
    // Columns 1 and 2 are invalid (any move of p1 or p2); column 0 is valid.
    const Residuals f = [&](std::span<const double> p, std::span<double> out) {
      if (p[1] != p0[1] || p[2] != p0[2]) return false;
      quadratic(p, out);
      return true;
    };
    const LmResult r = levenberg_marquardt(f, 3, p0, free_bounds(3));
    CHECK(r.status == LmStatus::Failed);
    REQUIRE(r.failed_variable.has_value());
    CHECK(*r.failed_variable == 1);
    CHECK(r.p == p0);
    CHECK(r.iterations == 0);
  }
  SECTION("an exception in the residual function counts as invalid") {
    const Residuals f = [&](std::span<const double> p, std::span<double> out) {
      if (p[0] > p0[0]) throw std::runtime_error("invalid");
      quadratic(p, out);
      return true;
    };
    const LmResult r = levenberg_marquardt(f, 3, p0, free_bounds(3));
    CHECK(r.status == LmStatus::Failed);
    REQUIRE(r.failed_variable.has_value());
    CHECK(*r.failed_variable == 0);
    CHECK(r.p == p0);
  }
}

TEST_CASE("LM: without an accepted step the values stay bitwise (T11)", "[optim][lm]") {
  const std::vector<double> p0{-1.2, 1.0};
  SECTION("k_max = 0") {
    LmOptions options;
    options.max_iterations = 0;
    const LmResult r = levenberg_marquardt(rosenbrock, 2, p0, free_bounds(2), options);
    CHECK(r.status == LmStatus::MaxIterations);
    CHECK(r.iterations == 0);
    CHECK(r.history.empty());
    CHECK(r.p == p0);
    std::vector<double> f0(2);
    REQUIRE(rosenbrock(p0, f0));
    CHECK(r.f == f0);
  }
  SECTION("every trial rejected") {
    // 1 start evaluation and 2n = 4 differences are valid, every later call is invalid.
    std::atomic<int> calls{0};
    const Residuals f = [&calls](std::span<const double> x, std::span<double> out) {
      if (++calls > 5) return false;
      return rosenbrock(x, out);
    };
    const LmResult r = levenberg_marquardt(f, 2, p0, free_bounds(2));
    CHECK(r.status == LmStatus::ConvergedStep);
    CHECK(r.p == p0);
    CHECK(r.changed == std::vector<std::uint8_t>{0, 0});
    CHECK(r.failed_evaluations == r.evaluations - 5);
  }
  SECTION("bounded variables keep their input value bitwise") {
    LmOptions options;
    options.max_iterations = 0;
    const std::vector<double> q0{0.3, 0.7};
    const LmResult r = levenberg_marquardt(rosenbrock, 2, q0,
                                           {Bounds{0.1, std::nullopt}, Bounds{-1.0, 1.0}}, options);
    CHECK(r.p == q0);
  }
}

TEST_CASE("LM: input errors are std::invalid_argument (T12)", "[optim][lm]") {
  const auto call = [](const Residuals& f, std::size_t m, std::vector<double> p0,
                       std::vector<Bounds> b, LmOptions o = {}) {
    return levenberg_marquardt(f, m, std::move(p0), std::move(b), o);
  };
  CHECK_THROWS_AS(call(Residuals{}, 2, {1.0, 1.0}, free_bounds(2)), std::invalid_argument);
  CHECK_THROWS_AS(call(rosenbrock, 0, {1.0, 1.0}, free_bounds(2)), std::invalid_argument);
  CHECK_THROWS_AS(call(rosenbrock, 2, {}, {}), std::invalid_argument);
  CHECK_THROWS_AS(call(rosenbrock, 2, {1.0, 1.0}, free_bounds(1)), std::invalid_argument);
  CHECK_THROWS_AS(call(rosenbrock, 2, {kNaN, 1.0}, free_bounds(2)), std::invalid_argument);
  CHECK_THROWS_AS(call(rosenbrock, 2, {1.0, 1.0}, {Bounds{2.0, 2.0}, Bounds{}}),
                  std::invalid_argument);
  CHECK_THROWS_AS(call(rosenbrock, 2, {1.0, 1.0}, {Bounds{kNaN, std::nullopt}, Bounds{}}),
                  std::invalid_argument);
  CHECK_THROWS_AS(call(rosenbrock, 2, {1.0, 1.0}, {Bounds{1.5, std::nullopt}, Bounds{}}),
                  std::invalid_argument);
  CHECK_THROWS_AS(call(rosenbrock, 2, {1.0, 1.0}, {Bounds{}, Bounds{std::nullopt, 0.5}}),
                  std::invalid_argument);
  const auto with = [&](auto&& set) {
    LmOptions o;
    set(o);
    return call(rosenbrock, 2, {-1.2, 1.0}, free_bounds(2), o);
  };
  CHECK_THROWS_AS(with([](LmOptions& o) { o.max_iterations = -1; }), std::invalid_argument);
  CHECK_THROWS_AS(with([](LmOptions& o) { o.ftol = -1.0; }), std::invalid_argument);
  CHECK_THROWS_AS(with([](LmOptions& o) { o.xtol = kNaN; }), std::invalid_argument);
  CHECK_THROWS_AS(with([](LmOptions& o) { o.gtol = -1.0; }), std::invalid_argument);
  CHECK_THROWS_AS(with([](LmOptions& o) { o.tau = 0.0; }), std::invalid_argument);
  CHECK_THROWS_AS(with([](LmOptions& o) { o.function_precision = 0.0; }), std::invalid_argument);
}

TEST_CASE("LM: zero residuals at the start end without differences", "[optim][lm]") {
  // ADR 0030 point 9: f(theta_0) = 0 gives g = 0 for every J; the run ends at the start.
  const Residuals zero = [](std::span<const double>, std::span<double> out) {
    for (double& v : out) v = 0.0;
    return true;
  };
  const LmResult r = levenberg_marquardt(zero, 3, {1.0, -2.0}, free_bounds(2));
  CHECK(r.status == LmStatus::ConvergedGradient);
  CHECK(r.evaluations == 1);
  CHECK(r.iterations == 0);
  CHECK(r.history.empty());
  CHECK(r.p == std::vector<double>{1.0, -2.0});
  CHECK(r.F == 0.0);
}

TEST_CASE("LM: a problem with F* > 0 stops by the merit test (MINPACK-1 F-convergence)",
          "[optim][lm]") {
  // The linear problem with default options (ftol = sqrt(eps_M)). A priori bound: f is linear,
  // so F(theta) = F* + 1/2 e_z^T A_z e_z with the error e_z = D (theta - theta*), and the damped
  // step h_z = -(A_z + mu I)^-1 A_z e_z predicts, in the eigenbasis of A_z (eigenvalues l),
  // L(0) - L(h) = 1/2 sum l^2 e^2 (l + 2 mu) / (l + mu)^2 >= (1 - mu^2/(l + mu)^2) (F - F*).
  // With mu <= lambda_min (checked) the factor is >= 3/4. The test stops when
  // L(0) - L(h) <= ftol F_old, so F_old - F* <= 4/3 ftol F_old and
  // ||e_z|| <= sqrt(8/3 ftol F_old / lambda_min) for the state before that solve; the result is
  // that state or an accepted step from it, whose error is smaller for a linear problem.
  // In theta: ||e|| <= ||e_z|| / min d.
  const LmOptions options;
  const LmResult r = levenberg_marquardt(linear, 4, {0.0, 0.0}, free_bounds(2), options);
  REQUIRE(r.status == LmStatus::ConvergedMerit);
  REQUIRE(r.history.size() >= 2);
  const double f_old = r.history[r.history.size() - 2].F;
  Matrix j{4, 2, {}};
  for (std::size_t i = 0; i < 4; ++i) {
    j.values.push_back(-kLinearA[i][0]);
    j.values.push_back(-kLinearA[i][1]);
  }
  const std::vector<double> d = column_norms(j);
  const double lambda_min = scaled_lambda_min(d, j);
  INFO("lambda_min = " << lambda_min << ", mu = " << r.history.back().mu);
  REQUIRE(r.history.back().mu <= lambda_min);
  const double tol = std::sqrt(8.0 / 3.0 * options.ftol * f_old / lambda_min) /
                     *std::min_element(d.begin(), d.end());
  const std::array<double, 2> x_star = linear_solution();
  INFO("tolerance = " << tol);
  CHECK(std::hypot(r.p[0] - x_star[0], r.p[1] - x_star[1]) <= tol);
}

TEST_CASE("LM: the merit test may end the run on a rejected trial", "[optim][lm]") {
  // Reference run: ends with ConvergedMerit after solve K, whose trial was the only evaluation
  // of that solve (call number E + 1 with E the evaluations after solve K - 1). The second run
  // returns at that call the residuals scaled to F_new = F_old (1 + 1e-9): rho < 0, so the
  // trial is rejected, but |F_old - F_new| / F_old = 1e-9 <= ftol and the predicted decrease
  // is that of the reference run, so the merit test holds on the rejected trial. The result
  // must be the last accepted state, the same as a run with k_max = K - 1.
  const LmResult ref = levenberg_marquardt(linear, 4, {0.0, 0.0}, free_bounds(2));
  REQUIRE(ref.status == LmStatus::ConvergedMerit);
  const std::size_t k = ref.history.size();
  REQUIRE(k >= 2);
  REQUIRE(ref.history[k - 1].accepted);
  const std::size_t e = ref.history[k - 2].evaluations;
  REQUIRE(ref.history[k - 1].evaluations == e + 1);
  const double f_old = ref.history[k - 2].F;
  std::atomic<std::size_t> calls{0};
  const Residuals f = [&](std::span<const double> x, std::span<double> out) {
    const std::size_t call = ++calls;
    if (!linear(x, out)) return false;
    if (call == e + 1) {
      double half = 0.0;
      for (const double v : out) half += 0.5 * v * v;
      const double s = std::sqrt(f_old * (1.0 + 1e-9) / half);
      for (double& v : out) v *= s;
    }
    return true;
  };
  const LmResult r = levenberg_marquardt(f, 4, {0.0, 0.0}, free_bounds(2));
  LmOptions before;
  before.max_iterations = static_cast<int>(k - 1);
  const LmResult b = levenberg_marquardt(linear, 4, {0.0, 0.0}, free_bounds(2), before);
  CHECK(r.status == LmStatus::ConvergedMerit);
  REQUIRE(r.history.size() == k);
  CHECK_FALSE(r.history.back().accepted);
  CHECK(r.history.back().rho < 0.0);
  CHECK(r.p == b.p);
  CHECK(r.f == b.f);
}

TEST_CASE("LM: a callback error is not hidden by a simultaneous cancellation", "[optim][lm]") {
  rtt::trace::RunControl control;
  control.cancel = rtt::trace::CancelToken{};
  control.min_interval = std::chrono::milliseconds(0);
  rtt::trace::CancelToken token = *control.cancel;
  control.progress = [token](const rtt::trace::Progress& progress) mutable {
    if (progress.stage == std::string_view("jacobian") && progress.done >= 1) {
      token.request_cancel();
      throw std::runtime_error("callback");
    }
  };
  CHECK_THROWS_AS(levenberg_marquardt(rosenbrock, 2, {-1.2, 1.0}, free_bounds(2), {}, control),
                  std::runtime_error);
}

TEST_CASE("LM: a start value exactly on a bound gives a zero column", "[optim][lm]") {
  // ADR 0030 point 5: at the bound theta = 0 (one-sided) or pi/2 (two-sided), p is symmetric
  // about it, the central difference is zero, g = 0, and the variable does not move.
  const Residuals f = [](std::span<const double> p, std::span<double> out) {
    out[0] = p[0] - 1.0;  // the minimum p = 1 lies inside
    return true;
  };
  SECTION("one-sided") {
    const LmResult r = levenberg_marquardt(f, 1, {0.0}, {Bounds{0.0, std::nullopt}});
    CHECK(r.status == LmStatus::ConvergedGradient);
    CHECK(r.p == std::vector<double>{0.0});
    CHECK(r.changed == std::vector<std::uint8_t>{0});
    CHECK(r.at_bound == std::vector<std::uint8_t>{1});
  }
  SECTION("two-sided") {
    const LmResult r = levenberg_marquardt(f, 1, {3.0}, {Bounds{-2.0, 3.0}});
    CHECK(r.status == LmStatus::ConvergedGradient);
    CHECK(r.p == std::vector<double>{3.0});
    CHECK(r.changed == std::vector<std::uint8_t>{0});
    CHECK(r.at_bound == std::vector<std::uint8_t>{1});
  }
}
