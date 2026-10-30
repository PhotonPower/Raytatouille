// The continuation fallback of real aiming (#287, ADR 0007 addendum) driven by a synthetic linear
// model instead of a traced system: hit(a, b) = A (a, b) + c on the stop, so the Newton run
// solves A (a, b) = g - c in one step and the Jacobian is A everywhere. The solution of the
// target s (tx, ty) moves on a straight line from x0 = A^-1 (-c) (chief ray, s = 0) to
// x1 = A^-1 ((tx, ty) - c) (s = 1). No real system was found where the bundle reserve decides
// (scan in #287), so the rule is tested here on both sides of 3 r_ep.

#include <algorithm>
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <optional>
#include <utility>
#include <vector>

#include "aim_continuation.hpp"
#include "rtt/trace/sources.hpp"

using rtt::trace::detail::continuation;
using rtt::trace::detail::kAimContinuationMinStep;
using rtt::trace::detail::kBundleReserve;
using rtt::trace::detail::NewtonRun;

namespace {

constexpr double kREp = 2.0;  // mm, paraxial EP radius

/// hit(a, b) = A (a, b) + c, A = {a11, a12, a21, a22}. The Newton run of the model fails for
/// targets beyond `fails_above` (in s) or more than `max_step` beyond the last accepted s, to
/// exercise the step rule. Like aim()'s run it takes no step (and computes no Jacobian) if it
/// starts on the target.
struct Linear {
  double a11 = 1.5, a12 = 0.25, a21 = -0.5, a22 = 2.0;
  double cx = 0.0, cy = 0.0;
  double tx = 0.0, ty = 0.0;
  double fails_above = 2.0;
  double max_step = 2.0;
  double last_ok = 0.0;                           ///< s of the last converged run
  std::vector<double> targets;                    ///< s of every Newton run, in call order
  std::vector<std::pair<double, double>> starts;  ///< start (u, v) of every Newton run
  int converged = 0;                              ///< converged runs
  int jacobian_calls = 0;                         ///< calls of jacobian_at

  [[nodiscard]] double det() const { return a11 * a22 - a12 * a21; }

  /// A^-1 (gx - cx, gy - cy).
  [[nodiscard]] std::pair<double, double> solve(double gx, double gy) const {
    const double rx = gx - cx;
    const double ry = gy - cy;
    return {(a22 * rx - a12 * ry) / det(), (-a21 * rx + a11 * ry) / det()};
  }

  /// Chooses c so that the chief-ray solution x0 lies at `x0`, i.e. c = -A x0.
  void put_chief_at(double x, double y) {
    cx = -(a11 * x + a12 * y);
    cy = -(a21 * x + a22 * y);
  }

  /// Chooses (tx, ty) so that the solution at s = 1 lies at `x1`, i.e. t = A x1 + c.
  void put_target_at(double x, double y) {
    tx = a11 * x + a12 * y + cx;
    ty = a21 * x + a22 * y + cy;
  }

  /// s of the target (gx, gy) = s (tx, ty).
  [[nodiscard]] double s_of(double gx, double gy) const {
    if (std::abs(tx) >= std::abs(ty)) return tx == 0.0 ? 0.0 : gx / tx;
    return gy / ty;
  }
};

std::optional<rtt::trace::detail::ContinuationResult> run(Linear& m, bool infinite) {
  const auto newton = [&m](double& u, double& v, double gx, double gy, NewtonRun& r) -> bool {
    const double s = m.s_of(gx, gy);
    m.targets.push_back(s);
    m.starts.emplace_back(u, v);
    if (s > m.fails_above || s - m.last_ok > m.max_step) return false;
    const double rx = m.a11 * u + m.a12 * v + m.cx - gx;
    const double ry = m.a21 * u + m.a22 * v + m.cy - gy;
    r.residual = std::hypot(rx, ry);
    if (r.residual >= rtt::trace::kAimTolerance) {
      const auto [x, y] = m.solve(gx, gy);
      u = x;
      v = y;
      r.steps = 1;
      r.residual = 0.0;
      r.jacobian = std::array<double, 4>{m.a11, m.a12, m.a21, m.a22};
    }
    m.last_ok = s;
    ++m.converged;
    return true;
  };
  const auto jacobian_at = [&m](double, double, double,
                                double) -> std::optional<std::array<double, 4>> {
    ++m.jacobian_calls;
    return std::array<double, 4>{m.a11, m.a12, m.a21, m.a22};
  };
  return continuation(newton, jacobian_at, infinite, kREp, m.tx, m.ty);
}
}  // namespace

TEST_CASE("#287 continuation: a linear model is followed in the documented steps",
          "[sources][aiming][continuation]") {
  // Every step is accepted: s = 1/16, then the step doubles, 3/16, 7/16, 15/16, and the last one
  // is cut to end exactly at s = 1. For a linear H with the exact Jacobian the Euler predictor
  // J^-1 (s_next - s) (tx, ty) lands on the solution, so every run after the chief ray starts on
  // its target and takes no step: one Newton step in all, and the Jacobian comes from
  // jacobian_at() after each of the five accepted steps.
  Linear m;
  m.put_chief_at(0.3 * kREp, -0.1 * kREp);
  m.put_target_at(1.2 * kREp, 0.8 * kREp);
  const auto r = run(m, true);
  REQUIRE(r.has_value());
  const auto result = *r;  // NOLINT(bugprone-unchecked-optional-access): REQUIRE above
  const std::vector<double> expected{0.0, 1.0 / 16, 3.0 / 16, 7.0 / 16, 15.0 / 16, 1.0};
  REQUIRE(m.targets.size() == expected.size());
  REQUIRE(m.starts.size() == expected.size());
  for (std::size_t i = 0; i < expected.size(); ++i) {
    INFO("run " << i);
    CHECK(std::abs(m.targets[i] - expected[i]) < 1e-15);
    if (i == 0) continue;  // the chief ray starts at the EP centre
    const auto [x, y] = m.solve(expected[i] * m.tx, expected[i] * m.ty);
    CHECK(std::abs(m.starts[i].first - x) < 1e-12);
    CHECK(std::abs(m.starts[i].second - y) < 1e-12);
  }
  CHECK(result.steps == 1);
  CHECK(m.jacobian_calls == 5);
  const auto [x1, y1] = m.solve(m.tx, m.ty);
  CHECK(std::abs(result.a - x1) < 1e-12);
  CHECK(std::abs(result.b - y1) < 1e-12);
}

TEST_CASE("#287 continuation: below the smallest step it gives up",
          "[sources][aiming][continuation]") {
  // The Newton run fails for every target beyond s = 1/2: s = 1/16, 3/16, 7/16 are accepted, then
  // the step halves down to 1/16 (s = 1/2 accepted) and again from 1/8 down to 1/1024, the
  // smallest step that is still tried; half of it is below the limit, so the result is none.
  Linear m;
  m.put_target_at(kREp, 0.0);
  m.fails_above = 0.5;
  REQUIRE_FALSE(run(m, true).has_value());
  double smallest = 1.0;
  for (std::size_t i = 1; i < m.targets.size(); ++i) {
    const double step = m.targets[i] - 0.5;
    if (step > 0.0) smallest = std::min(smallest, step);
  }
  CHECK(smallest == kAimContinuationMinStep);
  CHECK(m.targets.back() == 0.5 + kAimContinuationMinStep);
}

TEST_CASE("#287 continuation: a step cut to the rest is not tried twice",
          "[sources][aiming][continuation]") {
  // Fails at s = 1 only once it is reached from 15/16: the step 1 (cut to the rest 1/16) fails,
  // then half the rest, 1/32, is tried (31/32), not the same step to 1 again (old rule: ds 1/2,
  // 1/4, 1/8 and 1/16 all cut to the same try at s = 1).
  Linear m;
  m.put_target_at(kREp, 0.0);
  m.fails_above = 0.99;
  REQUIRE_FALSE(run(m, true).has_value());
  REQUIRE(m.targets.size() >= 7);
  CHECK(m.targets[4] == 15.0 / 16);
  CHECK(m.targets[5] == 1.0);
  CHECK(m.targets[6] == 31.0 / 32);
  // Later s = 1 is tried again, but from a later accepted point; never twice in a row.
  CHECK(std::adjacent_find(m.targets.begin(), m.targets.end()) == m.targets.end());
}

TEST_CASE("#287 continuation: at most 64 accepted steps", "[sources][aiming][continuation]") {
  // Only steps of at most 1/1024 converge: every accepted step is 1/1024 (each after one failed
  // try of 1/512), so s = 1 would need 1024 of them; the fallback stops after 64.
  Linear m;
  m.put_target_at(kREp, 0.0);
  m.max_step = rtt::trace::detail::kAimContinuationMinStep;
  REQUIRE_FALSE(run(m, true).has_value());
  CHECK(m.converged == 1 + rtt::trace::detail::kMaxAimContinuationSteps);  // with the chief ray
  CHECK(m.last_ok == 64.0 / 1024);
}
TEST_CASE("#287 continuation: the bundle reserve 3 r_ep, on both sides",
          "[sources][aiming][continuation]") {
  REQUIRE(kBundleReserve == 3.0);
  SECTION("chief-ray solution at 2.9 r_ep: accepted") {
    Linear m;
    m.put_chief_at(0.0, 2.9 * kREp);
    m.put_target_at(0.0, 2.5 * kREp);
    CHECK(run(m, true).has_value());
  }
  SECTION("chief-ray solution at 3.1 r_ep: none, before any step") {
    Linear m;
    m.put_chief_at(0.0, 3.1 * kREp);
    m.put_target_at(0.0, 2.5 * kREp);
    CHECK_FALSE(run(m, true).has_value());
    CHECK(m.targets.size() == 1);  // only the chief-ray run
  }
  SECTION("solution at s = 1 at 2.9 r_ep (diagonal): accepted") {
    Linear m;
    const double d = 2.9 * kREp / std::sqrt(2.0);
    m.put_target_at(d, d);
    const auto r = run(m, true);
    REQUIRE(r.has_value());
    const auto result = *r;  // NOLINT(bugprone-unchecked-optional-access): REQUIRE above
    CHECK(std::abs(std::hypot(result.a, result.b) - 2.9 * kREp) < 1e-12);
  }
  SECTION("solution at s = 1 at 3.1 r_ep: none") {
    Linear m;
    const double d = 3.1 * kREp / std::sqrt(2.0);
    m.put_target_at(d, d);
    CHECK_FALSE(run(m, true).has_value());
  }
  SECTION("finite object: no reserve, 3.1 r_ep is accepted") {
    Linear m;
    m.put_chief_at(0.0, 3.1 * kREp);
    m.put_target_at(0.0, 3.5 * kREp);
    CHECK(run(m, false).has_value());
  }
}
