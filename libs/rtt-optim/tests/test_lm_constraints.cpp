// Equality constraints in levenberg_marquardt (ADR 0056; #196): the damped KKT step of S. Boyd,
// L. Vandenberghe, "Introduction to Applied Linear Algebra" (VMLS), eq. (16.4), and the
// projection with the least-norm correction, eq. (16.12), on problems with a closed-form solution.
// Tolerances are derived beforehand in the comments, not measured.

#include <atomic>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <stdexcept>
#include <vector>

#include "rtt/optim/bounds.hpp"
#include "rtt/optim/levenberg_marquardt.hpp"

using rtt::optim::Bounds;
using rtt::optim::levenberg_marquardt;
using rtt::optim::LmOptions;
using rtt::optim::LmResult;
using rtt::optim::LmStatus;
using rtt::optim::Residuals;

namespace {

// Minimize ||x - a||^2, a = (1, 2, 3), subject to x1 + x2 + x3 = 1 and x1 - x2 = 0 (two linear
// constraints, p = 2, n = 3). The KKT equations (VMLS (16.4) with A = I, b = a) give
// x* = a + C^T (C C^T)^-1 (d - C a); here C C^T = diag(3, 2), d - C a = (-5, 1), so
// x* = (-1/6, -1/6, 4/3).
Residuals linear_problem() {
  return [](std::span<const double> p, std::span<double> f) {
    f[0] = p[0] - 1.0;
    f[1] = p[1] - 2.0;
    f[2] = p[2] - 3.0;
    f[3] = p[0] + p[1] + p[2] - 1.0;  // g1
    f[4] = p[0] - p[1];               // g2
    return true;
  };
}

// Minimize (x - 0.75)^2 + (y - 1)^2 subject to x^2 + y^2 = 1: the point of the unit circle nearest
// to P = (0.75, 1), |P| = 1.25, is P / |P| = (0.6, 0.8).
Residuals circle_problem() {
  return [](std::span<const double> p, std::span<double> f) {
    f[0] = p[0] - 0.75;
    f[1] = p[1] - 1.0;
    f[2] = p[0] * p[0] + p[1] * p[1] - 1.0;
    return true;
  };
}

}  // namespace

// Tolerance of the solution (derivation for both cases below). The error e of the result in p has
// three parts:
// - step test: the run ends when ||D h|| <= xtol (||D theta|| + xtol) (ADR 0030, point 8); near the
//   solution the remaining error in z = D theta is at most twice the last step (MNT p. 25, as in
//   test_levenberg_marquardt.cpp), so e_step <= 2 xtol (||D theta|| + xtol) / min d_j;
// - central differences: each residual is evaluated with an absolute rounding delta_i, so each
//   Jacobian entry is off by up to delta_i / h with h >= eps_M^(1/3) = 6.06e-6 (ADR 0030, point
//   9); the stationary point of the KKT steps then lies off by about ||delta J|| ||f*|| / lambda
//   plus |w| ||delta C|| / lambda (lambda the curvature of the Lagrangian F + w^T g along the
//   constraints, w its multipliers);
// - acceptance: the gain ratio compares the Lagrangian of two states (ADR 0056, point 4), known
//   only to eta = eta_g + eta_F: eta_g = 2 |w| delta_g from the rounding of g, and
//   eta_F = sum_i |f_i| |delta_i - delta_new_i| <= 4 eps_M F from the rounding of f (each f_i to
//   eps_M |f_i|); a step from the error e with the contraction q gains 1/2 lambda (1 - q^2) e^2,
//   so steps are accepted reliably down to about e_noise = sqrt(2 eta / (lambda (1 - q^2))), with
//   q <= 1/2 here (the undamped rates are 0 and 1/4, and mu falls by a factor of at least 1/3
//   per good step, Nielsen eq. (2.5)). This is an estimate of the order of magnitude, not a
//   strict bound, so the tolerances below round it up generously.
// The last part dominates; the runs themselves end far below these bounds.

TEST_CASE("constraints: linear residuals and two linear constraints", "[optim][196]") {
  // Linear f and g: the KKT step with mu -> 0 is the exact solution, the projection reaches g = 0
  // in one correction. ftol = 0, so the run ends by the step test (xtol = 1e-14) and not at a
  // relative decrease of F. Tolerance (derivation above): lambda = 1 (F = 1/2 ||x - a||^2 has the
  // Hessian I, g is linear); w solves C^T w = -(x* - a) = (7/6, 13/6, 5/3), so w = (5/3, -1/2),
  // |w| = 1.74; delta_g <= 3 eps_M 4/3 = 8.9e-16 (three additions of values <= 4/3), so
  // eta_g = 3.1e-15; F* = 4.42, so eta_F <= 4 eps_M F* = 3.9e-15; eta = 7.0e-15 and
  // e_noise = sqrt(2 * 7.0e-15 / 0.75) = 1.4e-7. Differences: |f| <= 3, so
  // |delta J| <= eps_M 3 / h = 1.1e-10 per entry, ||delta J|| <= 3.3e-10, times ||f*|| = 2.97:
  // 9.8e-10; ||delta C|| <= sqrt(6) 8.9e-16 / h = 3.6e-10, times |w|: 6.3e-10. Step test: D =
  // column norms <= sqrt(3), ||D theta|| <= 3: 2e-14 * 3 = 6e-14. Sum 1.4e-7: 2e-7.
  LmOptions o;
  o.xtol = 1e-14;
  o.ftol = 0.0;
  o.constraint_tolerances = {1e-14, 1e-14};
  const double x1 = -1.0 / 6.0;
  const double x3 = 4.0 / 3.0;
  constexpr double kTol = 2e-7;
  SECTION("free variables") {
    const LmResult r =
        levenberg_marquardt(linear_problem(), 5, {0.0, 0.0, 0.0}, std::vector<Bounds>(3), o);
    INFO("status " << static_cast<int>(r.status) << ", iterations " << r.iterations);
    REQUIRE(r.status != LmStatus::Failed);
    CHECK(r.status == LmStatus::ConvergedStep);  // ftol = 0: the step test ends the run
    CHECK(std::abs(r.p[0] - x1) <= kTol);
    CHECK(std::abs(r.p[1] - x1) <= kTol);
    CHECK(std::abs(r.p[2] - x3) <= kTol);
    CHECK(std::abs(r.f[3]) <= 1e-14);  // feasible
    CHECK(std::abs(r.f[4]) <= 1e-14);
    // F is 1/2 ||f||^2 over the residuals of F only, not the constraints.
    CHECK(r.F == 0.5 * (r.f[0] * r.f[0] + r.f[1] * r.f[1] + r.f[2] * r.f[2]));
    CHECK(r.F > 1.0);  // content guard: the residual at the solution is not zero
    CHECK(!r.dependent_constraints);
    CHECK(!r.infeasible_constraint);
  }
  SECTION("bounds around the solution (MINUIT transform)") {
    std::vector<Bounds> bounds(3);
    bounds[2].min = 0.0;
    bounds[2].max = 2.0;
    const LmResult r = levenberg_marquardt(linear_problem(), 5, {0.0, 0.0, 1.0}, bounds, o);
    REQUIRE(r.status != LmStatus::Failed);
    // The transform p = a + (b - a)(sin theta + 1) / 2 has |dp/dtheta| <= (b - a)/2 = 1, so the
    // bound in theta also holds in p.
    CHECK(std::abs(r.p[0] - x1) <= kTol);
    CHECK(std::abs(r.p[1] - x1) <= kTol);
    CHECK(std::abs(r.p[2] - x3) <= kTol);
    CHECK(std::abs(r.f[3]) <= 1e-14);
    CHECK(std::abs(r.f[4]) <= 1e-14);
    CHECK(r.at_bound[2] == 0);
  }
}

TEST_CASE("constraints: nearest point of a circle from an infeasible start", "[optim][196]") {
  // The start P = (0.75, 1) has g = 0.5625 and is projected first, by Gauss-Newton with a new
  // Jacobian per correction (ADR 0056, point 5). Along the circle (angle phi) the Gauss-Newton
  // curvature of F is |df/dphi|^2 = 1 and that of the Lagrangian lambda = |P| = 1.25 at the
  // solution, so the undamped step maps the error e to (1 - 1.25) e (rate 1/4, MNT section 3.2).
  // Tolerance (derivation above): grad F + w grad g = 0 with grad F = f* = -0.25 P/|P| and
  // grad g = 2 x* gives w = 0.125; delta_g <= 3 eps_M = 6.7e-16 (x^2 + y^2 - 1 near 1), so
  // eta_g = 1.7e-16; f = (x - 0.75, y - 1) is exact near the solution (Sterbenz lemma: x in
  // [0.375, 1.5], y in [0.5, 2]), so eta_F = 0 and e_noise = sqrt(2 * 1.7e-16 / (1.25 * 0.75)) =
  // 1.9e-8. Differences: |f| <=
  // 0.3 near the solution: ||delta J|| <= 2 eps_M 0.3 / h = 2.2e-11, times ||f*|| /
  // lambda: 4.4e-12;
  // ||delta C|| <= 2 delta_g / h = 2.2e-10, times |w| / lambda: 2.2e-11. Step test: xtol = 1e-13,
  // D <= 2.3, ||D theta|| <= 3: 6e-13. Sum 1.9e-8: 2e-8.
  LmOptions o;
  o.xtol = 1e-13;
  o.ftol = 0.0;
  o.constraint_tolerances = {1e-14};
  const LmResult r =
      levenberg_marquardt(circle_problem(), 3, {0.75, 1.0}, std::vector<Bounds>(2), o);
  INFO("status " << static_cast<int>(r.status) << ", iterations " << r.iterations);
  REQUIRE(r.status != LmStatus::Failed);
  CHECK(r.status == LmStatus::ConvergedStep);  // ftol = 0: the step test ends the run
  CHECK(std::abs(r.p[0] - 0.6) <= 2e-8);
  CHECK(std::abs(r.p[1] - 0.8) <= 2e-8);
  CHECK(std::abs(r.f[2]) <= 1e-14);
  CHECK(r.iterations > 1);  // content guard: not solved by the projection of the start alone
}

TEST_CASE("constraints: as many constraints as variables (p = n)", "[optim][196]") {
  // x + y = 1 and x - y = 0.5 fix (0.75, 0.25); F = 1/2 ((x - 2)^2 + y^2) does not matter. The
  // start projection (Newton, exact for linear constraints) reaches the point in one correction;
  // the KKT step is then h = 0 (C invertible, g = 0), so the run ends by the step test with the
  // values reached, up to the rounding of one correction (a few eps_M): 1e-15.
  const Residuals residuals = [](std::span<const double> p, std::span<double> f) {
    f[0] = p[0] - 2.0;
    f[1] = p[1];
    f[2] = p[0] + p[1] - 1.0;
    f[3] = p[0] - p[1] - 0.5;
    return true;
  };
  LmOptions o;
  o.constraint_tolerances = {1e-14, 1e-14};
  const LmResult r = levenberg_marquardt(residuals, 4, {0.0, 0.0}, std::vector<Bounds>(2), o);
  INFO("status " << static_cast<int>(r.status) << ", iterations " << r.iterations);
  CHECK(r.status == LmStatus::ConvergedStep);
  CHECK(std::abs(r.p[0] - 0.75) <= 1e-15);
  CHECK(std::abs(r.p[1] - 0.25) <= 1e-15);
  CHECK(r.iterations == 1);
}

TEST_CASE("constraints: the constrained minimum on a bound", "[optim][196]") {
  // Minimize 1/2 ((x - 2)^2 + y^2) subject to x + y = 1 with x in [0, 1]: along the constraint
  // F = 1/2 ((x - 2)^2 + (1 - x)^2) is smallest at x = 1.5, outside the bound, so the solution is
  // the bound x = 1, y = 0, reached from inside (theta -> pi/2 of the MINUIT transform, ADR 0030,
  // point 5) as T4 of test_levenberg_marquardt.cpp; the same tolerance 1e-6 (near the bound
  // p = 1 - (b - a)/4 dtheta^2, so the step test in theta leaves 1 - p tiny).
  std::atomic<bool> outside{false};
  const Residuals residuals = [&outside](std::span<const double> p, std::span<double> f) {
    if (p[0] < 0.0 || p[0] > 1.0) outside = true;
    f[0] = p[0] - 2.0;
    f[1] = p[1];
    f[2] = p[0] + p[1] - 1.0;
    return true;
  };
  LmOptions o;
  o.ftol = 0.0;
  o.constraint_tolerances = {1e-14};
  std::vector<Bounds> bounds(2);
  bounds[0].min = 0.0;
  bounds[0].max = 1.0;
  const LmResult r = levenberg_marquardt(residuals, 3, {0.5, 0.5}, bounds, o);
  INFO("status " << static_cast<int>(r.status) << ", p = " << r.p[0] << ", " << r.p[1]);
  REQUIRE(r.status != LmStatus::Failed);
  CHECK(1.0 - r.p[0] <= 1e-6);
  CHECK(std::abs(r.p[1]) <= 1e-6);
  CHECK(std::abs(r.f[2]) <= 1e-14);
  CHECK(r.at_bound[0] == 1);
  CHECK_FALSE(outside.load());
}

TEST_CASE("constraints: a constraint satisfiable only on a bound (known limit)", "[optim][196]") {
  // g = x^2 - 0.25 with x in [0.5, 1] holds only at the bound x = 0.5. In the internal variable
  // the MINUIT transform makes g quadratic in theta there (dp/dtheta = 0 on the bound), so the
  // Newton projection of the start converges only linearly and does not reach the tolerance in
  // 10 corrections (ADR 0056, point 5): the run fails at the start with that constraint. This
  // states the known limit; such a condition belongs in the parameter table (a fixed value).
  const Residuals residuals = [](std::span<const double> p, std::span<double> f) {
    f[0] = p[0] - 2.0;
    f[1] = p[1];
    f[2] = p[0] * p[0] - 0.25;
    return true;
  };
  LmOptions o;
  o.constraint_tolerances = {1e-12};
  std::vector<Bounds> bounds(2);
  bounds[0].min = 0.5;
  bounds[0].max = 1.0;
  const LmResult r = levenberg_marquardt(residuals, 3, {0.75, 0.0}, bounds, o);
  CHECK(r.status == LmStatus::Failed);
  CHECK(r.infeasible_constraint == std::optional<std::size_t>(0));
  CHECK(r.p == std::vector<double>{0.75, 0.0});  // the start
  CHECK(r.evaluations > 10);                     // content guard: the projection ran
}

TEST_CASE("constraints: dependent constraints end the run as Failed", "[optim][196]") {
  // g1 = x + y - 1 and g2 = 2 g1: the rows of C are parallel, so C D^-2 C^T and the Schur
  // complement of the KKT matrix are singular (VMLS (16.5)).
  const Residuals residuals = [](std::span<const double> p, std::span<double> f) {
    f[0] = p[0] - 2.0;
    f[1] = p[1] - 1.0;
    f[2] = p[0] + p[1] - 1.0;
    f[3] = 2.0 * (p[0] + p[1] - 1.0);
    return true;
  };
  LmOptions o;
  o.constraint_tolerances = {1e-12, 1e-12};
  SECTION("infeasible start: in the projection") {
    const LmResult r = levenberg_marquardt(residuals, 4, {0.0, 0.0}, std::vector<Bounds>(2), o);
    CHECK(r.status == LmStatus::Failed);
    CHECK(r.dependent_constraints);
    CHECK(r.p == std::vector<double>{0.0, 0.0});  // the start
  }
  SECTION("feasible start: in the KKT step") {
    const LmResult r = levenberg_marquardt(residuals, 4, {0.5, 0.5}, std::vector<Bounds>(2), o);
    CHECK(r.status == LmStatus::Failed);
    CHECK(r.dependent_constraints);
    CHECK(r.iterations == 1);
  }
}

TEST_CASE("constraints: a constraint that cannot be satisfied", "[optim][196]") {
  // g1 = x - y is satisfiable, g2 = x^2 + y^2 + 1 >= 1 is not: the projection stops when the
  // violation no longer decreases, and the run fails at the start with the index of the largest
  // |g_i| / tolerance_i at the start (g2: 1.1 against g1: 0.2).
  const Residuals residuals = [](std::span<const double> p, std::span<double> f) {
    f[0] = p[0] - 2.0;
    f[1] = p[1] - 1.0;
    f[2] = p[0] - p[1];
    f[3] = p[0] * p[0] + p[1] * p[1] + 1.0;
    return true;
  };
  LmOptions o;
  o.constraint_tolerances = {1e-12, 1e-12};
  const LmResult r = levenberg_marquardt(residuals, 4, {0.3, 0.1}, std::vector<Bounds>(2), o);
  CHECK(r.status == LmStatus::Failed);
  CHECK(r.infeasible_constraint == std::optional<std::size_t>(1));
  CHECK(!r.dependent_constraints);
  CHECK(r.p == std::vector<double>{0.3, 0.1});  // the start
  CHECK(r.evaluations > 2);                     // content guard: the projection ran
}

TEST_CASE("constraints: invalid input", "[optim][196]") {
  const auto run = [](std::size_t m, std::vector<double> p0, std::vector<double> tolerances) {
    LmOptions o;
    o.constraint_tolerances = std::move(tolerances);
    const std::size_t n = p0.size();
    const Residuals residuals = [m](std::span<const double>, std::span<double> f) {
      for (std::size_t i = 0; i < m; ++i) f[i] = 1.0;
      return true;
    };
    return levenberg_marquardt(residuals, m, std::move(p0), std::vector<Bounds>(n), o);
  };
  CHECK_THROWS_AS(run(2, {0.0, 0.0}, {0.0}), std::invalid_argument);           // tolerance 0
  CHECK_THROWS_AS(run(2, {0.0, 0.0}, {-1e-12}), std::invalid_argument);        // negative
  CHECK_THROWS_AS(run(2, {0.0, 0.0}, {std::nan("")}), std::invalid_argument);  // not finite
  CHECK_THROWS_AS(run(2, {0.0, 0.0}, {1e-12, 1e-12}), std::invalid_argument);  // p >= m
  CHECK_THROWS_AS(run(3, {0.0}, {1e-12, 1e-12}), std::invalid_argument);       // p > n
}

TEST_CASE("constraints: on_accepted reports every accepted state", "[optim][196]") {
  // As without constraints (ADR 0030, addendum #378): one call per accepted solve, the last with
  // the values of the result; every reported state is feasible (it passed the projection).
  LmOptions o;
  o.constraint_tolerances = {1e-14};
  std::size_t calls = 0;
  std::vector<double> last;
  bool feasible = true;
  const rtt::optim::LmAcceptedCallback on = [&](const rtt::optim::LmIteration&,
                                                std::span<const double> p,
                                                std::span<const std::uint8_t>) {
    ++calls;
    last.assign(p.begin(), p.end());
    feasible = feasible && std::abs(p[0] * p[0] + p[1] * p[1] - 1.0) <= 1e-14;
  };
  const LmResult r =
      levenberg_marquardt(circle_problem(), 3, {0.75, 1.0}, std::vector<Bounds>(2), o, {}, on);
  std::size_t accepted = 0;
  for (const auto& it : r.history) accepted += it.accepted ? 1 : 0;
  CHECK(calls == accepted);
  CHECK(calls > 1);  // content guard
  CHECK(last == r.p);
  CHECK(feasible);
}
