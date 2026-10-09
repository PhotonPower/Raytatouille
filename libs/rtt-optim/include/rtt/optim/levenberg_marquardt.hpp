#pragma once

/// @file levenberg_marquardt.hpp
/// Damped Levenberg-Marquardt for nonlinear least squares (ADR 0030, points 1 and 5-11).
///
/// Minimizes F = 1/2 f^T f over n variables with optional bounds. The solver is Algorithm 3.16
/// of K. Madsen, H. B. Nielsen, O. Tingleff, "Methods for Non-Linear Least Squares Problems",
/// 2nd ed., IMM DTU 2004 (MNT), applied to the scaled internal variables z = D theta, with the
/// damping update of H. B. Nielsen, IMM-REP-1999-05, eq. (2.5), the column scaling and the
/// F-convergence test of MINPACK-1 (Moré, Garbow, Hillstrom, ANL-80-74, sections 2.5 and 2.3)
/// and the bound transformations of MINUIT (bounds.hpp). The Jacobian is a central difference
/// in theta, evaluated in parallel over the variables. Sources and conventions:
/// docs/quellen.md.
///
/// Results are bitwise reproducible: they do not depend on the number of threads, the
/// partitioning or the RunControl (ADR 0004, ADR 0030 point 11).

#include <cstddef>
#include <cstdint>
#include <functional>
#include <optional>
#include <span>
#include <vector>

#include "rtt/optim/bounds.hpp"
#include "rtt/trace/run_control.hpp"

namespace rtt::optim {

/// Options of the solver (ADR 0030, points 7-9). Not part of the model or the file.
struct LmOptions {
  /// k_max: maximum number of solves of the damped normal equations (MNT eq. (3.15c)). Every
  /// solve counts, accepted or not.
  int max_iterations = 100;
  /// Relative decrease of F, actual and predicted (MINPACK-1 F-convergence, section 2.3).
  /// Default sqrt(eps_M) with eps_M = DBL_EPSILON.
  double ftol = 1.4901161193847656e-8;
  /// Scaled step: ||D h|| <= xtol (||D theta|| + xtol) (MNT eq. (3.15b)). Default sqrt(eps_M).
  double xtol = 1.4901161193847656e-8;
  /// Scaled gradient: ||D^-1 g||_inf <= gtol (MNT eq. (3.15a)). Default 0: only g = 0.
  double gtol = 0.0;
  /// mu_0 = tau max_i (J_z^T J_z)_ii (MNT eq. (3.14)).
  double tau = 1e-3;
  /// eps_f: absolute noise of the residuals relative to the fixed scale ||f(theta_0)||_inf;
  /// sets the central-difference step h_j = eps_f^(1/3) max(|theta_j|, 1). Default eps_M.
  double function_precision = 2.220446049250313e-16;
};

/// How the run ended (ADR 0030, point 8).
enum class LmStatus : std::uint8_t {
  ConvergedGradient,  ///< ||D^-1 g||_inf <= gtol
  ConvergedStep,      ///< ||D h|| <= xtol (||D theta|| + xtol); h is not applied
  ConvergedMerit,     ///< small relative actual and predicted decrease of F
  MaxIterations,      ///< k >= k_max
  Cancelled,          ///< cancellation through the RunControl; last accepted state
  Failed,             ///< a Jacobian column without valid evaluations; last accepted state
};

/// One solve of the damped normal equations.
struct LmIteration {
  int k = 0;                    ///< number of this solve (1, 2, ...)
  double F = 0.0;               ///< F of the current (last accepted) state after this solve
  double mu = 0.0;              ///< damping mu used in this solve
  double rho = 0.0;             ///< gain ratio; NaN without a valid trial
  double step_norm = 0.0;       ///< ||D h||; NaN if the decomposition failed
  bool accepted = false;        ///< the trial theta + h was accepted
  std::size_t evaluations = 0;  ///< evaluations of f so far (cumulative)
};

/// Result of levenberg_marquardt().
struct LmResult {
  LmStatus status = LmStatus::MaxIterations;
  /// External values of the last accepted state. A variable whose internal value did not
  /// change keeps its input value bitwise (ADR 0030, point 12).
  std::vector<double> p;
  std::vector<std::uint8_t> changed;   ///< 1 if theta_j differs from the start
  std::vector<std::uint8_t> at_bound;  ///< 1 if p_j lies at a bound (bounds.hpp, at_bound)
  /// Residuals at p; empty if the run was cancelled before the start evaluation.
  std::vector<double> f;
  double F = 0.0;                      ///< 1/2 f^T f; NaN if f is empty
  int iterations = 0;                  ///< k, the number of solves
  std::size_t evaluations = 0;         ///< evaluations of f
  std::size_t failed_evaluations = 0;  ///< evaluations that returned false or threw
  /// Status Failed: the smallest variable index whose Jacobian column had no valid evaluation.
  std::optional<std::size_t> failed_variable;
  std::vector<LmIteration> history;  ///< one entry per solve, accepted or not
};

/// Residuals f(p) in R^m for the external values p (n values, in the units of the variables).
/// Returns false if the evaluation is invalid; an exception counts as false, and so do
/// non-finite residuals. Called concurrently from several threads while the Jacobian is formed,
/// so it must be thread-safe, and it must be deterministic (same p, same f).
using Residuals = std::function<bool(std::span<const double> p, std::span<double> f)>;

/// Minimizes 1/2 f(p)^T f(p) from p0 within the bounds (ADR 0030).
///
/// @param residuals the residual function, see Residuals
/// @param m         number of residuals (>= 1), fixed for the run
/// @param p0        start values (n >= 1), finite and within the bounds
/// @param bounds    one entry per variable
/// @param options   tolerances and steps, see LmOptions
/// @param control   cancellation and progress (stages "jacobian" with 2n evaluations per
///                  Jacobian and "optimize" with k_max solves). A cancellation ends the run
///                  with status Cancelled and the last accepted state, without an exception; an
///                  exception of the progress callback is rethrown.
/// @throws std::invalid_argument for invalid input or if f is invalid at p0
[[nodiscard]] LmResult levenberg_marquardt(const Residuals& residuals,
                                           std::size_t m,
                                           std::vector<double> p0,
                                           std::vector<Bounds> bounds,
                                           const LmOptions& options = {},
                                           const trace::RunControl& control = {});

}  // namespace rtt::optim
