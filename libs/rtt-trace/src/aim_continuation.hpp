#pragma once

/// @file aim_continuation.hpp
/// Fallback of real aiming by continuation (#287, ADR 0007 addendum), separate from sources.cpp
/// so that the tests can drive it with a synthetic Newton run (tests/test_aim_continuation.cpp).
/// Not part of the API.

#include <array>
#include <cmath>
#include <limits>
#include <optional>

namespace rtt::trace::detail {

/// Half-width of the bundle of one field in units of the EP radius, around the chief ray
/// through the EP centre: reserve for the pupil aberration found by real aiming.
inline constexpr double kBundleReserve = 3.0;

/// Fallback of real aiming by continuation (#287, ADR 0007 addendum): first step of the
/// parameter s in [0, 1], smallest step before the fallback gives up, and largest number of
/// accepted steps.
inline constexpr double kAimContinuationFirstStep = 1.0 / 16.0;
inline constexpr double kAimContinuationMinStep = 1.0 / 1024.0;
inline constexpr int kMaxAimContinuationSteps = 64;

/// One Newton run of real aiming (aim()).
struct NewtonRun {
  int steps = 0;  ///< Newton steps taken
  /// Distance to the target on the stop surface at the last iterate, mm; +infinity if the stop
  /// was never reached.
  double residual = std::numeric_limits<double>::infinity();
  /// Last Jacobian computed, {j11, j12, j21, j22} (none if the run took no step).
  std::optional<std::array<double, 4>> jacobian;
};

/// Pupil point (a, b) found by continuation(), with its residual and all Newton steps taken.
struct ContinuationResult {
  double a = 0.0;
  double b = 0.0;
  double residual = 0.0;
  int steps = 0;
};

/// Fallback of real aiming by continuation (#287, ADR 0007 addendum), called by aim() only after
/// the Newton run from the paraxial solution failed.
///
/// Embedding method with the explicit parameter s in [0, 1] (E. L. Allgower, K. Georg,
/// "Continuation and path following", Acta Numerica 1993, Sec. 2, Eq. (2.10) and Algorithm 2.6,
/// with the predictor tangent to the curve of p. 6; docs/quellen.md): the residual
/// H(a, b, s) = hit(a, b) - s (tx, ty) = 0 is followed from the real chief ray (s = 0, target
/// the stop centre) to the target (tx, ty) (s = 1). Since dH/d(a, b) = J and dH/ds = -(tx, ty),
/// the tangent is d(a, b)/ds = J^-1 (tx, ty); the Euler predictor (their Eq. (3.2)) steps along
/// it with the Jacobian of the last corrector run, an approximation of H' as in their
/// Algorithm 3.1. Their t(H'(u)) is the normalized tangent in R^(N+1) and h an arc length; here
/// the same direction is parametrized by s. Corrector: the Newton run of aim() (their Eq. (3.1),
/// here with a square J). Steplength: their Sec. 3.3 asks for an automatic rule; this one is ours:
/// start at 1/16, double after a success, halve the step just tried after a failure (so a step cut
/// to the rest 1 - s is not tried again), give up below 1/1024 or after 64 accepted steps. For an
/// object at infinity every accepted point must stay inside the bundle reserve kBundleReserve r_ep,
/// for which make_field() places the start plane of the field.
/// @param newton      Newton run (u, v, gx, gy, NewtonRun&) -> bool of aim(): from (u, v) on the
///                    target (gx, gy), (u, v) hold the last iterate; true on convergence
/// @param jacobian_at central-difference Jacobian (u, v, gx, gy) of aim()
/// @param infinite    object at infinity (bundle reserve applies)
/// @param r_ep        paraxial EP radius, mm
/// @param tx, ty      target on the stop surface, local stop coordinates, mm
/// @return the pupil point at s = 1, or nullopt if the chief ray or the continuation fails
template <class Newton, class JacobianAt>
std::optional<ContinuationResult> continuation(const Newton& newton,
                                               const JacobianAt& jacobian_at,
                                               bool infinite,
                                               double r_ep,
                                               double tx,
                                               double ty) {
  const auto inside_reserve = [&](double u, double v) {
    return !infinite || std::hypot(u, v) <= kBundleReserve * r_ep;
  };
  ContinuationResult out;
  NewtonRun chief;
  const bool chief_converged = newton(out.a, out.b, 0.0, 0.0, chief);
  out.steps = chief.steps;
  if (!chief_converged || !inside_reserve(out.a, out.b)) return std::nullopt;
  out.residual = chief.residual;
  std::optional<std::array<double, 4>> jacobian =
      chief.jacobian ? chief.jacobian : jacobian_at(out.a, out.b, 0.0, 0.0);
  double s = 0.0;
  double ds = kAimContinuationFirstStep;
  int accepted = 0;
  while (s < 1.0) {
    if (!jacobian || accepted == kMaxAimContinuationSteps) return std::nullopt;
    const auto [j11, j12, j21, j22] = *jacobian;
    const double det = j11 * j22 - j12 * j21;
    if (!(det != 0.0) || !std::isfinite(det)) return std::nullopt;
    const double s_next = ds >= 1.0 - s ? 1.0 : s + ds;  // exactly 1 at the end
    const double dtx = (s_next - s) * tx;
    const double dty = (s_next - s) * ty;
    double u = out.a + (j22 * dtx - j12 * dty) / det;
    double v = out.b + (-j21 * dtx + j11 * dty) / det;
    NewtonRun step;
    const bool converged = newton(u, v, s_next * tx, s_next * ty, step);
    out.steps += step.steps;
    if (converged && inside_reserve(u, v)) {
      out.a = u;
      out.b = v;
      out.residual = step.residual;
      s = s_next;
      ++accepted;
      ds *= 2.0;
      // No predictor follows the step to s = 1, so no Jacobian is needed there.
      if (s < 1.0) jacobian = step.jacobian ? step.jacobian : jacobian_at(u, v, s * tx, s * ty);
    } else {
      ds = 0.5 * (s_next - s);  // the step just tried, possibly cut to 1 - s
      if (ds < kAimContinuationMinStep) return std::nullopt;
    }
  }
  return out;
}

}  // namespace rtt::trace::detail
