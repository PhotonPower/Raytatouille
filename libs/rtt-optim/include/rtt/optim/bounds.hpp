#pragma once

/// @file bounds.hpp
/// Bounds of optimization variables as variable transformations (ADR 0030, point 5).
///
/// The solver works with unbounded internal values theta; the model sees external values p in
/// the units of the parameter (mm, 1/mm, ...). Source: F. James, M. Winkler, "MINUIT User's
/// Guide", CERN 2004, section 1.3.1, eqs. (1.1)-(1.6) (docs/quellen.md):
/// - two bounds a < b:  p = a + (b - a)/2 (sin theta + 1),  theta = arcsin(2 (p - a)/(b - a) - 1)
///   (eqs. (1.2), (1.1));
/// - lower bound a:     p = a - 1 + sqrt(theta^2 + 1),     theta = sqrt((p - a + 1)^2 - 1)
///   (eqs. (1.4), (1.3));
/// - upper bound b:     p = b + 1 - sqrt(theta^2 + 1),     theta = sqrt((b - p + 1)^2 - 1)
///   (eqs. (1.6), (1.5));
/// - no bound:          p = theta.
/// The source leaves the sign of theta open (+-) for one bound; to_internal() returns theta >= 0,
/// and theta and -theta give the same p. The scale of the one-sided transformation is the 1 of
/// the source, in the units of the parameter (ADR 0030, point 5). At a bound dp/dtheta = 0.

#include <optional>

namespace rtt::optim {

/// Bounds of one variable, in the units of the variable. Either side may be missing; if both
/// are set, min < max (checked by levenberg_marquardt()).
struct Bounds {
  std::optional<double> min;
  std::optional<double> max;
};

/// Internal value theta of the external value p (p within the bounds). For one bound the
/// branch theta >= 0; for two bounds theta in [-pi/2, pi/2].
[[nodiscard]] double to_internal(double p, const Bounds& bounds);

/// External value p of the internal value theta; always within the bounds (clamped against
/// rounding at the bound).
[[nodiscard]] double to_external(double theta, const Bounds& bounds);

/// True if p lies at a bound for the warning optim.parameter_at_bound (ADR 0030, point 5):
/// within 1e-6 (b - a) of a bound for two bounds, within 1e-6 max(1, |a|) of the bound for one.
/// Always false without bounds.
[[nodiscard]] bool at_bound(double p, const Bounds& bounds);

}  // namespace rtt::optim
