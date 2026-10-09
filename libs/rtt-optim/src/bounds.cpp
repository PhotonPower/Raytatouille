#include "rtt/optim/bounds.hpp"

#include <algorithm>
#include <cmath>

namespace rtt::optim {

namespace {

// Relative threshold of at_bound() (ADR 0030, point 5).
constexpr double kAtBound = 1e-6;

// MINUIT User's Guide eqs. (1.3)/(1.5) for the distance x >= 0 from the bound, written without
// cancellation: (x + 1)^2 - 1 = x (x + 2). Same equation; the direct form loses the relative
// accuracy of small x (natural size 1e-4, ADR 0030 point 14).
double one_sided_theta(double x) {
  const double d = std::max(x, 0.0);  // rounding of p - a may give -0
  return std::sqrt(d * (d + 2.0));
}

// MINUIT User's Guide eqs. (1.4)/(1.6): distance from the bound sqrt(theta^2 + 1) - 1, written
// without cancellation as theta^2 / (sqrt(theta^2 + 1) + 1). Same equation; >= 0.
double one_sided_distance(double theta) {
  const double t2 = theta * theta;
  return t2 / (std::sqrt(t2 + 1.0) + 1.0);
}

}  // namespace

double to_internal(double p, const Bounds& bounds) {
  if (bounds.min && bounds.max) {
    // MINUIT eq. (1.1); the argument is clamped against rounding at the bounds.
    const double a = *bounds.min;
    const double b = *bounds.max;
    return std::asin(std::clamp(2.0 * (p - a) / (b - a) - 1.0, -1.0, 1.0));
  }
  if (bounds.min) return one_sided_theta(p - *bounds.min);
  if (bounds.max) return one_sided_theta(*bounds.max - p);
  return p;
}

double to_external(double theta, const Bounds& bounds) {
  if (bounds.min && bounds.max) {
    // MINUIT eq. (1.2); a + (b - a) may round past b, so the result is clamped to [a, b].
    const double a = *bounds.min;
    const double b = *bounds.max;
    return std::clamp(a + (b - a) / 2.0 * (std::sin(theta) + 1.0), a, b);
  }
  if (bounds.min) return *bounds.min + one_sided_distance(theta);
  if (bounds.max) return *bounds.max - one_sided_distance(theta);
  return theta;
}

bool at_bound(double p, const Bounds& bounds) {
  if (bounds.min && bounds.max) {
    const double tolerance = kAtBound * (*bounds.max - *bounds.min);
    return p - *bounds.min <= tolerance || *bounds.max - p <= tolerance;
  }
  if (bounds.min) return p - *bounds.min <= kAtBound * std::max(1.0, std::abs(*bounds.min));
  if (bounds.max) return *bounds.max - p <= kAtBound * std::max(1.0, std::abs(*bounds.max));
  return false;
}

}  // namespace rtt::optim
