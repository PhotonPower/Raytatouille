#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <limits>
#include <optional>
#include <utility>
#include <vector>

#include "rtt/geom/asphere.hpp"
#include "rtt/geom/conic.hpp"
#include "rtt/geom/intersect.hpp"
#include "rtt/geom/shape.hpp"

using rtt::geom::Conic;
using rtt::geom::EvenAsphere;
using rtt::geom::HitStatus;
using rtt::geom::Intersection;
using Vec3 = rtt::math::Vec3T<double>;

// General intersection (docs/architecture.md, Physik-Module, Geometrie): start at the analytic
// hit of the base conic, then Newton on F(t) = z(t) - sag(x(t), y(t)) with tolerance 1e-12 mm and
// at most 30 iterations, otherwise NoConvergence. Newton's method: e.g. Press et al.,
// Numerical Recipes, 3rd ed., Sec. 9.4.

namespace {
bool near(double a, double b, double tol) {
  return std::abs(a - b) <= tol;
}

bool near(const Vec3& a, const Vec3& b, double tol) {
  return (a - b).cwiseAbs().maxCoeff() <= tol;
}

bool all_finite(const Intersection<double>& hit) {
  return std::isfinite(hit.t) && hit.point.allFinite() && hit.normal.allFinite();
}

double residual(const rtt::geom::Shape<double>& shape, const Vec3& p) {
  return p.z() - shape.sag(p.x(), p.y());
}

/// Test shape z = (b + 2) u - u^3 - 2 with u = x / a. For the ray p(t) = (a (t - 1), 0, b (t - 1))
/// (|d| = 1 for a^2 + b^2 = 1) and tau = t - 1 this gives F = tau^3 - 2 tau + 2 with
/// F' = 3 tau^2 - 2. Newton cycles: tau = 0 -> 0 - 2 / (-2) = 1 -> 1 - 1 / 1 = 0. The cycle is
/// superattracting (F''(0) = 0), so rounding cannot break it. The plane start z = 0 lies at
/// tau = 0, so Newton never converges.
class NewtonCycleShape final : public rtt::geom::Shape<double> {
 public:
  static constexpr double kA = 0.6;
  static constexpr double kB = 0.8;

  [[nodiscard]] double sag(double x, double /*y*/) const override {
    const double u = x / kA;
    return (kB + 2.0) * u - u * u * u - 2.0;
  }
  [[nodiscard]] std::pair<double, double> grad(double x, double /*y*/) const override {
    const double u = x / kA;
    return {((kB + 2.0) - 3.0 * u * u) / kA, 0.0};
  }
  [[nodiscard]] std::pair<double, double> base_conic() const override { return {0.0, 0.0}; }
  [[nodiscard]] std::optional<double> max_radius() const override { return std::nullopt; }
};
/// Tilted plane z = slope * x with an optional domain limit, for the failure branches of Newton.
/// F is linear in t, so Newton lands on the exact crossing in one step.
class TiltedShape final : public rtt::geom::Shape<double> {
 public:
  TiltedShape(double slope, std::optional<double> max_radius)
      : slope_(slope), max_radius_(max_radius) {}

  [[nodiscard]] double sag(double x, double /*y*/) const override { return slope_ * x; }
  [[nodiscard]] std::pair<double, double> grad(double /*x*/, double /*y*/) const override {
    return {slope_, 0.0};
  }
  [[nodiscard]] std::pair<double, double> base_conic() const override { return {0.0, 0.0}; }
  [[nodiscard]] std::optional<double> max_radius() const override { return max_radius_; }

 private:
  double slope_;
  std::optional<double> max_radius_;
};

/// Independent reference: first sign change of F(t) = z(t) - sag(x(t), y(t)) on a fine scan
/// of (t_lo, t_hi], refined by bisection.
double first_crossing(
    const rtt::geom::Shape<double>& shape, const Vec3& o, const Vec3& d, double t_lo, double t_hi) {
  const auto f = [&](double t) { return residual(shape, o + t * d); };
  const int steps = 200000;
  double prev = t_lo;
  for (int i = 1; i <= steps; ++i) {
    const double next = t_lo + (t_hi - t_lo) * i / steps;
    if ((f(prev) < 0.0) != (f(next) < 0.0)) {
      double lo = prev;
      double hi = next;
      for (int j = 0; j < 200 && hi - lo > 1e-15; ++j) {
        const double mid = 0.5 * (lo + hi);
        if ((f(lo) < 0.0) == (f(mid) < 0.0)) {
          lo = mid;
        } else {
          hi = mid;
        }
      }
      return 0.5 * (lo + hi);
    }
    prev = next;
  }
  return std::numeric_limits<double>::quiet_NaN();
}
}  // namespace

TEST_CASE("asphere without coefficients intersects like the conic", "[newton]") {
  // Issue #4: identical intersection, absolute 1e-15 mm.
  for (const auto& [c, k] : {std::pair{0.02, 0.0}, std::pair{-0.02, -1.0}, std::pair{0.01, -2.0}}) {
    const Conic<double> conic(c, k);
    const EvenAsphere<double> asphere(c, k, {});
    for (const auto& [o, d] : {
             std::pair{Vec3(0.0, 5.0, -10.0), Vec3(0.0, 0.0, 1.0)},
             std::pair{Vec3(1.0, -2.0, -20.0), Vec3(0.05, 0.1, 1.0).normalized()},
             std::pair{Vec3(-3.0, 4.0, 15.0), Vec3(0.02, -0.03, -1.0).normalized()},
         }) {
      const auto expected = rtt::geom::intersect(conic, o, d);
      const auto hit = rtt::geom::intersect(asphere, o, d);
      REQUIRE(expected.status == HitStatus::Hit);
      REQUIRE(hit.status == HitStatus::Hit);
      REQUIRE(near(hit.t, expected.t, 1e-15));
      REQUIRE(near(hit.point, expected.point, 1e-15));
      REQUIRE(near(hit.normal, expected.normal, 1e-15));
    }
  }
}

TEST_CASE("asphere hits lie on the surface for a ray grid up to 0.9 aperture", "[newton]") {
  // Strong asphere on a prolate base; aperture radius 20 mm.
  const double aperture = 20.0;
  for (const auto& asphere : {
           EvenAsphere<double>(1.0 / 40.0, -0.5, {-2e-6, 1e-9, -3e-12}),
           EvenAsphere<double>(-1.0 / 30.0, 0.0, {5e-5, -1e-7}),
           EvenAsphere<double>(0.0, 0.0, {1e-4}),
       }) {
    int converged_with_steps = 0;
    for (const Vec3& d : {Vec3(0.0, 0.0, 1.0), Vec3(0.05, -0.1, 1.0).normalized(),
                          Vec3(-0.15, 0.0, 1.0).normalized()}) {
      for (int i = -4; i <= 4; ++i) {
        for (int j = -4; j <= 4; ++j) {
          const double x = 0.9 * aperture * i / 4.0;
          const double y = 0.9 * aperture * j / 4.0;
          if (x * x + y * y > 0.81 * aperture * aperture) {
            continue;
          }
          // Start 30 mm before the vertex plane, aimed at (x, y, 0).
          const Vec3 o = Vec3(x, y, 0.0) - 30.0 / d.z() * d;
          const auto hit = rtt::geom::intersect(asphere, o, d);
          REQUIRE(hit.status == HitStatus::Hit);
          REQUIRE(std::abs(residual(asphere, hit.point)) < 1e-12);
          REQUIRE(near(hit.point, o + hit.t * d, 1e-12));
          REQUIRE(hit.iterations <= rtt::geom::kMaxNewtonIterations);
          // Normal from the gradient (docs/architecture.md).
          const auto [gx, gy] = asphere.grad(hit.point.x(), hit.point.y());
          REQUIRE(near(hit.normal, Vec3(-gx, -gy, 1.0).normalized(), 1e-12));
          if (hit.iterations > 0) {
            ++converged_with_steps;
          }
        }
      }
    }
    // The polynomial terms move the hit away from the base conic, so Newton has to work.
    REQUIRE(converged_with_steps > 0);
  }
}

TEST_CASE("analytic intersections report zero Newton iterations", "[newton]") {
  const auto hit =
      rtt::geom::intersect(Conic<double>(0.02, 0.0), Vec3(0.0, 1.0, -5.0), Vec3(0.0, 0.0, 1.0));
  REQUIRE(hit.status == HitStatus::Hit);
  REQUIRE(hit.iterations == 0);
}

TEST_CASE("Newton that does not converge reports NoConvergence", "[newton]") {
  const NewtonCycleShape cycle;
  const rtt::geom::Shape<double>& shape = cycle;
  const double a = NewtonCycleShape::kA;
  const double b = NewtonCycleShape::kB;
  const Vec3 o(-a, 0.0, -b);
  const Vec3 d(a, 0.0, b);
  const auto hit = rtt::geom::intersect(shape, o, d);
  REQUIRE(hit.status == HitStatus::NoConvergence);
  REQUIRE(hit.iterations == rtt::geom::kMaxNewtonIterations);
  REQUIRE(all_finite(hit));
}

TEST_CASE("asphere: ray beyond the domain gives Missed without NaN", "[newton]") {
  // Spherical base R = 50 bounds the asphere at r = 50 mm.
  const EvenAsphere<double> asphere(0.02, 0.0, {1e-7});
  const auto hit = rtt::geom::intersect(asphere, Vec3(0.0, 60.0, -10.0), Vec3(0.0, 0.0, 1.0));
  REQUIRE(hit.status == HitStatus::Missed);
  REQUIRE(all_finite(hit));
}

TEST_CASE("asphere: surface behind the ray gives Missed", "[newton]") {
  const EvenAsphere<double> asphere(0.02, -1.0, {1e-6});
  const auto hit = rtt::geom::intersect(asphere, Vec3(0.0, 3.0, -10.0), Vec3(0.0, 0.0, -1.0));
  REQUIRE(hit.status == HitStatus::Missed);
  REQUIRE(all_finite(hit));
}

TEST_CASE("asphere on a flat base: start from the plane", "[newton]") {
  // c = 0: the base conic is the plane z = 0. z = A4 r^4 with A4 = 1e-4 /mm^3; a ray parallel to
  // the axis at r = 10 hits at z = 1 mm.
  const EvenAsphere<double> asphere(0.0, 0.0, {1e-4});
  const auto hit = rtt::geom::intersect(asphere, Vec3(10.0, 0.0, -5.0), Vec3(0.0, 0.0, 1.0));
  REQUIRE(hit.status == HitStatus::Hit);
  REQUIRE(near(hit.point.z(), 1.0, 1e-12));
  REQUIRE(near(hit.t, 6.0, 1e-12));
}

TEST_CASE("asphere: Newton finds the first crossing along the ray", "[newton]") {
  // Strong quartic bowl z = 0.01 r^4 and an oblique ray; the reference is a scan plus bisection
  // that does not use the gradient.
  const EvenAsphere<double> asphere(0.0, 0.0, {0.01});
  for (const auto& [o, d] : {
           std::pair{Vec3(-2.0, 0.0, -1.0), Vec3(0.3, 0.0, 1.0).normalized()},
           std::pair{Vec3(1.0, 1.0, -2.0), Vec3(-0.2, 0.1, 1.0).normalized()},
           std::pair{Vec3(0.5, -3.0, 2.0), Vec3(0.0, 0.4, -1.0).normalized()},
       }) {
    const auto hit = rtt::geom::intersect(asphere, o, d);
    REQUIRE(hit.status == HitStatus::Hit);
    REQUIRE(near(hit.t, first_crossing(asphere, o, d, rtt::geom::kDefaultTMin, 20.0), 1e-9));
  }
}

TEST_CASE("Newton leaving the domain of the shape reports NoConvergence", "[newton]") {
  // z = x / 2 limited to r <= 1. Ray o = (0, 0, -1), d = (0.6, 0, 0.8): the plane start is at
  // t = 1.25 (x = 0.75, inside); the crossing -1 + 0.8 t = 0.3 t lies at t = 2, x = 1.2 (outside).
  const TiltedShape tilted(0.5, 1.0);
  const auto hit = rtt::geom::intersect<double>(tilted, Vec3(0.0, 0.0, -1.0), Vec3(0.6, 0.0, 0.8));
  REQUIRE(hit.status == HitStatus::NoConvergence);
  REQUIRE(hit.iterations == 1);
  REQUIRE(all_finite(hit));
}

TEST_CASE("Newton with F' = 0 (ray parallel to the surface) reports NoConvergence", "[newton]") {
  // z = x / 2 and d ~ (2, 0, 1): F'(t) = d_z - d_x / 2 = 0, F = -0.1 mm everywhere.
  const TiltedShape tilted(0.5, 1.0);
  const auto hit =
      rtt::geom::intersect<double>(tilted, Vec3(0.0, 0.0, -0.1), Vec3(2.0, 0.0, 1.0).normalized());
  REQUIRE(hit.status == HitStatus::NoConvergence);
  REQUIRE(hit.iterations == 0);
  REQUIRE(all_finite(hit));
}

TEST_CASE("Newton converging behind the ray reports NoConvergence", "[newton]") {
  // z = 2 x, o = (0, 0, -1), d = (0.6, 0, 0.8): plane start at t = 1.25, but the only crossing
  // -1 + 0.8 t = 1.2 t is at t = -2.5. Newton cannot rule out a crossing ahead, so the result
  // is NoConvergence rather than Missed.
  const TiltedShape tilted(2.0, std::nullopt);
  const auto hit = rtt::geom::intersect<double>(tilted, Vec3(0.0, 0.0, -1.0), Vec3(0.6, 0.0, 0.8));
  REQUIRE(hit.status == HitStatus::NoConvergence);
  REQUIRE(hit.iterations == 1);
  REQUIRE(all_finite(hit));
}

TEST_CASE("Newton converges at the rounding limit for far ray origins", "[newton]") {
  // |F| < 1e-12 mm is unreachable when o + t d carries rounding errors of that size. The second
  // criterion |dt| <= 8 eps (1 + |t| + |o|) accepts the hit at the rounding limit.
  const Conic<double> sphere(0.02, 0.0);
  const rtt::geom::Shape<double>& shape = sphere;
  // Ray parallel to the axis: t = 1e5 + sag(r) exactly, rounding of the reference ~ ulp(1e5).
  const auto axial = rtt::geom::intersect(shape, Vec3(0.0, 3.0, -1e5), Vec3(0.0, 0.0, 1.0));
  REQUIRE(axial.status == HitStatus::Hit);
  REQUIRE(near(axial.t, 1e5 + sphere.sag(0.0, 3.0), 1e-10));
  REQUIRE(axial.iterations <= rtt::geom::kMaxNewtonIterations);
  // Oblique ray against the analytic intersection. Its quadratic coefficient g ~ c o_z^2 ~ 2e8 is
  // itself rounded, so the comparison is relative (1e-12, i.e. 1e-7 mm at t = 1e5 mm).
  const Vec3 o(-2.0, 4.0, -1e5);
  const Vec3 d = Vec3(1e-5, -2e-5, 1.0).normalized();
  const auto expected = rtt::geom::intersect(sphere, o, d);
  const auto hit = rtt::geom::intersect(shape, o, d);
  REQUIRE(expected.status == HitStatus::Hit);
  REQUIRE(hit.status == HitStatus::Hit);
  REQUIRE(std::abs(hit.t - expected.t) <= 1e-12 * expected.t);
}
