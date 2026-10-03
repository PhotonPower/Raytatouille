#include <catch2/catch_test_macros.hpp>
#include <cmath>
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
/// (|d| = 1 for a^2 + b^2 = 1) and tau = t - 1 this gives F = tau^3 - 2 tau + 2, the textbook case
/// in which Newton's method cycles 0 -> 1 -> 0 (e.g. Burden & Faires, Numerical Analysis,
/// Sec. 2.3). The plane start z = 0 lies at tau = 0, so Newton never converges.
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
