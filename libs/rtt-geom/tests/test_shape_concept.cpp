#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <optional>
#include <type_traits>
#include <utility>

#include "rtt/geom/asphere.hpp"
#include "rtt/geom/conic.hpp"
#include "rtt/geom/intersect.hpp"
#include "rtt/geom/plane.hpp"
#include "rtt/geom/shape.hpp"

// ADR 0041, point 4: the Newton intersection is a template over the concrete shape, constrained
// by IntersectableShape; there is no virtual base class, and sag, grad and the intersection are
// noexcept.

using rtt::geom::Conic;
using rtt::geom::EvenAsphere;
using rtt::geom::HitStatus;
using rtt::geom::IntersectableShape;
using rtt::geom::Plane;
using Vec3 = rtt::math::Vec3T<double>;

namespace {
/// Satisfies the concept: z = x^2 / 4 on r <= 3.
struct Bowl {
  [[nodiscard]] double sag(double x, double /*y*/) const noexcept { return 0.25 * x * x; }
  [[nodiscard]] std::pair<double, double> grad(double x, double /*y*/) const noexcept {
    return {0.5 * x, 0.0};
  }
  [[nodiscard]] std::pair<double, double> base_conic() const { return {0.0, 0.0}; }
  [[nodiscard]] std::optional<double> max_radius() const { return 3.0; }
};

/// Same as Bowl, but sag may throw.
struct ThrowingSag {
  [[nodiscard]] double sag(double x, double /*y*/) const { return 0.25 * x * x; }
  [[nodiscard]] std::pair<double, double> grad(double x, double /*y*/) const noexcept {
    return {0.5 * x, 0.0};
  }
  [[nodiscard]] std::pair<double, double> base_conic() const { return {0.0, 0.0}; }
  [[nodiscard]] std::optional<double> max_radius() const { return std::nullopt; }
};

/// Same as Bowl, but without max_radius.
struct NoDomain {
  [[nodiscard]] double sag(double x, double /*y*/) const noexcept { return 0.25 * x * x; }
  [[nodiscard]] std::pair<double, double> grad(double x, double /*y*/) const noexcept {
    return {0.5 * x, 0.0};
  }
  [[nodiscard]] std::pair<double, double> base_conic() const { return {0.0, 0.0}; }
};
}  // namespace

static_assert(IntersectableShape<Plane<double>, double>);
static_assert(IntersectableShape<Conic<double>, double>);
static_assert(IntersectableShape<EvenAsphere<double>, double>);
static_assert(IntersectableShape<Bowl, double>);
static_assert(!IntersectableShape<ThrowingSag, double>);
static_assert(!IntersectableShape<NoDomain, double>);

// No virtual base class: the shapes are plain final classes.
static_assert(!std::is_polymorphic_v<Plane<double>>);
static_assert(!std::is_polymorphic_v<Conic<double>>);
static_assert(!std::is_polymorphic_v<EvenAsphere<double>>);

// The intersection itself is noexcept (declval: only the call, not building its arguments).
static_assert(noexcept(rtt::geom::intersect(std::declval<const EvenAsphere<double>&>(),
                                            std::declval<const Vec3&>(),
                                            std::declval<const Vec3&>())));
static_assert(noexcept(rtt::geom::intersect_newton(std::declval<const Bowl&>(),
                                                   std::declval<const Vec3&>(),
                                                   std::declval<const Vec3&>())));

TEST_CASE("intersect forwards a shape without its own overload to intersect_newton", "[shape]") {
  const EvenAsphere<double> asphere(0.02, -0.5, {1e-5, -2e-8});
  int with_steps = 0;
  for (const auto& [o, d] : {
           std::pair{Vec3(0.0, 8.0, -10.0), Vec3(0.0, 0.0, 1.0)},
           std::pair{Vec3(3.0, -4.0, -20.0), Vec3(0.05, 0.1, 1.0).normalized()},
           std::pair{Vec3(-6.0, 2.0, -5.0), Vec3(0.1, -0.02, 1.0).normalized()},
       }) {
    const auto hit = rtt::geom::intersect(asphere, o, d);
    const auto newton = rtt::geom::intersect_newton(asphere, o, d);
    // Content guard: real hits off the base conic, so Newton has to iterate.
    REQUIRE(hit.status == HitStatus::Hit);
    with_steps += hit.iterations > 0 ? 1 : 0;
    // Bitwise identical: the same function.
    REQUIRE(hit.t == newton.t);
    REQUIRE(hit.point == newton.point);
    REQUIRE(hit.normal == newton.normal);
    REQUIRE(hit.iterations == newton.iterations);
  }
  REQUIRE(with_steps > 0);
}

TEST_CASE("intersect keeps the analytic overloads for plane and conic", "[shape]") {
  // The analytic paths report 0 Newton steps; a constrained template must not take them over.
  // Far ray origin: |F| at the analytic start is limited by ulp(z) ~ 1.5e-11 mm > 1e-12 mm, so
  // Newton takes at least one step there (precondition), while the analytic overload takes none.
  const Conic<double> sphere(0.02, 0.0);
  const Vec3 o(0.0, 3.0, -1e5);
  const Vec3 d = Vec3::UnitZ();
  const auto newton = rtt::geom::intersect_newton(sphere, o, d);
  REQUIRE(newton.status == HitStatus::Hit);
  REQUIRE(newton.iterations > 0);
  const auto conic = rtt::geom::intersect(sphere, o, d);
  REQUIRE(conic.status == HitStatus::Hit);
  REQUIRE(conic.iterations == 0);
  const auto plane = rtt::geom::intersect(Plane<double>(), Vec3(0.0, 1.0, -5.0), d);
  REQUIRE(plane.status == HitStatus::Hit);
  REQUIRE(plane.iterations == 0);
}

TEST_CASE("intersect_newton runs on a shape that only satisfies the concept", "[shape]") {
  // z = x^2 / 4 and the axis-parallel ray at x = 2: hit at z = 1, t = 6, normal ~ (-1, 0, 1).
  const auto hit = rtt::geom::intersect_newton(Bowl{}, Vec3(2.0, 0.0, -5.0), Vec3::UnitZ().eval());
  REQUIRE(hit.status == HitStatus::Hit);
  REQUIRE(std::abs(hit.t - 6.0) <= 1e-12);
  REQUIRE((hit.normal - Vec3(-1.0, 0.0, 1.0).normalized()).cwiseAbs().maxCoeff() <= 1e-12);
  // Outside r <= 3 there is no surface.
  const auto outside =
      rtt::geom::intersect_newton(Bowl{}, Vec3(4.0, 0.0, -5.0), Vec3::UnitZ().eval());
  REQUIRE(outside.status == HitStatus::Missed);
}
