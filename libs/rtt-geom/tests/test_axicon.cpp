#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <random>
#include <utility>
#include <vector>

#include "rtt/geom/axicon.hpp"
#include "rtt/geom/intersect.hpp"
#include "rtt/geom/odd_asphere.hpp"
#include "rtt/geom/plane.hpp"
#include "rtt/geom/shape.hpp"

// Axicon of ADR 0043, point 1 (derived there): z = r tan(alpha), and its analytic intersection of
// point 2 (cone as the quadric s^2 (x^2 + y^2) - z^2 = 0 on the sheet z s >= 0, discriminant
// without cancellation, shift into the vertex plane for directions inside the cone). Acceptance
// 9 and 10 of ADR 0043.

using rtt::geom::Axicon;
using rtt::geom::DomainShape;
using rtt::geom::HitStatus;
using rtt::geom::IntersectableShape;
using rtt::geom::Intersection;
using rtt::geom::OddAsphere;
using rtt::geom::ParaxialShape;
using rtt::geom::Plane;
using Vec3 = rtt::math::Vec3T<double>;

static_assert(IntersectableShape<Axicon<double>, double>);
static_assert(ParaxialShape<Axicon<double>, double>);
static_assert(!DomainShape<Axicon<double>, double>);

namespace {

const double kPi = std::acos(-1.0);

double rad(double deg) {
  return deg * kPi / 180.0;
}

/// F(t) = z(t) - s r(t) along the ray, written out independently of the library.
double cone_f(const Vec3& o, const Vec3& d, double s, double t) {
  return (o.z() + t * d.z()) - s * std::hypot(o.x() + t * d.x(), o.y() + t * d.y());
}

/// Root of cone_f in [lo, hi] (sign change required) by bisection down to adjacent doubles.
double bisect(const Vec3& o, const Vec3& d, double s, double lo, double hi) {
  double flo = cone_f(o, d, s, lo);
  for (int i = 0; i < 200; ++i) {
    const double mid = lo + (hi - lo) / 2.0;
    if (mid == lo || mid == hi) break;
    const double fm = cone_f(o, d, s, mid);
    if (fm == 0.0) return mid;
    if ((fm < 0.0) == (flo < 0.0)) {
      lo = mid;
      flo = fm;
    } else {
      hi = mid;
    }
  }
  return lo + (hi - lo) / 2.0;
}

/// Reference for flat rays: the first sign change of cone_f on a grid of t in [1e-9, 200] mm with
/// spacing 1e-3 mm, refined by bisection; nullopt without a sign change.
std::optional<double> first_crossing(const Vec3& o, const Vec3& d, double s) {
  constexpr int kSteps = 200000;
  const double t_lo = 1e-9;
  const double t_hi = 200.0;
  double previous_t = t_lo;
  double previous_f = cone_f(o, d, s, t_lo);
  for (int i = 1; i <= kSteps; ++i) {
    const double t = t_lo + (t_hi - t_lo) * i / kSteps;
    const double f = cone_f(o, d, s, t);
    if (f == 0.0) return t;
    if ((f < 0.0) != (previous_f < 0.0)) return bisect(o, d, s, previous_t, t);
    previous_t = t;
    previous_f = f;
  }
  return std::nullopt;
}

}  // namespace

TEST_CASE("axicon sag, slope and gradient", "[geom][axicon]") {
  const Axicon<double> axicon(rad(-5.0));
  CHECK(axicon.angle() == rad(-5.0));
  CHECK(axicon.slope() == std::tan(axicon.angle()));
  const double s = axicon.slope();
  CHECK(axicon.sag(3.0, -4.0) == s * 5.0);
  CHECK(axicon.sag(0.0, 0.0) == 0.0);
  // F8: the gradient at the apex is (0, 0), so an axial ray has the normal +z.
  CHECK(axicon.grad(0.0, 0.0) == std::pair{0.0, 0.0});
  for (const auto& [x, y] : {std::pair{3.0, -4.0}, std::pair{-0.25, 0.5}, std::pair{7.0, 0.0}}) {
    const double h = 1e-6;
    const auto [gx, gy] = axicon.grad(x, y);
    CHECK(std::abs(gx - (axicon.sag(x + h, y) - axicon.sag(x - h, y)) / (2.0 * h)) <= 1e-9);
    CHECK(std::abs(gy - (axicon.sag(x, y + h) - axicon.sag(x, y - h)) / (2.0 * h)) <= 1e-9);
    // The slope along r is tan(alpha) everywhere off the apex.
    CHECK(std::abs(std::hypot(gx, gy) - std::abs(s)) <= 4e-16);
  }
  // The domain is unbounded; the start value (0, 0) only serves the concept (ADR 0043).
  CHECK_FALSE(axicon.max_radius().has_value());
  CHECK(axicon.base_conic() == std::pair{0.0, 0.0});
}

TEST_CASE("axicon vertex expansion: order 0 for alpha != 0, the plane for alpha = 0",
          "[geom][axicon][vertex]") {
  // ADR 0043, point 3: no paraxial model at the apex of a cone (F6).
  const auto cone = Axicon<double>(rad(20.0)).vertex_expansion();
  CHECK(cone.order == 0);
  CHECK(cone.rotational4);
  CHECK(cone.value == 0.0);
  const auto flat = Axicon<double>(0.0).vertex_expansion();
  const auto plane = Plane<double>().vertex_expansion();
  CHECK(flat.order == 4);
  CHECK(flat.value == plane.value);
  CHECK(flat.slope == plane.slope);
  CHECK(flat.hessian == plane.hessian);
  CHECK(flat.rotational4);
  CHECK(flat.conic == plane.conic);
  CHECK(flat.a4 == plane.a4);
}

TEST_CASE("axicon intersection for rays passing near the apex", "[geom][axicon][intersect]") {
  // ADR 0043, acceptance 9: 20 000 rays per distance from the apex, tilted up to 0.9 of the
  // angle between axis and generatrix, so that a hit is certain; no miss, t against bisection
  // within 1e-12 mm. The naive discriminant misses about half of them at 1e-9 mm (ADR 0043,
  // point 2). Fixed seed (AGENTS.md, rule 7).
  //
  // Bisection bound: along such a ray F' = cos(theta + |alpha|) / cos(alpha) >= sin(0.1 (90 deg -
  // |alpha|)) >= sin(3 deg) for |alpha| <= 60 deg, and F carries a rounding error of a few eps
  // (|o| + t) <= 50 eps, so the reference root is good to ~2e-13 mm.
  for (const double off : {0.0, 1e-9, 1e-6, 1e-3}) {
    INFO("distance from the apex " << off << " mm");
    std::mt19937_64 rng(11);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    int misses = 0;
    double worst_t = 0.0;
    double worst_z = 0.0;
    for (int i = 0; i < 20000; ++i) {
      const double alpha = rad(-60.0 + 120.0 * unit(rng));
      const Axicon<double> axicon(alpha);
      const double theta = (kPi / 2.0 - std::abs(alpha)) * 0.9 * unit(rng);
      const double phi = 2.0 * kPi * unit(rng);
      Vec3 d(std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta));
      if (axicon.slope() < 0.0) d.z() = -d.z();  // a cone opening to -z: the ray comes from +z
      const Vec3 pass = Vec3(2.0 * unit(rng) - 1.0, 2.0 * unit(rng) - 1.0, 0.0) * off;
      const double length = 5.0 + 45.0 * unit(rng);
      const Vec3 o = pass - length * d;
      const Intersection<double> hit = rtt::geom::intersect(axicon, o, d);
      if (hit.status != HitStatus::Hit) {
        ++misses;
        continue;
      }
      const double sl = axicon.slope();
      const double t_ref = off > 0.0 ? bisect(o, d, sl, length - 1.0, length + 1.0) : length;
      worst_t = std::max(worst_t, std::abs(hit.t - t_ref));
      worst_z =
          std::max(worst_z, std::abs(hit.point.z() - axicon.sag(hit.point.x(), hit.point.y())));
    }
    CHECK(misses == 0);
    CHECK(worst_t <= 1e-12);
    CHECK(worst_z <= 1e-12);
  }
}

TEST_CASE("axicon intersection of flat rays against a grid and bisection",
          "[geom][axicon][intersect]") {
  // ADR 0043, acceptance 9: d_z = 0 and up to 2 deg from the horizontal; directions outside the
  // cone take the unshifted form of the discriminant. Same hits, t within 1e-11 mm.
  std::mt19937_64 rng(55);
  std::uniform_real_distribution<double> unit(0.0, 1.0);
  int disagreements = 0;
  int hits = 0;
  double worst = 0.0;
  for (int i = 0; i < 600; ++i) {
    const double alpha = rad(5.0 + 55.0 * unit(rng)) * (i % 2 == 1 ? 1.0 : -1.0);
    const Axicon<double> axicon(alpha);
    const double s = axicon.slope();
    const double phi = 2.0 * kPi * unit(rng);
    const double elevation = i % 3 == 0 ? 0.0 : rad(-2.0 + 4.0 * unit(rng));
    const Vec3 d(std::cos(elevation) * std::cos(phi), std::cos(elevation) * std::sin(phi),
                 std::sin(elevation));
    const Vec3 aim(-20.0 + 40.0 * unit(rng), -20.0 + 40.0 * unit(rng),
                   s * (0.5 + 14.5 * unit(rng)));
    const Vec3 o = aim - 60.0 * d;
    const Intersection<double> hit = rtt::geom::intersect(axicon, o, d);
    const std::optional<double> reference = first_crossing(o, d, s);
    if ((hit.status == HitStatus::Hit) != reference.has_value()) {
      ++disagreements;
      continue;
    }
    if (reference) {
      ++hits;
      worst = std::max(worst, std::abs(hit.t - *reference));
    }
  }
  CHECK(disagreements == 0);
  CHECK(worst <= 1e-11);
  CHECK(hits > 100);  // content guard: most of these rays hit
}

TEST_CASE("axicon intersection against intersect_newton", "[geom][axicon][intersect]") {
  // ADR 0043, acceptance 9: away from the apex the Newton iteration on the same shape agrees
  // within 1e-12 mm.
  std::mt19937_64 rng(314);
  std::uniform_real_distribution<double> unit(0.0, 1.0);
  int compared = 0;
  double worst = 0.0;
  for (int i = 0; i < 2000; ++i) {
    const Axicon<double> axicon(rad(-60.0 + 120.0 * unit(rng)));
    // Start 2 to 5 mm off the axis, tilt below 0.02: the ray stays more than 1 mm off the apex.
    const double r0 = 2.0 + 3.0 * unit(rng);
    const double phi = 2.0 * kPi * unit(rng);
    const Vec3 o(r0 * std::cos(phi), r0 * std::sin(phi), -12.0);
    const Vec3 d = Vec3(0.04 * unit(rng) - 0.02, 0.04 * unit(rng) - 0.02, 1.0).normalized();
    const Intersection<double> analytic = rtt::geom::intersect(axicon, o, d);
    const Intersection<double> newton = rtt::geom::intersect_newton(axicon, o, d);
    REQUIRE(analytic.status == HitStatus::Hit);
    REQUIRE(newton.status == HitStatus::Hit);
    worst = std::max(worst, std::abs(analytic.t - newton.t));
    CHECK((analytic.normal - newton.normal).cwiseAbs().maxCoeff() <= 1e-12);
    ++compared;
  }
  CHECK(worst <= 1e-12);
  CHECK(compared == 2000);
}

TEST_CASE("axicon with alpha = 0 intersects as the plane bit for bit",
          "[geom][axicon][intersect]") {
  const Axicon<double> axicon(0.0);
  const Plane<double> plane;
  for (const auto& [o, d] : {std::pair{Vec3(0.0, 0.0, -5.0), Vec3(0.0, 0.0, 1.0)},
                             std::pair{Vec3(1.0, -2.0, -3.0), Vec3(0.1, 0.2, 1.0).normalized()},
                             std::pair{Vec3(1.0, -2.0, 3.0), Vec3(0.1, 0.2, -1.0).normalized()},
                             std::pair{Vec3(1.0, -2.0, -3.0), Vec3(1.0, 0.0, 0.0)},
                             std::pair{Vec3(1.0, -2.0, 3.0), Vec3(0.1, 0.2, 1.0).normalized()}}) {
    const Intersection<double> a = rtt::geom::intersect(axicon, o, d);
    const Intersection<double> b = rtt::geom::intersect(plane, o, d);
    CHECK(a.status == b.status);
    CHECK(a.t == b.t);
    CHECK(a.point == b.point);
    CHECK(a.normal == b.normal);
    CHECK(a.iterations == b.iterations);
  }
}

TEST_CASE("a ray parallel to a generatrix hits the axicon once", "[geom][axicon][intersect]") {
  // ADR 0043, point 2: a = 0, the linear solution. alpha = 30 deg, the generatrix x > 0 has the
  // direction (cos alpha, 0, sin alpha). The ray from (-10, 0, -10 s + 1) along it stays above
  // that generatrix (z - s x = 1) and crosses the other half x < 0 where
  // z + s x = 0: t = (20 s - 1) / (2 sin alpha).
  const double alpha = rad(30.0);
  const Axicon<double> axicon(alpha);
  const double s = axicon.slope();
  const Vec3 d(std::cos(alpha), 0.0, std::sin(alpha));
  const Vec3 o(-10.0, 0.0, -10.0 * s + 1.0);
  const Intersection<double> hit = rtt::geom::intersect(axicon, o, d);
  REQUIRE(hit.status == HitStatus::Hit);
  const double t = (20.0 * s - 1.0) / (2.0 * std::sin(alpha));
  CHECK(std::abs(hit.t - t) <= 1e-12);
  REQUIRE(hit.point.x() < 0.0);
  CHECK(std::abs(hit.point.z() - axicon.sag(hit.point.x(), hit.point.y())) <= 1e-12);
  // No second hit further along the ray.
  const Intersection<double> after = rtt::geom::intersect(axicon, hit.point + 1e-6 * d, d);
  CHECK(after.status == HitStatus::Missed);
}

TEST_CASE("axicon normal and the axial ray through the apex", "[geom][axicon][intersect]") {
  const Axicon<double> axicon(rad(-5.0));
  const Intersection<double> axial =
      rtt::geom::intersect(axicon, Vec3(0.0, 0.0, -10.0), Vec3(0.0, 0.0, 1.0));
  REQUIRE(axial.status == HitStatus::Hit);
  CHECK(axial.t == 10.0);
  CHECK(axial.normal == Vec3(0.0, 0.0, 1.0));
  // Off the apex the normal is (-s x / r, -s y / r, 1) normalized.
  const double s = axicon.slope();
  const Intersection<double> hit =
      rtt::geom::intersect(axicon, Vec3(0.0, 2.0, -10.0), Vec3(0.0, 0.0, 1.0));
  REQUIRE(hit.status == HitStatus::Hit);
  CHECK(std::abs(hit.point.z() - 2.0 * s) <= 1e-15);
  const Vec3 n = Vec3(0.0, -s, 1.0).normalized();
  CHECK((hit.normal - n).cwiseAbs().maxCoeff() <= 1e-16);
}

TEST_CASE("axicon equals the odd asphere with c = 0 and A1 = tan(alpha)",
          "[geom][axicon][odd_asphere]") {
  // ADR 0043, acceptance 10: sag bit for bit, intersections within 1e-12 mm (analytic against
  // Newton on the odd asphere).
  for (const double alpha_deg : {-40.0, -5.0, 5.0, 25.0}) {
    INFO("alpha " << alpha_deg << " deg");
    const Axicon<double> axicon(rad(alpha_deg));
    const OddAsphere<double> odd(0.0, 0.0, {axicon.slope()});
    for (int i = 0; i <= 40; ++i) {
      const double x = -10.0 + 0.5 * i;
      const double y = 0.25 * i - 3.0;
      CHECK(axicon.sag(x, y) == odd.sag(x, y));
    }
    double worst = 0.0;
    for (int i = -5; i <= 5; ++i) {
      for (int j = -5; j <= 5; ++j) {
        if (i == 0 && j == 0) continue;
        const Vec3 o(1.1 * i, 0.9 * j, -20.0);
        const Vec3 d = Vec3(0.01 * j, -0.02 * i, 1.0).normalized();
        const Intersection<double> a = rtt::geom::intersect(axicon, o, d);
        const Intersection<double> b = rtt::geom::intersect(odd, o, d);
        REQUIRE(a.status == HitStatus::Hit);
        REQUIRE(b.status == HitStatus::Hit);
        worst = std::max(worst, (a.point - b.point).cwiseAbs().maxCoeff());
      }
    }
    CHECK(worst <= 1e-12);
  }
}
