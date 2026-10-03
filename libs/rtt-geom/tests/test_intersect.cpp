#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <limits>
#include <utility>

#include "rtt/geom/conic.hpp"
#include "rtt/geom/intersect.hpp"
#include "rtt/geom/plane.hpp"

using rtt::geom::Conic;
using rtt::geom::HitStatus;
using rtt::geom::Intersection;
using rtt::geom::Plane;
using Vec3 = rtt::math::Vec3T<double>;

// Reference: sag formula of the conic (W. T. Welford, Aberrations of Optical Systems, Ch. 2).
// Sphere intersections are checked against the elementary ray-sphere solution
// |o + t d - C|^2 = R^2 with centre C = (0, 0, R), which does not use the conic code.

namespace {
constexpr double kTol = 1e-12;

bool near(double a, double b, double tol = kTol) {
  return std::abs(a - b) <= tol;
}

bool near(const Vec3& a, const Vec3& b, double tol = kTol) {
  return (a - b).cwiseAbs().maxCoeff() <= tol;
}

bool all_finite(const Intersection<double>& hit) {
  return std::isfinite(hit.t) && hit.point.allFinite() && hit.normal.allFinite();
}

/// Elementary ray-sphere intersection: smallest t > 0 on the cap that contains the vertex.
double sphere_reference_t(double radius, const Vec3& o, const Vec3& d) {
  const Vec3 centre(0.0, 0.0, radius);
  const Vec3 oc = o - centre;
  const double b = oc.dot(d);
  const double disc = b * b - (oc.squaredNorm() - radius * radius);
  const double t1 = -b - std::sqrt(disc);
  const double t2 = -b + std::sqrt(disc);
  for (const double t : {t1, t2}) {
    const Vec3 p = o + t * d;
    // The vertex cap is the hemisphere on the vertex side of the centre.
    const bool on_cap = radius > 0 ? p.z() <= radius : p.z() >= radius;
    if (t > 0 && on_cap) {
      return t;
    }
  }
  return std::numeric_limits<double>::quiet_NaN();
}
}  // namespace

TEST_CASE("axial ray hits the vertex at t = distance with normal +z", "[intersect]") {
  const Vec3 dir(0.0, 0.0, 1.0);
  for (const auto& [c, k] : {std::pair{0.02, 0.0}, std::pair{-0.02, 0.0}, std::pair{0.01, -1.0},
                             std::pair{0.02, -2.0}, std::pair{0.02, 0.5}, std::pair{0.0, 0.0}}) {
    const Conic<double> conic(c, k);
    for (const double distance : {1.0, 37.5}) {
      const auto hit = rtt::geom::intersect(conic, Vec3(0.0, 0.0, -distance), dir);
      REQUIRE(hit.status == HitStatus::Hit);
      REQUIRE(near(hit.t, distance));
      REQUIRE(near(hit.point, Vec3::Zero()));
      REQUIRE(near(hit.normal, Vec3::UnitZ()));
    }
  }
  const auto hit = rtt::geom::intersect(Plane<double>(), Vec3(0.0, 0.0, -12.0), dir);
  REQUIRE(hit.status == HitStatus::Hit);
  REQUIRE(near(hit.t, 12.0));
  REQUIRE(near(hit.normal, Vec3::UnitZ()));
}

TEST_CASE("normal points to +z at the vertex for a ray travelling in -z", "[intersect]") {
  const auto hit =
      rtt::geom::intersect(Conic<double>(0.02, 0.0), Vec3(0.0, 0.0, 5.0), Vec3(0.0, 0.0, -1.0));
  REQUIRE(hit.status == HitStatus::Hit);
  REQUIRE(near(hit.t, 5.0));
  REQUIRE(near(hit.normal, Vec3::UnitZ()));
}

TEST_CASE("sphere: hit lies on the surface and is the near solution on the vertex cap",
          "[intersect][sphere]") {
  for (const double radius : {50.0, -50.0, 8.0}) {
    const Conic<double> sphere(1.0 / radius, 0.0);
    for (const auto& [o, d] : {
             std::pair{Vec3(0.0, 5.0, -10.0), Vec3(0.0, 0.0, 1.0)},
             std::pair{Vec3(1.0, -2.0, -20.0), Vec3(0.05, 0.1, 1.0).normalized()},
             std::pair{Vec3(-3.0, 4.0, 15.0), Vec3(0.02, -0.03, -1.0).normalized()},
             std::pair{Vec3(0.0, -6.0, -1.0), Vec3(0.0, 0.3, 1.0).normalized()},
         }) {
      const auto hit = rtt::geom::intersect(sphere, o, d);
      REQUIRE(hit.status == HitStatus::Hit);
      REQUIRE(near(hit.point.z(), sphere.sag(hit.point.x(), hit.point.y())));
      REQUIRE(near(hit.t, sphere_reference_t(radius, o, d), 1e-10));
      REQUIRE(near(hit.point, o + hit.t * d, 1e-12));
      // Normal from the gradient: (-dz/dx, -dz/dy, 1) normalised (docs/architecture.md).
      const auto [gx, gy] = sphere.grad(hit.point.x(), hit.point.y());
      REQUIRE(near(hit.normal, Vec3(-gx, -gy, 1.0).normalized()));
      REQUIRE(near(hit.normal.norm(), 1.0));
    }
  }
}

TEST_CASE("sphere: the far sheet behind the centre is not part of the surface",
          "[intersect][sphere]") {
  const Conic<double> sphere(0.02, 0.0);  // R = 50, vertex cap z <= 50
  // Starts inside the sphere behind the vertex and travels +z: the only crossing of the full
  // sphere is at z = 100 on the far sheet, which the sag does not describe.
  const auto miss = rtt::geom::intersect(sphere, Vec3(0.0, 0.0, 10.0), Vec3(0.0, 0.0, 1.0));
  REQUIRE(miss.status == HitStatus::Missed);
  REQUIRE(all_finite(miss));
}

TEST_CASE("sphere: of two valid hits on the vertex cap the nearer one is taken",
          "[intersect][sphere]") {
  // R = 50, centre (0, 0, 50). The line z = 10, x = 0 crosses the cap at y = +-sqrt(50^2 - 40^2)
  // = +-30, both with z <= R. From y = -60 travelling +y: t = 30 (not 90).
  const Conic<double> sphere(0.02, 0.0);
  const auto hit = rtt::geom::intersect(sphere, Vec3(0.0, -60.0, 10.0), Vec3(0.0, 1.0, 0.0));
  REQUIRE(hit.status == HitStatus::Hit);
  REQUIRE(near(hit.t, 30.0));
  REQUIRE(near(hit.point, Vec3(0.0, -30.0, 10.0)));
  // Same line travelling -y from y = +60.
  const auto back = rtt::geom::intersect(sphere, Vec3(0.0, 60.0, 10.0), Vec3(0.0, -1.0, 0.0));
  REQUIRE(back.status == HitStatus::Hit);
  REQUIRE(near(back.t, 30.0));
  REQUIRE(near(back.point, Vec3(0.0, 30.0, 10.0)));
}

TEST_CASE("only solutions with t > epsilon in propagation direction count", "[intersect]") {
  const Conic<double> sphere(0.02, 0.0);
  // Surface lies behind the ray.
  const auto behind = rtt::geom::intersect(sphere, Vec3(0.0, 0.0, -10.0), Vec3(0.0, 0.0, -1.0));
  REQUIRE(behind.status == HitStatus::Missed);
  // Ray starts on the vertex of a concave-to-(-z) sphere: t = 0 is excluded, the other root of
  // the full sphere (z = -100) lies behind the ray.
  const auto on_surface =
      rtt::geom::intersect(Conic<double>(-0.02, 0.0), Vec3(0.0, 0.0, 0.0), Vec3(0.0, 0.0, 1.0));
  REQUIRE(on_surface.status == HitStatus::Missed);
  // Ray starts on the plane: t = 0 is excluded and there is no further crossing.
  const auto plane =
      rtt::geom::intersect(Plane<double>(), Vec3(0.0, 0.0, 0.0), Vec3(0.0, 0.6, 0.8));
  REQUIRE(plane.status == HitStatus::Missed);
}

TEST_CASE("ray missing the sphere gives Missed without NaN", "[intersect][sphere]") {
  const Conic<double> sphere(0.02, 0.0);
  const auto hit = rtt::geom::intersect(sphere, Vec3(0.0, 60.0, -10.0), Vec3(0.0, 0.0, 1.0));
  REQUIRE(hit.status == HitStatus::Missed);
  REQUIRE(all_finite(hit));
}

TEST_CASE("ray outside the domain of the conic gives Missed without NaN", "[intersect]") {
  // Ellipsoid: r_max = 1 / (c sqrt(1 + k)) = 40.82 mm.
  const Conic<double> ellipsoid(0.02, 0.5);
  const auto hit = rtt::geom::intersect(ellipsoid, Vec3(45.0, 0.0, -10.0), Vec3(0.0, 0.0, 1.0));
  REQUIRE(hit.status == HitStatus::Missed);
  REQUIRE(all_finite(hit));
  // Sphere: beyond r = R.
  const auto sphere_hit =
      rtt::geom::intersect(Conic<double>(0.1, 0.0), Vec3(0.0, 10.5, -1.0), Vec3(0.0, 0.0, 1.0));
  REQUIRE(sphere_hit.status == HitStatus::Missed);
  REQUIRE(all_finite(sphere_hit));
}

TEST_CASE("paraboloid: parallel ray hits at z = y^2 / (2 R)", "[intersect][paraboloid]") {
  const double radius = 100.0;
  const double c = 1.0 / radius;
  const Conic<double> parabola(c, -1.0);
  for (const double y : {0.5, 20.0, 80.0, 300.0}) {
    const Vec3 o(0.0, y, -5.0);
    const auto hit = rtt::geom::intersect(parabola, o, Vec3(0.0, 0.0, 1.0));
    REQUIRE(hit.status == HitStatus::Hit);
    const double z = y * y / (2.0 * radius);
    REQUIRE(near(hit.t, z + 5.0, 1e-12 * (1.0 + z)));
    REQUIRE(near(hit.point, Vec3(0.0, y, z), 1e-12 * (1.0 + z)));
    // n ~ (-dz/dx, -dz/dy, 1) = (0, -c y, 1)
    REQUIRE(near(hit.normal, Vec3(0.0, -c * y, 1.0).normalized()));
  }
}

TEST_CASE("paraboloid: oblique ray hit lies on the surface", "[intersect][paraboloid]") {
  const Conic<double> parabola(-0.01, -1.0);
  const Vec3 o(3.0, -2.0, 10.0);
  const Vec3 d = Vec3(0.2, 0.1, -1.0).normalized();
  const auto hit = rtt::geom::intersect(parabola, o, d);
  REQUIRE(hit.status == HitStatus::Hit);
  REQUIRE(hit.t > 0.0);
  REQUIRE(near(hit.point.z(), parabola.sag(hit.point.x(), hit.point.y())));
}

TEST_CASE("hyperboloid: parallel ray hits at the analytic sag", "[intersect][hyperboloid]") {
  // k = -2: z = (sqrt(1 + c^2 r^2) - 1) / c, see test_conic.cpp.
  const double c = 0.02;
  const Conic<double> hyperbola(c, -2.0);
  for (const double y : {1.0, 30.0, 250.0}) {
    const auto hit = rtt::geom::intersect(hyperbola, Vec3(0.0, y, -1.0), Vec3(0.0, 0.0, 1.0));
    REQUIRE(hit.status == HitStatus::Hit);
    const double z = (std::sqrt(1.0 + c * c * y * y) - 1.0) / c;
    REQUIRE(near(hit.point.z(), z, 1e-12 * (1.0 + z)));
    const auto [gx, gy] = hyperbola.grad(0.0, y);
    REQUIRE(near(hit.normal, Vec3(-gx, -gy, 1.0).normalized()));
  }
}

TEST_CASE("hyperboloid: of two valid hits on the vertex sheet the nearer one is taken",
          "[intersect][hyperboloid]") {
  // k = -2, c = 0.02: z = (sqrt(1 + c^2 y^2) - 1) / c = 10 gives y^2 = 1100. Both crossings of
  // the line z = 10 lie on the vertex sheet (w = 1 + c z > 0). From y = -60: t = 60 - sqrt(1100).
  const Conic<double> hyperbola(0.02, -2.0);
  const auto hit = rtt::geom::intersect(hyperbola, Vec3(0.0, -60.0, 10.0), Vec3(0.0, 1.0, 0.0));
  REQUIRE(hit.status == HitStatus::Hit);
  REQUIRE(near(hit.t, 60.0 - std::sqrt(1100.0)));
  REQUIRE(near(hit.point.z(), hyperbola.sag(hit.point.x(), hit.point.y())));
}

TEST_CASE("hyperboloid: ray parallel to the asymptote (linear case)", "[intersect][hyperboloid]") {
  // For |d| = 1 the quadratic coefficient is c (1 + k d_z^2); it vanishes for d_z^2 = -1 / k.
  // Here it is zero only up to rounding, so the root g / q carries the hit. The exact branch
  // a == 0 is reached by the axial paraboloid case (d = +z, k = -1) and the flat conic (c = 0).
  const double c = 0.02;
  const double k = -2.0;
  const Conic<double> hyperbola(c, k);
  const double dz = std::sqrt(-1.0 / k);
  const Vec3 d(0.0, std::sqrt(1.0 - dz * dz), dz);
  const Vec3 o(0.0, -10.0, -20.0);
  const auto hit = rtt::geom::intersect(hyperbola, o, d);
  REQUIRE(hit.status == HitStatus::Hit);
  REQUIRE(all_finite(hit));
  REQUIRE(near(hit.point.z(), hyperbola.sag(hit.point.x(), hit.point.y())));
  REQUIRE(near(hit.point, o + hit.t * d));
}

TEST_CASE("hyperboloid: oblique rays hit on the vertex sheet", "[intersect][hyperboloid]") {
  const Conic<double> hyperbola(-0.05, -3.0);
  for (const auto& [o, d] : {
           std::pair{Vec3(1.0, 2.0, 30.0), Vec3(0.1, -0.2, -1.0).normalized()},
           std::pair{Vec3(-5.0, 0.0, -30.0), Vec3(0.0, 0.0, 1.0)},
       }) {
    const auto hit = rtt::geom::intersect(hyperbola, o, d);
    REQUIRE(hit.status == HitStatus::Hit);
    REQUIRE(near(hit.point.z(), hyperbola.sag(hit.point.x(), hit.point.y())));
    REQUIRE(near(hit.normal.norm(), 1.0));
  }
}

TEST_CASE("plane: oblique ray, parallel ray and ray moving away", "[intersect][plane]") {
  const Plane<double> plane;
  const Vec3 d = Vec3(0.3, -0.4, 1.0).normalized();
  const auto hit = rtt::geom::intersect(plane, Vec3(1.0, 2.0, -4.0), d);
  REQUIRE(hit.status == HitStatus::Hit);
  REQUIRE(near(hit.t, 4.0 / d.z()));
  REQUIRE(near(hit.point.z(), 0.0));
  REQUIRE(near(hit.normal, Vec3::UnitZ()));

  const auto parallel = rtt::geom::intersect(plane, Vec3(0.0, 0.0, -1.0), Vec3(1.0, 0.0, 0.0));
  REQUIRE(parallel.status == HitStatus::Missed);
  REQUIRE(all_finite(parallel));

  const auto away = rtt::geom::intersect(plane, Vec3(0.0, 0.0, -1.0), Vec3(0.0, 0.0, -1.0));
  REQUIRE(away.status == HitStatus::Missed);
  REQUIRE(all_finite(away));
}

TEST_CASE("plane shape has zero sag and gradient and no domain limit", "[plane]") {
  const Plane<double> plane;
  REQUIRE(plane.sag(3.0, -7.0) == 0.0);
  REQUIRE(plane.grad(3.0, -7.0) == std::pair{0.0, 0.0});
  REQUIRE(plane.base_conic() == std::pair{0.0, 0.0});
  REQUIRE_FALSE(plane.max_radius().has_value());
}

TEST_CASE("shapes are usable through the Shape interface", "[shape]") {
  const Conic<double> sphere(0.02, 0.0);
  const rtt::geom::Shape<double>& shape = sphere;
  REQUIRE(near(shape.sag(0.0, 10.0), 50.0 - std::sqrt(2400.0)));
}

TEST_CASE("intersection works for float", "[intersect]") {
  const auto hit =
      rtt::geom::intersect(Conic<float>(0.02F, 0.0F), rtt::math::Vec3T<float>(0, 0, -3),
                           rtt::math::Vec3T<float>(0, 0, 1));
  REQUIRE(hit.status == HitStatus::Hit);
  REQUIRE(std::abs(hit.t - 3.0F) <= 1e-6F);
}
