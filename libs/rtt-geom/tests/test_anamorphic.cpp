#include <algorithm>
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>

#include "rtt/geom/anamorphic.hpp"
#include "rtt/geom/asphere.hpp"
#include "rtt/geom/conic.hpp"
#include "rtt/geom/intersect.hpp"
#include "rtt/geom/shape.hpp"

// Anamorphic shapes of ADR 0042, point 1: biconic (Ansys Zemax OpticStudio User Guide 2025 R1,
// "Biconic"), toroid (same guide, "Toroidal": a y-z profile rotated about an axis parallel to y;
// the surface formula is derived in ADR 0042) and cylinder (the toroid profile without rotation).

using rtt::geom::Biconic;
using rtt::geom::Conic;
using rtt::geom::Cylinder;
using rtt::geom::DomainShape;
using rtt::geom::EvenAsphere;
using rtt::geom::HitStatus;
using rtt::geom::IntersectableShape;
using rtt::geom::ParaxialShape;
using rtt::geom::Toroid;
using Vec3 = rtt::math::Vec3T<double>;

static_assert(IntersectableShape<Biconic<double>, double>);
static_assert(IntersectableShape<Toroid<double>, double>);
static_assert(IntersectableShape<Cylinder<double>, double>);
static_assert(ParaxialShape<Biconic<double>, double>);
static_assert(ParaxialShape<Toroid<double>, double>);
static_assert(ParaxialShape<Cylinder<double>, double>);
static_assert(DomainShape<Biconic<double>, double>);
static_assert(DomainShape<Toroid<double>, double>);
static_assert(DomainShape<Cylinder<double>, double>);
// The existing shapes keep the circular domain of max_radius() (ADR 0042, point 1).
static_assert(!DomainShape<Conic<double>, double>);
static_assert(!DomainShape<EvenAsphere<double>, double>);

namespace {

constexpr double kEps = std::numeric_limits<double>::epsilon();

bool near(double a, double b, double tol) {
  return std::abs(a - b) <= tol;
}

/// |a - b| <= ulps * eps * max(|a|, |b|, floor): a few units in the last place.
bool near_ulps(double a, double b, double ulps, double floor = 0.0) {
  return std::abs(a - b) <= ulps * kEps * std::max({std::abs(a), std::abs(b), floor});
}

/// Central-difference gradient of `shape` at (x, y), step h.
template <class S>
std::pair<double, double> numeric_grad(const S& shape, double x, double y, double h) {
  return {(shape.sag(x + h, y) - shape.sag(x - h, y)) / (2.0 * h),
          (shape.sag(x, y + h) - shape.sag(x, y - h)) / (2.0 * h)};
}

/// The y-z profile of the toroid and cylinder, written out independently of the library:
/// c y^2 / (1 + sqrt(1 - (1 + k) c^2 y^2)) + A4 y^4 + A6 y^6 + ...
double profile(double c, double k, const std::vector<double>& a, double y) {
  double z = c * y * y / (1.0 + std::sqrt(1.0 - (1.0 + k) * c * c * y * y));
  double power = y * y * y * y;
  for (const double ai : a) {
    z += ai * power;
    power *= y * y;
  }
  return z;
}

bool all_finite(const rtt::geom::Intersection<double>& hit) {
  return std::isfinite(hit.t) && hit.point.allFinite() && hit.normal.allFinite();
}

}  // namespace

TEST_CASE("biconic with equal axes equals the conic", "[anamorphic]") {
  // ADR 0042, acceptance 2: sag and gradient up to a few ulp. The biconic forms
  // c x^2 + c y^2 and (1 + k) c^2 x^2 + (1 + k) c^2 y^2 where the conic forms c (x^2 + y^2):
  // each extra rounding is at most 1/2 ulp of a term, so 8 ulp of the larger of value and
  // 1e-3 (a scale below every tested sag) bounds the difference.
  for (const auto& [c, k] : {std::pair{0.02, 0.0}, std::pair{-0.03, -1.0}, std::pair{0.01, -2.5},
                             std::pair{0.025, 0.6}}) {
    const Conic<double> conic(c, k);
    const Biconic<double> biconic(c, c, k, k);
    for (const auto& [x, y] : {std::pair{0.0, 0.0}, std::pair{3.0, -4.0}, std::pair{-12.0, 7.5},
                               std::pair{1e-3, 2e-3}}) {
      CHECK(near_ulps(biconic.sag(x, y), conic.sag(x, y), 8.0, 1e-3));
      const auto [bx, by] = biconic.grad(x, y);
      const auto [cx, cy] = conic.grad(x, y);
      CHECK(near_ulps(bx, cx, 8.0, 1e-3));
      CHECK(near_ulps(by, cy, 8.0, 1e-3));
    }
    // The vertex expansion is the conic's bit for bit (ADR 0042, point 1).
    const auto e = biconic.vertex_expansion();
    const auto f = conic.vertex_expansion();
    CHECK(e.value == f.value);
    CHECK(e.slope == f.slope);
    CHECK(e.hessian == f.hessian);
    CHECK(e.order == f.order);
    CHECK(e.rotational4);
    CHECK(e.conic == f.conic);
    CHECK(e.a4 == f.a4);
    // Same start value of the intersection and same domain.
    CHECK(biconic.base_conic() == conic.base_conic());
    CHECK(biconic.max_radius() == conic.max_radius());
  }
}

TEST_CASE("biconic of an ellipsoid satisfies the ellipsoid equation", "[anamorphic]") {
  // OpticStudio guide, "Biconic", section "Modeling an ellipsoid": the ellipsoid
  // x^2/a^2 + y^2/b^2 + (z - c)^2/c^2 = 1 with the vertex at the origin is the biconic with
  // c_x = c/a^2, c_y = c/b^2, 1 + k_x = a^2/c^2, 1 + k_y = b^2/c^2. Independent check of the
  // biconic sag on the lower half of the ellipsoid.
  const double a = 30.0;
  const double b = 12.0;
  const double c = 20.0;
  const Biconic<double> biconic(c / (a * a), c / (b * b), a * a / (c * c) - 1.0,
                                b * b / (c * c) - 1.0);
  int checked = 0;
  for (const double x : {0.0, 5.0, -17.0, 25.0}) {
    for (const double y : {0.0, 3.0, -7.0, 9.5}) {
      if (x * x / (a * a) + y * y / (b * b) > 0.95) continue;
      const double z = biconic.sag(x, y);
      REQUIRE(std::isfinite(z));
      const double lhs = x * x / (a * a) + y * y / (b * b) + (z - c) * (z - c) / (c * c);
      CHECK(near(lhs, 1.0, 1e-14));
      ++checked;
    }
  }
  REQUIRE(checked >= 12);
}

TEST_CASE("biconic gradient matches central differences", "[anamorphic]") {
  // ADR 0042, point 1: dz/dx = 2 c_x x / (1 + S) + a x N / (S (1 + S)^2), likewise in y.
  const Biconic<double> biconic(1.0 / 30.0, -1.0 / 80.0, -2.0, 0.7);
  for (const auto& [x, y] : {std::pair{3.0, 4.0}, std::pair{-7.0, 2.0}, std::pair{10.0, -9.0}}) {
    const auto [gx, gy] = biconic.grad(x, y);
    const auto [nx, ny] = numeric_grad(biconic, x, y, 1e-5);
    CHECK(near(gx, nx, 1e-9));
    CHECK(near(gy, ny, 1e-9));
    // Content guard: the two sections differ.
    CHECK(std::abs(gx) > 0.05);
  }
}

TEST_CASE("biconic domain is the ellipse, max_radius its major semi-axis", "[anamorphic]") {
  // a x^2 + b y^2 <= 1 with a = (1 + k_x) c_x^2, b = (1 + k_y) c_y^2 (ADR 0042, point 1).
  const double cx = 1.0 / 10.0;
  const double cy = 1.0 / 50.0;
  const Biconic<double> biconic(cx, cy, 0.0, 0.0);
  const double ax = 1.0 / cx;  // semi-axis in x: 10 mm
  const double ay = 1.0 / cy;  // semi-axis in y: 50 mm
  REQUIRE(biconic.max_radius().has_value());
  CHECK(near(*biconic.max_radius(), ay, 1e-12));
  CHECK(biconic.in_domain(0.0, 0.99 * ay));
  CHECK(biconic.in_domain(0.99 * ax, 0.0));
  CHECK_FALSE(biconic.in_domain(1.01 * ax, 0.0));
  CHECK_FALSE(biconic.in_domain(0.8 * ax, 0.8 * ay));
  // Every domain point lies within max_radius (sampled on the boundary).
  for (int i = 0; i < 64; ++i) {
    const double phi = 2.0 * 3.141592653589793 * i / 64.0;
    CHECK(std::hypot(ax * std::cos(phi), ay * std::sin(phi)) <= *biconic.max_radius() * 1.0000001);
  }
  // Inside the domain sag and gradient are finite; outside NaN.
  CHECK(std::isfinite(biconic.sag(0.9 * ax, 0.0)));
  CHECK(std::isnan(biconic.sag(0.8 * ax, 0.8 * ay)));
  // Unbounded when a section is not closed (hyperbolic in x).
  CHECK_FALSE(Biconic<double>(cx, cy, -2.0, 0.0).max_radius().has_value());
  // Start value: the section with the larger (1 + k) c^2, here x.
  CHECK(biconic.base_conic() == std::pair{cx, 0.0});
}

TEST_CASE("biconic vertex expansion", "[anamorphic]") {
  const Biconic<double> biconic(0.04, -0.01, -1.0, 2.0);
  const auto e = biconic.vertex_expansion();
  CHECK(e.value == 0.0);
  CHECK(e.slope == std::array<double, 2>{0.0, 0.0});
  CHECK(e.hessian == std::array<double, 3>{0.04, 0.0, -0.01});
  CHECK(e.order == 4);
  CHECK_FALSE(e.rotational4);
  // Equal curvatures but different conics: not rotational to fourth order.
  CHECK_FALSE(Biconic<double>(0.04, 0.04, -1.0, 0.0).vertex_expansion().rotational4);
}

TEST_CASE("toroid sag lies on the surface of revolution", "[anamorphic]") {
  // ADR 0042, acceptance 3: the point (x, y, z) has the distance |R - z_y(y)| from the axis of
  // rotation, the line parallel to y through (x, z) = (0, R).
  const double c = 1.0 / 40.0;
  const double k = -0.5;
  const std::vector<double> a{1e-6, -2e-9};
  for (const double r : {-25.0, 60.0, 15.0}) {
    const Toroid<double> toroid(c, k, a, 1.0 / r);
    int checked = 0;
    for (const auto& [x, y] : {std::pair{3.0, 4.0}, std::pair{-7.0, 2.0}, std::pair{10.0, -9.0},
                               std::pair{0.0, 6.0}, std::pair{5.0, 0.0}}) {
      const double z = toroid.sag(x, y);
      REQUIRE(std::isfinite(z));
      const double distance = std::hypot(x, z - r);
      const double expected = std::abs(r - profile(c, k, a, y));
      CHECK(std::abs(distance - expected) <= 1e-12 * expected);
      ++checked;
    }
    REQUIRE(checked == 5);
    // In the y-z plane the toroid is the profile.
    CHECK(near_ulps(toroid.sag(0.0, 7.0), profile(c, k, a, 7.0), 4.0));
  }
}

TEST_CASE("toroid with R = radius and k = 0 is the sphere", "[anamorphic]") {
  const double c = 1.0 / 30.0;
  const Toroid<double> toroid(c, 0.0, {}, c);
  const Conic<double> sphere(c, 0.0);
  for (const auto& [x, y] : {std::pair{3.0, 4.0}, std::pair{-12.0, 7.5}, std::pair{20.0, -15.0}}) {
    CHECK(near(toroid.sag(x, y), sphere.sag(x, y), 1e-12));
    const auto [tx, ty] = toroid.grad(x, y);
    const auto [sx, sy] = sphere.grad(x, y);
    CHECK(near(tx, sx, 1e-12));
    CHECK(near(ty, sy, 1e-12));
  }
  const auto e = toroid.vertex_expansion();
  CHECK(e.hessian == std::array<double, 3>{c, 0.0, c});
  CHECK(e.rotational4);
  CHECK(e.conic == 0.0);
  CHECK(e.a4 == 0.0);
  // max_radius is the outer extent: |R| for the sphere, up to the margin of ADR 0042.
  REQUIRE(toroid.max_radius().has_value());
  CHECK(*toroid.max_radius() >= 30.0);
  CHECK(*toroid.max_radius() <= 30.0 * (1.0 + 2e-6));
}

TEST_CASE("toroid gradient matches central differences", "[anamorphic]") {
  // ADR 0042, point 1: dz/dx = c' x / sqrt(1 - c'^2 x^2), dz/dy = z_y'(y) / sqrt(1 - c'^2 x^2).
  const Toroid<double> toroid(1.0 / 40.0, -0.5, {1e-6, -2e-9}, -1.0 / 25.0);
  for (const auto& [x, y] : {std::pair{3.0, 4.0}, std::pair{-7.0, 2.0}, std::pair{10.0, -9.0}}) {
    const auto [gx, gy] = toroid.grad(x, y);
    const auto [nx, ny] = numeric_grad(toroid, x, y, 1e-5);
    CHECK(near(gx, nx, 1e-9));
    CHECK(near(gy, ny, 1e-9));
  }
}

TEST_CASE("toroid vertex expansion, domain and start value", "[anamorphic]") {
  const double c = 1.0 / 40.0;
  const double cr = -1.0 / 25.0;
  const Toroid<double> toroid(c, -0.5, {1e-6}, cr);
  const auto e = toroid.vertex_expansion();
  CHECK(e.hessian == std::array<double, 3>{cr, 0.0, c});
  CHECK(e.order == 4);
  CHECK_FALSE(e.rotational4);
  // Sphere radius but a conic profile: not rotational.
  CHECK_FALSE(Toroid<double>(c, -0.5, {}, c).vertex_expansion().rotational4);
  // Domain in x: |x| <= |R - z_y(y)|.
  CHECK(toroid.in_domain(24.0, 0.0));
  CHECK_FALSE(toroid.in_domain(26.0, 0.0));
  CHECK(std::isnan(toroid.sag(26.0, 0.0)));
  // Start value: the section with the larger (1 + k) c^2, here the circle of rotation.
  CHECK(toroid.base_conic() == std::pair{cr, 0.0});
  CHECK(Toroid<double>(c, -0.5, {}, 1.0 / 400.0).base_conic() == std::pair{c, -0.5});
}

TEST_CASE("toroid domain ends where the profile meets the axis of rotation", "[anamorphic]") {
  // R = 20 and the parabolic profile z_y = y^2 / 20 towards +z: z_y reaches R at |y| = 20, where
  // the domain ends. The outer extent is the largest sqrt(y^2 + d^2) with d = R - z_y > 0; with
  // f(y^2) = y^2 + (20 - y^2/20)^2 it is 20 at both ends (y = 0 and y = 20) and smaller between.
  const Toroid<double> toroid(0.1, -1.0, {}, 1.0 / 20.0);
  REQUIRE(toroid.max_radius().has_value());
  CHECK(*toroid.max_radius() >= 20.0);
  CHECK(*toroid.max_radius() <= 20.0 * (1.0 + 2e-6) + 1e-9);
  // Beyond the crossing (z_y(y) > R) there is no surface.
  CHECK_FALSE(toroid.in_domain(0.0, 25.0));
}

TEST_CASE("cylinder is the profile, flat in x", "[anamorphic]") {
  const double c = -1.0 / 35.0;
  const double k = 0.3;
  const std::vector<double> a{2e-6, 1e-9};
  const Cylinder<double> cylinder(c, k, a);
  const EvenAsphere<double> asphere(c, k, a);
  for (const double x : {0.0, 17.0, -120.0}) {
    for (const double y : {0.0, 4.0, -11.0}) {
      CHECK(cylinder.sag(x, y) == asphere.sag(0.0, y));
      const auto [gx, gy] = cylinder.grad(x, y);
      CHECK(gx == 0.0);
      CHECK(gy == asphere.grad(0.0, y).second);
    }
  }
  CHECK_FALSE(cylinder.max_radius().has_value());
  CHECK(cylinder.in_domain(1e6, 0.0));
  CHECK_FALSE(cylinder.in_domain(0.0, 40.0));
  const auto e = cylinder.vertex_expansion();
  CHECK(e.hessian == std::array<double, 3>{0.0, 0.0, c});
  CHECK_FALSE(e.rotational4);
  CHECK(cylinder.base_conic() == std::pair{c, k});
  // The limit c_R -> 0 of the toroid.
  const Toroid<double> flat(c, k, a, 1e-12);
  CHECK(near(flat.sag(9.0, 6.0), cylinder.sag(9.0, 6.0), 1e-9));
}

TEST_CASE("Newton intersects strongly anamorphic biconics", "[anamorphic][newton]") {
  // R_x / R_y = 20 and crossed: hits lie on the surface (ADR 0042, point 1).
  for (const auto& biconic : {Biconic<double>(1.0 / 10.0, 1.0 / 200.0, 0.0, -0.5),
                              Biconic<double>(-1.0 / 200.0, 1.0 / 10.0, 1.0, 0.0)}) {
    int hits = 0;
    for (const double px : {-8.0, -3.0, 0.0, 2.5, 7.0}) {
      for (const double py : {-8.0, -1.0, 0.0, 4.0, 9.0}) {
        const Vec3 o(px, py, -20.0);
        const Vec3 d = Vec3(0.002 * px, -0.02, 1.0).normalized();
        const auto hit = rtt::geom::intersect(biconic, o, d);
        REQUIRE(hit.status == HitStatus::Hit);
        CHECK(near(hit.point.z(), biconic.sag(hit.point.x(), hit.point.y()), 1e-10));
        ++hits;
      }
    }
    REQUIRE(hits == 25);
  }
}

TEST_CASE("Newton uses the domain of the shape, not the disk of max_radius", "[anamorphic][newton]") {
  // Semi-axes 10 (x) and 50 (y). The point (0, 40) lies outside the inscribed disk of radius 10
  // but inside the ellipse: the ray hits (ADR 0042, acceptance 5).
  const Biconic<double> biconic(1.0 / 10.0, 1.0 / 50.0, 0.0, 0.0);
  const auto hit = rtt::geom::intersect(biconic, Vec3(0.0, 40.0, -10.0), Vec3(0.0, 0.0, 1.0));
  REQUIRE(hit.status == HitStatus::Hit);
  CHECK(near(hit.point.z(), biconic.sag(0.0, 40.0), 1e-12));
  // Outside the ellipse: Missed, without NaN.
  const auto miss = rtt::geom::intersect(biconic, Vec3(12.0, 0.0, -10.0), Vec3(0.0, 0.0, 1.0));
  CHECK(miss.status == HitStatus::Missed);
  CHECK(all_finite(miss));
}

TEST_CASE("Newton on a toroid strongly curved in x, steep ray", "[anamorphic][newton]") {
  // R = 10 in x, nearly flat profile. The ray crosses z = 0 at x = 12, outside the domain
  // |x| <= 10, but meets the surface at (6, 0, 2) on its way: the start from the circle of
  // rotation (the larger section) finds it.
  const Toroid<double> toroid(1e-3, 0.0, {}, 1.0 / 10.0);
  const Vec3 d = Vec3(-3.0, 0.0, 1.0).normalized();
  const Vec3 o = Vec3(12.0, 0.0, 0.0) - 5.0 * d;
  const auto hit = rtt::geom::intersect(toroid, o, d);
  REQUIRE(hit.status == HitStatus::Hit);
  CHECK(near(hit.point.x(), 6.0, 1e-9));
  CHECK(near(hit.point.z(), 2.0, 1e-9));
}

TEST_CASE("Newton on a cylinder and a toroid, grid of rays", "[anamorphic][newton]") {
  const Cylinder<double> cylinder(1.0 / 30.0, -1.0, {1e-6});
  const Toroid<double> toroid(1.0 / 30.0, -0.5, {1e-6}, -1.0 / 60.0);
  int hits = 0;
  for (const double px : {-10.0, 0.0, 6.0}) {
    for (const double py : {-10.0, 0.0, 8.0}) {
      const Vec3 o(px, py, -15.0);
      const Vec3 d = Vec3(0.03, -0.01, 1.0).normalized();
      const auto hc = rtt::geom::intersect(cylinder, o, d);
      const auto ht = rtt::geom::intersect(toroid, o, d);
      REQUIRE(hc.status == HitStatus::Hit);
      REQUIRE(ht.status == HitStatus::Hit);
      CHECK(near(hc.point.z(), cylinder.sag(hc.point.x(), hc.point.y()), 1e-10));
      CHECK(near(ht.point.z(), toroid.sag(ht.point.x(), ht.point.y()), 1e-10));
      ++hits;
    }
  }
  REQUIRE(hits == 9);
}
