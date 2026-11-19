// Regular polygon and polygon apertures in the tracer (#320, ADR 0050, point 2 and point 4,
// PR 3): transmitted fraction of a hexagon, sign of the rotation, equality with the vertex form
// and with the rectangle, a concave polygon, and the rim tolerance.

#include <algorithm>
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numbers>
#include <random>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/geom/conic.hpp"
#include "rtt/math/isometry.hpp"
#include "rtt/trace/apply_event.hpp"

using rtt::compile::CompiledSurface;
using rtt::math::Vec3;
using rtt::model::Aperture;
using rtt::model::Param;
using rtt::model::PolygonAperture;
using rtt::model::RectangularAperture;
using rtt::model::RegularPolygonAperture;
using rtt::trace::inside_aperture;

namespace {

constexpr double kPi = std::numbers::pi;

/// inside_aperture() for a hit at the local point (x, y) of a plane surface with `aperture`.
bool inside(const Aperture& aperture, double x, double y) {
  CompiledSurface s;
  s.shape = rtt::geom::Plane<double>{};
  s.aperture = aperture;
  const rtt::trace::SurfaceHit hit{rtt::geom::HitStatus::Hit, 0.0, Vec3(x, y, 0.0), Vec3::UnitZ(),
                                   Vec3::UnitZ()};
  return inside_aperture(s, hit);
}

/// A uniform double in [0, 1) from 53 random bits (ADR 0050, point 4: not
/// std::uniform_real_distribution, whose result depends on the standard library).
double unit(std::mt19937_64& g) {
  return static_cast<double>(g() >> 11U) * 0x1.0p-53;
}

/// Fraction of N points, uniform in the disc of radius `r_disc`, that pass `aperture`.
/// Hergeleitet in #320: the radius has the distribution F(r) = r^2 / R^2, so r = R sqrt(u).
double transmitted(const Aperture& aperture, double r_disc, std::size_t n) {
  std::mt19937_64 g(20261118U);
  std::size_t pass = 0;
  for (std::size_t i = 0; i < n; ++i) {
    const double r = r_disc * std::sqrt(unit(g));
    const double phi = 2.0 * kPi * unit(g);
    if (inside(aperture, r * std::cos(phi), r * std::sin(phi))) ++pass;
  }
  return static_cast<double>(pass) / static_cast<double>(n);
}

/// Distance of (x, y) from the boundary of a polygon given by its vertices.
double boundary_distance(const std::vector<std::array<double, 2>>& v, double x, double y) {
  double best = std::numeric_limits<double>::infinity();
  for (std::size_t i = 0, j = v.size() - 1; i < v.size(); j = i++) {
    const double ex = v[i][0] - v[j][0];
    const double ey = v[i][1] - v[j][1];
    const double px = x - v[j][0];
    const double py = y - v[j][1];
    const double s = std::clamp((px * ex + py * ey) / (ex * ex + ey * ey), 0.0, 1.0);
    best = std::min(best, std::hypot(px - s * ex, py - s * ey));
  }
  return best;
}

}  // namespace

TEST_CASE("hexagon: uniform disc of the circumradius passes 3 sqrt(3) / (2 pi) (ADR 0050)",
          "[aperture][polygon]") {
  // A_6 = (1/2) 6 R^2 sin 60 deg (MathWorld "Regular Polygon", eq. 9), so p = A_6 / (pi R^2) =
  // 3 sqrt(3) / (2 pi). The hits are binomial: sigma = sqrt(p (1 - p) / N) = 3.78e-4 for
  // N = 1e6, and the bound is 5 sigma (hergeleitet in #320).
  constexpr double kR = 12.0;
  constexpr std::size_t kN = 1000000;
  const double p = 3.0 * std::sqrt(3.0) / (2.0 * kPi);
  const double bound = 5.0 * std::sqrt(p * (1.0 - p) / static_cast<double>(kN));
  CHECK(bound < 1.9e-3);
  for (const double theta0 : {0.0, 30.0}) {
    INFO("rotation " << theta0 << " deg");
    const double f = transmitted(RegularPolygonAperture{6, Param(kR), theta0}, kR, kN);
    INFO("fraction " << f << ", expected " << p);
    CHECK(std::abs(f - p) <= bound);
  }
  // Content guard: a circle of the same radius passes everything.
  CHECK(transmitted(rtt::model::CircularAperture{Param(kR), Param(0.0)}, kR, 10000) == 1.0);
}

TEST_CASE("pentagon at 10 deg: the sign of the rotation (ADR 0050)", "[aperture][polygon]") {
  // Vertices at 10, 82, 154, 226 and 298 deg. A = 0.999 R (cos 10, sin 10) lies just inside the
  // first vertex; its mirror image B = 0.999 R (cos 10, -sin 10) lies outside: the edge from 298
  // to 10 deg has the normal psi = -26 deg, and B projects to 0.999 R cos 16 deg = 0.960 R >
  // r = R cos 36 deg = 0.809 R. With the opposite sign of rotation it would be the other way.
  constexpr double kR = 5.0;
  const RegularPolygonAperture pentagon{5, Param(kR), 10.0};
  const double a = 10.0 * kPi / 180.0;
  CHECK(inside(pentagon, 0.999 * kR * std::cos(a), 0.999 * kR * std::sin(a)));
  CHECK_FALSE(inside(pentagon, 0.999 * kR * std::cos(a), -0.999 * kR * std::sin(a)));
  // The hand numbers: the projection of B on psi = -26 deg against the inradius.
  CHECK(0.999 * std::cos(16.0 * kPi / 180.0) > std::cos(36.0 * kPi / 180.0));
  // The vertices of the model helper follow the same convention.
  const std::vector<std::array<double, 2>> v = rtt::model::aperture_vertices(pentagon);
  REQUIRE(v.size() == 5);
  CHECK(std::abs(v[0][0] - kR * std::cos(a)) < 1e-15);
  CHECK(std::abs(v[0][1] - kR * std::sin(a)) < 1e-15);
  CHECK(v[1][1] > 0.0);  // 82 deg: counterclockwise
}

TEST_CASE("regular hexagon = polygon of its vertices; 4-vertex polygon = rectangle (ADR 0050)",
          "[aperture][polygon]") {
  // A 201 x 201 grid over [-1.1 R, 1.1 R]^2 (and the rectangle's box), leaving out a band of
  // 1e-6 mm about the edges, where the tolerance shapes differ at the corners.
  const RegularPolygonAperture hexagon{6, Param(10.0), 17.0};
  const std::vector<std::array<double, 2>> hv = rtt::model::aperture_vertices(hexagon);
  const PolygonAperture hexagon_vertices{hv};
  const RectangularAperture rectangle{Param(4.0), Param(2.5)};
  const std::vector<std::array<double, 2>> rv = {
      {-4.0, -2.5}, {4.0, -2.5}, {4.0, 2.5}, {-4.0, 2.5}};
  const PolygonAperture rectangle_vertices{rv};
  std::size_t hex_in = 0;
  std::size_t rect_in = 0;
  for (int i = 0; i <= 200; ++i) {
    for (int j = 0; j <= 200; ++j) {
      const double u = -1.1 + 2.2 * i / 200.0;
      const double w = -1.1 + 2.2 * j / 200.0;
      if (boundary_distance(hv, 10.0 * u, 10.0 * w) > 1e-6) {
        const bool a = inside(hexagon, 10.0 * u, 10.0 * w);
        CHECK(a == inside(hexagon_vertices, 10.0 * u, 10.0 * w));
        hex_in += a ? 1 : 0;
      }
      if (boundary_distance(rv, 4.4 * u, 2.75 * w) > 1e-6) {
        const bool a = inside(rectangle, 4.4 * u, 2.75 * w);
        CHECK(a == inside(rectangle_vertices, 4.4 * u, 2.75 * w));
        rect_in += a ? 1 : 0;
      }
    }
  }
  // Content guard: both sets are neither empty nor everything.
  CHECK(hex_in > 1000);
  CHECK(hex_in < 201 * 201);
  CHECK(rect_in > 1000);
  CHECK(rect_in < 201 * 201);
}

TEST_CASE("a concave L-shaped polygon blocks its missing corner (ADR 0050)",
          "[aperture][polygon]") {
  // The L: the square [0, 4]^2 without [2, 4]^2, counterclockwise and clockwise.
  std::vector<std::array<double, 2>> l = {{0.0, 0.0}, {4.0, 0.0}, {4.0, 2.0},
                                          {2.0, 2.0}, {2.0, 4.0}, {0.0, 4.0}};
  for (int pass = 0; pass < 2; ++pass) {
    const PolygonAperture shape{l};
    CHECK(inside(shape, 1.0, 1.0));
    CHECK(inside(shape, 3.0, 1.0));
    CHECK(inside(shape, 1.0, 3.0));
    CHECK_FALSE(inside(shape, 3.0, 3.0));  // the missing corner
    CHECK_FALSE(inside(shape, 5.0, 1.0));
    // The rims belong to the aperture, also at the reflex vertex and inside the tolerance.
    CHECK(inside(shape, 2.0, 2.0));
    CHECK(inside(shape, 3.0, 2.0 + 0.5e-9));
    CHECK_FALSE(inside(shape, 3.0, 2.0 + 2e-9));
    CHECK(inside(shape, 4.0 + 0.5e-9, 1.0));
    CHECK_FALSE(inside(shape, 4.0 + 2e-9, 1.0));
    std::reverse(l.begin(), l.end());
  }
}

TEST_CASE("regular polygon rims include the tolerance as a distance outwards (ADR 0050)",
          "[aperture][polygon]") {
  // Hexagon at 0 deg: horizontal edges at y = +-r, r = R cos 30 deg.
  const RegularPolygonAperture hexagon{6, Param(2.0), 0.0};
  const double r = 2.0 * std::cos(kPi / 6.0);
  CHECK(inside(hexagon, 0.0, r + 0.5 * rtt::trace::kApertureTolerance));
  CHECK_FALSE(inside(hexagon, 0.0, r + 2.0 * rtt::trace::kApertureTolerance));
  CHECK(inside(hexagon, 2.0, 0.0));  // a vertex
  CHECK_FALSE(inside(hexagon, 2.0 + 1e-6, 0.0));
}
