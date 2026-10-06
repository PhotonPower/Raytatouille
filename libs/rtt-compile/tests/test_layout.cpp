// Geometry export for the layout (#81, rtt/compile/layout.hpp).

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/compile/layout.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"

using rtt::compile::compile;
using rtt::compile::CompiledSystem;
using rtt::compile::Frame;
using rtt::compile::Polyline;
using rtt::compile::SectionPlane;
using rtt::material::MaterialLibrary;
using rtt::math::Vec3;
using rtt::model::Element;
using rtt::model::ElementKind;
using rtt::model::Param;
using rtt::model::Pose;
using rtt::model::Surface;
using rtt::model::SurfaceId;
using rtt::model::System;

namespace {

constexpr double kPi = std::numbers::pi;

System base_system() {
  System s;
  s.name = "layout";
  s.environment.medium = "VACUUM";
  s.wavelengths = {{0.5876, 1.0, true}};
  s.aperture = {rtt::model::SystemApertureType::EntrancePupilDiameter, Param(10.0)};
  s.fields = {rtt::model::FieldType::AngleDeg, {{0.0, 0.0, 1.0}}};
  s.root.name = "root";
  s.paths = {{"main", true, {}}};
  return s;
}

Surface surface(const std::string& id,
                double z_mm,
                rtt::model::BaseShape shape,
                std::optional<double> radius = 12.0) {
  Surface s;
  s.id = SurfaceId(id);
  s.pose = Pose::along_z(z_mm);
  s.shape.base = std::move(shape);
  if (radius) s.aperture = rtt::model::CircularAperture{*radius, 0.0};
  return s;
}

rtt::model::Conic conic(double radius, double k = 0.0) {
  return {Param(radius), Param(k)};
}

/// Closed form of the conic sag (Forbes 2011, Eq. (2.1)) with the even polynomial of
/// ISO 10110-12, as in rtt-geom (docs/quellen.md).
double sag_formula(double c, double k, const std::vector<double>& a, double r) {
  double z = c * r * r / (1.0 + std::sqrt(1.0 - (1.0 + k) * c * c * r * r));
  for (std::size_t i = 0; i < a.size(); ++i) z += a[i] * std::pow(r, 2.0 * i + 4.0);
  return z;
}

/// dz/dr of sag_formula().
double slope_formula(double c, double k, const std::vector<double>& a, double r) {
  double s = c * r / std::sqrt(1.0 - (1.0 + k) * c * c * r * r);
  for (std::size_t i = 0; i < a.size(); ++i) {
    s += (2.0 * i + 4.0) * a[i] * std::pow(r, 2.0 * i + 3.0);
  }
  return s;
}

std::uint32_t index_of(const CompiledSystem& cs, const std::string& id) {
  const auto i = cs.find_surface(SurfaceId(id));
  REQUIRE(i.has_value());
  return *i;
}

/// The y-z plane x = x0 (normal along x).
SectionPlane yz_plane(double x0 = 0.0) {
  return {Vec3(x0, 0.0, 0.0), Vec3(1.0, 0.0, 0.0)};
}

/// Signed area of a closed polygon projected onto the y-z plane (shoelace formula).
double yz_area(const Polyline& p) {
  double a = 0.0;
  for (std::size_t i = 0; i + 1 < p.size(); ++i) {
    a += p[i].y() * p[i + 1].z() - p[i + 1].y() * p[i].z();
  }
  return 0.5 * a;
}

/// True if two segments (a, b) and (c, d) in the y-z plane cross in their interiors.
bool segments_cross(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d) {
  const auto orient = [](const Vec3& p, const Vec3& q, const Vec3& r) {
    return (q.y() - p.y()) * (r.z() - p.z()) - (q.z() - p.z()) * (r.y() - p.y());
  };
  const double d1 = orient(a, b, c);
  const double d2 = orient(a, b, d);
  const double d3 = orient(c, d, a);
  const double d4 = orient(c, d, b);
  return ((d1 > 0.0 && d2 < 0.0) || (d1 < 0.0 && d2 > 0.0)) &&
         ((d3 > 0.0 && d4 < 0.0) || (d3 < 0.0 && d4 > 0.0));
}

/// Closed (last point = first) and simple in the y-z plane: no two non-adjacent edges cross, no
/// edge of zero length, and no edge runs back along the previous one (collinear overlap).
void require_closed_simple(const Polyline& p) {
  REQUIRE(p.size() >= 4);
  REQUIRE(p.front() == p.back());
  const std::size_t n = p.size() - 1;  // edges
  std::size_t degenerate = 0;
  for (std::size_t i = 0; i < n; ++i) {
    const Vec3 e1 = p[i + 1] - p[i];
    const Vec3 e2 = p[(i + 2) % n == 0 ? 1 : (i + 2) % n] - p[i + 1];
    const double cross = e1.y() * e2.z() - e1.z() * e2.y();
    if (e1.norm() == 0.0) ++degenerate;
    if (std::abs(cross) <= 1e-12 * e1.norm() * e2.norm() && e1.dot(e2) < 0.0) ++degenerate;
  }
  REQUIRE(degenerate == 0);
  std::size_t crossings = 0;
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = i + 2; j < n; ++j) {
      if (i == 0 && j == n - 1) continue;  // adjacent through the closing point
      if (segments_cross(p[i], p[i + 1], p[j], p[j + 1])) ++crossings;
    }
  }
  REQUIRE(crossings == 0);
}

/// Integral of the sphere sag R - sqrt(R^2 - y^2) (R > 0) over [-a, a]. With the primitive
/// F(y) = (y sqrt(R^2 - y^2) + R^2 asin(y / R)) / 2 of sqrt(R^2 - y^2):
/// I = 2 a R - (F(a) - F(-a)) = 2 a R - a sqrt(R^2 - a^2) - R^2 asin(a / R).
double sphere_sag_integral(double radius, double a) {
  return 2.0 * a * radius - a * std::sqrt(radius * radius - a * a) -
         radius * radius * std::asin(a / radius);
}

}  // namespace

TEST_CASE("sag and normals against the closed formulas", "[layout]") {
  // Sphere R = 50, conic R = -80 with k = -0.6 and even asphere R = 40, k = -1.2, A4 = 1e-5,
  // A6 = -2e-8 (Forbes 2011, Eq. (2.1); polynomial ISO 10110-12; docs/quellen.md), on a grid
  // of local points within
  // r <= 11.3 mm. The normal is (-dz/dx, -dz/dy, 1) / |...| with dz/dx = (dz/dr) x / r.
  // Tolerance 1e-12 mm for the sag, 1e-12 for the normal components.
  System s = base_system();
  Element lens{"L1",
               ElementKind::Lens,
               Pose::along_z(10.0),
               "CONST:1.5",
               {surface("S1", 0.0, conic(50.0)), surface("S2", 5.0, conic(-80.0, -0.6))}};
  Element asphere{
      "L2",
      ElementKind::Lens,
      Pose::along_z(30.0),
      "CONST:1.5",
      {surface("A1", 0.0,
               rtt::model::EvenAsphere{Param(40.0), Param(-1.2), {Param(1e-5), Param(-2e-8)}}),
       surface("A2", 4.0, rtt::model::Plane{})}};
  s.root.children = {{lens}, {asphere}};
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(s, lib);

  struct Case {
    std::string id;
    double c;
    double k;
    std::vector<double> a;
  };
  const std::vector<Case> cases = {{"S1", 1.0 / 50.0, 0.0, {}},
                                   {"S2", -1.0 / 80.0, -0.6, {}},
                                   {"A1", 1.0 / 40.0, -1.2, {1e-5, -2e-8}}};
  std::vector<double> x;
  std::vector<double> y;
  for (int i = -4; i <= 4; ++i) {
    for (int j = -4; j <= 4; ++j) {
      x.push_back(2.0 * i);
      y.push_back(2.0 * j);
    }
  }
  for (const Case& c : cases) {
    INFO("surface " << c.id);
    const std::uint32_t s_index = index_of(cs, c.id);
    std::vector<double> z(x.size());
    std::vector<double> n(3 * x.size());
    rtt::compile::surface_sag(cs, s_index, x, y, z);
    rtt::compile::surface_normal(cs, s_index, x, y, Frame::Local, n);
    for (std::size_t i = 0; i < x.size(); ++i) {
      const double r = std::hypot(x[i], y[i]);
      REQUIRE(std::abs(z[i] - sag_formula(c.c, c.k, c.a, r)) <= 1e-12);
      const double slope = slope_formula(c.c, c.k, c.a, r);
      const double gx = r > 0.0 ? slope * x[i] / r : 0.0;
      const double gy = r > 0.0 ? slope * y[i] / r : 0.0;
      const Vec3 expected = Vec3(-gx, -gy, 1.0).normalized();
      REQUIRE((Vec3(n[3 * i], n[3 * i + 1], n[3 * i + 2]) - expected).norm() <= 1e-12);
    }
  }
  std::vector<double> short_z(2);
  REQUIRE_THROWS_AS(rtt::compile::surface_sag(cs, 0, x, y, short_z), std::invalid_argument);
  REQUIRE_THROWS_AS(rtt::compile::surface_sag(cs, 99, x, y, short_z), std::out_of_range);
}

TEST_CASE("fold mirror and tilted element: global vertex and normal", "[layout]") {
  // Fold mirror at z = 20 rotated by 45 deg about x, and a lens at (0, 2, 30) rotated by 10 deg
  // about x with its second surface at local z = 5. With the right-handed R_x of
  // docs/architecture.md (Transformationen): the global vertex of the second lens surface is
  // (0, 2, 30) + R_x(10 deg) (0, 0, 5) = (0, 2 - 5 sin 10, 30 + 5 cos 10), and the global
  // normals at the vertices are R_x(a) z = (0, -sin a, cos a). Tolerance 1e-12.
  System s = base_system();
  Pose mirror_pose = Pose::along_z(20.0);
  mirror_pose.rotation_deg[0] = Param(45.0);
  Element mirror{"M",
                 ElementKind::Mirror,
                 mirror_pose,
                 std::nullopt,
                 {surface("M1", 0.0, rtt::model::Plane{})}};
  Pose lens_pose;
  lens_pose.position = {Param(0.0), Param(2.0), Param(30.0)};
  lens_pose.rotation_deg[0] = Param(10.0);
  Element lens{"L",
               ElementKind::Lens,
               lens_pose,
               "CONST:1.5",
               {surface("T1", 0.0, conic(60.0)), surface("T2", 5.0, conic(-60.0))}};
  s.root.children = {{mirror}, {lens}};
  s.paths = {{"main", false, {{SurfaceId("M1"), rtt::model::EventKind::Reflect, 0}}}};
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(s, lib);

  const auto vertex_normal = [&](std::uint32_t i) {
    std::vector<double> n(3);
    const std::vector<double> zero = {0.0};
    rtt::compile::surface_normal(cs, i, zero, zero, Frame::Global, n);
    return Vec3(n[0], n[1], n[2]);
  };
  const double a = 45.0 * kPi / 180.0;
  const double b = 10.0 * kPi / 180.0;
  const auto m1 = index_of(cs, "M1");
  const auto t2 = index_of(cs, "T2");
  REQUIRE((cs.surfaces()[m1].to_global.apply_point(Vec3::Zero()) - Vec3(0.0, 0.0, 20.0)).norm() <=
          1e-12);
  REQUIRE((vertex_normal(m1) - Vec3(0.0, -std::sin(a), std::cos(a))).norm() <= 1e-12);
  const Vec3 vertex(0.0, 2.0 - 5.0 * std::sin(b), 30.0 + 5.0 * std::cos(b));
  REQUIRE((cs.surfaces()[t2].to_global.apply_point(Vec3::Zero()) - vertex).norm() <= 1e-12);
  REQUIRE((vertex_normal(t2) - Vec3(0.0, -std::sin(b), std::cos(b))).norm() <= 1e-12);
  // The y-z plane contains the local z axes of both: profiles exist and lie in x = 0.
  for (const auto i : {m1, t2}) {
    const auto pieces = rtt::compile::surface_profile(cs, i, yz_plane(), 11);
    REQUIRE(pieces.size() == 1);
    for (const Vec3& p : pieces[0]) REQUIRE(std::abs(p.x()) <= 1e-12);
  }
}

TEST_CASE("singlet profile in the y-z plane", "[layout]") {
  // m1/singlet_const: L1.S1 conic R = 51.68 at z = 5, L1.S2 plane at z = 9, apertures r = 12.7.
  // In the plane x = 0 the profile of S1 is (0, y, 5 + sag(y)) for y from -12.7 to 12.7 in 51
  // equal steps (direction t = z x n = +y), S2 is (0, y, 9). Tolerance 1e-12.
  const System s =
      rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m1/singlet_const.rtt.json");
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(s, lib);
  const auto s1 = rtt::compile::surface_profile(cs, index_of(cs, "L1.S1"), yz_plane(), 51);
  const auto s2 = rtt::compile::surface_profile(cs, index_of(cs, "L1.S2"), yz_plane(), 51);
  REQUIRE(s1.size() == 1);
  REQUIRE(s2.size() == 1);
  REQUIRE(s1[0].size() == 51);
  for (std::size_t i = 0; i < 51; ++i) {
    const double y = -12.7 + 25.4 * static_cast<double>(i) / 50.0;
    INFO("point " << i);
    REQUIRE(
        (s1[0][i] - Vec3(0.0, y, 5.0 + sag_formula(1.0 / 51.68, 0.0, {}, std::abs(y)))).norm() <=
        1e-12);
    REQUIRE((s2[0][i] - Vec3(0.0, y, 9.0)).norm() <= 1e-12);
  }
}

TEST_CASE("offset section plane through a decentred surface, annulus and errors", "[layout]") {
  // A sphere R = 40 decentred by x = 2 mm with aperture r = 10. The plane x = 2 contains its
  // axis: y from -10 to 10. The plane x = 1 is 1 mm off the axis: the line meets the aperture at
  // y = +-sqrt(100 - 1), and z = sag(r) with r^2 = 1 + y^2. An annulus 3 <= r <= 10 cut through
  // the axis gives the two pieces [-10, -3] and [3, 10]. Tolerance 1e-12.
  System s = base_system();
  Surface decentred = surface("D", 0.0, conic(40.0), 10.0);
  decentred.pose.position[0] = Param(2.0);
  Surface ring = surface("R", 3.0, rtt::model::Plane{}, std::nullopt);
  ring.aperture = rtt::model::CircularAperture{10.0, 3.0};
  Element lens{"L", ElementKind::Lens, Pose::along_z(10.0), "CONST:1.5", {decentred, ring}};
  Element open{"O",
               ElementKind::ThinElement,
               Pose::along_z(40.0),
               std::nullopt,
               {surface("P", 0.0, rtt::model::Plane{}, std::nullopt)}};
  s.root.children = {{lens}, {open}};
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(s, lib);
  const auto d = index_of(cs, "D");
  const double c = 1.0 / 40.0;

  const auto through = rtt::compile::surface_profile(cs, d, yz_plane(2.0), 5);
  REQUIRE(through.size() == 1);
  REQUIRE((through[0].front() - Vec3(2.0, -10.0, 10.0 + sag_formula(c, 0.0, {}, 10.0))).norm() <=
          1e-12);
  REQUIRE((through[0][2] - Vec3(2.0, 0.0, 10.0)).norm() <= 1e-12);

  const auto offset = rtt::compile::surface_profile(cs, d, yz_plane(1.0), 5);
  REQUIRE(offset.size() == 1);
  const double edge = std::sqrt(99.0);
  REQUIRE((offset[0].back() - Vec3(1.0, edge, 10.0 + sag_formula(c, 0.0, {}, 10.0))).norm() <=
          1e-12);
  REQUIRE((offset[0][2] - Vec3(1.0, 0.0, 10.0 + sag_formula(c, 0.0, {}, 1.0))).norm() <= 1e-12);

  const auto pieces = rtt::compile::surface_profile(cs, index_of(cs, "R"), yz_plane(), 3);
  REQUIRE(pieces.size() == 2);
  REQUIRE((pieces[0].front() - Vec3(0.0, -10.0, 13.0)).norm() <= 1e-12);
  REQUIRE((pieces[0].back() - Vec3(0.0, -3.0, 13.0)).norm() <= 1e-12);
  REQUIRE((pieces[1].front() - Vec3(0.0, 3.0, 13.0)).norm() <= 1e-12);
  REQUIRE((pieces[1].back() - Vec3(0.0, 10.0, 13.0)).norm() <= 1e-12);

  // A plane beside the aperture gives no piece.
  REQUIRE(rtt::compile::surface_profile(cs, d, yz_plane(20.0), 5).empty());

  // Errors: a plane not parallel to the local z axis (with a hint to M10), too few samples,
  // an unbounded surface (plane without aperture), a zero normal, an unknown surface.
  const SectionPlane across{Vec3(0.0, 0.0, 10.0), Vec3(0.0, 0.0, 1.0)};
  REQUIRE_THROWS_WITH(rtt::compile::surface_profile(cs, d, across, 5),
                      Catch::Matchers::ContainsSubstring("M10"));
  REQUIRE_THROWS_AS(rtt::compile::surface_profile(cs, d, yz_plane(), 1), std::invalid_argument);
  REQUIRE_THROWS_AS(rtt::compile::surface_profile(cs, index_of(cs, "P"), yz_plane(), 5),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(rtt::compile::surface_profile(cs, d, {Vec3::Zero(), Vec3::Zero()}, 5),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(rtt::compile::surface_profile(cs, 99, yz_plane(), 5), std::out_of_range);
}

TEST_CASE("outlines of a singlet and a meniscus: closed, simple, analytic area", "[layout]") {
  // Singlet m1/singlet_const (S1 sphere R = 51.68 at z = 5, S2 plane at z = 9, a = 12.7): the
  // section between the curves has the area int (4 - sag(y)) dy = 8 a - I(R, a) with
  // I = int_{-a}^{a} (R - sqrt(R^2 - y^2)) dy. Meniscus S1 R = 30, S2 R = 40 at 3 mm, a = 10:
  // area = 6 a + I(40, a) - I(30, a). The polygon replaces each curve by chords (trapezoid
  // rule with h = 2a/2000 for 2001 samples), whose error is h^2/12 (z'(a) - z'(-a)) to leading
  // order (Euler-Maclaurin) with z'(a) = c a / sqrt(1 - c^2 a^2): 6.8e-6 mm^2 for the singlet,
  // 1.6e-6 mm^2 for the meniscus (the two curves partly cancel). Tolerance 1e-5 mm^2.
  {
    const System s =
        rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m1/singlet_const.rtt.json");
    const MaterialLibrary lib;
    const CompiledSystem cs = compile(s, lib);
    const auto outlines = rtt::compile::element_outlines(cs, 1, yz_plane(), 2001);
    REQUIRE(cs.elements()[1].name == "L1");
    REQUIRE(outlines.size() == 1);
    require_closed_simple(outlines[0]);
    const double a = 12.7;
    REQUIRE(std::abs(std::abs(yz_area(outlines[0])) - (8.0 * a - sphere_sag_integral(51.68, a))) <=
            1e-5);
  }
  {
    System s = base_system();
    s.root.children = {
        {Element{"M",
                 ElementKind::Lens,
                 Pose::along_z(0.0),
                 "CONST:1.5",
                 {surface("M1", 0.0, conic(30.0), 10.0), surface("M2", 3.0, conic(40.0), 10.0)}}}};
    const MaterialLibrary lib;
    const CompiledSystem cs = compile(s, lib);
    const auto outlines = rtt::compile::element_outlines(cs, 0, yz_plane(), 2001);
    REQUIRE(outlines.size() == 1);
    require_closed_simple(outlines[0]);
    const double a = 10.0;
    const double area = 6.0 * a + sphere_sag_integral(40.0, a) - sphere_sag_integral(30.0, a);
    REQUIRE(std::abs(std::abs(yz_area(outlines[0])) - area) <= 1e-5);
  }
}

TEST_CASE("outlines of a cemented doublet share the cemented profile", "[layout]") {
  // m2/achromat: one Lens element with three surfaces and two segments (N-BK7, F2): two closed
  // simple polygons; the profile of the cemented surface L1.S2 appears in both, bitwise.
  const System s = rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m2/achromat.rtt.json");
  MaterialLibrary lib;
  lib.add_catalog(std::string(RTT_CATALOG_DIR) + "/schott.agf");
  const CompiledSystem cs = compile(s, lib);
  const auto outlines = rtt::compile::element_outlines(cs, 1, yz_plane(), 101);
  REQUIRE(outlines.size() == 2);
  for (const Polyline& p : outlines) require_closed_simple(p);
  const auto cemented = rtt::compile::surface_profile(cs, index_of(cs, "L1.S2"), yz_plane(), 101);
  REQUIRE(cemented.size() == 1);
  for (const Vec3& point : cemented[0]) {
    for (const Polyline& p : outlines) {
      bool found = false;
      for (const Vec3& q : p) found = found || q == point;
      REQUIRE(found);
    }
  }
}

TEST_CASE("outline step between rims of different radius, annulus, non-segmented elements",
          "[layout]") {
  // Plano lens: S1 plane at z = 0 with r = 10, S2 plane at z = 4 with r = 12. Per side the edge
  // runs from the larger rim (+-12, 4) parallel to the axis of S1 to (+-12, 0), then along the
  // section line to the smaller rim (+-10, 0): the polygon contains the corners (+-12, 0) and
  // is the rectangle 24 x 4 mm, area 96 mm^2. An annulus 3 <= r <= 10 on both surfaces of a lens
  // of thickness 2 gives two outlines of 7 x 2 mm. A plate of one material (prism) and a mirror
  // give none. Tolerance 1e-12.
  System s = base_system();
  Element step{"S",
               ElementKind::Lens,
               Pose::along_z(0.0),
               "CONST:1.5",
               {surface("P1", 0.0, rtt::model::Plane{}, 10.0),
                surface("P2", 4.0, rtt::model::Plane{}, 12.0)}};
  Surface r1 = surface("R1", 0.0, rtt::model::Plane{}, std::nullopt);
  Surface r2 = surface("R2", 2.0, rtt::model::Plane{}, std::nullopt);
  r1.aperture = rtt::model::CircularAperture{10.0, 3.0};
  r2.aperture = rtt::model::CircularAperture{10.0, 3.0};
  Element ring{"R", ElementKind::Lens, Pose::along_z(10.0), "CONST:1.5", {r1, r2}};
  Element prism{"P",
                ElementKind::Plate,
                Pose::along_z(20.0),
                "CONST:1.5",
                {surface("Q1", 0.0, rtt::model::Plane{}), surface("Q2", 5.0, rtt::model::Plane{}),
                 surface("Q3", 10.0, rtt::model::Plane{})}};
  Element mirror{"M",
                 ElementKind::Mirror,
                 Pose::along_z(40.0),
                 "CONST:1.5",
                 {surface("M1", 0.0, rtt::model::Plane{})}};
  s.root.children = {{step}, {ring}, {prism}, {mirror}};
  s.paths = {{"main", false, {{SurfaceId("P1"), rtt::model::EventKind::Refract, 0}}}};
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(s, lib);

  const auto outlines = rtt::compile::element_outlines(cs, 0, yz_plane(), 5);
  REQUIRE(outlines.size() == 1);
  require_closed_simple(outlines[0]);
  REQUIRE(std::abs(std::abs(yz_area(outlines[0])) - 96.0) <= 1e-12);
  for (const double y : {-12.0, 12.0}) {
    bool corner = false;
    for (const Vec3& p : outlines[0]) corner = corner || (p - Vec3(0.0, y, 0.0)).norm() <= 1e-12;
    INFO("corner at y = " << y);
    REQUIRE(corner);
  }
  const auto rings = rtt::compile::element_outlines(cs, 1, yz_plane(), 5);
  REQUIRE(rings.size() == 2);
  for (const Polyline& p : rings) {
    require_closed_simple(p);
    REQUIRE(std::abs(std::abs(yz_area(p)) - 7.0 * 2.0) <= 1e-12);
  }
  REQUIRE(rtt::compile::element_outlines(cs, 2, yz_plane(), 5).empty());
  REQUIRE(rtt::compile::element_outlines(cs, 3, yz_plane(), 5).empty());
  REQUIRE_THROWS_AS(rtt::compile::element_outlines(cs, 9, yz_plane(), 5), std::out_of_range);
}

TEST_CASE("profiles of rectangular and elliptical apertures and at the shape domain", "[layout]") {
  // Plane surfaces with a rectangle (half widths 6 in x, 8 in y) and an ellipse (semi-axes 5 in
  // x, 9 in y), and a sphere R = 40 without aperture (domain r <= 40, the hemisphere).
  // - Rectangle: in the y-z plane y = -8 ... 8 (direction +y); in the plane y = 0 with normal
  //   (0, 1, 0) the direction is t = z x n = -x, so x runs from 6 to -6.
  // - Ellipse: y = -+9 through the axis; in the plane x = 3, y = -+9 sqrt(1 - 9/25) = -+7.2.
  // - Sphere: the profile ends 1e-12 R inside the domain (layout.hpp), at |y| = 40 (1 - 1e-12),
  //   with the sag of the formula there and within 40 sqrt(2e-12) = 5.7e-5 mm of z = 40.
  // Tolerance 1e-12. For the sag at the domain edge: dz/dr = r / sqrt(R^2 - r^2) = 7.1e5 there,
  // and |y| from the global point carries about one ulp of 40 (7.1e-15), so z = sag(|y|) agrees
  // to about 5e-9: tolerance 1e-8. (First set to 1e-9 without this product; the measured
  // 1.6e-9 showed the missing factor.)
  System s = base_system();
  Surface rect = surface("RE", 0.0, rtt::model::Plane{}, std::nullopt);
  rect.aperture = rtt::model::RectangularAperture{6.0, 8.0};
  Surface ellipse = surface("EL", 5.0, rtt::model::Plane{}, std::nullopt);
  ellipse.aperture = rtt::model::EllipticalAperture{5.0, 9.0};
  Element lens{"L",
               ElementKind::Lens,
               Pose::along_z(0.0),
               "CONST:1.5",
               {rect, ellipse, surface("SP", 50.0, conic(40.0), std::nullopt)}};
  s.root.children = {{lens}};
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(s, lib);

  const auto re_yz = rtt::compile::surface_profile(cs, 0, yz_plane(), 3);
  REQUIRE(re_yz.size() == 1);
  REQUIRE((re_yz[0].front() - Vec3(0.0, -8.0, 0.0)).norm() <= 1e-12);
  REQUIRE((re_yz[0].back() - Vec3(0.0, 8.0, 0.0)).norm() <= 1e-12);
  const SectionPlane xz{Vec3::Zero(), Vec3(0.0, 1.0, 0.0)};
  const auto re_xz = rtt::compile::surface_profile(cs, 0, xz, 3);
  REQUIRE(re_xz.size() == 1);
  REQUIRE((re_xz[0].front() - Vec3(6.0, 0.0, 0.0)).norm() <= 1e-12);
  REQUIRE((re_xz[0].back() - Vec3(-6.0, 0.0, 0.0)).norm() <= 1e-12);

  const auto el_axis = rtt::compile::surface_profile(cs, 1, yz_plane(), 3);
  REQUIRE((el_axis[0].front() - Vec3(0.0, -9.0, 5.0)).norm() <= 1e-12);
  REQUIRE((el_axis[0].back() - Vec3(0.0, 9.0, 5.0)).norm() <= 1e-12);
  const auto el_off = rtt::compile::surface_profile(cs, 1, yz_plane(3.0), 3);
  REQUIRE((el_off[0].front() - Vec3(3.0, -7.2, 5.0)).norm() <= 1e-12);
  REQUIRE((el_off[0].back() - Vec3(3.0, 7.2, 5.0)).norm() <= 1e-12);

  const auto sphere = rtt::compile::surface_profile(cs, 2, yz_plane(), 5);
  REQUIRE(sphere.size() == 1);
  const double edge = 40.0 * (1.0 - 1e-12);
  for (const Vec3& end : {sphere[0].front(), sphere[0].back()}) {
    REQUIRE(end.allFinite());
    REQUIRE(std::abs(std::abs(end.y()) - edge) <= 1e-12 * 40.0);
    REQUIRE(std::abs(end.z() - 50.0 - sag_formula(1.0 / 40.0, 0.0, {}, std::abs(end.y()))) <= 1e-8);
    REQUIRE(std::abs(end.z() - 90.0) <= 1e-4);
  }
}

TEST_CASE("outline edges: unequal bores, a turned-over surface, a tilted element, errors",
          "[layout]") {
  // 1. Annulus lens, thickness 4: front 3 <= r <= 10, back 2 <= r <= 10. The glass fills
  //    2 <= |y| <= 10 for 0 <= z <= 4 (bore at the smaller radius, flat shoulder at the front):
  //    two outlines of 8 x 4 = 32 mm^2 with the corners (+-2, 0).
  // 2. A singlet (sphere R = 50 at z = 0, plane at z = 4, r = 10) and the same with the plane
  //    turned over (rotated by 180 deg about x): identical geometry, equal areas.
  // 3. The step lens of the previous test (P1 r = 10 at 0, P2 r = 12 at 4) at (0, 2, 50),
  //    rotated by 10 deg about x: still 96 mm^2 (a rotation within the section plane), and the
  //    corners (0, +-12, 0) of P1's frame at (0, 2 +- 12 cos 10, 50 +- 12 sin 10).
  // 4. An annulus on the front surface only gives a different number of pieces: an error.
  // Tolerance 1e-12.
  System s = base_system();
  Surface b1 = surface("B1", 0.0, rtt::model::Plane{}, std::nullopt);
  Surface b2 = surface("B2", 4.0, rtt::model::Plane{}, std::nullopt);
  b1.aperture = rtt::model::CircularAperture{10.0, 3.0};
  b2.aperture = rtt::model::CircularAperture{10.0, 2.0};
  Element bores{"B", ElementKind::Lens, Pose::along_z(0.0), "CONST:1.5", {b1, b2}};
  Element plain{
      "S",
      ElementKind::Lens,
      Pose::along_z(20.0),
      "CONST:1.5",
      {surface("S1", 0.0, conic(50.0), 10.0), surface("S2", 4.0, rtt::model::Plane{}, 10.0)}};
  Surface turned_plane = surface("T2", 4.0, rtt::model::Plane{}, 10.0);
  turned_plane.pose.rotation_deg[0] = Param(180.0);
  Element turned{"T",
                 ElementKind::Lens,
                 Pose::along_z(30.0),
                 "CONST:1.5",
                 {surface("T1", 0.0, conic(50.0), 10.0), turned_plane}};
  Pose tilt;
  tilt.position = {Param(0.0), Param(2.0), Param(50.0)};
  tilt.rotation_deg[0] = Param(10.0);
  Element step{"P",
               ElementKind::Lens,
               tilt,
               "CONST:1.5",
               {surface("P1", 0.0, rtt::model::Plane{}, 10.0),
                surface("P2", 4.0, rtt::model::Plane{}, 12.0)}};
  Surface m1 = surface("M1", 0.0, rtt::model::Plane{}, std::nullopt);
  m1.aperture = rtt::model::CircularAperture{10.0, 3.0};
  Element mismatch{"X",
                   ElementKind::Lens,
                   Pose::along_z(80.0),
                   "CONST:1.5",
                   {m1, surface("M2", 4.0, rtt::model::Plane{}, 10.0)}};
  s.root.children = {{bores}, {plain}, {turned}, {step}, {mismatch}};
  s.paths = {{"main", false, {{SurfaceId("B1"), rtt::model::EventKind::Refract, 0}}}};
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(s, lib);
  const auto has_point = [](const Polyline& p, const Vec3& q) {
    for (const Vec3& x : p) {
      if ((x - q).norm() <= 1e-12) return true;
    }
    return false;
  };

  const auto bore_outlines = rtt::compile::element_outlines(cs, 0, yz_plane(), 5);
  REQUIRE(bore_outlines.size() == 2);
  for (const Polyline& o : bore_outlines) {
    require_closed_simple(o);
    REQUIRE(std::abs(std::abs(yz_area(o)) - 32.0) <= 1e-12);
  }
  REQUIRE(has_point(bore_outlines[0], Vec3(0.0, -2.0, 0.0)));
  REQUIRE(has_point(bore_outlines[1], Vec3(0.0, 2.0, 0.0)));

  const auto plain_outline = rtt::compile::element_outlines(cs, 1, yz_plane(), 101);
  const auto turned_outline = rtt::compile::element_outlines(cs, 2, yz_plane(), 101);
  REQUIRE(plain_outline.size() == 1);
  REQUIRE(turned_outline.size() == 1);
  require_closed_simple(turned_outline[0]);
  REQUIRE(std::abs(std::abs(yz_area(turned_outline[0])) - std::abs(yz_area(plain_outline[0]))) <=
          1e-12);

  const auto tilted = rtt::compile::element_outlines(cs, 3, yz_plane(), 5);
  REQUIRE(tilted.size() == 1);
  require_closed_simple(tilted[0]);
  REQUIRE(std::abs(std::abs(yz_area(tilted[0])) - 96.0) <= 1e-12);
  const double b = 10.0 * kPi / 180.0;
  for (const double sign : {-1.0, 1.0}) {
    const Vec3 corner(0.0, 2.0 + sign * 12.0 * std::cos(b), 50.0 + sign * 12.0 * std::sin(b));
    INFO("corner " << sign);
    REQUIRE(has_point(tilted[0], corner));
  }

  REQUIRE_THROWS_AS(rtt::compile::element_outlines(cs, 4, yz_plane(), 5), std::invalid_argument);
}

TEST_CASE("elements and the media in front of and behind each surface", "[layout]") {
  // #81: in the surface order of the element, from the environment, with the rules of ADR 0017:
  // lens segments, inside/environment toggle for a one-material plate (prism: env|glass,
  // glass|env, env|glass) and a mirror with substrate (env|substrate), environment on both
  // sides of a stop.
  System s = base_system();
  Element stop{"STOP",
               ElementKind::Stop,
               Pose::along_z(0.0),
               std::nullopt,
               {surface("STO", 0.0, rtt::model::Plane{})}};
  Element doublet{"D",
                  ElementKind::Lens,
                  Pose::along_z(5.0),
                  std::nullopt,
                  {surface("D1", 0.0, conic(50.0)), surface("D2", 4.0, conic(-40.0)),
                   surface("D3", 6.0, rtt::model::Plane{})}};
  doublet.segment_materials = {"CONST:1.5", "CONST:1.7"};
  Element prism{"P",
                ElementKind::Plate,
                Pose::along_z(20.0),
                "CONST:1.6",
                {surface("Q1", 0.0, rtt::model::Plane{}), surface("Q2", 5.0, rtt::model::Plane{}),
                 surface("Q3", 10.0, rtt::model::Plane{})}};
  Element mirror{"M",
                 ElementKind::Mirror,
                 Pose::along_z(40.0),
                 "CONST:1.8",
                 {surface("M1", 0.0, rtt::model::Plane{})}};
  s.root.children = {{stop}, {doublet}, {prism}, {mirror}};
  s.paths = {{"main", false, {{SurfaceId("STO"), rtt::model::EventKind::Transmit, 0}}}};
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(s, lib);
  const auto medium = [&](const std::string& reference) {
    for (std::uint32_t i = 0; i < cs.media().size(); ++i) {
      if (cs.media()[i].reference == reference) return i;
    }
    FAIL("no medium " << reference);
    return std::uint32_t{0};
  };
  const std::uint32_t env = cs.environment_medium();
  const std::uint32_t g15 = medium("CONST:1.5");
  const std::uint32_t g17 = medium("CONST:1.7");
  const std::uint32_t g16 = medium("CONST:1.6");
  const std::uint32_t g18 = medium("CONST:1.8");
  struct Expected {
    std::string id;
    std::uint32_t front;
    std::uint32_t back;
  };
  for (const Expected& e : std::vector<Expected>{{"STO", env, env},
                                                 {"D1", env, g15},
                                                 {"D2", g15, g17},
                                                 {"D3", g17, env},
                                                 {"Q1", env, g16},
                                                 {"Q2", g16, env},
                                                 {"Q3", env, g16},
                                                 {"M1", env, g18}}) {
    INFO("surface " << e.id);
    const auto& surf = cs.surfaces()[index_of(cs, e.id)];
    REQUIRE(surf.medium_front == e.front);
    REQUIRE(surf.medium_back == e.back);
  }
  REQUIRE(cs.elements().size() == 4);
  const auto& d = cs.elements()[1];
  REQUIRE(d.name == "D");
  REQUIRE(d.kind == ElementKind::Lens);
  REQUIRE(d.first_surface == index_of(cs, "D1"));
  REQUIRE(d.surface_count == 3);
  REQUIRE(d.media == std::vector<std::uint32_t>{g15, g17});
  REQUIRE(d.segmented);
  REQUIRE_FALSE(cs.elements()[2].segmented);
  REQUIRE(cs.elements()[0].media.empty());
  REQUIRE(cs.surfaces()[index_of(cs, "D2")].element == 1);
}
