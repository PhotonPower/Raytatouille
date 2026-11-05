// The anamorphic base shapes in rtt-compile (#315, ADR 0042, points 1 and 5): the mapping of the
// model to the rtt-geom shapes and the profiles of the layout at a domain that is not a disk.

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <optional>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/compile/layout.hpp"
#include "rtt/geom/anamorphic.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"

using rtt::compile::compile;
using rtt::compile::CompiledSystem;
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

/// A lens of CONST:1.5 in vacuum: `front` at z = 0 with `aperture`, a plane back at z = 30.
System lens(const rtt::model::BaseShape& front, std::optional<rtt::model::Aperture> aperture) {
  System s;
  s.name = "anamorphic";
  s.environment.medium = "VACUUM";
  s.wavelengths = {{0.5876, 1.0, true}};
  s.aperture = {rtt::model::SystemApertureType::EntrancePupilDiameter, Param(4.0)};
  s.fields = {rtt::model::FieldType::AngleDeg, {{0.0, 0.0, 1.0}}};
  s.root.name = "root";
  s.paths = {{"main", true, {}}};
  Surface s1;
  s1.id = SurfaceId("A.S1");
  s1.shape.base = front;
  s1.aperture = std::move(aperture);
  Surface s2;
  s2.id = SurfaceId("A.S2");
  s2.pose = Pose::along_z(30.0);
  s2.aperture = rtt::model::CircularAperture{80.0, 0.0};
  s.root.children = {{Element{"A", ElementKind::Lens, Pose::along_z(0.0), "CONST:1.5", {s1, s2}}}};
  return s;
}

CompiledSystem compiled(const System& s) {
  const MaterialLibrary lib;
  return compile(s, lib);
}

void require_finite(const std::vector<Polyline>& pieces) {
  for (const Polyline& p : pieces) {
    for (const Vec3& v : p) REQUIRE(v.allFinite());
  }
}

}  // namespace

TEST_CASE("compile maps the anamorphic shapes with c = 1/R", "[compile][anamorphic]") {
  SECTION("biconic") {
    const CompiledSystem cs = compiled(lens(
        rtt::model::Biconic{Param(40.0), Param(-60.0), Param(-0.5), Param(0.25)}, std::nullopt));
    const auto& b = std::get<rtt::geom::Biconic<double>>(cs.surfaces()[0].shape);
    CHECK(b.curvature_x() == 1.0 / 40.0);
    CHECK(b.curvature_y() == 1.0 / -60.0);
    CHECK(b.conic_x() == -0.5);
    CHECK(b.conic_y() == 0.25);
  }
  SECTION("toroid") {
    const CompiledSystem cs = compiled(lens(
        rtt::model::Toroid{Param(48.0), Param(-1.2), {Param(2e-6), Param(-3e-9)}, Param(-55.0)},
        rtt::model::CircularAperture{10.0, 0.0}));
    const auto& t = std::get<rtt::geom::Toroid<double>>(cs.surfaces()[0].shape);
    CHECK(t.curvature_of_rotation() == 1.0 / -55.0);
    CHECK(t.profile().base_conic() == std::pair{1.0 / 48.0, -1.2});
    CHECK(t.profile().coefficients() == std::vector<double>{2e-6, -3e-9});
  }
  SECTION("cylinder") {
    const CompiledSystem cs =
        compiled(lens(rtt::model::Cylinder{Param(-75.0), Param(0.3), {Param(1e-6)}},
                      rtt::model::CircularAperture{10.0, 0.0}));
    const auto& c = std::get<rtt::geom::Cylinder<double>>(cs.surfaces()[0].shape);
    CHECK(c.profile().base_conic() == std::pair{1.0 / -75.0, 0.3});
    CHECK(c.profile().coefficients() == std::vector<double>{1e-6});
  }
}

TEST_CASE("layout: the profile of a biconic ends on its elliptic domain without NaN",
          "[compile][anamorphic][layout]") {
  // Biconic R_x = 10, R_y = 200, k = 0, no aperture: domain x^2/100 + y^2/40000 <= 1, outer
  // extent (max_radius) 200 (ADR 0042, point 1). The disk of max_radius alone would let the
  // profile in the x-z plane run to |x| = 200, where the sag is NaN beyond |x| = 10; the layout
  // splits the piece where in_domain changes and bisects the boundary to 1e-12 relative
  // ("Layout"). So the profile ends at |x| = 10 (1 - O(1e-12)), relative tolerance 1e-11.
  // In the plane x = 5: |y| <= 200 sqrt(1 - 25/100) = 100 sqrt(3).
  // With 5 samples only the centre sample of [-200, 200] lies inside the domain; with 64 more.
  const CompiledSystem cs = compiled(
      lens(rtt::model::Biconic{Param(10.0), Param(200.0), Param(0.0), Param(0.0)}, std::nullopt));
  for (const std::size_t samples : {std::size_t{5}, std::size_t{64}}) {
    INFO(samples);
    const SectionPlane xz{Vec3::Zero(), Vec3(0.0, 1.0, 0.0)};
    const std::vector<Polyline> x_cut = rtt::compile::surface_profile(cs, 0, xz, samples);
    REQUIRE(x_cut.size() == 1);
    REQUIRE(x_cut[0].size() == samples);
    require_finite(x_cut);
    for (const Vec3& end : {x_cut[0].front(), x_cut[0].back()}) {
      CHECK(std::abs(std::abs(end.x()) - 10.0) <= 1e-11 * 10.0);
      CHECK(std::abs(end.x()) <= 10.0);
      // On the rim of the domain the surface is vertical at z = R_x = 10.
      CHECK(std::abs(end.z() - 10.0) <= 1e-4);
    }
    CHECK(x_cut[0].front().x() * x_cut[0].back().x() < 0.0);

    const SectionPlane offset{Vec3(5.0, 0.0, 0.0), Vec3(1.0, 0.0, 0.0)};
    const std::vector<Polyline> y_cut = rtt::compile::surface_profile(cs, 0, offset, samples);
    REQUIRE(y_cut.size() == 1);
    require_finite(y_cut);
    const double edge = 100.0 * std::sqrt(3.0);
    for (const Vec3& end : {y_cut[0].front(), y_cut[0].back()}) {
      CHECK(std::abs(std::abs(end.y()) - edge) <= 1e-11 * edge);
      CHECK(end.x() == 5.0);
    }
  }
}

TEST_CASE("layout: a domain that cuts a hole into the section gives two pieces",
          "[compile][anamorphic][layout]") {
  // Biconic R_x = 10, k_x = -3 (a = (1 + k_x) c_x^2 = -0.02), R_y = 20, k_y = 0
  // (b = 1/400), circular aperture r = 60. Domain 1 - a x^2 - b y^2 >= 0. On the line y = 25:
  // b y^2 = 1.5625, so the domain is x^2 >= 0.5625 / 0.02 = 28.125, a hole of half width
  // sqrt(28.125) around x = 0; the aperture allows |x| <= sqrt(3600 - 625). The profile has
  // two pieces whose inner ends lie on the hole's rim (relative tolerance 1e-11), no NaN.
  const CompiledSystem cs =
      compiled(lens(rtt::model::Biconic{Param(10.0), Param(20.0), Param(-3.0), Param(0.0)},
                    rtt::model::CircularAperture{60.0, 0.0}));
  const SectionPlane cut{Vec3(0.0, 25.0, 0.0), Vec3(0.0, 1.0, 0.0)};
  const std::vector<Polyline> pieces = rtt::compile::surface_profile(cs, 0, cut, 64);
  REQUIRE(pieces.size() == 2);
  require_finite(pieces);
  const double hole = std::sqrt(28.125);
  const double rim = std::sqrt(3600.0 - 625.0);
  for (const Polyline& p : pieces) {
    REQUIRE(p.size() == 64);
    const double a = std::abs(p.front().x());
    const double b = std::abs(p.back().x());
    const double inner = std::min(a, b);
    const double outer = std::max(a, b);
    CHECK(std::abs(inner - hole) <= 1e-11 * hole);
    CHECK(inner >= hole * (1.0 - 1e-11));
    CHECK(std::abs(outer - rim) <= 1e-12 * rim);
    for (const Vec3& v : p) CHECK(std::abs(v.y() - 25.0) <= 1e-12);
  }
}
