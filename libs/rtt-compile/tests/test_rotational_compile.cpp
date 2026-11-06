// The odd asphere and the axicon in rtt-compile (#314, ADR 0043, point 4): the mapping of the
// model to the rtt-geom shapes and their profiles in the layout.

#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <optional>
#include <utility>
#include <variant>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/compile/layout.hpp"
#include "rtt/geom/axicon.hpp"
#include "rtt/geom/odd_asphere.hpp"
#include "rtt/material/material.hpp"
#include "rtt/math/units.hpp"
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
  s.name = "rotational";
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

}  // namespace

TEST_CASE("compile maps the odd asphere with its curvature and the axicon with its angle",
          "[compile][rotational]") {
  SECTION("odd asphere: c is the Param itself, not 1/R") {
    const CompiledSystem cs = compiled(lens(
        rtt::model::OddAsphere{Param(0.02), Param(-0.5), {Param(1e-3), Param(0.0), Param(2e-6)}},
        std::nullopt));
    const auto& a = std::get<rtt::geom::OddAsphere<double>>(cs.surfaces()[0].shape);
    CHECK(a.base_conic() == std::pair{0.02, -0.5});
    CHECK(a.coefficients() == std::vector<double>{1e-3, 0.0, 2e-6});
  }
  SECTION("odd asphere with c = 0: a plane base") {
    const CompiledSystem cs =
        compiled(lens(rtt::model::OddAsphere{Param(0.0), Param(0.0), {Param(0.1)}}, std::nullopt));
    const auto& a = std::get<rtt::geom::OddAsphere<double>>(cs.surfaces()[0].shape);
    CHECK(a.base_conic() == std::pair{0.0, 0.0});
    CHECK_FALSE(a.max_radius().has_value());
  }
  SECTION("axicon: angle_deg in rad") {
    const CompiledSystem cs = compiled(lens(rtt::model::Axicon{Param(-5.0)}, std::nullopt));
    const auto& a = std::get<rtt::geom::Axicon<double>>(cs.surfaces()[0].shape);
    CHECK(a.angle() == rtt::math::deg_to_rad(-5.0));
    CHECK(std::abs(a.angle() + 5.0 * std::numbers::pi / 180.0) <= 1e-17);
  }
  SECTION("a bound curvature takes the value of its row") {
    System s =
        lens(rtt::model::OddAsphere{Param::bound("C"), Param(0.0), {Param(1e-3)}}, std::nullopt);
    rtt::model::ParameterRow row;
    row.name = "C";
    row.form = -0.0125;
    s.parameters = {row};
    const CompiledSystem cs = compiled(s);
    const auto& a = std::get<rtt::geom::OddAsphere<double>>(cs.surfaces()[0].shape);
    CHECK(a.base_conic().first == -0.0125);
  }
}

TEST_CASE("layout: profiles of an odd asphere and an axicon are cut at the aperture",
          "[compile][rotational][layout]") {
  // Both domains are unbounded at c = 0 (axicon, plane base): the aperture ends the profile.
  const SectionPlane yz{Vec3::Zero(), Vec3(1.0, 0.0, 0.0)};
  {
    const CompiledSystem cs =
        compiled(lens(rtt::model::Axicon{Param(10.0)}, rtt::model::CircularAperture{12.0, 0.0}));
    const std::vector<Polyline> cut = rtt::compile::surface_profile(cs, 0, yz, 49);
    REQUIRE(cut.size() == 1);
    REQUIRE(cut[0].size() == 49);
    const double s = std::tan(rtt::math::deg_to_rad(10.0));
    for (const Vec3& p : cut[0]) {
      REQUIRE(p.allFinite());
      CHECK(std::abs(p.z() - s * std::abs(p.y())) <= 1e-14);
    }
    CHECK(std::abs(std::abs(cut[0].front().y()) - 12.0) <= 1e-12);
    CHECK(cut[0][24].z() == 0.0);  // the apex lies on the axis sample
  }
  {
    const CompiledSystem cs = compiled(lens(
        rtt::model::OddAsphere{Param(0.02), Param(0.0), {Param(1e-3), Param(0.0), Param(2e-6)}},
        rtt::model::CircularAperture{10.0, 0.0}));
    const std::vector<Polyline> cut = rtt::compile::surface_profile(cs, 0, yz, 21);
    REQUIRE(cut.size() == 1);
    const Vec3& rim = cut[0].back();
    CHECK(std::abs(rim.y() - 10.0) <= 1e-12);
    // z(10) = 0.02 * 100 / (1 + sqrt(1 - 0.04)) + 1e-3 * 10 + 2e-6 * 1000.
    const double z = 2.0 / (1.0 + std::sqrt(0.96)) + 0.01 + 0.002;
    CHECK(std::abs(rim.z() - z) <= 1e-14);
  }
}
