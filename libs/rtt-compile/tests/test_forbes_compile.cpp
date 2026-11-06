// Q-con and Q-bfs in rtt-compile (#314, ADR 0043, point 4): the mapping of the model to the
// rtt-geom shapes and their profiles in the layout (acceptance 11).

#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <optional>
#include <utility>
#include <variant>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/compile/layout.hpp"
#include "rtt/geom/forbes.hpp"
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
  s.name = "forbes";
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

TEST_CASE("compile maps Q-con and Q-bfs with their curvature and normalization radius",
          "[compile][forbes]") {
  const CompiledSystem a = compiled(
      lens(rtt::model::QCon{Param(0.02), Param(-0.5), Param(12.5), {Param(1e-3), Param(-2e-4)}},
           std::nullopt));
  const auto& q = std::get<rtt::geom::QCon<double>>(a.surfaces()[0].shape);
  CHECK(q.base_conic() == std::pair{0.02, -0.5});
  CHECK(q.normalization_radius() == 12.5);
  CHECK(q.coefficients() == std::vector<double>{1e-3, -2e-4});
  const CompiledSystem b = compiled(
      lens(rtt::model::QBfs{Param(-0.01), Param(0.25), Param(8.0), {Param(3e-3)}}, std::nullopt));
  const auto& r = std::get<rtt::geom::QBfs<double>>(b.surfaces()[0].shape);
  CHECK(r.base_conic() == std::pair{-0.01, 0.25});
  CHECK(r.normalization_radius() == 8.0);
  CHECK(r.coefficients() == std::vector<double>{3e-3});
}

TEST_CASE("layout: a Q-bfs profile ends inside phi > 0 without NaN", "[compile][forbes][layout]") {
  // ADR 0043, acceptance 11: no aperture, c = 1/20, k = 0: the domain is r < 20 mm, where 1/sigma
  // diverges; the profile ends 1e-12 relative inside (max_radius and in_domain, ADR 0042).
  const CompiledSystem cs = compiled(
      lens(rtt::model::QBfs{Param(0.05), Param(0.0), Param(15.0), {Param(1e-3), Param(-4e-4)}},
           std::nullopt));
  const SectionPlane yz{Vec3::Zero(), Vec3(1.0, 0.0, 0.0)};
  for (const std::size_t samples : {std::size_t{5}, std::size_t{64}}) {
    const std::vector<Polyline> cut = rtt::compile::surface_profile(cs, 0, yz, samples);
    REQUIRE(cut.size() == 1);
    for (const Vec3& p : cut[0]) REQUIRE(p.allFinite());
    for (const Vec3& end : {cut[0].front(), cut[0].back()}) {
      CHECK(std::abs(end.y()) < 20.0);
      CHECK(std::abs(std::abs(end.y()) - 20.0) <= 1e-10 * 20.0);
    }
  }
}

TEST_CASE("layout: Q-con beyond rho extrapolates; c = 0 with a conic is flat",
          "[compile][forbes][layout]") {
  // ADR 0043, point 2 (F3) and acceptance 11: aperture 10 mm, rho = 5 mm.
  const std::vector<double> a = {1e-3, -2e-4, 5e-5};
  const CompiledSystem cs = compiled(lens(
      rtt::model::QCon{
          Param(0.02), Param(0.0), Param(5.0), {Param(a[0]), Param(a[1]), Param(a[2])}},
      rtt::model::CircularAperture{10.0, 0.0}));
  const SectionPlane yz{Vec3::Zero(), Vec3(1.0, 0.0, 0.0)};
  const std::vector<Polyline> cut = rtt::compile::surface_profile(cs, 0, yz, 21);
  REQUIRE(cut.size() == 1);
  const rtt::geom::QCon<double> shape(0.02, 0.0, 5.0, a);
  for (const Vec3& p : cut[0]) {
    REQUIRE(p.allFinite());
    CHECK(std::abs(p.z() - shape.sag(0.0, p.y())) <= 1e-14);
  }
  CHECK(std::abs(cut[0].back().y() - 10.0) <= 1e-12);  // u = 2
  // c = 0 with k = -1 and zero coefficients: a plane.
  const CompiledSystem flat =
      compiled(lens(rtt::model::QCon{Param(0.0), Param(-1.0), Param(5.0), {Param(0.0)}},
                    rtt::model::CircularAperture{10.0, 0.0}));
  for (const Vec3& p : rtt::compile::surface_profile(flat, 0, yz, 11)[0]) CHECK(p.z() == 0.0);
}
