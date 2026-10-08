#include <array>
#include <catch2/catch_test_macros.hpp>
#include <limits>
#include <string>
#include <variant>
#include <vector>

#include "test_support.hpp"

using namespace rtt::model;
using rtt::model::test::element;
using rtt::model::test::has_error_at;
using rtt::model::test::make_singlet;

TEST_CASE("reference singlet is valid", "[validate]") {
  const auto d = validate(make_singlet());
  INFO((d.empty() ? std::string() : to_string(d.front())));
  REQUIRE(d.empty());
}

TEST_CASE("wavelength checks", "[validate]") {
  System s = make_singlet();
  s.wavelengths.push_back({0.4861, 1.0, true});
  REQUIRE(has_error_at(validate(s), "/wavelengths"));  // two references

  s = make_singlet();
  s.wavelengths[0].um = -1.0;
  REQUIRE(has_error_at(validate(s), "/wavelengths/0/um"));

  s = make_singlet();
  s.wavelengths.clear();
  REQUIRE(has_error_at(validate(s), "/wavelengths"));
}

TEST_CASE("aperture, field and object checks", "[validate]") {
  System s = make_singlet();
  s.aperture.value = Param(0.0);
  REQUIRE(has_error_at(validate(s), "/aperture/value"));

  s = make_singlet();
  s.fields.points.clear();
  REQUIRE(has_error_at(validate(s), "/fields/points"));

  s = make_singlet();
  s.object.at_infinity = false;
  s.object.distance = Param(-5.0);
  REQUIRE(has_error_at(validate(s), "/object/distance"));
}

TEST_CASE("element kind constraints", "[validate]") {
  System s = make_singlet();
  element(s, 1).material.reset();
  REQUIRE(has_error_at(validate(s), "/root/children/1/material"));

  s = make_singlet();
  element(s, 1).surfaces.pop_back();
  REQUIRE(has_error_at(validate(s), "/root/children/1/surfaces"));

  s = make_singlet();
  element(s, 0).surfaces[0].aperture.reset();
  REQUIRE(has_error_at(validate(s), "/root/children/0/surfaces/0/aperture"));

  s = make_singlet();
  element(s, 2).material = "SCHOTT:N-BK7";
  REQUIRE(has_error_at(validate(s), "/root/children/2/material"));
}

TEST_CASE("surface checks", "[validate]") {
  System s = make_singlet();
  element(s, 1).surfaces[1].id = SurfaceId("L1.S1");
  REQUIRE(has_error_at(validate(s), "/root/children/1/surfaces/1/id"));

  s = make_singlet();
  element(s, 1).surfaces[0].shape.base = Conic{Param(0.0), Param(0.0)};
  REQUIRE(has_error_at(validate(s), "/root/children/1/surfaces/0/shape/base/radius"));

  s = make_singlet();
  element(s, 1).surfaces[0].shape.base =
      Conic{Param(std::numeric_limits<double>::infinity()), Param(0.0)};
  REQUIRE(has_error_at(validate(s), "/root/children/1/surfaces/0/shape/base/radius"));

  s = make_singlet();
  element(s, 1).surfaces[0].aperture = CircularAperture{5.0, 6.0};
  REQUIRE(has_error_at(validate(s), "/root/children/1/surfaces/0/aperture/inner_radius"));

  s = make_singlet();
  element(s, 1).surfaces[0].interaction = IdealBeamSplitter{1.5, 0.5};
  REQUIRE(has_error_at(validate(s), "/root/children/1/surfaces/0/interaction"));

  s = make_singlet();
  element(s, 1).surfaces[0].interaction = IdealPolarizer{{0.0, 0.0, 0.0}, 0.0};
  REQUIRE(has_error_at(validate(s), "/root/children/1/surfaces/0/interaction/transmission_axis"));
}

TEST_CASE("only one stop", "[validate]") {
  System s = make_singlet();
  Element second = element(s, 0);
  second.name = "stop2";
  second.surfaces[0].id = SurfaceId("STO2");
  s.root.children.push_back({second});
  REQUIRE(has_error_at(validate(s), "/root/children/3"));
}

TEST_CASE("node names are unique", "[validate]") {
  System s = make_singlet();
  element(s, 2).name = "L1";
  REQUIRE(has_error_at(validate(s), "/root/children/2/name"));
}

TEST_CASE("path checks", "[validate]") {
  System s = make_singlet();
  s.paths = {
      {"main",
       false,
       {{SurfaceId("L1.S1"), EventKind::Refract, 0}, {SurfaceId("nope"), EventKind::Refract, 0}}}};
  REQUIRE(has_error_at(validate(s), "/paths/0/events/1/surface"));

  s = make_singlet();
  s.paths = {{"main", false, {{SurfaceId("L1.S1"), EventKind::Reflect, 2}}}};
  REQUIRE(has_error_at(validate(s), "/paths/0/events/0/order"));

  s = make_singlet();
  s.paths = {{"main", false, {}}};
  REQUIRE(has_error_at(validate(s), "/paths/0/events"));

  s = make_singlet();
  s.paths.push_back(s.paths[0]);
  REQUIRE(has_error_at(validate(s), "/paths/1/name"));
}

TEST_CASE("nested assemblies are traversed", "[validate]") {
  System s = make_singlet();
  Assembly group;
  group.name = "group";
  group.children.push_back(s.root.children[1]);
  s.root.children[1] = {group};
  REQUIRE(validate(s).empty());
  std::get<Element>(std::get<Assembly>(s.root.children[1].value).children[0].value)
      .material.reset();
  REQUIRE(has_error_at(validate(s), "/root/children/1/children/0/material"));
}

TEST_CASE("an order is allowed at every event of a surface with a phase layer (ADR 0025)",
          "[validate]") {
  System s = make_singlet();
  element(s, 1).surfaces[0].phases.push_back(RadialPhase{Param(10.0), {Param(1.0)}});
  for (const EventKind kind :
       {EventKind::Refract, EventKind::Reflect, EventKind::Transmit, EventKind::Ordinary}) {
    s.paths = {{"main", false, {{SurfaceId("L1.S1"), kind, -2}}}};
    const auto d = validate(s);
    INFO((d.empty() ? std::string() : to_string(d.front())));
    REQUIRE(d.empty());
  }
  // The other surface of the lens has no phase layer.
  s.paths = {{"main", false, {{SurfaceId("L1.S2"), EventKind::Transmit, 1}}}};
  REQUIRE(has_error_at(validate(s), "/paths/0/events/0/order"));
}

TEST_CASE("diffraction efficiency checks (ADR 0025)", "[validate]") {
  const std::string loc = "/root/children/1/surfaces/0/diffraction_efficiency";
  System s = make_singlet();
  Surface& surface = element(s, 1).surfaces[0];
  surface.phases.push_back(LinearGrating{Param(300.0), 0.0});
  surface.diffraction_efficiency = std::vector<DiffractionEfficiency>{{1, 0.8}, {0, 0.0}};
  REQUIRE(validate(s).empty());

  surface.diffraction_efficiency = std::vector<DiffractionEfficiency>{};
  REQUIRE(has_error_at(validate(s), loc));  // empty: it would block every order

  surface.diffraction_efficiency = std::vector<DiffractionEfficiency>{{1, 0.8}, {1, 0.1}};
  REQUIRE(has_error_at(validate(s), loc + "/1/order"));

  for (const double eta : {-0.1, 1.1, std::numeric_limits<double>::quiet_NaN()}) {
    surface.diffraction_efficiency = std::vector<DiffractionEfficiency>{{1, eta}};
    REQUIRE(has_error_at(validate(s), loc + "/0/efficiency"));
  }

  surface.diffraction_efficiency = std::vector<DiffractionEfficiency>{{1, 1.0}};
  surface.phases.clear();
  REQUIRE(has_error_at(validate(s), loc));  // no phase layer
}

TEST_CASE("crystal checks (ADR 0026)", "[validate]") {
  System s = make_singlet();
  Element& lens = element(s, 1);
  lens.material.reset();
  lens.crystal = CrystalMaterial{"BIREFRINGENT:CALCITE", "BIREFRINGENT:CALCITE-E"};
  lens.optic_axis = std::array<double, 3>{0.0, 1.0, 1.0};
  REQUIRE(validate(s).empty());
  REQUIRE_FALSE(lens.segment_material(0).has_value());  // no isotropic material

  lens.kind = ElementKind::Plate;
  lens.surfaces[0].shape.base = Plane{};
  REQUIRE(validate(s).empty());

  lens.crystal->extraordinary.clear();
  REQUIRE(has_error_at(validate(s), "/root/children/1/material/extraordinary"));
  lens.crystal->extraordinary = "BIREFRINGENT:CALCITE-E";

  lens.optic_axis = std::array<double, 3>{0.0, std::numeric_limits<double>::infinity(), 1.0};
  REQUIRE(has_error_at(validate(s), "/root/children/1/optic_axis"));
  lens.optic_axis.reset();
  REQUIRE(has_error_at(validate(s), "/root/children/1/optic_axis"));

  // Without crystal, the plate without material is still an error, and an axis is not allowed.
  lens.crystal.reset();
  REQUIRE(has_error_at(validate(s), "/root/children/1/material"));
  lens.material = "SCHOTT:N-BK7";
  lens.optic_axis = std::array<double, 3>{0.0, 0.0, 1.0};
  REQUIRE(has_error_at(validate(s), "/root/children/1/optic_axis"));
}
