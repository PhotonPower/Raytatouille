#include <catch2/catch_test_macros.hpp>
#include <limits>

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
