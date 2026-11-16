// Ranges of Params bound to the parameter table (#353, ADR 0029): validate() skips bound values,
// so compile checks them after the resolution with the same codes and JSON pointers as unbound
// values. Valid bound values compile as the same literal values would.

#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/compile/errors.hpp"
#include "rtt/geom/axicon.hpp"
#include "rtt/geom/conic.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"

using rtt::compile::CompiledSystem;
using rtt::compile::CompileError;
using rtt::material::MaterialLibrary;
using rtt::model::Element;
using rtt::model::ElementKind;
using rtt::model::Param;
using rtt::model::Pose;
using rtt::model::Surface;
using rtt::model::SurfaceId;
using rtt::model::System;

namespace {

const std::string kFront = "/root/children/0/surfaces/0";

/// A lens of CONST:1.5: `front` at z = 0, a plane back at z = 5; EPD 4 mm, field 0 deg.
System lens(const rtt::model::BaseShape& front) {
  System s;
  s.name = "bound ranges";
  s.wavelengths = {{0.5876, 1.0, true}};
  s.aperture = {rtt::model::SystemApertureType::EntrancePupilDiameter, Param(4.0)};
  s.fields = {rtt::model::FieldType::AngleDeg, {{0.0, 0.0, 1.0}}};
  s.root.name = "root";
  s.paths = {{"main", true, {}}};
  Surface s1;
  s1.id = SurfaceId("A.S1");
  s1.shape.base = front;
  Surface s2;
  s2.id = SurfaceId("A.S2");
  s2.pose = Pose::along_z(5.0);
  s.root.children = {{Element{"A", ElementKind::Lens, Pose::along_z(0.0), "CONST:1.5", {s1, s2}}}};
  return s;
}

/// Two configurations "a" and "b" and the row `name` with the values a and b.
void table(System& s, const std::string& name, double a, double b) {
  s.configurations = {{"a"}, {"b"}};
  rtt::model::ParameterRow row;
  row.name = name;
  row.form = std::vector<double>{a, b};
  s.parameters.push_back(row);
}

/// The errors of compile() in `configuration`; fails if compile() does not throw.
std::vector<rtt::model::Diagnostic> errors(const System& s, std::size_t configuration) {
  const MaterialLibrary lib;
  try {
    (void)rtt::compile::compile(s, lib, configuration);
  } catch (const CompileError& e) {
    return e.diagnostics();
  }
  FAIL("compile() did not throw");
  return {};
}

/// Exactly one diagnostic, an error with `code` at `pointer`.
bool only(const std::vector<rtt::model::Diagnostic>& d,
          const std::string& code,
          const std::string& pointer) {
  return d.size() == 1 && d[0].code == code && d[0].location == pointer &&
         d[0].severity == rtt::model::Severity::Error;
}

CompiledSystem compiled(const System& s, std::size_t configuration) {
  const MaterialLibrary lib;
  return rtt::compile::compile(s, lib, configuration);
}

}  // namespace

TEST_CASE("a bound axicon angle of 95 deg is shape.axicon_angle_invalid (#353)", "[bound]") {
  System s = lens(rtt::model::Axicon{Param::bound("A")});
  table(s, "A", 5.0, 95.0);
  CHECK_NOTHROW(compiled(s, 0));
  CHECK(only(errors(s, 1), "shape.axicon_angle_invalid", kFront + "/shape/base/angle_deg"));
}

TEST_CASE("a bound radius 0 is shape.radius_invalid (#353)", "[bound]") {
  System s = lens(rtt::model::Conic{Param::bound("R"), Param(0.0)});
  table(s, "R", 0.0, 50.0);
  CHECK(only(errors(s, 0), "shape.radius_invalid", kFront + "/shape/base/radius"));
  CHECK_NOTHROW(compiled(s, 1));
}

TEST_CASE("a rule over several Params checks the resolved values (#353)", "[bound]") {
  // Q-bfs needs (1 + k) c^2 rho^2 < 1; the curvature is bound, the error stays at the unbound
  // normalization radius, as without the table. c = 0.2, rho = 10: 4 >= 1.
  System s = lens(rtt::model::QBfs{Param::bound("C"), Param(0.0), Param(10.0), {Param(1e-3)}});
  table(s, "C", 0.01, 0.2);
  CHECK_NOTHROW(compiled(s, 0));
  CHECK(only(errors(s, 1), "shape.normalization_radius_invalid",
             kFront + "/shape/base/normalization_radius"));
}

TEST_CASE("bound values outside the shapes are checked as well (#353)", "[bound]") {
  SECTION("lines per mm of a grating") {
    System s = lens(rtt::model::Plane{});
    auto& front = std::get<Element>(s.root.children[0].value).surfaces[0];
    front.phases.emplace_back(rtt::model::LinearGrating{Param::bound("G"), 0.0});
    table(s, "G", 100.0, -100.0);
    CHECK_NOTHROW(compiled(s, 0));
    CHECK(only(errors(s, 1), "phase.lines_per_mm_invalid", kFront + "/phases/0/lines_per_mm"));
  }
  SECTION("value of the system aperture") {
    System s = lens(rtt::model::Plane{});
    s.aperture.value = Param::bound("EPD");
    table(s, "EPD", 4.0, 0.0);
    CHECK_NOTHROW(compiled(s, 0));
    CHECK(only(errors(s, 1), "aperture.value_invalid", "/aperture/value"));
  }
  SECTION("finite object distance") {
    System s = lens(rtt::model::Plane{});
    s.object.at_infinity = false;
    s.object.distance = Param::bound("OBJ");
    table(s, "OBJ", -10.0, 100.0);
    CHECK(only(errors(s, 0), "object.distance_invalid", "/object/distance"));
    CHECK_NOTHROW(compiled(s, 1));
  }
}

TEST_CASE("valid bound shape values compile as the literal values (#353)", "[bound]") {
  System bound_axicon = lens(rtt::model::Axicon{Param::bound("A")});
  table(bound_axicon, "A", 5.0, -30.0);
  System bound_conic = lens(rtt::model::Conic{Param::bound("R"), Param(-0.5)});
  table(bound_conic, "R", 40.0, -25.0);
  for (const std::size_t k : {std::size_t{0}, std::size_t{1}}) {
    const System literal_axicon = lens(rtt::model::Axicon{Param(k == 0 ? 5.0 : -30.0)});
    const System literal_conic = lens(rtt::model::Conic{Param(k == 0 ? 40.0 : -25.0), Param(-0.5)});
    const CompiledSystem ba = compiled(bound_axicon, k);
    const CompiledSystem la = compiled(literal_axicon, 0);
    const CompiledSystem bc = compiled(bound_conic, k);
    const CompiledSystem lc = compiled(literal_conic, 0);
    const auto& a = std::get<rtt::geom::Axicon<double>>(ba.surfaces()[0].shape);
    const auto& b = std::get<rtt::geom::Axicon<double>>(la.surfaces()[0].shape);
    CHECK(a.slope() == b.slope());
    CHECK(a.sag(1.25, -0.5) == b.sag(1.25, -0.5));
    const auto& c = std::get<rtt::geom::Conic<double>>(bc.surfaces()[0].shape);
    const auto& d = std::get<rtt::geom::Conic<double>>(lc.surfaces()[0].shape);
    CHECK(c.sag(1.25, -0.5) == d.sag(1.25, -0.5));
    CHECK(c.grad(1.25, -0.5) == d.grad(1.25, -0.5));
  }
}
