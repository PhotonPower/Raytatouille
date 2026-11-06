// The anamorphic base shapes in the model (#315, ADR 0042, points 1 and 5): validate of their
// parameters, the class skips of the tolerance classes and the explicit tolerance kinds.

#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "rtt/model/model.hpp"
#include "rtt/model/tolerance_classes.hpp"
#include "test_support.hpp"
#include "tolerance_fixture.hpp"

using namespace rtt::model;
using namespace rtt::model::test;

namespace {

using K = ToleranceKind;
using A = ToleranceAnchor;

/// Pointer of the base shape of L1.S1 in make_singlet().
const std::string kBase = "/root/children/1/surfaces/0/shape/base";

System with_base(const BaseShape& base) {
  System s = make_singlet();
  element(s, 1).surfaces[0].shape.base = base;
  return s;
}

Biconic biconic() {
  return {Param(40.0), Param(-60.0), Param(-0.5), Param(0.25)};
}
Toroid toroid() {
  return {Param(48.0), Param(-1.2), {Param(2e-6)}, Param(55.0)};
}
Cylinder cylinder() {
  return {Param(-75.0), Param(0.0), {}};
}

std::vector<Diagnostic> errors(const System& s) {
  std::vector<Diagnostic> out;
  for (const Diagnostic& d : validate(s)) {
    if (d.severity == Severity::Error) out.push_back(d);
  }
  return out;
}

bool has(const std::vector<Diagnostic>& d, const std::string& code, const std::string& location) {
  for (const Diagnostic& x : d) {
    if (x.code == code && x.location == location) return true;
  }
  return false;
}

/// The class item with this id, by value (-Wdangling-reference, see test_tolerance_classes.cpp).
ClassTolerance find(const ClassTolerances& r, const std::string& id) {
  for (const ClassTolerance& c : r.items) {
    if (c.id == id) return c;
  }
  FAIL("no generated item " << id);
  return {};
}

}  // namespace

TEST_CASE("anamorphic shapes with valid parameters validate without errors",
          "[model][anamorphic]") {
  for (const BaseShape& base : {BaseShape{biconic()}, BaseShape{toroid()}, BaseShape{cylinder()},
                                BaseShape{Toroid{Param(48.0), Param(0.0), {}, Param(-20.0)}}}) {
    const std::vector<Diagnostic> d = errors(with_base(base));
    for (const Diagnostic& x : d) INFO(x.code << " at " << x.location << ": " << x.message);
    // A toroid or cylinder without coefficients is a plain toroid or cylinder, not an asphere
    // without coefficients.
    CHECK(d.empty());
  }
}

TEST_CASE("every radius of an anamorphic shape must be finite and non-zero",
          "[model][anamorphic][validate]") {
  const double inf = std::numeric_limits<double>::infinity();
  const double nan = std::numeric_limits<double>::quiet_NaN();
  for (const double bad : {0.0, -0.0, inf, -inf, nan}) {
    INFO(bad);
    Biconic b = biconic();
    b.radius_x = Param(bad);
    CHECK(has(errors(with_base(b)), "shape.radius_invalid", kBase + "/radius_x"));
    b = biconic();
    b.radius_y = Param(bad);
    CHECK(has(errors(with_base(b)), "shape.radius_invalid", kBase + "/radius_y"));
    Toroid t = toroid();
    t.radius = Param(bad);
    CHECK(has(errors(with_base(t)), "shape.radius_invalid", kBase + "/radius"));
    t = toroid();
    t.radius_of_rotation = Param(bad);
    CHECK(has(errors(with_base(t)), "shape.radius_invalid", kBase + "/radius_of_rotation"));
    Cylinder c = cylinder();
    c.radius = Param(bad);
    CHECK(has(errors(with_base(c)), "shape.radius_invalid", kBase + "/radius"));
  }
}

TEST_CASE("conics and coefficients of anamorphic shapes are checked as Params",
          "[model][anamorphic][validate]") {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  Biconic b = biconic();
  b.conic_y = Param(nan);
  CHECK(has(errors(with_base(b)), "value.not_finite", kBase + "/conic_y/value"));
  Toroid t = toroid();
  t.coefficients.push_back(Param(nan));
  CHECK(has(errors(with_base(t)), "value.not_finite", kBase + "/coefficients/1/value"));
  Cylinder c = cylinder();
  c.conic = Param(nan);
  CHECK(has(errors(with_base(c)), "value.not_finite", kBase + "/conic/value"));
  // A bound radius takes its value from the table: an unknown row is reported at the reference.
  b = biconic();
  b.radius_x = Param::bound("NOPE");
  CHECK(has(errors(with_base(b)), "param.unknown_parameter", kBase + "/radius_x/param"));
}

TEST_CASE("tolerance classes skip curvature and conic at anamorphic shapes with ShapeType",
          "[model][anamorphic][tolerances]") {
  // ADR 0042, point 5 (decision F2 of S0): the kinds are not defined for these types; the
  // reason is its own, not Plane.
  for (const BaseShape& base : {BaseShape{biconic()}, BaseShape{toroid()}, BaseShape{cylinder()}}) {
    System s = with_base(base);
    s.tolerances.tolerance_class = ToleranceClass::Normal;
    const ClassTolerances r = class_tolerances(s);
    CHECK(find(r, "class:surface:L1.S1:curvature_relative").skip == ClassSkip::ShapeType);
    CHECK(find(r, "class:surface:L1.S1:conic").skip == ClassSkip::ShapeType);
    // The other surfaces keep their reasons.
    CHECK(find(r, "class:surface:L1.S2:curvature_relative").skip == ClassSkip::Plane);
  }
  // Also when the radius is bound to a row: the kind is not defined at all.
  System s = with_base(Biconic{Param::bound("RX"), Param(-60.0), Param(0.0), Param(0.0)});
  s.parameters = {row("RX", 40.0)};
  s.tolerances.tolerance_class = ToleranceClass::Normal;
  const ClassTolerances r = class_tolerances(s);
  CHECK(find(r, "class:surface:L1.S1:curvature_relative").skip == ClassSkip::ShapeType);
}

TEST_CASE("a single value for a kind skipped by the shape type is reported as Skipped",
          "[model][anamorphic][tolerances]") {
  System s = with_base(toroid());
  s.tolerances.tolerance_class = ToleranceClass::Normal;
  s.tolerances.overrides = {value(A::Surface, "L1.S1", K::CurvatureRelative)};
  const ClassTolerances r = class_tolerances(s);
  REQUIRE(r.value_overrides.size() == 1);
  CHECK(r.value_overrides[0].effect == OverrideEffect::Skipped);
  CHECK(r.value_overrides[0].skip == ClassSkip::ShapeType);
}

TEST_CASE("explicit radius, curvature and conic items at anamorphic shapes are target_invalid",
          "[model][anamorphic][tolerances]") {
  // ADR 0042, point 4: tolerance kinds on R_x, R_y and the conic per axis come later (F5); the
  // message names the base type (ADR 0041, point 6).
  struct Case {
    BaseShape base;
    std::string name;
  };
  for (const Case& c :
       {Case{biconic(), "biconic"}, Case{toroid(), "toroid"}, Case{cylinder(), "cylinder"}}) {
    for (const K kind : {K::Radius, K::CurvatureRelative, K::Conic}) {
      INFO(c.name << " " << tolerance_kind_name(kind));
      System s = with_base(c.base);
      s.tolerances.items = {item("t", at(A::Surface, "L1.S1"), kind)};
      bool found = false;
      for (const Diagnostic& d : validate(s)) {
        if (d.code != "tolerances.target_invalid") continue;
        found = true;
        CHECK(d.location == "/tolerances/items/0/kind");
        CHECK(d.message.find("base type \"" + c.name + "\"") != std::string::npos);
      }
      CHECK(found);
    }
  }
}
