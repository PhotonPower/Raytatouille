// The odd asphere and the axicon in the model (#314, ADR 0043, point 4): validate of their
// parameters, the class skips of the tolerance classes and the explicit tolerance kinds (ADR
// 0035, addendum #314).

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

OddAsphere odd_asphere() {
  return {Param(1.0 / 50.0), Param(-0.5), {Param(1e-3), Param(0.0), Param(2e-6)}};
}
Axicon axicon() {
  return {Param(-5.0)};
}

std::vector<Diagnostic> with_severity(const System& s, Severity severity) {
  std::vector<Diagnostic> out;
  for (const Diagnostic& d : validate(s)) {
    if (d.severity == severity) out.push_back(d);
  }
  return out;
}

std::vector<Diagnostic> errors(const System& s) {
  return with_severity(s, Severity::Error);
}

std::vector<Diagnostic> warnings(const System& s) {
  return with_severity(s, Severity::Warning);
}

bool has(const std::vector<Diagnostic>& d, const std::string& code, const std::string& location) {
  for (const Diagnostic& x : d) {
    if (x.code == code && x.location == location) return true;
  }
  return false;
}

bool has_code(const std::vector<Diagnostic>& d, const std::string& code) {
  for (const Diagnostic& x : d) {
    if (x.code == code) return true;
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

TEST_CASE("odd asphere and axicon with valid parameters validate without errors",
          "[model][rotational]") {
  for (const BaseShape& base :
       {BaseShape{odd_asphere()}, BaseShape{axicon()}, BaseShape{Axicon{Param(0.0)}},
        BaseShape{Axicon{Param(89.9)}}, BaseShape{Axicon{Param(-89.9)}},
        BaseShape{OddAsphere{Param(0.0), Param(0.0), {Param(0.1)}}}}) {
    const System s = with_base(base);
    const std::vector<Diagnostic> d = validate(s);
    for (const Diagnostic& x : d) INFO(x.code << " at " << x.location << ": " << x.message);
    CHECK(d.empty());
  }
}

TEST_CASE("the curvature of an odd asphere must be finite; 0 is a plane base",
          "[model][rotational][validate]") {
  // ADR 0043, point 4: shape.curvature_invalid (error) for a curvature that is not finite.
  const double inf = std::numeric_limits<double>::infinity();
  const double nan = std::numeric_limits<double>::quiet_NaN();
  for (const double bad : {inf, -inf, nan}) {
    INFO(bad);
    OddAsphere a = odd_asphere();
    a.curvature = Param(bad);
    const std::vector<Diagnostic> e = errors(with_base(a));
    CHECK(has(e, "shape.curvature_invalid", kBase + "/curvature"));
    CHECK_FALSE(has_code(e, "value.not_finite"));  // one report per Param
  }
  for (const double zero : {0.0, -0.0}) {
    OddAsphere a = odd_asphere();
    a.curvature = Param(zero);
    a.conic = Param(0.0);
    CHECK(validate(with_base(a)).empty());
  }
  // A bound curvature takes its value from the table: an unknown row is reported at the
  // reference.
  OddAsphere a = odd_asphere();
  a.curvature = Param::bound("NOPE");
  CHECK(has(errors(with_base(a)), "param.unknown_parameter", kBase + "/curvature/param"));
}

TEST_CASE("a conic without curvature has no effect and is a warning",
          "[model][rotational][validate]") {
  // ADR 0043, point 4: shape.conic_without_curvature (warning) for conic != 0 at c = 0.
  OddAsphere a = odd_asphere();
  a.curvature = Param(0.0);
  a.conic = Param(-1.0);
  const System s = with_base(a);
  CHECK(errors(s).empty());
  CHECK(has(warnings(s), "shape.conic_without_curvature", kBase + "/conic"));
  // c != 0 with a conic, or c = 0 without one: no warning.
  CHECK_FALSE(has_code(warnings(with_base(odd_asphere())), "shape.conic_without_curvature"));
  a.conic = Param(0.0);
  CHECK_FALSE(has_code(warnings(with_base(a)), "shape.conic_without_curvature"));
}

TEST_CASE("an odd asphere without coefficients is an asphere without coefficients",
          "[model][rotational][validate]") {
  // ADR 0043, point 4: the warning also holds for an empty list of the odd asphere.
  const System s = with_base(OddAsphere{Param(0.02), Param(0.0), {}});
  CHECK(errors(s).empty());
  CHECK(has(warnings(s), "shape.asphere_without_coefficients", kBase + "/coefficients"));
}

TEST_CASE("conic and coefficients of an odd asphere are checked as Params",
          "[model][rotational][validate]") {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  OddAsphere a = odd_asphere();
  a.conic = Param(nan);
  CHECK(has(errors(with_base(a)), "value.not_finite", kBase + "/conic/value"));
  a = odd_asphere();
  a.coefficients[1] = Param(nan);
  CHECK(has(errors(with_base(a)), "value.not_finite", kBase + "/coefficients/1/value"));
}

TEST_CASE("the axicon angle must be finite and below 90 deg in magnitude",
          "[model][rotational][validate]") {
  // ADR 0043, point 4 (F5): shape.axicon_angle_invalid (error).
  const double inf = std::numeric_limits<double>::infinity();
  const double nan = std::numeric_limits<double>::quiet_NaN();
  for (const double bad : {90.0, -90.0, 120.0, -180.0, inf, -inf, nan}) {
    INFO(bad);
    const std::vector<Diagnostic> e = errors(with_base(Axicon{Param(bad)}));
    CHECK(has(e, "shape.axicon_angle_invalid", kBase + "/angle_deg"));
    CHECK_FALSE(has_code(e, "value.not_finite"));
  }
  CHECK(has(errors(with_base(Axicon{Param::bound("NOPE")})), "param.unknown_parameter",
            kBase + "/angle_deg/param"));
}

TEST_CASE("tolerance classes generate curvature_relative and conic at an odd asphere",
          "[model][rotational][tolerances]") {
  // ADR 0035, addendum #314: as for even_asphere, conic always (no "sphere" case).
  System s = with_base(OddAsphere{Param(0.02), Param(0.0), {Param(1e-3)}});
  s.tolerances.tolerance_class = ToleranceClass::Normal;
  ClassTolerances r = class_tolerances(s);
  CHECK(find(r, "class:surface:L1.S1:curvature_relative").skip == ClassSkip::None);
  CHECK(find(r, "class:surface:L1.S1:conic").skip == ClassSkip::None);
  // c = 0: a plane base, both skipped with Plane.
  s = with_base(OddAsphere{Param(0.0), Param(-1.0), {Param(1e-3)}});
  s.tolerances.tolerance_class = ToleranceClass::Normal;
  r = class_tolerances(s);
  CHECK(find(r, "class:surface:L1.S1:curvature_relative").skip == ClassSkip::Plane);
  CHECK(find(r, "class:surface:L1.S1:conic").skip == ClassSkip::Plane);
  // Bound Params are skipped with Bound.
  s = with_base(OddAsphere{Param::bound("C"), Param::bound("K"), {Param(1e-3)}});
  s.parameters = {row("C", 0.02), row("K", -0.5)};
  s.tolerances.tolerance_class = ToleranceClass::Normal;
  r = class_tolerances(s);
  CHECK(find(r, "class:surface:L1.S1:curvature_relative").skip == ClassSkip::Bound);
  CHECK(find(r, "class:surface:L1.S1:conic").skip == ClassSkip::Bound);
}

TEST_CASE("tolerance classes skip curvature and conic at an axicon with ShapeType",
          "[model][rotational][tolerances]") {
  System s = with_base(axicon());
  s.tolerances.tolerance_class = ToleranceClass::Normal;
  const ClassTolerances r = class_tolerances(s);
  CHECK(find(r, "class:surface:L1.S1:curvature_relative").skip == ClassSkip::ShapeType);
  CHECK(find(r, "class:surface:L1.S1:conic").skip == ClassSkip::ShapeType);
}

TEST_CASE("explicit tolerance kinds at an odd asphere and an axicon",
          "[model][rotational][tolerances]") {
  // ADR 0035, addendum #314: curvature_relative and conic at the odd asphere; radius there is
  // target_invalid, as every kind of shape at the axicon, with the base type in the message.
  struct Case {
    BaseShape base;
    std::string name;
    K kind;
    bool valid;
  };
  for (const Case& c :
       {Case{odd_asphere(), "odd_asphere", K::CurvatureRelative, true},
        Case{odd_asphere(), "odd_asphere", K::Conic, true},
        Case{odd_asphere(), "odd_asphere", K::Radius, false},
        Case{axicon(), "axicon", K::CurvatureRelative, false},
        Case{axicon(), "axicon", K::Conic, false}, Case{axicon(), "axicon", K::Radius, false}}) {
    INFO(c.name << " " << tolerance_kind_name(c.kind));
    System s = with_base(c.base);
    s.tolerances.items = {item("t", at(A::Surface, "L1.S1"), c.kind)};
    bool found = false;
    for (const Diagnostic& d : validate(s)) {
      if (d.code != "tolerances.target_invalid") continue;
      found = true;
      CHECK(d.location == "/tolerances/items/0/kind");
      CHECK(d.message.find("base type \"" + c.name + "\"") != std::string::npos);
    }
    CHECK(found == !c.valid);
    if (c.valid) CHECK(errors(s).empty());
  }
}

TEST_CASE("curvature_relative perturbs the curvature of an odd asphere",
          "[model][rotational][tolerances]") {
  // The perturbed quantity is /shape/base/curvature (ADR 0035, addendum #314): a param item on
  // it and a curvature_relative item at the same surface perturb the same quantity.
  System s = with_base(odd_asphere());
  ToleranceItem p = item("p", at(A::Surface, "L1.S1"), K::Param);
  p.target.param = "/shape/base/curvature";
  s.tolerances.items = {item("c", at(A::Surface, "L1.S1"), K::CurvatureRelative), p};
  CHECK(has(errors(s), "tolerances.quantity_duplicate", "/tolerances/items/1"));
  // conic and a param item on the conic, too.
  ToleranceItem q = item("q", at(A::Surface, "L1.S1"), K::Param);
  q.target.param = "/shape/base/conic";
  s.tolerances.items = {item("k", at(A::Surface, "L1.S1"), K::Conic), q};
  CHECK(has(errors(s), "tolerances.quantity_duplicate", "/tolerances/items/1"));
  // A bound curvature is tolerated through its row.
  s = with_base(OddAsphere{Param::bound("C"), Param(0.0), {Param(1e-3)}});
  s.parameters = {row("C", 0.02)};
  s.tolerances.items = {item("c", at(A::Surface, "L1.S1"), K::CurvatureRelative)};
  CHECK(has(errors(s), "tolerances.bound_param", "/tolerances/items/0"));
  // The angle of the axicon goes through the kind param.
  s = with_base(axicon());
  ToleranceItem angle = item("a", at(A::Surface, "L1.S1"), K::Param);
  angle.target.param = "/shape/base/angle_deg";
  s.tolerances.items = {angle};
  CHECK(errors(s).empty());
}
