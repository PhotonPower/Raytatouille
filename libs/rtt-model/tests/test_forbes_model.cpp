// Q-con and Q-bfs in the model (#314, ADR 0043, points 2 and 4): validate of their parameters,
// the warning for an aperture beyond the normalization radius, the tolerance classes and kinds
// (ADR 0035, addendum #314), and the limit of curvature_relative on a curvature (ADR 0043,
// addendum PR C).

#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <limits>
#include <optional>
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

/// Pointer of the base shape of L1.S1 in make_singlet() (aperture radius 12.7 mm).
const std::string kBase = "/root/children/1/surfaces/0/shape/base";

System with_base(const BaseShape& base) {
  System s = make_singlet();
  element(s, 1).surfaces[0].shape.base = base;
  return s;
}

QCon qcon() {
  return {Param(1.0 / 50.0), Param(-0.5), Param(12.7), {Param(1e-3), Param(-2e-4)}};
}
QBfs qbfs() {
  return {Param(1.0 / 50.0), Param(-0.5), Param(12.7), {Param(1e-3), Param(-2e-4)}};
}

std::vector<Diagnostic> with_severity(const System& s, Severity severity) {
  std::vector<Diagnostic> out;
  for (const Diagnostic& d : validate(s)) {
    if (d.severity == severity) out.push_back(d);
  }
  return out;
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

ClassTolerance find(const ClassTolerances& r, const std::string& id) {
  for (const ClassTolerance& c : r.items) {
    if (c.id == id) return c;
  }
  FAIL("no generated item " << id);
  return {};
}

/// Both Forbes types with the same parameters.
std::vector<BaseShape> both(Param c, Param k, Param rho, std::vector<Param> a) {
  return {QCon{c, k, rho, a}, QBfs{c, k, rho, a}};
}

}  // namespace

TEST_CASE("Q-con and Q-bfs with valid parameters validate without diagnostics", "[model][forbes]") {
  for (const BaseShape& base :
       {BaseShape{qcon()}, BaseShape{qbfs()},
        BaseShape{QCon{Param(0.0), Param(0.0), Param(5.0), {Param(1.0)}}},
        BaseShape{QBfs{Param(0.0), Param(0.0), Param(5.0), {Param(1.0)}}}}) {
    // Aperture radius 5 mm, within every normalization radius here.
    System s = with_base(base);
    element(s, 1).surfaces[0].aperture = CircularAperture{5.0, 0.0};
    const std::vector<Diagnostic> d = validate(s);
    for (const Diagnostic& x : d) INFO(x.code << " at " << x.location << ": " << x.message);
    CHECK(d.empty());
  }
}

TEST_CASE("curvature, conic and coefficients of Q-con and Q-bfs are checked as for the odd asphere",
          "[model][forbes][validate]") {
  const double inf = std::numeric_limits<double>::infinity();
  const double nan = std::numeric_limits<double>::quiet_NaN();
  for (const double bad : {inf, nan}) {
    for (const BaseShape& base : both(Param(bad), Param(0.0), Param(12.7), {Param(1e-3)})) {
      CHECK(has(with_severity(with_base(base), Severity::Error), "shape.curvature_invalid",
                kBase + "/curvature"));
    }
  }
  for (const BaseShape& base : both(Param(0.0), Param(-1.0), Param(12.7), {Param(1e-3)})) {
    const System s = with_base(base);
    CHECK(with_severity(s, Severity::Error).empty());
    CHECK(has(with_severity(s, Severity::Warning), "shape.conic_without_curvature",
              kBase + "/conic"));
  }
  for (const BaseShape& base : both(Param(0.02), Param(0.0), Param(12.7), {})) {
    CHECK(has(with_severity(with_base(base), Severity::Warning),
              "shape.asphere_without_coefficients", kBase + "/coefficients"));
  }
  for (const BaseShape& base : both(Param(0.02), Param(nan), Param(12.7), {Param(1e-3)})) {
    CHECK(has(with_severity(with_base(base), Severity::Error), "value.not_finite",
              kBase + "/conic/value"));
  }
}

TEST_CASE("the normalization radius must be finite and positive", "[model][forbes][validate]") {
  // ADR 0043, point 4: shape.normalization_radius_invalid (error).
  const double inf = std::numeric_limits<double>::infinity();
  const double nan = std::numeric_limits<double>::quiet_NaN();
  for (const double bad : {0.0, -0.0, -3.0, inf, nan}) {
    INFO(bad);
    for (const BaseShape& base : both(Param(0.02), Param(0.0), Param(bad), {Param(1e-3)})) {
      const std::vector<Diagnostic> e = with_severity(with_base(base), Severity::Error);
      CHECK(has(e, "shape.normalization_radius_invalid", kBase + "/normalization_radius"));
      CHECK_FALSE(has_code(e, "value.not_finite"));
    }
  }
  for (const BaseShape& base : both(Param(0.02), Param(0.0), Param::bound("NOPE"), {Param(1)})) {
    CHECK(has(with_severity(with_base(base), Severity::Error), "param.unknown_parameter",
              kBase + "/normalization_radius/param"));
  }
}

TEST_CASE("Q-bfs needs phi > 0 at the normalization radius", "[model][forbes][validate]") {
  // ADR 0043, point 4: (1 + k) c^2 rho^2 >= 1 is shape.normalization_radius_invalid at Q-bfs
  // (1/sigma diverges at phi = 0), not at Q-con.
  // c = 1/10, rho = 10: c^2 rho^2 = 1.
  System s = with_base(QBfs{Param(0.1), Param(0.0), Param(10.0), {Param(1e-3)}});
  CHECK(has(with_severity(s, Severity::Error), "shape.normalization_radius_invalid",
            kBase + "/normalization_radius"));
  s = with_base(QBfs{Param(0.1), Param(-0.5), Param(10.0), {Param(1e-3)}});  // 0.5 < 1
  CHECK_FALSE(has_code(with_severity(s, Severity::Error), "shape.normalization_radius_invalid"));
  s = with_base(QCon{Param(0.1), Param(0.0), Param(10.0), {Param(1e-3)}});
  CHECK_FALSE(has_code(with_severity(s, Severity::Error), "shape.normalization_radius_invalid"));
}

TEST_CASE("an aperture beyond the normalization radius is a warning", "[model][forbes][validate]") {
  // ADR 0043, point 2 (F3): the largest distance of the aperture from the axis exceeds rho by
  // more than 1e-9 relative: radius, half the diagonal, the major semi-axis.
  struct Case {
    std::optional<Aperture> aperture;
    double rho;
    bool warned;
  };
  const double diagonal = std::hypot(3.0, 4.0);  // 5
  for (const Case& c :
       {Case{CircularAperture{12.7, 0.0}, 12.7, false},
        Case{CircularAperture{12.7, 0.0}, 12.7 * (1.0 - 1e-12), false},
        Case{CircularAperture{12.7, 0.0}, 12.6, true},
        Case{RectangularAperture{3.0, 4.0}, diagonal, false},
        Case{RectangularAperture{3.0, 4.0}, 4.5, true},
        Case{EllipticalAperture{2.0, 6.0}, 6.0, false},
        Case{EllipticalAperture{2.0, 6.0}, 5.0, true}, Case{std::nullopt, 1.0, false}}) {
    INFO("rho " << c.rho << ", warned " << c.warned);
    for (const BaseShape& base : both(Param(0.02), Param(0.0), Param(c.rho), {Param(1e-3)})) {
      System s = with_base(base);
      element(s, 1).surfaces[0].aperture = c.aperture;
      const std::vector<Diagnostic> w = with_severity(s, Severity::Warning);
      CHECK(has(w, "shape.normalization_radius_exceeded", kBase + "/normalization_radius") ==
            c.warned);
      CHECK(with_severity(s, Severity::Error).empty());
    }
  }
}

TEST_CASE("tolerance classes and kinds at Q-con and Q-bfs", "[model][forbes][tolerances]") {
  // ADR 0035, addendum #314: as at the odd asphere.
  for (const BaseShape& base : {BaseShape{qcon()}, BaseShape{qbfs()}}) {
    System s = with_base(base);
    s.tolerances.tolerance_class = ToleranceClass::Normal;
    ClassTolerances r = class_tolerances(s);
    CHECK(find(r, "class:surface:L1.S1:curvature_relative").skip == ClassSkip::None);
    CHECK(find(r, "class:surface:L1.S1:conic").skip == ClassSkip::None);
  }
  for (const BaseShape& base : both(Param(0.0), Param(0.0), Param(12.7), {Param(1e-3)})) {
    System s = with_base(base);
    s.tolerances.tolerance_class = ToleranceClass::Normal;
    const ClassTolerances r = class_tolerances(s);
    CHECK(find(r, "class:surface:L1.S1:curvature_relative").skip == ClassSkip::Plane);
    CHECK(find(r, "class:surface:L1.S1:conic").skip == ClassSkip::Plane);
  }
  struct Case {
    BaseShape base;
    std::string name;
  };
  for (const Case& c : {Case{qcon(), "qcon"}, Case{qbfs(), "qbfs"}}) {
    for (const K kind : {K::CurvatureRelative, K::Conic, K::Radius}) {
      INFO(c.name << " " << tolerance_kind_name(kind));
      System s = with_base(c.base);
      s.tolerances.items = {item("t", at(A::Surface, "L1.S1"), kind)};
      bool found = false;
      for (const Diagnostic& d : validate(s)) {
        if (d.code != "tolerances.target_invalid") continue;
        found = true;
        CHECK(d.message.find("base type \"" + c.name + "\"") != std::string::npos);
      }
      CHECK(found == (kind == K::Radius));
    }
    // curvature_relative and a param item on the curvature perturb the same quantity.
    System s = with_base(c.base);
    ToleranceItem p = item("p", at(A::Surface, "L1.S1"), K::Param);
    p.target.param = "/shape/base/curvature";
    s.tolerances.items = {item("c", at(A::Surface, "L1.S1"), K::CurvatureRelative), p};
    CHECK(has(with_severity(s, Severity::Error), "tolerances.quantity_duplicate",
              "/tolerances/items/1"));
  }
}

TEST_CASE("curvature_relative >= 1 is allowed on a curvature, not on a radius",
          "[model][forbes][tolerances]") {
  // ADR 0043, addendum PR C: on odd_asphere, qcon and qbfs C' = C (1 + delta) may reach 0 or
  // change its sign as a valid trial (ADR 0035, addendum #314), so a limit >= 1 is no error
  // there; on a radius (conic, even asphere) it stays tolerances.limit_invalid.
  for (const BaseShape& base : {BaseShape{qcon()}, BaseShape{qbfs()},
                                BaseShape{OddAsphere{Param(0.02), Param(0.0), {Param(1e-3)}}}}) {
    System s = with_base(base);
    ToleranceItem i = item("c", at(A::Surface, "L1.S1"), K::CurvatureRelative);
    i.limit = 1.5;
    s.tolerances.items = {i};
    CHECK_FALSE(has_code(with_severity(s, Severity::Error), "tolerances.limit_invalid"));
    // The single value of a class, too.
    System v = with_base(base);
    v.tolerances.tolerance_class = ToleranceClass::Normal;
    ToleranceValueOverride o = value(A::Surface, "L1.S1", K::CurvatureRelative);
    o.limit = 2.0;
    v.tolerances.overrides = {o};
    CHECK_FALSE(has_code(with_severity(v, Severity::Error), "tolerances.limit_invalid"));
  }
  for (const BaseShape& base : {BaseShape{Conic{Param(51.68), Param(0.0)}},
                                BaseShape{EvenAsphere{Param(51.68), Param(-0.5), {Param(1e-6)}}}}) {
    System s = with_base(base);
    ToleranceItem i = item("c", at(A::Surface, "L1.S1"), K::CurvatureRelative);
    i.limit = 1.5;
    s.tolerances.items = {i};
    CHECK(has(with_severity(s, Severity::Error), "tolerances.limit_invalid",
              "/tolerances/items/0/limit"));
    System v = with_base(base);
    v.tolerances.tolerance_class = ToleranceClass::Normal;
    ToleranceValueOverride o = value(A::Surface, "L1.S1", K::CurvatureRelative);
    o.limit = 2.0;
    v.tolerances.overrides = {o};
    CHECK(has(with_severity(v, Severity::Error), "tolerances.limit_invalid",
              "/tolerances/overrides/0/limit"));
  }
}
