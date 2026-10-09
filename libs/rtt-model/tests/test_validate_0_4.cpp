// validate() for the model of schema 0.4.0 (#162): relative placement (ADR 0028, point 5) and the
// parameter table without evaluating expressions (ADR 0029, point 5; the evaluating codes come
// with #164).

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "rtt/model/validate.hpp"
#include "test_support.hpp"

using namespace rtt::model;
using rtt::model::test::element;
using rtt::model::test::make_singlet;

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

/// True if `d` has exactly one diagnostic with this code, and it is at `location`.
bool only_at(const std::vector<Diagnostic>& d, const std::string& code, const std::string& location) {
  const auto n = std::count_if(d.begin(), d.end(), [&](const Diagnostic& x) { return x.code == code; });
  return n == 1 && std::any_of(d.begin(), d.end(), [&](const Diagnostic& x) {
           return x.code == code && x.location == location;
         });
}

bool none_with(const std::vector<Diagnostic>& d, const std::string& code) {
  return std::none_of(d.begin(), d.end(), [&](const Diagnostic& x) { return x.code == code; });
}

ParameterRow row(std::string name, ParameterForm form) {
  ParameterRow r;
  r.name = std::move(name);
  r.form = std::move(form);
  return r;
}

}  // namespace

// ------------------------------------------------------------ relative placement -----

TEST_CASE("relative placement: valid references raise nothing (ADR 0028)", "[validate][0.4]") {
  System s = make_singlet();
  // L1 relative to the stop surface, L1.S2 relative to its sibling L1.S1, the detector relative
  // to its sibling L1, all with a rotate_first order somewhere.
  element(s, 1).pose.reference = PoseReference::RelativeToPreceding;
  element(s, 1).surfaces[1].pose.reference = PoseReference::RelativeToSibling;
  element(s, 2).pose.reference = PoseReference::RelativeToSibling;
  element(s, 2).pose.order = PoseOrder::RotateFirst;
  REQUIRE(validate(s).empty());
}

TEST_CASE("relative placement: the first surface of an element is absolute", "[validate][0.4]") {
  for (const PoseReference r : {PoseReference::RelativeToPreceding, PoseReference::RelativeToSibling}) {
    System s = make_singlet();
    element(s, 1).surfaces[0].pose.reference = r;
    const auto d = validate(s);
    INFO(static_cast<int>(r));
    REQUIRE(only_at(d, "pose.relative_first_surface", "/root/children/1/surfaces/0/pose/reference"));
    // One diagnostic per location (ADR 0028, point 5).
    REQUIRE(none_with(d, "pose.no_sibling"));
    REQUIRE(none_with(d, "pose.no_preceding"));
  }
}

TEST_CASE("relative placement: no preceding surface, no sibling", "[validate][0.4]") {
  System s = make_singlet();
  element(s, 0).pose.reference = PoseReference::RelativeToPreceding;
  REQUIRE(only_at(validate(s), "pose.no_preceding", "/root/children/0/pose/reference"));

  s = make_singlet();
  s.root.pose.reference = PoseReference::RelativeToPreceding;
  REQUIRE(only_at(validate(s), "pose.no_preceding", "/root/pose/reference"));

  s = make_singlet();
  element(s, 0).pose.reference = PoseReference::RelativeToSibling;
  REQUIRE(only_at(validate(s), "pose.no_sibling", "/root/children/0/pose/reference"));

  s = make_singlet();
  s.root.pose.reference = PoseReference::RelativeToSibling;
  REQUIRE(only_at(validate(s), "pose.no_sibling", "/root/pose/reference"));
}

TEST_CASE("relative placement: the preceding surface may lie in an earlier assembly",
          "[validate][0.4]") {
  // Pre-order (ADR 0028, point 2): the first element of a nested assembly refers to the last
  // surface before the assembly.
  System s = make_singlet();
  Assembly group;
  group.name = "group";
  Element lens = element(s, 1);
  lens.name = "L2";
  lens.surfaces[0].id = SurfaceId("L2.S1");
  lens.surfaces[1].id = SurfaceId("L2.S2");
  lens.pose.reference = PoseReference::RelativeToPreceding;
  group.children.push_back({lens});
  s.root.children.insert(s.root.children.begin() + 2, {group});
  REQUIRE(validate(s).empty());
}

// --------------------------------------------------------------- parameter table -----

TEST_CASE("parameter table: a valid table with configurations raises nothing", "[validate][0.4]") {
  System s = make_singlet();
  s.configurations = {{"wide"}, {"tele"}};
  ParameterRow g = row("G", std::vector<double>{20.0, 5.0});
  g.variable = true;
  g.min = 2.0;
  g.max = 30.0;
  s.parameters = {row("TOTAL", 40.0), g, row("B", ParameterExpression{"TOTAL - G"})};
  element(s, 2).pose.position[2] = Param::bound("B");
  REQUIRE(validate(s).empty());
}

TEST_CASE("parameter table: names", "[validate][0.4]") {
  for (const char* bad : {"", "1A", "A-B", "A B", "\xC3\x84"}) {  // last: a non-ASCII letter
    INFO(bad);
    System s = make_singlet();
    s.parameters = {row(bad, 1.0)};
    REQUIRE(only_at(validate(s), "parameters.name_invalid", "/parameters/0/name"));
  }
  System s = make_singlet();
  s.parameters = {row("_a1", 1.0), row("B", 2.0), row("_a1", 3.0)};
  REQUIRE(only_at(validate(s), "parameters.name_duplicate", "/parameters/2/name"));
  s.parameters = {row("a", 1.0), row("A", 2.0)};  // case matters
  REQUIRE(validate(s).empty());
}

TEST_CASE("parameter table: one value per configuration", "[validate][0.4]") {
  System s = make_singlet();
  s.parameters = {row("G", std::vector<double>{1.0, 2.0})};  // one configuration
  REQUIRE(only_at(validate(s), "parameters.values_count", "/parameters/0/values"));
  s.configurations = {{"a"}, {"b"}};
  REQUIRE(validate(s).empty());
  s.configurations = {{"a"}, {"b"}, {"c"}};
  REQUIRE(only_at(validate(s), "parameters.values_count", "/parameters/0/values"));
}

TEST_CASE("parameter table: a derived row is not variable and has no bounds", "[validate][0.4]") {
  System s = make_singlet();
  s.parameters = {row("A", 1.0), row("B", ParameterExpression{"2 * A"})};
  s.parameters[1].variable = true;
  REQUIRE(only_at(validate(s), "parameters.variable_expression", "/parameters/1/variable"));
  s.parameters[1].variable = false;
  s.parameters[1].min = 0.0;  // only over the API: a file cannot carry it (ADR 0029, point 1)
  REQUIRE(only_at(validate(s), "bounds.invalid", "/parameters/1/min"));
}

TEST_CASE("configuration names", "[validate][0.4]") {
  System s = make_singlet();
  s.configurations = {{"wide"}, {"  "}};
  REQUIRE(only_at(validate(s), "configurations.name_invalid", "/configurations/1/name"));
  s.configurations = {{""}, {"tele"}};
  REQUIRE(only_at(validate(s), "configurations.name_invalid", "/configurations/0/name"));
  s.configurations = {{"wide"}, {"tele"}, {"wide"}};
  REQUIRE(only_at(validate(s), "configurations.name_duplicate", "/configurations/2/name"));
}

// ----------------------------------------------------------- Param binding, bounds -----

TEST_CASE("a bound Param names an existing row and carries nothing else", "[validate][0.4]") {
  System s = make_singlet();
  s.parameters = {row("D", 106.363)};
  element(s, 2).pose.position[2] = Param::bound("D");
  REQUIRE(validate(s).empty());

  element(s, 2).pose.position[2] = Param::bound("X");
  REQUIRE(only_at(validate(s), "param.unknown_parameter", "/root/children/2/pose/position/2/param"));

  Param conflict = Param::bound("D");
  conflict.variable = true;
  element(s, 2).pose.position[2] = conflict;
  REQUIRE(only_at(validate(s), "param.bound_conflict", "/root/children/2/pose/position/2"));
  conflict.variable = false;
  conflict.max = 1.0;
  element(s, 2).pose.position[2] = conflict;
  REQUIRE(only_at(validate(s), "param.bound_conflict", "/root/children/2/pose/position/2"));
}

TEST_CASE("the value of a bound Param is not read (ADR 0029, point 3)", "[validate][0.4]") {
  // Its value is 0 until resolve_parameters; checks of the value must skip it.
  System s = make_singlet();
  s.parameters = {row("R", 51.68), row("EPD", 20.0), row("DIST", 100.0)};
  std::get<Conic>(element(s, 1).surfaces[0].shape.base).radius = Param::bound("R");
  s.aperture.value = Param::bound("EPD");
  s.object.at_infinity = false;
  s.object.distance = Param::bound("DIST");
  REQUIRE(validate(s).empty());
}

TEST_CASE("bounds of Params and rows", "[validate][0.4]") {
  System s = make_singlet();
  Param& conic = std::get<Conic>(element(s, 1).surfaces[0].shape.base).conic;
  conic.min = 1.0;
  conic.max = 1.0;
  REQUIRE(only_at(validate(s), "bounds.invalid", "/root/children/1/surfaces/0/shape/base/conic/min"));
  conic.min = -1.0;
  conic.max = -0.5;  // value 0 lies above max: a warning only
  const auto d = validate(s);
  REQUIRE(only_at(d, "bounds.value_outside", "/root/children/1/surfaces/0/shape/base/conic/value"));
  REQUIRE_FALSE(has_errors(d));

  s = make_singlet();
  s.configurations = {{"a"}, {"b"}};
  ParameterRow g = row("G", std::vector<double>{20.0, 5.0});
  g.min = 10.0;
  s.parameters = {g};
  REQUIRE(only_at(validate(s), "bounds.value_outside", "/parameters/0/values/1"));
  s.parameters[0] = row("H", 3.0);
  s.parameters[0].max = 2.0;
  REQUIRE(only_at(validate(s), "bounds.value_outside", "/parameters/0/value"));
  s.parameters[0].min = 5.0;
  REQUIRE(only_at(validate(s), "bounds.invalid", "/parameters/0/min"));
}

TEST_CASE("value.not_finite for numbers without their own check", "[validate][0.4]") {
  System s = make_singlet();
  std::get<Conic>(element(s, 1).surfaces[0].shape.base).conic = Param(kNaN);
  REQUIRE(only_at(validate(s), "value.not_finite", "/root/children/1/surfaces/0/shape/base/conic/value"));

  s = make_singlet();
  element(s, 1).pose.pivot[1] = kNaN;
  REQUIRE(only_at(validate(s), "value.not_finite", "/root/children/1/pose/pivot/1"));

  s = make_singlet();
  element(s, 1).pose.rotation_deg[0] = Param(std::numeric_limits<double>::infinity());
  REQUIRE(only_at(validate(s), "value.not_finite", "/root/children/1/pose/rotation_deg/0/value"));

  s = make_singlet();
  s.parameters = {row("A", kNaN)};
  REQUIRE(only_at(validate(s), "value.not_finite", "/parameters/0/value"));

  s = make_singlet();
  Param& bounded = std::get<Conic>(element(s, 1).surfaces[0].shape.base).conic;
  bounded.max = kNaN;
  REQUIRE(only_at(validate(s), "value.not_finite", "/root/children/1/surfaces/0/shape/base/conic/max"));

  // A number with its own check keeps it: a NaN radius stays shape.radius_invalid.
  s = make_singlet();
  std::get<Conic>(element(s, 1).surfaces[0].shape.base).radius = Param(kNaN);
  const auto d = validate(s);
  REQUIRE(none_with(d, "value.not_finite"));
  REQUIRE(only_at(d, "shape.radius_invalid", "/root/children/1/surfaces/0/shape/base/radius"));
}
