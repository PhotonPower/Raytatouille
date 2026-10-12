// validate() on the tolerance set (ADR 0035, point 9, and its addendum; #253 part B): every code
// tolerances.* at its pointer, the merit codes in criteria and compensator merit, and a valid set
// without diagnostics.

#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "rtt/model/model.hpp"
#include "test_support.hpp"
#include "tolerance_fixture.hpp"

using namespace rtt::model;
using namespace rtt::model::test;

namespace {

using K = ToleranceKind;
using A = ToleranceAnchor;

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

ToleranceTarget with_param(ToleranceTarget t, std::string pointer) {
  t.param = std::move(pointer);
  return t;
}

ToleranceTarget with_segment(ToleranceTarget t, std::uint32_t segment) {
  t.segment = segment;
  return t;
}

SpotRmsOperand spot(std::uint16_t field = 0) {
  SpotRmsOperand o;
  o.path = "main";
  o.field = field;
  return o;
}

/// The fixture with a valid set that uses every part of the section: no tolerances.* or merit.*
/// diagnostic.
System valid() {
  System s = toleranced();
  ParameterRow expr;
  expr.name = "HALF";
  expr.form = ParameterExpression{"BFD / 2"};
  s.parameters.push_back(expr);
  ToleranceSet& t = s.tolerances;
  t.overrides = {ToleranceClassOverride{"objective", ToleranceClass::Precise},
                 ToleranceValueOverride{with_segment(at(A::Node, "L1"), 1), K::Index, 0.001,
                                        ToleranceDistribution::Uniform},
                 TolerancePivotOverride{at(A::Node, "objective"), {0.0, 0.0, 9.5}},
                 ToleranceClassOverride{"image", ToleranceClass::None}};
  t.items = {
      item("S3-conic", with_param(at(A::Surface, "L1.S3"), "/shape/base/conic"), K::Param),
      item("S1-radius", at(A::Surface, "L1.S1"), K::Radius),
      item("S2-curv", at(A::Surface, "L1.S2"), K::CurvatureRelative),
      item("L1-glass", with_segment(at(A::Node, "L1"), 0), K::Abbe),
      item("L2-index", at(A::Node, "L2"), K::Index),
      item("asph", with_param(at(A::Surface, "L2.S1"), "/shape/base/coefficients/0"), K::Param),
      item("bfd", at(A::Parameter, "BFD"), K::Parameter),
      item("roll", at(A::Node, "fold"), K::TiltZDeg),
      item("shift", at(A::Surface, "FOLD"), K::ShiftZ)};
  t.compensators = {Compensator{at(A::Parameter, "BFD"), -1.0, 1.0},
                    Compensator{with_param(at(A::Node, "fold"), "/pose/position/2"), -0.5, 0.5}};
  Optimization merit;
  merit.operands = {spot()};
  t.compensator_merit = merit;
  FirstOrderOperand efl;
  efl.path = "main";
  t.criteria = {Criterion{"spot", spot(1), std::nullopt, 0.03, std::nullopt},
                Criterion{"efl", efl, std::nullopt, std::nullopt, 0.5}};
  t.monte_carlo = MonteCarloDefaults{(std::uint64_t{1} << 53) - 1, std::uint64_t{1} << 32};
  return s;
}

/// The diagnostics of validate() whose code starts with tolerances. or merit.
std::vector<Diagnostic> tolerance_diagnostics(const System& s) {
  std::vector<Diagnostic> out;
  for (const Diagnostic& d : validate(s)) {
    if (d.code.starts_with("tolerances.") || d.code.starts_with("merit.")) out.push_back(d);
  }
  return out;
}

bool reports(const System& s, const std::string& code, const std::string& location) {
  for (const Diagnostic& d : validate(s)) {
    if (d.code == code && d.location == location) return true;
  }
  return false;
}

struct Case {
  std::string name;
  std::function<void(System&)> change;
  std::string code;
  std::string location;
};

}  // namespace

TEST_CASE("tolerances: a valid set has no diagnostics", "[model][tolerances][validate]") {
  const std::vector<Diagnostic> d = tolerance_diagnostics(valid());
  for (const Diagnostic& x : d) UNSCOPED_INFO(to_string(x));
  CHECK(d.empty());
  // The fixture without tolerances neither.
  CHECK(tolerance_diagnostics(toleranced()).empty());
}

TEST_CASE("tolerances: every code at its pointer", "[model][tolerances][validate]") {
  const std::string items = "/tolerances/items/";
  const std::string over = "/tolerances/overrides/";
  const std::string comp = "/tolerances/compensators/";
  const std::string crit = "/tolerances/criteria/";
  const std::vector<Case> cases = {
      // unknown_target
      {"unknown node", [](System& s) { s.tolerances.items[0].target = at(A::Node, "L9"); },
       "tolerances.unknown_target", items + "0/node"},
      {"unknown surface", [](System& s) { s.tolerances.items[1].target = at(A::Surface, "L9.S1"); },
       "tolerances.unknown_target", items + "1/surface"},
      {"unknown row", [](System& s) { s.tolerances.items[6].target = at(A::Parameter, "X"); },
       "tolerances.unknown_target", items + "6/parameter"},
      {"unknown class node",
       [](System& s) { std::get<ToleranceClassOverride>(s.tolerances.overrides[0]).node = "L9"; },
       "tolerances.unknown_target", over + "0/node"},
      {"unknown pivot surface",
       [](System& s) {
         std::get<TolerancePivotOverride>(s.tolerances.overrides[2]).target = at(A::Surface, "X");
       },
       "tolerances.unknown_target", over + "2/surface"},
      {"unknown compensator row",
       [](System& s) { s.tolerances.compensators[0].target = at(A::Parameter, "X"); },
       "tolerances.unknown_target", comp + "0/parameter"},
      // target_invalid
      {"radius at a plane",
       [](System& s) { s.tolerances.items[1].target = at(A::Surface, "L2.S2"); },
       "tolerances.target_invalid", items + "1/kind"},
      {"radius at a node", [](System& s) { s.tolerances.items[1].target = at(A::Node, "L1"); },
       "tolerances.target_invalid", items + "1/kind"},
      {"index at an assembly",
       [](System& s) { s.tolerances.items[4].target = at(A::Node, "objective"); },
       "tolerances.target_invalid", items + "4/kind"},
      {"rigid at a row", [](System& s) { s.tolerances.items[7].target = at(A::Parameter, "BFD"); },
       "tolerances.target_invalid", items + "7/kind"},
      {"the root", [](System& s) { s.tolerances.items[7].target = at(A::Node, "system"); },
       "tolerances.target_invalid", items + "7/node"},
      {"param missing", [](System& s) { s.tolerances.items[0].target.param.reset(); },
       "tolerances.target_invalid", items + "0"},
      {"param not a Param",
       [](System& s) { s.tolerances.items[0].target.param = "/shape/base/conic/value"; },
       "tolerances.target_invalid", items + "0/param"},
      {"param at another kind",
       [](System& s) { s.tolerances.items[7].target.param = "/pose/position/0"; },
       "tolerances.target_invalid", items + "7/param"},
      {"segment missing", [](System& s) { s.tolerances.items[3].target.segment.reset(); },
       "tolerances.target_invalid", items + "3"},
      {"segment out of range", [](System& s) { s.tolerances.items[3].target.segment = 2; },
       "tolerances.target_invalid", items + "3/segment"},
      {"segment at one material", [](System& s) { s.tolerances.items[4].target.segment = 0; },
       "tolerances.target_invalid", items + "4/segment"},
      {"segment at another kind", [](System& s) { s.tolerances.items[7].target.segment = 0; },
       "tolerances.target_invalid", items + "7/segment"},
      {"two pivots",
       [](System& s) {
         s.tolerances.overrides.push_back(
             TolerancePivotOverride{at(A::Node, "objective"), {0.0, 0.0, 1.0}});
       },
       "tolerances.target_invalid", over + "4"},
      {"single value of a kind no class generates there",
       [](System& s) {
         s.tolerances.overrides.push_back(ToleranceValueOverride{at(A::Node, "L1"), K::ShiftZ, 0.1,
                                                                 ToleranceDistribution::Normal});
       },
       "tolerances.target_invalid", over + "4/kind"},
      // bound_param
      {"thickness at a bound position",
       [](System& s) { s.tolerances.items[7] = item("gap", at(A::Node, "image"), K::Thickness); },
       "tolerances.bound_param", items + "7"},
      {"param at a bound Param",
       [](System& s) {
         s.tolerances.items[0].target = with_param(at(A::Surface, "L2.S2"), "/pose/position/2");
       },
       "tolerances.bound_param", items + "0/param"},
      {"expression row", [](System& s) { s.tolerances.items[6].target = at(A::Parameter, "HALF"); },
       "tolerances.bound_param", items + "6/parameter"},
      {"compensator at an expression row",
       [](System& s) { s.tolerances.compensators[0].target = at(A::Parameter, "HALF"); },
       "tolerances.bound_param", comp + "0/parameter"},
      // limit_invalid
      {"limit zero", [](System& s) { s.tolerances.items[7].limit = 0.0; },
       "tolerances.limit_invalid", items + "7/limit"},
      {"limit NaN", [](System& s) { s.tolerances.items[7].limit = kNaN; },
       "tolerances.limit_invalid", items + "7/limit"},
      {"curvature_relative 1", [](System& s) { s.tolerances.items[2].limit = 1.0; },
       "tolerances.limit_invalid", items + "2/limit"},
      {"radius |R|", [](System& s) { s.tolerances.items[1].limit = 61.5; },
       "tolerances.limit_invalid", items + "1/limit"},
      {"single value negative",
       [](System& s) { std::get<ToleranceValueOverride>(s.tolerances.overrides[1]).limit = -1.0; },
       "tolerances.limit_invalid", over + "1/limit"},
      {"pivot NaN",
       [](System& s) {
         std::get<TolerancePivotOverride>(s.tolerances.overrides[2]).pivot[1] = kNaN;
       },
       "tolerances.limit_invalid", over + "2/pivot"},
      // ids
      {"empty id", [](System& s) { s.tolerances.items[7].id = ""; }, "tolerances.id_invalid",
       items + "7/id"},
      {"class: id", [](System& s) { s.tolerances.items[7].id = "class:node:fold:tilt_z_deg"; },
       "tolerances.id_invalid", items + "7/id"},
      {"empty criterion id", [](System& s) { s.tolerances.criteria[1].id = " "; },
       "tolerances.id_invalid", crit + "1/id"},
      {"duplicate id", [](System& s) { s.tolerances.items[8].id = "roll"; },
       "tolerances.id_duplicate", items + "8/id"},
      {"duplicate criterion id", [](System& s) { s.tolerances.criteria[1].id = "spot"; },
       "tolerances.id_duplicate", crit + "1/id"},
      {"hash collision",
       [](System& s) {
         // FNV-1a 64 of both is 0x4EAC0C95540867E4 (computed with RFC 9923, Sec. 2 and 5).
         s.tolerances.items[7].id = "8yn0iYCKYHlIj4-BwPqk";
         s.tolerances.items[8].id = "GReLUrM4wMqfg9yzV3KQ";
       },
       "tolerances.id_hash_collision", items + "8/id"},
      // quantity_duplicate
      {"radius and curvature_relative",
       [](System& s) { s.tolerances.items[2].target = at(A::Surface, "L1.S1"); },
       "tolerances.quantity_duplicate", items + "2"},
      {"the same Param through another anchor",
       [](System& s) {
         s.tolerances.items.push_back(item(
             "via-node", with_param(at(A::Node, "L1"), "/surfaces/2/shape/base/conic"), K::Param));
       },
       "tolerances.quantity_duplicate", items + "9"},
      // class_unknown
      {"class out of range",
       [](System& s) { s.tolerances.tolerance_class = static_cast<ToleranceClass>(7); },
       "tolerances.class_unknown", "/tolerances/class"},
      // materials
      {"abbe at CONST", [](System& s) { s.tolerances.items[4].kind = K::Abbe; },
       "tolerances.abbe_on_const", items + "4/kind"},
      {"index at a metal", [](System& s) { s.tolerances.items[4].target = at(A::Node, "fold"); },
       "tolerances.material_unsupported", items + "4/node"},
      // compensators
      {"compensator without param",
       [](System& s) { s.tolerances.compensators[1].target.param.reset(); },
       "tolerances.compensator_target", comp + "1"},
      {"compensator at a values row",
       [](System& s) {
         s.configurations = {{"a"}, {"b"}};
         s.parameters[0].form = std::vector<double>{4.0, 5.0};
         s.tolerances.compensators[0].target = at(A::Parameter, "T2");
       },
       "tolerances.compensator_target", comp + "0/parameter"},
      {"two compensators on one Param",
       [](System& s) { s.tolerances.compensators.push_back(s.tolerances.compensators[1]); },
       "tolerances.compensator_target", comp + "2"},
      {"travel_min >= 0", [](System& s) { s.tolerances.compensators[0].travel_min = 0.0; },
       "tolerances.compensator_bounds", comp + "0/travel_min"},
      {"travel_max NaN", [](System& s) { s.tolerances.compensators[0].travel_max = kNaN; },
       "tolerances.compensator_bounds", comp + "0/travel_max"},
      {"no compensator merit and no optimization",
       [](System& s) { s.tolerances.compensator_merit.reset(); },
       "tolerances.compensator_merit_missing", "/tolerances/compensators"},
      {"an empty compensator merit (addendum)",
       [](System& s) { s.tolerances.compensator_merit = Optimization{}; },
       "tolerances.compensator_merit_missing", "/tolerances/compensator_merit"},
      // criteria
      {"criterion without bound", [](System& s) { s.tolerances.criteria[0].max.reset(); },
       "tolerances.criterion_bounds", crit + "0"},
      {"min > max", [](System& s) { s.tolerances.criteria[0].min = 0.05; },
       "tolerances.criterion_bounds", crit + "0/min"},
      {"max_change 0", [](System& s) { s.tolerances.criteria[1].max_change = 0.0; },
       "tolerances.criterion_bounds", crit + "1/max_change"},
      {"max NaN", [](System& s) { s.tolerances.criteria[0].max = kNaN; },
       "tolerances.criterion_bounds", crit + "0/max"},
      {"criterion target",
       [](System& s) {
         std::get<SpotRmsOperand>(s.tolerances.criteria[0].operand).common.target = 1.0;
       },
       "tolerances.criterion_bounds", crit + "0/operand/target"},
      {"criterion weight",
       [](System& s) {
         std::get<SpotRmsOperand>(s.tolerances.criteria[0].operand).common.weight = 2.0;
       },
       "tolerances.criterion_bounds", crit + "0/operand/weight"},
      {"merit code in a criterion",
       [](System& s) { std::get<SpotRmsOperand>(s.tolerances.criteria[0].operand).path = "x"; },
       "merit.unknown_path", crit + "0/operand/path"},
      {"merit code in the compensator merit",
       [](System& s) {
         std::get<SpotRmsOperand>(s.tolerances.compensator_merit->operands[0]).field = 9;
       },
       "merit.index_out_of_range", "/tolerances/compensator_merit/operands/0/field"},
      // Monte Carlo
      {"seed 2^53", [](System& s) { s.tolerances.monte_carlo->seed = std::uint64_t{1} << 53; },
       "tolerances.monte_carlo_invalid", "/tolerances/monte_carlo/seed"},
      {"trials 0", [](System& s) { s.tolerances.monte_carlo->trials = 0; },
       "tolerances.monte_carlo_invalid", "/tolerances/monte_carlo/trials"},
      {"trials 2^32 + 1",
       [](System& s) { s.tolerances.monte_carlo->trials = (std::uint64_t{1} << 32) + 1; },
       "tolerances.monte_carlo_invalid", "/tolerances/monte_carlo/trials"},
  };
  for (const Case& c : cases) {
    INFO(c.name);
    System s = valid();
    c.change(s);
    const bool found = reports(s, c.code, c.location);
    if (!found) {
      for (const Diagnostic& d : tolerance_diagnostics(s)) UNSCOPED_INFO(to_string(d));
    }
    CHECK(found);
  }
}

TEST_CASE("tolerances: the fallback on optimization and sets without compensators",
          "[model][tolerances][validate]") {
  System s = valid();
  s.tolerances.compensator_merit.reset();
  SECTION("a non-empty optimization is the compensator merit") {
    s.optimization.operands = {spot()};
    CHECK(tolerance_diagnostics(s).empty());
  }
  SECTION("without compensators an empty merit is fine") {
    s.tolerances.compensators.clear();
    CHECK(tolerance_diagnostics(s).empty());
    s.tolerances.compensator_merit = Optimization{};
    CHECK(tolerance_diagnostics(s).empty());
  }
}

TEST_CASE("tolerances: a class alone never invalidates a valid file",
          "[model][tolerances][validate]") {
  // ADR 0035, point 5: skips and single values without effect are not errors.
  for (const ToleranceClass c :
       {ToleranceClass::Coarse, ToleranceClass::Normal, ToleranceClass::Precise}) {
    System s = valid();
    s.tolerances.tolerance_class = c;
    s.tolerances.overrides.push_back(ToleranceValueOverride{
        at(A::Surface, "L1.S1"), K::Conic, 0.1, ToleranceDistribution::Normal});  // a sphere
    s.tolerances.overrides.push_back(ToleranceValueOverride{
        at(A::Node, "image"), K::Thickness, 0.1, ToleranceDistribution::Normal});  // class none
    CHECK(tolerance_diagnostics(s).empty());
  }
}
