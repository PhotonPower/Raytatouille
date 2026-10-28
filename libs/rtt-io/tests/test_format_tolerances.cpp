// The section "tolerances" of schema 0.7.0 (ADR 0035; #253 part A): reader, canonical writer,
// edit form, read errors with their pointers and the migration 0.5 -> 0.6. The semantic checks
// (codes tolerances.*) are part B.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cstdint>
#include <functional>
#include <nlohmann/json.hpp>
#include <optional>
#include <string>
#include <utility>
#include <variant>

#include "rtt/io/edit.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/model/system.hpp"
#include "test_support.hpp"

using Catch::Matchers::ContainsSubstring;
using nlohmann::json;
using namespace rtt::model;

namespace {

/// Pointer and message of the ParseError for `j`, or "<no error>".
std::pair<std::string, std::string> error_of(const json& j) {
  try {
    [[maybe_unused]] const System s = rtt::io::parse_system(j.dump());
  } catch (const rtt::io::ParseError& e) {
    return {e.pointer(), e.what()};
  }
  return {"<no error>", ""};
}

ToleranceTarget at(ToleranceAnchor anchor, std::string name) {
  ToleranceTarget t;
  t.anchor = anchor;
  t.name = std::move(name);
  return t;
}

/// The singlet with a row "BFD" and a tolerance set that uses every form of ADR 0035 once, with
/// defaults and non-defaults mixed. Only the format is tested here: validate (part B) would
/// reject some of it (e.g. a segment at an element with one material).
System toleranced_singlet() {
  System s = test::make_singlet();
  ParameterRow bfd;
  bfd.name = "BFD";
  bfd.form = 100.0;
  s.parameters = {bfd};
  ToleranceSet& t = s.tolerances;
  t.tolerance_class = ToleranceClass::Normal;

  ToleranceValueOverride curvature{at(ToleranceAnchor::Surface, "L1.S1"),
                                   ToleranceKind::CurvatureRelative, 0.002,
                                   ToleranceDistribution::Uniform};
  ToleranceValueOverride index{at(ToleranceAnchor::Node, "L1"), ToleranceKind::Index, 0.0005,
                               ToleranceDistribution::Normal};
  index.target.segment = 0;
  t.overrides = {ToleranceClassOverride{"L1", ToleranceClass::Precise}, curvature, index,
                 TolerancePivotOverride{at(ToleranceAnchor::Node, "L1"), {0.0, 0.0, 2.0}},
                 ToleranceClassOverride{"image", ToleranceClass::None}};

  ToleranceItem dx{"L1-dx", at(ToleranceAnchor::Node, "L1"), ToleranceKind::DecenterX, 0.02,
                   ToleranceDistribution::Uniform};
  ToleranceItem row{"bfd", at(ToleranceAnchor::Parameter, "BFD"), ToleranceKind::Parameter, 0.05,
                    ToleranceDistribution::Normal};
  ToleranceItem conic{"conic", at(ToleranceAnchor::Surface, "L1.S1"), ToleranceKind::Param, 0.01,
                      ToleranceDistribution::Endpoint};
  conic.target.param = "/shape/base/conic";
  ToleranceItem tilt{"tilt", at(ToleranceAnchor::Surface, "L1.S2"), ToleranceKind::TiltYDeg, 0.05,
                     ToleranceDistribution::Normal};
  t.items = {dx, row, conic, tilt};

  Compensator focus{at(ToleranceAnchor::Node, "image"), -0.5, 0.5};
  focus.target.param = "/pose/position/2";
  t.compensators = {Compensator{at(ToleranceAnchor::Parameter, "BFD"), -1.0, 2.0}, focus};

  Optimization merit;
  SpotRmsOperand spot;
  spot.path = "main";
  merit.operands = {spot};
  t.compensator_merit = merit;

  SpotRmsOperand criterion_spot;
  criterion_spot.path = "main";
  criterion_spot.field = 1;
  FirstOrderOperand efl;
  efl.path = "main";
  RayOperand ray;
  ray.path = "main";
  ray.surface = SurfaceId("IMG");
  ray.py = 1.0;
  t.criteria = {Criterion{"spot", criterion_spot, std::nullopt, 0.01, std::nullopt},
                Criterion{"efl", efl, std::nullopt, std::nullopt, 0.5},
                Criterion{"y", ray, -0.1, 0.1, 0.02}};

  t.monte_carlo = MonteCarloDefaults{(std::uint64_t{1} << 53) - 1, 200};
  return s;
}

/// The canonical file of toleranced_singlet() as JSON, changed by `change`.
json changed(const std::function<void(json&)>& change) {
  json j = json::parse(rtt::io::to_json(toleranced_singlet()));
  change(j);
  return j;
}

}  // namespace

TEST_CASE("schema version is 0.7.0", "[io][tolerances]") {
  REQUIRE(kSchemaVersion == "0.7.0");
}

TEST_CASE("tolerances: the set round-trips in both forms", "[io][tolerances]") {
  const System s = toleranced_singlet();
  const std::string text = rtt::io::to_json(s);
  INFO(text);
  REQUIRE(rtt::io::parse_system(text) == s);
  REQUIRE(rtt::io::to_json(rtt::io::parse_system(text)) == text);
  const std::string edit = rtt::io::to_edit_json(s);
  REQUIRE(rtt::io::parse_system(edit) == s);
  REQUIRE(rtt::io::to_edit_json(rtt::io::parse_system(edit)) == edit);
  // A system without tolerances has an empty set.
  REQUIRE(rtt::io::parse_system(rtt::io::to_json(test::make_singlet())).tolerances.empty());
}

TEST_CASE("tolerances: canonical form", "[io][tolerances]") {
  const std::string text = rtt::io::to_json(toleranced_singlet());
  const json j = json::parse(text);
  REQUIRE(j.contains("tolerances"));  // guard: const access below
  // The section comes last, after optimization (ADR 0035, point 1).
  REQUIRE(text.find("\"paths\"") < text.find("\"tolerances\""));
  // Keys in the order of ADR 0035, point 1; defaults omitted (distribution normal).
  REQUIRE(j["tolerances"] == json::parse(R"({
    "class": "normal",
    "overrides": [
      {"node": "L1", "class": "precise"},
      {"surface": "L1.S1", "kind": "curvature_relative", "limit": 0.002,
       "distribution": "uniform"},
      {"node": "L1", "segment": 0, "kind": "index", "limit": 0.0005},
      {"node": "L1", "pivot": [0.0, 0.0, 2.0]},
      {"node": "image", "class": "none"}
    ],
    "items": [
      {"id": "L1-dx", "node": "L1", "kind": "decenter_x", "limit": 0.02,
       "distribution": "uniform"},
      {"id": "bfd", "parameter": "BFD", "kind": "parameter", "limit": 0.05},
      {"id": "conic", "surface": "L1.S1", "param": "/shape/base/conic", "kind": "param",
       "limit": 0.01, "distribution": "endpoint"},
      {"id": "tilt", "surface": "L1.S2", "kind": "tilt_y_deg", "limit": 0.05}
    ],
    "compensators": [
      {"parameter": "BFD", "travel_min": -1.0, "travel_max": 2.0},
      {"node": "image", "param": "/pose/position/2", "travel_min": -0.5, "travel_max": 0.5}
    ],
    "compensator_merit": {"operands": [{"type": "spot_rms", "path": "main", "target": 0.0}]},
    "criteria": [
      {"id": "spot", "operand": {"type": "spot_rms", "path": "main", "field": 1}, "max": 0.01},
      {"id": "efl", "operand": {"type": "efl", "path": "main"}, "max_change": 0.5},
      {"id": "y", "operand": {"type": "ray_y", "path": "main", "surface": "IMG", "py": 1.0},
       "min": -0.1, "max": 0.1, "max_change": 0.02}
    ],
    "monte_carlo": {"seed": 9007199254740991, "trials": 200}
  })"));
  // Without tolerances the section is missing; an empty set writes nothing either.
  REQUIRE_FALSE(json::parse(rtt::io::to_json(test::make_singlet())).contains("tolerances"));
  // Only a class: the lists are omitted.
  System only_class = test::make_singlet();
  only_class.tolerances.tolerance_class = ToleranceClass::Coarse;
  REQUIRE(json::parse(rtt::io::to_json(only_class))["tolerances"] == json{{"class", "coarse"}});
}

TEST_CASE("tolerances: \"class\": \"none\" reads as no class", "[io][tolerances]") {
  const System s = rtt::io::parse_system(
      changed([](json& j) { j["tolerances"] = json{{"class", "none"}}; }).dump());
  REQUIRE(s.tolerances.tolerance_class == ToleranceClass::None);
  REQUIRE(s.tolerances.empty());
  REQUIRE_FALSE(json::parse(rtt::io::to_json(s)).contains("tolerances"));
}

TEST_CASE("tolerances: a criterion from the file has target 0 and weight 1", "[io][tolerances]") {
  // ADR 0035, point 7: the operand of a criterion has neither key; the model keeps the defaults.
  const System s = rtt::io::parse_system(rtt::io::to_json(toleranced_singlet()));
  REQUIRE(s.tolerances.criteria.size() == 3);
  for (const Criterion& c : s.tolerances.criteria) {
    // By value: GCC 13 takes a reference to the result of std::visit for a dangling temporary.
    const OperandCommon common = std::visit([](const auto& op) { return op.common; }, c.operand);
    CHECK(common.target == 0.0);
    CHECK(common.weight == 1.0);
  }
}

TEST_CASE("tolerances: edit form (ADR 0024)", "[io][tolerances]") {
  SECTION("without tolerances: the class and every list, nothing optional") {
    const json e = json::parse(rtt::io::to_edit_json(test::make_singlet()));
    REQUIRE(e.contains("tolerances"));
    REQUIRE(e["tolerances"] == json::parse(R"({"class": "none", "overrides": [], "items": [],
                                              "compensators": [], "criteria": []})"));
  }
  SECTION("every value written, criteria without target and weight") {
    const json e = json::parse(rtt::io::to_edit_json(toleranced_singlet()));
    REQUIRE(e.contains("tolerances"));
    const json& t = e["tolerances"];
    CHECK(t["overrides"][1]["distribution"] == "uniform");
    CHECK(t["overrides"][2]["distribution"] == "normal");
    CHECK(t["items"][1]["distribution"] == "normal");
    CHECK(t["criteria"][0]["operand"] ==
          json::parse(R"({"type": "spot_rms", "path": "main", "field": 1,
                          "polychromatic": false, "reference": "centroid", "rings": 6})"));
    CHECK(t["compensator_merit"]["generators"] == json::array());
    CHECK(t["monte_carlo"] == json{{"seed", 9007199254740991ULL}, {"trials", 200}});
  }
  SECTION("a patch adds the first item and its inverse removes it") {
    const System s = test::make_singlet();
    const rtt::io::PatchResult r = rtt::io::apply_patch(s, R"([{"op": "add",
        "path": "/tolerances/items/-",
        "value": {"id": "gap", "node": "image", "kind": "thickness", "limit": 0.1}}])");
    REQUIRE(r.system.tolerances.items.size() == 1);
    CHECK(r.system.tolerances.items[0].kind == ToleranceKind::Thickness);
    CHECK(r.system.tolerances.items[0].distribution == ToleranceDistribution::Normal);
    const rtt::io::PatchResult back =
        rtt::io::apply_patch(r.system, r.inverse, rtt::io::EditCheck::StructureOnly);
    CHECK(back.system == s);
  }
}

TEST_CASE("tolerances: every kind, distribution and class name reads back", "[io][tolerances]") {
  const std::pair<const char*, ToleranceKind> kinds[] = {
      {"radius", ToleranceKind::Radius},
      {"curvature_relative", ToleranceKind::CurvatureRelative},
      {"conic", ToleranceKind::Conic},
      {"thickness", ToleranceKind::Thickness},
      {"param", ToleranceKind::Param},
      {"parameter", ToleranceKind::Parameter},
      {"decenter_x", ToleranceKind::DecenterX},
      {"decenter_y", ToleranceKind::DecenterY},
      {"shift_z", ToleranceKind::ShiftZ},
      {"tilt_x_deg", ToleranceKind::TiltXDeg},
      {"tilt_y_deg", ToleranceKind::TiltYDeg},
      {"tilt_z_deg", ToleranceKind::TiltZDeg},
      {"index", ToleranceKind::Index},
      {"abbe", ToleranceKind::Abbe}};
  for (const auto& entry : kinds) {
    const char* const name = entry.first;  // a local, not a structured binding, for the lambda
    const ToleranceKind kind = entry.second;
    INFO(name);
    const System s = rtt::io::parse_system(
        changed([&](json& j) { j["tolerances"]["items"][0]["kind"] = name; }).dump());
    CHECK(s.tolerances.items[0].kind == kind);
    CHECK(json::parse(rtt::io::to_json(s))["tolerances"]["items"][0]["kind"] == name);
  }
  const std::pair<const char*, ToleranceDistribution> distributions[] = {
      {"normal", ToleranceDistribution::Normal},
      {"uniform", ToleranceDistribution::Uniform},
      {"endpoint", ToleranceDistribution::Endpoint}};
  for (const auto& entry : distributions) {
    const char* const name = entry.first;
    const ToleranceDistribution d = entry.second;
    INFO(name);
    const System s = rtt::io::parse_system(
        changed([&](json& j) { j["tolerances"]["items"][0]["distribution"] = name; }).dump());
    CHECK(s.tolerances.items[0].distribution == d);
  }
  const std::pair<const char*, ToleranceClass> classes[] = {{"coarse", ToleranceClass::Coarse},
                                                            {"normal", ToleranceClass::Normal},
                                                            {"precise", ToleranceClass::Precise}};
  for (const auto& entry : classes) {
    const char* const name = entry.first;
    const ToleranceClass c = entry.second;
    INFO(name);
    const System s =
        rtt::io::parse_system(changed([&](json& j) { j["tolerances"]["class"] = name; }).dump());
    CHECK(s.tolerances.tolerance_class == c);
  }
}

TEST_CASE("tolerances: Monte Carlo seed and trials take the full uint64 range when read",
          "[io][tolerances]") {
  // The ranges of ADR 0035, point 8 are validate's (monte_carlo_invalid, part B).
  const System s = rtt::io::parse_system(changed([](json& j) {
                                           j["tolerances"]["monte_carlo"] = json{
                                               {"seed", 18446744073709551615ULL}, {"trials", 0}};
                                         }).dump());
  REQUIRE(s.tolerances.monte_carlo);
  CHECK(s.tolerances.monte_carlo->seed == 18446744073709551615ULL);
  CHECK(s.tolerances.monte_carlo->trials == 0U);
  const System empty = rtt::io::parse_system(
      changed([](json& j) { j["tolerances"]["monte_carlo"] = json::object(); }).dump());
  REQUIRE(empty.tolerances.monte_carlo);
  CHECK_FALSE(empty.tolerances.monte_carlo->seed);
  CHECK_FALSE(empty.tolerances.monte_carlo->trials);
}

TEST_CASE("tolerances: malformed entries are errors at their pointer", "[io][tolerances]") {
  // The reader checks keys and types (ADR 0008); whether kind, anchor, param and segment fit
  // together is validate's (tolerances.target_invalid, part B).
  using Change = std::function<void(json&)>;
  const std::pair<Change, std::string> cases[] = {
      // section
      {[](json& j) { j["tolerances"] = json::array(); }, "/tolerances"},
      {[](json& j) { j["tolerances"]["extra"] = 1; }, "/tolerances/extra"},
      {[](json& j) { j["tolerances"]["class"] = "fine"; }, "/tolerances/class"},
      {[](json& j) { j["tolerances"]["items"] = json::object(); }, "/tolerances/items"},
      // items: anchors
      {[](json& j) { j["tolerances"]["items"][0].erase("node"); }, "/tolerances/items/0"},
      {[](json& j) { j["tolerances"]["items"][0]["surface"] = "L1.S1"; }, "/tolerances/items/0"},
      {[](json& j) { j["tolerances"]["items"][0]["node"] = 3; }, "/tolerances/items/0/node"},
      {[](json& j) { j["tolerances"]["items"][1]["parameter"] = 5; },
       "/tolerances/items/1/parameter"},
      {[](json& j) { j["tolerances"]["items"][2]["param"] = 0; }, "/tolerances/items/2/param"},
      // items: other keys
      {[](json& j) { j["tolerances"]["items"][0].erase("id"); }, "/tolerances/items/0"},
      {[](json& j) { j["tolerances"]["items"][0].erase("limit"); }, "/tolerances/items/0"},
      {[](json& j) { j["tolerances"]["items"][0]["kind"] = "wedge"; }, "/tolerances/items/0/kind"},
      {[](json& j) { j["tolerances"]["items"][0]["limit"] = "0.1"; }, "/tolerances/items/0/limit"},
      {[](json& j) { j["tolerances"]["items"][0]["distribution"] = "gauss"; },
       "/tolerances/items/0/distribution"},
      {[](json& j) { j["tolerances"]["items"][0]["segment"] = -1; }, "/tolerances/items/0/segment"},
      {[](json& j) { j["tolerances"]["items"][0]["segment"] = 0.5; },
       "/tolerances/items/0/segment"},
      {[](json& j) { j["tolerances"]["items"][0]["pivot"] = json::array({0, 0, 1}); },
       "/tolerances/items/0/pivot"},
      // overrides: exactly one of the three forms
      {[](json& j) { j["tolerances"]["overrides"][0]["kind"] = "index"; },
       "/tolerances/overrides/0"},
      {[](json& j) { j["tolerances"]["overrides"][3]["kind"] = "tilt_x_deg"; },
       "/tolerances/overrides/3"},
      {[](json& j) { j["tolerances"]["overrides"][1]["pivot"] = json::array({0, 0, 0}); },
       "/tolerances/overrides/1"},
      {[](json& j) { j["tolerances"]["overrides"][0] = json{{"node", "L1"}}; },
       "/tolerances/overrides/0"},
      {[](json& j) {
         j["tolerances"]["overrides"][0] = json{{"surface", "L1.S1"}, {"class", "coarse"}};
       },
       "/tolerances/overrides/0/surface"},
      {[](json& j) { j["tolerances"]["overrides"][1]["parameter"] = "BFD"; },
       "/tolerances/overrides/1/parameter"},
      {[](json& j) { j["tolerances"]["overrides"][1]["param"] = "/shape/base/radius"; },
       "/tolerances/overrides/1/param"},
      {[](json& j) { j["tolerances"]["overrides"][3]["pivot"] = json::array({0, 0}); },
       "/tolerances/overrides/3/pivot"},
      {[](json& j) { j["tolerances"]["overrides"][3]["segment"] = 0; },
       "/tolerances/overrides/3/segment"},
      // compensators
      {[](json& j) { j["tolerances"]["compensators"][0].erase("travel_min"); },
       "/tolerances/compensators/0"},
      {[](json& j) { j["tolerances"]["compensators"][0]["segment"] = 0; },
       "/tolerances/compensators/0/segment"},
      {[](json& j) { j["tolerances"]["compensators"][0]["limit"] = 1; },
       "/tolerances/compensators/0/limit"},
      {[](json& j) { j["tolerances"]["compensators"][1]["travel_max"] = true; },
       "/tolerances/compensators/1/travel_max"},
      // compensator merit: as optimization
      {[](json& j) { j["tolerances"]["compensator_merit"]["operands"][0].erase("target"); },
       "/tolerances/compensator_merit/operands/0"},
      // criteria: an operand without target and weight, no generator
      {[](json& j) { j["tolerances"]["criteria"][0]["operand"]["target"] = 0.0; },
       "/tolerances/criteria/0/operand/target"},
      {[](json& j) { j["tolerances"]["criteria"][0]["operand"]["weight"] = 1.0; },
       "/tolerances/criteria/0/operand/weight"},
      {[](json& j) { j["tolerances"]["criteria"][0]["operand"]["type"] = "rms_spot"; },
       "/tolerances/criteria/0/operand/type"},
      {[](json& j) { j["tolerances"]["criteria"][0].erase("operand"); }, "/tolerances/criteria/0"},
      {[](json& j) { j["tolerances"]["criteria"][0].erase("id"); }, "/tolerances/criteria/0"},
      {[](json& j) { j["tolerances"]["criteria"][0]["max"] = "1"; }, "/tolerances/criteria/0/max"},
      {[](json& j) { j["tolerances"]["criteria"][0]["target"] = 0.0; },
       "/tolerances/criteria/0/target"},
      // Monte Carlo
      {[](json& j) { j["tolerances"]["monte_carlo"]["seed"] = -1; },
       "/tolerances/monte_carlo/seed"},
      {[](json& j) { j["tolerances"]["monte_carlo"]["seed"] = 1.5; },
       "/tolerances/monte_carlo/seed"},
      {[](json& j) { j["tolerances"]["monte_carlo"]["trials"] = "200"; },
       "/tolerances/monte_carlo/trials"},
      {[](json& j) { j["tolerances"]["monte_carlo"]["threads"] = 4; },
       "/tolerances/monte_carlo/threads"},
  };
  for (std::size_t i = 0; i < std::size(cases); ++i) {
    const auto& [change, pointer] = cases[i];
    INFO("case " << i << ": " << pointer);
    CHECK(error_of(changed(change)).first == pointer);
  }
}

TEST_CASE("migration 0.5 -> 0.6 changes only the version", "[io][tolerances][migration]") {
  json j = json::parse(rtt::io::to_json(test::make_singlet()));
  j["schema_version"] = "0.5.0";
  const System s = rtt::io::parse_system(j.dump());
  CHECK(s.schema_version == "0.7.0");
  CHECK(s.tolerances.empty());
  json written = json::parse(rtt::io::to_json(s));
  written["schema_version"] = "0.5.0";
  CHECK(written == j);
}

TEST_CASE("a tolerance section in an older file is an error", "[io][tolerances][migration]") {
  for (const char* version : {"0.5.0", "0.4.0", "0.2.0"}) {
    INFO(version);
    json j = json::parse(rtt::io::to_json(toleranced_singlet()));
    j["schema_version"] = version;
    const auto [pointer, message] = error_of(j);
    CHECK(pointer == "/tolerances");
    CHECK_THAT(message, ContainsSubstring("0.6"));
  }
}
