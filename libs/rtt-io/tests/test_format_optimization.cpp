// The section "optimization" of schema 0.4.0 (ADR 0030, points 2-4; #162 part B): reader,
// canonical writer, edit form and read errors with their pointers.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cstdint>
#include <functional>
#include <nlohmann/json.hpp>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "rtt/io/edit.hpp"
#include "rtt/io/json_io.hpp"
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

OperandCommon target(double t) {
  OperandCommon c;
  c.target = t;
  return c;
}

/// The singlet (path "main", 2 fields, 1 wavelength) with configurations, a row "D" and a merit
/// function with every type: defaults and non-defaults mixed.
System merit_singlet() {
  System s = test::make_singlet();
  s.configurations = {{"near"}, {"far"}};
  ParameterRow d;
  d.name = "D";
  d.form = 5.0;
  s.parameters = {d};
  FirstOrderOperand efl{target(100.0), FirstOrderQuantity::Efl, "main", std::nullopt};
  efl.common.configuration = "far";
  FirstOrderOperand fno{target(4.0), FirstOrderQuantity::ImageFNumber, "main", 0};
  fno.common.weight = 0.5;
  RayOperand ray{target(-0.25), RayCoordinate::X, "main", SurfaceId("IMG"), 0, 1, 0.5,
                 -1.0,          std::nullopt};
  SpotRmsOperand spot{target(0.0), "main", 1, std::nullopt, true, SpotReference::Chief, 4};
  SpotRmsOperand plain_spot{
      target(0.0), "main", 0, std::nullopt, false, SpotReference::Centroid, 6};
  OpdRmsOperand opd{target(0.0), "main", 0, 0, 17};
  ParamValueOperand value{target(5.0), "D"};
  value.common.configuration = "near";
  s.optimization.operands = {efl, fno, ray, spot, plain_spot, opd, value};
  SpotGenerator spots;
  spots.path = "main";
  spots.fields = std::vector<std::uint16_t>{0, 1};
  spots.wavelengths = std::vector<std::uint16_t>{0};
  spots.reference = SpotReference::Chief;
  spots.rings = 4;
  spots.arms = 8;
  spots.weight = 2.0;
  WavefrontGenerator waves;
  waves.path = "main";
  waves.configuration = "near";
  s.optimization.generators = {spots, waves};
  return s;
}

/// The canonical file of merit_singlet() as JSON, changed by `change`.
json changed(const std::function<void(json&)>& change) {
  json j = json::parse(rtt::io::to_json(merit_singlet()));
  change(j);
  return j;
}

}  // namespace

TEST_CASE("optimization: the merit function round-trips in both forms", "[io][optimization]") {
  const System s = merit_singlet();
  const std::string text = rtt::io::to_json(s);
  INFO(text);
  REQUIRE(rtt::io::parse_system(text) == s);
  REQUIRE(rtt::io::to_json(rtt::io::parse_system(text)) == text);
  REQUIRE(rtt::io::parse_system(rtt::io::to_edit_json(s)) == s);
}

TEST_CASE("optimization: canonical form", "[io][optimization]") {
  const std::string text = rtt::io::to_json(merit_singlet());
  const json j = json::parse(text);
  // The section comes last, after the paths it refers to.
  REQUIRE(text.find("\"paths\"") < text.find("\"optimization\""));
  // Defaults are omitted; the target is always written (ADR 0030, point 2).
  REQUIRE(j["optimization"]["operands"] == json::parse(R"([
    {"type": "efl", "path": "main", "configuration": "far", "target": 100.0},
    {"type": "image_fnumber", "path": "main", "wavelength": 0, "target": 4.0, "weight": 0.5},
    {"type": "ray_x", "path": "main", "surface": "IMG", "occurrence": 0, "field": 1, "px": 0.5,
     "py": -1.0, "target": -0.25},
    {"type": "spot_rms", "path": "main", "field": 1, "polychromatic": true, "reference": "chief",
     "rings": 4, "target": 0.0},
    {"type": "spot_rms", "path": "main", "target": 0.0},
    {"type": "opd_rms", "path": "main", "wavelength": 0, "grid": 17, "target": 0.0},
    {"type": "param_value", "parameter": "D", "configuration": "near", "target": 5.0}
  ])"));
  REQUIRE(j["optimization"]["generators"] == json::parse(R"([
    {"type": "rms_spot", "path": "main", "fields": [0, 1], "wavelengths": [0],
     "reference": "chief", "rings": 4, "arms": 8, "weight": 2.0},
    {"type": "rms_wavefront", "path": "main", "configuration": "near"}
  ])"));
  // Without a merit function the section is missing; with only operands there is no generator
  // list.
  REQUIRE_FALSE(json::parse(rtt::io::to_json(test::make_singlet())).contains("optimization"));
  System only_operands = merit_singlet();
  only_operands.optimization.generators.clear();
  const json o = json::parse(rtt::io::to_json(only_operands))["optimization"];
  REQUIRE(o.contains("operands"));
  REQUIRE_FALSE(o.contains("generators"));
}

TEST_CASE("optimization: edit form writes every value (ADR 0024)", "[io][optimization]") {
  const json j = json::parse(rtt::io::to_edit_json(merit_singlet()));
  const json& ops = j["optimization"]["operands"];
  REQUIRE(ops[0] == json::parse(R"({"type": "efl", "path": "main", "configuration": "far",
                                     "target": 100.0, "weight": 1.0})"));
  REQUIRE(ops[4] == json::parse(R"({"type": "spot_rms", "path": "main", "field": 0,
                                     "polychromatic": false, "reference": "centroid", "rings": 6,
                                     "target": 0.0, "weight": 1.0})"));
  REQUIRE(ops[2] == json::parse(R"({"type": "ray_x", "path": "main", "surface": "IMG",
                                     "occurrence": 0, "field": 1, "px": 0.5, "py": -1.0,
                                     "target": -0.25, "weight": 1.0})"));
  REQUIRE(j["optimization"]["generators"][1] ==
          json::parse(R"({"type": "rms_wavefront", "path": "main", "configuration": "near",
                          "rings": 3, "arms": 6, "weight": 1.0})"));
  // The edit form writes the section also without a merit function, so that a patch can add the
  // first operand with "add /optimization/operands/-".
  const json plain = json::parse(rtt::io::to_edit_json(test::make_singlet()));
  REQUIRE(plain["optimization"] == json::parse(R"({"operands": [], "generators": []})"));
  REQUIRE(rtt::io::parse_system(plain.dump()) == test::make_singlet());
}

TEST_CASE("optimization: read errors point at the offending value", "[io][optimization]") {
  const auto at = [](const std::function<void(json&)>& change) {
    return error_of(changed(change));
  };
  const std::string op = "/optimization/operands";
  const std::string gen = "/optimization/generators";
  CHECK(at([](json& j) { j["optimization"]["operands"][0]["type"] = "focus"; }).first ==
        op + "/0/type");
  CHECK(at([](json& j) { j["optimization"]["operands"][0].erase("target"); }).first == op + "/0");
  CHECK(at([](json& j) { j["optimization"]["operands"][0].erase("path"); }).first == op + "/0");
  CHECK(at([](json& j) { j["optimization"]["operands"][6]["path"] = "main"; }).first ==
        op + "/6/path");
  CHECK(at([](json& j) { j["optimization"]["operands"][0]["grid"] = 5; }).first == op + "/0/grid");
  CHECK(at([](json& j) { j["optimization"]["generators"][0]["target"] = 0.0; }).first ==
        gen + "/0/target");
  CHECK(at([](json& j) { j["optimization"]["generators"][1]["reference"] = "chief"; }).first ==
        gen + "/1/reference");
  CHECK(at([](json& j) { j["optimization"]["generators"][0]["type"] = "rms_ghost"; }).first ==
        gen + "/0/type");
  // A polychromatic spot has no wavelength.
  const auto poly = at([](json& j) { j["optimization"]["operands"][3]["wavelength"] = 0; });
  CHECK(poly.first == op + "/3/wavelength");
  CHECK_THAT(poly.second, ContainsSubstring("polychromatic"));
  // Selections are not empty; missing means all.
  CHECK(at([](json& j) { j["optimization"]["generators"][0]["fields"] = json::array(); }).first ==
        gen + "/0/fields");
  // Indices are integers 0 ... 65535.
  CHECK(at([](json& j) { j["optimization"]["operands"][2]["field"] = -1; }).first ==
        op + "/2/field");
  CHECK(at([](json& j) { j["optimization"]["operands"][2]["field"] = 65536; }).first ==
        op + "/2/field");
  CHECK(at([](json& j) { j["optimization"]["operands"][1]["wavelength"] = 0.5; }).first ==
        op + "/1/wavelength");
  CHECK(at([](json& j) { j["optimization"]["generators"][0]["wavelengths"][0] = 70000; }).first ==
        gen + "/0/wavelengths/0");
  CHECK(at([](json& j) { j["optimization"]["operands"][2]["occurrence"] = -1; }).first ==
        op + "/2/occurrence");
  CHECK(at([](json& j) { j["optimization"]["operands"][3]["rings"] = 2.5; }).first ==
        op + "/3/rings");
  CHECK(at([](json& j) { j["optimization"]["operands"][3]["reference"] = "best"; }).first ==
        op + "/3/reference");
  CHECK(at([](json& j) { j["optimization"]["unknown"] = 1; }).first == "/optimization/unknown");
  CHECK(at([](json& j) { j["optimization"]["operands"] = json::object(); }).first ==
        "/optimization/operands");
  // The values themselves (weight < 0, an unknown path, ...) are validate's (merit.*).
  CHECK(at([](json& j) { j["optimization"]["operands"][0]["weight"] = -1.0; }).first ==
        "<no error>");
}

TEST_CASE("optimization: a merit function needs schema 0.4", "[io][optimization]") {
  json j = changed([](json&) {});
  j["schema_version"] = "0.3.0";
  j.erase("configurations");
  j.erase("parameters");
  const auto [pointer, message] = error_of(j);
  CHECK(pointer == "/optimization");
  CHECK_THAT(message, ContainsSubstring("0.4"));
}
