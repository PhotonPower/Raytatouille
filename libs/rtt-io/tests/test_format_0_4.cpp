// Schema 0.4 (#162; ADR 0028: Pose.reference and Pose.order; ADR 0029: parameter table,
// configurations, bound Params with bounds; pickup dropped with the warning io.pickup_dropped).

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <filesystem>
#include <fstream>
#include <functional>
#include <limits>
#include <nlohmann/json.hpp>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "rtt/io/edit.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/model/validate.hpp"
#include "test_support.hpp"

using Catch::Matchers::ContainsSubstring;
using nlohmann::json;
using rtt::model::Diagnostic;
using rtt::model::Element;
using rtt::model::Param;
using rtt::model::ParameterExpression;
using rtt::model::ParameterRow;
using rtt::model::PoseOrder;
using rtt::model::PoseReference;
using rtt::model::Severity;
using rtt::model::System;
using rtt::model::test::element;

namespace {

/// Pointer of the ParseError for `j`, or "<no error>".
std::string error_pointer(const json& j) {
  try {
    [[maybe_unused]] const System s = rtt::io::parse_system(j.dump());
  } catch (const rtt::io::ParseError& e) {
    return e.pointer();
  }
  return "<no error>";
}

std::string error_message(const json& j) {
  try {
    [[maybe_unused]] const System s = rtt::io::parse_system(j.dump());
  } catch (const rtt::io::ParseError& e) {
    return e.what();
  }
  return "<no error>";
}

/// The canonical file of `s` as JSON, changed by `change`.
json changed(const System& s, const std::function<void(json&)>& change) {
  json j = json::parse(rtt::io::to_json(s));
  change(j);
  return j;
}

const Element& element_of(const System& s, std::size_t i) {
  return std::get<Element>(s.root.children[i].value);
}

ParameterRow row(std::string name, rtt::model::ParameterForm form) {
  ParameterRow r;
  r.name = std::move(name);
  r.form = std::move(form);
  return r;
}

/// Zoom of ADR 0029, point 8, on the singlet: two configurations, three rows, the detector
/// relative to its sibling at the distance B, L1 relative to the stop with rotate_first, and a
/// variable conic with bounds.
System zoom_singlet() {
  System s = rtt::model::test::make_singlet();
  s.configurations = {{"wide"}, {"tele"}};
  ParameterRow g = row("G", std::vector<double>{20.0, 5.0});
  g.variable = true;
  g.min = 2.0;
  g.max = 30.0;
  s.parameters = {row("TOTAL", 40.0), g, row("B", ParameterExpression{"TOTAL - G"})};
  element(s, 1).pose.reference = PoseReference::RelativeToPreceding;
  element(s, 1).pose.order = PoseOrder::RotateFirst;
  element(s, 2).pose.reference = PoseReference::RelativeToSibling;
  element(s, 2).pose.position[2] = Param::bound("B");
  Param conic(-0.5);
  conic.variable = true;
  conic.min = -1.0;
  conic.max = 0.0;
  std::get<rtt::model::Conic>(element(s, 1).surfaces[0].shape.base).conic = conic;
  return s;
}

/// A 0.3 file: stop and detector; `detector_pose` and `extra` (top-level members) are inserted.
json file_0_3(const json& detector_pose, const json& extra = json::object()) {
  json j = json::parse(R"({
    "schema_version": "0.3.0",
    "units": {"length": "mm", "wavelength": "um"},
    "wavelengths": [{"um": 0.55, "reference": true}],
    "aperture": {"type": "epd", "value": 10.0},
    "fields": {"points": [{}]},
    "root": {"type": "assembly", "name": "r", "children": [
      {"type": "stop", "name": "stop",
       "surfaces": [{"id": "STO", "aperture": {"type": "circular", "radius": 5.0}}]},
      {"type": "detector", "name": "image", "surfaces": [{"id": "IMG"}]}
    ]},
    "paths": [{"name": "main", "events": "auto"}]
  })");
  if (!detector_pose.empty()) j["root"]["children"][1]["pose"] = detector_pose;
  for (const auto& item : extra.items()) j[item.key()] = item.value();
  return j;
}

}  // namespace

TEST_CASE("schema version is 0.4.0", "[io][0.4]") {
  REQUIRE(rtt::model::kSchemaVersion == "0.4.0");
}

TEST_CASE("relative poses, bounds, bound Params and the table round-trip", "[io][0.4]") {
  const System s = zoom_singlet();
  const std::string text = rtt::io::to_json(s);
  INFO(text);
  REQUIRE(rtt::io::parse_system(text) == s);
  REQUIRE(rtt::io::to_json(rtt::io::parse_system(text)) == text);
  // The edit form reads back to the same model as well (ADR 0024).
  REQUIRE(rtt::io::parse_system(rtt::io::to_edit_json(s)) == s);
}

TEST_CASE("canonical form of the 0.4 fields", "[io][0.4]") {
  const std::string text = rtt::io::to_json(zoom_singlet());
  const json j = json::parse(text);
  // Both sections before root, configurations first (ADR 0029, point 1).
  REQUIRE(text.find("\"configurations\"") < text.find("\"parameters\""));
  REQUIRE(text.find("\"parameters\"") < text.find("\"root\""));
  REQUIRE(j["configurations"] == json::parse(R"([{"name": "wide"}, {"name": "tele"}])"));
  REQUIRE(j["parameters"] == json::parse(R"([
    {"name": "TOTAL", "value": 40.0},
    {"name": "G", "values": [20.0, 5.0], "variable": true, "min": 2.0, "max": 30.0},
    {"name": "B", "expression": "TOTAL - G"}
  ])"));
  const json& children = j["root"]["children"];
  REQUIRE(children[1]["pose"] == json::parse(R"({"reference": "relative_to_preceding",
    "order": "rotate_first", "position": [0.0, 0.0, 5.0]})"));
  REQUIRE(children[2]["pose"] == json::parse(R"({"reference": "relative_to_sibling",
    "position": [0.0, 0.0, {"param": "B"}]})"));
  REQUIRE(children[1]["surfaces"][0]["shape"]["base"]["conic"] ==
          json::parse(R"({"value": -0.5, "variable": true, "min": -1.0, "max": 0.0})"));
  // Defaults are omitted: an absolute translate_first pose has neither key.
  REQUIRE_FALSE(children[0].contains("pose"));
  REQUIRE_FALSE(children[1]["surfaces"][1]["pose"].contains("reference"));
  REQUIRE_FALSE(children[1]["surfaces"][1]["pose"].contains("order"));
  // Without a table or configurations both sections are missing.
  const json plain = json::parse(rtt::io::to_json(rtt::model::test::make_singlet()));
  REQUIRE_FALSE(plain.contains("configurations"));
  REQUIRE_FALSE(plain.contains("parameters"));
}

TEST_CASE("edit form of the 0.4 fields (ADR 0024, ADR 0029 point 3)", "[io][0.4]") {
  const json j = json::parse(rtt::io::to_edit_json(zoom_singlet()));
  // Every row carries variable, also a derived one; bounds only where set.
  REQUIRE(j["parameters"][0] ==
          json::parse(R"({"name": "TOTAL", "value": 40.0, "variable": false})"));
  REQUIRE(j["parameters"][2] ==
          json::parse(R"({"name": "B", "expression": "TOTAL - G", "variable": false})"));
  // A bound Param is only the reference, without a value.
  REQUIRE(j["root"]["children"][2]["pose"]["position"][2] == json::parse(R"({"param": "B"})"));
  // An unbound Param: value and variable, bounds if set; poses with reference and order.
  REQUIRE(j["root"]["children"][2]["pose"]["position"][0] ==
          json::parse(R"({"value": 0.0, "variable": false})"));
  REQUIRE(j["root"]["children"][0]["pose"]["reference"] == "absolute");
  REQUIRE(j["root"]["children"][0]["pose"]["order"] == "translate_first");
  // The edit form writes the table also when empty; configurations only if there are any,
  // because an empty section is not readable (ADR 0029, point 1).
  const json plain = json::parse(rtt::io::to_edit_json(rtt::model::test::make_singlet()));
  REQUIRE_FALSE(plain.contains("configurations"));
  REQUIRE(plain["parameters"] == json::array());
  REQUIRE(rtt::io::parse_system(plain.dump()) == rtt::model::test::make_singlet());
}

TEST_CASE("malformed 0.4 values are errors at their pointer", "[io][0.4]") {
  const System s = zoom_singlet();
  const auto bound_d = [](json& j) -> json& {
    return j["root"]["children"][2]["pose"]["position"][2];
  };
  // A bound Param carries nothing but the name (ADR 0029, point 3).
  CHECK(error_pointer(changed(s, [&](json& j) { bound_d(j)["value"] = 1.0; })) ==
        "/root/children/2/pose/position/2/value");
  CHECK(error_pointer(changed(s, [&](json& j) { bound_d(j)["param"] = 3; })) ==
        "/root/children/2/pose/position/2/param");
  // pickup is gone in 0.4.
  CHECK(error_pointer(changed(s, [&](json& j) {
          bound_d(j) = json::parse(R"({"value": 1.0, "pickup": "2 * 3"})");
        })) == "/root/children/2/pose/position/2/pickup");
  CHECK(error_pointer(changed(s, [](json& j) {
          j["root"]["children"][2]["pose"]["reference"] = "sibling";
        })) == "/root/children/2/pose/reference");
  CHECK(error_pointer(changed(s, [](json& j) {
          j["root"]["children"][1]["pose"]["order"] = "zyx";
        })) == "/root/children/1/pose/order");
  // A row has exactly one of value, values, expression; bounds only at independent rows.
  CHECK(error_pointer(changed(s, [](json& j) { j["parameters"][0].erase("value"); })) ==
        "/parameters/0");
  CHECK(error_pointer(changed(s, [](json& j) { j["parameters"][0]["expression"] = "1"; })) ==
        "/parameters/0");
  CHECK(error_pointer(changed(s, [](json& j) { j["parameters"][2]["min"] = 0.0; })) ==
        "/parameters/2/min");
  CHECK(error_pointer(changed(s, [](json& j) { j["parameters"][1]["values"] = 20.0; })) ==
        "/parameters/1/values");
  CHECK(error_pointer(changed(s, [](json& j) { j["parameters"][0]["name"] = 1; })) ==
        "/parameters/0/name");
  // Configurations are objects with a name only; a section that is present is not empty.
  CHECK(error_pointer(changed(s, [](json& j) { j["configurations"][1]["weight"] = 1.0; })) ==
        "/configurations/1/weight");
  CHECK(error_pointer(changed(s, [](json& j) { j["configurations"][1] = "tele"; })) ==
        "/configurations/1");
  CHECK(error_pointer(changed(s, [](json& j) { j["configurations"] = json::array(); })) ==
        "/configurations");
}

TEST_CASE("migration 0.3 -> 0.4 drops pickups with io.pickup_dropped (ADR 0029, point 6)",
          "[io][0.4][migration]") {
  json j = file_0_3(json::parse(R"({"position": [0.0, 0.0, {"value": 6.0, "pickup": "2 * 3"}]})"));
  j["root"]["children"][1]["surfaces"][0]["pose"] =
      json::parse(R"({"position": [0.0, {"value": 1.0, "variable": true, "pickup": "x"}, 0.0]})");
  const std::string text = j.dump();
  std::vector<Diagnostic> warnings;
  const System s = rtt::io::parse_system(text, warnings);
  REQUIRE(s.schema_version == "0.4.0");
  // The value stays, the pickup is gone, variable is kept.
  REQUIRE(element_of(s, 1).pose.position[2] == Param(6.0));
  const Param& y = element_of(s, 1).surfaces[0].pose.position[1];
  REQUIRE(y.value == 1.0);
  REQUIRE(y.variable);
  REQUIRE_FALSE(y.is_bound());
  REQUIRE(warnings.size() == 2);
  for (const Diagnostic& w : warnings) {
    REQUIRE(w.code == "io.pickup_dropped");
    REQUIRE(w.severity == Severity::Warning);
  }
  // In file order, with the pointer into the file that was read and the dropped text.
  REQUIRE(warnings[0].location == "/root/children/1/pose/position/2/pickup");
  REQUIRE_THAT(warnings[0].message, ContainsSubstring("2 * 3"));
  REQUIRE(warnings[1].location == "/root/children/1/surfaces/0/pose/position/1/pickup");
  REQUIRE_THAT(warnings[1].message, ContainsSubstring("'x'"));
  // The form without warnings reads the same model; it drops the warnings (documented).
  REQUIRE(rtt::io::parse_system(text) == s);
  // A 0.4 file has no warnings.
  std::vector<Diagnostic> none;
  REQUIRE(rtt::io::parse_system(rtt::io::to_json(s), none) == s);
  REQUIRE(none.empty());
}

TEST_CASE("0.4 forms in an older file are errors, the version is never upgraded silently",
          "[io][0.4][migration]") {
  CHECK(error_pointer(file_0_3(json::parse(R"({"reference": "relative_to_sibling"})"))) ==
        "/root/children/1/pose/reference");
  CHECK(error_pointer(file_0_3(json::parse(R"({"order": "rotate_first"})"))) ==
        "/root/children/1/pose/order");
  CHECK(error_pointer(file_0_3(json::parse(R"({"position": [0.0, 0.0, {"param": "D"}]})"))) ==
        "/root/children/1/pose/position/2/param");
  CHECK(error_pointer(
            file_0_3(json::parse(R"({"position": [0.0, 0.0, {"value": 1.0, "min": 0.0}]})"))) ==
        "/root/children/1/pose/position/2/min");
  CHECK(error_pointer(file_0_3(json::object(),
                               json::parse(R"({"parameters": [{"name": "D", "value": 1.0}]})"))) ==
        "/parameters");
  CHECK(error_pointer(
            file_0_3(json::object(), json::parse(R"({"configurations": [{"name": "a"}]})"))) ==
        "/configurations");
  CHECK_THAT(error_message(file_0_3(json::parse(R"({"order": "rotate_first"})"))),
             ContainsSubstring("0.4"));
  // The same file as 0.4 is fine.
  json ok = file_0_3(json::parse(R"({"order": "rotate_first"})"));
  ok["schema_version"] = "0.4.0";
  CHECK(error_pointer(ok) == "<no error>");
}

TEST_CASE("load_system reports the warnings of the file as well", "[io][0.4][migration]") {
  const std::string text =
      file_0_3(json::parse(R"({"position": [0.0, 0.0, {"value": 6.0, "pickup": "2 * 3"}]})"))
          .dump();
  const auto path = std::filesystem::temp_directory_path() / "rtt_io_test_pickup_0_3.rtt.json";
  {
    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    out << text;
  }
  std::vector<Diagnostic> warnings;
  const System s = rtt::io::load_system(path, warnings);
  std::filesystem::remove(path);
  REQUIRE(warnings.size() == 1);
  REQUIRE(warnings[0].location == "/root/children/1/pose/position/2/pickup");
  REQUIRE(element_of(s, 1).pose.position[2] == Param(6.0));
}

TEST_CASE("the writer does not write the value of a bound Param", "[io][0.4]") {
  System s = rtt::model::test::make_singlet();
  s.parameters = {row("D", 106.363)};
  Param bound = Param::bound("D");
  bound.value = std::numeric_limits<double>::quiet_NaN();  // meaningless before evaluation
  element(s, 2).pose.position[2] = bound;
  const json j = json::parse(rtt::io::to_json(s));
  REQUIRE(j["root"]["children"][2]["pose"]["position"][2] == json::parse(R"({"param": "D"})"));
  // Read back, the value is 0: the reader sets it (ADR 0029, point 3).
  const System back = rtt::io::parse_system(j.dump());
  REQUIRE(element_of(back, 2).pose.position[2].value == 0.0);
  REQUIRE(element_of(back, 2).pose.position[2].param == "D");
}
