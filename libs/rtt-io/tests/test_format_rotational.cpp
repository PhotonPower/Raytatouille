// Schema 0.9 (#314, ADR 0043, point 4): the base shapes odd_asphere and axicon. The schema stays
// 0.9.0 (F1); the types are additive, and a file before 0.9 that holds them is a read error.

#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <set>
#include <string>
#include <string_view>
#include <variant>
#include <vector>

#include "rtt/io/edit.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/model/parameters.hpp"
#include "rtt/model/system.hpp"
#include "test_support.hpp"

using nlohmann::json;
using ojson = nlohmann::ordered_json;
using rtt::model::Axicon;
using rtt::model::BaseShape;
using rtt::model::OddAsphere;
using rtt::model::Param;
using rtt::model::System;

namespace {

/// Pointer of the base shape of L1.S1 in the singlet of make_singlet().
constexpr std::string_view kBase = "/root/children/1/surfaces/0/shape/base";

System with_base(const BaseShape& base) {
  System s = rtt::model::test::make_singlet();
  rtt::model::test::element(s, 1).surfaces[0].shape.base = base;
  return s;
}

OddAsphere odd_asphere() {
  Param variable(0.02);
  variable.variable = true;
  return {variable, Param(-0.5), {Param(1e-3), Param(0.0), Param(-2e-6)}};
}

Axicon axicon() {
  return {Param(-5.0)};
}

json base_json(const std::string& text) {
  return json::parse(text)[json::json_pointer(std::string(kBase))];
}

std::vector<std::string> keys(const ojson& o) {
  std::vector<std::string> out;
  for (const auto& [key, value] : o.items()) out.push_back(key);
  return out;
}

std::string error_pointer(const std::string& text) {
  try {
    [[maybe_unused]] const System s = rtt::io::parse_system(text);
  } catch (const rtt::io::ParseError& e) {
    return e.pointer();
  }
  return "<no error>";
}

std::string with_version(const System& s, const std::string& version) {
  json j = json::parse(rtt::io::to_json(s));
  j["schema_version"] = version;
  return j.dump();
}

std::string with_base_json(const System& s, const json& base) {
  json j = json::parse(rtt::io::to_json(s));
  j[json::json_pointer(std::string(kBase))] = base;
  return j.dump();
}

/// A Param object of the edit form (as in test_format_0_9.cpp).
bool param_object(const ojson& j) {
  if (!j.is_object()) return false;
  if (j.contains("param")) return j.size() == 1;
  if (!j.contains("value")) return false;
  for (const auto& [key, value] : j.items()) {
    if (key != "value" && key != "variable" && key != "min" && key != "max") return false;
  }
  return true;
}

void params_in(const ojson& j, const std::string& pointer, std::vector<std::string>& out) {
  if (j.is_object()) {
    if (param_object(j) && !pointer.starts_with("/parameters/")) {
      out.push_back(pointer);
      return;
    }
    for (const auto& [key, value] : j.items()) params_in(value, pointer + "/" + key, out);
  } else if (j.is_array()) {
    for (std::size_t i = 0; i < j.size(); ++i) {
      params_in(j[i], pointer + "/" + std::to_string(i), out);
    }
  }
}

}  // namespace

TEST_CASE("odd asphere and axicon: canonical form and round trip", "[io][0.9][rotational]") {
  SECTION("odd asphere: every key, in the order of the struct") {
    const System s = with_base(odd_asphere());
    const std::string text = rtt::io::to_json(s);
    CHECK(base_json(text) == json{{"type", "odd_asphere"},
                                  {"curvature", {{"value", 0.02}, {"variable", true}}},
                                  {"conic", -0.5},
                                  {"coefficients", {1e-3, 0.0, -2e-6}}});
    CHECK(keys(ojson::parse(text)[ojson::json_pointer(std::string(kBase))]) ==
          std::vector<std::string>{"type", "curvature", "conic", "coefficients"});
    CHECK(rtt::io::parse_system(text) == s);
    CHECK(rtt::io::to_json(rtt::io::parse_system(text)) == text);
  }
  SECTION("odd asphere: zero conic and no coefficients are not written; c = 0 is") {
    const System s = with_base(OddAsphere{Param(0.0), Param(0.0), {}});
    const std::string text = rtt::io::to_json(s);
    CHECK(base_json(text) == json{{"type", "odd_asphere"}, {"curvature", 0.0}});
    CHECK(rtt::io::parse_system(text) == s);
  }
  SECTION("axicon") {
    const System s = with_base(axicon());
    const std::string text = rtt::io::to_json(s);
    CHECK(base_json(text) == json{{"type", "axicon"}, {"angle_deg", -5.0}});
    CHECK(rtt::io::parse_system(text) == s);
    CHECK(rtt::io::to_json(rtt::io::parse_system(text)) == text);
  }
  SECTION("bound curvature and angle read back as references") {
    for (const BaseShape& base :
         {BaseShape{OddAsphere{Param::bound("P"), Param(0.0), {Param(1e-3)}}},
          BaseShape{Axicon{Param::bound("P")}}}) {
      System s = with_base(base);
      rtt::model::ParameterRow row;
      row.name = "P";
      row.form = 0.01;
      s.parameters = {row};
      const std::string text = rtt::io::to_json(s);
      const std::string key = std::holds_alternative<Axicon>(base) ? "angle_deg" : "curvature";
      CHECK(base_json(text)[key] == json{{"param", "P"}});
      CHECK(rtt::io::parse_system(text) == s);
    }
  }
}

TEST_CASE("odd asphere and axicon: the edit form has every key as a Param",
          "[io][0.9][rotational]") {
  struct Case {
    BaseShape base;
    std::vector<std::string> keys;
    std::vector<std::string> params;  ///< relative to the base
  };
  const std::vector<Case> cases = {
      {OddAsphere{Param(0.0), Param(0.0), {}},
       {"type", "curvature", "conic", "coefficients"},
       {"/curvature", "/conic"}},
      {odd_asphere(),
       {"type", "curvature", "conic", "coefficients"},
       {"/curvature", "/conic", "/coefficients/0", "/coefficients/1", "/coefficients/2"}},
      {axicon(), {"type", "angle_deg"}, {"/angle_deg"}},
  };
  for (const Case& c : cases) {
    const System s = with_base(c.base);
    const std::string text = rtt::io::to_edit_json(s);
    const ojson base = ojson::parse(text)[ojson::json_pointer(std::string(kBase))];
    INFO(base.dump());
    CHECK(keys(base) == c.keys);
    std::vector<std::string> found;
    params_in(base, "", found);
    CHECK(found == c.params);
    CHECK(rtt::io::parse_system(text) == s);
  }
}

TEST_CASE(
    "for_each_param visits the Params of odd asphere and axicon as the edit form writes "
    "them",
    "[io][0.9][rotational][params]") {
  for (const BaseShape& base : {BaseShape{odd_asphere()}, BaseShape{axicon()}}) {
    const System s = with_base(base);
    std::vector<std::string> expected;
    params_in(ojson::parse(rtt::io::to_edit_json(s)), "", expected);
    std::vector<std::string> visited;
    rtt::model::for_each_param(
        s, [&](std::string_view pointer, const Param&) { visited.emplace_back(pointer); });
    CHECK(visited == expected);
    CHECK(std::set<std::string>(visited.begin(), visited.end()).size() == visited.size());
    std::size_t in_base = 0;
    for (const std::string& p : visited) in_base += p.starts_with(kBase) ? 1 : 0;
    CHECK(in_base == (std::holds_alternative<Axicon>(base) ? 1U : 5U));
  }
}

TEST_CASE("odd asphere and axicon: missing, unknown and malformed keys are read errors",
          "[io][0.9][rotational]") {
  const System s = rtt::model::test::make_singlet();
  const std::string base(kBase);
  struct Bad {
    json value;
    std::string pointer;
  };
  const std::vector<Bad> bad = {
      {{{"type", "odd_asphere"}}, base},
      {{{"type", "odd_asphere"}, {"radius", 10.0}}, base},
      {{{"type", "odd_asphere"}, {"curvature", 0.1}, {"radius", 10.0}}, base},
      {{{"type", "odd_asphere"}, {"curvature", "0.1"}}, base + "/curvature"},
      {{{"type", "odd_asphere"}, {"curvature", 0.1}, {"coefficients", 1.0}},
       base + "/coefficients"},
      {{{"type", "axicon"}}, base},
      {{{"type", "axicon"}, {"angle_deg", 5.0}, {"curvature", 0.1}}, base},
      {{{"type", "axicon"}, {"angle", 5.0}}, base},
      {{{"type", "axicon"}, {"angle_deg", true}}, base + "/angle_deg"},
  };
  for (const Bad& b : bad) {
    INFO(b.value.dump());
    CHECK(error_pointer(with_base_json(s, b.value)).starts_with(b.pointer));
  }
  try {
    (void)rtt::io::parse_system(with_base_json(s, json{{"type", "cone"}, {"angle_deg", 1.0}}));
    FAIL("no ParseError");
  } catch (const rtt::io::ParseError& e) {
    const std::string what = e.what();
    for (const char* name : {"odd_asphere", "axicon"}) {
      INFO(name);
      CHECK(what.find(name) != std::string::npos);
    }
  }
}

TEST_CASE("odd asphere and axicon in a file before 0.9 are an error at the type",
          "[io][0.9][rotational]") {
  // ADR 0043, point 4: as check_0_9_forms for the shapes of S2.
  for (const BaseShape& base : {BaseShape{odd_asphere()}, BaseShape{axicon()}}) {
    const System s = with_base(base);
    INFO(rtt::io::to_json(s));
    for (const std::string version : {"0.8.0", "0.7.0", "0.5.0", "0.2.0"}) {
      INFO(version);
      CHECK(error_pointer(with_version(s, version)) == std::string(kBase) + "/type");
    }
    CHECK(error_pointer(with_version(s, "0.9.0")) == "<no error>");
  }
  // Inside an assembly.
  System s = rtt::model::test::make_singlet();
  rtt::model::Assembly group;
  group.name = "group";
  group.children = {s.root.children[1]};
  std::get<rtt::model::Element>(group.children[0].value).surfaces[0].shape.base = axicon();
  s.root.children[1] = rtt::model::Node{group};
  CHECK(error_pointer(with_version(s, "0.8.0")) ==
        "/root/children/1/children/0/surfaces/0/shape/base/type");
  CHECK(error_pointer(with_version(s, "0.9.0")) == "<no error>");
}
