// Schema 0.9 (#315, ADR 0042, point 5): the anamorphic base shapes biconic, toroid and cylinder;
// the migration 0.8 -> 0.9 changes only the version.

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
using rtt::model::BaseShape;
using rtt::model::Biconic;
using rtt::model::Cylinder;
using rtt::model::Param;
using rtt::model::System;
using rtt::model::Toroid;

namespace {

/// Pointer of the base shape of L1.S1 in the singlet of make_singlet().
constexpr std::string_view kBase = "/root/children/1/surfaces/0/shape/base";

/// The singlet with `base` at L1.S1.
System with_base(const BaseShape& base) {
  System s = rtt::model::test::make_singlet();
  rtt::model::test::element(s, 1).surfaces[0].shape.base = base;
  return s;
}

Biconic biconic() {
  return {Param(40.0), Param(-60.0), Param(-0.5), Param(0.25)};
}

Toroid toroid() {
  Param variable(55.0);
  variable.variable = true;
  return {Param(48.0), Param(-1.2), {Param(2e-6), Param(-3e-9)}, variable};
}

Cylinder cylinder() {
  return {Param(-75.0), Param(0.0), {Param(1e-6)}};
}

/// The base shape of L1.S1 in the file `text`.
json base_json(const std::string& text) {
  return json::parse(text)[json::json_pointer(std::string(kBase))];
}

/// Keys of `o` in written order.
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

/// The canonical file of `s` with the schema version `version`.
std::string with_version(const System& s, const std::string& version) {
  json j = json::parse(rtt::io::to_json(s));
  j["schema_version"] = version;
  return j.dump();
}

/// The canonical file of `s` with the base shape of L1.S1 replaced by `base`.
std::string with_base_json(const System& s, const json& base) {
  json j = json::parse(rtt::io::to_json(s));
  j[json::json_pointer(std::string(kBase))] = base;
  return j.dump();
}

/// A Param object of the edit form: {"param": …} alone, or "value" with keys only from
/// {value, variable, min, max} (as in test_for_each_param.cpp).
bool param_object(const ojson& j) {
  if (!j.is_object()) return false;
  if (j.contains("param")) return j.size() == 1;
  if (!j.contains("value")) return false;
  for (const auto& [key, value] : j.items()) {
    if (key != "value" && key != "variable" && key != "min" && key != "max") return false;
  }
  return true;
}

/// Pointers of the Param objects under `j`, in document order.
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

TEST_CASE("schema version is 0.9.0", "[io][0.9]") {
  REQUIRE(rtt::model::kSchemaVersion == "0.9.0");
}

TEST_CASE("migration 0.8 -> 0.9: content unchanged, the version becomes 0.9.0", "[io][0.9]") {
  const System s = rtt::model::test::make_singlet();
  const System read = rtt::io::parse_system(with_version(s, "0.8.0"));
  CHECK(read.schema_version == "0.9.0");
  CHECK(read == s);
  CHECK(json::parse(rtt::io::to_json(read))["schema_version"] == "0.9.0");
}

TEST_CASE("biconic, toroid and cylinder: canonical form and round trip", "[io][0.9]") {
  SECTION("biconic: every key, in the order of the struct") {
    const System s = with_base(biconic());
    const std::string text = rtt::io::to_json(s);
    const json b = base_json(text);
    CHECK(b == json{{"type", "biconic"},
                    {"radius_x", 40.0},
                    {"radius_y", -60.0},
                    {"conic_x", -0.5},
                    {"conic_y", 0.25}});
    CHECK(keys(ojson::parse(text)[ojson::json_pointer(std::string(kBase))]) ==
          std::vector<std::string>{"type", "radius_x", "radius_y", "conic_x", "conic_y"});
    CHECK(rtt::io::parse_system(text) == s);
    CHECK(rtt::io::to_json(rtt::io::parse_system(text)) == text);
  }
  SECTION("biconic: zero conics are not written and read back as 0") {
    const System s = with_base(Biconic{Param(40.0), Param(-60.0), Param(0.0), Param(0.0)});
    const std::string text = rtt::io::to_json(s);
    CHECK(base_json(text) == json{{"type", "biconic"}, {"radius_x", 40.0}, {"radius_y", -60.0}});
    CHECK(rtt::io::parse_system(text) == s);
  }
  SECTION("toroid: every key, a variable radius of rotation") {
    const System s = with_base(toroid());
    const std::string text = rtt::io::to_json(s);
    CHECK(base_json(text) == json{{"type", "toroid"},
                                  {"radius", 48.0},
                                  {"conic", -1.2},
                                  {"coefficients", {2e-6, -3e-9}},
                                  {"radius_of_rotation", {{"value", 55.0}, {"variable", true}}}});
    CHECK(
        keys(ojson::parse(text)[ojson::json_pointer(std::string(kBase))]) ==
        std::vector<std::string>{"type", "radius", "conic", "coefficients", "radius_of_rotation"});
    CHECK(rtt::io::parse_system(text) == s);
    CHECK(rtt::io::to_json(rtt::io::parse_system(text)) == text);
  }
  SECTION("toroid: without conic and coefficients") {
    const System s = with_base(Toroid{Param(48.0), Param(0.0), {}, Param(-30.0)});
    const std::string text = rtt::io::to_json(s);
    CHECK(base_json(text) ==
          json{{"type", "toroid"}, {"radius", 48.0}, {"radius_of_rotation", -30.0}});
    CHECK(rtt::io::parse_system(text) == s);
  }
  SECTION("cylinder: every key; without conic and coefficients") {
    const System s = with_base(cylinder());
    const std::string text = rtt::io::to_json(s);
    CHECK(base_json(text) ==
          json{{"type", "cylinder"}, {"radius", -75.0}, {"coefficients", {1e-6}}});
    CHECK(rtt::io::parse_system(text) == s);
    const System plain = with_base(Cylinder{Param(-75.0), Param(0.0), {}});
    CHECK(base_json(rtt::io::to_json(plain)) == json{{"type", "cylinder"}, {"radius", -75.0}});
    CHECK(rtt::io::parse_system(rtt::io::to_json(plain)) == plain);
  }
  SECTION("bound radii read back as references") {
    System s = with_base(Biconic{Param::bound("RX"), Param(-60.0), Param(0.0), Param(0.0)});
    rtt::model::ParameterRow row;
    row.name = "RX";
    row.form = 40.0;
    s.parameters = {row};
    const std::string text = rtt::io::to_json(s);
    CHECK(base_json(text)["radius_x"] == json{{"param", "RX"}});
    CHECK(rtt::io::parse_system(text) == s);
  }
}

TEST_CASE("biconic, toroid and cylinder: the edit form has every key as a Param", "[io][0.9]") {
  struct Case {
    BaseShape base;
    std::vector<std::string> keys;
    std::vector<std::string> params;  ///< relative to the base
  };
  const std::vector<Case> cases = {
      {Biconic{Param(40.0), Param(-60.0), Param(0.0), Param(0.0)},
       {"type", "radius_x", "radius_y", "conic_x", "conic_y"},
       {"/radius_x", "/radius_y", "/conic_x", "/conic_y"}},
      {Toroid{Param(48.0), Param(0.0), {}, Param(55.0)},
       {"type", "radius", "conic", "coefficients", "radius_of_rotation"},
       {"/radius", "/conic", "/radius_of_rotation"}},
      {toroid(),
       {"type", "radius", "conic", "coefficients", "radius_of_rotation"},
       {"/radius", "/conic", "/coefficients/0", "/coefficients/1", "/radius_of_rotation"}},
      {Cylinder{Param(-75.0), Param(0.0), {}},
       {"type", "radius", "conic", "coefficients"},
       {"/radius", "/conic"}},
      {cylinder(),
       {"type", "radius", "conic", "coefficients"},
       {"/radius", "/conic", "/coefficients/0"}},
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
    // The edit form reads back to the same system.
    CHECK(rtt::io::parse_system(text) == s);
  }
}

TEST_CASE("for_each_param visits the Params of the new shapes as the edit form writes them",
          "[io][0.9][params]") {
  for (const BaseShape& base : {BaseShape{biconic()}, BaseShape{toroid()}, BaseShape{cylinder()}}) {
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
    CHECK(in_base == (std::holds_alternative<Biconic>(base)  ? 4U
                      : std::holds_alternative<Toroid>(base) ? 5U
                                                             : 3U));
  }
}

TEST_CASE("the new shapes: missing, unknown and malformed keys are read errors", "[io][0.9]") {
  const System s = rtt::model::test::make_singlet();
  const std::string base(kBase);
  struct Bad {
    json value;
    std::string pointer;
  };
  const std::vector<Bad> bad = {
      {{{"type", "biconic"}, {"radius_y", 10.0}}, base},
      {{{"type", "biconic"}, {"radius_x", 10.0}}, base},
      {{{"type", "biconic"}, {"radius_x", 10.0}, {"radius_y", 10.0}, {"radius", 10.0}}, base},
      {{{"type", "biconic"}, {"radius_x", "10"}, {"radius_y", 10.0}}, base + "/radius_x"},
      {{{"type", "toroid"}, {"radius", 10.0}}, base},
      {{{"type", "toroid"}, {"radius_of_rotation", 10.0}}, base},
      {{{"type", "toroid"}, {"radius", 10.0}, {"radius_of_rotation", 5.0}, {"conic_x", 0.0}}, base},
      {{{"type", "toroid"}, {"radius", 10.0}, {"radius_of_rotation", 5.0}, {"coefficients", 1.0}},
       base + "/coefficients"},
      {{{"type", "cylinder"}}, base},
      {{{"type", "cylinder"}, {"radius", 10.0}, {"radius_of_rotation", 5.0}}, base},
  };
  for (const Bad& b : bad) {
    INFO(b.value.dump());
    CHECK(error_pointer(with_base_json(s, b.value)).starts_with(b.pointer));
  }
  // The message of an unknown type lists the new types.
  try {
    (void)rtt::io::parse_system(with_base_json(s, json{{"type", "torus"}, {"radius", 1.0}}));
    FAIL("no ParseError");
  } catch (const rtt::io::ParseError& e) {
    const std::string what = e.what();
    for (const char* name : {"biconic", "toroid", "cylinder"}) {
      INFO(name);
      CHECK(what.find(name) != std::string::npos);
    }
  }
}

TEST_CASE("the new shapes in a file before 0.9 are an error at the type", "[io][0.9]") {
  // A version is never silently upgraded (as check_0_5_forms): the pointer is the type of the
  // base shape, also inside an assembly.
  for (const BaseShape& base : {BaseShape{biconic()}, BaseShape{toroid()}, BaseShape{cylinder()}}) {
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
  std::get<rtt::model::Element>(group.children[0].value).surfaces[0].shape.base = cylinder();
  s.root.children[1] = rtt::model::Node{group};
  CHECK(error_pointer(with_version(s, "0.8.0")) ==
        "/root/children/1/children/0/surfaces/0/shape/base/type");
  CHECK(error_pointer(with_version(s, "0.9.0")) == "<no error>");
}
