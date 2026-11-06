// Schema 0.9 (#314, ADR 0043, point 4): the base shapes qcon and qbfs. The schema stays 0.9.0
// (F1); the types are additive, and a file before 0.9 that holds them is a read error.

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
using rtt::model::Param;
using rtt::model::QBfs;
using rtt::model::QCon;
using rtt::model::System;

namespace {

/// Pointer of the base shape of L1.S1 in the singlet of make_singlet().
constexpr std::string_view kBase = "/root/children/1/surfaces/0/shape/base";

System with_base(const BaseShape& base) {
  System s = rtt::model::test::make_singlet();
  rtt::model::test::element(s, 1).surfaces[0].shape.base = base;
  return s;
}

QCon qcon() {
  Param variable(0.02);
  variable.variable = true;
  return {variable, Param(-0.5), Param(12.5), {Param(1e-3), Param(0.0), Param(-2e-6)}};
}

QBfs qbfs() {
  return {Param(-0.01), Param(0.25), Param(8.0), {Param(3e-3)}};
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

TEST_CASE("qcon and qbfs: canonical form and round trip", "[io][0.9][forbes]") {
  SECTION("qcon: every key, in the order of the struct") {
    const System s = with_base(qcon());
    const std::string text = rtt::io::to_json(s);
    CHECK(base_json(text) == json{{"type", "qcon"},
                                  {"curvature", {{"value", 0.02}, {"variable", true}}},
                                  {"conic", -0.5},
                                  {"normalization_radius", 12.5},
                                  {"coefficients", {1e-3, 0.0, -2e-6}}});
    CHECK(keys(ojson::parse(text)[ojson::json_pointer(std::string(kBase))]) ==
          std::vector<std::string>{"type", "curvature", "conic", "normalization_radius",
                                   "coefficients"});
    CHECK(rtt::io::parse_system(text) == s);
    CHECK(rtt::io::to_json(rtt::io::parse_system(text)) == text);
  }
  SECTION("qbfs: every key") {
    const System s = with_base(qbfs());
    const std::string text = rtt::io::to_json(s);
    CHECK(base_json(text) == json{{"type", "qbfs"},
                                  {"curvature", -0.01},
                                  {"conic", 0.25},
                                  {"normalization_radius", 8.0},
                                  {"coefficients", {3e-3}}});
    CHECK(rtt::io::parse_system(text) == s);
    CHECK(rtt::io::to_json(rtt::io::parse_system(text)) == text);
  }
  SECTION("zero conic and no coefficients are not written; c = 0 is") {
    for (const BaseShape& base : {BaseShape{QCon{Param(0.0), Param(0.0), Param(5.0), {}}},
                                  BaseShape{QBfs{Param(0.0), Param(0.0), Param(5.0), {}}}}) {
      const System s = with_base(base);
      const std::string text = rtt::io::to_json(s);
      const std::string type = std::holds_alternative<QCon>(base) ? "qcon" : "qbfs";
      CHECK(base_json(text) ==
            json{{"type", type}, {"curvature", 0.0}, {"normalization_radius", 5.0}});
      CHECK(rtt::io::parse_system(text) == s);
    }
  }
}

TEST_CASE("qcon and qbfs: the edit form has every key as a Param", "[io][0.9][forbes]") {
  for (const BaseShape& base : {BaseShape{qcon()}, BaseShape{qbfs()}}) {
    const System s = with_base(base);
    const std::string text = rtt::io::to_edit_json(s);
    const ojson b = ojson::parse(text)[ojson::json_pointer(std::string(kBase))];
    INFO(b.dump());
    CHECK(keys(b) == std::vector<std::string>{"type", "curvature", "conic", "normalization_radius",
                                              "coefficients"});
    std::vector<std::string> found;
    params_in(b, "", found);
    std::vector<std::string> expected = {"/curvature", "/conic", "/normalization_radius"};
    const std::size_t n = std::holds_alternative<QCon>(base) ? 3 : 1;
    for (std::size_t i = 0; i < n; ++i) expected.push_back("/coefficients/" + std::to_string(i));
    CHECK(found == expected);
    CHECK(rtt::io::parse_system(text) == s);
    // for_each_param visits them as the edit form writes them.
    std::vector<std::string> all;
    params_in(ojson::parse(text), "", all);
    std::vector<std::string> visited;
    rtt::model::for_each_param(
        s, [&](std::string_view pointer, const Param&) { visited.emplace_back(pointer); });
    CHECK(visited == all);
  }
}

TEST_CASE("qcon and qbfs: missing, unknown and malformed keys are read errors",
          "[io][0.9][forbes]") {
  const System s = rtt::model::test::make_singlet();
  const std::string base(kBase);
  struct Bad {
    json value;
    std::string pointer;
  };
  const std::vector<Bad> bad = {
      {{{"type", "qcon"}, {"curvature", 0.1}}, base},
      {{{"type", "qcon"}, {"normalization_radius", 5.0}}, base},
      {{{"type", "qbfs"}, {"curvature", 0.1}, {"normalization_radius", 5.0}, {"radius", 1.0}},
       base},
      {{{"type", "qbfs"}, {"curvature", 0.1}, {"normalization_radius", "5"}},
       base + "/normalization_radius"},
      {{{"type", "qcon"}, {"curvature", 0.1}, {"normalization_radius", 5.0}, {"coefficients", 1}},
       base + "/coefficients"},
  };
  for (const Bad& b : bad) {
    INFO(b.value.dump());
    CHECK(error_pointer(with_base_json(s, b.value)).starts_with(b.pointer));
  }
  try {
    (void)rtt::io::parse_system(with_base_json(s, json{{"type", "q"}, {"curvature", 1.0}}));
    FAIL("no ParseError");
  } catch (const rtt::io::ParseError& e) {
    const std::string what = e.what();
    for (const char* name : {"qcon", "qbfs"}) {
      INFO(name);
      CHECK(what.find(name) != std::string::npos);
    }
  }
}

TEST_CASE("qcon and qbfs in a file before 0.9 are an error at the type", "[io][0.9][forbes]") {
  for (const BaseShape& base : {BaseShape{qcon()}, BaseShape{qbfs()}}) {
    const System s = with_base(base);
    for (const std::string version : {"0.8.0", "0.7.0", "0.2.0"}) {
      INFO(version);
      CHECK(error_pointer(with_version(s, version)) == std::string(kBase) + "/type");
    }
    CHECK(error_pointer(with_version(s, "0.9.0")) == "<no error>");
  }
}
