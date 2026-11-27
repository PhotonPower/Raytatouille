// Schema 0.12 (#358, ADR 0055): object.slope = [sx, sy], two Params, the tilt of a finite object
// plane. The migration 0.11 -> 0.12 changes only the version; a slope in an older file is a read
// error at its pointer.

#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <string>

#include "rtt/io/edit.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/model/parameters.hpp"
#include "rtt/model/system.hpp"
#include "rtt/model/validate.hpp"
#include "test_support.hpp"

using nlohmann::json;
using rtt::model::Param;
using rtt::model::System;

namespace {

std::string error_pointer(const std::string& text) {
  try {
    [[maybe_unused]] const System s = rtt::io::parse_system(text);
  } catch (const rtt::io::ParseError& e) {
    return e.pointer();
  }
  return "<no error>";
}

/// The singlet with a finite object 100 mm in front, tilted by sy = -1.
System tilted() {
  System s = rtt::model::test::make_singlet();
  s.object.at_infinity = false;
  s.object.distance = Param(100.0);
  s.object.slope = {Param(0.25), Param(-1.0)};
  return s;
}

}  // namespace

TEST_CASE("schema version is 0.13.0", "[io][0.12]") {
  REQUIRE(rtt::model::kSchemaVersion == "0.13.0");
}

TEST_CASE("object.slope round-trips; default slopes are not written", "[io][0.12]") {
  const System s = tilted();
  const json j = json::parse(rtt::io::to_json(s));
  REQUIRE(j["object"]["slope"] == json::array({0.25, -1.0}));
  CHECK(rtt::io::parse_system(rtt::io::to_json(s)) == s);
  // without a tilt the canonical form has no slope, the edit form always (a Param, as distance)
  System plane = s;
  plane.object.slope = {};
  CHECK_FALSE(json::parse(rtt::io::to_json(plane))["object"].contains("slope"));
  CHECK(json::parse(rtt::io::to_edit_json(plane))["object"]["slope"].size() == 2);
  // a slope bound to a row of the parameter table (ADR 0029)
  System bound = s;
  rtt::model::ParameterRow row;
  row.name = "tilt";
  row.form = 0.5;
  bound.parameters.push_back(row);
  bound.object.slope[1] = Param::bound("tilt");
  CHECK(rtt::io::parse_system(rtt::io::to_json(bound)) == bound);
  CHECK(rtt::model::resolve_parameters(bound, 0).object.slope[1].value == 0.5);
  CHECK(rtt::model::validate(bound).empty());
}

TEST_CASE("migration 0.11 -> 0.12: content unchanged, a slope needs 0.12", "[io][0.12]") {
  System plane = tilted();
  plane.object.slope = {};
  json j = json::parse(rtt::io::to_json(plane));
  j["schema_version"] = "0.11.0";
  const System read = rtt::io::parse_system(j.dump());
  CHECK(read.schema_version == "0.13.0");
  CHECK(read == plane);
  json old = json::parse(rtt::io::to_json(tilted()));
  old["schema_version"] = "0.11.0";
  CHECK(error_pointer(old.dump()) == "/object/slope");
  // the form of the key
  json bad = json::parse(rtt::io::to_json(tilted()));
  bad["object"]["slope"] = json::array({1.0});
  CHECK(error_pointer(bad.dump()) == "/object/slope");
  bad["object"]["slope"] = json::array({1.0, "a"});
  CHECK(error_pointer(bad.dump()) == "/object/slope/1");
}
