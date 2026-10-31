// Schema 0.8 (#288, ADR 0038): the reference of the object distance, object.reference; the
// migration 0.7 -> 0.8 changes only the version.

#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <optional>
#include <string>

#include "rtt/io/edit.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/model/system.hpp"
#include "test_support.hpp"

using nlohmann::json;
using rtt::model::Param;
using rtt::model::SurfaceId;
using rtt::model::System;

namespace {

/// The singlet with a finite object 50 mm in front of the origin or of `reference`.
System finite_object(const std::optional<SurfaceId>& reference) {
  System s = rtt::model::test::make_singlet();
  s.object.at_infinity = false;
  s.object.distance = Param(50.0);
  s.object.reference = reference;
  return s;
}

std::string error_pointer(const std::string& text) {
  try {
    [[maybe_unused]] const System s = rtt::io::parse_system(text);
  } catch (const rtt::io::ParseError& e) {
    return e.pointer();
  }
  return "<no error>";
}

/// The canonical file of `s` with the schema version and `object` replaced.
std::string with_version_and_object(const System& s,
                                    const std::string& version,
                                    const json& object) {
  json j = json::parse(rtt::io::to_json(s));
  j["schema_version"] = version;
  j["object"] = object;
  return j.dump();
}

const json kFinite = {{"at_infinity", false}, {"distance", 50.0}};

}  // namespace

TEST_CASE("schema version is 0.8.0", "[io][0.8]") {
  REQUIRE(rtt::model::kSchemaVersion == "0.8.0");
}

TEST_CASE("migration 0.7 -> 0.8: content unchanged, the version becomes 0.8.0", "[io][0.8]") {
  const System s = finite_object(std::nullopt);
  const System read = rtt::io::parse_system(with_version_and_object(s, "0.7.0", kFinite));
  CHECK(read.schema_version == "0.8.0");
  CHECK(read == s);
  CHECK(json::parse(rtt::io::to_json(read))["schema_version"] == "0.8.0");
}

TEST_CASE("object.reference: the reference form round trip", "[io][0.8]") {
  const System s = finite_object(SurfaceId("L1.S1"));
  const json j = json::parse(rtt::io::to_json(s));
  CHECK(j["object"]["reference"] == json{{"surface", "L1.S1"}});
  // Written after medium (ADR 0038, point 8): at_infinity, distance, reference.
  json keys = json::array();
  for (const auto& [key, value] : j["object"].items()) keys.push_back(key);
  CHECK(keys == json{"at_infinity", "distance", "reference"});
  const System back = rtt::io::parse_system(rtt::io::to_json(s));
  CHECK(back == s);
  REQUIRE(back.object.reference.has_value());
  CHECK(back.object.reference->str() == "L1.S1");
}

TEST_CASE("object.reference: \"origin\" is the default and is not written", "[io][0.8]") {
  const System s = finite_object(std::nullopt);
  json object = kFinite;
  object["reference"] = "origin";
  const System read = rtt::io::parse_system(with_version_and_object(s, "0.8.0", object));
  CHECK(!read.object.reference.has_value());
  CHECK(read == s);
  CHECK(!json::parse(rtt::io::to_json(read))["object"].contains("reference"));
  // Nor in the edit form, like the medium (ADR 0038, point 8).
  CHECK(!json::parse(rtt::io::to_edit_json(read))["object"].contains("reference"));
}

TEST_CASE("object.reference: other forms are read errors at its pointer", "[io][0.8]") {
  const System s = finite_object(std::nullopt);
  for (const json& bad : {json("vertex"), json(1.0), json::object(), json{{"node", "L1"}},
                          json{{"surface", "L1.S1"}, {"node", "L1"}}, json{{"surface", 3}}}) {
    INFO(bad.dump());
    json object = kFinite;
    object["reference"] = bad;
    const std::string pointer = error_pointer(with_version_and_object(s, "0.8.0", object));
    CHECK(pointer.starts_with("/object/reference"));
  }
}

TEST_CASE("object.reference in a file before 0.8 is an error at its pointer", "[io][0.8]") {
  const System s = finite_object(std::nullopt);
  json object = kFinite;
  object["reference"] = json{{"surface", "L1.S1"}};
  CHECK(error_pointer(with_version_and_object(s, "0.7.0", object)) == "/object/reference");
  CHECK(error_pointer(with_version_and_object(s, "0.5.0", object)) == "/object/reference");
  object["reference"] = "origin";  // also the default form: a version is never silently upgraded
  CHECK(error_pointer(with_version_and_object(s, "0.7.0", object)) == "/object/reference");
  // the same file with the current version is fine
  object["reference"] = json{{"surface", "L1.S1"}};
  CHECK(error_pointer(with_version_and_object(s, "0.8.0", object)) == "<no error>");
}
