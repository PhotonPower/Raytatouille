// Schema 0.7 (#282, ADR 0037): the medium of the object space, object.medium; the migration
// 0.6 -> 0.7 changes only the version.

#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <optional>
#include <string>

#include "rtt/io/json_io.hpp"
#include "rtt/model/system.hpp"
#include "test_support.hpp"

using nlohmann::json;
using rtt::model::Param;
using rtt::model::System;

namespace {

/// The singlet with a finite object at 50 mm, with or without an object medium.
System finite_object(const std::optional<std::string>& medium) {
  System s = rtt::model::test::make_singlet();
  s.object.at_infinity = false;
  s.object.distance = Param(50.0);
  s.object.medium = medium;
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

}  // namespace

TEST_CASE("schema version is 0.8.0", "[io][0.7]") {
  REQUIRE(rtt::model::kSchemaVersion == "0.8.0");
}

TEST_CASE("object.medium: canonical form and round trip", "[io][0.7]") {
  SECTION("without a medium the file has no key, as before") {
    const System s = finite_object(std::nullopt);
    const json j = json::parse(rtt::io::to_json(s));
    CHECK_FALSE(j["object"].contains("medium"));
    CHECK(rtt::io::parse_system(rtt::io::to_json(s)) == s);
  }
  SECTION("with a medium") {
    const System s = finite_object("CONST:1.53");
    const json j = json::parse(rtt::io::to_json(s));
    CHECK(j["object"]["medium"] == "CONST:1.53");
    const System back = rtt::io::parse_system(rtt::io::to_json(s));
    CHECK(back == s);
    REQUIRE(back.object.medium.has_value());
    CHECK(*back.object.medium == "CONST:1.53");
  }
  SECTION("a medium makes the system differ from one without") {
    CHECK_FALSE(finite_object("CONST:1.53") == finite_object(std::nullopt));
  }
}

TEST_CASE("migration 0.6 -> 0.7: content unchanged, the version becomes 0.8.0", "[io][0.7]") {
  const System s = finite_object(std::nullopt);
  const System read = rtt::io::parse_system(
      with_version_and_object(s, "0.6.0", json{{"at_infinity", false}, {"distance", 50.0}}));
  CHECK(read.schema_version == "0.8.0");
  CHECK(read == s);
  CHECK(json::parse(rtt::io::to_json(read))["schema_version"] == "0.8.0");
}

TEST_CASE("object.medium in a file before 0.7 is an error at its pointer", "[io][0.7]") {
  const System s = finite_object(std::nullopt);
  const json object = {{"at_infinity", false}, {"distance", 50.0}, {"medium", "CONST:1.53"}};
  CHECK(error_pointer(with_version_and_object(s, "0.6.0", object)) == "/object/medium");
  CHECK(error_pointer(with_version_and_object(s, "0.5.0", object)) == "/object/medium");
  // the same file with the current version is fine
  CHECK(error_pointer(with_version_and_object(s, "0.7.0", object)) == "<no error>");
}

TEST_CASE("object.medium must be a string", "[io][0.7]") {
  const System s = finite_object(std::nullopt);
  const json object = {{"at_infinity", false}, {"distance", 50.0}, {"medium", 1.53}};
  CHECK(error_pointer(with_version_and_object(s, "0.7.0", object)) == "/object/medium");
}
