// aperture.reading (#296, ADR 0040, points 1 and 7): "paraxial" (default) or "real", only for
// object_na and only from schema 0.8; written only if real.

#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <string>

#include "rtt/io/edit.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/model/system.hpp"
#include "test_support.hpp"

using nlohmann::json;
using rtt::model::ApertureReading;
using rtt::model::Param;
using rtt::model::System;
using rtt::model::SystemApertureType;

namespace {

/// The singlet with a finite object and object_na 0.1, read as `reading`.
System na_system(ApertureReading reading) {
  System s = rtt::model::test::make_singlet();
  s.object.at_infinity = false;
  s.object.distance = Param(100.0);
  s.fields = {rtt::model::FieldType::ObjectHeight, {{0.0, 0.0, 1.0}}};
  s.aperture = {SystemApertureType::ObjectSpaceNA, Param(0.1)};
  s.aperture.reading = reading;
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

}  // namespace

TEST_CASE("aperture.reading: real is written and read back", "[io][0.8][aperture-reading]") {
  const System s = na_system(ApertureReading::Real);
  const nlohmann::ordered_json o = nlohmann::ordered_json::parse(rtt::io::to_json(s));
  CHECK(o["aperture"]["reading"] == "real");
  json keys = json::array();
  for (const auto& [key, value] : o["aperture"].items()) keys.push_back(key);
  CHECK(keys == json{"type", "value", "reading"});
  const System back = rtt::io::parse_system(rtt::io::to_json(s));
  CHECK(back == s);
  CHECK(back.aperture.reading == ApertureReading::Real);
}

TEST_CASE("aperture.reading: paraxial is the default and is not written",
          "[io][0.8][aperture-reading]") {
  const System s = na_system(ApertureReading::Paraxial);
  CHECK_FALSE(json::parse(rtt::io::to_json(s))["aperture"].contains("reading"));
  CHECK_FALSE(json::parse(rtt::io::to_edit_json(s))["aperture"].contains("reading"));
  // An explicit "paraxial" is read and dropped by rtt format.
  json j = json::parse(rtt::io::to_json(s));
  j["aperture"]["reading"] = "paraxial";
  const System read = rtt::io::parse_system(j.dump());
  CHECK(read == s);
  CHECK_FALSE(json::parse(rtt::io::to_json(read))["aperture"].contains("reading"));
}

TEST_CASE("aperture.reading: read errors", "[io][0.8][aperture-reading]") {
  json j = json::parse(rtt::io::to_json(na_system(ApertureReading::Real)));
  // An unknown value.
  j["aperture"]["reading"] = "exact";
  CHECK(error_pointer(j.dump()) == "/aperture/reading");
  // Another aperture type, also with "paraxial".
  for (const char* value : {"real", "paraxial"}) {
    INFO(value);
    json e = json::parse(rtt::io::to_json(rtt::model::test::make_singlet()));
    e["aperture"]["reading"] = value;
    CHECK(error_pointer(e.dump()) == "/aperture/reading");
  }
  // A file before 0.8.
  json old = json::parse(rtt::io::to_json(na_system(ApertureReading::Real)));
  old["schema_version"] = "0.7.0";
  CHECK(error_pointer(old.dump()) == "/aperture/reading");
  // Content guard: the same 0.7 file without reading is read.
  old["aperture"].erase("reading");
  CHECK(error_pointer(old.dump()) == "<no error>");
}
