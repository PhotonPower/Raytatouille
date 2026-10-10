#include <catch2/catch_test_macros.hpp>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "rtt/io/json_io.hpp"
#include "test_support.hpp"

using rtt::model::Element;

namespace {

/// Canonical file with one cemented lens of three surfaces. `material` is the JSON text of the
/// material value, `version` the schema version.
std::string doublet_file(const std::string& version, const std::string& material) {
  return R"({
  "schema_version": ")" +
         version + R"(",
  "units": {"length": "mm", "wavelength": "um"},
  "wavelengths": [
    {"um": 0.5876, "reference": true}
  ],
  "aperture": {"type": "epd", "value": 10.0},
  "fields": {
    "points": [
      {}
    ]
  },
  "root": {
    "type": "assembly",
    "name": "system",
    "children": [
      {
        "type": "lens",
        "name": "L1",
        "material": )" +
         material + R"(,
        "surfaces": [
          {"id": "L1.S1"},
          {"id": "L1.S2"},
          {"id": "L1.S3"}
        ]
      }
    ]
  },
  "paths": [
    {"name": "main", "events": "auto"}
  ]
}
)";
}

const Element& lens(const rtt::model::System& s) {
  return std::get<Element>(s.root.children[0].value);
}

std::string error_pointer(const std::string& text) {
  try {
    (void)rtt::io::parse_system(text);
  } catch (const rtt::io::ParseError& e) {
    return e.pointer();
  }
  return "<no error>";
}

constexpr const char* kList = R"(["SCHOTT:N-BK7", "SCHOTT:F2"])";
constexpr const char* kShorthand = R"("SCHOTT:N-BK7")";

}  // namespace

TEST_CASE("schema version is 0.5.0", "[io][segments]") {
  REQUIRE(rtt::model::kSchemaVersion == "0.5.0");
}

TEST_CASE("material list: file -> model -> file is byte-identical", "[io][segments][roundtrip]") {
  const std::string text = doublet_file("0.5.0", kList);
  const rtt::model::System s = rtt::io::parse_system(text);
  REQUIRE_FALSE(lens(s).material.has_value());
  REQUIRE(lens(s).segment_materials == std::vector<std::string>{"SCHOTT:N-BK7", "SCHOTT:F2"});
  REQUIRE(rtt::io::to_json(s) == text);
}

TEST_CASE("material shorthand is kept on write", "[io][segments][roundtrip]") {
  const std::string text = doublet_file("0.5.0", kShorthand);
  const rtt::model::System s = rtt::io::parse_system(text);
  REQUIRE(lens(s).material == std::optional<std::string>("SCHOTT:N-BK7"));
  REQUIRE(lens(s).segment_materials.empty());
  REQUIRE(rtt::io::to_json(s) == text);
}

TEST_CASE("a one-entry list stays a list", "[io][segments][roundtrip]") {
  const std::string text = doublet_file("0.5.0", R"(["SCHOTT:N-BK7"])");
  const rtt::model::System s = rtt::io::parse_system(text);
  REQUIRE_FALSE(lens(s).material.has_value());
  REQUIRE(lens(s).segment_materials == std::vector<std::string>{"SCHOTT:N-BK7"});
  REQUIRE(rtt::io::to_json(s) == text);
}

TEST_CASE("migration: a 0.1 file with one material is written in the current version",
          "[io][segments]") {
  const std::string old_text = doublet_file("0.1.0", kShorthand);
  const rtt::model::System s = rtt::io::parse_system(old_text);
  REQUIRE(s.schema_version == "0.5.0");
  REQUIRE(lens(s).material == std::optional<std::string>("SCHOTT:N-BK7"));
  REQUIRE(lens(s).segment_materials.empty());
  REQUIRE(rtt::io::to_json(s) == doublet_file("0.5.0", kShorthand));
}

TEST_CASE("migration: a 0.1 file must not contain material lists", "[io][segments]") {
  REQUIRE(error_pointer(doublet_file("0.1.0", kList)) == "/root/children/0/material");
}

TEST_CASE("material list structural errors carry a JSON pointer", "[io][segments]") {
  REQUIRE(error_pointer(doublet_file("0.2.0", "[]")) == "/root/children/0/material");
  REQUIRE(error_pointer(doublet_file("0.2.0", R"(["SCHOTT:N-BK7", 2])")) ==
          "/root/children/0/material/1");
  REQUIRE(error_pointer(doublet_file("0.2.0", R"([["SCHOTT:N-BK7"]])")) ==
          "/root/children/0/material/0");
  REQUIRE(error_pointer(doublet_file("0.2.0", R"({"name": "SCHOTT:N-BK7"})")) ==
          "/root/children/0/material");
  REQUIRE(error_pointer(doublet_file("0.2.0", "true")) == "/root/children/0/material");
}

TEST_CASE("model -> file -> model preserves segment materials", "[io][segments][roundtrip]") {
  rtt::model::System s = rtt::model::test::make_singlet();
  Element& l = rtt::model::test::element(s, 1);
  l.material.reset();
  l.segment_materials = {"SCHOTT:N-BK7"};
  const std::string text = rtt::io::to_json(s);
  const rtt::model::System back = rtt::io::parse_system(text);
  REQUIRE(back == s);
  REQUIRE(rtt::io::to_json(back) == text);
}

TEST_CASE("shorthand and list together cannot be written", "[io][segments]") {
  rtt::model::System s = rtt::model::test::make_singlet();
  rtt::model::test::element(s, 1).segment_materials = {"SCHOTT:N-BK7"};
  REQUIRE_THROWS_AS(rtt::io::to_json(s), std::invalid_argument);
}
