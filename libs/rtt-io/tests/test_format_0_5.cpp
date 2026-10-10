// Schema 0.5 (#178, ADR 0031, point 1): the interactions ideal_lens and ideal_cylinder_lens of a
// thin_element; the migration 0.4 -> 0.5 changes only the version.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <functional>
#include <nlohmann/json.hpp>
#include <optional>
#include <string>
#include <variant>

#include "rtt/io/edit.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/model/system.hpp"
#include "test_support.hpp"

using Catch::Matchers::ContainsSubstring;
using nlohmann::json;
using rtt::model::Element;
using rtt::model::ElementKind;
using rtt::model::IdealCylinderLens;
using rtt::model::IdealLens;
using rtt::model::Param;
using rtt::model::System;

namespace {

/// The singlet with a thin element "IL" at z = 50 (between the lens and the image) whose only
/// surface carries `interaction`.
System with_thin(const rtt::model::Interaction& interaction) {
  System s = rtt::model::test::make_singlet();
  Element thin;
  thin.name = "IL";
  thin.kind = ElementKind::ThinElement;
  thin.pose.position[2] = Param(50.0);
  rtt::model::Surface surface;
  surface.id = rtt::model::SurfaceId{"IL.S1"};
  surface.aperture = rtt::model::CircularAperture{10.0, 0.0};
  surface.interaction = interaction;
  thin.surfaces = {surface};
  s.root.children.insert(s.root.children.begin() + 2, rtt::model::Node{thin});
  return s;
}

/// The interaction of the thin element in the canonical file of `s`.
json interaction_json(const System& s) {
  const json j = json::parse(rtt::io::to_json(s));
  return j["root"]["children"][2]["surfaces"][0]["interaction"];
}

std::string error_pointer(const std::string& text) {
  try {
    [[maybe_unused]] const System s = rtt::io::parse_system(text);
  } catch (const rtt::io::ParseError& e) {
    return e.pointer();
  }
  return "<no error>";
}

/// The canonical file of `s` with its interaction replaced by `interaction`.
std::string with_interaction_json(const System& s, const json& interaction) {
  json j = json::parse(rtt::io::to_json(s));
  j["root"]["children"][2]["surfaces"][0]["interaction"] = interaction;
  return j.dump();
}

}  // namespace

TEST_CASE("schema version is 0.5.0", "[io][0.5]") {
  REQUIRE(rtt::model::kSchemaVersion == "0.5.0");
}

TEST_CASE("ideal lens: canonical form and round trip", "[io][0.5]") {
  SECTION("object at infinity: no object_distance in the file") {
    const System s = with_thin(IdealLens{Param(50.0), std::nullopt});
    CHECK(interaction_json(s) == json{{"type", "ideal_lens"}, {"focal_length", 50.0}});
    CHECK(rtt::io::parse_system(rtt::io::to_json(s)) == s);
  }
  SECTION("object_distance set, focal length variable with bounds") {
    Param f(50.0);
    f.variable = true;
    f.min = 20.0;
    f.max = 80.0;
    const System s = with_thin(IdealLens{f, Param(75.0)});
    const json i = interaction_json(s);
    CHECK(i["type"] == "ideal_lens");
    CHECK(i["focal_length"] ==
          json{{"value", 50.0}, {"variable", true}, {"min", 20.0}, {"max", 80.0}});
    CHECK(i["object_distance"] == 75.0);
    CHECK(rtt::io::parse_system(rtt::io::to_json(s)) == s);
  }
  SECTION("negative values: a diverging lens and a virtual object") {
    const System s = with_thin(IdealLens{Param(-50.0), Param(-30.0)});
    CHECK(rtt::io::parse_system(rtt::io::to_json(s)) == s);
  }
}

TEST_CASE("ideal cylinder lens: axis_deg is left out at 0", "[io][0.5]") {
  const System zero = with_thin(IdealCylinderLens{Param(50.0), 0.0, std::nullopt});
  CHECK(interaction_json(zero) == json{{"type", "ideal_cylinder_lens"}, {"focal_length", 50.0}});
  CHECK(rtt::io::parse_system(rtt::io::to_json(zero)) == zero);
  const System turned = with_thin(IdealCylinderLens{Param(50.0), 30.0, Param(75.0)});
  const json i = interaction_json(turned);
  CHECK(i["axis_deg"] == 30.0);
  CHECK(i["object_distance"] == 75.0);
  CHECK(rtt::io::parse_system(rtt::io::to_json(turned)) == turned);
}

TEST_CASE("ideal lens: a bound focal length and object distance (ADR 0029)", "[io][0.5]") {
  System s = with_thin(IdealLens{Param::bound("F"), Param::bound("OBJ")});
  rtt::model::ParameterRow f;
  f.name = "F";
  f.form = 50.0;
  rtt::model::ParameterRow obj;
  obj.name = "OBJ";
  obj.form = 75.0;
  s.parameters = {f, obj};
  const json i = interaction_json(s);
  CHECK(i["focal_length"] == json{{"param", "F"}});
  CHECK(i["object_distance"] == json{{"param", "OBJ"}});
  CHECK(rtt::io::parse_system(rtt::io::to_json(s)) == s);
}

TEST_CASE("edit form: axis_deg always, object_distance only when set (ADR 0024)", "[io][0.5]") {
  const System cyl = with_thin(IdealCylinderLens{Param(50.0), 0.0, std::nullopt});
  const json e = json::parse(rtt::io::to_edit_json(cyl));
  const json& i = e["root"]["children"][2]["surfaces"][0]["interaction"];
  CHECK(i["axis_deg"] == 0.0);
  CHECK_FALSE(i.contains("object_distance"));
  CHECK(i["focal_length"]["value"] == 50.0);
  CHECK(rtt::io::parse_system(rtt::io::to_edit_json(cyl)) == cyl);
  const System lens = with_thin(IdealLens{Param(50.0), Param(75.0)});
  const json l = json::parse(rtt::io::to_edit_json(lens));
  CHECK(l["root"]["children"][2]["surfaces"][0]["interaction"]["object_distance"]["value"] == 75.0);
  CHECK(rtt::io::parse_system(rtt::io::to_edit_json(lens)) == lens);
}

TEST_CASE("malformed ideal lenses are errors at their pointer", "[io][0.5]") {
  const System s = with_thin(IdealLens{Param(50.0), std::nullopt});
  const std::string at = "/root/children/2/surfaces/0/interaction";
  CHECK(error_pointer(with_interaction_json(s, json{{"type", "ideal_lens"}})) == at);
  CHECK(error_pointer(with_interaction_json(
            s, json{{"type", "ideal_lens"}, {"focal_length", "fifty"}})) == at + "/focal_length");
  CHECK(error_pointer(with_interaction_json(
            s, json{{"type", "ideal_lens"}, {"focal_length", 50.0}, {"axis_deg", 30.0}})) ==
        at + "/axis_deg");  // only the cylinder lens has an axis
  CHECK(error_pointer(with_interaction_json(
            s, json{{"type", "ideal_cylinder_lens"}, {"focal_length", 50.0}, {"axis_deg", "x"}})) ==
        at + "/axis_deg");
  CHECK(error_pointer(with_interaction_json(
            s, json{{"type", "ideal_lens"}, {"focal_length", 50.0}, {"object_distance", true}})) ==
        at + "/object_distance");
}

TEST_CASE("migration 0.4 -> 0.5 changes only the version", "[io][0.5][migration]") {
  json j = json::parse(rtt::io::to_json(rtt::model::test::make_singlet()));
  j["schema_version"] = "0.4.0";
  const System s = rtt::io::parse_system(j.dump());
  CHECK(s.schema_version == "0.5.0");
  json written = json::parse(rtt::io::to_json(s));
  written["schema_version"] = "0.4.0";
  CHECK(written == j);
}

TEST_CASE("the 0.5 interactions in an older file are errors", "[io][0.5][migration]") {
  const System s = with_thin(IdealLens{Param(50.0), std::nullopt});
  json j = json::parse(rtt::io::to_json(s));
  j["schema_version"] = "0.4.0";
  try {
    [[maybe_unused]] const System old = rtt::io::parse_system(j.dump());
    FAIL("an ideal lens in a 0.4 file was accepted");
  } catch (const rtt::io::ParseError& e) {
    CHECK(e.pointer() == "/root/children/2/surfaces/0/interaction/type");
    CHECK_THAT(std::string(e.what()), ContainsSubstring("0.5"));
  }
}
