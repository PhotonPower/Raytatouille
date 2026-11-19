// Schema 0.11 (#320, ADR 0050, point 1): aperture dimensions are Params; the migration
// 0.10 -> 0.11 changes only the version, and a Param object in an older file is a read error.

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <nlohmann/json.hpp>
#include <string>
#include <string_view>
#include <variant>
#include <vector>

#include "rtt/io/edit.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/model/parameters.hpp"
#include "rtt/model/system.hpp"
#include "rtt/model/validate.hpp"
#include "test_support.hpp"

using nlohmann::json;
using rtt::model::CircularAperture;
using rtt::model::EllipticalAperture;
using rtt::model::Param;
using rtt::model::ParameterRow;
using rtt::model::RectangularAperture;
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

/// The canonical file of `s` with the schema version `version`.
std::string with_version(const System& s, const std::string& version) {
  json j = json::parse(rtt::io::to_json(s));
  j["schema_version"] = version;
  return j.dump();
}

CircularAperture& stop_aperture(System& s) {
  return std::get<CircularAperture>(*rtt::model::test::element(s, 0).surfaces[0].aperture);
}

constexpr std::string_view kStop = "/root/children/0/surfaces/0/aperture";

}  // namespace

TEST_CASE("schema version is 0.11.0", "[io][0.11]") {
  REQUIRE(rtt::model::kSchemaVersion == "0.11.0");
}

TEST_CASE("aperture dimensions: plain Params stay numbers in the canonical form", "[io][0.11]") {
  const System s = rtt::model::test::make_singlet();
  const json j = json::parse(rtt::io::to_json(s));
  const json& a = j["root"]["children"][0]["surfaces"][0]["aperture"];
  CHECK(a["radius"] == 10.0);
  CHECK_FALSE(a.contains("inner_radius"));
  CHECK(rtt::io::parse_system(j.dump()) == s);
}

TEST_CASE("aperture dimensions: object and bound forms round trip", "[io][0.11]") {
  System s = rtt::model::test::make_singlet();
  Param r(10.0);
  r.variable = true;
  r.min = 5.0;
  r.max = 12.0;
  stop_aperture(s).radius = r;
  ParameterRow row;
  row.name = "R_IN";
  row.form = 2.0;
  s.parameters.push_back(row);
  stop_aperture(s).inner_radius = Param::bound("R_IN");
  auto& lens = rtt::model::test::element(s, 1);
  lens.surfaces[0].aperture = RectangularAperture{Param(9.0), Param::bound("R_IN")};
  lens.surfaces[1].aperture = EllipticalAperture{Param::bound("R_IN"), Param(8.0)};

  const json j = json::parse(rtt::io::to_json(s));
  const json& a = j["root"]["children"][0]["surfaces"][0]["aperture"];
  CHECK(a["radius"] ==
        json::parse(R"({"value": 10.0, "variable": true, "min": 5.0, "max": 12.0})"));
  CHECK(a["inner_radius"] == json::parse(R"({"param": "R_IN"})"));
  const System back = rtt::io::parse_system(j.dump());
  CHECK(back == s);
  CHECK_FALSE(rtt::model::has_errors(rtt::model::validate(back)));
}

TEST_CASE("aperture dimensions: the edit form writes objects; a patch sets the value",
          "[io][0.11]") {
  const System s = rtt::model::test::make_singlet();
  const json edit = json::parse(rtt::io::to_edit_json(s));
  CHECK(edit["root"]["children"][0]["surfaces"][0]["aperture"]["radius"]["value"] == 10.0);
  const json patch = json::array(
      {{{"op", "replace"}, {"path", std::string(kStop) + "/radius/value"}, {"value", 8.0}}});
  const rtt::io::PatchResult r =
      rtt::io::apply_patch(s, patch.dump(), rtt::io::EditCheck::NoNewErrors);
  System expected = s;
  stop_aperture(expected).radius = Param(8.0);
  CHECK(r.system == expected);
}

TEST_CASE("for_each_param visits the aperture after the shape and before the phases",
          "[io][0.11]") {
  const System s = rtt::model::test::make_singlet();
  std::vector<std::string> pointers;
  rtt::model::for_each_param(s,
                             [&](std::string_view p, const Param&) { pointers.emplace_back(p); });
  const std::string lens = "/root/children/1/surfaces/0";
  const auto pos = [&](const std::string& p) {
    return std::find(pointers.begin(), pointers.end(), p) - pointers.begin();
  };
  const auto end = static_cast<std::ptrdiff_t>(pointers.size());
  REQUIRE(pos(std::string(kStop) + "/radius") < end);
  REQUIRE(pos(std::string(kStop) + "/inner_radius") == pos(std::string(kStop) + "/radius") + 1);
  REQUIRE(pos(lens + "/aperture/radius") < end);
  CHECK(pos(lens + "/shape/base/conic") < pos(lens + "/aperture/radius"));
  CHECK(pos(lens + "/aperture/inner_radius") < pos("/root/children/1/surfaces/1/pose/position/0"));
}

TEST_CASE("validate: aperture dimensions are checked as Params; bound ones after resolving",
          "[io][0.11]") {
  System s = rtt::model::test::make_singlet();
  stop_aperture(s).radius = Param(-1.0);
  CHECK(rtt::model::test::has_error_at(rtt::model::validate(s), std::string(kStop) + "/radius"));
  // A bound radius with a bound inner radius: no cross check before resolving.
  ParameterRow r;
  r.name = "R";
  r.form = 10.0;
  s.parameters.push_back(r);
  stop_aperture(s).radius = Param::bound("R");
  stop_aperture(s).inner_radius = Param(20.0);
  CHECK_FALSE(rtt::model::has_errors(rtt::model::validate(s)));
  // The other way round: a bound inner radius with a plain radius.
  stop_aperture(s).radius = Param(10.0);
  stop_aperture(s).inner_radius = Param::bound("R");
  CHECK_FALSE(rtt::model::has_errors(rtt::model::validate(s)));
  Param bad = Param::bound("R");
  bad.variable = true;
  stop_aperture(s).inner_radius = bad;
  CHECK(rtt::model::test::has_error_at(rtt::model::validate(s),
                                       std::string(kStop) + "/inner_radius"));
}

TEST_CASE("migration 0.10 -> 0.11: content unchanged, the version becomes 0.11.0", "[io][0.11]") {
  const System s = rtt::model::test::make_singlet();
  const System read = rtt::io::parse_system(with_version(s, "0.10.0"));
  CHECK(read == s);
  CHECK(read.schema_version == "0.11.0");
  CHECK(json::parse(rtt::io::to_json(read))["schema_version"] == "0.11.0");
}

TEST_CASE("migration: a Param object at an aperture dimension before 0.11 is an error there",
          "[io][0.11]") {
  System s = rtt::model::test::make_singlet();
  Param r(10.0);
  r.variable = true;
  stop_aperture(s).radius = r;
  CHECK(error_pointer(with_version(s, "0.10.0")) == std::string(kStop) + "/radius");
  CHECK(error_pointer(with_version(s, "0.11.0")) == "<no error>");
  // Also in a file before 0.4 that goes through migrate_0_4 first.
  CHECK(error_pointer(with_version(s, "0.3.0")) == std::string(kStop) + "/radius");
  // In a nested assembly and at the other dimensions too.
  System nested = rtt::model::test::make_singlet();
  rtt::model::Assembly group;
  group.name = "group";
  group.children.push_back(nested.root.children[1]);
  nested.root.children[1] = rtt::model::Node{group};
  auto& lens = std::get<rtt::model::Element>(
      std::get<rtt::model::Assembly>(nested.root.children[1].value).children[0].value);
  Param b(4.0);
  b.variable = true;
  lens.surfaces[0].aperture = RectangularAperture{b, Param(3.0)};
  lens.surfaces[1].aperture = EllipticalAperture{Param(3.0), b};
  CHECK(error_pointer(with_version(nested, "0.10.0")) ==
        "/root/children/1/children/0/surfaces/0/aperture/half_width_x");
  lens.surfaces[0].aperture = RectangularAperture{Param(4.0), Param(3.0)};
  CHECK(error_pointer(with_version(nested, "0.10.0")) ==
        "/root/children/1/children/0/surfaces/1/aperture/semi_axis_y");
  // A plain number reads in every version.
  CHECK(error_pointer(with_version(rtt::model::test::make_singlet(), "0.10.0")) == "<no error>");
}

TEST_CASE("regular_polygon and polygon: canonical form, round trip, edit form (ADR 0050)",
          "[io][0.11]") {
  System s = rtt::model::test::make_singlet();
  auto& lens = rtt::model::test::element(s, 1);
  Param r(12.7);
  r.variable = true;
  lens.surfaces[0].aperture = rtt::model::RegularPolygonAperture{6, r, 0.0};
  lens.surfaces[1].aperture =
      rtt::model::PolygonAperture{{{-6.0, -6.0}, {6.0, -6.0}, {6.0, 6.0}, {-6.0, 6.0}}};
  const json j = json::parse(rtt::io::to_json(s));
  const json& a0 = j["root"]["children"][1]["surfaces"][0]["aperture"];
  CHECK(a0 == json::parse(
                  R"({"type": "regular_polygon", "sides": 6,
                      "circumradius": {"value": 12.7, "variable": true}})"));
  const json& a1 = j["root"]["children"][1]["surfaces"][1]["aperture"];
  CHECK(a1 == json::parse(R"({"type": "polygon",
                              "vertices": [[-6.0, -6.0], [6.0, -6.0], [6.0, 6.0], [-6.0, 6.0]]})"));
  CHECK(rtt::io::parse_system(j.dump()) == s);
  CHECK_FALSE(rtt::model::has_errors(rtt::model::validate(s)));
  // A rotation is written when it is not 0; the edit form writes it always.
  lens.surfaces[0].aperture = rtt::model::RegularPolygonAperture{5, Param(12.7), 18.0};
  CHECK(json::parse(rtt::io::to_json(
            s))["root"]["children"][1]["surfaces"][0]["aperture"]["rotation_deg"] == 18.0);
  const json edit = json::parse(rtt::io::to_edit_json(rtt::model::test::make_singlet()));
  System hex = rtt::model::test::make_singlet();
  rtt::model::test::element(hex, 1).surfaces[0].aperture =
      rtt::model::RegularPolygonAperture{6, Param(12.7), 0.0};
  const json edit_hex = json::parse(rtt::io::to_edit_json(hex));
  CHECK(edit_hex["root"]["children"][1]["surfaces"][0]["aperture"]["rotation_deg"] == 0.0);
  CHECK(edit["root"]["children"][1]["surfaces"][0]["aperture"]["type"] == "circular");
  // for_each_param visits the circumradius, not the vertices.
  std::vector<std::string> pointers;
  rtt::model::for_each_param(hex,
                             [&](std::string_view p, const Param&) { pointers.emplace_back(p); });
  CHECK(std::find(pointers.begin(), pointers.end(),
                  "/root/children/1/surfaces/0/aperture/circumradius") != pointers.end());
}

TEST_CASE("regular_polygon and polygon: read errors (ADR 0050)", "[io][0.11]") {
  const auto with_aperture = [](const std::string& aperture, const std::string& version) {
    json j = json::parse(rtt::io::to_json(rtt::model::test::make_singlet()));
    j["schema_version"] = version;
    j["root"]["children"][1]["surfaces"][0]["aperture"] = json::parse(aperture);
    return j.dump();
  };
  const std::string lens = "/root/children/1/surfaces/0/aperture";
  // The new types before 0.11 are an error at the type.
  CHECK(error_pointer(with_aperture(R"({"type": "regular_polygon", "sides": 6,
                                        "circumradius": 5.0})",
                                    "0.10.0")) == lens + "/type");
  CHECK(error_pointer(with_aperture(R"({"type": "polygon", "vertices": [[0, 0], [1, 0], [0, 1]]})",
                                    "0.9.0")) == lens + "/type");
  // Structural errors in 0.11.
  CHECK(error_pointer(with_aperture(R"({"type": "regular_polygon", "circumradius": 5.0})",
                                    "0.11.0")) == lens);
  CHECK(error_pointer(with_aperture(R"({"type": "regular_polygon", "sides": -1,
                                        "circumradius": 5.0})",
                                    "0.11.0")) == lens + "/sides");
  CHECK(error_pointer(with_aperture(R"({"type": "polygon", "vertices": [[0, 0, 1]]})", "0.11.0")) ==
        lens + "/vertices/0");
  CHECK(error_pointer(with_aperture(R"({"type": "polygon", "vertices": [[0, "a"]]})", "0.11.0")) ==
        lens + "/vertices/0/1");
  CHECK(error_pointer(with_aperture(R"({"type": "hexagonal", "circumradius": 5.0})", "0.11.0")) ==
        lens + "/type");
}

TEST_CASE("obscurations: canonical form, round trip, edit form, for_each_param (ADR 0050)",
          "[io][0.11]") {
  System s = rtt::model::test::make_singlet();
  auto& lens = rtt::model::test::element(s, 1);
  Param dx(1.5);
  dx.variable = true;
  lens.surfaces[0].obscurations = {
      rtt::model::CircularObscuration{Param(2.0), Param(0.0), Param(0.0)},
      rtt::model::RectangularObscuration{Param(12.0), Param(0.5), dx, Param(-0.25), 45.0},
      rtt::model::EllipticalObscuration{Param(1.0), Param(0.5), Param(0.0), Param(0.0), 0.0}};
  const json j = json::parse(rtt::io::to_json(s));
  const json& o = j["root"]["children"][1]["surfaces"][0]["obscurations"];
  // Defaults (decenter 0, rotation 0) are omitted in the canonical form.
  CHECK(o == json::parse(R"([
      {"type": "circular", "radius": 2.0},
      {"type": "rectangular", "half_width_x": 12.0, "half_width_y": 0.5,
       "decenter_x": {"value": 1.5, "variable": true}, "decenter_y": -0.25,
       "rotation_deg": 45.0},
      {"type": "elliptical", "semi_axis_x": 1.0, "semi_axis_y": 0.5}])"));
  CHECK(rtt::io::parse_system(j.dump()) == s);
  CHECK_FALSE(rtt::model::has_errors(rtt::model::validate(s)));
  // Without obscurations the key is missing in the canonical form, present in the edit form.
  CHECK_FALSE(j["root"]["children"][1]["surfaces"][1].contains("obscurations"));
  const json edit = json::parse(rtt::io::to_edit_json(s));
  CHECK(edit["root"]["children"][1]["surfaces"][1]["obscurations"] == json::array());
  const json& e0 = edit["root"]["children"][1]["surfaces"][0]["obscurations"][0];
  CHECK(e0["decenter_x"] == json::parse(R"({"value": 0.0, "variable": false})"));
  // A patch adds the first obscuration to a surface.
  const json patch = json::parse(R"([{"op": "add",
      "path": "/root/children/1/surfaces/1/obscurations/0",
      "value": {"type": "circular", "radius": 1.0}}])");
  const rtt::io::PatchResult r =
      rtt::io::apply_patch(s, patch.dump(), rtt::io::EditCheck::NoNewErrors);
  System added = r.system;
  REQUIRE(rtt::model::test::element(added, 1).surfaces[1].obscurations.size() == 1);
  // for_each_param: after the aperture, in list order and field order, before the next surface.
  std::vector<std::string> pointers;
  rtt::model::for_each_param(s,
                             [&](std::string_view p, const Param&) { pointers.emplace_back(p); });
  const auto pos = [&](const std::string& p) {
    return std::find(pointers.begin(), pointers.end(), p) - pointers.begin();
  };
  const std::string ob = "/root/children/1/surfaces/0/obscurations/";
  const auto end = static_cast<std::ptrdiff_t>(pointers.size());
  REQUIRE(pos(ob + "2/decenter_y") < end);
  CHECK(pos("/root/children/1/surfaces/0/aperture/inner_radius") < pos(ob + "0/radius"));
  CHECK(pos(ob + "0/radius") + 1 == pos(ob + "0/decenter_x"));
  CHECK(pos(ob + "0/decenter_y") + 1 == pos(ob + "1/half_width_x"));
  CHECK(pos(ob + "1/half_width_y") + 1 == pos(ob + "1/decenter_x"));
  CHECK(pos(ob + "2/decenter_y") < pos("/root/children/1/surfaces/1/pose/position/0"));
}

TEST_CASE("obscurations: read errors and the version before 0.11 (ADR 0050)", "[io][0.11]") {
  const auto with_obscurations = [](const std::string& list, const std::string& version) {
    json j = json::parse(rtt::io::to_json(rtt::model::test::make_singlet()));
    j["schema_version"] = version;
    j["root"]["children"][1]["surfaces"][0]["obscurations"] = json::parse(list);
    return j.dump();
  };
  const std::string at = "/root/children/1/surfaces/0/obscurations";
  // The key before 0.11 is an error at the key.
  CHECK(error_pointer(with_obscurations(R"([{"type": "circular", "radius": 1.0}])", "0.10.0")) ==
        at);
  CHECK(error_pointer(with_obscurations(R"([])", "0.9.0")) == at);
  // Structural errors in 0.11.
  CHECK(error_pointer(with_obscurations(R"({"type": "circular"})", "0.11.0")) == at);
  CHECK(error_pointer(with_obscurations(R"([{"type": "circular"}])", "0.11.0")) == at + "/0");
  CHECK(error_pointer(with_obscurations(R"([{"type": "circular", "radius": 1.0,
                                             "rotation_deg": 5.0}])",
                                        "0.11.0")) == at + "/0/rotation_deg");
  CHECK(error_pointer(with_obscurations(R"([{"type": "spider", "radius": 1.0}])", "0.11.0")) ==
        at + "/0/type");
  CHECK(error_pointer(with_obscurations(R"([{"type": "rectangular", "half_width_x": 1.0,
                                             "half_width_y": 1.0, "rotation_deg": "a"}])",
                                        "0.11.0")) == at + "/0/rotation_deg");
  CHECK(error_pointer(with_obscurations(R"([{"type": "circular", "radius": 1.0}])", "0.11.0")) ==
        "<no error>");
}

TEST_CASE("zernike_fringe_sag: canonical form, round trip, edit form, version (ADR 0051)",
          "[io][0.11]") {
  System s = rtt::model::test::make_singlet();
  auto& lens = rtt::model::test::element(s, 1);
  Param a(1e-4);
  a.variable = true;
  lens.surfaces[0].shape.terms = {
      rtt::model::ZernikeSag{Param(12.0), {Param(0.0), Param(0.0), Param(0.0), Param(2e-4)}},
      rtt::model::ZernikeFringeSag{Param(12.0), {Param(0.0), Param(0.0), Param(0.0), a}}};
  const json j = json::parse(rtt::io::to_json(s));
  const json& t = j["root"]["children"][1]["surfaces"][0]["shape"]["terms"];
  CHECK(t[1] == json::parse(R"({"type": "zernike_fringe_sag", "normalization_radius": 12.0,
      "coefficients": [0.0, 0.0, 0.0, {"value": 0.0001, "variable": true}]})"));
  CHECK(t[0]["type"] == "zernike_sag");
  CHECK(rtt::io::parse_system(j.dump()) == s);
  CHECK_FALSE(rtt::model::has_errors(rtt::model::validate(s)));
  // for_each_param visits the Fringe coefficients after the Noll term.
  std::vector<std::string> pointers;
  rtt::model::for_each_param(s,
                             [&](std::string_view p, const Param&) { pointers.emplace_back(p); });
  const std::string term = "/root/children/1/surfaces/0/shape/terms/1";
  CHECK(std::find(pointers.begin(), pointers.end(), term + "/coefficients/3") != pointers.end());
  // Before 0.11 the new term type is a read error at its type; zernike_sag stays readable.
  json old = j;
  old["schema_version"] = "0.10.0";
  CHECK(error_pointer(old.dump()) == term + "/type");
  old["root"]["children"][1]["surfaces"][0]["shape"]["terms"].erase(1);
  CHECK(error_pointer(old.dump()) == "<no error>");
  // Too many Fringe coefficients are a validation error, not a read error.
  lens.surfaces[0].shape.terms = {
      rtt::model::ZernikeFringeSag{Param(12.0), std::vector<Param>(38, Param(0.0))}};
  CHECK(rtt::model::has_errors(rtt::model::validate(s)));
}
