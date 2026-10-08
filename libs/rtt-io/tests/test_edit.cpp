// Changing a system with JSON Patch on its edit form, with inverse patches (ADR 0024, #82).

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <limits>
#include <map>
#include <nlohmann/json.hpp>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <variant>
#include <vector>

#include "rtt/diagnostics/codes.hpp"
#include "rtt/io/edit.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/model/validate.hpp"
#include "test_support.hpp"

namespace fs = std::filesystem;
using nlohmann::json;
using rtt::io::apply_patch;
using rtt::io::EditCheck;
using rtt::io::EditError;
using rtt::io::PatchResult;
using rtt::model::System;

namespace {

std::vector<fs::path> reference_files() {
  std::vector<fs::path> files;
  for (const auto& entry : fs::recursive_directory_iterator(RTT_REFERENCE_DIR)) {
    const std::string name = entry.path().filename().string();
    if (entry.is_regular_file() && name.ends_with(".rtt.json")) files.push_back(entry.path());
  }
  std::sort(files.begin(), files.end());
  return files;
}

std::string read_text(const fs::path& p) {
  std::ifstream in(p, std::ios::binary);
  std::ostringstream s;
  s << in.rdbuf();
  return s.str();
}

System load(const std::string& relative) {
  return rtt::io::load_system(fs::path(RTT_REFERENCE_DIR) / relative);
}

/// Applies the patch, undoes it with the inverse (StructureOnly) and checks that the canonical
/// JSON is bit for bit the original (acceptance criterion of #82). `changes`: whether the patch
/// changes the model; if it does, the result must differ from the input, so that an empty or
/// wrong inverse cannot pass.
PatchResult applied_and_undone(const System& s, const std::string& patch, bool changes = true) {
  INFO("patch " << patch);
  const PatchResult r = apply_patch(s, patch);
  INFO("inverse " << r.inverse);
  CHECK((rtt::io::to_json(r.system) != rtt::io::to_json(s)) == changes);
  CHECK(json::parse(r.inverse) != json::array());
  const PatchResult back = apply_patch(r.system, r.inverse, EditCheck::StructureOnly);
  REQUIRE(rtt::io::to_json(back.system) == rtt::io::to_json(s));
  return r;
}

bool identity_pose(const System& s, std::size_t child) {
  return std::visit([](const auto& n) { return n.pose.is_identity(); },
                    s.root.children[child].value);
}

// ---- Keys of the edit form (ADR 0024 point 3): every value written, every Param an object ----

using Keys = std::set<std::string>;

/// `o` has exactly the keys `keys` plus those of `optional` it has.
void has_keys(const json& o, Keys keys, const Keys& optional, const std::string& where) {
  INFO(where);
  REQUIRE(o.is_object());
  Keys have;
  for (const auto& item : o.items()) have.insert(item.key());
  for (const std::string& k : optional) {
    if (have.contains(k)) keys.insert(k);
  }
  CHECK(have == keys);
}

void param(const json& p, const std::string& where) {
  has_keys(p, {"value", "variable"}, {"pickup"}, where);
}

void params(const json& list, const std::string& where) {
  REQUIRE(list.is_array());
  for (std::size_t i = 0; i < list.size(); ++i) param(list[i], where + "/" + std::to_string(i));
}

void pose(const json& p, const std::string& where) {
  has_keys(p, {"position", "rotation_deg", "pivot"}, {}, where);
  params(p["position"], where + "/position");
  params(p["rotation_deg"], where + "/rotation_deg");
  CHECK(p["pivot"].size() == 3);
}

/// A variant object: its keys by "type". Where `params_inside` is set (shapes and phases; not
/// apertures, whose values are plain numbers), radius, conic, normalization_radius and
/// lines_per_mm hold a Param and coefficients a list of Params.
void variant(const json& o,
             const std::map<std::string, Keys>& by_type,
             const std::string& where,
             bool params_inside = true) {
  INFO(where);
  REQUIRE(o.contains("type"));
  const auto it = by_type.find(o["type"].get<std::string>());
  REQUIRE(it != by_type.end());
  has_keys(o, it->second, {}, where);
  if (!params_inside) return;
  for (const char* k : {"radius", "conic", "normalization_radius", "lines_per_mm"}) {
    if (o.contains(k)) param(o[k], where + "/" + k);
  }
  if (o.contains("coefficients")) params(o["coefficients"], where + "/coefficients");
}

void surface_keys(const json& s, const std::string& at) {
  has_keys(s, {"id", "pose", "shape", "phases", "interaction"}, {"aperture"}, at);
  pose(s["pose"], at + "/pose");
  has_keys(s["shape"], {"base", "terms"}, {}, at + "/shape");
  variant(s["shape"]["base"],
          {{"plane", {"type"}},
           {"conic", {"type", "radius", "conic"}},
           {"even_asphere", {"type", "radius", "conic", "coefficients"}}},
          at + "/shape/base");
  for (const json& t : s["shape"]["terms"]) {
    variant(t, {{"zernike_sag", {"type", "normalization_radius", "coefficients"}}},
            at + "/shape/terms");
  }
  if (s.contains("aperture")) {
    variant(s["aperture"],
            {{"circular", {"type", "radius", "inner_radius"}},
             {"rectangular", {"type", "half_width_x", "half_width_y"}},
             {"elliptical", {"type", "semi_axis_x", "semi_axis_y"}}},
            at + "/aperture", false);
  }
  for (const json& ph : s["phases"]) {
    variant(ph,
            {{"linear_grating", {"type", "lines_per_mm", "orientation_deg"}},
             {"radial_phase", {"type", "normalization_radius", "coefficients"}}},
            at + "/phases");
  }
  variant(s["interaction"],
          {{"fresnel", {"type"}},
           {"ideal_mirror", {"type"}},
           {"ideal_anti_reflection", {"type"}},
           {"absorber", {"type"}},
           {"ideal_beam_splitter", {"type", "reflectance_s", "reflectance_p"}},
           {"coating", {"type", "name"}},
           {"ideal_polarizer", {"type", "transmission_axis", "extinction_ratio"}},
           {"ideal_retarder", {"type", "fast_axis", "retardance_waves"}}},
          at + "/interaction");
}

void node_keys(const json& n, const std::string& where) {
  if (n["type"] == "assembly") {
    has_keys(n, {"type", "name", "pose", "children"}, {}, where);
    for (std::size_t i = 0; i < n["children"].size(); ++i) {
      node_keys(n["children"][i], where + "/children/" + std::to_string(i));
    }
  } else {
    has_keys(n, {"type", "name", "pose", "surfaces"}, {"material", "optic_axis"}, where);
    for (std::size_t i = 0; i < n["surfaces"].size(); ++i) {
      surface_keys(n["surfaces"][i], where + "/surfaces/" + std::to_string(i));
    }
  }
  pose(n["pose"], where + "/pose");
}

void edit_form_keys(const json& e) {
  has_keys(e,
           {"schema_version", "name", "units", "environment", "object", "wavelengths", "aperture",
            "fields", "root", "paths"},
           {}, "");
  has_keys(e["environment"], {"temperature_c", "pressure_atm", "medium"}, {}, "/environment");
  has_keys(e["object"], {"at_infinity", "distance"}, {}, "/object");
  param(e["object"]["distance"], "/object/distance");
  for (const json& w : e["wavelengths"]) {
    has_keys(w, {"um", "weight", "reference"}, {}, "/wavelengths");
  }
  has_keys(e["aperture"], {"type", "value"}, {}, "/aperture");
  param(e["aperture"]["value"], "/aperture/value");
  has_keys(e["fields"], {"type", "points"}, {}, "/fields");
  for (const json& f : e["fields"]["points"]) has_keys(f, {"x", "y", "weight"}, {}, "/points");
  node_keys(e["root"], "/root");
  for (const json& p : e["paths"]) {
    has_keys(p, {"name", "events"}, {}, "/paths");
    if (p["events"].is_array()) {
      for (const json& ev : p["events"]) has_keys(ev, {"surface", "kind", "order"}, {}, "/events");
    } else {
      CHECK(p["events"] == "auto");
    }
  }
}

EditError edit_error(const System& s,
                     const std::string& patch,
                     EditCheck check = EditCheck::NoNewErrors) {
  const System before = s;
  try {
    (void)apply_patch(s, patch, check);
  } catch (const EditError& e) {
    // The input is never changed (a NaN temperature compares unequal to itself).
    if (!std::isnan(s.environment.temperature_c)) REQUIRE(s == before);
    return e;
  }
  FAIL("no EditError for " << patch);
  return EditError("", "", std::nullopt, {}, "");
}

/// Every JSON object below `j` (depth first) with its pointer.
void objects(const json& j,
             const std::string& ptr,
             std::vector<std::pair<std::string, json>>& out) {
  if (j.is_object()) {
    out.emplace_back(ptr, j);
    for (const auto& item : j.items()) objects(item.value(), ptr + "/" + item.key(), out);
  } else if (j.is_array()) {
    for (std::size_t i = 0; i < j.size(); ++i) objects(j[i], ptr + "/" + std::to_string(i), out);
  }
}

bool is_param_object(const json& j) {
  return j.is_object() && j.contains("value") && j["value"].is_number() && j.contains("variable") &&
         j["variable"].is_boolean();
}

}  // namespace

TEST_CASE("the edit form reads back to the same system for every reference system", "[io][edit]") {
  for (const fs::path& file : reference_files()) {
    INFO(file.string());
    const std::string text = read_text(file);
    const System s = rtt::io::parse_system(text);
    const std::string edit = rtt::io::to_edit_json(s);
    const System again = rtt::io::parse_system(edit);
    REQUIRE(again == s);
    REQUIRE(rtt::io::to_json(again) == text);
    // Every pointer into the canonical file exists in the edit form.
    std::vector<std::pair<std::string, json>> canonical;
    objects(json::parse(text), "", canonical);
    const json e = json::parse(edit);
    for (const auto& [ptr, obj] : canonical) {
      for (const auto& item : obj.items()) {
        INFO(ptr << "/" << item.key());
        CHECK(e.contains(json::json_pointer(ptr + "/" + item.key())));
      }
    }
  }
}

TEST_CASE("the edit form has every key of every object in every reference system", "[io][edit]") {
  for (const fs::path& file : reference_files()) {
    INFO(file.string());
    edit_form_keys(json::parse(rtt::io::to_edit_json(rtt::io::load_system(file))));
  }
}

TEST_CASE("the edit form writes every value", "[io][edit]") {
  const json e = json::parse(rtt::io::to_edit_json(load("m0/singlet.rtt.json")));
  for (const char* key :
       {"name", "environment", "object", "wavelengths", "aperture", "fields", "root", "paths"}) {
    CHECK(e.contains(key));
  }
  CHECK(e["environment"].size() == 3);
  CHECK(e["object"].contains("distance"));
  CHECK(is_param_object(e["object"]["distance"]));
  CHECK(is_param_object(e["aperture"]["value"]));
  CHECK(e["fields"].contains("type"));
  for (const json& w : e["wavelengths"]) CHECK((w.contains("weight") && w.contains("reference")));
  for (const json& p : e["fields"]["points"]) CHECK(p.size() == 3);
  CHECK(e["root"].contains("pose"));
  for (const json& node : e["root"]["children"]) {
    CHECK(node.contains("pose"));
    for (const json& p : node["pose"]["position"]) CHECK(is_param_object(p));
    CHECK(node["pose"]["pivot"].size() == 3);
    for (const json& s : node["surfaces"]) {
      for (const char* key : {"pose", "shape", "phases", "interaction"}) CHECK(s.contains(key));
      CHECK(s["shape"].contains("base"));
      CHECK(s["shape"].contains("terms"));
    }
  }
  const json lens = e["root"]["children"][1];
  REQUIRE(lens["name"] == "L1");
  CHECK(lens["surfaces"][0]["shape"]["base"].contains("conic"));
  CHECK(is_param_object(lens["surfaces"][0]["shape"]["base"]["radius"]));
  const json events = e["paths"][0]["events"];
  CHECK((events == "auto" || (events[0].contains("kind") && events[0].contains("order"))));
  // A pickup stays; an unset optional member is missing.
  const json tour = json::parse(rtt::io::to_edit_json(load("m0/feature_tour.rtt.json")));
  const json z = tour["root"]["children"][2]["surfaces"][1]["pose"]["position"][2];
  CHECK(z == json::parse(R"({"value": 6.0, "variable": false, "pickup": "2 * 3"})"));
  CHECK(!tour["root"]["children"][5]["surfaces"][0].contains("aperture"));
  CHECK(!tour["root"]["children"][5].contains("material"));
}

TEST_CASE("the acceptance patches are undone bit for bit on every reference system", "[io][edit]") {
  for (const fs::path& file : reference_files()) {
    INFO(file.string());
    const System s = rtt::io::load_system(file);
    const json e = json::parse(rtt::io::to_edit_json(s));
    const std::string value = json(e["aperture"]["value"]["value"].get<double>() + 1.0).dump();
    // A Param object replaced by a number (drops variable and pickup).
    applied_and_undone(s,
                       R"([{"op": "replace", "path": "/aperture/value", "value": )" + value + "}]");
    // A move into an ancestor of its source.
    applied_and_undone(
        s, R"([{"op": "replace", "path": "/aperture/value/value", "value": )" + value + R"(},
                              {"op": "move", "from": "/aperture/value/value",
                               "path": "/aperture/value"}])");
    // A member with its default value removed (the reader restores the default: no change).
    applied_and_undone(s, R"([{"op": "remove", "path": "/root/pose"}])",
                       !s.root.pose.is_identity());
    // Append with "-", and a move to "-" of the same array.
    applied_and_undone(s, R"([{"op": "add", "path": "/fields/points/-",
                              "value": {"x": 0.0, "y": 0.5}}])");
    applied_and_undone(s, R"([{"op": "move", "from": "/wavelengths/0", "path": "/wavelengths/-"}])",
                       s.wavelengths.size() >= 2);
    // copy, and several operations in one patch.
    applied_and_undone(
        s, R"([{"op": "copy", "from": "/fields/points/0", "path": "/fields/points/0"}])");
    applied_and_undone(s, R"([{"op": "replace", "path": "/name", "value": "changed"},
                              {"op": "add", "path": "/fields/points/0", "value": {"y": 0.25}},
                              {"op": "remove", "path": "/root/pose"},
                              {"op": "test", "path": "/name", "value": "changed"}])");
    // A move onto an existing member.
    if (s.root.children.size() >= 2) {
      applied_and_undone(s, R"([{"op": "move", "from": "/root/children/0/pose",
                                "path": "/root/children/1/pose"}])",
                         !(identity_pose(s, 0) && identity_pose(s, 1)));
    }
    // A move onto an existing member whose way back shifts array indices (second review of
    // ADR 0024): after removing point 0, point 1 is the former point 2.
    System three = s;
    while (three.fields.points.size() < 3) {
      three =
          apply_patch(three, R"([{"op": "add", "path": "/fields/points/-", "value": {"y": 0.1}}])")
              .system;
    }
    applied_and_undone(three, R"([{"op": "replace", "path": "/fields/points/0", "value": 0.5},
                                  {"op": "move", "from": "/fields/points/0",
                                   "path": "/fields/points/1/y"}])");
  }
}

TEST_CASE("a repair is accepted and undone; a new error is rejected", "[io][edit]") {
  const System singlet = load("m0/singlet.rtt.json");
  // A broken state, reached with StructureOnly (as an undo could).
  const System broken =
      apply_patch(singlet, R"([{"op": "replace", "path": "/wavelengths/0/um", "value": -1.0}])",
                  EditCheck::StructureOnly)
          .system;
  REQUIRE(rtt::model::has_errors(rtt::model::validate(broken)));
  // The repair is accepted; its undo restores the broken state bit for bit.
  const PatchResult repaired = applied_and_undone(
      broken, R"([{"op": "replace", "path": "/wavelengths/0/um", "value": 0.4861}])");
  CHECK_FALSE(rtt::model::has_errors(repaired.diagnostics));
  // Its inverse adds the error back, so it needs StructureOnly (why undo does not use the check).
  CHECK_THROWS_AS(apply_patch(repaired.system, repaired.inverse), EditError);
  // Another error of the same code instead of the repaired one is accepted (ADR 0024 point 4,
  // limit): the number of errors per code does not grow.
  CHECK_NOTHROW(apply_patch(broken, R"([
      {"op": "replace", "path": "/wavelengths/0/um", "value": 0.5},
      {"op": "replace", "path": "/wavelengths/1/um", "value": -2.0}])"));
  // A new error is rejected with the code and place of validate().
  const std::string duplicate = R"([{"op": "replace", "path": "/root/children/1/surfaces/1/id",
                                     "value": "L1.S1"}])";
  const EditError e = edit_error(singlet, duplicate);
  CHECK(e.code() == "surface.id_duplicate");
  CHECK(!e.op_index().has_value());
  const System dup = apply_patch(singlet, duplicate, EditCheck::StructureOnly).system;
  std::vector<rtt::model::Diagnostic> expected;
  for (const auto& d : rtt::model::validate(dup)) {
    if (d.code == "surface.id_duplicate") expected.push_back(d);
  }
  REQUIRE(!expected.empty());
  CHECK(e.location() == expected.front().location);
  CHECK(e.diagnostics().size() == expected.size());
}

TEST_CASE("an empty System is editable", "[io][edit]") {
  const System empty{};
  REQUIRE(rtt::model::has_errors(rtt::model::validate(empty)));
  const PatchResult r = applied_and_undone(
      empty, R"([{"op": "add", "path": "/wavelengths/-", "value": {"um": 0.55, "reference": true}},
                 {"op": "replace", "path": "/root/name", "value": "system"}])");
  CHECK(r.system.wavelengths.size() == 1);
  CHECK(r.system.root.name == "system");
}

TEST_CASE("every edit error has its code, place and operation", "[io][edit]") {
  const System s = load("m0/singlet.rtt.json");
  std::set<std::string> seen;
  const auto check = [&](const EditError& e, const std::string& code, const std::string& location,
                         std::optional<std::size_t> op) {
    INFO(e.what());
    CHECK(e.code() == code);
    CHECK(e.location() == location);
    CHECK(e.op_index() == op);
    REQUIRE(rtt::diagnostics::find_code(e.code()) != nullptr);
    seen.insert(e.code());
  };
  check(edit_error(s, R"({"op": "remove", "path": "/name"})"), "edit.patch_invalid", "",
        std::nullopt);
  check(edit_error(s, R"([{"op": "add", "path": "/name", "value": 1, "op": "remove"}])"),
        "edit.patch_invalid", "", std::nullopt);  // RFC 6902 A.13: duplicate key
  check(
      edit_error(s, R"([{"op": "test", "path": "/name"}, {"op": "frobnicate", "path": "/name"}])"),
      "edit.patch_invalid", "/name", 0);  // test without value
  check(edit_error(s, R"([{"op": "remove", "path": "/root/children/9"}])"), "edit.path_not_found",
        "/root/children/9", 0);
  check(edit_error(s, R"([{"op": "replace", "path": "/name", "value": "x"},
                          {"op": "replace", "path": "/schema_version", "value": "0.1.0"}])"),
        "edit.read_only", "/schema_version", 1);
  check(edit_error(s, R"([{"op": "remove", "path": "/units/length"}])"), "edit.read_only",
        "/units/length", 0);
  check(edit_error(s, R"([{"op": "move", "from": "/units", "path": "/name"}])"), "edit.read_only",
        "/units", 0);
  check(edit_error(s, R"([{"op": "replace", "path": "", "value": {}}])"), "edit.read_only", "", 0);
  check(edit_error(s, R"([{"op": "test", "path": "/name", "value": "other"}])"), "edit.test_failed",
        "/name", 0);
  // The strict reader rejects the result: a wrong type, a key that does not fit the variant.
  check(edit_error(s, R"([{"op": "replace", "path": "/name", "value": 5}])"), "edit.invalid_value",
        "/name", std::nullopt);
  check(edit_error(s, R"([{"op": "replace", "path": "/root/children/1/surfaces/0/shape/base/type",
                           "value": "plane"}])"),
        "edit.invalid_value", "/root/children/1/surfaces/0/shape/base/conic",
        std::nullopt);  // the first key a plane does not have, in key order
  // A key with "~" and "/" is masked in the reader's pointer (RFC 6901).
  check(edit_error(s, R"([{"op": "add", "path": "/environment/a~1b~0", "value": 1}])"),
        "edit.invalid_value", "/environment/a~1b~0", std::nullopt);
  // The input system cannot be written.
  System nan = s;
  nan.environment.temperature_c = std::numeric_limits<double>::quiet_NaN();
  check(edit_error(nan, R"([{"op": "test", "path": "/name", "value": "x"}])"),
        "edit.base_not_representable", "", std::nullopt);
  // Reading /schema_version with test is allowed.
  // (Built outside the macro: MSVC's preprocessor rejected this raw string followed by a string
  // with an escaped quote inside a macro argument.)
  const std::string read_version =
      R"([{"op": "test", "path": "/schema_version", "value": ")" + s.schema_version + "\"}]";
  CHECK_NOTHROW(apply_patch(s, read_version));

  // Every registered code of the producer "edit" has a case here.
  for (const auto& info : rtt::diagnostics::kCodes) {
    if (info.producer != "edit") continue;
    INFO(info.code);
    CHECK(seen.contains(std::string(info.code)));
    CHECK(info.severity == rtt::diagnostics::Severity::Error);
  }
}

TEST_CASE("an integer outside int in a patch is an invalid value (#35)", "[io][edit]") {
  // The strict reader rejects it (read_int); before #35 E it became 0 without an error.
  const System s = load("m2/cooke_triplet.rtt.json");
  const EditError e = edit_error(
      s, R"([{"op": "replace", "path": "/paths/0/events/0/order", "value": 4294967296}])");
  CHECK(e.code() == "edit.invalid_value");
  CHECK(e.location() == "/paths/0/events/0/order");
  CHECK_THAT(std::string(e.what()), Catch::Matchers::ContainsSubstring("out of range"));
}

TEST_CASE("a failing operation leaves the system unchanged", "[io][edit]") {
  const System s = load("m2/cooke_triplet.rtt.json");
  const EditError e = edit_error(s, R"([{"op": "replace", "path": "/name", "value": "x"},
                                        {"op": "remove", "path": "/paths/7"}])");
  CHECK(e.op_index() == 1);
  CHECK(e.code() == "edit.path_not_found");
}

TEST_CASE("cost of a patch on large systems", "[io][edit][.benchmark]") {
  const System lens = load("m2/cooke_triplet.rtt.json");
  for (const std::size_t copies : {500U, 2500U}) {
    System big = lens;
    const auto& first = std::get<rtt::model::Element>(lens.root.children[0].value);
    big.root.children.clear();
    big.paths = {rtt::model::Path{"main", true, {}}};
    for (std::size_t i = 0; i < copies; ++i) {
      rtt::model::Element e = first;
      e.name = "L" + std::to_string(i);
      for (std::size_t k = 0; k < e.surfaces.size(); ++k) {
        e.surfaces[k].id = rtt::model::SurfaceId(e.name + ".S" + std::to_string(k));
      }
      big.root.children.push_back(rtt::model::Node{e});
    }
    const auto t0 = std::chrono::steady_clock::now();
    const std::string edit = rtt::io::to_edit_json(big);
    const auto t1 = std::chrono::steady_clock::now();
    const PatchResult r =
        apply_patch(big, R"([{"op": "replace", "path": "/root/children/0/name", "value": "X"}])");
    const auto t2 = std::chrono::steady_clock::now();
    const auto ms = [](auto a, auto b) {
      return std::chrono::duration<double, std::milli>(b - a).count();
    };
    WARN(2 * copies << " surfaces: edit form " << edit.size() / 1024 << " KiB, to_edit_json "
                    << ms(t0, t1) << " ms, apply_patch " << ms(t1, t2) << " ms");
    CHECK(r.system.root.children.size() == copies);
  }
}
