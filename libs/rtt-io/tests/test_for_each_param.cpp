// Guard of rtt::model::for_each_param (#164, ADR 0029, point 4): for every reference system it
// visits exactly the Param objects of the edit form (rtt::io::to_edit_json), in the order in
// which the edit form writes them. A new Param field that for_each_param misses, or a writer
// change of the order, turns this red.

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <filesystem>
#include <nlohmann/json.hpp>
#include <regex>
#include <set>
#include <string>
#include <string_view>
#include <vector>

#include "rtt/io/edit.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/model/parameters.hpp"

namespace fs = std::filesystem;
using ojson = nlohmann::ordered_json;

namespace {

std::vector<fs::path> reference_files() {
  std::vector<fs::path> files;
  for (const auto& entry : fs::recursive_directory_iterator(RTT_REFERENCE_DIR)) {
    if (entry.is_regular_file() && entry.path().string().ends_with(".rtt.json")) {
      files.push_back(entry.path());
    }
  }
  return files;
}

/// A Param object of the edit form: {"param": …} alone, or "value" with keys only from
/// {value, variable, min, max}.
bool param_object(const ojson& j) {
  if (!j.is_object()) return false;
  if (j.contains("param")) return j.size() == 1;
  if (!j.contains("value")) return false;
  for (const auto& [key, value] : j.items()) {
    if (key != "value" && key != "variable" && key != "min" && key != "max") return false;
  }
  return true;
}

struct Walk {
  std::vector<std::string> params;    ///< pointers of Param objects, in document order
  std::vector<std::string> excluded;  ///< other objects with a "value" key
};

void walk(const ojson& j, const std::string& pointer, Walk& out) {
  if (j.is_object()) {
    if (param_object(j) && !pointer.starts_with("/parameters/")) {
      out.params.push_back(pointer);
      return;
    }
    if (j.contains("value")) out.excluded.push_back(pointer);
    for (const auto& [key, value] : j.items()) walk(value, pointer + "/" + key, out);
  } else if (j.is_array()) {
    for (std::size_t i = 0; i < j.size(); ++i) walk(j[i], pointer + "/" + std::to_string(i), out);
  }
}

}  // namespace

TEST_CASE("for_each_param visits exactly the Params of the edit form, in order", "[io][params]") {
  // Objects with "value" that are not Params, by design: the system aperture {"type", "value"}
  // and the rows of the parameter table. Anything else with a "value" key is unexpected here.
  const std::regex excluded_ok(R"(^/aperture$|^/parameters/[0-9]+$)");
  const std::vector<fs::path> files = reference_files();
  REQUIRE(files.size() >= 20U);
  std::size_t total = 0;
  bool tour_seen = false;
  for (const fs::path& file : files) {
    INFO(file.string());
    const rtt::model::System system = rtt::io::load_system(file);
    Walk w;
    walk(ojson::parse(rtt::io::to_edit_json(system)), "", w);
    for (const std::string& p : w.excluded) {
      INFO("object with value at " << p);
      CHECK(std::regex_match(p, excluded_ok));
    }
    std::vector<std::string> visited;
    rtt::model::for_each_param(system, [&](std::string_view pointer, const rtt::model::Param&) {
      visited.emplace_back(pointer);
    });
    CHECK(visited == w.params);
    CHECK(std::set<std::string>(visited.begin(), visited.end()).size() == visited.size());
    total += visited.size();
    tour_seen = tour_seen || file.filename() == "feature_tour.rtt.json";
  }
  CHECK(tour_seen);
  CHECK(total > 0);
}

TEST_CASE("for_each_param reaches bound Params and every field kind", "[io][params]") {
  // The feature tour has asphere coefficients, Zernike terms, both phase layers and bounds;
  // a binding is added here (the tour itself binds only from #165 on).
  rtt::model::System s =
      rtt::io::load_system(fs::path(RTT_REFERENCE_DIR) / "m0" / "feature_tour.rtt.json");
  REQUIRE_FALSE(s.parameters.empty());
  s.aperture.value = rtt::model::Param::bound(s.parameters[0].name);
  std::vector<std::string> kinds;
  bool bound_seen = false;
  rtt::model::for_each_param(s, [&](std::string_view pointer, const rtt::model::Param& p) {
    bound_seen = bound_seen || (pointer == "/aperture/value" && p.is_bound());
    for (const char* kind : {"/coefficients/", "/terms/", "/lines_per_mm", "/normalization_radius",
                             "/conic", "/radius"}) {
      if (pointer.find(kind) != std::string_view::npos) kinds.emplace_back(kind);
    }
  });
  CHECK(bound_seen);
  for (const char* kind : {"/coefficients/", "/terms/", "/lines_per_mm", "/normalization_radius",
                           "/conic", "/radius"}) {
    INFO(kind);
    CHECK(std::find(kinds.begin(), kinds.end(), std::string(kind)) != kinds.end());
  }
}
