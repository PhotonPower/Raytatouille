// Variables of an optimization run (ADR 0030, point 5; #167): order, pointers, bounds, and
// setting values without touching anything else.

#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <filesystem>
#include <nlohmann/json.hpp>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "rtt/io/edit.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/optim/variables.hpp"

namespace fs = std::filesystem;
using ojson = nlohmann::ordered_json;
using rtt::optim::collect_variables;
using rtt::optim::Variable;
using rtt::optim::with_values;

namespace {

rtt::model::System tour() {
  return rtt::io::load_system(fs::path(RTT_REFERENCE_DIR) / "m0" / "feature_tour.rtt.json");
}

/// Pointers of the variable model Params in the edit form, in document order, independently of
/// for_each_param: objects with "variable": true outside /parameters.
void variable_params(const ojson& j, const std::string& pointer, std::vector<std::string>& out) {
  if (j.is_object()) {
    if (j.contains("value") && j.contains("variable") && j["variable"] == true &&
        !pointer.starts_with("/parameters/")) {
      out.push_back(pointer + "/value");
      return;
    }
    for (const auto& [key, value] : j.items()) variable_params(value, pointer + "/" + key, out);
  } else if (j.is_array()) {
    for (std::size_t i = 0; i < j.size(); ++i) {
      variable_params(j[i], pointer + "/" + std::to_string(i), out);
    }
  }
}

}  // namespace

TEST_CASE("variables: table rows first, then model Params in edit-form order", "[optim][vars]") {
  const rtt::model::System s = tour();
  const std::vector<Variable> vars = collect_variables(s);
  // Feature tour: A_S2_Z is an expression (never variable), FOCUS a `values` row over the two
  // configurations, SCALE a `value` row.
  REQUIRE(vars.size() >= 3);
  CHECK(vars[0].pointer == "/parameters/1/values/0");
  CHECK(vars[0].row == "FOCUS");
  CHECK(vars[0].configuration == std::optional<std::size_t>(0));
  CHECK(vars[0].start == 0.0);
  CHECK(vars[0].bounds.min == std::optional<double>(-1.0));
  CHECK(vars[0].bounds.max == std::optional<double>(1.0));
  CHECK(vars[1].pointer == "/parameters/1/values/1");
  CHECK(vars[1].configuration == std::optional<std::size_t>(1));
  CHECK(vars[1].start == 0.5);
  CHECK(vars[2].pointer == "/parameters/2/value");
  CHECK(vars[2].row == "SCALE");
  CHECK_FALSE(vars[2].configuration.has_value());
  CHECK(vars[2].start == 1.0);
  // The model Params: exactly the variable Param objects of the edit form, in order.
  std::vector<std::string> expected;
  variable_params(ojson::parse(rtt::io::to_edit_json(s)), "", expected);
  REQUIRE(expected.size() >= 3);  // object distance, a radius with bounds, a coefficient
  std::vector<std::string> got;
  for (std::size_t v = 3; v < vars.size(); ++v) {
    got.push_back(vars[v].pointer);
    CHECK(vars[v].row.empty());
    CHECK(vars[v].param_index.has_value());
  }
  CHECK(got == expected);
  CHECK(vars[3].pointer == "/object/distance/value");
  CHECK(vars[3].start == 250.0);
}

TEST_CASE("variables: with_values sets exactly the variables", "[optim][vars]") {
  const rtt::model::System s = tour();
  const std::vector<Variable> vars = collect_variables(s);
  // Content guard (red check #193): three table variables and at least three model Params.
  REQUIRE(vars.size() >= 6);
  std::vector<double> start;
  for (const Variable& v : vars) start.push_back(v.start);
  // The start values give the input back, bitwise in the edit form.
  CHECK(rtt::io::to_edit_json(with_values(s, vars, start)) == rtt::io::to_edit_json(s));
  // New values appear at the pointers; nothing else changes.
  std::vector<double> next;
  for (std::size_t v = 0; v < vars.size(); ++v)
    next.push_back(start[v] + 0.25 * static_cast<double>(v + 1));
  ojson changed = ojson::parse(rtt::io::to_edit_json(with_values(s, vars, next)));
  ojson original = ojson::parse(rtt::io::to_edit_json(s));
  CHECK(changed != original);  // the new values arrived
  for (std::size_t v = 0; v < vars.size(); ++v) {
    INFO(vars[v].pointer);
    const ojson::json_pointer at(vars[v].pointer);
    CHECK(changed[at].get<double>() == next[v]);
    changed[at] = original[at];  // undo this one change
  }
  CHECK(changed == original);
}

TEST_CASE("variables: wrong number of values", "[optim][vars]") {
  const rtt::model::System s = tour();
  const std::vector<Variable> vars = collect_variables(s);
  const std::vector<double> too_few(vars.size() - 1, 0.0);
  CHECK_THROWS_AS(with_values(s, vars, too_few), std::invalid_argument);
}
