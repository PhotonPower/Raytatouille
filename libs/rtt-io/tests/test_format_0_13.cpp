// Schema 0.13 (#196, ADR 0056, point 6): an operand with "constraint": "equal" is an exact
// constraint; it has no weight. The migration 0.12 -> 0.13 changes only the version; a constraint
// in an older file is a read error at its pointer.

#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <string>
#include <vector>

#include "rtt/io/edit.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/model/optimization.hpp"
#include "rtt/model/system.hpp"
#include "rtt/model/validate.hpp"
#include "test_support.hpp"

using nlohmann::json;
using rtt::model::ConstraintKind;
using rtt::model::FirstOrderOperand;
using rtt::model::FirstOrderQuantity;
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

FirstOrderOperand efl(double target, bool constraint) {
  FirstOrderOperand o;
  o.path = "main";
  o.quantity = FirstOrderQuantity::Efl;
  o.common.target = target;
  if (constraint) o.common.constraint = ConstraintKind::Equal;
  return o;
}

/// The singlet with efl = 100 as a constraint and a weighted bfl operand.
System constrained() {
  System s = rtt::model::test::make_singlet();
  s.optimization.operands.push_back(efl(100.0, true));
  FirstOrderOperand bfl = efl(95.0, false);
  bfl.quantity = FirstOrderQuantity::Bfl;
  bfl.common.weight = 2.0;
  s.optimization.operands.push_back(bfl);
  return s;
}

std::vector<std::string> codes(const System& s) {
  std::vector<std::string> out;
  for (const auto& d : rtt::model::validate(s)) out.push_back(d.code);
  return out;
}

}  // namespace

TEST_CASE("schema version is 0.13.0", "[io][0.13]") {
  REQUIRE(rtt::model::kSchemaVersion == "0.13.0");
}

TEST_CASE("a constraint round-trips without a weight", "[io][0.13]") {
  const System s = constrained();
  REQUIRE(codes(s).empty());
  const json j = json::parse(rtt::io::to_json(s));
  const json& op = j["optimization"]["operands"][0];
  CHECK(op["constraint"] == "equal");
  CHECK_FALSE(op.contains("weight"));
  CHECK_FALSE(j["optimization"]["operands"][1].contains("constraint"));
  CHECK(j["optimization"]["operands"][1]["weight"] == 2.0);
  CHECK(rtt::io::parse_system(rtt::io::to_json(s)) == s);
  // the edit form writes every weight, but none for a constraint
  const json e = json::parse(rtt::io::to_edit_json(s));
  CHECK(e["optimization"]["operands"][0]["constraint"] == "equal");
  CHECK_FALSE(e["optimization"]["operands"][0].contains("weight"));
  CHECK(e["optimization"]["operands"][1].contains("weight"));
  CHECK(rtt::io::parse_system(rtt::io::to_edit_json(s)) == s);
}

TEST_CASE("reading a constraint: weight and value are checked", "[io][0.13]") {
  const json j = json::parse(rtt::io::to_json(constrained()));
  json with_weight = j;
  with_weight["optimization"]["operands"][0]["weight"] = 1.0;
  CHECK(error_pointer(with_weight.dump()) == "/optimization/operands/0/weight");
  json unknown = j;
  unknown["optimization"]["operands"][0]["constraint"] = "less";
  CHECK(error_pointer(unknown.dump()) == "/optimization/operands/0/constraint");
}

TEST_CASE("validate: a constraint with a weight; a constraint in a criterion", "[io][0.13]") {
  System s = constrained();
  std::get<FirstOrderOperand>(s.optimization.operands[0]).common.weight = 3.0;
  CHECK(codes(s) == std::vector<std::string>{"merit.constraint_weight"});
  System t = rtt::model::test::make_singlet();
  rtt::model::Criterion k;
  k.id = "efl";
  k.operand = efl(0.0, true);
  k.max = 101.0;
  t.tolerances.criteria.push_back(k);
  const auto diagnostics = rtt::model::validate(t);
  bool found = false;
  for (const auto& d : diagnostics) {
    found = found || (d.code == "tolerances.criterion_bounds" &&
                      d.location == "/tolerances/criteria/0/operand/constraint");
  }
  CHECK(found);
}

TEST_CASE("migration 0.12 -> 0.13: content unchanged, a constraint needs 0.13", "[io][0.13]") {
  System plain = constrained();
  plain.optimization.operands.erase(plain.optimization.operands.begin());
  json j = json::parse(rtt::io::to_json(plain));
  j["schema_version"] = "0.12.0";
  const System read = rtt::io::parse_system(j.dump());
  CHECK(read.schema_version == "0.13.0");
  CHECK(read == plain);
  json old = json::parse(rtt::io::to_json(constrained()));
  old["schema_version"] = "0.12.0";
  CHECK(error_pointer(old.dump()) == "/optimization/operands/0/constraint");
  // also in the compensator merit and in the operand of a criterion
  json merit = json::parse(rtt::io::to_json(plain));
  merit["schema_version"] = "0.12.0";
  merit["tolerances"]["compensator_merit"]["operands"] = json::array(
      {json{{"type", "efl"}, {"path", "main"}, {"target", 100.0}, {"constraint", "equal"}}});
  CHECK(error_pointer(merit.dump()) == "/tolerances/compensator_merit/operands/0/constraint");
  json crit = json::parse(rtt::io::to_json(plain));
  crit["schema_version"] = "0.12.0";
  crit["tolerances"]["criteria"] = json::array(
      {json{{"id", "efl"},
            {"operand", json{{"type", "efl"}, {"path", "main"}, {"constraint", "equal"}}},
            {"max", 101.0}}});
  CHECK(error_pointer(crit.dump()) == "/tolerances/criteria/0/operand/constraint");
}
