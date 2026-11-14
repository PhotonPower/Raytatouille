// Schema 0.10 (#380, ADR 0049, point 1): the flag "active" of assemblies and elements; the
// migration 0.9 -> 0.10 changes only the version.

#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <string>
#include <variant>
#include <vector>

#include "rtt/io/edit.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/model/parameters.hpp"
#include "rtt/model/system.hpp"
#include "rtt/model/validate.hpp"
#include "test_support.hpp"

using nlohmann::json;
using rtt::model::Assembly;
using rtt::model::Node;
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

/// The singlet with its lens L1 (child 1) inside an assembly "group" at child 1.
System grouped_singlet() {
  System s = rtt::model::test::make_singlet();
  Assembly group;
  group.name = "group";
  group.children.push_back(s.root.children[1]);
  s.root.children[1] = Node{group};
  return s;
}

}  // namespace

TEST_CASE("schema version is 0.10.0", "[io][0.10]") {
  REQUIRE(rtt::model::kSchemaVersion == "0.10.0");
}

TEST_CASE("active: the canonical form writes it only if false, after the name", "[io][0.10]") {
  System s = rtt::model::test::make_singlet();
  const json all_active = json::parse(rtt::io::to_json(s));
  CHECK_FALSE(all_active["root"].contains("active"));
  CHECK_FALSE(all_active["root"]["children"][1].contains("active"));

  rtt::model::test::element(s, 1).active = false;
  const std::string text = rtt::io::to_json(s);
  const json j = json::parse(text);
  CHECK(j["root"]["children"][1]["active"] == false);
  CHECK_FALSE(j["root"]["children"][2].contains("active"));
  // Key order: type, name, active, then the rest.
  CHECK(text.find(R"("name": "L1",)") < text.find(R"("active": false)"));
  // Round trip, bitwise and as a model.
  const System back = rtt::io::parse_system(text);
  CHECK(back == s);
  CHECK(rtt::io::to_json(back) == text);
  CHECK_FALSE(std::get<rtt::model::Element>(back.root.children[1].value).active);
}

TEST_CASE("active: an inactive assembly round-trips and compares", "[io][0.10]") {
  System s = grouped_singlet();
  std::get<Assembly>(s.root.children[1].value).active = false;
  const std::string text = rtt::io::to_json(s);
  CHECK(json::parse(text)["root"]["children"][1]["active"] == false);
  const System back = rtt::io::parse_system(text);
  CHECK(back == s);
  // Assembly::operator== is written by hand (element.cpp): it compares active too.
  System on = s;
  std::get<Assembly>(on.root.children[1].value).active = true;
  CHECK_FALSE(on == s);
}

TEST_CASE("active: the edit form writes it always; a patch toggles it with undo", "[io][0.10]") {
  const System s = grouped_singlet();
  const json edit = json::parse(rtt::io::to_edit_json(s));
  CHECK(edit["root"]["active"] == true);
  CHECK(edit["root"]["children"][1]["active"] == true);                 // the assembly
  CHECK(edit["root"]["children"][1]["children"][0]["active"] == true);  // L1 inside it

  const rtt::io::PatchResult off = rtt::io::apply_patch(
      s, R"([{"op": "replace", "path": "/root/children/1/children/0/active", "value": false}])");
  const auto& group = std::get<Assembly>(off.system.root.children[1].value);
  CHECK_FALSE(std::get<rtt::model::Element>(group.children[0].value).active);
  // No new errors: switching a part off is always accepted (ADR 0049, point 8).
  CHECK_FALSE(rtt::model::has_errors(rtt::model::validate(off.system)));
  const rtt::io::PatchResult back =
      rtt::io::apply_patch(off.system, off.inverse, rtt::io::EditCheck::StructureOnly);
  CHECK(back.system == s);
}

TEST_CASE("active: a wrong type is a read error at the key", "[io][0.10]") {
  json j = json::parse(rtt::io::to_json(rtt::model::test::make_singlet()));
  j["root"]["children"][1]["active"] = "no";
  CHECK(error_pointer(j.dump()) == "/root/children/1/active");
}

TEST_CASE("migration 0.9 -> 0.10: content unchanged, the version becomes 0.10.0", "[io][0.10]") {
  const System s = rtt::model::test::make_singlet();
  const System read = rtt::io::parse_system(with_version(s, "0.9.0"));
  CHECK(read == s);
  CHECK(read.schema_version == "0.10.0");
  CHECK(json::parse(rtt::io::to_json(read))["schema_version"] == "0.10.0");
}

TEST_CASE("migration: the key active in a file before 0.10 is an error at the key", "[io][0.10]") {
  System s = grouped_singlet();
  std::get<Assembly>(s.root.children[1].value).active = false;
  CHECK(error_pointer(with_version(s, "0.9.0")) == "/root/children/1/active");
  System e = rtt::model::test::make_singlet();
  rtt::model::test::element(e, 2).active = false;
  CHECK(error_pointer(with_version(e, "0.8.0")) == "/root/children/2/active");
  CHECK(error_pointer(with_version(e, "0.10.0")) == "<no error>");
  // At the root too, and in a file before 0.4 that goes through migrate_0_4 first.
  json root = json::parse(with_version(rtt::model::test::make_singlet(), "0.9.0"));
  root["root"]["active"] = true;
  CHECK(error_pointer(root.dump()) == "/root/active");
  CHECK(error_pointer(with_version(e, "0.3.0")) == "/root/children/2/active");
}

TEST_CASE("inactive_nodes lists the topmost inactive nodes; in_inactive_node tests pointers",
          "[io][0.10]") {
  System s = grouped_singlet();
  CHECK(rtt::model::inactive_nodes(s).empty());
  auto& group = std::get<Assembly>(s.root.children[1].value);
  std::get<rtt::model::Element>(group.children[0].value).active = false;  // L1 in the group
  CHECK(rtt::model::inactive_nodes(s) == std::vector<std::string>{"/root/children/1/children/0"});
  group.active = false;  // now the whole group; L1 below it is not listed again
  const std::vector<std::string> inactive = rtt::model::inactive_nodes(s);
  CHECK(inactive == std::vector<std::string>{"/root/children/1"});
  CHECK(rtt::model::in_inactive_node(inactive, "/root/children/1"));
  CHECK(rtt::model::in_inactive_node(inactive, "/root/children/1/children/0/pose/position/2"));
  CHECK_FALSE(rtt::model::in_inactive_node(inactive, "/root/children/10/pose/position/2"));
  CHECK_FALSE(rtt::model::in_inactive_node(inactive, "/root/children/2"));
  CHECK_FALSE(rtt::model::in_inactive_node(inactive, "/aperture/value"));
  // for_each_param still visits every Param, inactive or not (ADR 0049, point 1).
  std::size_t all = 0;
  std::size_t in_inactive = 0;
  rtt::model::for_each_param(s, [&](std::string_view pointer, const rtt::model::Param&) {
    ++all;
    if (rtt::model::in_inactive_node(inactive, pointer)) ++in_inactive;
  });
  CHECK(in_inactive > 0);
  CHECK(in_inactive < all);
}
