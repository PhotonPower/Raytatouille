// aperture.stop (#53, ADR 0039, points 1 and 9): the stop on a surface, in the reference form of
// ADR 0038 ({"surface": "<id>"}), within schema 0.8; in a file before 0.8 it is a read error.

#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <string>

#include "rtt/io/edit.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/model/system.hpp"
#include "test_support.hpp"

using nlohmann::json;
using rtt::model::SurfaceId;
using rtt::model::System;

namespace {

/// The singlet without its stop element, the stop on L1.S1 (which has an aperture).
System surface_stop() {
  System s = rtt::model::test::make_singlet();
  s.root.children.erase(s.root.children.begin());  // the stop element
  s.aperture.stop = SurfaceId("L1.S1");
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

TEST_CASE("aperture.stop: the reference form round trip", "[io][0.8][aperture-stop]") {
  const System s = surface_stop();
  const json j = json::parse(rtt::io::to_json(s));
  CHECK(j["aperture"]["stop"] == json{{"surface", "L1.S1"}});
  // Written after value (ADR 0039, point 9): type, value, stop.
  json keys = json::array();
  for (const auto& [key, value] : j["aperture"].items()) keys.push_back(key);
  CHECK(keys == json{"type", "value", "stop"});
  const System back = rtt::io::parse_system(rtt::io::to_json(s));
  CHECK(back == s);
  REQUIRE(back.aperture.stop.has_value());
  CHECK(back.aperture.stop->str() == "L1.S1");
  // Canonical: written again, the same bytes.
  CHECK(rtt::io::to_json(back) == rtt::io::to_json(s));
}

TEST_CASE("aperture.stop: not written without a stop surface", "[io][0.8][aperture-stop]") {
  const System s = rtt::model::test::make_singlet();
  const json j = json::parse(rtt::io::to_json(s));
  CHECK_FALSE(j["aperture"].contains("stop"));
  // Also not in the edit form (like object.reference): a patch adds /aperture/stop.
  const json e = json::parse(rtt::io::to_edit_json(s));
  CHECK_FALSE(e["aperture"].contains("stop"));
  const json es = json::parse(rtt::io::to_edit_json(surface_stop()));
  CHECK(es["aperture"]["stop"] == json{{"surface", "L1.S1"}});
}

TEST_CASE("aperture.stop: other forms are read errors at its pointer", "[io][0.8][aperture-stop]") {
  // Only the reference form; "origin" is a form of object.reference, not of aperture.stop.
  for (const json& bad :
       {json("L1.S1"), json("origin"), json(1.0), json::object(), json{{"node", "L1"}},
        json{{"surface", "L1.S1"}, {"node", "L1"}}, json{{"surface", 3}}}) {
    INFO(bad.dump());
    json j = json::parse(rtt::io::to_json(surface_stop()));
    j["aperture"]["stop"] = bad;
    CHECK(error_pointer(j.dump()).starts_with("/aperture/stop"));
  }
}

TEST_CASE("aperture.stop: a read error in a file before 0.8", "[io][0.8][aperture-stop]") {
  json j = json::parse(rtt::io::to_json(surface_stop()));
  j["schema_version"] = "0.7.0";
  CHECK(error_pointer(j.dump()) == "/aperture/stop");
  // Content guard: the same 0.7 file without aperture.stop is read.
  j["aperture"].erase("stop");
  CHECK(error_pointer(j.dump()) == "<no error>");
}

TEST_CASE("aperture.stop: an edit patch adds and removes it", "[io][0.8][aperture-stop]") {
  System s = surface_stop();
  s.aperture.stop.reset();
  const rtt::io::PatchResult added = rtt::io::apply_patch(
      s, json::parse(R"([{"op": "add", "path": "/aperture/stop", "value": {"surface": "L1.S1"}}])")
             .dump());
  REQUIRE(added.system.aperture.stop.has_value());
  CHECK(added.system.aperture.stop->str() == "L1.S1");
  const rtt::io::PatchResult removed = rtt::io::apply_patch(
      added.system, json::parse(R"([{"op": "remove", "path": "/aperture/stop"}])").dump());
  CHECK_FALSE(removed.system.aperture.stop.has_value());
}
