// Migration 0.2 -> 0.3 (ADR 0025, point 4): "diffract" becomes "transmit" with the same order.
// Both kept the medium, so the compiled media of the path do not change.

#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <string>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"

using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::material::MaterialLibrary;
using rtt::model::Element;
using rtt::model::ElementKind;
using rtt::model::EventKind;
using rtt::model::Param;
using rtt::model::Pose;
using rtt::model::Surface;
using rtt::model::SurfaceId;
using rtt::model::System;

namespace {

Surface plane(const std::string& id, double z_mm = 0.0) {
  Surface s;
  s.id = SurfaceId(id);
  s.pose = Pose::along_z(z_mm);
  return s;
}

/// Lens L (CONST:1.5) with a grating on its exit surface and plate P (CONST:1.6) with a grating
/// on its entrance surface; the path diffracts at L.S2 from inside the lens and at P.S1 from
/// outside the plate, written with `kind` at both.
System system_with(EventKind kind) {
  System s;
  s.wavelengths = {{0.5876, 1.0, true}};
  s.aperture = {rtt::model::SystemApertureType::EntrancePupilDiameter, Param(10.0)};
  s.fields = {rtt::model::FieldType::AngleDeg, {{0.0, 0.0, 1.0}}};
  s.root.name = "root";
  Surface exit = plane("L.S2", 5.0);
  exit.phases.push_back(rtt::model::LinearGrating{Param(300.0), 0.0});
  s.root.children.push_back(
      {Element{"L", ElementKind::Lens, Pose::along_z(10.0), "CONST:1.5", {plane("L.S1"), exit}}});
  Surface entrance = plane("P.S1");
  entrance.phases.push_back(rtt::model::LinearGrating{Param(600.0), 0.0});
  s.root.children.push_back({Element{
      "P", ElementKind::Plate, Pose::along_z(30.0), "CONST:1.6", {entrance, plane("P.S2", 3.0)}}});
  s.paths = {{"main",
              false,
              {{SurfaceId("L.S1"), EventKind::Refract, 0},
               {SurfaceId("L.S2"), kind, 1},
               {SurfaceId("L.S2"), EventKind::Refract, 0},
               {SurfaceId("P.S1"), kind, -1},
               {SurfaceId("P.S1"), EventKind::Refract, 0},
               {SurfaceId("P.S2"), EventKind::Refract, 0}}}};
  return s;
}

std::string replace_all(std::string text, const std::string& from, const std::string& to) {
  std::size_t count = 0;
  for (std::size_t at = text.find(from); at != std::string::npos;
       at = text.find(from, at + to.size())) {
    text.replace(at, from.size(), to);
    ++count;
  }
  REQUIRE(count > 0);
  return text;
}

struct Media {
  std::string before;
  std::string after;
  std::string beyond;
  bool from_inside = false;
};

std::vector<Media> media(const CompiledSystem& cs) {
  std::vector<Media> out;
  for (const auto& e : cs.path(PathId{0}).events) {
    out.push_back({cs.media()[e.medium_before].reference, cs.media()[e.medium_after].reference,
                   cs.media()[e.medium_beyond].reference, e.from_inside});
  }
  return out;
}

}  // namespace

TEST_CASE("migration: diffract -> transmit keeps the compiled media of the path",
          "[compile][migration]") {
  // The 0.2 file as the 0.2 writer wrote it: "diffract" instead of "transmit".
  const std::string current = rtt::io::to_json(system_with(EventKind::Transmit));
  const std::string old = replace_all(
      replace_all(current, R"("schema_version": "0.4.0")", R"("schema_version": "0.2.0")"),
      R"("kind": "transmit")", R"("kind": "diffract")");
  const System migrated = rtt::io::parse_system(old);
  REQUIRE(migrated == system_with(EventKind::Transmit));

  const MaterialLibrary lib;
  const std::vector<Media> got = media(rtt::compile::compile(migrated, lib));
  // Media of the 0.2 rules (decided in #5, docs/architecture.md before 0.3): Diffract, like
  // Reflect and Transmit, kept the medium; beyond is the medium a Refract would enter.
  const std::vector<Media> expected = {
      {"AIR", "CONST:1.5", "CONST:1.5", false},  // L.S1 refract: enter the lens
      {"CONST:1.5", "CONST:1.5", "AIR", true},   // L.S2 diffract, inside the lens
      {"CONST:1.5", "AIR", "AIR", true},         // L.S2 refract: leave
      {"AIR", "AIR", "CONST:1.6", false},        // P.S1 diffract, outside the plate
      {"AIR", "CONST:1.6", "CONST:1.6", false},  // P.S1 refract: enter the plate
      {"CONST:1.6", "AIR", "AIR", true},         // P.S2 refract: leave
  };
  REQUIRE(got.size() == expected.size());
  for (std::size_t k = 0; k < got.size(); ++k) {
    INFO("event " << k);
    REQUIRE(got[k].before == expected[k].before);
    REQUIRE(got[k].after == expected[k].after);
    REQUIRE(got[k].beyond == expected[k].beyond);
    REQUIRE(got[k].from_inside == expected[k].from_inside);
  }
}
