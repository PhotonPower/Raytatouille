#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "rtt/coating/catalog.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/compile/ghosts.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"

// Ghost generator (#123, ADR 0027): two-reflection ghosts of a path as explicit model paths.

using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::material::MaterialLibrary;
using rtt::model::Element;
using rtt::model::ElementKind;
using rtt::model::Event;
using rtt::model::EventKind;
using rtt::model::Param;
using rtt::model::Pose;
using rtt::model::Surface;
using rtt::model::SurfaceId;
using rtt::model::System;

namespace {

System load(const std::string& relative) {
  return rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/" + relative);
}

PathId path(const CompiledSystem& cs, const std::string& name) {
  const std::optional<PathId> id = cs.find_path(name);
  REQUIRE(id.has_value());
  return *id;
}

Event event(const std::string& surface, EventKind kind) {
  return Event{SurfaceId(surface), kind, 0};
}

Surface plane(const std::string& id, double z) {
  Surface s;
  s.id = SurfaceId(id);
  s.pose = Pose::along_z(z);
  return s;
}

/// Number of Refract events of a compiled path.
std::size_t refract_events(const CompiledSystem& cs, PathId id) {
  std::size_t n = 0;
  for (const auto& e : cs.path(id).events) n += e.kind == EventKind::Refract ? 1 : 0;
  return n;
}

}  // namespace

TEST_CASE("ghosts: a plane plate has exactly one ghost (#123)", "[ghosts]") {
  // tests/reference/m3/fresnel_bk7.rtt.json, automatic path "main": STO (Transmit), P.S1 and
  // P.S2 (Refract), IMG (Transmit). The one ghost reflects at P.S2 back and at P.S1 forward.
  const MaterialLibrary lib;
  const CompiledSystem cs = rtt::compile::compile(load("m3/fresnel_bk7.rtt.json"), lib);
  const auto ghosts = rtt::compile::ghost_paths(cs, path(cs, "main"));
  REQUIRE(ghosts.size() == 1);
  REQUIRE(ghosts[0].name == "main ghost P.S2/P.S1");
  REQUIRE_FALSE(ghosts[0].automatic);
  const std::vector<Event> expected{
      event("STO", EventKind::Transmit), event("P.S1", EventKind::Refract),
      event("P.S2", EventKind::Reflect), event("P.S1", EventKind::Reflect),
      event("P.S2", EventKind::Refract), event("IMG", EventKind::Transmit)};
  REQUIRE(ghosts[0].events == expected);
}

TEST_CASE("ghosts: N refracting events give N (N - 1) / 2 ghosts in a fixed order (#123)",
          "[ghosts]") {
  MaterialLibrary lib;
  lib.add_catalog(std::string(RTT_CATALOG_DIR) + "/schott.agf");
  MaterialLibrary lib_m2;
  lib_m2.add_catalog(std::string(RTT_CATALOG_DIR) + "/m2/schott.agf");

  SECTION("cemented achromat: 3 refracting surfaces, 3 ghosts") {
    // m2/achromat: STO, L1.S1, L1.S2 (cemented), L1.S3, IMG. Order j ascending, then i.
    const CompiledSystem cs = rtt::compile::compile(load("m2/achromat.rtt.json"), lib);
    const PathId main = path(cs, "main");
    REQUIRE(refract_events(cs, main) == 3);
    const auto ghosts = rtt::compile::ghost_paths(cs, main);
    REQUIRE(ghosts.size() == 3);
    REQUIRE(ghosts[0].name == "main ghost L1.S2/L1.S1");
    REQUIRE(ghosts[1].name == "main ghost L1.S3/L1.S1");
    REQUIRE(ghosts[2].name == "main ghost L1.S3/L1.S2");
    // j = L1.S3, i = L1.S1: the way back crosses the cemented surface L1.S2 (Refract).
    const std::vector<Event> expected{
        event("STO", EventKind::Transmit),  event("L1.S1", EventKind::Refract),
        event("L1.S2", EventKind::Refract), event("L1.S3", EventKind::Reflect),
        event("L1.S2", EventKind::Refract), event("L1.S1", EventKind::Reflect),
        event("L1.S2", EventKind::Refract), event("L1.S3", EventKind::Refract),
        event("IMG", EventKind::Transmit)};
    REQUIRE(ghosts[1].events == expected);
  }
  SECTION("Cooke triplet: 6 refracting surfaces, 15 ghosts") {
    const CompiledSystem cs = rtt::compile::compile(load("m2/cooke_triplet.rtt.json"), lib_m2);
    const PathId main = path(cs, "main");
    REQUIRE(refract_events(cs, main) == 6);
    REQUIRE(rtt::compile::ghost_paths(cs, main).size() == 15);
  }
}

TEST_CASE("ghosts: a mirror on the base path stays a reflection, repeated surfaces (#123)",
          "[ghosts]") {
  // Plate P (CONST:1.5) passed twice: P.S1, P.S2, mirror M, P.S2, P.S1 (explicit path).
  // Refract events 0, 1, 3, 4 give 6 ghosts. The ghost j = 3 (P.S2, second pass), i = 1
  // (P.S2, first pass) runs back over the mirror: base[0..2], Reflect at P.S2, M (Reflect),
  // Reflect at P.S2, base[2..4]. P.S2 occurs twice, so the event indices are in the name.
  System s;
  s.name = "double pass";
  s.wavelengths = {{0.5876, 1.0, true}};
  s.aperture = {rtt::model::SystemApertureType::EntrancePupilDiameter, Param(5.0)};
  s.fields = {rtt::model::FieldType::AngleDeg, {{0.0, 0.0, 1.0}}};
  s.root.name = "root";
  s.root.children = {{Element{"P",
                              ElementKind::Plate,
                              Pose::along_z(0.0),
                              "CONST:1.5",
                              {plane("P.S1", 0.0), plane("P.S2", 5.0)}}},
                     {Element{"M", ElementKind::Mirror, Pose::along_z(20.0), std::nullopt, {[] {
                                Surface m = plane("M.S", 0.0);
                                m.interaction = rtt::model::IdealMirror{};
                                return m;
                              }()}}}};
  s.paths = {{"main",
              false,
              {event("P.S1", EventKind::Refract), event("P.S2", EventKind::Refract),
               event("M.S", EventKind::Reflect), event("P.S2", EventKind::Refract),
               event("P.S1", EventKind::Refract)}}};
  const MaterialLibrary lib;
  const CompiledSystem cs = rtt::compile::compile(s, lib);
  const auto ghosts = rtt::compile::ghost_paths(cs, PathId{0});
  REQUIRE(ghosts.size() == 6);
  // Order (j, i): (1, 0), (3, 0), (3, 1), (4, 0), (4, 1), (4, 3).
  REQUIRE(ghosts[2].name == "main ghost P.S2#3/P.S2#1");
  const std::vector<Event> expected{
      event("P.S1", EventKind::Refract), event("P.S2", EventKind::Refract),
      event("M.S", EventKind::Reflect),  event("P.S2", EventKind::Reflect),
      event("M.S", EventKind::Reflect),  event("P.S2", EventKind::Reflect),
      event("M.S", EventKind::Reflect),  event("P.S2", EventKind::Refract),
      event("P.S1", EventKind::Refract)};
  REQUIRE(ghosts[2].events == expected);
  REQUIRE(ghosts[0].name == "main ghost P.S2#1/P.S1#0");
}

TEST_CASE("ghosts: limit, name conflicts and invalid base paths (#123)", "[ghosts]") {
  const MaterialLibrary lib;
  System s = load("m3/fresnel_bk7.rtt.json");
  const CompiledSystem cs = rtt::compile::compile(s, lib);
  rtt::compile::GhostOptions none;
  none.max_paths = 0;
  REQUIRE_THROWS_AS(rtt::compile::ghost_paths(cs, path(cs, "main"), none), std::invalid_argument);
  rtt::compile::GhostOptions one;
  one.max_paths = 1;
  REQUIRE(rtt::compile::ghost_paths(cs, path(cs, "main"), one).size() == 1);
  REQUIRE_THROWS_AS(rtt::compile::ghost_paths(cs, PathId{7}), std::invalid_argument);
  // M4: no ghosts for paths with a diffraction order or a crystal mode (ADR 0027).
  System modes = s;
  modes.paths.push_back({"crystal",
                         false,
                         {event("STO", EventKind::Transmit), event("P.S1", EventKind::Ordinary),
                          event("P.S2", EventKind::Refract), event("IMG", EventKind::Transmit)}});
  const CompiledSystem with_mode = rtt::compile::compile(modes, lib);
  REQUIRE_THROWS_AS(rtt::compile::ghost_paths(with_mode, path(with_mode, "crystal")),
                    std::invalid_argument);
  // The generator checks order != 0, not the event kind. Before schema 0.3.0 an order needs a
  // `diffract` event; with ADR 0025 (format PR) this case becomes an order at an event of a
  // surface with a phase layer, and the test is adapted there.
  System grating = s;
  grating.paths.push_back(
      {"grating",
       false,
       {event("STO", EventKind::Transmit), event("P.S1", EventKind::Refract),
        event("P.S2", EventKind::Refract), Event{SurfaceId("IMG"), EventKind::Diffract, 1}}});
  const CompiledSystem with_order = rtt::compile::compile(grating, lib);
  REQUIRE_THROWS_AS(rtt::compile::ghost_paths(with_order, path(with_order, "grating")),
                    std::invalid_argument);
  // A path that already has the ghost's name.
  rtt::model::Path taken = s.paths[1];
  taken.name = "main ghost P.S2/P.S1";
  s.paths.push_back(taken);
  const CompiledSystem conflict = rtt::compile::compile(s, lib);
  REQUIRE_THROWS_AS(rtt::compile::ghost_paths(conflict, path(conflict, "main")),
                    std::invalid_argument);
}

TEST_CASE("compile_with_ghosts: compile determines the media of the ghost (#123)", "[ghosts]") {
  // The plate ghost: glass after the refraction at P.S1, glass on both sides of the way back
  // (reflection at P.S2 from inside, reflection at P.S1 from inside), glass before P.S2 again,
  // vacuum after it. The model is not changed.
  const MaterialLibrary lib;
  const rtt::coating::CoatingLibrary coatings;
  const System s = load("m3/fresnel_bk7.rtt.json");
  const auto g = rtt::compile::compile_with_ghosts(s, "main", lib, coatings);
  REQUIRE(g.ghosts.size() == 1);
  REQUIRE(g.system.paths().size() == s.paths.size() + 1);
  const auto& ghost = g.ghosts[0];
  REQUIRE(ghost.base == *g.system.find_path("main"));
  REQUIRE(ghost.path == *g.system.find_path("main ghost P.S2/P.S1"));
  REQUIRE(ghost.surface_j == *g.system.find_surface(SurfaceId("P.S2")));
  REQUIRE(ghost.surface_i == *g.system.find_surface(SurfaceId("P.S1")));
  REQUIRE(ghost.event_j == 2);
  REQUIRE(ghost.event_i == 1);
  const auto& events = g.system.path(ghost.path).events;
  REQUIRE(events.size() == 6);
  const auto ref = [&](std::uint32_t medium) { return g.system.media()[medium].reference; };
  REQUIRE(ref(events[1].medium_after) == "CONST:1.5168");   // into the plate at P.S1
  REQUIRE(ref(events[2].medium_before) == "CONST:1.5168");  // reflection at P.S2 from inside
  REQUIRE(events[2].from_inside);
  REQUIRE(ref(events[2].medium_after) == "CONST:1.5168");
  REQUIRE(ref(events[3].medium_before) == "CONST:1.5168");  // reflection at P.S1 from inside
  REQUIRE(events[3].from_inside);
  REQUIRE(ref(events[4].medium_before) == "CONST:1.5168");  // out of the plate at P.S2
  REQUIRE(ref(events[4].medium_after) == "VACUUM");
  REQUIRE(s.paths.size() == 2);  // the model itself is unchanged

  // Cemented achromat, ghost L1.S3/L1.S1: on the way back the cemented surface L1.S2 is
  // crossed from F2 into N-BK7 (segment rule backwards), then L1.S1 reflects inside N-BK7.
  MaterialLibrary schott;
  schott.add_catalog(std::string(RTT_CATALOG_DIR) + "/schott.agf");
  const auto a =
      rtt::compile::compile_with_ghosts(load("m2/achromat.rtt.json"), "main", schott, coatings);
  REQUIRE(a.ghosts.size() == 3);
  const auto& back = a.system.path(a.ghosts[1].path).events;
  REQUIRE(back.size() == 9);
  const auto mat = [&](std::uint32_t medium) { return a.system.media()[medium].reference; };
  REQUIRE(mat(back[3].medium_before) == "SCHOTT:F2");     // reflection at L1.S3 inside F2
  REQUIRE(mat(back[4].medium_before) == "SCHOTT:F2");     // L1.S2 backwards ...
  REQUIRE(mat(back[4].medium_after) == "SCHOTT:N-BK7");   // ... into N-BK7
  REQUIRE(mat(back[5].medium_before) == "SCHOTT:N-BK7");  // reflection at L1.S1 inside N-BK7
  REQUIRE(back[5].from_inside);
  REQUIRE_THROWS_AS(rtt::compile::compile_with_ghosts(s, "no such path", lib, coatings),
                    std::invalid_argument);
}
