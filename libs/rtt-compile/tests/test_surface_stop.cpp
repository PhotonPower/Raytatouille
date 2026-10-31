// The stop on a surface (#53, ADR 0039, points 3 and 6): aperture.stop makes that surface the
// stop surface; its event keeps its kind (the lens surface refracts), and a refracting stop
// surface is a reflection candidate of the ghost generator.

#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>

#include "rtt/coating/catalog.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/compile/errors.hpp"
#include "rtt/compile/ghosts.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"

using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::material::MaterialLibrary;
using rtt::model::EventKind;
using rtt::model::SurfaceId;
using rtt::model::System;

namespace {

System load(const std::string& relative) {
  return rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/" + relative);
}

const MaterialLibrary& schott_m2() {
  struct Lib {  // MaterialLibrary is not movable
    MaterialLibrary lib;
    Lib() { lib.add_catalog(std::string(RTT_CATALOG_DIR) + "/m2/schott.agf"); }
  };
  static const Lib l;
  return l.lib;
}

}  // namespace

TEST_CASE("surface stop: L2.S1 of the triplet is the stop surface and refracts",
          "[compile][stop][aperture-stop]") {
  // m2/cooke_triplet_surface_stop: path auto = L1.S1, L1.S2, L2.S1, L2.S2, L3.S1, L3.S2, IMG.
  const CompiledSystem cs =
      rtt::compile::compile(load("m2/cooke_triplet_surface_stop.rtt.json"), schott_m2());
  REQUIRE(cs.stop_surface().has_value());
  CHECK(cs.stop_surface() == cs.find_surface(SurfaceId("L2.S1")));
  CHECK(cs.surfaces()[*cs.stop_surface()].element_kind == rtt::model::ElementKind::Lens);
  const PathId main{0};
  REQUIRE(rtt::compile::stop_event(cs, main) == std::optional<std::size_t>{2});
  CHECK(cs.path(main).events[2].kind == EventKind::Refract);
  CHECK_NOTHROW(rtt::compile::require_stop(cs, main));
  CHECK(cs.diagnostics().empty());
}

TEST_CASE("surface stop: a path that does not visit it", "[compile][stop][aperture-stop]") {
  System s = load("m2/cooke_triplet_surface_stop.rtt.json");
  rtt::model::Path bypass;
  bypass.name = "bypass";
  bypass.automatic = false;
  for (const char* id : {"L3.S1", "L3.S2", "IMG"}) {
    bypass.events.push_back(
        {SurfaceId(id), std::string(id) == "IMG" ? EventKind::Transmit : EventKind::Refract, 0});
  }
  s.paths.push_back(bypass);
  const CompiledSystem cs = rtt::compile::compile(s, schott_m2());
  const PathId id{1};
  CHECK(!rtt::compile::stop_event(cs, id).has_value());
  CHECK_THROWS_AS(rtt::compile::require_stop(cs, id), rtt::compile::NoStopError);
  bool warned = false;
  for (const auto& d : cs.diagnostics()) {
    warned = warned || (d.code == "stop.not_on_path" && d.location == "/paths/1");
  }
  CHECK(warned);
}

TEST_CASE("surface stop: 15 ghosts, the stop surface is a reflection candidate",
          "[compile][stop][aperture-stop][ghosts]") {
  // Six refracting surfaces as with the stop element (STO is no candidate there, L2.S1 refracts
  // in both): C(6, 2) = 15 ghosts (ADR 0039, point 6). The ghost (j, i) reflects at j, then at
  // i < j; with the surfaces 0..5 = L1.S1 .. L3.S2 and the stop s = 2 (L2.S1) the path meets s
  // [s <= j] + [i <= s < j] + [s > i] times: more than once for j >= 3, i <= 2 (9 ghosts) and
  // for j = 2, i in {0, 1} (2 ghosts, which reflect at the stop first).
  const rtt::coating::CoatingLibrary coatings;
  const rtt::compile::GhostSystem g = rtt::compile::compile_with_ghosts(
      load("m2/cooke_triplet_surface_stop.rtt.json"), "main", schott_m2(), coatings);
  const CompiledSystem& cs = g.system;
  REQUIRE(cs.paths().size() == 16);  // main and 15 ghosts
  REQUIRE(cs.stop_surface().has_value());
  const std::uint32_t stop = *cs.stop_surface();
  std::size_t multiple = 0;
  std::size_t first_reflects = 0;
  for (std::size_t p = 0; p < cs.paths().size(); ++p) {
    const PathId id{static_cast<std::uint32_t>(p)};
    const auto& events = cs.path(id).events;
    std::size_t passes = 0;
    for (const auto& e : events) passes += e.surface == stop ? 1 : 0;
    if (passes > 1) ++multiple;
    const std::optional<std::size_t> e = rtt::compile::stop_event(cs, id);
    INFO("path " << cs.path(id).name);
    REQUIRE(e.has_value());
    CHECK(events[*e].surface == stop);
    CHECK(events[*e].kind != EventKind::Transmit);
    if (events[*e].kind == EventKind::Reflect) ++first_reflects;
  }
  CHECK(multiple == 11);
  CHECK(first_reflects == 2);
}
