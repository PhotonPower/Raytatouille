// The stop surface and the stop event of a path (ADR 0039, point 3): CompiledSystem::stop_surface()
// is the surface of the Stop element, compile::stop_event() the first event of a path at it,
// whatever its kind. Every user of the stop (aiming, paraxial pupils, parabasal operands, the
// loss counter, the obscuration, the system report, require_stop and stop.not_on_path) reads
// these two; the cases below compare them with the lookup the users had before (the first event
// at a surface of a Stop element).

#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

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

/// The lookup of the users before ADR 0039: the first event at a surface of a Stop element.
std::optional<std::size_t> old_stop_event(const CompiledSystem& cs, PathId path) {
  const auto& events = cs.path(path).events;
  for (std::size_t e = 0; e < events.size(); ++e) {
    if (cs.surfaces()[events[e].surface].element_kind == rtt::model::ElementKind::Stop) return e;
  }
  return std::nullopt;
}

}  // namespace

TEST_CASE("stop: the Stop element of the Cooke triplet and its event", "[compile][stop]") {
  // m2/cooke_triplet: path main = L1.S1, L1.S2, L2.S1, STO (transmit), L2.S2, L3.S1, L3.S2, IMG.
  const CompiledSystem cs = rtt::compile::compile(load("m2/cooke_triplet.rtt.json"), schott_m2());
  REQUIRE(cs.stop_surface().has_value());
  CHECK(cs.stop_surface() == cs.find_surface(rtt::model::SurfaceId("STO")));
  CHECK(cs.surfaces()[*cs.stop_surface()].element_kind == rtt::model::ElementKind::Stop);
  const PathId main{0};
  CHECK(rtt::compile::stop_event(cs, main) == std::optional<std::size_t>{3});
  CHECK(rtt::compile::stop_event(cs, main) == old_stop_event(cs, main));
  CHECK_NOTHROW(rtt::compile::require_stop(cs, main));
  CHECK_THROWS_AS((void)rtt::compile::stop_event(cs, PathId{7}), std::out_of_range);
}

TEST_CASE("stop: a system without a Stop element has neither surface nor event",
          "[compile][stop]") {
  const CompiledSystem cs =
      rtt::compile::compile(load("m1/paraboloid_mirror.rtt.json"), MaterialLibrary{});
  CHECK(!cs.stop_surface().has_value());
  for (std::size_t p = 0; p < cs.paths().size(); ++p) {
    const PathId id{static_cast<std::uint32_t>(p)};
    CHECK(!rtt::compile::stop_event(cs, id).has_value());
    CHECK(!old_stop_event(cs, id).has_value());
    CHECK_THROWS_AS(rtt::compile::require_stop(cs, id), rtt::compile::NoStopError);
  }
}

TEST_CASE("stop: a path that does not visit the stop", "[compile][stop]") {
  // A second explicit path of the triplet without STO: the system has a stop surface, the path
  // no stop event; compile warns stop.not_on_path for it (ADR 0023), require_stop throws.
  System s = load("m2/cooke_triplet.rtt.json");
  rtt::model::Path bypass = s.paths.at(0);
  bypass.name = "bypass";
  bypass.automatic = false;
  std::erase_if(bypass.events, [](const rtt::model::Event& e) { return e.surface.str() == "STO"; });
  s.paths.push_back(bypass);
  const CompiledSystem cs = rtt::compile::compile(s, schott_m2());
  REQUIRE(cs.stop_surface().has_value());
  const PathId id{1};
  CHECK(!rtt::compile::stop_event(cs, id).has_value());
  CHECK(!old_stop_event(cs, id).has_value());
  CHECK_THROWS_AS(rtt::compile::require_stop(cs, id), rtt::compile::NoStopError);
  bool warned = false;
  for (const auto& d : cs.diagnostics()) {
    warned = warned || (d.code == "stop.not_on_path" && d.location == "/paths/1");
  }
  CHECK(warned);
}

TEST_CASE("stop: the first of several passes in the ghost paths", "[compile][stop]") {
  // Every ghost of the triplet whose reflections enclose STO passes it three times (forward,
  // back, forward, ADR 0027); the stop event is the first pass, as before. Content guard: at
  // least one ghost passes STO more than once.
  const rtt::coating::CoatingLibrary coatings;
  const rtt::compile::GhostSystem g = rtt::compile::compile_with_ghosts(
      load("m2/cooke_triplet.rtt.json"), "main", schott_m2(), coatings);
  const CompiledSystem& cs = g.system;
  REQUIRE(cs.paths().size() == 16);  // main and 15 ghosts
  const std::uint32_t stop = cs.stop_surface().value();
  std::size_t multiple = 0;
  for (std::size_t p = 0; p < cs.paths().size(); ++p) {
    const PathId id{static_cast<std::uint32_t>(p)};
    const auto& events = cs.path(id).events;
    std::size_t passes = 0;
    for (const auto& e : events) passes += e.surface == stop ? 1 : 0;
    if (passes > 1) ++multiple;
    const std::optional<std::size_t> e = rtt::compile::stop_event(cs, id);
    INFO("path " << cs.path(id).name);
    CHECK(e == old_stop_event(cs, id));
    REQUIRE(e.has_value());
    CHECK(events[*e].surface == stop);
    CHECK(events[*e].kind == EventKind::Transmit);
  }
  CHECK(multiple > 0);
}
