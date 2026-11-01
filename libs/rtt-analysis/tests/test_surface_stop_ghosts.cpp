// Ghost ranking with the stop on a surface (#53, ADR 0039, point 6): the refracting stop surface
// L2.S1 is a reflection candidate, the triplet has 15 ghosts as with the workaround, and the
// ranking runs, since it aims only on the base path.

#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <set>
#include <string>
#include <utility>

#include "rtt/analysis/ghosts.hpp"
#include "rtt/coating/catalog.hpp"
#include "rtt/compile/ghosts.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"

namespace {

rtt::compile::GhostSystem ghosts(const std::string& relative) {
  static const rtt::coating::CoatingLibrary coatings;
  struct Lib {  // MaterialLibrary is not movable
    rtt::material::MaterialLibrary lib;
    Lib() { lib.add_catalog(std::string(RTT_CATALOG_DIR) + "/m2/schott.agf"); }
  };
  static const Lib l;
  return rtt::compile::compile_with_ghosts(
      rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/" + relative), "main", l.lib,
      coatings);
}

}  // namespace

TEST_CASE("surface stop: the ghost ranking of the triplet runs",
          "[analysis][ghosts][aperture-stop]") {
  const rtt::compile::GhostSystem g = ghosts("m2/cooke_triplet_surface_stop.rtt.json");
  const std::uint16_t wl = g.system.reference_wavelength();
  const rtt::analysis::GhostRanking r = rtt::analysis::ghost_ranking(g, 2, wl);
  REQUIRE(r.entries.size() == 15);
  CHECK(r.base_power > 0.0);
  CHECK(std::isfinite(r.base_rms_radius));
  // Every pair (j, i) of the six refracting surfaces once, the stop surface L2.S1 among them.
  const std::uint32_t stop = g.system.stop_surface().value();
  std::set<std::pair<std::uint32_t, std::uint32_t>> pairs;
  std::size_t at_stop = 0;
  for (const auto& e : r.entries) {
    INFO(g.system.path(e.path).name);
    CHECK(std::isfinite(e.relative_power));
    CHECK(e.relative_power >= 0.0);
    pairs.emplace(e.surface_j, e.surface_i);
    if (e.surface_j == stop || e.surface_i == stop) ++at_stop;
  }
  CHECK(pairs.size() == 15);
  CHECK(at_stop == 5);  // L2.S1 pairs with each of the other five refracting surfaces
  // Content guard: the workaround has as many ghosts (STO is no candidate there).
  CHECK(rtt::analysis::ghost_ranking(ghosts("m2/cooke_triplet.rtt.json"), 2, wl).entries.size() ==
        15);
}
