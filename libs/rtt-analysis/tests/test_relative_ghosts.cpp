// Ghost ranking of relatively placed systems (ADR 0028, acceptance; #163): the M4 ghost
// reference systems rewritten to relative poses rank their ghosts as the absolute files do.

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <map>
#include <string>

#include "relative_placement.hpp"
#include "rtt/analysis/ghosts.hpp"
#include "rtt/coating/catalog.hpp"
#include "rtt/compile/ghosts.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"

using rtt::analysis::GhostRanking;
using rtt::compile::GhostSystem;

namespace {

GhostRanking rank(const rtt::model::System& s) {
  const rtt::material::MaterialLibrary lib;
  const rtt::coating::CoatingLibrary coatings;
  const GhostSystem g = rtt::compile::compile_with_ghosts(s, "main", lib, coatings);
  return rtt::analysis::ghost_ranking(g, 0, 0);
}

/// Ranking entries by ghost path name.
std::map<std::string, const rtt::analysis::GhostEntry*> by_name(const GhostRanking& r,
                                                                const rtt::model::System& s) {
  const rtt::material::MaterialLibrary lib;
  const rtt::coating::CoatingLibrary coatings;
  const GhostSystem g = rtt::compile::compile_with_ghosts(s, "main", lib, coatings);
  std::map<std::string, const rtt::analysis::GhostEntry*> out;
  for (const auto& e : r.entries) out[g.system.path(e.path).name] = &e;
  return out;
}

bool close_rel(double a, double b, double rel) {
  return std::abs(a - b) <= rel * std::max(std::abs(a), std::abs(b));
}

}  // namespace

TEST_CASE("relative ghosts: ranking of a relatively placed system equals the absolute one",
          "[ghosts][ranking][relative]") {
  // The rewrite (relative_placement.hpp: to_relative) keeps every surface transform to
  // delta_t = 1e-10 mm and delta_r = 1e-12 rad (test_relative_placement.cpp). Over at most
  // 150 mm a ray then moves by at most delta = 1e-10 + 150 * 1e-12 mm < 1e-9 mm on the image.
  // Hence: RMS radii to 1e-9 mm (absolute); powers (Fresnel factors, smooth in the angles of
  // incidence) to 1e-9 relative; rho = P (r_b^2 + r0^2) / (r_g^2 + r0^2) to
  // 2 delta (r_b / (r_b^2 + r0^2) + r_g / (r_g^2 + r0^2)) <= 2 delta / r0 = 4e-7 relative with
  // r0 = 0.005 mm, so 1e-6.
  for (const char* file : {"m4/ghost_plates.rtt.json", "m4/ghost_singlet.rtt.json"}) {
    INFO(file);
    const rtt::model::System absolute =
        rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/" + file);
    const rtt::model::test::RelativeSystem relative = rtt::model::test::to_relative(absolute);
    REQUIRE(relative.preceding + relative.sibling > 0);
    const GhostRanking a = rank(absolute);
    const GhostRanking b = rank(relative.system);
    REQUIRE(!a.entries.empty());
    REQUIRE(a.entries.size() == b.entries.size());
    CHECK(close_rel(a.base_power, b.base_power, 1e-9));
    CHECK(std::abs(a.base_rms_radius - b.base_rms_radius) <= 1e-9);
    const auto ea = by_name(a, absolute);
    const auto eb = by_name(b, relative.system);
    for (const auto& [name, x] : ea) {
      INFO(name);
      REQUIRE(eb.count(name) == 1);
      const auto* y = eb.at(name);
      CHECK(close_rel(x->relative_power, y->relative_power, 1e-9));
      CHECK(std::abs(x->rms_radius - y->rms_radius) <= 1e-9);
      CHECK(close_rel(x->relative_irradiance, y->relative_irradiance, 1e-6));
      CHECK(x->rays_arrived == y->rays_arrived);
    }
    // The order is the same wherever the rank values differ by more than the tolerance.
    for (std::size_t k = 0; k + 1 < a.entries.size(); ++k) {
      const double r0 = a.entries[k].relative_irradiance;
      const double r1 = a.entries[k + 1].relative_irradiance;
      if (close_rel(r0, r1, 2e-6)) continue;
      CHECK(b.entries[k].path.index == a.entries[k].path.index);
    }
  }
}
