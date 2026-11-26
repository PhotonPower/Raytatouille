// The optimization on folded paths (#357, ADR 0054, point 5, PR B1): first_order has data for a
// plane fold mirror, so the merit function runs on such a path as on its unfolded counterpart,
// without and with the overlay of a tolerance trial (the nominal basis skips the fold as a
// vertex for the BFL of a trial).

#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>

#include "rtt/io/json_io.hpp"
#include "rtt/model/system.hpp"
#include "rtt/optim/merit.hpp"
#include "rtt/optim/overlay.hpp"

namespace {

/// The reference file `name` with an EFL and a BFL operand and an RMS spot generator.
rtt::model::System with_merit(const std::string& name) {
  const std::ifstream in(std::string(RTT_REFERENCE_DIR) + "/" + name);
  REQUIRE(in);
  std::ostringstream buffer;
  buffer << in.rdbuf();
  std::string text = buffer.str();
  const std::size_t end = text.rfind('}');
  REQUIRE(end != std::string::npos);
  text.insert(end, R"(, "optimization": {
    "operands": [{"type": "efl", "path": "main", "target": 100.0, "weight": 1.0},
                 {"type": "bfl", "path": "main", "target": 90.0, "weight": 1.0}],
    "generators": [{"type": "rms_spot", "path": "main"}]})");
  return rtt::io::parse_system(text);
}

rtt::optim::MeritEvaluation evaluate(const rtt::model::System& s, bool overlay) {
  const rtt::material::MaterialLibrary lib;
  std::optional<rtt::optim::TrialOverlay> o;
  if (overlay) o.emplace(rtt::optim::TrialOverlay{{}, rtt::optim::nominal_basis(s, lib, nullptr)});
  const rtt::optim::MeritFunction merit(s, lib, nullptr, o ? &*o : nullptr);
  return merit.evaluate(merit.start());
}

}  // namespace

TEST_CASE("operands on a folded path equal those of the unfolded one", "[optim][fold][357]") {
  const rtt::model::System folded = with_merit("m6/singlet_fold45.rtt.json");
  const rtt::model::System flat = with_merit("m1/singlet_const.rtt.json");
  for (const bool overlay : {false, true}) {
    INFO("overlay " << overlay);
    const auto f = evaluate(folded, overlay);
    const auto u = evaluate(flat, overlay);
    REQUIRE(f.values.size() == u.values.size());
    REQUIRE(f.residuals.size() == u.residuals.size());
    for (std::size_t i = 0; i < f.values.size(); ++i) {
      INFO("value " << i);
      CHECK(std::abs(f.values[i] - u.values[i]) <= 1e-9 * std::max(1.0, std::abs(u.values[i])));
    }
    // The spot residuals are the x and y of each ray in the local frame of the detector, which
    // the fold turns (the detector faces +y): the same magnitudes, signs may differ.
    for (std::size_t i = 0; i < f.residuals.size(); ++i) {
      INFO("residual " << i);
      CHECK(std::abs(std::abs(f.residuals[i]) - std::abs(u.residuals[i])) <= 1e-9);
    }
    // content guard: EFL near 100 mm, BFL near 95 mm (the singlet of m1)
    CHECK(std::abs(u.values[0] - 100.0) < 2.0);
    CHECK(std::abs(u.values[1] - 95.0) < 3.0);
  }
}
