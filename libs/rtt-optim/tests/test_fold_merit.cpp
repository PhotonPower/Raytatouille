// The optimization does not take folded paths yet (#357, ADR 0054, point 4): first_order has
// data for a plane fold mirror since PR A, but the nominal basis of an overlay needs the
// prescription of the path, which follows in PR B. Every operand and generator on such a path is
// merit.operand_unsupported, with and without overlay alike (test_trial_overlay.cpp).

#include <catch2/catch_test_macros.hpp>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "rtt/io/json_io.hpp"
#include "rtt/model/system.hpp"
#include "rtt/optim/merit.hpp"
#include "rtt/optim/optimize.hpp"

TEST_CASE("operands on a folded path are not supported yet", "[optim][fold][357]") {
  const std::ifstream in(std::string(RTT_REFERENCE_DIR) + "/m6/singlet_fold45.rtt.json");
  REQUIRE(in);
  std::ostringstream buffer;
  buffer << in.rdbuf();
  std::string text = buffer.str();
  const std::size_t end = text.rfind('}');
  REQUIRE(end != std::string::npos);
  text.insert(end, R"(, "optimization": {
    "operands": [{"type": "efl", "path": "main", "target": 100.0, "weight": 1.0}],
    "generators": [{"type": "rms_spot", "path": "main"}]})");
  const rtt::model::System s = rtt::io::parse_system(text);
  const rtt::material::MaterialLibrary lib;
  try {
    const rtt::optim::MeritFunction merit(s, lib, nullptr, nullptr);
    FAIL("the merit function was constructed");
  } catch (const rtt::optim::OptimError& e) {
    const std::vector<rtt::model::Diagnostic> d = e.diagnostics();
    REQUIRE(d.size() == 2);
    for (const auto& x : d) {
      CHECK(x.code == "merit.operand_unsupported");
      CHECK(x.message.find("folded path") != std::string::npos);
    }
  }
}
