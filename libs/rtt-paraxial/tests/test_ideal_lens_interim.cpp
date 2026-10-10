// Interim state of #178 until rtt-paraxial knows the ideal lens (ADR 0031, point 8, PR 1 -> PR 3):
// first_order throws ParaxialError at the lens surface, so no EFL without the lens arises. PR 3
// keeps this only for the cylinder lens.

#include <catch2/catch_test_macros.hpp>
#include <optional>
#include <string>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/paraxial/paraxial.hpp"

using rtt::compile::CompiledSystem;

TEST_CASE("interim: first_order rejects a path through an ideal lens", "[ideal][interim]") {
  const CompiledSystem cs = rtt::compile::compile(
      rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/r2/ideal_lens.rtt.json"),
      rtt::material::MaterialLibrary{});
  for (const char* const name : {"ideal lens", "cylinder lens"}) {
    INFO(name);
    try {
      static_cast<void>(rtt::paraxial::first_order(cs, *cs.find_path(name), 0));
      FAIL("first_order accepted an ideal lens");
    } catch (const rtt::paraxial::ParaxialError& e) {
      const std::string lens = std::string(name) == "ideal lens" ? "IL" : "CL";
      REQUIRE(e.surface().has_value());
      CHECK(e.surface()->str() == lens);
    }
  }
}
