// The ideal cylinder lens has no paraxial model (#178, ADR 0031, point 6): it is not
// rotationally symmetric, so first_order throws ParaxialError at its surface. (Until PR 3 the
// ideal lens did so as well; that interim state ended with PR 3.)

#include <catch2/catch_test_macros.hpp>
#include <optional>
#include <string>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/paraxial/paraxial.hpp"

using rtt::compile::CompiledSystem;

TEST_CASE("first_order rejects a path through an ideal cylinder lens", "[ideal][interim]") {
  const CompiledSystem cs = rtt::compile::compile(
      rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/r2/ideal_lens.rtt.json"),
      rtt::material::MaterialLibrary{});
  try {
    static_cast<void>(rtt::paraxial::first_order(cs, *cs.find_path("cylinder lens"), 0));
    FAIL("first_order accepted an ideal cylinder lens");
  } catch (const rtt::paraxial::ParaxialError& e) {
    REQUIRE(e.surface().has_value());
    CHECK(e.surface()->str() == "CL");
    CHECK(e.location() == std::optional<std::string>{"/root/children/2/surfaces/0"});
  }
  // The ideal lens of the same file has a paraxial model.
  CHECK_NOTHROW(rtt::paraxial::first_order(cs, *cs.find_path("ideal lens"), 0));
}
