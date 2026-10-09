// Interim state of #178 until the tracer applies the ideal lens (ADR 0031, point 8, PR 1 -> PR 2):
// a ray reaching an ideal lens or ideal cylinder lens ends EventImpossible at that surface, never
// silently undeflected. PR 2 replaces this test with the reference tests of the lens.

#include <catch2/catch_test_macros.hpp>
#include <string>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/sequential.hpp"

using rtt::compile::CompiledSystem;
using rtt::trace::RayBatch;
using rtt::trace::RayStatus;

TEST_CASE("interim: a ray at an ideal lens ends EventImpossible there", "[ideal][interim]") {
  const CompiledSystem cs = rtt::compile::compile(
      rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/r2/ideal_lens.rtt.json"),
      rtt::material::MaterialLibrary{});
  for (const char* const name : {"ideal lens", "cylinder lens"}) {
    INFO(name);
    // From the object point on the axis (z = -75), in front of the stop plane z = -1.
    RayBatch rays(1);
    rays.pos_z()[0] = -75.0;
    const rtt::compile::PathId path = *cs.find_path(name);
    static_cast<void>(rtt::trace::SequentialTracer{}.trace(cs, path, rays));
    CHECK(rays.status()[0] == RayStatus::EventImpossible);
    const std::string lens = std::string(name) == "ideal lens" ? "IL" : "CL";
    CHECK(cs.surfaces()[rays.last_surface()[0]].id.str() == lens);
  }
}
