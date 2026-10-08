// "The path has no stop" is one error type everywhere (ADR 0022, #86): rtt::compile::NoStopError,
// thrown by rtt::compile::require_stop before any computation.

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <string>

#include "rtt/analysis/chromatic.hpp"
#include "rtt/analysis/field.hpp"
#include "rtt/analysis/opd.hpp"
#include "rtt/analysis/spot.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/compile/errors.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/paraxial/paraxial.hpp"
#include "rtt/paraxial/seidel.hpp"
#include "rtt/trace/sources.hpp"

using rtt::compile::CompiledSystem;
using rtt::compile::NoStopError;
using rtt::compile::PathId;

namespace {

CompiledSystem compiled(const std::string& relative) {
  const rtt::material::MaterialLibrary lib;
  return rtt::compile::compile(
      rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/" + relative), lib);
}

}  // namespace

TEST_CASE("every analysis that needs the stop throws NoStopError", "[analysis][no-stop]") {
  // m1/paraboloid_mirror has no Stop element; m0/michelson has none on either arm.
  for (const char* file : {"m1/paraboloid_mirror.rtt.json", "m0/michelson.rtt.json"}) {
    INFO(file);
    const CompiledSystem cs = compiled(file);
    const PathId path{0};
    const std::array<std::uint16_t, 1> fields{0};
    using rtt::trace::Aiming;
    REQUIRE_THROWS_AS(rtt::trace::aim_ray(cs, path, 0, 0, 0.0, 0.0, Aiming::Real), NoStopError);
    REQUIRE_THROWS_AS(rtt::trace::aim_ray(cs, path, 0, 0, 0.0, 0.0, Aiming::Paraxial), NoStopError);
    REQUIRE_THROWS_AS(rtt::trace::make_rays(cs, path, fields, 0, rtt::trace::HexapolarPupil{}),
                      NoStopError);
    REQUIRE_THROWS_AS(rtt::analysis::spot(cs, path, 0, std::nullopt), NoStopError);
    REQUIRE_THROWS_AS(rtt::analysis::ray_fan(cs, path, 0, 0), NoStopError);
    REQUIRE_THROWS_AS(rtt::analysis::opd_map(cs, path, 0, 0), NoStopError);
    REQUIRE_THROWS_AS(rtt::analysis::opd_fan(cs, path, 0, 0), NoStopError);
    // A NoStopError is also a std::invalid_argument, the type sources threw before.
    REQUIRE_THROWS_AS(rtt::analysis::spot(cs, path, 0, std::nullopt), std::invalid_argument);
  }
  // Paraxial optics needs a rotationally symmetric path; of the two, only the paraboloid mirror
  // is one. first_order needs no stop (its pupils are empty), Seidel sums do.
  const CompiledSystem mirror = compiled("m1/paraboloid_mirror.rtt.json");
  REQUIRE_FALSE(rtt::paraxial::first_order(mirror, PathId{0}, 0).entrance_pupil.has_value());
  REQUIRE_THROWS_AS(rtt::paraxial::seidel(mirror, PathId{0}, 0), NoStopError);
}

TEST_CASE("NoStopError names the path and points at it", "[analysis][no-stop]") {
  const CompiledSystem cs = compiled("m0/michelson.rtt.json");
  try {
    (void)rtt::analysis::spot(cs, PathId{1}, 0, std::nullopt);
    FAIL("no NoStopError");
  } catch (const NoStopError& e) {
    REQUIRE(e.path_name() == cs.path(PathId{1}).name);
    REQUIRE(e.location() == "/paths/1");
  }
}

TEST_CASE("the colour and field analyses throw NoStopError as well", "[analysis][no-stop]") {
  // Second review of #98 (H2): they aim through rtt::trace::aim_ray. Argument checks of their
  // own (e.g. an on-axis field for distortion) may report first; no ray is traced before
  // require_stop.
  const CompiledSystem cs = compiled("m1/paraboloid_mirror.rtt.json");
  const rtt::model::Field off_axis{0.0, 0.1, 1.0};
  REQUIRE_THROWS_AS(rtt::analysis::distortion_at(cs, PathId{0}, off_axis, 0), NoStopError);
  REQUIRE_THROWS_AS(rtt::analysis::field_curvature_at(cs, PathId{0}, off_axis, 0), NoStopError);
  REQUIRE_THROWS_AS(rtt::analysis::lateral_colour(cs, PathId{0}, 0), NoStopError);
  REQUIRE_THROWS_AS(rtt::analysis::longitudinal_colour(cs, PathId{0}), NoStopError);
}
