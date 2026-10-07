// Error locations and the diagnostics channel of rtt-compile (ADR 0022, #86).

#include <catch2/catch_test_macros.hpp>
#include <stdexcept>
#include <string>
#include <type_traits>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/compile/errors.hpp"
#include "rtt/io/json_io.hpp"
#include "test_support.hpp"

using rtt::compile::compile;
using rtt::compile::CompiledSystem;
using rtt::compile::NoStopError;
using rtt::compile::PathId;
using rtt::compile::require_stop;
using rtt::material::MaterialLibrary;
using rtt::model::System;
using rtt::model::test::element;

static_assert(std::is_base_of_v<std::invalid_argument, NoStopError>);

namespace {

System load(const std::string& relative) {
  return rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/" + relative);
}

/// Singlet of rtt-model's test support in CONST:1.5168 (stop, L1, detector).
System singlet() {
  System s = rtt::model::test::make_singlet();
  element(s, 1).material = "CONST:1.5168";
  return s;
}

}  // namespace

TEST_CASE("every compiled surface knows its JSON pointer", "[compile][errors]") {
  const MaterialLibrary lib;
  System s = singlet();
  // One level of assembly: the detector moves into /root/children/2/children/0.
  rtt::model::Assembly group;
  group.name = "group";
  group.children.push_back(s.root.children[2]);
  s.root.children[2] = {group};
  const CompiledSystem cs = compile(s, lib);
  REQUIRE(cs.surfaces().size() == 4);
  REQUIRE(cs.surfaces()[0].location == "/root/children/0/surfaces/0");
  REQUIRE(cs.surfaces()[1].location == "/root/children/1/surfaces/0");
  REQUIRE(cs.surfaces()[2].location == "/root/children/1/surfaces/1");
  REQUIRE(cs.surfaces()[3].location == "/root/children/2/children/0/surfaces/0");
}

TEST_CASE("compile keeps the warnings of validate", "[compile][errors]") {
  const MaterialLibrary lib;
  REQUIRE(compile(singlet(), lib).diagnostics().empty());

  System s = singlet();
  element(s, 1).surfaces[0].shape.base =
      rtt::model::EvenAsphere{rtt::model::Param(51.68), rtt::model::Param(0.0), {}};
  const CompiledSystem cs = compile(s, lib);
  REQUIRE(cs.diagnostics().size() == 1);
  const rtt::model::Diagnostic& w = cs.diagnostics()[0];
  REQUIRE(w.severity == rtt::model::Severity::Warning);
  REQUIRE(w.code == "shape.asphere_without_coefficients");
  REQUIRE(w.location == "/root/children/1/surfaces/0/shape/base/coefficients");
}

TEST_CASE("require_stop accepts a path with a stop and names the path without one",
          "[compile][errors]") {
  const MaterialLibrary lib;
  REQUIRE_NOTHROW(require_stop(compile(singlet(), lib), PathId{0}));

  SECTION("no stop element in the system") {
    const CompiledSystem cs = compile(load("m1/paraboloid_mirror.rtt.json"), lib);
    try {
      require_stop(cs, PathId{0});
      FAIL("no NoStopError");
    } catch (const NoStopError& e) {
      REQUIRE(e.path_name() == cs.path(PathId{0}).name);
      REQUIRE(e.location() == "/paths/0");
      REQUIRE(std::string(e.what()).find("'" + cs.path(PathId{0}).name + "'") != std::string::npos);
    }
  }
  SECTION("a stop element that the path does not visit") {
    System s = singlet();
    s.paths.push_back({"lens only",
                       false,
                       {{rtt::model::SurfaceId("L1.S1"), rtt::model::EventKind::Refract, 0},
                        {rtt::model::SurfaceId("L1.S2"), rtt::model::EventKind::Refract, 0}}});
    const CompiledSystem cs = compile(s, lib);
    REQUIRE_NOTHROW(require_stop(cs, PathId{0}));
    try {
      require_stop(cs, PathId{1});
      FAIL("no NoStopError");
    } catch (const NoStopError& e) {
      REQUIRE(e.path_name() == "lens only");
      REQUIRE(e.location() == "/paths/1");
    }
  }
  SECTION("an unknown path id is out of range") {
    REQUIRE_THROWS_AS(require_stop(compile(singlet(), lib), PathId{7}), std::out_of_range);
  }
}
