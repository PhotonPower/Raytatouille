// object.medium (#282, ADR 0037): the medium of the object space is the medium before the first
// event of every path; the media of the elements (layout) still start in the environment.

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <string>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/material/material.hpp"
#include "test_support.hpp"

using rtt::compile::CompiledSystem;
using rtt::compile::CompileError;
using rtt::compile::PathId;
using rtt::material::MaterialLibrary;
using rtt::model::Param;
using rtt::model::System;

namespace {

/// Singlet with a finite object at 100 mm; the lens is made of CONST:1.6.
System finite_singlet() {
  System s = rtt::model::test::make_singlet();
  rtt::model::test::element(s, 1).material = "CONST:1.6";
  s.object.at_infinity = false;
  s.object.distance = Param(100.0);
  s.fields = {rtt::model::FieldType::ObjectHeight, {{0.0, 0.0, 1.0}}};
  return s;
}

std::uint32_t medium_index(const CompiledSystem& cs, const std::string& reference) {
  for (std::size_t i = 0; i < cs.media().size(); ++i) {
    if (cs.media()[i].reference == reference) return static_cast<std::uint32_t>(i);
  }
  FAIL("medium " << reference << " not found");
  return 0;
}

bool has_code(const CompiledSystem& cs, const std::string& code) {
  return std::any_of(cs.diagnostics().begin(), cs.diagnostics().end(),
                     [&](const rtt::model::Diagnostic& d) { return d.code == code; });
}

}  // namespace

TEST_CASE("without object.medium the object space is the environment", "[compile][object-medium]") {
  const MaterialLibrary lib;
  const CompiledSystem cs = rtt::compile::compile(finite_singlet(), lib);
  const auto& events = cs.path(PathId{0}).events;
  REQUIRE(events.front().medium_before == cs.environment_medium());
}

TEST_CASE("object.medium is the medium before the first event of every path",
          "[compile][object-medium]") {
  System s = finite_singlet();
  s.object.medium = "CONST:1.5";
  const MaterialLibrary lib;
  const CompiledSystem cs = rtt::compile::compile(s, lib);
  const std::uint32_t object = medium_index(cs, "CONST:1.5");
  const std::uint32_t lens = medium_index(cs, "CONST:1.6");
  REQUIRE(object != cs.environment_medium());
  // The environment stays media()[0].
  REQUIRE(cs.media()[cs.environment_medium()].reference == "AIR");
  const auto& events = cs.path(PathId{0}).events;
  // events: stop (transmit), L1.S1, L1.S2, detector
  REQUIRE(events.size() == 4);
  REQUIRE(events[0].medium_before == object);
  REQUIRE(events[0].medium_after == object);   // a stop keeps the medium
  REQUIRE(events[1].medium_before == object);  // object medium -> lens
  REQUIRE(events[1].medium_after == lens);
  REQUIRE(events[2].medium_before == lens);  // lens -> environment
  REQUIRE(events[2].medium_after == cs.environment_medium());
  REQUIRE(events[3].medium_before == cs.environment_medium());
}

TEST_CASE("the media of the elements (layout) still start in the environment",
          "[compile][object-medium]") {
  // medium_front/medium_back are per element and independent of the path: an object medium must
  // not change the surfaces of the elements that are not first on a path.
  System s = finite_singlet();
  s.object.medium = "CONST:1.5";
  const MaterialLibrary lib;
  const CompiledSystem cs = rtt::compile::compile(s, lib);
  REQUIRE(cs.surfaces()[0].medium_front == cs.environment_medium());  // STO
  REQUIRE(cs.surfaces()[1].medium_front == cs.environment_medium());  // L1.S1
  REQUIRE(cs.surfaces()[1].medium_back == medium_index(cs, "CONST:1.6"));
}

TEST_CASE("an unknown object medium is material.unknown at /object/medium",
          "[compile][object-medium]") {
  System s = finite_singlet();
  s.object.medium = "NOPE:X";
  const MaterialLibrary lib;
  try {
    (void)rtt::compile::compile(s, lib);
    FAIL("no CompileError");
  } catch (const CompileError& e) {
    REQUIRE(e.diagnostics().size() == 1);
    REQUIRE(e.diagnostics()[0].code == "material.unknown");
    REQUIRE(e.diagnostics()[0].location == "/object/medium");
  }
}

TEST_CASE("object-space NA at or above the index of the object medium is reported by compile",
          "[compile][object-medium]") {
  const MaterialLibrary lib;
  SECTION("NA 1.6 in glass of index 1.5: aperture.na_exceeds_index at the aperture value") {
    System s = finite_singlet();
    s.aperture = {rtt::model::SystemApertureType::ObjectSpaceNA, Param(1.6)};
    s.object.medium = "CONST:1.5";
    const CompiledSystem cs = rtt::compile::compile(s, lib);
    bool found = false;
    for (const auto& d : cs.diagnostics()) {
      if (d.code == "aperture.na_exceeds_index") {
        found = true;
        REQUIRE(d.location == "/aperture/value");
        REQUIRE(d.severity == rtt::model::Severity::Warning);
      }
    }
    REQUIRE(found);
    REQUIRE_FALSE(has_code(cs, "aperture.na_not_physical"));
  }
  SECTION("NA 0.9 in glass of index 1.5 is physical: no warning") {
    System s = finite_singlet();
    s.aperture = {rtt::model::SystemApertureType::ObjectSpaceNA, Param(0.9)};
    s.object.medium = "CONST:1.5";
    const CompiledSystem cs = rtt::compile::compile(s, lib);
    REQUIRE_FALSE(has_code(cs, "aperture.na_exceeds_index"));
    REQUIRE_FALSE(has_code(cs, "aperture.na_not_physical"));
  }
  SECTION("air: only the warning of validate, not a second one from compile") {
    System s = finite_singlet();
    s.aperture = {rtt::model::SystemApertureType::ObjectSpaceNA, Param(1.2)};
    const CompiledSystem cs = rtt::compile::compile(s, lib);
    REQUIRE(has_code(cs, "aperture.na_not_physical"));
    REQUIRE_FALSE(has_code(cs, "aperture.na_exceeds_index"));
  }
}
