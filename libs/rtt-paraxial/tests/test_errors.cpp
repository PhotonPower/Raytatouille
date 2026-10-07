// ParaxialError carries the surface and the JSON pointer where they are known (ADR 0022, #86).

#include <catch2/catch_test_macros.hpp>
#include <string>
#include <variant>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/paraxial/paraxial.hpp"
#include "rtt/paraxial/seidel.hpp"

using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::model::Element;
using rtt::model::Param;
using rtt::model::System;
using rtt::paraxial::ParaxialError;

namespace {

/// m1/singlet_const: stop STO (/root/children/0), lens L1 with L1.S1 and L1.S2
/// (/root/children/1), image IMG (/root/children/2).
System singlet() {
  return rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m1/singlet_const.rtt.json");
}

CompiledSystem compiled(const System& s) {
  const rtt::material::MaterialLibrary lib;
  return rtt::compile::compile(s, lib);
}

Element& element(System& s, std::size_t i) {
  return std::get<Element>(s.root.children[i].value);
}

/// Runs f and returns the ParaxialError it throws.
template <typename F>
ParaxialError paraxial_error(F&& f) {
  try {
    f();
  } catch (const ParaxialError& e) {
    return e;
  }
  FAIL("no ParaxialError");
  return ParaxialError("unreachable");
}

}  // namespace

TEST_CASE("ParaxialError names the surface that breaks the symmetry", "[paraxial][errors]") {
  System s = singlet();
  element(s, 1).pose.rotation_deg[0] = Param(1e-3);
  const CompiledSystem cs = compiled(s);
  const ParaxialError e =
      paraxial_error([&] { (void)rtt::paraxial::first_order(cs, PathId{0}, 0); });
  REQUIRE(e.surface() == rtt::model::SurfaceId("L1.S1"));
  REQUIRE(e.location() == "/root/children/1/surfaces/0");
}

TEST_CASE("ParaxialError names a stop without circular aperture", "[paraxial][errors]") {
  System s = singlet();
  element(s, 0).surfaces[0].aperture = rtt::model::RectangularAperture{5.0, 5.0};
  const CompiledSystem cs = compiled(s);
  const ParaxialError e =
      paraxial_error([&] { (void)rtt::paraxial::first_order(cs, PathId{0}, 0); });
  REQUIRE(e.surface() == rtt::model::SurfaceId("STO"));
  REQUIRE(e.location() == "/root/children/0/surfaces/0");
}

TEST_CASE("Seidel field errors point at the field data", "[paraxial][errors]") {
  SECTION("field angle outside (-90, 90) degree: the field point") {
    System s = singlet();
    s.fields.points.push_back({0.0, 95.0, 1.0});
    const CompiledSystem cs = compiled(s);
    const ParaxialError e = paraxial_error([&] { (void)rtt::paraxial::seidel(cs, PathId{0}, 0); });
    REQUIRE_FALSE(e.surface().has_value());
    REQUIRE(e.location() == "/fields/points/" + std::to_string(s.fields.points.size() - 1));
  }
  SECTION("object height with the object at infinity: the field type") {
    System s = singlet();
    s.fields = {rtt::model::FieldType::ObjectHeight, {{0.0, 1.0, 1.0}}};
    const CompiledSystem cs = compiled(s);
    const ParaxialError e = paraxial_error([&] { (void)rtt::paraxial::seidel(cs, PathId{0}, 0); });
    REQUIRE_FALSE(e.surface().has_value());
    REQUIRE(e.location() == "/fields/type");
  }
}

TEST_CASE("ParaxialError without a place in the file has neither", "[paraxial][errors]") {
  const CompiledSystem cs = compiled(singlet());
  const ParaxialError e =
      paraxial_error([&] { (void)rtt::paraxial::first_order(cs, PathId{5}, 0); });
  REQUIRE_FALSE(e.surface().has_value());
  REQUIRE_FALSE(e.location().has_value());
}
