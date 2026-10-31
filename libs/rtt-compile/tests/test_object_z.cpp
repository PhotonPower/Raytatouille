// CompiledSystem::object_z(): the one place that turns object.distance into a position (#288).
// The distance is measured from the global origin along -z (ADR 0031, addendum #279), so the
// axial object point lies at z = -distance, bitwise, also with the distance bound to a row in
// several configurations. With object.reference (ADR 0038) it lies the distance in front of the
// nominal vertex of that surface.

#include <bit>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#include "rtt/coating/catalog.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/compile/perturbation.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"

namespace {

bool same(double a, double b) {
  return std::bit_cast<std::uint64_t>(a) == std::bit_cast<std::uint64_t>(b);
}

rtt::model::System ideal_lens() {
  // Object 75 mm in front of the origin (ADR 0031, #178).
  return rtt::io::load_system(std::filesystem::path(RTT_REFERENCE_DIR) / "r2" /
                              "ideal_lens.rtt.json");
}

rtt::model::System referenced() {
  // Ideal lens at z = 20, object 75 mm in front of its vertex (ADR 0038): z_obj = 20 - 75 = -55.
  return rtt::io::load_system(std::filesystem::path(RTT_REFERENCE_DIR) / "r3" /
                              "object_reference.rtt.json");
}

/// The element "ideal lens" of `s` (root child 1).
rtt::model::Element& lens(rtt::model::System& s) {
  return std::get<rtt::model::Element>(s.root.children[1].value);
}

}  // namespace

TEST_CASE("object_z: the object point lies at -distance", "[compile][object]") {
  const rtt::material::MaterialLibrary materials;
  const rtt::model::System s = ideal_lens();
  REQUIRE(!s.object.at_infinity);
  REQUIRE(s.object.distance.value == 75.0);
  const rtt::compile::CompiledSystem cs = rtt::compile::compile(s, materials);
  CHECK(same(cs.object_z(), -75.0));
  CHECK(same(cs.object_z(), -cs.object().distance.value));
}

TEST_CASE("object_z: an object at infinity keeps -distance as stored", "[compile][object]") {
  // At infinity object_z() has no meaning and every caller checks at_infinity first; the value is
  // still the stored -distance, so nothing else can get in.
  const rtt::material::MaterialLibrary materials;
  const rtt::model::System s = rtt::io::load_system(std::filesystem::path(RTT_REFERENCE_DIR) /
                                                    "m7" / "decentred_lens.rtt.json");
  REQUIRE(s.object.at_infinity);
  const rtt::compile::CompiledSystem cs = rtt::compile::compile(s, materials);
  CHECK(same(cs.object_z(), -cs.object().distance.value));
}

TEST_CASE("object_z: a distance from a row follows each configuration", "[compile][object]") {
  // The distance bound to a row with one value per configuration: the resolved value of each
  // configuration (compile copies the object after resolving the table).
  const rtt::material::MaterialLibrary materials;
  rtt::model::System s = ideal_lens();
  s.configurations = {{"near"}, {"far"}};
  rtt::model::ParameterRow row;
  row.name = "OBJ";
  row.form = std::vector<double>{75.0, 112.5};
  s.parameters = {row};
  s.object.distance = rtt::model::Param::bound("OBJ");
  CHECK(same(rtt::compile::compile(s, materials, 0).object_z(), -75.0));
  CHECK(same(rtt::compile::compile(s, materials, 1).object_z(), -112.5));
}

TEST_CASE("object_z: in front of the nominal vertex of object.reference", "[compile][object]") {
  // ADR 0038, point 3: z_obj = z_S - distance, z_S the vertex z of the surface (exact here).
  const rtt::material::MaterialLibrary materials;
  const rtt::model::System s = referenced();
  REQUIRE(s.object.reference.has_value());
  const rtt::compile::CompiledSystem cs = rtt::compile::compile(s, materials);
  CHECK(same(cs.object_z(), -55.0));
  REQUIRE(cs.object_reference_z().has_value());
  CHECK(same(cs.object_reference_z().value(), 20.0));
  // The same object from the origin: distance 55, no reference; the same position, bitwise.
  rtt::model::System origin = s;
  origin.object.reference.reset();
  origin.object.distance = rtt::model::Param(55.0);
  const rtt::compile::CompiledSystem co = rtt::compile::compile(origin, materials);
  CHECK(same(co.object_z(), cs.object_z()));
  CHECK(!co.object_reference_z().has_value());
}

TEST_CASE("object_z: the object follows the surface in the model, not in the overlay",
          "[compile][object]") {
  const rtt::material::MaterialLibrary materials;
  const rtt::coating::CoatingLibrary coatings;
  // A model change (an optimizer variable, ADR 0038, point 5): the object follows.
  rtt::model::System moved = referenced();
  lens(moved).pose.position[2] = rtt::model::Param(25.0);
  CHECK(same(rtt::compile::compile(moved, materials).object_z(), 25.0 - 75.0));
  // The overlay of a tolerance trial (ADR 0036, point 1b): the nominal frame counts.
  const rtt::model::System s = referenced();
  rtt::compile::Perturbation p;
  p.nodes["ideal lens"].translation_mm = rtt::math::Vec3(0.5, 0.0, 3.0);
  const rtt::compile::CompiledSystem cs = rtt::compile::compile(s, materials, coatings, 0, p);
  CHECK(same(cs.object_z(), -55.0));
  // Content guard: the lens itself did move in the overlay.
  const std::optional<std::uint32_t> k = cs.find_surface(rtt::model::SurfaceId("IL"));
  REQUIRE(k.has_value());
  CHECK(cs.surfaces()[*k].to_global.translation().z() != 20.0);
}

TEST_CASE("object_z: Perturbation::object_reference_z replaces the vertex of the trial",
          "[compile][object]") {
  // ADR 0038, point 5: rtt-tolerance passes the nominal z_S; the trial's distance still counts.
  const rtt::material::MaterialLibrary materials;
  const rtt::coating::CoatingLibrary coatings;
  rtt::model::System s = referenced();
  lens(s).pose.position[2] = rtt::model::Param(25.0);  // e.g. a thickness item before the lens
  rtt::compile::Perturbation p;
  p.object_reference_z = {20.0};
  const rtt::compile::CompiledSystem cs = rtt::compile::compile(s, materials, coatings, 0, p);
  CHECK(same(cs.object_z(), -55.0));
  REQUIRE(cs.object_reference_z().has_value());
  CHECK(same(cs.object_reference_z().value(), 20.0));
  s.object.distance = rtt::model::Param(80.0);  // a tolerance on the distance moves the object
  CHECK(same(rtt::compile::compile(s, materials, coatings, 0, p).object_z(), 20.0 - 80.0));
}

TEST_CASE("object_z: input errors of object_reference_z", "[compile][object]") {
  const rtt::material::MaterialLibrary materials;
  const rtt::coating::CoatingLibrary coatings;
  rtt::compile::Perturbation p;
  p.object_reference_z = {20.0};
  CHECK_THROWS_AS(rtt::compile::compile(ideal_lens(), materials, coatings, 0, p),
                  std::invalid_argument);  // the system has no reference
  const rtt::model::System s = referenced();
  p.object_reference_z = {20.0, 20.0};  // one configuration
  CHECK_THROWS_AS(rtt::compile::compile(s, materials, coatings, 0, p), std::invalid_argument);
  p.object_reference_z = {std::numeric_limits<double>::infinity()};
  CHECK_THROWS_AS(rtt::compile::compile(s, materials, coatings, 0, p), std::invalid_argument);
}

TEST_CASE("object_z: the reference follows each configuration", "[compile][object]") {
  // The lens position from a row with one value per configuration: z_S per configuration.
  const rtt::material::MaterialLibrary materials;
  rtt::model::System s = referenced();
  s.configurations = {{"near"}, {"far"}};
  rtt::model::ParameterRow row;
  row.name = "ZL";
  row.form = std::vector<double>{20.0, 40.0};
  s.parameters = {row};
  lens(s).pose.position[2] = rtt::model::Param::bound("ZL");
  CHECK(same(rtt::compile::compile(s, materials, 0).object_z(), -55.0));
  CHECK(same(rtt::compile::compile(s, materials, 1).object_z(), -35.0));
}
