// CompiledSystem::object_z(): the one place that turns object.distance into a position (#288).
// The distance is measured from the global origin along -z (ADR 0031, addendum #279), so the
// axial object point lies at z = -distance, bitwise, also with the distance bound to a row in
// several configurations.

#include <bit>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
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
