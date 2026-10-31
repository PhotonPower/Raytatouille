// object.reference (#288, ADR 0038): an object given in front of the vertex of a surface is the
// same object as one given from the origin with the distance plus the vertex z. Rays, aiming and
// the matching of the aiming data (ADR 0007, addendum #255) only see object_z().

#include <bit>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <filesystem>
#include <stdexcept>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/sequential.hpp"
#include "rtt/trace/sources.hpp"

using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::model::Param;
using rtt::model::System;
using rtt::trace::Aiming;
using rtt::trace::RayBatch;

namespace {

bool same(double a, double b) {
  return std::bit_cast<std::uint64_t>(a) == std::bit_cast<std::uint64_t>(b);
}

/// Ideal lens at z = 20, object 75 mm in front of its vertex: z_obj = -55.
System referenced() {
  return rtt::io::load_system(std::filesystem::path(RTT_REFERENCE_DIR) / "r3" /
                              "object_reference.rtt.json");
}

/// The same object from the origin, `distance` in front of z = 0.
System from_origin(double distance) {
  System s = referenced();
  s.object.reference.reset();
  s.object.distance = Param(distance);
  return s;
}

RayBatch traced(const CompiledSystem& cs) {
  const std::vector<std::uint16_t> fields{0, 1};
  RayBatch rays =
      rtt::trace::make_rays(cs, PathId{0}, fields, 0, rtt::trace::HexapolarPupil{3}, Aiming::Real);
  [[maybe_unused]] const auto stats = rtt::trace::SequentialTracer().trace(cs, PathId{0}, rays);
  return rays;
}

}  // namespace

TEST_CASE("object.reference: the rays equal those from the origin, bitwise",
          "[sources][object-reference]") {
  const rtt::material::MaterialLibrary lib;
  const CompiledSystem a = rtt::compile::compile(referenced(), lib);
  const CompiledSystem b = rtt::compile::compile(from_origin(55.0), lib);
  REQUIRE(same(a.object_z(), b.object_z()));
  const RayBatch ra = traced(a);
  const RayBatch rb = traced(b);
  REQUIRE(ra.size() == rb.size());
  REQUIRE(ra.size() > 10);
  for (std::size_t i = 0; i < ra.size(); ++i) {
    INFO("ray " << i);
    CHECK(ra.status()[i] == rb.status()[i]);
    CHECK(same(ra.pos_x()[i], rb.pos_x()[i]));
    CHECK(same(ra.pos_y()[i], rb.pos_y()[i]));
    CHECK(same(ra.pos_z()[i], rb.pos_z()[i]));
    CHECK(same(ra.dir_y()[i], rb.dir_y()[i]));
    CHECK(same(ra.opl()[i], rb.opl()[i]));
  }
  // Content guard: the object lies at z = -55, not at -75 (the distance alone from the origin).
  const RayBatch rc = traced(rtt::compile::compile(from_origin(75.0), lib));
  CHECK(!same(rc.dir_y()[ra.size() - 1], ra.dir_y()[ra.size() - 1]));
}

TEST_CASE("object.reference: the aiming data match on object_z", "[sources][object-reference]") {
  // ADR 0038, point 6: the signature compares object_z(); a reference and the origin with the
  // same z match, another z does not.
  const rtt::material::MaterialLibrary lib;
  const CompiledSystem a = rtt::compile::compile(referenced(), lib);
  const CompiledSystem b = rtt::compile::compile(from_origin(55.0), lib);
  const rtt::trace::AimReference ref = rtt::trace::aim_reference(b, PathId{0}, 0);
  const rtt::trace::AimedRay with =
      rtt::trace::aim_ray(a, PathId{0}, 1, 0, 0.0, 1.0, Aiming::Real, ref);
  const rtt::trace::AimedRay plain = rtt::trace::aim_ray(a, PathId{0}, 1, 0, 0.0, 1.0);
  CHECK(same(with.ray.pos.y(), plain.ray.pos.y()));
  CHECK(same(with.ray.dir.y(), plain.ray.dir.y()));
  const CompiledSystem c = rtt::compile::compile(from_origin(60.0), lib);
  CHECK_THROWS_AS(rtt::trace::aim_ray(c, PathId{0}, 1, 0, 0.0, 1.0, Aiming::Real, ref),
                  std::invalid_argument);
}
