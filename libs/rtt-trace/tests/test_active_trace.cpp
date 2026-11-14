// Reference case 1 of #380 (ADR 0049, point 2): a system with an inactive plane-parallel window
// traces bit for bit like the same system without the window.

#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/ray_paths.hpp"
#include "rtt/trace/sequential.hpp"
#include "test_support.hpp"

using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::material::MaterialLibrary;
using rtt::math::Vec3;
using rtt::model::CircularAperture;
using rtt::model::Element;
using rtt::model::ElementKind;
using rtt::model::Node;
using rtt::model::Pose;
using rtt::model::Surface;
using rtt::model::SurfaceId;
using rtt::model::System;
using rtt::trace::RayBatch;
using rtt::trace::RayStatus;
using rtt::trace::SequentialTracer;

namespace {

/// The reference singlet in CONST:1.5168; every node is placed absolutely (make_singlet).
System singlet() {
  System s = rtt::model::test::make_singlet();
  rtt::model::test::element(s, 1).material = "CONST:1.5168";
  return s;
}

/// The singlet with a plane-parallel window (CONST:1.5, 2 mm) between stop and lens, at z = 1 mm,
/// switched off.
System with_inactive_window() {
  System s = singlet();
  Surface a;
  a.id = SurfaceId("W.S1");
  a.aperture = CircularAperture{12.7, 0.0};
  Surface b = a;
  b.id = SurfaceId("W.S2");
  b.pose = Pose::along_z(2.0);
  Element window{"window", ElementKind::Plate, Pose::along_z(1.0), "CONST:1.5", {a, b}};
  window.active = false;
  s.root.children.insert(s.root.children.begin() + 1, Node{window});
  return s;
}

/// Rays before the stop (z = -1 mm, AGENTS.md): on axis, in the pupil and inclined, one of them
/// outside the stop aperture (radius 10 mm), so that it ends vignetted.
RayBatch start_rays() {
  const std::vector<std::pair<Vec3, Vec3>> start{
      {Vec3(0.0, 0.0, -1.0), Vec3::UnitZ()},
      {Vec3(0.0, 7.0, -1.0), Vec3::UnitZ()},
      {Vec3(-4.0, 3.0, -1.0), Vec3(0.02, -0.05, 1.0).normalized()},
      {Vec3(2.5, -6.0, -1.0), Vec3(0.0, 0.0874, 1.0).normalized()},
      {Vec3(0.0, 11.0, -1.0), Vec3::UnitZ()}};
  RayBatch rays(start.size());
  for (std::size_t i = 0; i < start.size(); ++i) {
    rays.pos_x()[i] = start[i].first.x();
    rays.pos_y()[i] = start[i].first.y();
    rays.pos_z()[i] = start[i].first.z();
    rays.dir_x()[i] = start[i].second.x();
    rays.dir_y()[i] = start[i].second.y();
    rays.dir_z()[i] = start[i].second.z();
  }
  return rays;
}

}  // namespace

TEST_CASE("active: an inactive window traces bit for bit like no window (reference case 1)",
          "[active]") {
  const MaterialLibrary lib;
  const CompiledSystem a = rtt::compile::compile(with_inactive_window(), lib);
  const CompiledSystem b = rtt::compile::compile(singlet(), lib);
  // A keeps the surfaces of the window, so the indices differ; compare by id and name.
  REQUIRE(a.surfaces().size() == b.surfaces().size() + 2);
  const auto& ea = a.path(PathId{0}).events;
  const auto& eb = b.path(PathId{0}).events;
  REQUIRE(ea.size() == eb.size());
  for (std::size_t k = 0; k < ea.size(); ++k) {
    INFO("event " << k);
    CHECK(a.surfaces()[ea[k].surface].id == b.surfaces()[eb[k].surface].id);
    CHECK(ea[k].kind == eb[k].kind);
    CHECK(a.media()[ea[k].medium_before].reference == b.media()[eb[k].medium_before].reference);
    CHECK(a.media()[ea[k].medium_after].reference == b.media()[eb[k].medium_after].reference);
  }

  RayBatch ra = start_rays();
  RayBatch rb = start_rays();
  rtt::trace::RayPaths pa;
  rtt::trace::RayPaths pb;
  [[maybe_unused]] const auto sa = SequentialTracer().trace(a, PathId{0}, ra, pa);
  [[maybe_unused]] const auto sb = SequentialTracer().trace(b, PathId{0}, rb, pb);
  // Content guard (AGENTS.md): the rays arrive except the one outside the stop, and they differ.
  CHECK(ra.status()[0] == RayStatus::Alive);
  CHECK(ra.status()[3] == RayStatus::Alive);
  CHECK(ra.status()[4] != RayStatus::Alive);
  CHECK(ra.pos_y()[1] != ra.pos_y()[3]);
  for (std::size_t i = 0; i < ra.size(); ++i) {
    INFO("ray " << i);
    CHECK(ra.status()[i] == rb.status()[i]);
    CHECK(ra.pos_x()[i] == rb.pos_x()[i]);
    CHECK(ra.pos_y()[i] == rb.pos_y()[i]);
    CHECK(ra.pos_z()[i] == rb.pos_z()[i]);
    CHECK(ra.dir_x()[i] == rb.dir_x()[i]);
    CHECK(ra.dir_y()[i] == rb.dir_y()[i]);
    CHECK(ra.dir_z()[i] == rb.dir_z()[i]);
    CHECK(ra.opl()[i] == rb.opl()[i]);
    CHECK(ra.weight()[i] == rb.weight()[i]);
  }
  // The recorded paths agree at every valid slot too (the slots after a loss are NaN).
  REQUIRE(pa.slots == pb.slots);
  for (std::size_t i = 0; i < ra.size(); ++i) {
    REQUIRE(pa.count[i] == pb.count[i]);
    for (std::size_t k = 0; k < pa.count[i]; ++k) {
      INFO("ray " << i << ", slot " << k);
      CHECK(pa.position_at(i, k) == pb.position_at(i, k));
      CHECK(pa.direction_at(i, k) == pb.direction_at(i, k));
    }
  }
}
