// Tracing relatively placed systems (ADR 0028, #163): the fold mirror of ADR 0028, point 4, and
// the calcite plate of the M4 acceptance placed relative to a rotated element.

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "relative_placement.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/sequential.hpp"

using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::material::MaterialLibrary;
using rtt::math::Vec3;
using rtt::model::Element;
using rtt::model::ElementKind;
using rtt::model::Param;
using rtt::model::Pose;
using rtt::model::PoseOrder;
using rtt::model::PoseReference;
using rtt::model::Surface;
using rtt::model::SurfaceId;
using rtt::model::System;
using rtt::trace::RayBatch;
using rtt::trace::RayStatus;
using rtt::trace::SequentialTracer;

namespace {

/// Rays set by hand: position and direction (as test_crystal_trace.cpp).
RayBatch rays_at(const std::vector<std::pair<Vec3, Vec3>>& start) {
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

Vec3 pos(const RayBatch& r, std::size_t i) {
  return {r.pos_x()[i], r.pos_y()[i], r.pos_z()[i]};
}
Vec3 dir(const RayBatch& r, std::size_t i) {
  return {r.dir_x()[i], r.dir_y()[i], r.dir_z()[i]};
}

RayBatch trace(const CompiledSystem& cs, const std::string& path, RayBatch rays) {
  const auto id = cs.find_path(path);
  REQUIRE(id.has_value());
  [[maybe_unused]] const auto stats = SequentialTracer().trace(cs, *id, rays);
  return rays;
}

}  // namespace

TEST_CASE("relative trace: the fold mirror equals its absolute counterpart", "[relative]") {
  // ADR 0028, point 4 (test with #163): the same rays through the relatively placed fold mirror
  // and its absolute counterpart (relative_placement.hpp: fold_mirror). The surface transforms
  // agree to 1e-12 mm and 1e-12 rad (test_relative_placement.cpp). A ray meets a surface moved
  // by 1e-12 mm and tilted by 1e-12 rad at most a few 1e-12 mm away (heights below 2 mm,
  // incidence below 45 deg); the direction changes by a few 1e-12 rad at the plane mirror and
  // at the lens (R = 40 mm, n = 1.5), and over the 75 mm to the image that is below 1e-9 mm.
  // Tolerance 1e-9 mm for position and OPL, 1e-10 for the direction.
  const MaterialLibrary lib;
  const CompiledSystem rel = rtt::compile::compile(rtt::model::test::fold_mirror(true, 20.0), lib);
  const CompiledSystem abs = rtt::compile::compile(rtt::model::test::fold_mirror(false, 20.0), lib);
  const std::vector<std::pair<Vec3, Vec3>> start{{Vec3(0.0, 0.0, 0.0), Vec3::UnitZ()},
                                                 {Vec3(1.0, 0.0, 0.0), Vec3::UnitZ()},
                                                 {Vec3(0.0, 1.5, 0.0), Vec3::UnitZ()},
                                                 {Vec3(-0.7, -1.2, 0.0), Vec3::UnitZ()}};
  const RayBatch a = trace(rel, "main", rays_at(start));
  const RayBatch b = trace(abs, "main", rays_at(start));
  for (std::size_t i = 0; i < start.size(); ++i) {
    INFO("ray " << i);
    REQUIRE(a.status()[i] == RayStatus::Alive);
    REQUIRE(b.status()[i] == RayStatus::Alive);
    CHECK((pos(a, i) - pos(b, i)).norm() <= 1e-9);
    CHECK((dir(a, i) - dir(b, i)).norm() <= 1e-10);
    CHECK(std::abs(a.opl()[i] - b.opl()[i]) <= 1e-9);
  }
  // Content: the image lies on the reflected axis (y = 55 at d = 20), the chief ray ends there
  // along +y, and the off-axis rays are refracted (not parallel to the axis).
  REQUIRE((pos(a, 0) - Vec3(0.0, 55.0, 50.0)).norm() <= 1e-9);
  REQUIRE((dir(a, 0) - Vec3::UnitY()).norm() <= 1e-12);
  REQUIRE(std::abs(dir(a, 1).x()) > 1e-3);
}

TEST_CASE("relative trace: the calcite plate placed relative to a rotated element", "[relative]") {
  // ADR 0028, acceptance: tests/reference/m4/calcite_walkoff.rtt.json with the plate P placed
  // relative to an element F in front of it. F (thin, z = 4) is rotated by 90 deg about z; P
  // follows with relative_to_preceding, rotate_first, rotation (0, 0, -90 deg), position
  // (0, 0, 6): P lies at z = 10 with the orientation of the file again (Rz(90) Rz(-90) = I up to
  // ~1e-16). The optic axis is given in element coordinates (ADR 0026, point 2) and turns with
  // the frame of P: with the parent as reference P would stay rotated by -90 deg, the axis in
  // the y-z plane, and the plate would start at z = 6. The
  // expected values are those of the M4 acceptance (test_crystal_trace.cpp, 1e-10), and the
  // compiled optic axis is that of the absolute file to 1e-12.
  const System original =
      rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m4/calcite_walkoff.rtt.json");
  System s = original;
  Surface f_surface;
  f_surface.id = SurfaceId("F.S");
  Element f{"F", ElementKind::ThinElement, Pose::along_z(4.0), std::nullopt, {f_surface}};
  f.pose.rotation_deg[2] = Param(90.0);
  auto& plate = std::get<Element>(s.root.children[0].value);
  plate.pose = Pose{};
  plate.pose.reference = PoseReference::RelativeToPreceding;
  plate.pose.order = PoseOrder::RotateFirst;
  plate.pose.rotation_deg[2] = Param(-90.0);
  plate.pose.position[2] = Param(6.0);
  s.root.children.insert(s.root.children.begin(), {f});
  const MaterialLibrary lib;
  const CompiledSystem cs = rtt::compile::compile(s, lib);
  const CompiledSystem reference = rtt::compile::compile(original, lib);
  const auto crystal_axis = [](const CompiledSystem& c) {
    for (const auto& m : c.media()) {
      if (m.optic_axis) return *m.optic_axis;
    }
    FAIL("no crystal medium");
    return Vec3(0.0, 0.0, 0.0);
  };
  REQUIRE((crystal_axis(cs) - crystal_axis(reference)).norm() <= 1e-12);
  REQUIRE(std::abs(crystal_axis(cs).x()) > 0.5);  // (1, 0, 1)/sqrt(2), not turned onto y
  // n_e(45 deg) and the walk-off from the M4 acceptance (Lam, Eq. (2.39) and p. 107):
  constexpr double kNo = 1.6584;
  constexpr double kNe = 1.4864;
  constexpr double kThickness = 2.0;
  const double n_e45 = std::sqrt(2.0 / (1.0 / (kNo * kNo) + 1.0 / (kNe * kNe)));
  const double rho45 = std::atan((kNo * kNo - kNe * kNe) / (kNo * kNo + kNe * kNe));
  for (const char* name : {"o", "e"}) {
    INFO(name);
    const RayBatch r = trace(cs, name, rays_at({{Vec3(0.0, 0.3, 0.0), Vec3::UnitZ()}}));
    REQUIRE(r.status()[0] == RayStatus::Alive);
    const bool e = std::string(name) == "e";
    const double offset = e ? -kThickness * std::tan(rho45) : 0.0;
    REQUIRE((pos(r, 0) - Vec3(offset, 0.3, 20.0)).norm() <= 1e-10);
    REQUIRE(std::abs(r.opl()[0] - (18.0 + (e ? n_e45 : kNo) * kThickness)) <= 1e-10);
  }
}
