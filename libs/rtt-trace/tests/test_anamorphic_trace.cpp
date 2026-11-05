// A biconic with equal sections traces like the conic (#315, ADR 0042, acceptance 2): the
// intersection points of a bundle agree to 1e-12 mm on every surface. The biconic is
// intersected by the Newton iteration with its domain test (intersect_newton, in_domain), the
// conic by its closed form, so the points agree up to the convergence of the iteration.

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <string>
#include <variant>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/ray_paths.hpp"
#include "rtt/trace/sequential.hpp"

using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::math::Vec3;
using rtt::model::Element;
using rtt::model::Param;
using rtt::model::System;
using rtt::trace::RayBatch;
using rtt::trace::RayPaths;
using rtt::trace::RayStatus;
using rtt::trace::SequentialTracer;

namespace {

/// m1/singlet_const (stop at z = 0, L1.S1 at z = 5, image at z = 106.363) with `front` at
/// L1.S1.
System singlet(const rtt::model::BaseShape& front) {
  System s = rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m1/singlet_const.rtt.json");
  std::get<Element>(s.root.children[1].value).surfaces[0].shape.base = front;
  return s;
}

/// A grid of rays from z = -1 (in front of the stop, AGENTS.md) inside r <= 9.5, parallel to
/// the axis and in the oblique direction (0.05, 0.08, 1) normalized.
RayBatch bundle() {
  std::vector<Vec3> starts;
  for (int i = -6; i <= 6; ++i) {
    for (int j = -6; j <= 6; ++j) {
      const double x = 1.5 * i;
      const double y = 1.5 * j;
      if (x * x + y * y <= 9.5 * 9.5) starts.emplace_back(x, y, -1.0);
    }
  }
  const Vec3 oblique = Vec3(0.05, 0.08, 1.0).normalized();
  RayBatch rays(2 * starts.size());
  for (std::size_t k = 0; k < rays.size(); ++k) {
    const Vec3& p = starts[k % starts.size()];
    const Vec3 d = k < starts.size() ? Vec3::UnitZ() : oblique;
    rays.pos_x()[k] = p.x();
    rays.pos_y()[k] = p.y();
    rays.pos_z()[k] = p.z();
    rays.dir_x()[k] = d.x();
    rays.dir_y()[k] = d.y();
    rays.dir_z()[k] = d.z();
  }
  return rays;
}

RayPaths traced(const System& s) {
  const rtt::material::MaterialLibrary lib;
  const CompiledSystem cs = rtt::compile::compile(s, lib);
  RayBatch rays = bundle();
  RayPaths paths;
  const auto stats = SequentialTracer().trace(cs, PathId{0}, rays, paths);
  REQUIRE(stats.count(RayStatus::Alive) == rays.size());
  return paths;
}

}  // namespace

TEST_CASE("a biconic with equal sections hits the conic's points to 1e-12 mm",
          "[trace][anamorphic]") {
  for (const double k : {0.0, -0.6}) {
    INFO("k = " << k);
    const RayPaths a = traced(singlet(rtt::model::Conic{Param(51.68), Param(k)}));
    const RayPaths b =
        traced(singlet(rtt::model::Biconic{Param(51.68), Param(51.68), Param(k), Param(k)}));
    REQUIRE(a.ray_count() == b.ray_count());
    REQUIRE(a.slots == b.slots);
    REQUIRE(a.slots == 5);  // start, stop, L1.S1, L1.S2, image
    double spread = 0.0;    // content guard: the rays leave L1.S1 in different places
    for (std::size_t r = 0; r < a.ray_count(); ++r) {
      REQUIRE(a.count[r] == a.slots);
      REQUIRE(b.count[r] == b.slots);
      for (std::size_t s = 0; s < a.slots; ++s) {
        INFO("ray " << r << ", slot " << s);
        CHECK(b.status[r * b.slots + s] == RayStatus::Alive);
        CHECK((a.position_at(r, s) - b.position_at(r, s)).cwiseAbs().maxCoeff() <= 1e-12);
      }
      // On the curved surface the sag is not zero away from the axis.
      spread = std::max(spread, std::abs(a.position_at(r, 2).z() - 5.0));
    }
    CHECK(spread > 0.1);
  }
}
