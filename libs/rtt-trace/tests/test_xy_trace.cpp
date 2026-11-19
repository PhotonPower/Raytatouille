// Real rays through an XY polynomial (#316, ADR 0051, point 4; AGENTS.md, checklist point 6):
// c_x2 = c_y2 = R^2 / (2 R_p) on a plane is mathematically the paraboloid z = (x^2 + y^2) / (2
// R_p), the conic with radius R_p and k = -1 (numerically (x/R)^2 c and the conic formula round
// differently, so no bitwise agreement is expected). The stack is intersected by the Newton
// iteration, the conic by its closed form (an independent computation).
//
// Bounds (derived): Newton stops at |F| <= kNewtonTolerance = 1e-12 mm (rtt/geom/intersect.hpp),
// F(t) = z(t) - sag(x(t), y(t)), so the hit is off by at most 1e-12 / min |F'| along the ray, with
// F' = d_z - grad sag . d_xy >= 0.9955 - 0.19 * 0.095 > 0.97 here (|grad sag| <= r / R_p = 0.19,
// the oblique d_z above 0.9955 and |d_xy| below 0.095): eps_s = 1.1e-12 mm at L1.S1. The normal
// there tilts by at most eps_s / R_p (the curvature 1 / R_p bounds the change of the gradient), the
// refracted direction by at most 2 eps_s / R_p, and the plane L1.S2 enlarges a direction error by
// at most n = 1.5: eps_d = 3 eps_s / R_p = 7e-14. Over the at most 105 mm to the image the position
// error grows by eps_d L = 7e-12. Rounding (a few ulp of 100 mm, 1e-13) is below that. The test
// uses positions 1e-11 mm and directions 1e-13.

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
System singlet(const rtt::model::ShapeStack& front) {
  System s = rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m1/singlet_const.rtt.json");
  std::get<Element>(s.root.children[1].value).surfaces[0].shape = front;
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

TEST_CASE("an XY polynomial paraboloid traces like the conic k = -1", "[trace][xy]") {
  constexpr double kR = 10.0;   // normalization radius, mm
  constexpr double kRp = 50.0;  // paraboloid vertex radius, mm
  const double c = kR * kR / (2.0 * kRp);
  rtt::model::ShapeStack xy;
  xy.terms = {rtt::model::XyPolynomial{Param(kR),
                                       {Param(0.0), Param(0.0), Param(c), Param(0.0), Param(c)}}};
  rtt::model::ShapeStack conic;
  conic.base = rtt::model::Conic{Param(kRp), Param(-1.0)};
  const RayPaths a = traced(singlet(conic));
  const RayPaths b = traced(singlet(xy));
  REQUIRE(a.ray_count() == b.ray_count());
  REQUIRE(a.slots == b.slots);
  REQUIRE(a.slots == 5);  // start, stop, L1.S1, L1.S2, image
  double spread = 0.0;    // content guard: the rays meet L1.S1 off its vertex plane
  for (std::size_t r = 0; r < a.ray_count(); ++r) {
    REQUIRE(a.count[r] == a.slots);
    REQUIRE(b.count[r] == b.slots);
    for (std::size_t s = 0; s < a.slots; ++s) {
      INFO("ray " << r << ", slot " << s);
      CHECK(b.status[r * b.slots + s] == RayStatus::Alive);
      CHECK((a.position_at(r, s) - b.position_at(r, s)).cwiseAbs().maxCoeff() <= 1e-11);
      CHECK((a.direction_at(r, s) - b.direction_at(r, s)).cwiseAbs().maxCoeff() <= 1e-13);
    }
    spread = std::max(spread, std::abs(a.position_at(r, 2).z() - 5.0));
  }
  CHECK(spread > 0.5);
}
