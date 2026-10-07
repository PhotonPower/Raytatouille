// Edge cases of the loss counting behind RayLosses and the warnings (ADR 0023; second review
// of #104: T1-T3). They test the internal counter directly: exactly half lost or a tie between
// two surfaces cannot be arranged reliably through the geometry of a traced system.

#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <string>
#include <vector>

#include "../src/common.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"

using rtt::analysis::detail::LossCounter;
using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::trace::RayStatus;

namespace {

/// m1/singlet_const: surfaces STO (0), L1.S1 (1), L1.S2 (2), IMG (3).
CompiledSystem singlet() {
  const rtt::material::MaterialLibrary lib;
  return rtt::compile::compile(
      rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m1/singlet_const.rtt.json"), lib);
}

bool warns(const LossCounter& c, const std::string& code) {
  for (const auto& d : c.warnings()) {
    if (d.code == code) return true;
  }
  return false;
}

}  // namespace

TEST_CASE("T1: exactly the threshold fraction lost does not warn", "[analysis][losses]") {
  const CompiledSystem cs = singlet();
  // 2 of 4 rays lost: exactly one half.
  const auto counter = [&](double threshold) {
    LossCounter c(cs, PathId{0}, threshold);
    c.add_ray(RayStatus::Alive, 3);
    c.add_ray(RayStatus::Alive, 3);
    c.add_ray(RayStatus::Vignetted, 1);
    c.add_ray(RayStatus::Vignetted, 1);
    return c;
  };
  REQUIRE_FALSE(warns(counter(0.5), "rays.lost"));  // "more than", not "at least"
  REQUIRE(warns(counter(std::nextafter(0.5, 0.0)), "rays.lost"));
  // 270 of 381 (the fraction of the small-lens test): no rounding of threshold * launched.
  LossCounter c(cs, PathId{0}, 270.0 / 381.0);
  for (int i = 0; i < 111; ++i) c.add_ray(RayStatus::Alive, 3);
  for (int i = 0; i < 270; ++i) c.add_ray(RayStatus::Vignetted, 1);
  REQUIRE_FALSE(warns(c, "rays.lost"));
}

TEST_CASE("T2: a tie between two loss surfaces goes to the lowest index", "[analysis][losses]") {
  const CompiledSystem cs = singlet();
  LossCounter c(cs, PathId{0}, 0.5);
  c.add_ray(RayStatus::Vignetted, 2);
  c.add_ray(RayStatus::Missed, 2);
  c.add_ray(RayStatus::Vignetted, 1);
  c.add_ray(RayStatus::Tir, 1);
  const auto losses = c.result();
  REQUIRE(losses.worst_surface == 1u);
  REQUIRE(losses.worst_surface_count == 2);
}

TEST_CASE("T3: a ray that ends Alive elsewhere than on the image counts as Vignetted",
          "[analysis][losses]") {
  const CompiledSystem cs = singlet();
  rtt::trace::RayBatch rays(2);
  rays.status()[0] = RayStatus::Alive;
  rays.last_surface()[0] = 3;  // on the image surface: arrived
  rays.status()[1] = RayStatus::Alive;
  rays.last_surface()[1] = 2;  // the path ended elsewhere
  LossCounter c(cs, PathId{0}, 0.5);
  c.add(rays, 3);
  const auto losses = c.result();
  REQUIRE(losses.count(RayStatus::Alive) == 1);
  REQUIRE(losses.count(RayStatus::Vignetted) == 1);
  REQUIRE(losses.worst_surface == 2u);
}
