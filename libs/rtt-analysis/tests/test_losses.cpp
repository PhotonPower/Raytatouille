// Where the rays of an analysis went (ADR 0023, #86): counts per final status and the surface
// where most lost rays ended, always part of the result.

#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <variant>
#include <vector>

#include "rtt/analysis/opd.hpp"
#include "rtt/analysis/paths.hpp"
#include "rtt/analysis/spot.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/trace/sources.hpp"

using rtt::analysis::RayLosses;
using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::model::System;
using rtt::trace::RayStatus;

namespace {

std::size_t count(const RayLosses& l, RayStatus s) {
  return l.by_status[static_cast<std::size_t>(s)];
}

/// m1/singlet_const (STO, L1.S1, L1.S2, IMG; EPD 20 mm), optionally with the front lens
/// aperture shrunk to `lens_radius_mm`.
CompiledSystem singlet(std::optional<double> lens_radius_mm = std::nullopt) {
  System s = rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m1/singlet_const.rtt.json");
  if (lens_radius_mm) {
    std::get<rtt::model::Element>(s.root.children[1].value).surfaces[0].aperture =
        rtt::model::CircularAperture{*lens_radius_mm, 0.0};
  }
  const rtt::material::MaterialLibrary lib;
  return rtt::compile::compile(s, lib);
}

/// The counts add up and agree with the arrived rays of the result.
void require_consistent(const RayLosses& l, std::size_t launched, std::size_t arrived) {
  REQUIRE(l.launched == launched);
  std::size_t sum = 0;
  for (const std::size_t n : l.by_status) sum += n;
  REQUIRE(sum == launched);
  REQUIRE(count(l, RayStatus::Alive) == arrived);
}

}  // namespace

TEST_CASE("without losses only Alive is counted and there is no loss surface",
          "[analysis][losses]") {
  const CompiledSystem cs = singlet();
  const auto spot = rtt::analysis::spot(cs, PathId{0}, 0, std::nullopt);
  require_consistent(spot.losses, spot.rays_launched, spot.rays_arrived);
  REQUIRE(spot.rays_arrived == spot.rays_launched);
  REQUIRE_FALSE(spot.losses.worst_surface.has_value());
  REQUIRE(spot.losses.worst_surface_count == 0);
}

TEST_CASE("rays clipped by the lens aperture are counted at that surface", "[analysis][losses]") {
  // EPD 20 mm on axis: the marginal rays meet L1.S1 near 10 mm, outside a 6 mm aperture
  // (hit, then Vignetted: last_surface is L1.S1, index 1).
  const CompiledSystem cs = singlet(6.0);
  const std::uint32_t lens = *cs.find_surface(rtt::model::SurfaceId("L1.S1"));

  SECTION("spot over all wavelengths") {
    const auto spot = rtt::analysis::spot(cs, PathId{0}, 0, std::nullopt);
    require_consistent(spot.losses, spot.rays_launched, spot.rays_arrived);
    const std::size_t lost = spot.rays_launched - spot.rays_arrived;
    REQUIRE(lost > 0);
    REQUIRE(count(spot.losses, RayStatus::Vignetted) == lost);
    REQUIRE(spot.losses.worst_surface == lens);
    REQUIRE(spot.losses.worst_surface_count == lost);
  }
  SECTION("ray fans: both fans, not the chief ray") {
    const auto fan = rtt::analysis::ray_fan(cs, PathId{0}, 0, 1);
    std::size_t arrived = 0;
    for (const auto* points : {&fan.tangential, &fan.sagittal}) {
      for (const auto& p : *points) arrived += p.status == RayStatus::Alive ? 1 : 0;
    }
    require_consistent(fan.losses, fan.tangential.size() + fan.sagittal.size(), arrived);
    REQUIRE(fan.losses.worst_surface == lens);
  }
  SECTION("OPD map and fan") {
    const auto map = rtt::analysis::opd_map(cs, PathId{0}, 0, 1);
    require_consistent(map.losses, map.points.size(), map.arrived);
    REQUIRE(map.losses.worst_surface == lens);
    const auto opd_fan = rtt::analysis::opd_fan(cs, PathId{0}, 0, 1);
    std::size_t arrived = 0;
    for (const auto* points : {&opd_fan.tangential, &opd_fan.sagittal}) {
      for (const auto& p : *points) arrived += p.status == RayStatus::Alive ? 1 : 0;
    }
    require_consistent(opd_fan.losses, opd_fan.tangential.size() + opd_fan.sagittal.size(),
                       arrived);
    REQUIRE(opd_fan.losses.worst_surface == lens);
  }
}

TEST_CASE("an Evanescent ray counts as lost at the surface where it stopped (#127)",
          "[analysis][losses]") {
  // ADR 0025, point 7: Evanescent stops a ray at its surface like Tir (last_surface is that
  // surface). The losses count it under its own status, index 7, and at that surface.
  const CompiledSystem cs = singlet();
  const std::vector<std::uint16_t> field{0};
  rtt::trace::RayBatch start =
      rtt::trace::make_rays(cs, PathId{0}, field, 0, rtt::trace::HexapolarPupil{1});
  REQUIRE(start.size() == 7);
  start.status()[1] = RayStatus::Evanescent;
  start.last_surface()[1] = 1;  // L1.S1
  const auto t = rtt::analysis::path_transmission(cs, PathId{0}, start);
  REQUIRE(t.losses.by_status.size() == 8);
  REQUIRE(count(t.losses, RayStatus::Evanescent) == 1);
  REQUIRE(t.losses.worst_surface == std::optional<std::uint32_t>{1});
  require_consistent(t.losses, 7, 6);
}
