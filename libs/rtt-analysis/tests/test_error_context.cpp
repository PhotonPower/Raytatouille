// AnalysisError carries the lost ray's surface, JSON pointer, status, field and wavelength
// (ADR 0022, #86).

#include <catch2/catch_test_macros.hpp>
#include <string>
#include <variant>
#include <vector>

#include "rtt/analysis/opd.hpp"
#include "rtt/analysis/spot.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"

using rtt::analysis::AnalysisError;
using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::model::System;
using rtt::trace::RayStatus;

namespace {

template <typename F>
AnalysisError analysis_error(F&& f) {
  try {
    f();
  } catch (const AnalysisError& e) {
    return e;
  }
  FAIL("no AnalysisError");
  return AnalysisError("unreachable");
}

/// m1/singlet_const with a field point at 80 degree (index 3): its chief ray passes the stop
/// centre (STO, /root/children/0/surfaces/0) and then misses the lens surface L1.S1, so it is
/// Missed with STO as its last surface.
CompiledSystem steep_singlet() {
  System s = rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m1/singlet_const.rtt.json");
  s.fields.points.push_back({0.0, 80.0, 1.0});
  REQUIRE(s.fields.points.size() == 4);
  const rtt::material::MaterialLibrary lib;
  return rtt::compile::compile(s, lib);
}

void require_missed_chief(const AnalysisError& e, const CompiledSystem& cs) {
  INFO(e.what() << "; surface " << (e.surface() ? e.surface()->str() : std::string("-"))
                << ", status " << (e.ray_status() ? static_cast<int>(*e.ray_status()) : -1));
  REQUIRE(e.surface() == rtt::model::SurfaceId("STO"));
  REQUIRE(e.location() == "/root/children/0/surfaces/0");
  REQUIRE(e.ray_status() == RayStatus::Missed);
  REQUIRE(e.field() == 3);
  REQUIRE(e.wavelength() == cs.reference_wavelength());
}

}  // namespace

TEST_CASE("a lost chief ray gives surface, location, status, field and wavelength",
          "[analysis][errors]") {
  const CompiledSystem cs = steep_singlet();
  require_missed_chief(
      analysis_error([&] { (void)rtt::analysis::spot(cs, PathId{0}, 3, std::nullopt); }), cs);
  require_missed_chief(analysis_error([&] { (void)rtt::analysis::ray_fan(cs, PathId{0}, 3, 1); }),
                       cs);
  require_missed_chief(analysis_error([&] { (void)rtt::analysis::opd_map(cs, PathId{0}, 3, 1); }),
                       cs);
}

TEST_CASE("a vignetted chief ray names the surface where it stopped", "[analysis][errors]") {
  // Lens aperture shrunk to 1 mm: the chief ray of 20 degree (field index 3) passes the stop
  // centre and hits L1.S1 about 1.8 mm off axis, outside the aperture.
  System s = rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m1/singlet_const.rtt.json");
  s.fields.points.push_back({0.0, 20.0, 1.0});
  std::get<rtt::model::Element>(s.root.children[1].value).surfaces[0].aperture =
      rtt::model::CircularAperture{1.0, 0.0};
  const rtt::material::MaterialLibrary lib;
  const CompiledSystem cs = rtt::compile::compile(s, lib);
  const AnalysisError e =
      analysis_error([&] { (void)rtt::analysis::spot(cs, PathId{0}, 3, std::nullopt); });
  REQUIRE(e.surface() == rtt::model::SurfaceId("L1.S1"));
  REQUIRE(e.location() == "/root/children/1/surfaces/0");
  REQUIRE(e.ray_status() == RayStatus::Vignetted);
  REQUIRE(e.field() == 3);
}

TEST_CASE("an AnalysisError without a ray has no context", "[analysis][errors]") {
  const AnalysisError e = analysis_error(
      [] { (void)rtt::analysis::spot_statistics(std::vector<rtt::analysis::SpotPoint>{}, {}); });
  REQUIRE_FALSE(e.surface().has_value());
  REQUIRE_FALSE(e.location().has_value());
  REQUIRE_FALSE(e.ray_status().has_value());
  REQUIRE_FALSE(e.field().has_value());
  REQUIRE_FALSE(e.wavelength().has_value());
}
