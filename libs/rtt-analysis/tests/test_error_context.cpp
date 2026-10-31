// AnalysisError carries the lost ray's surface, JSON pointer, status, field and wavelength
// (ADR 0022, #86).

#include <catch2/catch_test_macros.hpp>
#include <string>
#include <variant>
#include <vector>

#include "rtt/analysis/chromatic.hpp"
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

namespace {

/// A singlet with R1 = 20 mm and EPD 50 mm: the marginal rays above 20 mm miss the sphere, so
/// real aiming cannot reach the outer targets on the stop 10 mm behind the lens (paraxial stop
/// radius about 20 mm). Stop radius `stop`: 3 mm puts those targets outside the stop, 50 mm
/// inside (#310).
CompiledSystem unreachable_rim(double stop) {
  const std::string json = R"json({
  "schema_version": "0.7.0",
  "name": "unreachable rim, issue 310",
  "units": {"length": "mm", "wavelength": "um"},
  "wavelengths": [{"um": 0.4861}, {"um": 0.5876, "reference": true}, {"um": 0.6563}],
  "aperture": {"type": "epd", "value": 50.0},
  "fields": {"points": [{}]},
  "root": {"type": "assembly", "name": "system", "children": [
    {"type": "lens", "name": "L1", "material": "CONST:1.5168",
     "surfaces": [{"id": "L1.S1", "shape": {"base": {"type": "conic", "radius": 20.0}},
                   "aperture": {"type": "circular", "radius": 19.9}},
                  {"id": "L1.S2", "pose": {"position": [0.0, 0.0, 5.0]},
                   "aperture": {"type": "circular", "radius": 19.9}}]},
    {"type": "stop", "name": "stop", "pose": {"position": [0.0, 0.0, 10.0]},
     "surfaces": [{"id": "STO", "aperture": {"type": "circular", "radius": )json" +
                           std::to_string(stop) + R"json(}}]},
    {"type": "detector", "name": "image", "pose": {"position": [0.0, 0.0, 40.0]},
     "surfaces": [{"id": "IMG", "aperture": {"type": "circular", "radius": 50.0}}]}]},
  "paths": [{"name": "main", "events": "auto"}]
})json";
  const rtt::material::MaterialLibrary lib;
  return rtt::compile::compile(rtt::io::parse_system(json), lib);
}

}  // namespace

TEST_CASE("a zone ray aimed outside the stop names the stop as Vignetted (#310)",
          "[analysis][errors]") {
  // ADR 0007, addendum #310: real aiming fails for the zone ray p = (0, 1), whose target lies
  // outside the 3 mm stop; the ray ends Vignetted at the stop, and the AnalysisError of the
  // longitudinal colour names it. With the 50 mm stop the target lies inside: NoConvergence
  // without a surface, as before (content guard: the aim fails in both systems).
  const CompiledSystem wide = unreachable_rim(50.0);
  const AnalysisError open =
      analysis_error([&] { (void)rtt::analysis::longitudinal_colour(wide, PathId{0}); });
  CHECK(open.ray_status() == RayStatus::NoConvergence);
  CHECK_FALSE(open.surface().has_value());

  const CompiledSystem narrow = unreachable_rim(3.0);
  const AnalysisError e =
      analysis_error([&] { (void)rtt::analysis::longitudinal_colour(narrow, PathId{0}); });
  INFO(e.what());
  CHECK(e.ray_status() == RayStatus::Vignetted);
  CHECK(e.surface() == rtt::model::SurfaceId("STO"));
  CHECK(e.location() == "/root/children/1/surfaces/0");

  // The spot through make_rays: the failed outer rays count as Vignetted at the stop, none as
  // NoConvergence (with the wide stop some do).
  const auto lost = [](const CompiledSystem& cs) {
    return rtt::analysis::spot(cs, PathId{0}, 0, std::nullopt).losses;
  };
  const rtt::analysis::RayLosses open_losses = lost(wide);
  CHECK(open_losses.count(RayStatus::NoConvergence) > 0);
  const rtt::analysis::RayLosses losses = lost(narrow);
  CHECK(losses.count(RayStatus::NoConvergence) == 0);
  CHECK(losses.count(RayStatus::Vignetted) >= open_losses.count(RayStatus::NoConvergence));
  REQUIRE(losses.worst_surface.has_value());
  CHECK(narrow.surfaces()[*losses.worst_surface].id == rtt::model::SurfaceId("STO"));
}
