// Warnings of the analyses (ADR 0023, #86): rays.lost above a threshold, stop.clips_beam when
// rays end vignetted at the stop. Every analysis code of the registry has a case here.

#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <limits>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#include "rtt/analysis/opd.hpp"
#include "rtt/analysis/reports.hpp"
#include "rtt/analysis/spot.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/diagnostics/codes.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"

using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::model::Diagnostic;
using rtt::model::System;

namespace {

/// m1/singlet_const: stop STO (radius 10 mm, /root/children/0/surfaces/0) in front of the lens
/// L1 (L1.S1 /root/children/1/surfaces/0, radius 12.7 mm), EPD 20 mm.
System singlet() {
  return rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m1/singlet_const.rtt.json");
}

CompiledSystem compiled(const System& s) {
  const rtt::material::MaterialLibrary lib;
  return rtt::compile::compile(s, lib);
}

void lens_radius(System& s, double radius_mm) {
  std::get<rtt::model::Element>(s.root.children[1].value).surfaces[0].aperture =
      rtt::model::CircularAperture{radius_mm, 0.0};
}

/// The warnings of all four analyses of field `field` with the given threshold.
std::vector<std::vector<Diagnostic>> all_warnings(const CompiledSystem& cs,
                                                  std::uint16_t field,
                                                  double threshold = 0.5) {
  rtt::analysis::SpotOptions spot;
  spot.lost_warning_fraction = threshold;
  rtt::analysis::FanOptions fan;
  fan.lost_warning_fraction = threshold;
  rtt::analysis::OpdOptions opd;
  opd.lost_warning_fraction = threshold;
  return {rtt::analysis::spot(cs, PathId{0}, field, std::nullopt, spot).warnings,
          rtt::analysis::ray_fan(cs, PathId{0}, field, 1, fan).warnings,
          rtt::analysis::opd_map(cs, PathId{0}, field, 1, opd).warnings,
          rtt::analysis::opd_fan(cs, PathId{0}, field, 1, opd).warnings};
}

bool has(const std::vector<Diagnostic>& d, const std::string& code, const std::string& location) {
  for (const Diagnostic& x : d) {
    if (x.code == code && x.location == location) return true;
  }
  return false;
}

bool has_code(const std::vector<Diagnostic>& d, const std::string& code) {
  for (const Diagnostic& x : d) {
    if (x.code == code) return true;
  }
  return false;
}

}  // namespace

TEST_CASE("the reference singlet has no analysis warnings", "[analysis][warnings]") {
  for (const auto& w : all_warnings(compiled(singlet()), 1)) REQUIRE(w.empty());
}

TEST_CASE("rays.lost above the threshold, at the surface where most rays ended",
          "[analysis][warnings]") {
  // Lens aperture 3 mm with EPD 20 mm on axis: most rays end Vignetted at L1.S1.
  System s = singlet();
  lens_radius(s, 3.0);
  const CompiledSystem cs = compiled(s);
  for (const auto& w : all_warnings(cs, 0)) {
    REQUIRE(has(w, "rays.lost", "/root/children/1/surfaces/0"));
    REQUIRE_FALSE(has_code(w, "stop.clips_beam"));  // the stop does not clip here
    for (const Diagnostic& x : w) REQUIRE(x.severity == rtt::model::Severity::Warning);
  }
  // Vignetting is wanted at the field edge: a threshold above the lost fraction stays silent.
  for (const auto& w : all_warnings(cs, 0, 0.99)) REQUIRE_FALSE(has_code(w, "rays.lost"));
  // Losses are data in any case.
  REQUIRE(rtt::analysis::spot(cs, PathId{0}, 0, std::nullopt).losses.worst_surface ==
          cs.find_surface(rtt::model::SurfaceId("L1.S1")));
}

TEST_CASE("threshold 0 warns about any loss, threshold 1 never", "[analysis][warnings]") {
  // Lens aperture 9 mm: a few marginal rays of EPD 20 mm are lost, far less than half.
  System s = singlet();
  lens_radius(s, 9.0);
  const CompiledSystem few = compiled(s);
  REQUIRE(rtt::analysis::spot(few, PathId{0}, 0, std::nullopt).losses.worst_surface_count > 0);
  for (const auto& w : all_warnings(few, 0, 0.0)) REQUIRE(has_code(w, "rays.lost"));
  for (const auto& w : all_warnings(few, 0)) REQUIRE_FALSE(has_code(w, "rays.lost"));
  // Lens aperture 3 mm: most rays are lost, yet "more than all" never holds.
  lens_radius(s, 3.0);
  for (const auto& w : all_warnings(compiled(s), 0, 1.0)) REQUIRE_FALSE(has_code(w, "rays.lost"));
}

TEST_CASE("the threshold must be a fraction", "[analysis][warnings]") {
  const CompiledSystem cs = compiled(singlet());
  for (const double bad : {-0.1, 1.5, std::numeric_limits<double>::quiet_NaN()}) {
    rtt::analysis::SpotOptions options;
    options.lost_warning_fraction = bad;
    REQUIRE_THROWS_AS(rtt::analysis::spot(cs, PathId{0}, 0, std::nullopt, options),
                      std::invalid_argument);
  }
}

TEST_CASE("stop.clips_beam when the entrance pupil is larger than the stop opening",
          "[analysis][warnings]") {
  // The stop is the first surface, so the entrance pupil is the stop: EPD 30 mm against a stop
  // radius of 10 mm leaves the rays beyond 10 mm Vignetted at STO.
  System s = singlet();
  s.aperture.value = rtt::model::Param(30.0);
  const CompiledSystem cs = compiled(s);
  for (const auto& w : all_warnings(cs, 0)) {
    REQUIRE(has(w, "stop.clips_beam", "/root/children/0/surfaces/0"));
  }
}

TEST_CASE("no stop.clips_beam for marginal rays aimed exactly at the stop rim",
          "[analysis][warnings]") {
  // Aperture type stop_size with real aiming: the rays with px^2 + py^2 = 1 (outer hexapolar
  // ring, fan ends, grid corners inside the unit circle) are aimed at the rim of the stop. They
  // must arrive and must not warn because of rounding in the aiming.
  System s = singlet();
  s.aperture = {rtt::model::SystemApertureType::StopSize, rtt::model::Param(0.0)};
  const CompiledSystem cs = compiled(s);
  for (const std::uint16_t field : {std::uint16_t{0}, std::uint16_t{1}, std::uint16_t{2}}) {
    for (const auto& w : all_warnings(cs, field)) REQUIRE_FALSE(has_code(w, "stop.clips_beam"));
    const auto spot = rtt::analysis::spot(cs, PathId{0}, field, std::nullopt);
    REQUIRE(spot.losses.count(rtt::trace::RayStatus::Vignetted) == 0);
  }
}

TEST_CASE("report.paraxial_unavailable for a path the prescription rejects",
          "[analysis][warnings]") {
  // The system report of a diffraction order (#177): no paraxial data, a warning at the path.
  const CompiledSystem cs = rtt::compile::compile(
      rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m4/grating_transmission.rtt.json"),
      rtt::material::MaterialLibrary{});
  const PathId order = *cs.find_path("order +1");
  const rtt::analysis::SystemReport r = rtt::analysis::system_report(cs, order, 0);
  REQUIRE(r.warnings.size() == 1);
  CHECK(r.warnings[0].code == "report.paraxial_unavailable");
  CHECK(r.warnings[0].location == "/paths/" + std::to_string(order.index));
  // A path with paraxial data has none.
  CHECK(rtt::analysis::system_report(cs, *cs.find_path("order 0"), 0).warnings.empty());
}

TEST_CASE("every analysis code has a case in this file", "[analysis][warnings]") {
  const std::set<std::string> covered = {"rays.lost", "stop.clips_beam",
                                         "report.paraxial_unavailable"};
  for (const rtt::diagnostics::CodeInfo& info : rtt::diagnostics::kCodes) {
    if (info.producer != "analysis") continue;
    INFO("no case for " << info.code);
    REQUIRE(covered.contains(std::string(info.code)));
  }
}
