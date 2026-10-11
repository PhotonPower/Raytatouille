// The ideal lens in the analyses (#178, ADR 0031, tests of point 9 for PR 3): make_rays aims
// through it (rtt-paraxial knows its power), the spot of an object point of the design plane is
// a point and its OPD map is 0; without the design conjugate the OPD shows the r^4 rest.
//
// Setup of tests/reference/r2/ideal_lens.rtt.json: f = 50 mm, object at 75 mm, image at 150 mm,
// fields (0, 0), (0, 5), (3, -4) mm (object heights). Tolerances as in test_ideal_lens_trace.cpp:
// 1e-9 mm and 1e-6 waves, both more than 100 times the rounding.

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstdint>
#include <optional>
#include <string>
#include <variant>

#include "rtt/analysis/opd.hpp"
#include "rtt/analysis/spot.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"

using rtt::model::Element;
using rtt::model::IdealLens;
using rtt::model::System;

namespace {

System reference() {
  return rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/r2/ideal_lens.rtt.json");
}

rtt::compile::CompiledSystem compiled(const System& s) {
  return rtt::compile::compile(s, rtt::material::MaterialLibrary{});
}

}  // namespace

TEST_CASE("ideal lens: the spot of every field is a point", "[ideal][analysis]") {
  const rtt::compile::CompiledSystem cs = compiled(reference());
  const rtt::compile::PathId path = *cs.find_path("ideal lens");
  // Object heights of the fields; the image lies at m (x0, y0) with m = s'/s = -2.
  const std::array<std::array<double, 2>, 3> object{{{0.0, 0.0}, {0.0, 5.0}, {3.0, -4.0}}};
  for (const std::uint16_t field : {std::uint16_t{0}, std::uint16_t{1}, std::uint16_t{2}}) {
    INFO("field " << field);
    const rtt::analysis::SpotDiagram spot =
        rtt::analysis::spot(cs, path, field, cs.reference_wavelength());
    CHECK(spot.rays_arrived > 100);  // content guard: hexapolar 6 rings
    CHECK(spot.stats.rms_centroid < 1e-9);
    CHECK(spot.stats.geo_chief < 1e-9);
    CHECK(std::abs(spot.stats.centroid.x + 2.0 * object[field][0]) < 1e-9);
    CHECK(std::abs(spot.stats.centroid.y + 2.0 * object[field][1]) < 1e-9);
  }
}

TEST_CASE("ideal lens: OPD map 0 for the design conjugate", "[ideal][analysis][opd]") {
  const rtt::compile::CompiledSystem cs = compiled(reference());
  const rtt::compile::PathId path = *cs.find_path("ideal lens");
  for (const std::uint16_t field : {std::uint16_t{0}, std::uint16_t{1}, std::uint16_t{2}}) {
    INFO("field " << field);
    const rtt::analysis::OpdMap map =
        rtt::analysis::opd_map(cs, path, field, cs.reference_wavelength());
    CHECK(map.arrived > 100);
    CHECK(map.rms < 1e-6);  // waves
    CHECK(map.pv < 1e-6);
  }
}

TEST_CASE("ideal lens: without the design conjugate the OPD shows the r^4 rest",
          "[ideal][analysis][opd]") {
  // Content guard (ADR 0031, point 3): design at infinity, object at 75 mm: about 11 waves PV.
  System s = reference();
  std::get<IdealLens>(std::get<Element>(s.root.children[1].value).surfaces[0].interaction)
      .object_distance.reset();
  const rtt::compile::CompiledSystem cs = compiled(s);
  const rtt::analysis::OpdMap map =
      rtt::analysis::opd_map(cs, *cs.find_path("ideal lens"), 0, cs.reference_wavelength());
  CHECK(map.pv > 1.0);
  // The spot stays a point: the directions do not depend on the design conjugate.
  CHECK(rtt::analysis::spot(cs, *cs.find_path("ideal lens"), 0, cs.reference_wavelength())
            .stats.rms_centroid < 1e-9);
}
