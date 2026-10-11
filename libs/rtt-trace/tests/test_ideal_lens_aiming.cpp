// Ray aiming through an ideal lens (#178, ADR 0031, point 6 and test of point 9 for PR 3): the
// lens 20 mm in front of the stop, a finite object, object height 5 mm. make_rays aims with the
// paraxial start values of first_order and traces with apply_event, so the aimed chief ray hits
// the centre of the stop and the rim rays (p_x, p_y) R_s, each within the aiming tolerance.

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstdint>
#include <string>
#include <utility>
#include <variant>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/sequential.hpp"
#include "rtt/trace/sources.hpp"

using rtt::model::Element;
using rtt::model::EventKind;
using rtt::model::Param;
using rtt::model::SurfaceId;
using rtt::model::System;

TEST_CASE("ray aiming through an ideal lens in front of the stop", "[ideal][aiming]") {
  // r2/ideal_lens with the stop (radius 5 mm, aperture stop_size) moved to z = 20, behind the
  // lens; a path that ends at the stop, so the last positions are the stop hits.
  System s = rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/r2/ideal_lens.rtt.json");
  Element& stop = std::get<Element>(s.root.children[0].value);
  stop.pose.position[2] = Param(20.0);
  stop.surfaces[0].aperture = rtt::model::CircularAperture{5.0, 0.0};
  s.paths.push_back(
      {"to stop",
       false,
       {{SurfaceId{"IL"}, EventKind::Transmit, 0}, {SurfaceId{"STO"}, EventKind::Transmit, 0}}});
  const rtt::compile::CompiledSystem cs =
      rtt::compile::compile(s, rtt::material::MaterialLibrary{});
  const rtt::compile::PathId path = *cs.find_path("to stop");
  const std::uint16_t field = 1;  // object height (0, 5) mm
  const std::array<std::uint16_t, 1> fields{field};
  const double r_stop = 5.0;
  for (const auto& [px, py] : {std::pair{0.0, 0.0}, std::pair{1.0, 0.0}, std::pair{0.0, 1.0},
                               std::pair{0.0, -1.0}, std::pair{-0.6, 0.8}}) {
    INFO("pupil (" << px << ", " << py << ")");
    rtt::trace::RayBatch rays =
        rtt::trace::make_rays(cs, path, fields, cs.reference_wavelength(),
                              rtt::trace::SinglePupilPoint{px, py}, rtt::trace::Aiming::Real);
    static_cast<void>(rtt::trace::SequentialTracer{}.trace(cs, path, rays));
    REQUIRE(rays.status()[0] == rtt::trace::RayStatus::Alive);
    // The stop lies at z = 20 without rotation: local x, y = global x, y.
    CHECK(std::abs(rays.pos_x()[0] - px * r_stop) < 1e-9);
    CHECK(std::abs(rays.pos_y()[0] - py * r_stop) < 1e-9);
    CHECK(std::abs(rays.pos_z()[0] - 20.0) < 1e-12);
  }
}
