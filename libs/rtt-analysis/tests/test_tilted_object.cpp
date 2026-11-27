// A tilted object plane (#358, ADR 0055): object.slope = [sx, sy], the plane z = z_obj + sx x +
// sy y through the axial object point. m1/tilted_object_ideal_lens: an ideal lens f = 50 mm, the
// object 100 mm in front tilted by 45 deg (sy = -1), stop 1 mm before the lens, detector at the
// paraxial image of the axial point. The ideal lens is a collinear map, exact for every ray
// (ADR 0031, point 2): an object point at the distance u in front of it is imaged at
// v = u f / (u - f) behind it with the height -h v / u, so the tilted object plane is imaged onto
// a tilted plane (Scheimpflug), here at 45 deg (m = -1 at the axis).

#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <string>

#include "rtt/analysis/field.hpp"
#include "rtt/analysis/nominal_reference.hpp"
#include "rtt/analysis/spot.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/model/system.hpp"
#include "rtt/paraxial/paraxial.hpp"
#include "rtt/paraxial/prescription.hpp"
#include "rtt/paraxial/seidel.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/sources.hpp"

using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::model::System;

namespace {

System tilted() {
  return rtt::io::load_system(std::string(RTT_REFERENCE_DIR) +
                              "/m1/tilted_object_ideal_lens.rtt.json");
}

CompiledSystem compiled(const System& s) {
  const rtt::material::MaterialLibrary lib;
  return rtt::compile::compile(s, lib);
}

/// The chief ray of `field` (pupil centre) as make_rays starts it: at the object point.
rtt::trace::RayBatch chief(const CompiledSystem& cs, std::uint16_t field) {
  return rtt::trace::make_rays(cs, PathId{0}, std::span<const std::uint16_t>(&field, 1), 0,
                               rtt::trace::SinglePupilPoint{}, rtt::trace::Aiming::Real);
}

rtt::analysis::SpotOptions real_aiming() {
  rtt::analysis::SpotOptions o;
  o.aiming = rtt::trace::Aiming::Real;
  return o;
}

}  // namespace

TEST_CASE("field points lie on the tilted object plane", "[analysis][object][358]") {
  // object heights keep x, y: z = -100 + 0 x - 1 y
  const System s = tilted();
  const CompiledSystem cs = compiled(s);
  for (std::uint16_t f = 0; f < 4; ++f) {
    INFO("field " << f);
    const auto r = chief(cs, f);
    REQUIRE(r.size() == 1);
    CHECK(r.pos_y()[0] == s.fields.points[f].y);
    CHECK(r.pos_z()[0] == -100.0 - s.fields.points[f].y);
  }
  // field angles: the chief ray through the EP centre (the stop at z = -1) meets the plane
  System a = s;
  a.fields.type = rtt::model::FieldType::AngleDeg;
  a.fields.points = {{}, {0.0, 2.0, 1.0}, {1.0, -3.0, 1.0}};
  const CompiledSystem ca = compiled(a);
  const double z_ep = *rtt::paraxial::first_order(ca, PathId{0}, 0).entrance_pupil->z;
  REQUIRE(std::abs(z_ep + 1.0) <= 1e-12);
  for (std::uint16_t f = 1; f < 3; ++f) {
    INFO("angle field " << f);
    const auto r = chief(ca, f);
    const double x = r.pos_x()[0];
    const double y = r.pos_y()[0];
    const double z = r.pos_z()[0];
    CHECK(std::abs(z - (-100.0 - y)) <= 1e-12);  // on the plane
    const double deg = 3.14159265358979323846 / 180.0;
    CHECK(std::abs(x / (z - z_ep) - std::tan(a.fields.points[f].x * deg)) <= 1e-14);
    CHECK(std::abs(y / (z - z_ep) - std::tan(a.fields.points[f].y * deg)) <= 1e-14);
    CHECK(std::abs(z + 100.0) > 1.0);  // content guard: off the plane z = -100
  }
}

TEST_CASE("an ideal lens images the tilted object sharply onto the Scheimpflug plane",
          "[analysis][object][358]") {
  // The image plane of the collinear map: through (0, 0, 100) with dz/dy = 1, the detector
  // turned by 45 deg about x. A tilted detector has no paraxial model (ADR 0044), so the rays
  // are aimed with the reference of the same system with the detector normal to the axis.
  const System nominal = tilted();
  System scheimpflug = nominal;
  std::get<rtt::model::Element>(scheimpflug.root.children[2].value).pose.rotation_deg[0] =
      rtt::model::Param(45.0);
  const CompiledSystem cn = compiled(nominal);
  const CompiledSystem cs = compiled(scheimpflug);
  rtt::analysis::SpotOptions o = real_aiming();
  o.reference = rtt::analysis::nominal_reference(cn, PathId{0});
  for (std::uint16_t f = 1; f < 4; ++f) {
    INFO("field " << f);
    const auto sharp = rtt::analysis::spot(cs, PathId{0}, f, 0, o);
    const auto flat = rtt::analysis::spot(cn, PathId{0}, f, 0, real_aiming());
    REQUIRE(sharp.rays_arrived == 127);
    CHECK(sharp.stats.rms_centroid <= 1e-9);
    CHECK(flat.stats.rms_centroid > 0.1);  // content guard: defocused on the normal detector
  }
}

TEST_CASE("distortion of a tilted object through an ideal lens is zero",
          "[analysis][object][358]") {
  // The ideal lens is exact for every ray, so the real chief ray from the object point equals
  // the paraxial one from the same point (field.cpp, paraxial_chief for a tilted plane).
  const CompiledSystem cs = compiled(tilted());
  for (std::uint16_t f = 1; f < 4; ++f) {
    INFO("field " << f);
    const auto d = rtt::analysis::distortion_at(cs, PathId{0}, cs.fields().points[f], 0);
    CHECK(std::abs(d.real_height - d.paraxial_height) <= 1e-12);
    CHECK(std::abs(d.real_height) > 1.0);  // content guard
  }
}

TEST_CASE("first order ignores the tilt, Seidel and prescription reject it",
          "[analysis][object][358]") {
  const System s = tilted();
  System plane = s;
  plane.object.slope = {};
  const auto a = rtt::paraxial::first_order(compiled(s), PathId{0}, 0);
  const auto b = rtt::paraxial::first_order(compiled(plane), PathId{0}, 0);
  CHECK(a.efl == b.efl);
  CHECK(a.image_z == b.image_z);
  CHECK(a.entrance_pupil->z == b.entrance_pupil->z);
  const CompiledSystem cs = compiled(s);
  for (const bool seidel : {true, false}) {
    INFO((seidel ? "seidel" : "prescription"));
    try {
      if (seidel) {
        (void)rtt::paraxial::seidel(cs, PathId{0}, 0);
      } else {
        (void)rtt::paraxial::prescription(cs, PathId{0}, 0);
      }
      FAIL("no ParaxialError");
    } catch (const rtt::paraxial::ParaxialError& e) {
      CHECK(std::string(e.what()).find("tilted") != std::string::npos);
      CHECK(e.location() == "/object/slope");
    }
  }
  CHECK_NOTHROW(rtt::paraxial::seidel(compiled(plane), PathId{0}, 0));
}
