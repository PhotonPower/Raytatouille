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
#include <numbers>
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
    const double deg = std::numbers::pi / 180.0;
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
        static_cast<void>(rtt::paraxial::seidel(cs, PathId{0}, 0));
      } else {
        static_cast<void>(rtt::paraxial::prescription(cs, PathId{0}, 0));
      }
      FAIL("no ParaxialError");
    } catch (const rtt::paraxial::ParaxialError& e) {
      CHECK(std::string(e.what()).find("tilted") != std::string::npos);
      CHECK(e.location() == "/object/slope");
    }
  }
  CHECK_NOTHROW(rtt::paraxial::seidel(compiled(plane), PathId{0}, 0));
}

TEST_CASE("an object point on or beyond the entrance pupil plane is an error",
          "[analysis][object][358]") {
  // review of #358: the direction of the rays and the pupil labels follow the side of z_obj
  System s = tilted();
  s.fields.points = {{}, {0.0, -99.0, 1.0}, {0.0, -120.0, 1.0}};  // P.z = -1 = z_EP, beyond
  const CompiledSystem cs = compiled(s);
  for (std::uint16_t f = 1; f < 3; ++f) {
    INFO("field " << f);
    CHECK_THROWS_AS(chief(cs, f), std::invalid_argument);
    CHECK_THROWS_AS(rtt::analysis::distortion_at(cs, PathId{0}, cs.fields().points[f], 0),
                    std::invalid_argument);
  }
  // a field angle whose chief ray meets the plane behind the EP: k = sy tan(theta) > 1
  System a = s;
  a.fields.type = rtt::model::FieldType::AngleDeg;
  a.fields.points = {{}, {0.0, -50.0, 1.0}};
  CHECK_THROWS_AS(chief(compiled(a), 1), std::invalid_argument);
}

TEST_CASE("a field angle on an anamorphic path meets the tilted plane per section",
          "[analysis][object][358]") {
  // examples/anamorphic/crossed_cylinders with the stop between the lenses (z = 20), so the
  // entrance pupils of the sections differ, and a tilted object 200 mm in front: the chief ray
  // x = (z - z_EP,x) tan(theta_x), y = (z - z_EP,y) tan(theta_y) meets z = z_obj + sx x + sy y
  // (ADR 0055, point 3).
  const CompiledSystem cs = compiled(rtt::io::parse_system(
      R"({"schema_version": "0.12.0", "name": "Crossed plano-convex cylinder lenses, n = 1.5168, t = 4 mm: L1 R = 51.68 mm (f_y = 100 mm) at z = 5, L2 R = 38.76 mm (f = 75 mm) at z = 30 rotated by 90 degree about z (power in x); x focus at z = 106.363, y focus at z = 107.726 (L2 acts as a plate in y); VACUUM", "units": {"length": "mm", "wavelength": "um"}, "environment": {"medium": "VACUUM"}, "wavelengths": [{"um": 0.5876, "reference": true}], "aperture": {"type": "epd", "value": 10.0}, "fields": {"type": "angle_deg", "points": [{}, {"x": 2.0, "y": 3.0}]}, "root": {"type": "assembly", "name": "system", "children": [{"type": "lens", "name": "L1", "pose": {"position": [0.0, 0.0, 5.0]}, "material": "CONST:1.5168", "surfaces": [{"id": "L1.S1", "shape": {"base": {"type": "cylinder", "radius": 51.68}}, "aperture": {"type": "rectangular", "half_width_x": 12.7, "half_width_y": 12.7}}, {"id": "L1.S2", "pose": {"position": [0.0, 0.0, 4.0]}, "aperture": {"type": "rectangular", "half_width_x": 12.7, "half_width_y": 12.7}}]}, {"type": "stop", "name": "stop", "surfaces": [{"id": "STO", "aperture": {"type": "circular", "radius": 10.0}}], "pose": {"position": [0.0, 0.0, 20.0]}}, {"type": "lens", "name": "L2", "pose": {"position": [0.0, 0.0, 30.0], "rotation_deg": [0.0, 0.0, 90.0]}, "material": "CONST:1.5168", "surfaces": [{"id": "L2.S1", "shape": {"base": {"type": "cylinder", "radius": 38.76}}, "aperture": {"type": "rectangular", "half_width_x": 12.7, "half_width_y": 12.7}}, {"id": "L2.S2", "pose": {"position": [0.0, 0.0, 4.0]}, "aperture": {"type": "rectangular", "half_width_x": 12.7, "half_width_y": 12.7}}]}, {"type": "detector", "name": "image", "pose": {"position": [0.0, 0.0, 106.36286919831225]}, "surfaces": [{"id": "IMG"}]}]}, "paths": [{"name": "main", "events": "auto"}], "object": {"at_infinity": false, "distance": 200.0, "slope": [0.3, -0.5]}})"));
  const auto zx = rtt::paraxial::first_order_section(cs, PathId{0}, cs.reference_wavelength(),
                                                     rtt::paraxial::Section::x)
                      .entrance_pupil->z;
  const auto zy = rtt::paraxial::first_order_section(cs, PathId{0}, cs.reference_wavelength(),
                                                     rtt::paraxial::Section::y)
                      .entrance_pupil->z;
  REQUIRE(zx);
  REQUIRE(zy);
  CHECK(std::abs(*zx - *zy) > 1.0);  // content guard: two pupils
  const std::uint16_t field = 1;
  const auto r = rtt::trace::make_rays(cs, PathId{0}, std::span<const std::uint16_t>(&field, 1),
                                       cs.reference_wavelength(), rtt::trace::SinglePupilPoint{},
                                       rtt::trace::Aiming::Paraxial);
  const double x = r.pos_x()[0];
  const double y = r.pos_y()[0];
  const double z = r.pos_z()[0];
  const double deg = std::numbers::pi / 180.0;
  CHECK(std::abs(z - (-200.0 + 0.3 * x - 0.5 * y)) <= 1e-12);
  CHECK(std::abs(x - (z - *zx) * std::tan(2.0 * deg)) <= 1e-12);
  CHECK(std::abs(y - (z - *zy) * std::tan(3.0 * deg)) <= 1e-12);
}

TEST_CASE("distortion with field angles and the place of the Scheimpflug image",
          "[analysis][object][358]") {
  // field angles on the tilted plane through the ideal lens: no distortion (paraxial_chief,
  // branch AngleDeg)
  System a = tilted();
  a.fields.type = rtt::model::FieldType::AngleDeg;
  a.fields.points = {{}, {0.0, 2.0, 1.0}, {1.5, -2.0, 1.0}};
  const CompiledSystem ca = compiled(a);
  for (std::uint16_t f = 1; f < 3; ++f) {
    const auto d = rtt::analysis::distortion_at(ca, PathId{0}, ca.fields().points[f], 0);
    CHECK(std::abs(d.real_height - d.paraxial_height) <= 1e-12);
    CHECK(std::abs(d.real_height) > 1.0);
  }
  // the chief ray of the object height h on the 45 deg detector: the collinear image point
  // (0, y', 100 + y') with y' = -50 h / (50 + h), at the local y = sqrt(2) |y'| from the vertex
  const System nominal = tilted();
  System scheimpflug = nominal;
  std::get<rtt::model::Element>(scheimpflug.root.children[2].value).pose.rotation_deg[0] =
      rtt::model::Param(45.0);
  rtt::analysis::SpotOptions o = real_aiming();
  o.reference = rtt::analysis::nominal_reference(compiled(nominal), PathId{0});
  const CompiledSystem cs = compiled(scheimpflug);
  for (std::uint16_t f = 1; f < 4; ++f) {
    const double h = nominal.fields.points[f].y;
    const double y_img = -50.0 * h / (50.0 + h);
    const auto s = rtt::analysis::spot(cs, PathId{0}, f, 0, o);
    CHECK(std::abs(s.chief.x) <= 1e-9);
    CHECK(std::abs(std::abs(s.chief.y) - std::sqrt(2.0) * std::abs(y_img)) <= 1e-9);
  }
}

TEST_CASE("the paraxial working F-number ignores the tilt", "[analysis][object][358]") {
  // decision on #358 after the review: it needs only the axial marginal ray
  System plane = tilted();
  plane.object.slope = {};
  const auto a = rtt::paraxial::paraxial_working_f_number(compiled(tilted()), PathId{0}, 0);
  const auto b = rtt::paraxial::prescription(compiled(plane), PathId{0}, 0);
  REQUIRE(a);
  CHECK(a == b.paraxial_working_f_number);
  CHECK(*a > 1.0);
}
