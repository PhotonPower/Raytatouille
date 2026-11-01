// object_na read as real NA (#296, ADR 0040, points 2, 3, 5 and 6): with reading real the marginal
// ray from the axial object point in the meridional plane has n_ref sin U = NA; r_EP is the radius
// of that cone in the EP plane, and the stop target R_s = sigma S comes from the stop hit S of
// that ray. Paraxial quantities stay at u = NA / n.
//
// Test system with analytic values: the object lies in immersion oil (n = 1.52) at z = -0.1 mm,
// an index-matched lens (plane front at z = 0) ends in a sphere of radius 1 mm around the object
// point (vertex at z = 0.9), a plane stop sits in vacuum at z = 2. Every ray from the object
// point leaves the sphere undeviated (normal incidence), so a ray of angle U meets the stop plane
// at y = 2.1 tan U. Paraxially the sphere is a lens of power (n - 1)/1 mm with the object at its
// centre: n u -> u, so the vacuum slope equals the slope u in the oil, and the paraxial stop
// height of the marginal ray is 2.1 u.

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <variant>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
#include "rtt/paraxial/paraxial.hpp"
#include "rtt/trace/sequential.hpp"
#include "rtt/trace/sources.hpp"

using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::model::ApertureReading;
using rtt::model::Element;
using rtt::model::ElementKind;
using rtt::model::FieldType;
using rtt::model::Param;
using rtt::model::Pose;
using rtt::model::Surface;
using rtt::model::SurfaceId;
using rtt::model::System;
using rtt::trace::Aiming;

namespace {

constexpr double kN = 1.52;        // immersion oil and lens
constexpr double kNa = 1.45;       // object-space NA
constexpr double kStopDist = 2.1;  // object to stop plane, mm

Surface surface(const std::string& id, double z, double radius = 0.0) {
  Surface s;
  s.id = SurfaceId(id);
  s.pose = Pose::along_z(z);
  if (radius != 0.0) s.shape.base = rtt::model::Conic{Param(radius), Param(0.0)};
  return s;
}

/// The immersion system of the file comment; `sphere` = false ends the lens in a plane instead,
/// where the steep marginal ray goes into total internal reflection (n sin U = NA > 1).
System immersion(double na, ApertureReading reading, bool sphere = true) {
  System s;
  s.name = "immersion objective front";
  s.environment.medium = "VACUUM";
  s.wavelengths = {{0.5876, 1.0, true}};
  s.object.at_infinity = false;
  s.object.distance = Param(0.1);
  s.object.medium = "CONST:1.52";
  s.aperture = {rtt::model::SystemApertureType::ObjectSpaceNA, Param(na)};
  s.aperture.reading = reading;
  s.fields = {FieldType::ObjectHeight, {{0.0, 0.0, 1.0}}};
  s.root.name = "root";
  Surface back = sphere ? surface("L.S2", 0.9, -1.0) : surface("L.S2", 0.9);
  back.aperture = rtt::model::CircularAperture{0.99, 0.0};
  Surface stop = surface("STO", 0.0);
  stop.aperture = rtt::model::CircularAperture{10.0, 0.0};
  s.root.children = {
      {Element{
          "L", ElementKind::Lens, Pose::along_z(0.0), "CONST:1.52", {surface("L.S1", 0.0), back}}},
      {Element{"stop", ElementKind::Stop, Pose::along_z(2.0), std::nullopt, {stop}}},
      {Element{"image",
               ElementKind::Detector,
               Pose::along_z(3.0),
               std::nullopt,
               {surface("IMG", 0.0)}}}};
  s.paths = {{"main", true, {}}};
  return s;
}

CompiledSystem compiled(const System& s) {
  return rtt::compile::compile(s, rtt::material::MaterialLibrary{});
}

/// The immersion system with N-BK7 as object medium and lens (index-matched at every wavelength)
/// and two wavelengths, the reference d line first.
CompiledSystem immersion_bk7() {
  System s = immersion(kNa, ApertureReading::Real);
  s.wavelengths = {{0.5876, 1.0, true}, {0.4861, 1.0, false}};
  s.object.medium = "SCHOTT:N-BK7";
  std::get<Element>(s.root.children[0].value).material = "SCHOTT:N-BK7";
  rtt::material::MaterialLibrary lib;
  lib.add_catalog(std::string(RTT_CATALOG_DIR) + "/schott.agf");
  return rtt::compile::compile(s, lib);
}

/// n sin U of a ray direction in the object medium (U against the axis, meridional).
double n_sin_u(const rtt::math::Vec3& dir) {
  return kN * std::hypot(dir.x(), dir.y()) / dir.norm();
}

}  // namespace

TEST_CASE("real NA: the marginal ray has n sin U = NA", "[trace][real-na]") {
  // ADR 0040, point 3: at the reference wavelength the real-aimed ray at p = (0, 1) is bitwise the
  // start ray through sigma r_EP, so Newton takes no step, and it has n_ref sin U = NA up to
  // rounding (1e-12). Content guard: it is alive (AimedRay status), and the stop target is
  // 2.1 tan U, far beyond the paraxial 2.1 NA / n.
  const CompiledSystem cs = compiled(immersion(kNa, ApertureReading::Real));
  for (const Aiming aiming : {Aiming::Real, Aiming::Paraxial}) {
    INFO("aiming " << static_cast<int>(aiming));
    const rtt::trace::AimedRay ray = rtt::trace::aim_ray(cs, PathId{0}, 0, 0, 0.0, 1.0, aiming);
    REQUIRE(ray.ray.status == rtt::trace::RayStatus::Alive);
    CHECK(ray.iterations == 0);
    CHECK(std::abs(n_sin_u(ray.ray.dir) - kNa) <= 1e-12);
    CHECK(ray.ray.dir.y() > 0.0);  // rising in object space
  }
}

TEST_CASE("real NA: the marginal ray is alive at the image", "[trace][real-na]") {
  // Content guard of ADR 0040 (acceptance): the real-aimed ray at p = (0, 1) passes the sphere
  // (height 0.954 mm < 0.99) and the stop (2.1 tan U = 6.63 mm < 10) and arrives at the image.
  const CompiledSystem cs = compiled(immersion(kNa, ApertureReading::Real));
  const std::array<std::uint16_t, 1> fields{0};
  rtt::trace::RayBatch rays =
      rtt::trace::make_rays(cs, PathId{0}, fields, 0, rtt::trace::SinglePupilPoint{0.0, 1.0});
  REQUIRE(rays.status()[0] == rtt::trace::RayStatus::Alive);
  [[maybe_unused]] const auto stats = rtt::trace::SequentialTracer().trace(cs, PathId{0}, rays);
  REQUIRE(rays.status()[0] == rtt::trace::RayStatus::Alive);
  // At the image (z = 3) the ray left the stop plane (z = 2) at 2.1 tan U with slope tan U
  // (undeviated at the sphere, file comment): y = 3.1 tan U, within the aim residual 1e-9 mm
  // (times 3.1 / 2.1).
  const double tan_u = std::tan(std::asin(kNa / kN));
  CHECK(std::abs(rays.pos_y()[0] - 3.1 * tan_u) <= 2e-9);
}

TEST_CASE("real NA: the same U and S at every wavelength", "[trace][real-na]") {
  // ADR 0040, points 2 and 3: U is formed with n_ref only, and R_s = sigma S is the same for all
  // wavelengths. N-BK7 in object space and lens keeps the system index-matched and the sphere
  // normal to every ray from the object point at every wavelength, so a ray meeting the stop at
  // S has the angle U of the reference wavelength: n(lambda) sin U = NA n(lambda) / n_ref.
  const CompiledSystem cs = immersion_bk7();
  const double n_ref = rtt::paraxial::first_order(cs, PathId{0}, 0).object_index;
  const double n_f = rtt::paraxial::first_order(cs, PathId{0}, 1).object_index;
  REQUIRE(n_f - n_ref > 5e-3);  // dispersion of N-BK7, d to F line: 0.0080
  const double tan_u = std::tan(std::asin(kNa / n_ref));
  // The stop target is bitwise the same (sigma = +1 at both wavelengths), the paraxial one not.
  CHECK(rtt::trace::stop_radius(cs, PathId{0}, 1) == rtt::trace::stop_radius(cs, PathId{0}, 0));
  CHECK(std::abs(rtt::trace::stop_radius(cs, PathId{0}, 0) - kStopDist * tan_u) <=
        1e-12 * kStopDist * tan_u);
  CHECK(rtt::trace::paraxial_stop_radius(cs, PathId{0}, 1) <
        rtt::trace::paraxial_stop_radius(cs, PathId{0}, 0));
  // Real aiming at the F line: the residual 1e-9 mm on the stop moves sin U by at most
  // cos^3 U 1e-9 / 2.1 < 1e-9.
  const rtt::trace::AimedRay ray = rtt::trace::aim_ray(cs, PathId{0}, 0, 1, 0.0, 1.0);
  REQUIRE(ray.ray.status == rtt::trace::RayStatus::Alive);
  const double sin_u = std::hypot(ray.ray.dir.x(), ray.ray.dir.y()) / ray.ray.dir.norm();
  CHECK(std::abs(sin_u - kNa / n_ref) <= 1e-9);
  CHECK(std::abs(n_f * sin_u - kNa) > 5e-3);  // n(lambda) sin U != NA
}

TEST_CASE("real NA: the paraxial reading gives a smaller real NA", "[trace][real-na]") {
  // reading paraxial, real aiming: the ray meets the stop at the paraxial R_s = 2.1 NA / n
  // (file comment), so tan U = NA / n and n sin U = NA / sqrt(1 + (NA / n)^2) (derived, not
  // given). The aim residual kAimTolerance on the stop moves tan U by 1e-9 / 2.1 and n sin U by
  // at most n cos^3 U 1e-9 / 2.1 < 1e-9.
  const CompiledSystem cs = compiled(immersion(kNa, ApertureReading::Paraxial));
  const rtt::trace::AimedRay ray = rtt::trace::aim_ray(cs, PathId{0}, 0, 0, 0.0, 1.0);
  REQUIRE(ray.ray.status == rtt::trace::RayStatus::Alive);
  const double t = kNa / kN;
  const double expected = kNa / std::sqrt(1.0 + t * t);
  CHECK(std::abs(n_sin_u(ray.ray.dir) - expected) <= 1e-9);
  CHECK(std::abs(n_sin_u(ray.ray.dir) - kNa) > 0.3);  // about 1.05 against 1.45
}

TEST_CASE("real NA: stop_radius and paraxial_stop_radius", "[trace][real-na]") {
  // ADR 0040, point 5: stop_radius() is |R_s| of the aim target, real with reading real;
  // paraxial_stop_radius() the paraxial value |(NA/n) |z_EP - z_obj| stop_scale|, bitwise equal to
  // stop_radius() with reading paraxial. Analytic: 2.1 tan U and 2.1 NA / n (file comment).
  const double tan_u = std::tan(std::asin(kNa / kN));
  const CompiledSystem real = compiled(immersion(kNa, ApertureReading::Real));
  const CompiledSystem par = compiled(immersion(kNa, ApertureReading::Paraxial));
  const double r_real = rtt::trace::stop_radius(real, PathId{0}, 0);
  const double r_par = rtt::trace::paraxial_stop_radius(real, PathId{0}, 0);
  CHECK(std::abs(r_real - kStopDist * tan_u) <= 1e-12 * kStopDist * tan_u);
  CHECK(std::abs(r_par - kStopDist * kNa / kN) <= 1e-12 * kStopDist);
  CHECK(rtt::trace::stop_radius(par, PathId{0}, 0) ==
        rtt::trace::paraxial_stop_radius(par, PathId{0}, 0));
  CHECK(rtt::trace::paraxial_stop_radius(par, PathId{0}, 0) == r_par);
  // The same through an AimReference.
  const rtt::trace::AimReference ref = rtt::trace::aim_reference(real, PathId{0}, 0);
  CHECK(rtt::trace::stop_radius(ref) == r_real);
  CHECK(rtt::trace::paraxial_stop_radius(ref) == r_par);
}

TEST_CASE("real NA: a virtual entrance pupil beyond the object", "[trace][real-na]") {
  // The singlet of #93 with the stop behind F' (z_EP < z_obj, sigma = -1), in vacuum, NA 0.02 read
  // real: the ray at p = (0, 1) rises in object space, Newton takes no step at lambda_ref, and
  // sin U = NA up to rounding.
  System s;
  s.name = "stop behind F'";
  s.environment.medium = "VACUUM";
  s.wavelengths = {{0.5876, 1.0, true}};
  s.object.at_infinity = false;
  s.object.distance = Param(200.0);
  s.aperture = {rtt::model::SystemApertureType::ObjectSpaceNA, Param(0.02)};
  s.aperture.reading = ApertureReading::Real;
  s.fields = {FieldType::ObjectHeight, {{0.0, 0.0, 1.0}}};
  s.root.name = "root";
  Surface s1 = surface("L1.S1", 0.0, 51.68);
  Surface s2 = surface("L1.S2", 4.0);
  s1.aperture = rtt::model::CircularAperture{25.0, 0.0};
  s2.aperture = rtt::model::CircularAperture{25.0, 0.0};
  Surface stop = surface("STO", 0.0);
  stop.aperture = rtt::model::CircularAperture{10.0, 0.0};
  s.root.children = {
      {Element{"L1", ElementKind::Lens, Pose::along_z(0.0), "CONST:1.5168", {s1, s2}}},
      {Element{"stop", ElementKind::Stop, Pose::along_z(150.0), std::nullopt, {stop}}},
      {Element{"image",
               ElementKind::Detector,
               Pose::along_z(205.0),
               std::nullopt,
               {surface("IMG", 0.0)}}}};
  s.paths = {{"main", true, {}}};
  const CompiledSystem cs = compiled(s);
  const rtt::trace::AimedRay ray = rtt::trace::aim_ray(cs, PathId{0}, 0, 0, 0.0, 1.0);
  REQUIRE(ray.ray.status == rtt::trace::RayStatus::Alive);
  CHECK(ray.iterations == 0);
  CHECK(ray.ray.dir.y() > 0.0);
  CHECK(std::abs(ray.ray.dir.y() / ray.ray.dir.norm() - 0.02) <= 1e-12);
}

TEST_CASE("real NA: aperture.na_unreachable", "[trace][real-na]") {
  // ADR 0040, point 6. The marginal ray in total internal reflection at a plane exit before the
  // stop: no S, so real aiming and stop_radius throw ApertureError, paraxial aiming does not.
  const CompiledSystem tir = compiled(immersion(kNa, ApertureReading::Real, false));
  try {
    (void)rtt::trace::aim_ray(tir, PathId{0}, 0, 0, 0.0, 1.0, Aiming::Real);
    FAIL("no ApertureError");
  } catch (const rtt::trace::ApertureError& e) {
    CHECK(e.code() == "aperture.na_unreachable");
    CHECK(e.location() == "/aperture/value");
  }
  CHECK_THROWS_AS(rtt::trace::stop_radius(tir, PathId{0}, 0), rtt::trace::ApertureError);
  const std::array<std::uint16_t, 1> fields{0};
  CHECK_THROWS_AS(rtt::trace::make_rays(tir, PathId{0}, fields, 0,
                                        rtt::trace::SinglePupilPoint{0.0, 1.0}, Aiming::Real),
                  rtt::trace::ApertureError);
  CHECK_NOTHROW(rtt::trace::aim_ray(tir, PathId{0}, 0, 0, 0.0, 1.0, Aiming::Paraxial));
  // NA >= n_ref: U is not defined, every aiming throws (ApertureError is an invalid_argument).
  const CompiledSystem over = compiled(immersion(1.6, ApertureReading::Real));
  for (const Aiming aiming : {Aiming::Real, Aiming::Paraxial}) {
    CHECK_THROWS_AS(rtt::trace::aim_ray(over, PathId{0}, 0, 0, 0.0, 0.0, aiming),
                    rtt::trace::ApertureError);
  }
  CHECK_THROWS_AS(rtt::trace::aim_ray(over, PathId{0}, 0, 0, 0.0, 0.0), std::invalid_argument);
  // Content guard: the same NA read paraxially is aimed (u = NA / n stays defined).
  CHECK_NOTHROW(rtt::trace::aim_ray(compiled(immersion(1.6, ApertureReading::Paraxial)), PathId{0},
                                    0, 0, 0.0, 0.0));
}

TEST_CASE("real NA: the reference compares the reading", "[trace][real-na][aim-reference]") {
  // ADR 0040, point 6, last item: aiming data with reading real do not fit a trial read paraxially.
  const CompiledSystem real = compiled(immersion(kNa, ApertureReading::Real));
  const CompiledSystem par = compiled(immersion(kNa, ApertureReading::Paraxial));
  const rtt::trace::AimReference ref = rtt::trace::aim_reference(real, PathId{0}, 0);
  CHECK_THROWS_WITH(rtt::trace::aim_ray(par, PathId{0}, 0, 0, 0.0, 1.0, Aiming::Real, ref),
                    Catch::Matchers::ContainsSubstring("system aperture"));
  // Content guard: with its own reference the ray is aimed, as without one.
  const rtt::trace::AimedRay a =
      rtt::trace::aim_ray(real, PathId{0}, 0, 0, 0.0, 1.0, Aiming::Real, ref);
  const rtt::trace::AimedRay b = rtt::trace::aim_ray(real, PathId{0}, 0, 0, 0.0, 1.0);
  CHECK(a.ray.dir == b.ray.dir);
}
