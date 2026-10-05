#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numbers>
#include <optional>
#include <random>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
#include "rtt/paraxial/paraxial.hpp"
#include "rtt/trace/apply_event.hpp"
#include "rtt/trace/sequential.hpp"
#include "rtt/trace/sources.hpp"

using rtt::compile::compile;
using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::material::MaterialLibrary;
using rtt::math::Vec3;
using rtt::model::Element;
using rtt::model::ElementKind;
using rtt::model::FieldType;
using rtt::model::Param;
using rtt::model::Pose;
using rtt::model::Surface;
using rtt::model::SurfaceId;
using rtt::model::System;
using rtt::trace::Aiming;
using rtt::trace::RayState;
using rtt::trace::RayStatus;

// Reference cases of issue #8. Paraxial reference values (entrance pupil, stop radius, focal
// points, image heights) come from rtt-paraxial (#7), whose formulas are verified there.

namespace {

constexpr double kDeg = std::numbers::pi / 180.0;

Surface surface(const std::string& id, double z_mm = 0.0, std::optional<double> radius = {}) {
  Surface s;
  s.id = SurfaceId(id);
  s.pose = Pose::along_z(z_mm);
  if (radius) s.shape.base = rtt::model::Conic{Param(*radius), Param(0.0)};
  return s;
}

enum class StopPlace { Before, Between, After };

/// Two singlets (CONST:1.5168) with the stop before, between or after them; EPD 10 mm; fields
/// as angles (object at infinity) or object heights (object at z = -100 mm).
System two_lenses(StopPlace place, bool finite_object = false) {
  System s;
  s.name = "two lenses";
  s.wavelengths = {{0.5876, 1.0, true}};
  s.aperture = {rtt::model::SystemApertureType::EntrancePupilDiameter, Param(10.0)};
  if (finite_object) {
    s.object.at_infinity = false;
    s.object.distance = Param(100.0);
    s.fields = {FieldType::ObjectHeight, {{0.0, 0.0, 1.0}, {0.0, 8.0, 1.0}, {-5.0, 3.0, 1.0}}};
  } else {
    s.fields = {FieldType::AngleDeg, {{0.0, 0.0, 1.0}, {0.0, 5.0, 1.0}, {-3.0, 4.0, 1.0}}};
  }
  s.root.name = "root";
  Surface stop = surface("STO");
  stop.aperture = rtt::model::CircularAperture{50.0, 0.0};  // never vignettes in these tests
  const Element stop_element{"stop", ElementKind::Stop, Pose::along_z(0.0), std::nullopt, {stop}};
  Element l1{"L1",
             ElementKind::Lens,
             Pose::along_z(10.0),
             "CONST:1.5168",
             {surface("L1.S1", 0.0, 60.0), surface("L1.S2", 5.0, -60.0)}};
  Element l2{"L2",
             ElementKind::Lens,
             Pose::along_z(25.0),
             "CONST:1.5168",
             {surface("L2.S1", 0.0, 45.0), surface("L2.S2", 4.0)}};
  Element stop_at = stop_element;
  switch (place) {
    case StopPlace::Before:
      stop_at.pose = Pose::along_z(0.0);
      s.root.children = {{stop_at}, {l1}, {l2}};
      break;
    case StopPlace::Between:
      stop_at.pose = Pose::along_z(20.0);
      s.root.children = {{l1}, {stop_at}, {l2}};
      break;
    case StopPlace::After:
      stop_at.pose = Pose::along_z(35.0);
      s.root.children = {{l1}, {l2}, {stop_at}};
      break;
  }
  s.root.children.push_back(
      {Element{"D", ElementKind::Detector, Pose::along_z(80.0), std::nullopt, {surface("IMG")}}});
  s.paths = {{"main", true, {}}};
  return s;
}

System singlet() {
  return rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m1/singlet_const.rtt.json");
}

/// Index of the first Stop event of path 0.
std::size_t stop_event(const CompiledSystem& cs) {
  const auto& events = cs.path(PathId{0}).events;
  for (std::size_t i = 0; i < events.size(); ++i) {
    if (cs.surfaces()[events[i].surface].element_kind == ElementKind::Stop) return i;
  }
  FAIL("no stop");
  return 0;
}

/// Traces the ray to the stop surface (apertures of the surfaces before the stop included) and
/// returns the hit in local stop coordinates; the stop aperture itself is not applied, so that
/// rays aimed exactly at its rim can be checked.
Vec3 hit_on_stop(const CompiledSystem& cs, RayState ray) {
  const auto& events = cs.path(PathId{0}).events;
  const std::size_t stop = stop_event(cs);
  for (std::size_t i = 0; i < stop; ++i) {
    const auto& e = events[i];
    ray = rtt::trace::sequential_step(ray, cs.surfaces()[e.surface], e.surface, e.kind,
                                      cs.media()[e.medium_before].index[0].real(),
                                      cs.media()[e.medium_after].index[0].real());
  }
  REQUIRE(ray.status == RayStatus::Alive);
  const auto hit = rtt::trace::intersect_surface(ray, cs.surfaces()[events[stop].surface]);
  REQUIRE(hit.status == rtt::geom::HitStatus::Hit);
  return hit.point;
}

/// Signed paraxial stop radius that belongs to the entrance pupil: height of the paraxial
/// marginal ray at the stop (object at infinity: y = EPD / 2 parallel to the axis; finite
/// object: from the axial object point through the EP rim). Negative if the stop sees the EP
/// inverted; pupil point (px, py) then aims at (px, py) * R_s.
double paraxial_stop_radius(const CompiledSystem& cs) {
  const auto fo = rtt::paraxial::first_order(cs, PathId{0}, 0);
  REQUIRE(fo.entrance_pupil);
  REQUIRE(fo.entrance_pupil->diameter);
  const double h = *fo.entrance_pupil->diameter / 2.0;
  std::vector<rtt::paraxial::RayAtEvent> ray;
  if (cs.object().at_infinity) {
    ray = rtt::paraxial::trace_ray(cs, PathId{0}, 0, -50.0, h, 0.0);
  } else {
    const double z_obj = -cs.object().distance.value;
    REQUIRE(fo.entrance_pupil->z);
    ray =
        rtt::paraxial::trace_ray(cs, PathId{0}, 0, z_obj, 0.0, h / (*fo.entrance_pupil->z - z_obj));
  }
  return ray[stop_event(cs)].y;
}

}  // namespace

TEST_CASE("pupil samplings follow the documented order", "[sources]") {
  using rtt::trace::pupil_points;
  SECTION("hexapolar: centre plus 6k rays on ring k, starting at +y") {
    const auto p = pupil_points(rtt::trace::HexapolarPupil{2});
    REQUIRE(p.size() == 1 + 6 + 12);
    REQUIRE((p[0].px == 0.0 && p[0].py == 0.0));
    REQUIRE((std::abs(p[1].px) < 1e-15 && std::abs(p[1].py - 0.5) < 1e-15));
    // Second ray of ring 1 at 60 deg from +y towards +x.
    REQUIRE(std::abs(p[2].px - 0.5 * std::sin(60.0 * kDeg)) < 1e-15);
    REQUIRE(std::abs(p[2].py - 0.5 * std::cos(60.0 * kDeg)) < 1e-15);
    REQUIRE((std::abs(p[7].px) < 1e-15 && std::abs(p[7].py - 1.0) < 1e-15));
    REQUIRE(pupil_points(rtt::trace::HexapolarPupil{0}).size() == 1);
  }
  SECTION("grid: n x n on [-1, 1]^2, rim included, corners outside") {
    const auto p = pupil_points(rtt::trace::GridPupil{3});
    REQUIRE(p.size() == 5);
    REQUIRE((p[0].px == 0.0 && p[0].py == -1.0));
    REQUIRE((p[1].px == -1.0 && p[1].py == 0.0));
    REQUIRE((p[2].px == 0.0 && p[2].py == 0.0));
    REQUIRE((p[4].px == 0.0 && p[4].py == 1.0));
    REQUIRE(pupil_points(rtt::trace::GridPupil{1}).size() == 1);
    // Step 0.2: lattice points (i, j) with i^2 + j^2 <= 25, Gauss circle N(5) = 81.
    REQUIRE(pupil_points(rtt::trace::GridPupil{11}).size() == 81);
  }
  SECTION("fans") {
    const auto y = pupil_points(rtt::trace::FanYPupil{5});
    REQUIRE(y.size() == 5);
    for (std::size_t i = 0; i < 5; ++i) {
      REQUIRE(y[i].px == 0.0);
      REQUIRE(std::abs(y[i].py - (-1.0 + 0.5 * static_cast<double>(i))) < 1e-15);
    }
    const auto x = pupil_points(rtt::trace::FanXPupil{3});
    REQUIRE((x[0].px == -1.0 && x[2].px == 1.0 && x[1].py == 0.0));
    const auto one = pupil_points(rtt::trace::FanXPupil{1});
    REQUIRE((one.size() == 1 && one[0].px == 0.0 && one[0].py == 0.0));
  }
  SECTION("single point") {
    const auto p = pupil_points(rtt::trace::SinglePupilPoint{0.3, -0.4});
    REQUIRE((p.size() == 1 && p[0].px == 0.3 && p[0].py == -0.4));
  }
  SECTION("invalid sizes") {
    REQUIRE_THROWS_AS(pupil_points(rtt::trace::GridPupil{0}), std::invalid_argument);
    REQUIRE_THROWS_AS(pupil_points(rtt::trace::FanYPupil{0}), std::invalid_argument);
    REQUIRE_THROWS_AS(pupil_points(rtt::trace::HexapolarPupil{-1}), std::invalid_argument);
  }
}

TEST_CASE("random pupil sampling is reproducible on every platform", "[sources]") {
  // The C++ standard fixes the 10000th output of a default-constructed std::mt19937_64
  // ([rand.predef]); this anchors the engine itself.
  std::mt19937_64 engine;
  engine.discard(9999);
  REQUIRE(engine() == 9981545732273789042ULL);

  const auto a = rtt::trace::pupil_points(rtt::trace::RandomPupil{1000, 42});
  const auto b = rtt::trace::pupil_points(rtt::trace::RandomPupil{1000, 42});
  const auto c = rtt::trace::pupil_points(rtt::trace::RandomPupil{1000, 43});
  REQUIRE(a.size() == 1000);
  for (std::size_t i = 0; i < a.size(); ++i) {
    REQUIRE((a[i].px == b[i].px && a[i].py == b[i].py));
    REQUIRE(a[i].px * a[i].px + a[i].py * a[i].py <= 1.0);
  }
  REQUIRE((a[0].px != c[0].px || a[0].py != c[0].py));

  // Independent re-computation of the documented mapping for the first ray.
  std::mt19937_64 gen(42);
  const double u1 = static_cast<double>(gen() >> 11) * 0x1.0p-53;
  const double u2 = static_cast<double>(gen() >> 11) * 0x1.0p-53;
  const double r = std::sqrt(u1);
  const double phi = 2.0 * std::numbers::pi * u2;
  REQUIRE(a[0].px == r * std::sin(phi));
  REQUIRE(a[0].py == r * std::cos(phi));

  // Fixed values, recorded once with GCC 16 (MinGW). They must not change across platforms or
  // library versions. The draws and sqrt are exact; std::sin and std::cos are not correctly
  // rounded and may differ by a few ulp between C libraries, hence 1e-15 instead of ==.
  const auto fixed = [&](std::size_t i, double px, double py) {
    INFO("ray " << i);
    REQUIRE(std::abs(a[i].px - px) <= 1e-15);
    REQUIRE(std::abs(a[i].py - py) <= 1e-15);
  };
  fixed(0, -0x1.5516de0aed9dp-1, -0x1.1db018d033d6fp-1);
  fixed(1, 0x1.4f6a72e3e6742p-1, 0x1.22f9a869a799ep-1);
  fixed(999, -0x1.843e6827c46b4p-1, -0x1.29c18995e5e9dp-2);
}

TEST_CASE("chief ray hits the stop centre for stop before, inside and after the group",
          "[sources][aiming]") {
  // Issue #8: < 1e-9 mm.
  const MaterialLibrary lib;
  for (const StopPlace place : {StopPlace::Before, StopPlace::Between, StopPlace::After}) {
    for (const bool finite : {false, true}) {
      const CompiledSystem cs = compile(two_lenses(place, finite), lib);
      for (std::uint16_t f = 0; f < 3; ++f) {
        INFO("stop " << static_cast<int>(place) << ", finite " << finite << ", field " << f);
        const auto aimed = rtt::trace::aim_ray(cs, PathId{0}, f, 0, 0.0, 0.0);
        REQUIRE(aimed.ray.status == RayStatus::Alive);
        const Vec3 hit = hit_on_stop(cs, aimed.ray);
        REQUIRE(std::hypot(hit.x(), hit.y()) < 1e-9);
      }
    }
  }
}

TEST_CASE("pupil rays hit their target on the stop", "[sources][aiming]") {
  // Target (px R_s, py R_s) with the signed paraxial stop radius R_s that belongs to the EP.
  const MaterialLibrary lib;
  for (const StopPlace place : {StopPlace::Before, StopPlace::Between, StopPlace::After}) {
    for (const bool finite : {false, true}) {
      const CompiledSystem cs = compile(two_lenses(place, finite), lib);
      const double r_s = paraxial_stop_radius(cs);
      for (const auto& [px, py] :
           {std::pair{0.0, 1.0}, std::pair{1.0, 0.0}, std::pair{-0.7, 0.7}, std::pair{0.2, -0.5}}) {
        INFO("stop " << static_cast<int>(place) << ", finite " << finite << ", p = (" << px << ", "
                     << py << ")");
        const auto aimed = rtt::trace::aim_ray(cs, PathId{0}, 2, 0, px, py);
        REQUIRE(aimed.ray.status == RayStatus::Alive);
        const Vec3 hit = hit_on_stop(cs, aimed.ray);
        REQUIRE(std::hypot(hit.x() - px * r_s, hit.y() - py * r_s) < 1e-9);
      }
    }
  }
}

TEST_CASE("field conventions: angle, object height, direction of +y", "[sources]") {
  const MaterialLibrary lib;
  SECTION("angle: d = (tan theta_x, tan theta_y, 1) normalised") {
    const CompiledSystem cs = compile(two_lenses(StopPlace::Between), lib);
    const auto aimed = rtt::trace::aim_ray(cs, PathId{0}, 2, 0, 0.3, -0.2);
    const Vec3 expected = Vec3(std::tan(-3.0 * kDeg), std::tan(4.0 * kDeg), 1.0).normalized();
    REQUIRE((aimed.ray.dir - expected).cwiseAbs().maxCoeff() <= 1e-15);
  }
  SECTION("object height: the ray starts in (x, y, -distance)") {
    const CompiledSystem cs = compile(two_lenses(StopPlace::Between, true), lib);
    const auto aimed = rtt::trace::aim_ray(cs, PathId{0}, 2, 0, 0.5, 0.5);
    REQUIRE(aimed.ray.pos == Vec3(-5.0, 3.0, -100.0));
    REQUIRE(aimed.ray.opl == 0.0);
  }
  SECTION("object height needs a finite object") {
    System s = two_lenses(StopPlace::Between);
    s.fields.type = FieldType::ObjectHeight;
    const CompiledSystem cs = compile(s, lib);
    REQUIRE_THROWS_AS(rtt::trace::aim_ray(cs, PathId{0}, 1, 0, 0.0, 0.0), std::invalid_argument);
  }
}

TEST_CASE("paraxial image height is converted with the paraxial chief ray", "[sources]") {
  const MaterialLibrary lib;
  SECTION("object at infinity: tan(theta) = y' / EFL") {
    // Paraxial image height of a distant object: y' = EFL tan(theta), independent of the stop.
    // That relation holds for object and image space with n = 1, so the file's AIR (Ciddor
    // since #25) is replaced by VACUUM here.
    System s = singlet();
    s.environment.medium = "VACUUM";
    s.fields = {FieldType::ParaxialImageHeight, {{0.0, 0.0, 1.0}, {1.0, 2.0, 1.0}}};
    const CompiledSystem cs = compile(s, lib);
    const double efl = *rtt::paraxial::first_order(cs, PathId{0}, 1).efl;
    const auto aimed = rtt::trace::aim_ray(cs, PathId{0}, 1, 1, 0.0, 0.0);
    REQUIRE(std::abs(aimed.ray.dir.x() / aimed.ray.dir.z() - 1.0 / efl) <= 1e-12);
    REQUIRE(std::abs(aimed.ray.dir.y() / aimed.ray.dir.z() - 2.0 / efl) <= 1e-12);
  }
  SECTION("finite object: object height = y' / m") {
    System s = two_lenses(StopPlace::Between, true);
    s.fields = {FieldType::ParaxialImageHeight, {{0.0, 0.0, 1.0}, {0.0, -3.0, 1.0}}};
    const CompiledSystem cs = compile(s, lib);
    const auto fo = rtt::paraxial::first_order(cs, PathId{0}, 0);
    REQUIRE(fo.lateral_magnification);
    const auto aimed = rtt::trace::aim_ray(cs, PathId{0}, 1, 0, 0.0, 0.0);
    REQUIRE(std::abs(aimed.ray.pos.y() - (-3.0 / *fo.lateral_magnification)) <= 1e-12);
    REQUIRE(aimed.ray.pos.z() == -100.0);
  }
}

TEST_CASE("object at infinity: rays of a field start on one plane wave before the system",
          "[sources]") {
  // Large field: the tilted start plane must still lie before the first surface for the whole
  // bundle, and d . pos is the same for all rays of a field (plane wave, needed for OPD in M2).
  // Stop between the lenses: with the stop behind the group, no ray of a 30 deg bundle can reach
  // its stop target (reachable stop heights 8.6 ... 18.4 mm, targets within +-2.88 mm; lower
  // rays would hit L1 outside its rim, where its edge thickness drops to zero at r = 17.1 mm),
  // which the NoConvergence tests cover.
  const MaterialLibrary lib;
  System s = two_lenses(StopPlace::Between);
  s.fields = {FieldType::AngleDeg, {{0.0, 0.0, 1.0}, {0.0, 30.0, 1.0}, {20.0, -20.0, 1.0}}};
  const CompiledSystem cs = compile(s, lib);
  const double z_ep = *rtt::paraxial::first_order(cs, PathId{0}, 0).entrance_pupil->z;
  const std::vector<std::uint16_t> fields{0, 1, 2};
  rtt::trace::RayBatch rays =
      rtt::trace::make_rays(cs, PathId{0}, fields, 0, rtt::trace::HexapolarPupil{4});
  const std::size_t per_field = rays.size() / 3;
  REQUIRE(per_field == 61);
  for (std::size_t f = 0; f < 3; ++f) {
    const std::size_t first = f * per_field;
    const Vec3 d0(rays.dir_x()[first], rays.dir_y()[first], rays.dir_z()[first]);
    const double plane =
        d0.dot(Vec3(rays.pos_x()[first], rays.pos_y()[first], rays.pos_z()[first]));
    for (std::size_t i = first; i < first + per_field; ++i) {
      INFO("field " << f << ", pupil (" << rays.pupil_x()[i] << ", " << rays.pupil_y()[i] << ")");
      REQUIRE(rays.field()[i] == f);
      REQUIRE(rays.status()[i] == RayStatus::Alive);
      const Vec3 p(rays.pos_x()[i], rays.pos_y()[i], rays.pos_z()[i]);
      REQUIRE(std::abs(d0.dot(p) - plane) <= 1e-9);
      // At least 1 mm before the first vertex (z = 10) and the entrance pupil.
      REQUIRE(p.z() <= 9.0);
      REQUIRE(p.z() <= z_ep - 1.0);
    }
  }
  // All rays reach the detector.
  const auto stats = rtt::trace::SequentialTracer().trace(cs, PathId{0}, rays);
  REQUIRE(stats.count(RayStatus::Alive) == rays.size());
}

TEST_CASE("paraxial aiming passes through the paraxial entrance pupil", "[sources]") {
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(two_lenses(StopPlace::After), lib);
  const auto fo = rtt::paraxial::first_order(cs, PathId{0}, 0);
  const double z_ep = *fo.entrance_pupil->z;
  const double r_ep = *fo.entrance_pupil->diameter / 2.0;
  const auto aimed = rtt::trace::aim_ray(cs, PathId{0}, 2, 0, 0.6, -0.8, Aiming::Paraxial);
  REQUIRE(aimed.iterations == 0);
  const RayState& r = aimed.ray;
  const double t = (z_ep - r.pos.z()) / r.dir.z();
  REQUIRE(std::abs(r.pos.x() + t * r.dir.x() - 0.6 * r_ep) <= 1e-12);
  REQUIRE(std::abs(r.pos.y() + t * r.dir.y() + 0.8 * r_ep) <= 1e-12);
}

TEST_CASE("make_rays fills the batch per field and pupil point", "[sources]") {
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(two_lenses(StopPlace::Between), lib);
  const std::vector<std::uint16_t> fields{2, 0};
  const auto rays = rtt::trace::make_rays(cs, PathId{0}, fields, 0, rtt::trace::FanYPupil{3});
  REQUIRE(rays.size() == 6);
  const std::vector<double> py{-1.0, 0.0, 1.0};
  for (std::size_t i = 0; i < 6; ++i) {
    REQUIRE(rays.field()[i] == fields[i / 3]);
    REQUIRE(rays.pupil_x()[i] == 0.0);
    REQUIRE(rays.pupil_y()[i] == py[i % 3]);
    REQUIRE(rays.wl()[i] == 0);
    REQUIRE(rays.opl()[i] == 0.0);
    REQUIRE(rays.weight()[i] == 1.0);
    REQUIRE(rays.last_surface()[i] == rtt::trace::kNoSurface);
    const auto aimed = rtt::trace::aim_ray(cs, PathId{0}, fields[i / 3], 0, 0.0, py[i % 3]);
    REQUIRE(rays.pos_y()[i] == aimed.ray.pos.y());
    REQUIRE(rays.dir_z()[i] == aimed.ray.dir.z());
  }
}

TEST_CASE("small aperture: real focus converges to the paraxial focus with NA^2", "[sources]") {
  // Third-order longitudinal spherical aberration grows with the square of the aperture, so
  // the axis crossing of a real marginal ray at relative pupil height rho deviates from the
  // paraxial focus F' by ~ rho^2 (architecture: "Fehler proportional NA^2", fit exponent
  // 2 +- 0.05, issue #8). Reference focus: rear_focal_z from rtt-paraxial (#7).
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(singlet(), lib);
  const auto fo = rtt::paraxial::first_order(cs, PathId{0}, 1);
  const double z_f = *fo.rear_focal_z;
  const std::vector<double> rho{0.02, 0.04, 0.08, 0.16};
  std::vector<double> deviation;
  for (const double r : rho) {
    rtt::trace::RayBatch rays(1);
    const auto aimed = rtt::trace::aim_ray(cs, PathId{0}, 0, 1, 0.0, r);
    REQUIRE(aimed.ray.status == RayStatus::Alive);
    rays.pos_x()[0] = aimed.ray.pos.x();
    rays.pos_y()[0] = aimed.ray.pos.y();
    rays.pos_z()[0] = aimed.ray.pos.z();
    rays.dir_x()[0] = aimed.ray.dir.x();
    rays.dir_y()[0] = aimed.ray.dir.y();
    rays.dir_z()[0] = aimed.ray.dir.z();
    rays.wl()[0] = 1;
    [[maybe_unused]] const auto stats = rtt::trace::SequentialTracer().trace(cs, PathId{0}, rays);
    REQUIRE(rays.status()[0] == RayStatus::Alive);
    // Axis crossing of the ray after the lens (meridional ray, x = 0).
    const double z_cross = rays.pos_z()[0] - rays.pos_y()[0] * rays.dir_z()[0] / rays.dir_y()[0];
    deviation.push_back(std::abs(z_cross - z_f));
    // Issue #8: focus matches the paraxial BFL, relative 1e-3, at small aperture.
    if (r <= 0.08) {
      REQUIRE(std::abs(z_cross - z_f) <= 1e-3 * *fo.bfl);
    }
  }
  // Least-squares slope of log(deviation) over log(rho).
  double mx = 0.0, my = 0.0;
  for (std::size_t i = 0; i < rho.size(); ++i) {
    mx += std::log(rho[i]);
    my += std::log(deviation[i]);
  }
  mx /= static_cast<double>(rho.size());
  my /= static_cast<double>(rho.size());
  double sxy = 0.0, sxx = 0.0;
  for (std::size_t i = 0; i < rho.size(); ++i) {
    sxy += (std::log(rho[i]) - mx) * (std::log(deviation[i]) - my);
    sxx += (std::log(rho[i]) - mx) * (std::log(rho[i]) - mx);
  }
  const double exponent = sxy / sxx;
  INFO("fit exponent " << exponent);
  REQUIRE(std::abs(exponent - 2.0) <= 0.05);
}

TEST_CASE("aiming that cannot reach the stop gives NoConvergence, not an exception",
          "[sources][aiming]") {
  // Stop behind the group and a field so steep that the rays cannot reach the stop through the
  // lenses.
  const MaterialLibrary lib;
  System s = two_lenses(StopPlace::After);
  s.fields = {FieldType::AngleDeg, {{0.0, 0.0, 1.0}, {0.0, 80.0, 1.0}}};
  const CompiledSystem cs = compile(s, lib);
  const auto aimed = rtt::trace::aim_ray(cs, PathId{0}, 1, 0, 0.0, 1.0);
  REQUIRE(aimed.ray.status == RayStatus::NoConvergence);
  REQUIRE(aimed.ray.pos.allFinite());
}

TEST_CASE("invalid source input throws at the API boundary", "[sources]") {
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(two_lenses(StopPlace::Between), lib);
  REQUIRE_THROWS_AS(rtt::trace::aim_ray(cs, PathId{0}, 3, 0, 0.0, 0.0), std::invalid_argument);
  REQUIRE_THROWS_AS(rtt::trace::aim_ray(cs, PathId{0}, 0, 2, 0.0, 0.0), std::invalid_argument);
  REQUIRE_THROWS_AS(rtt::trace::aim_ray(cs, PathId{4}, 0, 0, 0.0, 0.0), std::invalid_argument);
  // A path without a stop.
  System no_stop = two_lenses(StopPlace::Between);
  no_stop.root.children.erase(no_stop.root.children.begin() + 1);
  const CompiledSystem cs2 = compile(no_stop, lib);
  REQUIRE_THROWS_AS(rtt::trace::aim_ray(cs2, PathId{0}, 0, 0, 0.0, 0.0), std::invalid_argument);
}

TEST_CASE("field angle with a finite object: object point on the chief ray through the EP",
          "[sources]") {
  // Decided for #8: a field angle with a finite object places the object point on the chief
  // ray d ~ (tan theta_x, tan theta_y, 1) through the EP centre, i.e. at
  // ((z_obj - z_ep) tan theta_x, (z_obj - z_ep) tan theta_y, z_obj).
  const MaterialLibrary lib;
  for (const StopPlace place : {StopPlace::Before, StopPlace::Between, StopPlace::After}) {
    System s = two_lenses(place, true);
    s.fields = {FieldType::AngleDeg, {{0.0, 0.0, 1.0}, {-3.0, 4.0, 1.0}}};
    const CompiledSystem cs = compile(s, lib);
    const double z_ep = *rtt::paraxial::first_order(cs, PathId{0}, 0).entrance_pupil->z;
    const auto aimed = rtt::trace::aim_ray(cs, PathId{0}, 1, 0, 0.0, 0.0);
    INFO("stop " << static_cast<int>(place));
    REQUIRE(aimed.ray.status == RayStatus::Alive);
    const double dz = -100.0 - z_ep;
    REQUIRE(std::abs(aimed.ray.pos.x() - dz * std::tan(-3.0 * kDeg)) <= 1e-12);
    REQUIRE(std::abs(aimed.ray.pos.y() - dz * std::tan(4.0 * kDeg)) <= 1e-12);
    REQUIRE(aimed.ray.pos.z() == -100.0);
    const Vec3 hit = hit_on_stop(cs, aimed.ray);
    REQUIRE(std::hypot(hit.x(), hit.y()) < 1e-9);
  }
}

TEST_CASE("start plane lies before a surface that curves back upstream", "[sources]") {
  // Strongly concave first surface without aperture and a steep field: the bundle meets the
  // surface far upstream of its vertex (review of #8). Every ray must start before it.
  const MaterialLibrary lib;
  System s;
  s.name = "concave front";
  s.wavelengths = {{0.5876, 1.0, true}};
  s.aperture = {rtt::model::SystemApertureType::EntrancePupilDiameter, Param(1.0)};
  s.fields = {FieldType::AngleDeg, {{0.0, 0.0, 1.0}, {0.0, 40.0, 1.0}}};
  s.root.name = "root";
  s.root.children.push_back({Element{"L",
                                     ElementKind::Lens,
                                     Pose::along_z(10.0),
                                     "CONST:1.5168",
                                     {surface("L.S1", 0.0, -40.0), surface("L.S2", 5.0, -60.0)}}});
  Surface stop = surface("STO");
  stop.aperture = rtt::model::CircularAperture{0.6, 0.0};
  s.root.children.push_back(
      {Element{"S", ElementKind::Stop, Pose::along_z(30.0), std::nullopt, {stop}}});
  s.paths = {{"main", true, {}}};
  const CompiledSystem cs = compile(s, lib);
  const auto rays = rtt::trace::make_rays(cs, PathId{0}, std::vector<std::uint16_t>{1}, 0,
                                          rtt::trace::HexapolarPupil{2});
  for (std::size_t i = 0; i < rays.size(); ++i) {
    INFO("pupil (" << rays.pupil_x()[i] << ", " << rays.pupil_y()[i] << ")");
    REQUIRE(rays.status()[i] == RayStatus::Alive);
    RayState r;
    r.pos = Vec3(rays.pos_x()[i], rays.pos_y()[i], rays.pos_z()[i]);
    r.dir = Vec3(rays.dir_x()[i], rays.dir_y()[i], rays.dir_z()[i]);
    // The first surface lies ahead of the start point and is hit from upstream.
    const auto hit = rtt::trace::intersect_surface(r, cs.surfaces()[0]);
    REQUIRE(hit.status == rtt::geom::HitStatus::Hit);
    REQUIRE(hit.t > 0.0);
  }
  const auto chief = rtt::trace::aim_ray(cs, PathId{0}, 1, 0, 0.0, 0.0);
  REQUIRE(chief.ray.status == RayStatus::Alive);
  const Vec3 h = hit_on_stop(cs, chief.ray);
  REQUIRE(std::hypot(h.x(), h.y()) < 1e-9);
}

TEST_CASE("inverted pupil: pupil coordinates refer to the entrance pupil", "[sources]") {
  // Stop behind the focus of a positive lens (f = R / (n - 1) = 200 mm): the stop sees the
  // entrance pupil inverted, R_s < 0. Pupil point py = 1 must still start at +y on the EP side
  // and hit the stop at py R_s < 0 (decided for #8: R_s signed).
  const MaterialLibrary lib;
  System s;
  s.name = "stop behind focus";
  s.wavelengths = {{0.5876, 1.0, true}};
  s.aperture = {rtt::model::SystemApertureType::EntrancePupilDiameter, Param(4.0)};
  s.fields = {FieldType::AngleDeg, {{0.0, 0.0, 1.0}}};
  s.root.name = "root";
  s.root.children.push_back({Element{"L",
                                     ElementKind::Lens,
                                     Pose::along_z(10.0),
                                     "CONST:1.5168",
                                     {surface("L.S1", 0.0, 103.36), surface("L.S2", 3.0)}}});
  Surface stop = surface("STO");
  stop.aperture = rtt::model::CircularAperture{5.0, 0.0};
  s.root.children.push_back(
      {Element{"S", ElementKind::Stop, Pose::along_z(300.0), std::nullopt, {stop}}});
  s.paths = {{"main", true, {}}};
  const CompiledSystem cs = compile(s, lib);
  const double r_s = paraxial_stop_radius(cs);
  REQUIRE(r_s < 0.0);
  for (const Aiming aiming : {Aiming::Paraxial, Aiming::Real}) {
    const auto aimed = rtt::trace::aim_ray(cs, PathId{0}, 0, 0, 0.0, 1.0, aiming);
    REQUIRE(aimed.ray.status == RayStatus::Alive);
    REQUIRE(aimed.ray.pos.y() > 1.9);  // EP radius 2 mm at +y
  }
  const auto real = rtt::trace::aim_ray(cs, PathId{0}, 0, 0, 0.0, 1.0);
  const Vec3 hit = hit_on_stop(cs, real.ray);
  REQUIRE(std::abs(hit.y() - r_s) < 1e-9);
}

TEST_CASE("stop_size aperture: the pupil rim hits the stop at its aperture radius", "[sources]") {
  // With aperture type stop_size the stop radius defines the pupil, so (px, py) on the unit
  // circle lands on the stop rim: R_s = stop aperture radius, independent of the paraxial
  // pupil computation. Reference system tests/reference/m1/two_lenses_stop_between.rtt.json.
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(
      rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m1/two_lenses_stop_between.rtt.json"),
      lib);
  for (std::uint16_t f = 0; f < 3; ++f) {
    for (const auto& [px, py] : {std::pair{0.0, 1.0}, std::pair{1.0, 0.0}, std::pair{-0.6, -0.8}}) {
      INFO("field " << f << ", p = (" << px << ", " << py << ")");
      const auto aimed = rtt::trace::aim_ray(cs, PathId{0}, f, 0, px, py);
      REQUIRE(aimed.ray.status == RayStatus::Alive);
      const Vec3 hit = hit_on_stop(cs, aimed.ray);
      REQUIRE(std::abs(std::hypot(hit.x(), hit.y()) - 4.0) < 1e-9);
      REQUIRE(std::hypot(hit.x() - 4.0 * px, hit.y() - 4.0 * py) < 1e-9);
    }
  }
}

TEST_CASE("residual and pupil coordinates at the API boundary", "[sources]") {
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(two_lenses(StopPlace::Between), lib);
  REQUIRE_THROWS_AS(rtt::trace::aim_ray(cs, PathId{0}, 0, 0, std::nan(""), 0.0),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(rtt::trace::make_rays(
                        cs, PathId{0}, std::vector<std::uint16_t>{0}, 0,
                        rtt::trace::SinglePupilPoint{0.0, std::numeric_limits<double>::infinity()}),
                    std::invalid_argument);
  // A ray that never reaches the stop reports an infinite residual.
  System s = two_lenses(StopPlace::After);
  s.fields = {FieldType::AngleDeg, {{0.0, 0.0, 1.0}, {0.0, 80.0, 1.0}}};
  const CompiledSystem steep = compile(s, lib);
  const auto aimed = rtt::trace::aim_ray(steep, PathId{0}, 1, 0, 0.0, 1.0);
  REQUIRE(aimed.ray.status == RayStatus::NoConvergence);
  REQUIRE(std::isinf(aimed.residual));
}

TEST_CASE("no safe start plane before an unbounded surface curving back: error at the API",
          "[sources]") {
  // Concave paraboloid (k = -1) without aperture and a steep field: its sag grows with r^2
  // while the bundle bound grows linearly with the distance from the EP, so the bound keeps
  // dropping and no start plane exists (review of #8).
  const MaterialLibrary lib;
  System s;
  s.name = "concave paraboloid";
  s.wavelengths = {{0.5876, 1.0, true}};
  s.aperture = {rtt::model::SystemApertureType::EntrancePupilDiameter, Param(3.0)};
  s.fields = {FieldType::AngleDeg, {{0.0, 0.0, 1.0}, {0.0, 80.0, 1.0}}};
  s.root.name = "root";
  Surface front = surface("L.S1");
  front.shape.base = rtt::model::Conic{Param(-40.0), Param(-1.0)};
  s.root.children.push_back({Element{
      "L", ElementKind::Lens, Pose::along_z(10.0), "CONST:1.5168", {front, surface("L.S2", 5.0)}}});
  Surface stop = surface("STO");
  stop.aperture = rtt::model::CircularAperture{1.0, 0.0};
  s.root.children.push_back(
      {Element{"S", ElementKind::Stop, Pose::along_z(30.0), std::nullopt, {stop}}});
  s.paths = {{"main", true, {}}};
  const CompiledSystem cs = compile(s, lib);
  REQUIRE_THROWS_AS(rtt::trace::aim_ray(cs, PathId{0}, 1, 0, 0.0, 0.0), std::invalid_argument);
  // The on-axis field is fine.
  REQUIRE(rtt::trace::aim_ray(cs, PathId{0}, 0, 0, 0.0, 0.0).ray.status == RayStatus::Alive);
}

TEST_CASE("target beyond the reachable stop heights: NoConvergence with a finite residual",
          "[sources][aiming]") {
  // Stop behind the group, 30 deg: every iterate reaches the stop, but only at heights of about
  // 8.6 ... 18.4 mm, while the target of py = 1 is R_s = 2.88 mm. The damped Newton iteration
  // cannot reduce the residual any further and stops (backtracking exhausted).
  const MaterialLibrary lib;
  System s = two_lenses(StopPlace::After);
  s.fields = {FieldType::AngleDeg, {{0.0, 0.0, 1.0}, {0.0, 30.0, 1.0}}};
  const CompiledSystem cs = compile(s, lib);
  const double r_s = paraxial_stop_radius(cs);
  const auto aimed = rtt::trace::aim_ray(cs, PathId{0}, 1, 0, 0.0, 1.0);
  INFO("residual " << aimed.residual << ", R_s " << r_s);
  REQUIRE(aimed.ray.status == RayStatus::NoConvergence);
  REQUIRE(std::isfinite(aimed.residual));
  // The best reachable stop height lies above the target by more than 5 mm.
  REQUIRE(aimed.residual > 5.0);
  REQUIRE(aimed.ray.pos.allFinite());
}

TEST_CASE("aim_ray with a field value equals the field index and accepts other values",
          "[sources]") {
  // Overload for field values outside the model list (#31): same field type as the system.
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(two_lenses(StopPlace::Between), lib);
  const rtt::model::Field& f2 = cs.fields().points[2];
  for (const Aiming aiming : {Aiming::Real, Aiming::Paraxial}) {
    const auto by_index = rtt::trace::aim_ray(cs, PathId{0}, 2, 0, 0.3, -0.4, aiming);
    const auto by_value = rtt::trace::aim_ray(cs, PathId{0}, f2, 0, 0.3, -0.4, aiming);
    REQUIRE(by_value.ray.pos == by_index.ray.pos);
    REQUIRE(by_value.ray.dir == by_index.ray.dir);
    REQUIRE(by_value.ray.status == by_index.ray.status);
  }
  // Half of field 1 (5 deg): d ~ (0, tan 2.5 deg, 1), not in the model list.
  const auto half =
      rtt::trace::aim_ray(cs, PathId{0}, rtt::model::Field{0.0, 2.5, 1.0}, 0, 0.0, 0.0);
  REQUIRE(half.ray.status == RayStatus::Alive);
  REQUIRE(std::abs(half.ray.dir.y() / half.ray.dir.z() - std::tan(2.5 * kDeg)) <= 1e-15);
  REQUIRE(std::hypot(hit_on_stop(cs, half.ray).x(), hit_on_stop(cs, half.ray).y()) < 1e-9);
  REQUIRE_THROWS_AS(
      rtt::trace::aim_ray(cs, PathId{0}, rtt::model::Field{0.0, 95.0, 1.0}, 0, 0.0, 0.0),
      std::invalid_argument);
}

TEST_CASE("field values are converted at the reference wavelength", "[sources]") {
  // A field point is wavelength-independent: angle, object height and paraxial image height
  // describing the same point must give the same start (direction for an object at infinity,
  // object point otherwise) for every wavelength (fix of #8, found in #31). Dispersive N-BK7
  // lenses with the stop between them make entrance pupil and image heights depend on lambda.
  rtt::material::MaterialLibrary lib;
  lib.add_catalog(std::string(RTT_CATALOG_DIR) + "/schott.agf");
  const auto dispersive = [](System s) {
    s.wavelengths = {{0.4861, 1.0, false}, {0.5876, 1.0, true}, {0.6563, 1.0, false}};
    for (auto& child : s.root.children) {
      if (auto* e = std::get_if<Element>(&child.value); e && e->kind == ElementKind::Lens) {
        e->material = "SCHOTT:N-BK7";
      }
    }
    return s;
  };
  constexpr double kTheta = 4.0;  // degree
  SECTION("object at infinity: angle and paraxial image height") {
    // y' = EFL tan(theta) needs n = 1 in object and image space; Ciddor AIR (since #25) gives
    // y' = n_air EFL tan(theta), so the environment is VACUUM here.
    System angle = dispersive(two_lenses(StopPlace::Between));
    angle.environment.medium = "VACUUM";
    const CompiledSystem ca = compile(angle, lib);
    const auto fo = rtt::paraxial::first_order(ca, PathId{0}, ca.reference_wavelength());
    REQUIRE(fo.efl);
    // Paraxial image height of a distant object at the reference wavelength: EFL tan(theta).
    System image = angle;
    image.fields = {FieldType::ParaxialImageHeight,
                    {{0.0, 0.0, 1.0}, {0.0, *fo.efl * std::tan(kTheta * kDeg), 1.0}}};
    const CompiledSystem ci = compile(image, lib);
    const Vec3 expected = Vec3(0.0, std::tan(kTheta * kDeg), 1.0).normalized();
    for (std::uint16_t wl = 0; wl < 3; ++wl) {
      INFO("wavelength " << wl);
      const auto a =
          rtt::trace::aim_ray(ca, PathId{0}, rtt::model::Field{0.0, kTheta, 1.0}, wl, 0.0, 0.0);
      const auto i = rtt::trace::aim_ray(ci, PathId{0}, 1, wl, 0.0, 0.0);
      REQUIRE((a.ray.dir - expected).cwiseAbs().maxCoeff() <= 1e-15);
      REQUIRE((i.ray.dir - expected).cwiseAbs().maxCoeff() <= 1e-12);
    }
  }
  SECTION("finite object: angle, object height and paraxial image height") {
    System height = dispersive(two_lenses(StopPlace::Between, true));
    const CompiledSystem ch = compile(height, lib);
    const std::uint16_t ref = ch.reference_wavelength();
    const auto fo = rtt::paraxial::first_order(ch, PathId{0}, ref);
    REQUIRE(fo.lateral_magnification);
    const double z_ep = *fo.entrance_pupil->z;
    // Object point of the angle field: on the chief ray through the reference EP centre.
    const double h = (-100.0 - z_ep) * std::tan(kTheta * kDeg);
    height.fields = {FieldType::ObjectHeight, {{0.0, 0.0, 1.0}, {0.0, h, 1.0}}};
    const CompiledSystem c_height = compile(height, lib);
    System angle = height;
    angle.fields = {FieldType::AngleDeg, {{0.0, 0.0, 1.0}, {0.0, kTheta, 1.0}}};
    const CompiledSystem c_angle = compile(angle, lib);
    System image = height;
    image.fields = {FieldType::ParaxialImageHeight,
                    {{0.0, 0.0, 1.0}, {0.0, *fo.lateral_magnification * h, 1.0}}};
    const CompiledSystem c_image = compile(image, lib);
    for (std::uint16_t wl = 0; wl < 3; ++wl) {
      INFO("wavelength " << wl);
      const Vec3 expected(0.0, h, -100.0);
      REQUIRE((rtt::trace::aim_ray(c_height, PathId{0}, 1, wl, 0.0, 0.0).ray.pos - expected)
                  .cwiseAbs()
                  .maxCoeff() <= 1e-12);
      REQUIRE((rtt::trace::aim_ray(c_angle, PathId{0}, 1, wl, 0.0, 0.0).ray.pos - expected)
                  .cwiseAbs()
                  .maxCoeff() <= 1e-12);
      REQUIRE((rtt::trace::aim_ray(c_image, PathId{0}, 1, wl, 0.0, 0.0).ray.pos - expected)
                  .cwiseAbs()
                  .maxCoeff() <= 1e-12);
    }
  }
}

TEST_CASE("rim rays of the reference singlet pass the stop", "[sources]") {
  // Rays with |p| = 1 are aimed at the stop rim to within kAimTolerance; the aperture check
  // includes the same tolerance (#50), so none of them is vignetted at the stop (previously 7 of
  // the 36 rim rays per field were). The lens apertures (12.7 mm) do not vignette these fields.
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(singlet(), lib);
  const std::vector<std::uint16_t> fields{0, 1, 2};
  rtt::trace::RayBatch rays =
      rtt::trace::make_rays(cs, PathId{0}, fields, 1, rtt::trace::HexapolarPupil{6});
  [[maybe_unused]] const auto stats = rtt::trace::SequentialTracer().trace(cs, PathId{0}, rays);
  std::size_t rim = 0;
  for (std::size_t i = 0; i < rays.size(); ++i) {
    INFO("field " << rays.field()[i] << ", pupil (" << rays.pupil_x()[i] << ", "
                  << rays.pupil_y()[i] << ")");
    REQUIRE(rays.status()[i] == RayStatus::Alive);
    if (std::hypot(rays.pupil_x()[i], rays.pupil_y()[i]) > 1.0 - 1e-12) ++rim;
  }
  REQUIRE(rim == 3 * 36);
}

TEST_CASE("Cooke triplet: every field has a start plane before the first surface (#72)",
          "[sources]") {
  // tests/reference/m2/cooke_triplet.rtt.json (#34): object at infinity, fields 0, 14 and 20
  // degree, no apertures on the lens surfaces, stop 1e-8 mm behind the vertex of L2.S1. Before
  // #72 field 1 threw "no start plane before surface 'L2.S1'" (found in #61).
  MaterialLibrary lib;
  lib.add_catalog(std::string(RTT_CATALOG_DIR) + "/m2/schott.agf");
  const CompiledSystem cs = compile(
      rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m2/cooke_triplet.rtt.json"), lib);
  const auto vertex_z = [&](const std::string& id) {
    for (const auto& s : cs.surfaces()) {
      if (s.id.str() == id) return s.to_global.translation().z();
    }
    FAIL("no surface " << id);
    return 0.0;
  };
  // L1.S1 is convex towards -z, so its vertex is its lowest point.
  const double z_first = vertex_z("L1.S1");
  // Field 1: the bound for L2.S1 (sphere R = -24.456 mm, no aperture) converges only linearly
  // (factor tan(14 deg) |dsag/dr| ~ 0.71), so after 50 rounds the whole cap disk bounds the
  // bundle: its lowest point is the bottom of the hemisphere, vertex z - 24.456 mm.
  const double z_hemisphere = vertex_z("L2.S1") - 24.456;
  for (std::uint16_t field = 0; field < 3; ++field) {
    for (std::uint16_t wl = 0; wl < cs.wavelengths_um().size(); ++wl) {
      for (const auto& [px, py] :
           {std::pair{0.0, 0.0}, std::pair{0.0, 1.0}, std::pair{0.0, -1.0}, std::pair{1.0, 0.0}}) {
        INFO("field " << field << ", wavelength " << wl << ", pupil " << px << " " << py);
        const auto aimed = rtt::trace::aim_ray(cs, PathId{0}, field, wl, px, py);
        REQUIRE(aimed.ray.status == RayStatus::Alive);
        // At least 1 mm before L1.S1 (start plane rule, #8).
        REQUIRE(aimed.ray.pos.z() <= z_first - 1.0);
        if (field == 1) REQUIRE(aimed.ray.pos.z() <= z_hemisphere - 1.0 + 1e-9);
      }
    }
  }
  // The spot of fields 1 and 2 with real aiming (default sampling of spot()): every chief ray
  // reaches the detector, and all rays start at least 1 mm before L1.S1.
  const std::vector<std::uint16_t> fields{1, 2};
  for (std::uint16_t wl = 0; wl < cs.wavelengths_um().size(); ++wl) {
    INFO("wavelength " << wl);
    rtt::trace::RayBatch rays = rtt::trace::make_rays(cs, PathId{0}, fields, wl,
                                                      rtt::trace::HexapolarPupil{}, Aiming::Real);
    for (std::size_t i = 0; i < rays.size(); ++i) REQUIRE(rays.pos_z()[i] <= z_first - 1.0);
    [[maybe_unused]] const auto stats = rtt::trace::SequentialTracer().trace(cs, PathId{0}, rays);
    std::size_t chief = 0;
    for (std::size_t i = 0; i < rays.size(); ++i) {
      if (rays.pupil_x()[i] == 0.0 && rays.pupil_y()[i] == 0.0) {
        ++chief;
        REQUIRE(rays.status()[i] == RayStatus::Alive);
      }
    }
    REQUIRE(chief == 2);
  }
}
