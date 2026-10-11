// The ideal lens and the ideal cylinder lens in the tracer (#178, ADR 0031, points 2 to 5, tests
// of point 9 for PR 2): perfect imaging of every ray (collinear map), the optical path exact for
// the design conjugate, P = R(k_in -> k_out) with the weight unchanged, the backward passage.
//
// Setup of tests/reference/r2/ideal_lens.rtt.json: AIR, f = 50 mm at z = 0, object plane at
// z = -75 (s = -75), stop r = 10 mm at z = -1, detector at z = 150 = s' with 1/s' = 1/s + 1/f,
// magnification m = s'/s = -2. Start rays from the object point through points of the stop,
// so the tests do not depend on ray aiming (PR 3).
//
// Tolerances (ADR 0031, point 9): a step in double rounds by at most u |value|, u = 2^-53, so
// coordinates up to 300 mm carry some 1e-14 mm per step and a few dozen steps at most 1e-12 mm;
// 1e-9 mm is more than 1000 times that. OPL up to 250 mm: 3e-14 mm per step, at most 3e-12 mm
// = 5e-9 waves at 0.5876 um; 1e-6 waves is more than 100 times that. Directions: 1e-15.

#include <Eigen/Geometry>
#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <complex>
#include <cstddef>
#include <numbers>
#include <optional>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/math/types.hpp"
#include "rtt/model/model.hpp"
#include "rtt/trace/apply_event.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/sequential.hpp"

using rtt::compile::CompiledSystem;
using rtt::math::Vec3;
using rtt::model::Element;
using rtt::model::IdealCylinderLens;
using rtt::model::IdealLens;
using rtt::model::Param;
using rtt::model::System;
using rtt::trace::RayBatch;
using rtt::trace::RayStatus;

namespace {

constexpr double kLambdaMm = 0.5876e-3;  ///< reference wavelength in mm
constexpr double kOplTolerance = 1e-6 * kLambdaMm;
constexpr double kPositionTolerance = 1e-9;

System reference() {
  return rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/r2/ideal_lens.rtt.json");
}

/// Children of r2/ideal_lens: 0 stop, 1 ideal lens IL, 2 cylinder lens CL, 3 image.
Element& child(System& s, std::size_t i) {
  return std::get<Element>(s.root.children[i].value);
}
IdealLens& ideal(System& s) {
  return std::get<IdealLens>(child(s, 1).surfaces[0].interaction);
}
IdealCylinderLens& cylinder(System& s) {
  return std::get<IdealCylinderLens>(child(s, 2).surfaces[0].interaction);
}

CompiledSystem compiled(const System& s) {
  return rtt::compile::compile(s, rtt::material::MaterialLibrary{});
}

/// Points of the stop disk (radius 9.9 mm at z = -1): centre and 6 rings of 6 k points.
std::vector<Vec3> stop_points() {
  std::vector<Vec3> out{Vec3(0.0, 0.0, -1.0)};
  for (int k = 1; k <= 6; ++k) {
    const double r = 9.9 * k / 6.0;
    for (int j = 0; j < 6 * k; ++j) {
      const double phi = 2.0 * std::numbers::pi * j / (6.0 * k);
      out.emplace_back(r * std::cos(phi), r * std::sin(phi), -1.0);
    }
  }
  return out;
}

/// Rays from `object` (OPL 0) towards each point of the stop.
RayBatch rays_from(const Vec3& object) {
  const std::vector<Vec3> targets = stop_points();
  RayBatch rays(targets.size());
  for (std::size_t i = 0; i < targets.size(); ++i) {
    const Vec3 d = (targets[i] - object).normalized();
    rays.pos_x()[i] = object.x();
    rays.pos_y()[i] = object.y();
    rays.pos_z()[i] = object.z();
    rays.dir_x()[i] = d.x();
    rays.dir_y()[i] = d.y();
    rays.dir_z()[i] = d.z();
  }
  return rays;
}

void trace(const CompiledSystem& cs, const std::string& path, RayBatch& rays) {
  static_cast<void>(rtt::trace::SequentialTracer{}.trace(cs, *cs.find_path(path), rays));
}

Vec3 pos(const RayBatch& r, std::size_t i) {
  return {r.pos_x()[i], r.pos_y()[i], r.pos_z()[i]};
}
Vec3 dir(const RayBatch& r, std::size_t i) {
  return {r.dir_x()[i], r.dir_y()[i], r.dir_z()[i]};
}

/// The point of the ray (p, d) in the plane z = `z` (forwards or backwards).
Vec3 at_z(const Vec3& p, const Vec3& d, double z) {
  return p + ((z - p.z()) / d.z()) * d;
}

/// max - min of the OPL of all rays.
double opl_span(const RayBatch& rays) {
  double lo = rays.opl()[0];
  double hi = lo;
  for (std::size_t i = 0; i < rays.size(); ++i) {
    lo = std::min(lo, rays.opl()[i]);
    hi = std::max(hi, rays.opl()[i]);
  }
  return hi - lo;
}

void all_alive(const RayBatch& rays) {
  for (std::size_t i = 0; i < rays.size(); ++i) REQUIRE(rays.status()[i] == RayStatus::Alive);
}

/// Media of a transmit in the environment of `cs` at the reference wavelength.
rtt::trace::EventMedia environment(const CompiledSystem& cs) {
  rtt::trace::EventMedia m;
  m.before = m.after = m.beyond = cs.media()[0].index[cs.reference_wavelength()];
  m.wavelength_um = 0.5876;
  return m;
}

const rtt::compile::CompiledSurface& surface(const CompiledSystem& cs, const char* id) {
  return cs.surfaces()[*cs.find_surface(rtt::model::SurfaceId{id})];
}

}  // namespace

// ------------------------------------------------------------------ point 2: directions -----

TEST_CASE("ideal lens: every ray of an object point meets in its image, m = -2", "[ideal][trace]") {
  const CompiledSystem cs = compiled(reference());
  const double s = -75.0;
  const double s_image = 1.0 / (1.0 / 50.0 + 1.0 / s);  // 150 up to rounding
  const double m = s_image / s;
  for (const Vec3& object :
       {Vec3(0.0, 0.0, -75.0), Vec3(0.0, 5.0, -75.0), Vec3(3.0, -4.0, -75.0)}) {
    INFO("object " << object.transpose());
    RayBatch rays = rays_from(object);
    trace(cs, "ideal lens", rays);
    all_alive(rays);
    const Vec3 image(m * object.x(), m * object.y(), s_image);
    for (std::size_t i = 0; i < rays.size(); ++i) {
      CHECK((at_z(pos(rays, i), dir(rays, i), s_image) - image).norm() < kPositionTolerance);
    }
  }
}

TEST_CASE("ideal lens: f = -50 gives a virtual image at s' = -30", "[ideal][trace]") {
  System sys = reference();
  ideal(sys).focal_length = Param(-50.0);
  const CompiledSystem cs = compiled(sys);
  const double s_image = 1.0 / (1.0 / -50.0 + 1.0 / -75.0);  // -30
  const Vec3 object(3.0, -4.0, -75.0);
  const Vec3 image((s_image / -75.0) * object.x(), (s_image / -75.0) * object.y(), s_image);
  RayBatch rays = rays_from(object);
  trace(cs, "ideal lens", rays);
  all_alive(rays);
  for (std::size_t i = 0; i < rays.size(); ++i) {
    // The rays leave diverging; extended backwards they meet in the virtual image.
    CHECK((at_z(pos(rays, i), dir(rays, i), s_image) - image).norm() < kPositionTolerance);
  }
}

TEST_CASE("ideal lens: a plane wave at 10 deg meets in (0, f tan 10 deg, f)", "[ideal][trace]") {
  System sys = reference();
  ideal(sys).object_distance.reset();  // design conjugate at infinity
  child(sys, 3).pose.position[2] = Param(50.0);
  const CompiledSystem cs = compiled(sys);
  const double a = 10.0 * std::numbers::pi / 180.0;
  const Vec3 d(0.0, std::sin(a), std::cos(a));
  // Start points on the plane through (0, 0, -20) perpendicular to d, OPL 0 there (a plane
  // wave): offsets along x and along e = d x x.
  const Vec3 e = d.cross(Vec3::UnitX());
  std::vector<Vec3> starts;
  for (int i = -4; i <= 4; ++i) {
    for (int j = -4; j <= 4; ++j) {
      starts.push_back(Vec3(0.0, 0.0, -20.0) + 2.0 * i * Vec3::UnitX() + 2.0 * j * e);
    }
  }
  RayBatch rays(starts.size());
  for (std::size_t i = 0; i < starts.size(); ++i) {
    rays.pos_x()[i] = starts[i].x();
    rays.pos_y()[i] = starts[i].y();
    rays.pos_z()[i] = starts[i].z();
    rays.dir_x()[i] = d.x();
    rays.dir_y()[i] = d.y();
    rays.dir_z()[i] = d.z();
  }
  trace(cs, "ideal lens", rays);
  const Vec3 focus(0.0, 50.0 * std::tan(a), 50.0);
  std::size_t arrived = 0;
  for (std::size_t i = 0; i < rays.size(); ++i) {
    if (rays.status()[i] != RayStatus::Alive) continue;  // outside the stop
    ++arrived;
    CHECK((pos(rays, i) - focus).norm() < kPositionTolerance);
  }
  CHECK(arrived > 20);  // content guard
  // The plane wave arrives with equal OPL in F' (ADR 0031, point 3, Delta_infinity).
  RayBatch alive(arrived);
  std::size_t k = 0;
  for (std::size_t i = 0; i < rays.size(); ++i) {
    if (rays.status()[i] == RayStatus::Alive) alive.opl()[k++] = rays.opl()[i];
  }
  CHECK(opl_span(alive) < kOplTolerance);
}

// ---------------------------------------------------------------- point 3: optical path -----

TEST_CASE("ideal lens: equal OPL from an object point of the design plane to its image",
          "[ideal][trace][opl]") {
  const CompiledSystem cs = compiled(reference());
  for (const Vec3& object :
       {Vec3(0.0, 0.0, -75.0), Vec3(0.0, 5.0, -75.0), Vec3(3.0, -4.0, -75.0)}) {
    INFO("object " << object.transpose());
    RayBatch rays = rays_from(object);
    trace(cs, "ideal lens", rays);
    all_alive(rays);
    // The detector lies in the image plane: every ray ends in the image point.
    CHECK(opl_span(rays) < kOplTolerance);
  }
}

TEST_CASE("ideal lens: without the design conjugate the OPL differs by more than 1 wave",
          "[ideal][trace][opl]") {
  // Content guard (ADR 0031, point 3): design at infinity, object at 75 mm: about 11 waves at
  // the rim (n r^4 (2/s + 1/s') / (8 f s)), the directions stay perfect.
  System sys = reference();
  ideal(sys).object_distance.reset();
  const CompiledSystem cs = compiled(sys);
  RayBatch rays = rays_from(Vec3(0.0, 0.0, -75.0));
  trace(cs, "ideal lens", rays);
  all_alive(rays);
  CHECK(opl_span(rays) > kLambdaMm);
  for (std::size_t i = 0; i < rays.size(); ++i) {
    CHECK(at_z(pos(rays, i), dir(rays, i), 150.0).head<2>().norm() < kPositionTolerance);
  }
}

TEST_CASE("ideal lens: the OPL jump is the eikonal of the direction change",
          "[ideal][trace][opl]") {
  // For the rays of a point of the design plane Delta is a function W(q) of the hit point, and
  // grad W = n (t'_par - t_par) (ADR 0031, point 3). Central difference with h = 1e-3 mm:
  // truncation h^2 |W'''| / 6 with |W'''| <~ 3 n / f^2 is about 2e-10. delta() forms
  // after.opl - n t with |opl| about 75 mm, so the rounding is about u 75 / h = 1e-11; tolerance
  // 1e-8.
  const CompiledSystem cs = compiled(reference());
  const rtt::compile::CompiledSurface& lens = surface(cs, "IL");
  const std::uint32_t index = *cs.find_surface(rtt::model::SurfaceId{"IL"});
  const rtt::trace::EventMedia media = environment(cs);
  const double n = media.before.real();
  const Vec3 object(0.0, 5.0, -75.0);
  // Delta of the ray from `object` to the point q of the lens plane; out = state after.
  const auto delta = [&](double qx, double qy, rtt::trace::RayState* out) {
    rtt::trace::RayState ray;
    ray.pos = object;
    ray.dir = (Vec3(qx, qy, 0.0) - object).normalized();
    const rtt::trace::SurfaceHit hit = rtt::trace::intersect_surface(ray, lens);
    const rtt::trace::RayState after =
        rtt::trace::apply_event(ray, lens, hit, index, rtt::model::EventKind::Transmit, 0, media);
    REQUIRE(after.status == RayStatus::Alive);
    if (out != nullptr) *out = after;
    return after.opl - n * hit.t;
  };
  const double h = 1e-3;
  for (const auto& [qx, qy] : {std::pair{0.0, 0.0}, std::pair{3.0, -2.0}, std::pair{-6.0, 7.0}}) {
    INFO("q = (" << qx << ", " << qy << ")");
    rtt::trace::RayState after;
    static_cast<void>(delta(qx, qy, &after));
    const Vec3 t_in = (Vec3(qx, qy, 0.0) - object).normalized();
    const double gx = (delta(qx + h, qy, nullptr) - delta(qx - h, qy, nullptr)) / (2.0 * h);
    const double gy = (delta(qx, qy + h, nullptr) - delta(qx, qy - h, nullptr)) / (2.0 * h);
    CHECK(std::abs(gx - n * (after.dir.x() - t_in.x())) < 1e-8);
    CHECK(std::abs(gy - n * (after.dir.y() - t_in.y())) < 1e-8);
  }
}

// --------------------------------------------------------------- point 4: cylinder lens -----

TEST_CASE("ideal cylinder lens: focuses y only, x stays the straight line", "[ideal][trace]") {
  const CompiledSystem cs = compiled(reference());
  const Vec3 object(3.0, 2.0, -75.0);
  const RayBatch start = rays_from(object);
  RayBatch rays = start;
  trace(cs, "cylinder lens", rays);
  all_alive(rays);
  for (std::size_t i = 0; i < rays.size(); ++i) {
    // psi = 0: power along y, image line y = m y0 = -4 at z = 150.
    CHECK(std::abs(at_z(pos(rays, i), dir(rays, i), 150.0).y() + 4.0) < kPositionTolerance);
    // x is that of the straight line without the lens in every plane behind it.
    for (const double z : {50.0, 150.0, 300.0}) {
      const double free = at_z(object, dir(start, i), z).x();
      CHECK(std::abs(at_z(pos(rays, i), dir(rays, i), z).x() - free) < kPositionTolerance);
    }
  }
}

TEST_CASE("ideal cylinder lens: a ray with b . q = 0 keeps its direction", "[ideal][trace]") {
  const CompiledSystem cs = compiled(reference());
  const rtt::compile::CompiledSurface& lens = surface(cs, "CL");
  const std::uint32_t index = *cs.find_surface(rtt::model::SurfaceId{"CL"});
  rtt::trace::RayState ray;
  ray.pos = Vec3(0.0, 0.0, -75.0);
  ray.dir = (Vec3(7.0, 0.0, 0.0) - ray.pos).normalized();  // hits y = 0
  const rtt::trace::SurfaceHit hit = rtt::trace::intersect_surface(ray, lens);
  const rtt::trace::RayState out = rtt::trace::apply_event(
      ray, lens, hit, index, rtt::model::EventKind::Transmit, 0, environment(cs));
  REQUIRE(out.status == RayStatus::Alive);
  CHECK((out.dir - ray.dir).norm() < 1e-15);
}

TEST_CASE("ideal cylinder lens: psi = 90 deg focuses x; psi = 30 deg is psi = 0 turned",
          "[ideal][trace]") {
  const Vec3 object(3.0, 2.0, -75.0);
  SECTION("psi = 90 deg: image line x = -6") {
    System sys = reference();
    cylinder(sys).axis_deg = 90.0;
    const CompiledSystem cs = compiled(sys);
    RayBatch rays = rays_from(object);
    trace(cs, "cylinder lens", rays);
    all_alive(rays);
    for (std::size_t i = 0; i < rays.size(); ++i) {
      CHECK(std::abs(at_z(pos(rays, i), dir(rays, i), 150.0).x() + 6.0) < kPositionTolerance);
    }
  }
  SECTION("psi = 30 deg against psi = 0 in coordinates turned by 30 deg") {
    System turned = reference();
    cylinder(turned).axis_deg = 30.0;
    const CompiledSystem cs30 = compiled(turned);
    const CompiledSystem cs0 = compiled(reference());
    const double a = 30.0 * std::numbers::pi / 180.0;
    const rtt::math::Mat3 rot = Eigen::AngleAxisd(a, Vec3::UnitZ()).toRotationMatrix();
    RayBatch rays30 = rays_from(object);
    RayBatch rays0 = rays_from(rot.transpose() * object);
    trace(cs30, "cylinder lens", rays30);
    trace(cs0, "cylinder lens", rays0);
    all_alive(rays30);
    all_alive(rays0);
    for (std::size_t i = 0; i < rays30.size(); ++i) {
      // The stop points are not turned, so compare the image line of each system: b . X = m b . P
      // at z = 150 with b = (-sin psi, cos psi).
      const Vec3 b(-std::sin(a), std::cos(a), 0.0);
      const double bx = b.dot(at_z(pos(rays30, i), dir(rays30, i), 150.0));
      CHECK(std::abs(bx - (-2.0) * b.dot(object)) < kPositionTolerance);
      const double y0 = at_z(pos(rays0, i), dir(rays0, i), 150.0).y();
      CHECK(std::abs(y0 - (-2.0) * (rot.transpose() * object).y()) < kPositionTolerance);
    }
  }
}

TEST_CASE("ideal cylinder lens: crossed cylinder lenses act as the ideal lens", "[ideal][trace]") {
  // psi = 0 and psi = 90 deg at the same place: the slope changes add up (ADR 0031, point 4).
  System crossed = reference();
  crossed.root.children.push_back(crossed.root.children[2]);  // a second CL, renamed below
  Element& second = std::get<Element>(crossed.root.children.back().value);
  second.name = "cylinder lens 90";
  second.surfaces[0].id = rtt::model::SurfaceId{"CL90"};
  std::get<IdealCylinderLens>(second.surfaces[0].interaction).axis_deg = 90.0;
  const CompiledSystem cs = compiled(crossed);
  const CompiledSystem reference_cs = compiled(reference());
  const rtt::trace::EventMedia media = environment(cs);
  const auto& cl = surface(cs, "CL");
  const auto& cl90 = surface(cs, "CL90");
  const auto& il = surface(reference_cs, "IL");
  for (const Vec3& q : {Vec3(4.0, -3.0, 0.0), Vec3(-7.0, 2.0, 0.0), Vec3(0.5, 9.0, 0.0)}) {
    INFO("q = " << q.transpose());
    rtt::trace::RayState ray;
    ray.pos = Vec3(3.0, -4.0, -75.0);
    ray.dir = (q - ray.pos).normalized();
    const auto hit = rtt::trace::intersect_surface(ray, cl);
    const auto a =
        rtt::trace::apply_event(ray, cl, hit, 0, rtt::model::EventKind::Transmit, 0, media);
    // The second lens at the same place: the same hit point, the new direction.
    rtt::trace::SurfaceHit hit2 = hit;
    hit2.t = 0.0;
    hit2.direction = cl90.to_local.apply_vector(a.dir);
    const auto b =
        rtt::trace::apply_event(a, cl90, hit2, 1, rtt::model::EventKind::Transmit, 0, media);
    const auto hit_il = rtt::trace::intersect_surface(ray, il);
    const auto c =
        rtt::trace::apply_event(ray, il, hit_il, 0, rtt::model::EventKind::Transmit, 0, media);
    REQUIRE((a.status == RayStatus::Alive && b.status == RayStatus::Alive));
    // Same hit point and same direction as the ideal lens, whose rays meet in the image within
    // 1e-9 mm (first test): the spot of the crossed lenses follows (ADR 0031, point 9).
    CHECK((b.pos - c.pos).norm() < 1e-12);
    CHECK((b.dir - c.dir).norm() < 1e-15);
  }
}

TEST_CASE("ideal cylinder lens: OPL exact in the power section, zero along the axis",
          "[ideal][trace][opl]") {
  const CompiledSystem cs = compiled(reference());
  const rtt::compile::CompiledSurface& lens = surface(cs, "CL");
  const rtt::trace::EventMedia media = environment(cs);
  const double n = media.before.real();
  const Vec3 object(0.0, 0.0, -75.0);
  std::vector<double> b_fan;
  for (int k = -6; k <= 6; ++k) {
    rtt::trace::RayState ray;
    ray.pos = object;
    // b fan (along y): OPL from the object to the image line at z = 150 is the same for all.
    ray.dir = (Vec3(0.0, 1.5 * k, 0.0) - object).normalized();
    const auto hit = rtt::trace::intersect_surface(ray, lens);
    const auto out =
        rtt::trace::apply_event(ray, lens, hit, 0, rtt::model::EventKind::Transmit, 0, media);
    REQUIRE(out.status == RayStatus::Alive);
    const Vec3 end = at_z(out.pos, out.dir, 150.0);
    b_fan.push_back(out.opl + n * (end - out.pos).norm());
    // a fan (along x): Delta = 0, the lens does nothing there.
    rtt::trace::RayState ray_a;
    ray_a.pos = object;
    ray_a.dir = (Vec3(1.5 * k, 0.0, 0.0) - object).normalized();
    const auto hit_a = rtt::trace::intersect_surface(ray_a, lens);
    const auto out_a =
        rtt::trace::apply_event(ray_a, lens, hit_a, 0, rtt::model::EventKind::Transmit, 0, media);
    CHECK(std::abs(out_a.opl - n * hit_a.t) < 1e-15);
  }
  double lo = b_fan[0];
  double hi = lo;
  for (const double v : b_fan) {
    lo = std::min(lo, v);
    hi = std::max(hi, v);
  }
  CHECK(hi - lo < kOplTolerance);
}

// ------------------------------------------------------------ point 5: P and the weight -----

TEST_CASE("ideal lens: P is the rotation k_in -> k_out, the weight stays", "[ideal][trace]") {
  const CompiledSystem cs = compiled(reference());
  const rtt::compile::CompiledSurface& lens = surface(cs, "IL");
  const rtt::trace::EventMedia media = environment(cs);
  SECTION("off the centre") {
    for (const Vec3& q : {Vec3(4.0, -3.0, 0.0), Vec3(-9.0, 1.0, 0.0)}) {
      rtt::trace::RayState ray;
      ray.pos = Vec3(3.0, -4.0, -75.0);
      ray.dir = (q - ray.pos).normalized();
      const auto hit = rtt::trace::intersect_surface(ray, lens);
      const auto out =
          rtt::trace::apply_event(ray, lens, hit, 0, rtt::model::EventKind::Transmit, 0, media);
      REQUIRE(out.status == RayStatus::Alive);
      const rtt::math::CMat3& p = out.prt;  // P_event, since the incoming P is I
      CHECK((p * ray.dir.cast<rtt::math::Complex>() - out.dir.cast<rtt::math::Complex>()).norm() <
            1e-15);
      CHECK((p.transpose() * p - rtt::math::CMat3::Identity()).norm() < 1e-15);
      CHECK(std::abs(out.weight - 1.0) < 1e-15);
    }
  }
  SECTION("at the centre P = I bitwise") {
    rtt::trace::RayState ray;
    ray.pos = Vec3(0.0, 0.0, -75.0);
    const auto hit = rtt::trace::intersect_surface(ray, lens);
    const auto out =
        rtt::trace::apply_event(ray, lens, hit, 0, rtt::model::EventKind::Transmit, 0, media);
    CHECK(out.prt == rtt::math::CMat3::Identity());
    CHECK(out.weight == 1.0);
    CHECK(out.dir == ray.dir);
  }
  SECTION("normal incidence at (0, y): x polarization stays x") {
    rtt::trace::RayState ray;
    ray.pos = Vec3(0.0, 6.0, -75.0);
    const auto hit = rtt::trace::intersect_surface(ray, lens);
    const auto out =
        rtt::trace::apply_event(ray, lens, hit, 0, rtt::model::EventKind::Transmit, 0, media);
    const Eigen::Vector3cd ex = out.prt * Vec3::UnitX().cast<rtt::math::Complex>();
    CHECK((ex - Vec3::UnitX().cast<rtt::math::Complex>()).norm() < 1e-15);
  }
}

// --------------------------------------------------------------------- backward passage -----

TEST_CASE("ideal lens: the backward passage is the mirrored forward passage", "[ideal][trace]") {
  // A ray in -z sees the lens mirrored (z -> -z): the same slopes, the same q, the same f, so
  // the same deflection and the same Delta (ADR 0031, "Bezeichnungen").
  const CompiledSystem cs = compiled(reference());
  const rtt::compile::CompiledSurface& lens = surface(cs, "IL");
  const rtt::trace::EventMedia media = environment(cs);
  const Vec3 object(3.0, -4.0, -75.0);
  const Vec3 mirrored_object(3.0, -4.0, 75.0);
  for (const Vec3& q : {Vec3(4.0, -3.0, 0.0), Vec3(-9.0, 1.0, 0.0), Vec3(0.0, 0.0, 0.0)}) {
    INFO("q = " << q.transpose());
    rtt::trace::RayState forward;
    forward.pos = object;
    forward.dir = (q - object).normalized();
    rtt::trace::RayState backward;
    backward.pos = mirrored_object;
    backward.dir = (q - mirrored_object).normalized();
    const auto hf = rtt::trace::intersect_surface(forward, lens);
    const auto hb = rtt::trace::intersect_surface(backward, lens);
    const auto f =
        rtt::trace::apply_event(forward, lens, hf, 0, rtt::model::EventKind::Transmit, 0, media);
    const auto b =
        rtt::trace::apply_event(backward, lens, hb, 0, rtt::model::EventKind::Transmit, 0, media);
    REQUIRE((f.status == RayStatus::Alive && b.status == RayStatus::Alive));
    CHECK(std::abs(b.dir.x() - f.dir.x()) < 1e-15);
    CHECK(std::abs(b.dir.y() - f.dir.y()) < 1e-15);
    CHECK(std::abs(b.dir.z() + f.dir.z()) < 1e-15);
    CHECK(std::abs(b.opl - f.opl) < 1e-12);
    // The backward image lies mirrored at z = -150 (local), up to 1e-9 mm.
    const Vec3 image = at_z(b.pos, b.dir, -150.0);
    CHECK((image - Vec3(-6.0, 8.0, -150.0)).norm() < kPositionTolerance);
  }
}

TEST_CASE("ideal lens: lens -> mirror -> lens in the tracer", "[ideal][trace]") {
  // ADR 0031, point 9 (Rueckweg): a plane ideal mirror M at z = d = 25 behind IL, the path
  // STO, IL, M, IL. The forward image at z = s' = 150 is mirrored to z = 2 d - 150 = -100, a
  // virtual object s_b = 100 mm behind the lens along the backward propagation; its image lies
  // at s_b' = 1 / (1/f + 1/s_b) along -z with m_b = s_b'/s_b. The tracer hits IL from behind with
  // the reflected local direction; the second passage must equal the mirrored forward passage
  // of the same ray (apply_event at the lens in its basic pose), Delta included. The backward
  // object s_b = 100 is not the design conjugate (object_distance 75), so the OPL of the
  // backward image is not equal for all rays (the r^4 rest of point 3); hence the comparison
  // with the mirrored passage and no OPL span.
  const double d = 25.0;
  const double f = 50.0;
  System s = reference();
  rtt::model::Node mirror_node = s.root.children[1];  // a copy of the thin element IL
  Element& mirror = std::get<Element>(mirror_node.value);
  mirror.kind = rtt::model::ElementKind::Mirror;
  mirror.name = "mirror";
  mirror.pose.position[2] = Param(d);
  mirror.surfaces[0].id = rtt::model::SurfaceId{"M"};
  mirror.surfaces[0].interaction = rtt::model::IdealMirror{};
  s.root.children.push_back(mirror_node);
  const double s_back = 2.0 * d - 150.0;  // z of the mirrored forward image: -100
  const double sb = -s_back;              // virtual object 100 mm along the backward propagation
  const double sb_image = 1.0 / (1.0 / f + 1.0 / sb);
  rtt::model::Node detector_node = s.root.children[3];  // a copy of the detector IMG
  Element& detector = std::get<Element>(detector_node.value);
  detector.name = "back image";
  detector.pose.position[2] = Param(-sb_image);
  detector.surfaces[0].id = rtt::model::SurfaceId{"BIMG"};
  s.root.children.push_back(detector_node);
  using rtt::model::EventKind;
  using rtt::model::SurfaceId;
  const auto path = [](std::string name, std::vector<rtt::model::Event> events) {
    return rtt::model::Path{std::move(name), false, std::move(events)};
  };
  const rtt::model::Event sto{SurfaceId{"STO"}, EventKind::Transmit, 0};
  const rtt::model::Event il{SurfaceId{"IL"}, EventKind::Transmit, 0};
  const rtt::model::Event m{SurfaceId{"M"}, EventKind::Reflect, 0};
  s.paths.push_back(path("to mirror", {sto, il, m}));
  s.paths.push_back(path("back to lens", {sto, il, m, il}));
  s.paths.push_back(
      path("back image", {sto, il, m, il, {SurfaceId{"BIMG"}, EventKind::Transmit, 0}}));
  const CompiledSystem cs = compiled(s);
  const rtt::compile::CompiledSurface& lens = surface(cs, "IL");
  const rtt::trace::EventMedia media = environment(cs);

  const Vec3 object(3.0, -4.0, -75.0);
  RayBatch at_mirror = rays_from(object);
  RayBatch at_lens = rays_from(object);
  RayBatch at_image = rays_from(object);
  trace(cs, "to mirror", at_mirror);
  trace(cs, "back to lens", at_lens);
  trace(cs, "back image", at_image);
  all_alive(at_mirror);
  all_alive(at_lens);
  all_alive(at_image);
  const double m_back = sb_image / sb;
  const Vec3 image(-2.0 * m_back * object.x(), -2.0 * m_back * object.y(), -sb_image);
  for (std::size_t i = 0; i < at_mirror.size(); ++i) {
    INFO("ray " << i);
    // The mirrored forward passage: the ray after M with z -> -z, through IL in +z.
    rtt::trace::RayState mirrored;
    mirrored.pos = pos(at_mirror, i).cwiseProduct(Vec3(1.0, 1.0, -1.0));
    mirrored.dir = dir(at_mirror, i).cwiseProduct(Vec3(1.0, 1.0, -1.0));
    mirrored.opl = at_mirror.opl()[i];
    const auto hit = rtt::trace::intersect_surface(mirrored, lens);
    const auto forward =
        rtt::trace::apply_event(mirrored, lens, hit, 0, EventKind::Transmit, 0, media);
    REQUIRE(forward.status == RayStatus::Alive);
    CHECK(std::abs(at_lens.pos_x()[i] - forward.pos.x()) < 1e-12);
    CHECK(std::abs(at_lens.pos_y()[i] - forward.pos.y()) < 1e-12);
    CHECK(std::abs(at_lens.dir_x()[i] - forward.dir.x()) < 1e-15);
    CHECK(std::abs(at_lens.dir_y()[i] - forward.dir.y()) < 1e-15);
    CHECK(std::abs(at_lens.dir_z()[i] + forward.dir.z()) < 1e-15);
    CHECK(std::abs(at_lens.opl()[i] - forward.opl) < 1e-12);  // Delta of the backward passage
    // The backward image: every ray in m_back times the mirrored forward image.
    CHECK((pos(at_image, i) - image).norm() < kPositionTolerance);
  }
  // Content guard: the second passage deflects (without it the rays would meet at z = -100).
  CHECK(std::abs(at_lens.dir_y()[at_lens.size() - 1] - dir(at_mirror, at_mirror.size() - 1).y()) >
        1e-3);
}

TEST_CASE("ideal lens: refract and reflect at the lens are EventImpossible", "[ideal][trace]") {
  // ADR 0031, point 7: only transmit; compile rejects the other kinds at the path
  // (paths.ideal_lens_event), apply_event itself gives EventImpossible as before.
  const CompiledSystem cs = compiled(reference());
  for (const char* id : {"IL", "CL"}) {
    const rtt::compile::CompiledSurface& lens = surface(cs, id);
    rtt::trace::RayState ray;
    ray.pos = Vec3(0.0, 0.0, -75.0);
    ray.dir = (Vec3(2.0, 3.0, 0.0) - ray.pos).normalized();
    const auto hit = rtt::trace::intersect_surface(ray, lens);
    for (const auto kind : {rtt::model::EventKind::Refract, rtt::model::EventKind::Reflect}) {
      INFO(id << ", kind " << static_cast<int>(kind));
      CHECK(rtt::trace::apply_event(ray, lens, hit, 0, kind, 0, environment(cs)).status ==
            RayStatus::EventImpossible);
    }
  }
}

TEST_CASE("ideal lens: a ray in the lens plane is EventImpossible", "[ideal][trace]") {
  const CompiledSystem cs = compiled(reference());
  const rtt::compile::CompiledSurface& lens = surface(cs, "IL");
  rtt::trace::RayState ray;
  ray.pos = Vec3(-5.0, 0.0, 0.0);
  ray.dir = Vec3::UnitX();
  rtt::trace::SurfaceHit hit;  // a grazing hit as intersect_surface cannot give it
  hit.status = rtt::geom::HitStatus::Hit;
  hit.point = Vec3(0.0, 0.0, 0.0);
  hit.normal = Vec3::UnitZ();
  hit.direction = Vec3::UnitX();
  hit.t = 5.0;
  const auto out = rtt::trace::apply_event(ray, lens, hit, 0, rtt::model::EventKind::Transmit, 0,
                                           environment(cs));
  CHECK(out.status == RayStatus::EventImpossible);
}
