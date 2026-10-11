// Accumulated geometric transformation Q, P and Q per slot in RayPaths, and the ideal conductor
// flag (#207, ADR 0034, points 1, 3, 7 and 9; closes #76).
//
// Sources: W.-S. T. Lam, "Anisotropic Ray Trace", Eqs. (4.5)/(4.6) (geometric transformation,
// with the (3, 3) element 1; rtt::polar::geometric_transform), Sec. 4.5.1 (physical part Q^T P);
// S. J. Byrnes, arXiv:1603.02720v5, Eq. (6) (Fresnel amplitudes, rtt::polar::fresnel); the forms
// of Q without singularity are derived in ADR 0034, point 1 (docs/quellen.md).
//
// Tolerances: Q of one event is formed from unit vectors with a few dozen roundings of u = 2^-53
// per element, so two forms of the same matrix differ by a few 1e-16; 1e-14 is more than 10
// times that. Products of two events add as much again. The retardance compares phases from
// Fresnel amplitudes of a metal (|r| about 0.9): 1e-10 as in test_m3_acceptance.cpp.

#include <Eigen/Geometry>
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <numbers>
#include <optional>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/geom/phase.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
#include "rtt/polar/fresnel.hpp"
#include "rtt/polar/prt.hpp"
#include "rtt/polar/prt_analysis.hpp"
#include "rtt/trace/apply_event.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/ray_paths.hpp"
#include "rtt/trace/sequential.hpp"

using Cx = std::complex<double>;
using rtt::compile::compile;
using rtt::compile::CompiledSurface;
using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::material::MaterialLibrary;
using rtt::math::CMat3;
using rtt::math::Mat3;
using rtt::math::Vec3;
using rtt::model::Element;
using rtt::model::ElementKind;
using rtt::model::EventKind;
using rtt::model::Param;
using rtt::model::Pose;
using rtt::model::Surface;
using rtt::model::SurfaceId;
using rtt::model::System;
using rtt::trace::EventMedia;
using rtt::trace::RayBatch;
using rtt::trace::RayPaths;
using rtt::trace::RayState;
using rtt::trace::RayStatus;
using rtt::trace::SequentialTracer;

namespace {

constexpr double kPi = std::numbers::pi;
constexpr double kDeg = kPi / 180.0;
const Cx kMetal(1.2, 7.26);  // invented test metal of m3/metal_periscope (#56)

/// Plane surface at the global origin (local = global), normal +z, Fresnel interaction.
CompiledSurface plane_surface() {
  CompiledSurface s;
  s.shape = rtt::geom::Plane<double>{};
  return s;
}

/// Media of one event: real indices before, after (refraction) and beyond (reflection).
EventMedia media(double before, double after, double beyond, double wavelength_um = 0.55) {
  EventMedia m;
  m.before = before;
  m.after = after;
  m.beyond = beyond;
  m.wavelength_um = wavelength_um;
  return m;
}

/// Ray with unit direction `d` towards the global origin from z = -10 mm.
RayState ray_along(const Vec3& d) {
  RayState r;
  r.dir = d;
  r.pos = -10.0 / d.z() * d;
  return r;
}

RayState step(
    const RayState& ray, const CompiledSurface& s, EventKind kind, int order, const EventMedia& m) {
  return rtt::trace::sequential_step(ray, s, 0, kind, order, m);
}

double max_abs(const Mat3& m) {
  return m.cwiseAbs().maxCoeff();
}

System load(const std::string& file) {
  return rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/" + file);
}

PathId path_named(const CompiledSystem& cs, const std::string& name) {
  return *cs.find_path(name);
}

/// One ray from `start` along `dir` on `path`.
RayBatch trace_one(const CompiledSystem& cs, PathId path, const Vec3& start, const Vec3& dir) {
  RayBatch rays(1);
  rays.pos_x()[0] = start.x();
  rays.pos_y()[0] = start.y();
  rays.pos_z()[0] = start.z();
  rays.dir_x()[0] = dir.x();
  rays.dir_y()[0] = dir.y();
  rays.dir_z()[0] = dir.z();
  [[maybe_unused]] const auto stats = SequentialTracer().trace(cs, path, rays);
  return rays;
}

Vec3 direction(const RayBatch& r, std::size_t i) {
  return {r.dir_x()[i], r.dir_y()[i], r.dir_z()[i]};
}

/// Global normal of surface `id`: its local +z.
Vec3 normal_of(const CompiledSystem& cs, const std::string& id) {
  return cs.surfaces()[*cs.find_surface(SurfaceId(id))].to_global.apply_vector(Vec3::UnitZ());
}

Vec3 mirrored(const Vec3& k, const Vec3& n) {
  return k - 2.0 * k.dot(n) * n;
}

/// m3/metal_periscope with the second mirror folding +y into +x instead of +z: its plane of
/// incidence (x-y) is perpendicular to that of the first mirror (y-z), so s of the first
/// reflection is p of the second ("crossed fold"). R = Rx Ry Rz (docs/architecture.md): rotation
/// (90, 45, 0) deg turns the local z into (sin 45, -cos 45, 0); the detector at x = 20 faces +x.
System crossed_periscope() {
  System s = load("m3/metal_periscope.rtt.json");
  for (auto& child : s.root.children) {
    auto* e = std::get_if<Element>(&child.value);
    if (e == nullptr) continue;
    if (e->name == "mirror 2") {
      e->pose.rotation_deg = {Param(90.0), Param(45.0), Param(0.0)};
    } else if (e->name == "image") {
      e->pose.position = {Param(20.0), Param(20.0), Param(10.0)};
      e->pose.rotation_deg = {Param(0.0), Param(90.0), Param(0.0)};
    }
  }
  return s;
}

}  // namespace

// ------------------------------------------------------------ point 1: Q of one event -----

TEST_CASE("Q of one event equals Lam's geometric transformation (4.5)/(4.6)", "[polarization][q]") {
  // Refraction: Q = R(k_in -> k_out); reflection: Q = I - 2 N N^T (ADR 0034, point 1, derived
  // there); both against rtt::polar::geometric_transform with s = k_in x N, at many angles and
  // azimuths, outside and inside the glass (TIR as a Reflect event). 1e-14 (file comment).
  const CompiledSurface surface = plane_surface();
  const double critical = std::asin(1.0 / 1.5);
  std::size_t checked = 0;
  for (const double theta_deg : {1.0, 10.0, 30.0, 50.0, 70.0, 85.0}) {
    for (const double phi_deg : {0.0, 37.0, 90.0, 200.0}) {
      const double t = theta_deg * kDeg;
      const double p = phi_deg * kDeg;
      const Vec3 d(std::sin(t) * std::cos(p), std::sin(t) * std::sin(p), std::cos(t));
      struct Case {
        const char* name;
        EventKind kind;
        EventMedia m;
        bool possible;
      };
      const std::array<Case, 4> cases{{
          {"refract 1 -> 1.5", EventKind::Refract, media(1.0, 1.5, 1.5), true},
          {"refract 1.5 -> 1", EventKind::Refract, media(1.5, 1.0, 1.0), t < critical},
          {"reflect outside", EventKind::Reflect, media(1.0, 1.5, 1.5), true},
          {"reflect inside (TIR)", EventKind::Reflect, media(1.5, 1.0, 1.0), true},
      }};
      for (const Case& c : cases) {
        if (!c.possible) continue;
        INFO(c.name << ", theta " << theta_deg << " deg, phi " << phi_deg << " deg");
        const RayState in = ray_along(d);
        const RayState out = step(in, surface, c.kind, 0, c.m);
        REQUIRE(out.status == RayStatus::Alive);
        const Mat3 expected = rtt::polar::geometric_transform<double>(
            in.dir, out.dir, Vec3::UnitZ(), c.kind == EventKind::Reflect);
        CHECK(max_abs(out.q - expected) <= 1e-14);
        CHECK((out.q * in.dir - out.dir).norm() <= 1e-14);
        // Content guard: away from normal incidence Q is not the identity.
        if (theta_deg >= 10.0 || c.kind == EventKind::Reflect) {
          CHECK(max_abs(out.q - Mat3::Identity()) > 1e-3);
        }
        ++checked;
      }
    }
  }
  CHECK(checked == 6 * 4 * 3 + 2 * 4);  // refract 1.5 -> 1 only for theta < 41.8 deg
}

TEST_CASE("Q at normal incidence without a special case", "[polarization][q]") {
  // Lam's s = k x N is not defined there; the forms of ADR 0034 are: refraction R(k -> k) = I
  // (the refracted direction may differ from k by an ulp: 1e-15), reflection I - 2 N N^T =
  // diag(1, 1, -1) for N = z, bitwise.
  const CompiledSurface surface = plane_surface();
  const RayState in = ray_along(Vec3::UnitZ());
  const RayState refracted = step(in, surface, EventKind::Refract, 0, media(1.0, 1.5, 1.5));
  REQUIRE(refracted.status == RayStatus::Alive);
  CHECK(max_abs(refracted.q - Mat3::Identity()) <= 1e-15);
  const RayState reflected = step(in, surface, EventKind::Reflect, 0, media(1.0, 1.5, 1.5));
  REQUIRE(reflected.status == RayStatus::Alive);
  CHECK(reflected.q == Mat3(Eigen::Vector3d(1.0, 1.0, -1.0).asDiagonal()));
  CHECK((reflected.q * in.dir - reflected.dir).norm() <= 1e-15);
}

TEST_CASE("Q of a path is the product of the events: two mirrors", "[polarization][q]") {
  // The crossed periscope (two 45 deg metal mirrors with perpendicular planes of incidence):
  // Q = Q_2 Q_1 with Q_i = geometric_transform of each reflection (normals of the compiled
  // mirrors). 1e-14 (file comment). The transmit at the detector keeps Q (point 1: events
  // without deflection). Not the parallel periscope: there N_2 = -N_1, and the two Householder
  // reflections give Q = I.
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(crossed_periscope(), lib);
  const Vec3 k0 = Vec3::UnitZ();
  const RayBatch rays = trace_one(cs, path_named(cs, "periscope"), Vec3::Zero(), k0);
  REQUIRE(rays.status()[0] == RayStatus::Alive);
  const Vec3 n1 = normal_of(cs, "M1");
  const Vec3 n2 = normal_of(cs, "M2");
  const Vec3 k1 = mirrored(k0, n1);
  const Vec3 k2 = mirrored(k1, n2);
  const Mat3 expected = rtt::polar::geometric_transform<double>(k1, k2, n2, true) *
                        rtt::polar::geometric_transform<double>(k0, k1, n1, true);
  const Mat3 q = rays.q_matrix(0);
  CHECK(max_abs(q - expected) <= 1e-14);
  CHECK((q * k0 - direction(rays, 0)).norm() <= 1e-14);
  // Content guard: two reflections give a rotation (det +1) that is not I.
  CHECK(std::abs(q.determinant() - 1.0) <= 1e-14);
  CHECK(max_abs(q - Mat3::Identity()) > 0.5);
}

TEST_CASE("Q of a diffraction order: Q = R(k_0 -> k_m) Q_0, so Q^T P = Q_0^T P_0",
          "[polarization][q]") {
  // ADR 0034, point 1, with ADR 0025, point 6: the physical part of an order is that of order 0.
  // Linear grating of 300 lines/mm at the origin, refraction 1 -> 1.5 and reflection. 1e-14.
  CompiledSurface grating = plane_surface();
  grating.phase_functions = {rtt::geom::LinearGratingPhase<double>(300.0, 0.0)};
  const RayState in = ray_along(Vec3(0.2, 0.1, 1.0).normalized());
  for (const EventKind kind : {EventKind::Refract, EventKind::Reflect}) {
    INFO("kind " << static_cast<int>(kind));
    const RayState zero = step(in, grating, kind, 0, media(1.0, 1.5, 1.5));
    const RayState first = step(in, grating, kind, 1, media(1.0, 1.5, 1.5));
    REQUIRE(zero.status == RayStatus::Alive);
    REQUIRE(first.status == RayStatus::Alive);
    REQUIRE((first.dir - zero.dir).norm() > 1e-2);  // content guard: a real order
    const Mat3 r = rtt::trace::rotation_between(zero.dir, first.dir);
    CHECK(max_abs(first.q - r * zero.q) <= 1e-14);
    const CMat3 m_first = first.q.transpose().cast<Cx>() * first.prt;
    const CMat3 m_zero = zero.q.transpose().cast<Cx>() * zero.prt;
    CHECK((m_first - m_zero).cwiseAbs().maxCoeff() <= 1e-14);
  }
}

TEST_CASE("events without deflection keep Q bitwise", "[polarization][q]") {
  // ADR 0034, point 1: Q = I for polarizer, retarder and transmit (dummy passage). A ray that
  // starts with a Q != I (set in the batch) leaves with the same Q, bit for bit.
  System s;
  s.name = "dummy events";
  s.environment.medium = "VACUUM";
  s.wavelengths = {{0.55, 1.0, true}};
  s.aperture = {rtt::model::SystemApertureType::EntrancePupilDiameter, Param(10.0)};
  s.fields = {rtt::model::FieldType::AngleDeg, {{0.0, 0.0, 1.0}}};
  s.root.name = "root";
  const auto thin = [](const std::string& name, double z, rtt::model::Interaction interaction) {
    Surface surface;
    surface.id = SurfaceId(name);
    surface.interaction = std::move(interaction);
    return Element{name, ElementKind::ThinElement, Pose::along_z(z), std::nullopt, {surface}};
  };
  Surface plate_front;
  plate_front.id = SurfaceId("W.S1");
  Surface plate_back;
  plate_back.id = SurfaceId("W.S2");
  plate_back.pose = Pose::along_z(2.0);
  Surface image;
  image.id = SurfaceId("IMG");
  s.root.children = {
      {thin("POL", 0.0, rtt::model::IdealPolarizer{{1.0, 0.0, 0.0}, 0.0})},
      {thin("RET", 5.0, rtt::model::IdealRetarder{{1.0, 1.0, 0.0}, 0.25})},
      {Element{
          "W", ElementKind::Plate, Pose::along_z(10.0), "CONST:1.5", {plate_front, plate_back}}},
      {Element{"D", ElementKind::Detector, Pose::along_z(20.0), std::nullopt, {image}}}};
  s.paths = {{"dummy",
              false,
              {{SurfaceId("POL"), EventKind::Transmit, 0},
               {SurfaceId("RET"), EventKind::Transmit, 0},
               {SurfaceId("W.S1"), EventKind::Transmit, 0},
               {SurfaceId("W.S2"), EventKind::Transmit, 0},
               {SurfaceId("IMG"), EventKind::Transmit, 0}}}};
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(s, lib);
  RayBatch rays(1);
  rays.pos_z()[0] = -1.0;
  const Mat3 start = Eigen::AngleAxisd(0.3, Vec3::UnitZ()).toRotationMatrix();
  rays.set_q_matrix(0, start);
  [[maybe_unused]] const auto stats = SequentialTracer().trace(cs, PathId{0}, rays);
  REQUIRE(rays.status()[0] == RayStatus::Alive);
  CHECK(rays.q_matrix(0) == start);
  // Content guard: the polarizer acted on P.
  CHECK(std::abs(rays.weight()[0] - 0.5) <= 1e-12);
}

TEST_CASE("Q is NaN from the entry into a crystal on", "[polarization][q]") {
  // ADR 0034, point 1: Q of ordinary and extraordinary events is not defined in M6. Calcite plate
  // of m4/calcite_walkoff, both modes; P stays finite (content guard).
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(load("m4/calcite_walkoff.rtt.json"), lib);
  for (std::uint32_t p = 0; p < cs.paths().size(); ++p) {
    INFO("path " << cs.paths()[p].name);
    const RayBatch rays = trace_one(cs, PathId{p}, Vec3(0.0, 0.0, -1.0), Vec3::UnitZ());
    REQUIRE(rays.status()[0] == RayStatus::Alive);
    CHECK(rays.prt_matrix(0).allFinite());
    const Mat3 q = rays.q_matrix(0);
    for (int r = 0; r < 3; ++r) {
      for (int c = 0; c < 3; ++c) CHECK(std::isnan(q(r, c)));
    }
  }
}

// ------------------------------------------ point 4: physical retardance of a path -----

TEST_CASE("physical retardance of two metal mirrors: parallel folds add, crossed folds cancel",
          "[polarization][q]") {
  // ADR 0034, point 4: M = Q^T P maps k_0 onto itself; for one reflection M = diag(r_s, -r_p) in
  // (s, p_in). Two 45 deg reflections in the same plane (periscope): x gets r_s^2, the other
  // direction r_p^2, so delta = |arg((-r_p / r_s)^2)| reduced to [0, pi], the sum of the two
  // single retardances. Crossed folds: s of the first is p of the second, both directions get
  // -r_s r_p, delta = 0. Fresnel amplitudes from Byrnes, Eq. (6) (rtt::polar::fresnel), not from
  // the tracer. 1e-10 (file comment).
  const auto a = rtt::polar::fresnel_at_angle<double>(Cx(1.0), kMetal, std::cos(45.0 * kDeg));
  const Cx single = -a.rp / a.rs;
  const double sum = std::abs(std::arg(single * single));
  REQUIRE(sum > 0.1);  // content guard: a real retardance
  const MaterialLibrary lib;
  const Vec3 k0 = Vec3::UnitZ();
  {
    const CompiledSystem cs = compile(load("m3/metal_periscope.rtt.json"), lib);
    const RayBatch rays = trace_one(cs, path_named(cs, "periscope"), Vec3::Zero(), k0);
    REQUIRE(rays.status()[0] == RayStatus::Alive);
    REQUIRE((direction(rays, 0) - k0).norm() <= 1e-12);
    const auto r = rtt::polar::physical_retardance(rays.prt_matrix(0), rays.q_matrix(0), k0);
    CHECK(std::abs(r.value - sum) <= 1e-10);
  }
  {
    const CompiledSystem cs = compile(crossed_periscope(), lib);
    const RayBatch rays = trace_one(cs, path_named(cs, "periscope"), Vec3::Zero(), k0);
    REQUIRE(rays.status()[0] == RayStatus::Alive);
    REQUIRE((direction(rays, 0) - Vec3::UnitX()).norm() <= 1e-12);  // the crossed geometry
    const auto r = rtt::polar::physical_retardance(rays.prt_matrix(0), rays.q_matrix(0), k0);
    CHECK(r.value <= 1e-10);
  }
}

// -------------------------------------------------------- point 3: P and Q in RayPaths -----

TEST_CASE("RayPaths with polarization: P and Q per slot, the last valid slot is the RayBatch",
          "[polarization][paths]") {
  // ADR 0034, point 3: slot 0 is the input state, slot s the product after event s - 1; the last
  // valid slot (count - 1) is bitwise P and Q of the RayBatch, also for a ray lost on the way
  // (aperture of radius 3 mm on M2). Without the request prt and q stay empty and every other
  // array is bitwise the same.
  System s = load("m3/metal_periscope.rtt.json");
  for (auto& child : s.root.children) {
    auto* e = std::get_if<Element>(&child.value);
    if (e != nullptr && e->name == "mirror 2") {
      e->surfaces[0].aperture = rtt::model::CircularAperture{3.0, 0.0};
    }
  }
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(s, lib);
  const PathId path = path_named(cs, "periscope");
  RayBatch start(3);
  start.pos_x()[1] = 1.0;
  start.pos_x()[2] = 5.0;  // outside the aperture of M2: vignetted there
  // Ray 1 enters with P and Q != I: slot 0 must show them.
  const Mat3 q_in = Eigen::AngleAxisd(0.2, Vec3::UnitZ()).toRotationMatrix();
  start.set_q_matrix(1, q_in);
  start.set_prt_matrix(1, q_in.cast<Cx>() * Cx(0.9, 0.1));

  RayBatch with = start;
  RayPaths paths;
  paths.polarization = true;
  [[maybe_unused]] const auto s1 = SequentialTracer().trace(cs, path, with, paths);
  CHECK(paths.polarization);
  REQUIRE(paths.prt.size() == 3 * paths.slots * 9);
  REQUIRE(paths.q.size() == 3 * paths.slots * 9);
  REQUIRE(with.status()[0] == RayStatus::Alive);
  REQUIRE(with.status()[2] == RayStatus::Vignetted);  // content guard: a ray lost on the way
  for (std::size_t r = 0; r < 3; ++r) {
    INFO("ray " << r);
    CHECK(paths.prt_at(r, 0) == start.prt_matrix(r));
    CHECK(paths.q_at(r, 0) == start.q_matrix(r));
    const std::size_t last = paths.count[r] - 1;
    CHECK(paths.prt_at(r, last) == with.prt_matrix(r));
    CHECK(paths.q_at(r, last) == with.q_matrix(r));
    // Slot 1 is the product after the first event: Q_1 times the input (content guard).
    CHECK(max_abs(paths.q_at(r, 1) - paths.q_at(r, 0)) > 1e-3);
  }
  CHECK(std::isnan(paths.q_at(2, paths.slots - 1)(0, 0)));  // unused slot of the lost ray

  RayBatch without = start;
  RayPaths plain;
  [[maybe_unused]] const auto s2 = SequentialTracer().trace(cs, path, without, plain);
  CHECK_FALSE(plain.polarization);
  CHECK(plain.prt.empty());
  CHECK(plain.q.empty());
  CHECK(plain.count == paths.count);
  CHECK(plain.lost_at == paths.lost_at);
  CHECK(plain.status == paths.status);
  const auto same_bits = [](const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i) {
      if (std::isnan(a[i]) != std::isnan(b[i]) || (!std::isnan(a[i]) && a[i] != b[i])) {
        return false;
      }
    }
    return true;
  };
  CHECK(same_bits(plain.position, paths.position));
  CHECK(same_bits(plain.direction, paths.direction));
  CHECK(same_bits(plain.opl, paths.opl));
  CHECK(same_bits(plain.weight, paths.weight));
  for (std::size_t r = 0; r < 3; ++r) {
    CHECK(without.prt_matrix(r) == with.prt_matrix(r));
    CHECK(without.q_matrix(r) == with.q_matrix(r));
  }
}

// ------------------------------------------------------- point 7: the ideal conductor -----

TEST_CASE("a mirror with a material of the surrounding index is no ideal conductor",
          "[polarization][conductor]") {
  // ADR 0034, point 7: only a mirror without material is an ideal conductor. A mirror of
  // CONST:1.0 in VACUUM reflects with Fresnel against n = 1: r_s = r_p = 0, weight 0 (before
  // #207 it counted as an ideal conductor, weight 1). The mirror without material keeps r_s = -1,
  // r_p = +1, weight 1.
  const MaterialLibrary lib;
  for (const bool with_material : {true, false}) {
    INFO("material " << with_material);
    System s = load("m3/metal_periscope.rtt.json");
    for (auto& child : s.root.children) {
      auto* e = std::get_if<Element>(&child.value);
      if (e != nullptr && e->kind == ElementKind::Mirror) {
        e->material = with_material ? std::optional<std::string>("CONST:1.0") : std::nullopt;
      }
    }
    const CompiledSystem cs = compile(s, lib);
    const RayBatch rays = trace_one(cs, path_named(cs, "single"), Vec3::Zero(), Vec3::UnitZ());
    REQUIRE(rays.status()[0] == RayStatus::Alive);
    CHECK(std::abs(rays.weight()[0] - (with_material ? 0.0 : 1.0)) <= 1e-15);
  }
}
