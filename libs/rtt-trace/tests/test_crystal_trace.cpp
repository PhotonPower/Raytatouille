// Uniaxial crystals in the tracer (#132, ADR 0026, points 3-5): the ray follows S, the phase
// follows k (OPL = n l (k . S), Lam, Eq. (2.17)), entry with Ordinary/Extraordinary and exit with
// Refract through rtt/polar/birefringence.hpp (#130), modes from compile (#131), orders at the
// exit (ADR 0025, with #127). Expected values are closed forms, computed independently.

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
#include "rtt/polar/birefringence.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/sequential.hpp"

using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::material::MaterialLibrary;
using rtt::math::Complex;
using rtt::math::Vec3;
using rtt::model::CrystalMaterial;
using rtt::model::Element;
using rtt::model::ElementKind;
using rtt::model::EventKind;
using rtt::model::Param;
using rtt::model::Pose;
using rtt::model::Surface;
using rtt::model::SurfaceId;
using rtt::model::System;
using rtt::trace::RayBatch;
using rtt::trace::RayStatus;
using rtt::trace::SequentialTracer;

namespace {

constexpr double kNo = 1.6584;  // calcite, CONST indices of the ADR 0026 acceptance
constexpr double kNe = 1.4864;
constexpr double kThickness = 2.0;  // mm, plate from z = 10 to z = 12, detector at z = 20

/// Lam, Eq. (2.39) at 45 degree.
double n_e45() {
  return std::sqrt(2.0 / (1.0 / (kNo * kNo) + 1.0 / (kNe * kNe)));
}

/// Walk-off at 45 degree from S normal to the K-surface (Lam, p. 107; derivation in
/// docs/quellen.md): tan(rho) = (nO^2 - nE^2) / (nO^2 + nE^2).
double rho45() {
  return std::atan((kNo * kNo - kNe * kNe) / (kNo * kNo + kNe * kNe));
}

Surface plane(const std::string& id, double z_mm = 0.0) {
  Surface s;
  s.id = SurfaceId(id);
  s.pose = Pose::along_z(z_mm);
  return s;
}

/// Calcite plate P (z = 10 ... 12) with the given optic axis, detector at z = 20, in vacuum.
System crystal_system(const std::array<double, 3>& axis) {
  System s;
  s.name = "crystal";
  s.environment.medium = "VACUUM";
  s.wavelengths = {{0.5876, 1.0, true}};
  s.aperture = {rtt::model::SystemApertureType::EntrancePupilDiameter, Param(2.0)};
  s.fields = {rtt::model::FieldType::AngleDeg, {{0.0, 0.0, 1.0}}};
  s.root.name = "root";
  Element plate{"P",
                ElementKind::Plate,
                Pose::along_z(10.0),
                std::nullopt,
                {plane("P.S1"), plane("P.S2", kThickness)}};
  plate.crystal = CrystalMaterial{"CONST:1.6584", "CONST:1.4864"};
  plate.optic_axis = axis;
  s.root.children.push_back({plate});
  s.root.children.push_back(
      {Element{"image", ElementKind::Detector, Pose::along_z(20.0), std::nullopt, {plane("IMG")}}});
  return s;
}

/// Path through the plate in `mode` and to the detector.
rtt::model::Path through(EventKind mode, int exit_order = 0) {
  return {"through",
          false,
          {{SurfaceId("P.S1"), mode, 0},
           {SurfaceId("P.S2"), EventKind::Refract, exit_order},
           {SurfaceId("IMG"), EventKind::Transmit, 0}}};
}

/// One ray per entry, set by hand (position and direction only, ADR 0026 reading rule).
RayBatch rays_at(const std::vector<std::pair<Vec3, Vec3>>& start) {
  RayBatch rays(start.size());
  for (std::size_t i = 0; i < start.size(); ++i) {
    rays.pos_x()[i] = start[i].first.x();
    rays.pos_y()[i] = start[i].first.y();
    rays.pos_z()[i] = start[i].first.z();
    rays.dir_x()[i] = start[i].second.x();
    rays.dir_y()[i] = start[i].second.y();
    rays.dir_z()[i] = start[i].second.z();
  }
  return rays;
}

Vec3 pos(const RayBatch& r, std::size_t i) {
  return {r.pos_x()[i], r.pos_y()[i], r.pos_z()[i]};
}
Vec3 dir(const RayBatch& r, std::size_t i) {
  return {r.dir_x()[i], r.dir_y()[i], r.dir_z()[i]};
}
Vec3 wave(const RayBatch& r, std::size_t i) {
  return {r.wave_x()[i], r.wave_y()[i], r.wave_z()[i]};
}

RayBatch trace(const CompiledSystem& cs, std::uint32_t path, RayBatch rays) {
  [[maybe_unused]] const auto stats = SequentialTracer().trace(cs, PathId{path}, rays);
  return rays;
}

}  // namespace

TEST_CASE("M4 acceptance: calcite walk-off from tests/reference/m4/calcite_walkoff.rtt.json",
          "[crystal][acceptance]") {
  // Plate 2 mm, axis (1, 0, 1) at 45 degree in the x-z plane (ADR 0026, Abnahme M4), normal
  // incidence in vacuum, detector at z = 20. Expected (independently computed):
  // n_e(45) = 1.565356826060665,
  // rho = 6.224117602966 degree, offset t tan(rho) = 0.218121366133243 mm away from the axis,
  // OPL_o = 18 + nO t = 21.3168, OPL_e = 18 + n_e(45) t = 21.130713652121329 (in the plate
  // l = t / cos(rho) along S and k . S = cos(rho)).
  const System s =
      rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m4/calcite_walkoff.rtt.json");
  const MaterialLibrary lib;
  const CompiledSystem cs = rtt::compile::compile(s, lib);
  const Vec3 start(0.0, 0.3, 0.0);
  for (const char* name : {"o", "e"}) {
    INFO(name);
    const auto path = cs.find_path(name);
    REQUIRE(path.has_value());
    const RayBatch r = trace(cs, path->index, rays_at({{start, Vec3::UnitZ()}}));
    REQUIRE(r.status()[0] == RayStatus::Alive);
    REQUIRE((dir(r, 0) - Vec3::UnitZ()).norm() <= 1e-12);  // exit parallel to the incidence
    const bool e = std::string(name) == "e";
    const double offset = e ? -kThickness * std::tan(rho45()) : 0.0;  // away from a_x > 0
    REQUIRE((pos(r, 0) - Vec3(offset, 0.3, 20.0)).norm() <= 1e-10);
    REQUIRE(std::abs(r.opl()[0] - (18.0 + (e ? n_e45() : kNo) * kThickness)) <= 1e-10);
    REQUIRE(r.mode_index()[0] == 0.0);  // back in vacuum: no mode, wave = dir
    REQUIRE((wave(r, 0) - dir(r, 0)).norm() == 0.0);
  }
  REQUIRE(std::abs(kThickness * std::tan(rho45()) - 0.218121366133243) <= 1e-12);
  REQUIRE(std::abs(n_e45() - 1.565356826060665) <= 1e-12);
}

TEST_CASE("crystal: no walk-off with the axis along or across the normal", "[crystal]") {
  // Axis along N: k along the axis, n = n_O; axis in the surface: k across the axis, n = n_E
  // (Lam, Eq. (2.39)); in both cases S = k.
  const MaterialLibrary lib;
  for (const auto& [axis, n] :
       {std::pair{std::array{0.0, 0.0, 1.0}, kNo}, std::pair{std::array{1.0, 0.0, 0.0}, kNe}}) {
    System s = crystal_system(axis);
    s.paths = {through(EventKind::Extraordinary)};
    const CompiledSystem cs = rtt::compile::compile(s, lib);
    const RayBatch r = trace(cs, 0, rays_at({{Vec3(0.0, 0.3, 0.0), Vec3::UnitZ()}}));
    REQUIRE(r.status()[0] == RayStatus::Alive);
    REQUIRE((pos(r, 0) - Vec3(0.0, 0.3, 20.0)).norm() <= 1e-12);
    REQUIRE(std::abs(r.opl()[0] - (18.0 + n * kThickness)) <= 1e-12);
  }
}

TEST_CASE("crystal: oblique ordinary ray follows Snell with n_O", "[crystal]") {
  // As the isotropic plate test (test_sequential.cpp): offset d sin(a - a') / cos(a') across
  // the ray, exit parallel to the incidence, OPL = vacuum path + n_O * glass path.
  const MaterialLibrary lib;
  System s = crystal_system({0.2, -0.5, 0.8});
  s.paths = {through(EventKind::Ordinary)};
  const CompiledSystem cs = rtt::compile::compile(s, lib);
  const double a = 0.4;
  const Vec3 d0(0.0, std::sin(a), std::cos(a));
  const RayBatch r = trace(cs, 0, rays_at({{Vec3(0.0, 0.3, 0.0), d0}}));
  REQUIRE(r.status()[0] == RayStatus::Alive);
  REQUIRE((dir(r, 0) - d0).norm() <= 1e-12);
  const double at = std::asin(std::sin(a) / kNo);
  const double offset = kThickness * std::sin(a - at) / std::cos(at);
  const Vec3 straight = Vec3(0.0, 0.3, 0.0) + (20.0 / d0.z()) * d0;
  REQUIRE(std::abs(((pos(r, 0) - straight).cross(d0)).norm() - offset) <= 1e-12);
  REQUIRE(pos(r, 0).y() < straight.y());  // bent towards the normal inside the plate
  const double glass = kThickness / std::cos(at);
  const double vacuum = (20.0 - kThickness) / std::cos(a);
  REQUIRE(std::abs(r.opl()[0] - (vacuum + kNo * glass)) <= 1e-12);
}

TEST_CASE("crystal: a hand-filled batch with only dir set enters with the vector solution",
          "[crystal]") {
  // ADR 0026, point 3, reading rule: wave and mode_index count only for mode_index > 0; trace
  // sets wave := dir at the start. An oblique ray from a batch that sets only dir must refract
  // as rtt::polar::uniaxial_mode() says (here: ending inside the crystal at P.S2).
  const MaterialLibrary lib;
  const Vec3 axis = Vec3(0.2, -0.5, 0.8).normalized();
  System s = crystal_system({axis.x(), axis.y(), axis.z()});
  s.paths = {{"inside",
              false,
              {{SurfaceId("P.S1"), EventKind::Extraordinary, 0},
               {SurfaceId("P.S2"), EventKind::Transmit, 0}}}};
  const CompiledSystem cs = rtt::compile::compile(s, lib);
  const Vec3 d0 = Vec3(0.1, 0.45, std::sqrt(1.0 - 0.01 - 0.2025));
  const Vec3 start(0.0, 0.0, 0.0);
  const RayBatch r = trace(cs, 0, rays_at({{start, d0}}));
  REQUIRE(r.status()[0] == RayStatus::Alive);
  const auto m =
      rtt::polar::uniaxial_mode(rtt::polar::Mode::Extraordinary, Vec3(d0.x(), d0.y(), 0.0),
                                Vec3(Vec3::UnitZ()), kNo, kNe, axis);
  REQUIRE(m.has_value());
  REQUIRE((dir(r, 0) - m->s).norm() <= 1e-12);
  REQUIRE((wave(r, 0) - m->k).norm() <= 1e-12);
  REQUIRE(std::abs(r.mode_index()[0] - m->n) <= 1e-12);
  const Vec3 entry = start + (10.0 / d0.z()) * d0;
  const Vec3 end = entry + (kThickness / m->s.z()) * m->s;
  REQUIRE((pos(r, 0) - end).norm() <= 1e-12);
  // OPL (Lam, Eq. (2.17)): vacuum to the entry plus n l (k . S) in the crystal.
  const double l = (end - entry).norm();
  REQUIRE(std::abs(r.opl()[0] - ((entry - start).norm() + m->n * l * m->k.dot(m->s))) <= 1e-12);
}

TEST_CASE("crystal: P and weight of the projection model", "[crystal]") {
  // Normal incidence, axis at 45 degree in the x-z plane: y couples only to the o-mode, x only
  // to the e-mode (ADR 0026, point 5 and Abnahme M4); each mode carries half of unpolarized
  // light, and P S0 = S at the end (ADR 0021).
  const MaterialLibrary lib;
  System s = crystal_system({1.0, 0.0, 1.0});
  s.paths = {through(EventKind::Ordinary), through(EventKind::Extraordinary)};
  s.paths[1].name = "through e";
  const CompiledSystem cs = rtt::compile::compile(s, lib);
  const Vec3 x = Vec3::UnitX();
  const Vec3 y = Vec3::UnitY();
  for (std::uint32_t p : {0U, 1U}) {
    INFO("path " << p);
    const RayBatch r = trace(cs, p, rays_at({{Vec3(0.0, 0.3, 0.0), Vec3::UnitZ()}}));
    REQUIRE(r.status()[0] == RayStatus::Alive);
    REQUIRE(std::abs(r.weight()[0] - 0.5) <= 1e-12);
    const auto prt = r.prt_matrix(0);
    REQUIRE((prt * Vec3::UnitZ().cast<Complex>() - dir(r, 0).cast<Complex>()).norm() <= 1e-12);
    const Vec3& passes = p == 0 ? y : x;
    const Vec3& blocked = p == 0 ? x : y;
    REQUIRE(std::abs((prt * passes.cast<Complex>()).norm() - 1.0) <= 1e-12);
    REQUIRE((prt * blocked.cast<Complex>()).norm() <= 1e-12);
  }
}

TEST_CASE("crystal: unpolarized power of both modes sums to 1, also oblique", "[crystal]") {
  // Each mode carries s ||P_T||^2 / 2 = 1/2 of unpolarized light at the entry (P_T = E_v e_v^T,
  // ||P_T||_F = 1): the sum over both modes is exactly 1 (trace of u u^T + v v^T = 2), while
  // polarized powers may lie between 1 - |u . v| and 1 + |u . v| (ADR 0026, point 5; #130).
  const MaterialLibrary lib;
  System s = crystal_system({1.0, 0.0, 0.3});
  s.paths = {{"o", false, {{SurfaceId("P.S1"), EventKind::Ordinary, 0}}},
             {"e", false, {{SurfaceId("P.S1"), EventKind::Extraordinary, 0}}}};
  const CompiledSystem cs = rtt::compile::compile(s, lib);
  const double a = 1.2;
  const Vec3 d0(0.6 * std::sin(a), 0.8 * std::sin(a), std::cos(a));
  double sum = 0.0;
  for (std::uint32_t p : {0U, 1U}) {
    const RayBatch r = trace(cs, p, rays_at({{Vec3(0.0, 0.0, 9.0), d0}}));
    REQUIRE(r.status()[0] == RayStatus::Alive);
    sum += r.weight()[0];
  }
  REQUIRE(std::abs(sum - 1.0) <= 1e-12);
}

TEST_CASE("crystal: exit with a diffraction order (ADR 0025 with ADR 0026)", "[crystal]") {
  // Grating of G lines/mm (psi = 0) on the exit face P.S2. At normal incidence the wave normal
  // in the crystal is z for both modes (axis at 45 degree), so the tangential momentum at the
  // exit is n k_par + m lambda0 G x = m lambda0 G x, the same for o and e: the exit direction is
  // (lambda0 G, 0, sqrt(1 - (lambda0 G)^2)) in vacuum. lambda0 = 0.5876e-3 mm. The ray hits
  // at x = 0.4 (the walk-off runs along y), where the order adds m lambda0 G x to the OPL
  // (ADR 0025, point 3: m phi lambda0 / (2 pi) with phi = 2 pi G x).
  const MaterialLibrary lib;
  for (const double lines : {300.0, 2000.0}) {
    System s = crystal_system({0.0, 1.0, 1.0});
    auto& plate = std::get<Element>(s.root.children[0].value);
    plate.surfaces[1].phases.emplace_back(rtt::model::LinearGrating{Param(lines), 0.0});
    s.paths = {through(EventKind::Ordinary, 1), through(EventKind::Extraordinary, 1)};
    s.paths[1].name = "through e";
    const CompiledSystem cs = rtt::compile::compile(s, lib);
    const double sx = 0.5876e-3 * lines;
    for (std::uint32_t p : {0U, 1U}) {
      INFO("lines/mm " << lines << ", path " << p);
      const RayBatch r = trace(cs, p, rays_at({{Vec3(0.4, 0.3, 0.0), Vec3::UnitZ()}}));
      if (sx < 1.0) {
        REQUIRE(r.status()[0] == RayStatus::Alive);
        const double cz = std::sqrt(1.0 - sx * sx);
        REQUIRE((dir(r, 0) - Vec3(sx, 0.0, cz)).norm() <= 1e-12);
        const double n = p == 0 ? kNo : n_e45();
        REQUIRE(std::abs(r.opl()[0] - (10.0 + n * kThickness + sx * 0.4 + 8.0 / cz)) <= 1e-12);
      } else {
        // lambda0 G = 1.1752 > 1: order 0 exists, order 1 not (ADR 0025, point 7).
        REQUIRE(r.status()[0] == RayStatus::Evanescent);
        REQUIRE(r.last_surface()[0] == 1);
      }
    }
  }
}

TEST_CASE("crystal: no real mode at the entry gives Tir", "[crystal]") {
  // Environment n = 2: at 1.2 rad the tangential momentum 2 sin(1.2) = 1.864 exceeds n_O and
  // n_e(theta) <= n_O, so neither mode propagates (ADR 0026, point 4: Tir for order 0).
  const MaterialLibrary lib;
  System s = crystal_system({1.0, 0.0, 1.0});
  s.environment.medium = "CONST:2.0";
  s.paths = {through(EventKind::Ordinary), through(EventKind::Extraordinary)};
  s.paths[1].name = "through e";
  const CompiledSystem cs = rtt::compile::compile(s, lib);
  const Vec3 d0(0.0, std::sin(1.2), std::cos(1.2));
  for (std::uint32_t p : {0U, 1U}) {
    INFO("path " << p);
    const RayBatch r = trace(cs, p, rays_at({{Vec3(0.0, 0.0, 9.0), d0}}));
    REQUIRE(r.status()[0] == RayStatus::Tir);
    REQUIRE(r.last_surface()[0] == 0);
  }
}

TEST_CASE("isotropic systems keep wave = dir and mode_index = 0", "[crystal]") {
  // The OPL and every other value of isotropic systems stays bitwise as before #132 (the
  // existing tests); the new columns stay in their isotropic state.
  const MaterialLibrary lib;
  System s = crystal_system({0.0, 0.0, 1.0});
  auto& plate = std::get<Element>(s.root.children[0].value);
  plate.crystal.reset();
  plate.optic_axis.reset();
  plate.material = "CONST:1.5";
  s.paths = {{"main", true, {}}};
  const CompiledSystem cs = rtt::compile::compile(s, lib);
  const RayBatch r =
      trace(cs, 0,
            rays_at({{Vec3(0.0, 0.3, 0.0), Vec3::UnitZ()},
                     {Vec3(0.0, 0.0, 0.0), Vec3(Vec3(0.0, 0.3, 1.0).normalized())}}));
  for (std::size_t i = 0; i < r.size(); ++i) {
    REQUIRE(r.status()[i] == RayStatus::Alive);
    REQUIRE(r.mode_index()[i] == 0.0);
    REQUIRE((wave(r, i) - dir(r, i)).norm() == 0.0);
  }
}
