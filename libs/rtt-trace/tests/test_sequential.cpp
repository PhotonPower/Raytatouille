#include <oneapi/tbb/task_arena.h>

#include <Eigen/Geometry>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/sequential.hpp"

using rtt::compile::compile;
using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::material::MaterialLibrary;
using rtt::math::Vec3;
using rtt::model::Element;
using rtt::model::ElementKind;
using rtt::model::Param;
using rtt::model::Pose;
using rtt::model::Surface;
using rtt::model::SurfaceId;
using rtt::model::System;
using rtt::trace::RayBatch;
using rtt::trace::RayStatus;
using rtt::trace::SequentialTracer;

namespace {

/// Valid system with one wavelength, one field, EPD aperture, empty root and an automatic path.
System bare_system() {
  System s;
  s.name = "test";
  s.wavelengths = {{0.5876, 1.0, true}};
  s.aperture = {rtt::model::SystemApertureType::EntrancePupilDiameter, Param(10.0)};
  s.fields = {rtt::model::FieldType::AngleDeg, {{0.0, 0.0, 1.0}}};
  s.root.name = "root";
  s.paths = {{"main", true, {}}};
  return s;
}

Surface surface(const std::string& id, double z_mm = 0.0) {
  Surface s;
  s.id = SurfaceId(id);
  s.pose = Pose::along_z(z_mm);
  return s;
}

void set_ray(RayBatch& rays, std::size_t i, const Vec3& pos, const Vec3& dir) {
  rays.pos_x()[i] = pos.x();
  rays.pos_y()[i] = pos.y();
  rays.pos_z()[i] = pos.z();
  rays.dir_x()[i] = dir.x();
  rays.dir_y()[i] = dir.y();
  rays.dir_z()[i] = dir.z();
}

Vec3 pos(const RayBatch& rays, std::size_t i) {
  return {rays.pos_x()[i], rays.pos_y()[i], rays.pos_z()[i]};
}

Vec3 dir(const RayBatch& rays, std::size_t i) {
  return {rays.dir_x()[i], rays.dir_y()[i], rays.dir_z()[i]};
}

}  // namespace

TEST_CASE("paraboloid mirror focuses axial rays at R/2", "[sequential]") {
  // Concave paraboloid (k = -1) with R = -200 mm: centre of curvature on the -z side, focus at
  // z = R/2 = -100 mm in front of the mirror (focal property of the paraboloid; reference case
  // of docs/architecture.md, Validierung). Issue #6: RMS spot < 1e-9 mm.
  const System s =
      rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m1/paraboloid_mirror.rtt.json");
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(s, lib);

  std::vector<Vec3> starts;
  for (int i = -6; i <= 6; ++i) {
    for (int j = -6; j <= 6; ++j) {
      if (i * i + j * j <= 36) {
        starts.emplace_back(5.0 * i, 5.0 * j, -150.0);  // up to r = 30 mm (f/1.7)
      }
    }
  }
  RayBatch rays(starts.size());
  for (std::size_t i = 0; i < starts.size(); ++i) {
    set_ray(rays, i, starts[i], Vec3::UnitZ());
  }
  const auto stats = SequentialTracer().trace(cs, PathId{0}, rays);
  REQUIRE(stats.count(RayStatus::Alive) == starts.size());

  double sum_r2 = 0.0;
  for (std::size_t i = 0; i < rays.size(); ++i) {
    const Vec3 p = pos(rays, i);
    REQUIRE(std::abs(p.z() + 100.0) <= 1e-12);
    sum_r2 += p.x() * p.x() + p.y() * p.y();
    REQUIRE(rays.last_surface()[i] == 1);
  }
  const double rms = std::sqrt(sum_r2 / static_cast<double>(rays.size()));
  REQUIRE(rms < 1e-9);
  // Equal optical paths to the focus (docs/architecture.md: OPD < 1e-6 waves): from the start
  // plane z = -150 to the focus every ray travels 250 mm (directrix of the parabola at z = +100,
  // so 150 + 100 from start plane to directrix equals the path via the mirror to the focus).
  // 1e-6 waves at 0.5876 um is 5.9e-10 mm.
  for (std::size_t i = 0; i < rays.size(); ++i) {
    REQUIRE(std::abs(rays.opl()[i] - 250.0) <= 5e-10);
  }
}

TEST_CASE("plane-parallel plate shifts the ray by d sin(theta - theta') / cos(theta')",
          "[sequential]") {
  // Inside the plate the ray travels d / cos(theta') along a direction rotated by
  // theta - theta' against the incident one; the perpendicular offset after the plate is
  // therefore d sin(theta - theta') / cos(theta'), and the exit direction equals the incident
  // one (Snell's law: de Greve, Reflections and Refractions in Ray Tracing, 2006, Eq. (23)).
  // Issue #6: 1e-12 mm.
  const double thickness = 5.0;
  const double n = 1.5;
  System s = bare_system();
  s.root.children.push_back({Element{"P",
                                     ElementKind::Plate,
                                     Pose::along_z(10.0),
                                     "CONST:1.5",
                                     {surface("P.S1"), surface("P.S2", thickness)}}});
  s.root.children.push_back(
      {Element{"D", ElementKind::Detector, Pose::along_z(40.0), std::nullopt, {surface("D")}}});
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(s, lib);

  const std::vector<double> angles{0.0, 0.1, 0.3, 0.6, 1.0};
  RayBatch rays(angles.size());
  for (std::size_t i = 0; i < angles.size(); ++i) {
    set_ray(rays, i, Vec3(0.0, 1.0, 0.0), Vec3(0.0, std::sin(angles[i]), std::cos(angles[i])));
  }
  const RayBatch before = rays;
  [[maybe_unused]] const auto stats = SequentialTracer().trace(cs, PathId{0}, rays);
  for (std::size_t i = 0; i < angles.size(); ++i) {
    INFO("theta = " << angles[i] << " rad");
    REQUIRE(rays.status()[i] == RayStatus::Alive);
    const Vec3 d0 = dir(before, i);
    REQUIRE((dir(rays, i) - d0).cwiseAbs().maxCoeff() <= 1e-15);
    const double theta = angles[i];
    const double theta_t = std::asin(std::sin(theta) / n);
    const double expected = thickness * std::sin(theta - theta_t) / std::cos(theta_t);
    const double offset = (pos(rays, i) - pos(before, i)).cross(d0).norm();
    REQUIRE(std::abs(offset - expected) <= 1e-12);
    // OPL: air before and after the plate, glass inside.
    const double glass = thickness / std::cos(theta_t);
    const double air = (40.0 - thickness) / std::cos(theta);
    REQUIRE(std::abs(rays.opl()[i] - (air + n * glass)) <= 1e-12);
  }
}

TEST_CASE("trace is bitwise identical with one thread and with all threads", "[sequential]") {
  // ADR 0004: parallelisation must not change the result.
  const System s =
      rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m1/singlet_const.rtt.json");
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(s, lib);

  constexpr int kGrid = 61;
  RayBatch rays(static_cast<std::size_t>(kGrid) * kGrid * 3);
  std::size_t k = 0;
  for (int f = 0; f < 3; ++f) {
    for (int i = 0; i < kGrid; ++i) {
      for (int j = 0; j < kGrid; ++j, ++k) {
        const double x = -12.0 + 24.0 * i / (kGrid - 1);
        const double y = -12.0 + 24.0 * j / (kGrid - 1);
        set_ray(rays, k, Vec3(x, y, -5.0), Vec3(0.0, 0.05 * f, 1.0).normalized());
        rays.wl()[k] = static_cast<std::uint16_t>(f);
      }
    }
  }
  RayBatch serial = rays;
  RayBatch parallel = rays;
  oneapi::tbb::task_arena one(1);
  const auto serial_stats =
      one.execute([&] { return SequentialTracer().trace(cs, PathId{0}, serial); });
  const auto parallel_stats = SequentialTracer().trace(cs, PathId{0}, parallel);

  REQUIRE(serial_stats.rays == parallel_stats.rays);
  REQUIRE(serial_stats.count(RayStatus::Vignetted) > 0);  // the grid overfills the stop
  REQUIRE(serial_stats.count(RayStatus::Alive) > 0);
  for (std::size_t i = 0; i < rays.size(); ++i) {
    REQUIRE(serial.pos_x()[i] == parallel.pos_x()[i]);
    REQUIRE(serial.pos_y()[i] == parallel.pos_y()[i]);
    REQUIRE(serial.pos_z()[i] == parallel.pos_z()[i]);
    REQUIRE(serial.dir_x()[i] == parallel.dir_x()[i]);
    REQUIRE(serial.dir_y()[i] == parallel.dir_y()[i]);
    REQUIRE(serial.dir_z()[i] == parallel.dir_z()[i]);
    REQUIRE(serial.opl()[i] == parallel.opl()[i]);
    REQUIRE(serial.status()[i] == parallel.status()[i]);
    REQUIRE(serial.last_surface()[i] == parallel.last_surface()[i]);
  }
}

TEST_CASE("trace stats count the rays per status", "[sequential]") {
  System s = bare_system();
  Surface stop = surface("STO");
  stop.aperture = rtt::model::CircularAperture{5.0, 0.0};
  s.root.children.push_back({Element{"S", ElementKind::Stop, {}, std::nullopt, {stop}}});
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(s, lib);

  RayBatch rays(4);
  set_ray(rays, 0, Vec3(0.0, 1.0, -1.0), Vec3::UnitZ());   // passes
  set_ray(rays, 1, Vec3(0.0, 6.0, -1.0), Vec3::UnitZ());   // vignetted
  set_ray(rays, 2, Vec3(0.0, 1.0, -1.0), -Vec3::UnitZ());  // plane behind: missed
  set_ray(rays, 3, Vec3(0.0, 1.0, -1.0), Vec3::UnitZ());
  rays.status()[3] = RayStatus::Absorbed;  // already stopped: skipped
  const auto stats = SequentialTracer().trace(cs, PathId{0}, rays);
  REQUIRE(stats.count(RayStatus::Alive) == 1);
  REQUIRE(stats.count(RayStatus::Vignetted) == 1);
  REQUIRE(stats.count(RayStatus::Missed) == 1);
  REQUIRE(stats.count(RayStatus::Absorbed) == 1);
  REQUIRE(rays.last_surface()[3] == rtt::trace::kNoSurface);
}

TEST_CASE("invalid trace input throws before tracing", "[sequential]") {
  System s = bare_system();
  s.root.children.push_back(
      {Element{"D", ElementKind::Detector, {}, std::nullopt, {surface("D")}}});
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(s, lib);
  RayBatch rays(1);
  REQUIRE_THROWS_AS(SequentialTracer().trace(cs, PathId{5}, rays), std::out_of_range);
  rays.wl()[0] = 3;  // only one system wavelength
  REQUIRE_THROWS_AS(SequentialTracer().trace(cs, PathId{0}, rays), std::invalid_argument);
  rays.wl()[0] = 0;
  rays.status()[0] = static_cast<RayStatus>(9);  // e.g. garbage from a foreign caller
  REQUIRE_THROWS_AS(SequentialTracer().trace(cs, PathId{0}, rays), std::invalid_argument);
}
