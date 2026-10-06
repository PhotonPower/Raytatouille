// Recording ray paths during the sequential trace (#80, rtt/trace/ray_paths.hpp).

#include <oneapi/tbb/task_arena.h>

#include <bit>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/ray_paths.hpp"
#include "rtt/trace/sequential.hpp"
#include "rtt/trace/sources.hpp"

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
using rtt::trace::RayPaths;
using rtt::trace::RayStatus;
using rtt::trace::SequentialTracer;

namespace {

System singlet_in_vacuum() {
  System s = rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m1/singlet_const.rtt.json");
  s.environment.medium = "VACUUM";  // n = 1 exactly, so the OPL below is exact
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

/// Bitwise equality, also for NaN and signed zeros.
bool same_bits(double a, double b) {
  return std::bit_cast<std::uint64_t>(a) == std::bit_cast<std::uint64_t>(b);
}

/// Rays of m1/singlet_const with every status of a normal trace (as hand_filled_rays() of #32):
/// on axis and parallel up to y = 19.5 mm (ALIVE at IMG, VIGNETTED at the stop of radius 10 mm
/// above it), oblique, away from the system (MISSED) and two rays that do not start (ABSORBED).
RayBatch mixed_rays() {
  RayBatch rays(50);
  for (std::size_t i = 0; i < rays.size(); ++i) {
    const double k = static_cast<double>(i);
    rays.pos_z()[i] = -10.0;
    if (i < 40) {
      rays.pos_y()[i] = 0.5 * k;
      rays.wl()[i] = static_cast<std::uint16_t>(i % 3);
    } else if (i < 44) {
      rays.pos_y()[i] = 2.0 * (k - 40.0);
      rays.dir_x()[i] = 0.6;
      rays.dir_z()[i] = 0.8;
      rays.wl()[i] = 1;
    } else if (i < 48) {
      rays.pos_y()[i] = k - 44.0;
      rays.dir_z()[i] = -1.0;
    } else {
      rays.status()[i] = RayStatus::Absorbed;
    }
  }
  return rays;
}

/// Requires that the batch columns of a and b are bitwise equal.
void require_same_batch(const RayBatch& a, const RayBatch& b) {
  REQUIRE(a.size() == b.size());
  for (std::size_t i = 0; i < a.size(); ++i) {
    INFO("ray " << i);
    REQUIRE(same_bits(a.pos_x()[i], b.pos_x()[i]));
    REQUIRE(same_bits(a.pos_y()[i], b.pos_y()[i]));
    REQUIRE(same_bits(a.pos_z()[i], b.pos_z()[i]));
    REQUIRE(same_bits(a.dir_x()[i], b.dir_x()[i]));
    REQUIRE(same_bits(a.dir_y()[i], b.dir_y()[i]));
    REQUIRE(same_bits(a.dir_z()[i], b.dir_z()[i]));
    REQUIRE(same_bits(a.opl()[i], b.opl()[i]));
    REQUIRE(same_bits(a.weight()[i], b.weight()[i]));
    REQUIRE(a.status()[i] == b.status()[i]);
    REQUIRE(a.last_surface()[i] == b.last_surface()[i]);
    REQUIRE(a.prt_matrix(i) == b.prt_matrix(i));
  }
}

/// Requires that recorded ray r of a and recorded ray q of b have bitwise equal slots.
void require_same_record(const RayPaths& a, std::size_t r, const RayPaths& b, std::size_t q) {
  REQUIRE(a.slots == b.slots);
  REQUIRE(a.count[r] == b.count[q]);
  REQUIRE(a.lost_at[r] == b.lost_at[q]);
  for (std::size_t s = 0; s < a.slots; ++s) {
    const std::size_t i = r * a.slots + s;
    const std::size_t j = q * b.slots + s;
    for (std::size_t c = 0; c < 3; ++c) {
      REQUIRE(same_bits(a.position[3 * i + c], b.position[3 * j + c]));
      REQUIRE(same_bits(a.direction[3 * i + c], b.direction[3 * j + c]));
    }
    REQUIRE(same_bits(a.opl[i], b.opl[j]));
    REQUIRE(same_bits(a.weight[i], b.weight[j]));
    REQUIRE(a.status[i] == b.status[j]);
  }
}

/// Requires that two RayPaths are bitwise equal in every member.
void require_same_paths(const RayPaths& a, const RayPaths& b) {
  REQUIRE(a.slots == b.slots);
  REQUIRE(a.ray_indices == b.ray_indices);
  REQUIRE(a.event_surfaces == b.event_surfaces);
  REQUIRE(a.count == b.count);
  REQUIRE(a.lost_at == b.lost_at);
  REQUIRE(a.status == b.status);
  const auto same = [](const std::vector<double>& x, const std::vector<double>& y) {
    if (x.size() != y.size()) return false;
    for (std::size_t i = 0; i < x.size(); ++i) {
      if (!same_bits(x[i], y[i])) return false;
    }
    return true;
  };
  REQUIRE(same(a.position, b.position));
  REQUIRE(same(a.direction, b.direction));
  REQUIRE(same(a.opl, b.opl));
  REQUIRE(same(a.weight, b.weight));
}

}  // namespace

TEST_CASE("axial ray through the singlet: recorded points are the vertices", "[ray_paths]") {
  // Issue #80: the axial ray meets every surface at its vertex (0, 0, z_i), with z_i = 0 (stop),
  // 5 and 9 (lens) and 106.363 mm (image) from m1/singlet_const, and keeps the direction
  // (0, 0, 1). In vacuum with n = 1.5168 inside the lens and the start at z = -10 mm the OPL of
  // the slots is 10, 15, 15 + 4 n and 15 + 4 n + 97.363 mm. The weight of the slots is 1, 1
  // (the stop is a Transmit without effect), 1 - R, (1 - R)^2 and (1 - R)^2 (detector) with
  // R = ((n - 1)/(n + 1))^2 for both uncoated lens surfaces at normal incidence (Byrnes,
  // Eq. (6); ADR 0021). Tolerance 1e-12 mm, 1e-12 for the weight.
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(singlet_in_vacuum(), lib);
  RayBatch rays(1);
  set_ray(rays, 0, Vec3(0.0, 0.0, -10.0), Vec3(0.0, 0.0, 1.0));
  RayPaths paths;
  static_cast<void>(SequentialTracer().trace(cs, PathId{0}, rays, paths));

  const double n = 1.5168;
  const std::vector<double> z = {-10.0, 0.0, 5.0, 9.0, 106.363};
  const std::vector<double> opl = {0.0, 10.0, 15.0, 15.0 + 4.0 * n, 15.0 + 4.0 * n + 97.363};
  const double big_r = std::pow((n - 1.0) / (n + 1.0), 2);
  const std::vector<double> weight = {1.0, 1.0, 1.0 - big_r, std::pow(1.0 - big_r, 2),
                                      std::pow(1.0 - big_r, 2)};
  REQUIRE(paths.slots == 5);
  REQUIRE(paths.ray_count() == 1);
  REQUIRE(paths.ray_indices == std::vector<std::size_t>{0});
  REQUIRE(paths.count[0] == 5);
  REQUIRE(paths.lost_at[0] == -1);
  REQUIRE(paths.event_surfaces.size() == 4);
  for (std::size_t s = 0; s < 5; ++s) {
    INFO("slot " << s);
    REQUIRE((paths.position_at(0, s) - Vec3(0.0, 0.0, z[s])).norm() <= 1e-12);
    REQUIRE((paths.direction_at(0, s) - Vec3(0.0, 0.0, 1.0)).norm() <= 1e-12);
    REQUIRE(std::abs(paths.opl[s] - opl[s]) <= 1e-12);
    REQUIRE(std::abs(paths.weight[s] - weight[s]) <= 1e-12);
    REQUIRE(paths.status[s] == RayStatus::Alive);
    if (s > 0) {
      const auto surface = paths.event_surfaces[s - 1];
      const Vec3 vertex = cs.surfaces()[surface].to_global.apply_point(Vec3::Zero());
      REQUIRE((paths.position_at(0, s) - vertex).norm() <= 1e-12);
    }
  }
}

TEST_CASE("mirror tilted by 45 deg: the recorded direction after the reflection", "[ray_paths]") {
  // Plane mirror at z = 10 mm, rotated by 45 deg about x. Its normal is the rotated z axis
  // R_x(45 deg) z = (0, -sin 45, cos 45) with the right-handed R_x = [[1, 0, 0], [0, c, -s],
  // [0, s, c]] (docs/architecture.md, Transformationen; rtt/math/isometry.hpp), so the ray along
  // +z leaves along d - 2 (d . n) n = (0, 2 sin 45 cos 45, 1 - 2 cos^2 45) = (0, 1, 0)
  // (de Greve, Eq. (13)) from the vertex (0, 0, 10). Tolerance 1e-12.
  // With a mirror aperture of radius 1 mm a second ray at x = 2 mm is vignetted at the only
  // event, which is also the last one: count = S = 2, lost_at = S - 2 = 0.
  System s;
  s.name = "mirror";
  s.environment.medium = "VACUUM";
  s.wavelengths = {{0.55, 1.0, true}};
  s.aperture = {rtt::model::SystemApertureType::EntrancePupilDiameter, Param(10.0)};
  s.fields = {rtt::model::FieldType::AngleDeg, {{0.0, 0.0, 1.0}}};
  s.root.name = "root";
  s.paths = {{"main", true, {}}};
  Surface m;
  m.id = SurfaceId("M");
  m.aperture = rtt::model::CircularAperture{1.0, 0.0};
  Pose pose = Pose::along_z(10.0);
  pose.rotation_deg[0] = Param(45.0);
  s.root.children = {{Element{"M", ElementKind::Mirror, pose, std::nullopt, {m}}}};
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(s, lib);
  RayBatch rays(2);
  set_ray(rays, 1, Vec3(2.0, 0.0, 0.0), Vec3(0.0, 0.0, 1.0));
  RayPaths paths;
  static_cast<void>(SequentialTracer().trace(cs, PathId{0}, rays, paths));
  REQUIRE(paths.slots == 2);
  REQUIRE(paths.count[0] == 2);
  REQUIRE(paths.lost_at[0] == -1);
  REQUIRE(paths.count[1] == 2);
  REQUIRE(paths.lost_at[1] == 0);
  REQUIRE(paths.status[3] == RayStatus::Vignetted);
  REQUIRE((paths.position_at(1, 1) - Vec3(2.0, 0.0, 10.0)).norm() <= 1e-12);
  REQUIRE((paths.position_at(0, 0) - Vec3(0.0, 0.0, 0.0)).norm() <= 1e-12);
  REQUIRE((paths.direction_at(0, 0) - Vec3(0.0, 0.0, 1.0)).norm() <= 1e-12);
  REQUIRE((paths.position_at(0, 1) - Vec3(0.0, 0.0, 10.0)).norm() <= 1e-12);
  REQUIRE((paths.direction_at(0, 1) - Vec3(0.0, 1.0, 0.0)).norm() <= 1e-12);
  REQUIRE(std::abs(paths.opl[1] - 10.0) <= 1e-12);
}

TEST_CASE("recording does not change the trace; the last valid slot is the final state",
          "[ray_paths]") {
  // Issue #80: without recording the hot path is unchanged, with recording the batch ends
  // bitwise the same, and slot count - 1 is bitwise the final state of the ray.
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(singlet_in_vacuum(), lib);
  const RayBatch start = mixed_rays();
  RayBatch plain = start;
  RayBatch recorded = start;
  const auto plain_stats = SequentialTracer().trace(cs, PathId{0}, plain);
  RayPaths paths;
  const auto recorded_stats = SequentialTracer().trace(cs, PathId{0}, recorded, paths);
  REQUIRE(plain_stats.rays == recorded_stats.rays);
  require_same_batch(plain, recorded);
  REQUIRE(plain_stats.count(RayStatus::Alive) > 0);
  REQUIRE(plain_stats.count(RayStatus::Vignetted) > 0);
  REQUIRE(plain_stats.count(RayStatus::Missed) > 0);
  REQUIRE(plain_stats.count(RayStatus::Absorbed) == 2);

  REQUIRE(paths.ray_count() == start.size());
  for (std::size_t r = 0; r < paths.ray_count(); ++r) {
    INFO("ray " << r);
    REQUIRE(paths.ray_indices[r] == r);
    const std::size_t last = r * paths.slots + paths.count[r] - 1;
    REQUIRE(same_bits(paths.position[3 * last], plain.pos_x()[r]));
    REQUIRE(same_bits(paths.position[3 * last + 1], plain.pos_y()[r]));
    REQUIRE(same_bits(paths.position[3 * last + 2], plain.pos_z()[r]));
    REQUIRE(same_bits(paths.direction[3 * last], plain.dir_x()[r]));
    REQUIRE(same_bits(paths.direction[3 * last + 1], plain.dir_y()[r]));
    REQUIRE(same_bits(paths.direction[3 * last + 2], plain.dir_z()[r]));
    REQUIRE(same_bits(paths.opl[last], plain.opl()[r]));
    REQUIRE(same_bits(paths.weight[last], plain.weight()[r]));
    REQUIRE(paths.status[last] == plain.status()[r]);
    // Slot 0 is the start state.
    REQUIRE(same_bits(paths.position[3 * r * paths.slots + 1], start.pos_y()[r]));
    REQUIRE(paths.status[r * paths.slots] == start.status()[r]);
  }
}

TEST_CASE("lost rays: count, lost_at and NaN after the loss", "[ray_paths]") {
  // Events of m1/singlet_const: 0 stop, 1 and 2 lens, 3 image. A parallel ray at y = 15 mm is
  // vignetted at the stop (event 0): count 2, lost_at 0, slot 1 at the hit point (0, 15, 0)
  // with status Vignetted, slots 2 to 4 NaN with status Vignetted. A ray along -z misses the
  // stop: count 2, lost_at 0, slot 1 = start state with status Missed. A ray that is not Alive
  // at the start is not traced: count 1, lost_at -1, NaN from slot 1. A ray from (0, 2.4, -10)
  // along (0, 0.6, 0.8) passes the stop at y = 2.4 + 0.75 * 10 = 9.9 mm < 10 mm and reaches the
  // first lens surface (vertex z = 5, radius 51.68) at y > 2.4 + 0.75 * 15 = 13.65 mm > 12.7 mm:
  // vignetted at event 1, count 3, lost_at 1.
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(singlet_in_vacuum(), lib);
  RayBatch rays(4);
  set_ray(rays, 0, Vec3(0.0, 15.0, -10.0), Vec3(0.0, 0.0, 1.0));
  set_ray(rays, 3, Vec3(0.0, 2.4, -10.0), Vec3(0.0, 0.6, 0.8));
  set_ray(rays, 1, Vec3(0.0, 1.0, -10.0), Vec3(0.0, 0.0, -1.0));
  set_ray(rays, 2, Vec3(0.0, 0.0, -10.0), Vec3(0.0, 0.0, 1.0));
  rays.status()[2] = RayStatus::Absorbed;
  RayPaths paths;
  static_cast<void>(SequentialTracer().trace(cs, PathId{0}, rays, paths));
  REQUIRE(paths.slots == 5);

  REQUIRE(paths.count[0] == 2);
  REQUIRE(paths.lost_at[0] == 0);
  REQUIRE((paths.position_at(0, 1) - Vec3(0.0, 15.0, 0.0)).norm() <= 1e-12);
  REQUIRE(paths.status[1] == RayStatus::Vignetted);

  REQUIRE(paths.count[1] == 2);
  REQUIRE(paths.lost_at[1] == 0);
  REQUIRE(paths.position_at(1, 1) == Vec3(0.0, 1.0, -10.0));
  REQUIRE(paths.status[paths.slots + 1] == RayStatus::Missed);

  REQUIRE(paths.count[2] == 1);
  REQUIRE(paths.lost_at[2] == -1);
  REQUIRE(paths.status[2 * paths.slots] == RayStatus::Absorbed);

  REQUIRE(paths.count[3] == 3);
  REQUIRE(paths.lost_at[3] == 1);
  REQUIRE(paths.status[3 * paths.slots + 1] == RayStatus::Alive);
  REQUIRE(std::abs(paths.position_at(3, 1).y() - 9.9) <= 1e-12);
  REQUIRE(paths.status[3 * paths.slots + 2] == RayStatus::Vignetted);
  REQUIRE(paths.position_at(3, 2).y() > 13.65);

  for (std::size_t r = 0; r < 4; ++r) {
    for (std::size_t s = paths.count[r]; s < paths.slots; ++s) {
      INFO("ray " << r << ", slot " << s);
      const std::size_t i = r * paths.slots + s;
      REQUIRE(std::isnan(paths.position[3 * i]));
      REQUIRE(std::isnan(paths.position[3 * i + 1]));
      REQUIRE(std::isnan(paths.position[3 * i + 2]));
      REQUIRE(std::isnan(paths.direction[3 * i]));
      REQUIRE(std::isnan(paths.opl[i]));
      REQUIRE(std::isnan(paths.weight[i]));
      REQUIRE(paths.status[i] == paths.status[r * paths.slots + paths.count[r] - 1]);
    }
  }
}

TEST_CASE("record_rays selects rays of a hexapolar bundle, bitwise as in the full record",
          "[ray_paths]") {
  // Decided for #80: an explicit selection instead of truncation; the selected rays are
  // recorded in the given order and equal the full record of the same rays bitwise.
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(singlet_in_vacuum(), lib);
  const std::vector<std::uint16_t> fields = {0, 1, 2};
  const RayBatch start = rtt::trace::make_rays(
      cs, PathId{0}, fields, 1, rtt::trace::HexapolarPupil{4}, rtt::trace::Aiming::Real);
  const std::vector<std::size_t> selection = {start.size() - 1, 0, start.size() / 2, 7};
  RayBatch full_rays = start;
  RayBatch selected_rays = start;
  RayPaths full;
  RayPaths selected;
  static_cast<void>(SequentialTracer().trace(cs, PathId{0}, full_rays, full));
  static_cast<void>(SequentialTracer().trace(cs, PathId{0}, selected_rays, selected, selection));
  require_same_batch(full_rays, selected_rays);
  REQUIRE(selected.ray_indices == selection);
  REQUIRE(selected.ray_count() == selection.size());
  REQUIRE(selected.event_surfaces == full.event_surfaces);
  for (std::size_t r = 0; r < selection.size(); ++r) {
    INFO("selected " << r);
    require_same_record(selected, r, full, selection[r]);
  }
}

TEST_CASE("record_rays and the size limit are checked before tracing", "[ray_paths]") {
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(singlet_in_vacuum(), lib);
  RayBatch rays(5);
  RayPaths paths;
  const SequentialTracer tracer;
  // A previous record that the rejected calls below must leave untouched (strong guarantee).
  RayBatch earlier(3);
  static_cast<void>(tracer.trace(cs, PathId{0}, earlier, paths));
  const RayPaths snapshot = paths;
  const std::vector<std::size_t> duplicate = {1, 3, 1};
  const std::vector<std::size_t> outside = {0, 5};
  REQUIRE_THROWS_AS(static_cast<void>(tracer.trace(cs, PathId{0}, rays, paths, duplicate)),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(static_cast<void>(tracer.trace(cs, PathId{0}, rays, paths, outside)),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(static_cast<void>(tracer.trace(cs, PathId{0}, rays, paths, {}, 4)),
                    std::invalid_argument);
  // An invalid wavelength is found by the trace itself, after the record was allocated.
  RayBatch bad_wavelength(2);
  bad_wavelength.wl()[1] = 99;
  REQUIRE_THROWS_AS(static_cast<void>(tracer.trace(cs, PathId{0}, bad_wavelength, paths)),
                    std::invalid_argument);
  // Nothing was traced by the rejected calls, and the earlier record is unchanged.
  REQUIRE(rays.last_surface()[0] == rtt::trace::kNoSurface);
  REQUIRE(bad_wavelength.last_surface()[0] == rtt::trace::kNoSurface);
  require_same_paths(paths, snapshot);
  // A selection may exceed nothing: 2 rays of 5 with a limit of 4 are fine.
  const std::vector<std::size_t> two = {4, 2};
  static_cast<void>(tracer.trace(cs, PathId{0}, rays, paths, two, 4));
  REQUIRE(paths.ray_count() == 2);
  // The limit message names the memory need and record_rays.
  RayBatch more(5);
  try {
    static_cast<void>(tracer.trace(cs, PathId{0}, more, paths, {}, 4));
    FAIL("no exception");
  } catch (const std::invalid_argument& e) {
    const std::string message = e.what();
    REQUIRE(message.find("record_rays") != std::string::npos);
    REQUIRE(message.find("bytes") != std::string::npos);
  }
}

TEST_CASE("recorded paths are bitwise identical with one thread and with all threads",
          "[ray_paths]") {
  // ADR 0004: every ray writes only its own slots.
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(singlet_in_vacuum(), lib);
  const std::vector<std::uint16_t> fields = {0, 1, 2};
  const RayBatch start = rtt::trace::make_rays(cs, PathId{0}, fields, 1, rtt::trace::GridPupil{31},
                                               rtt::trace::Aiming::Real);
  RayBatch serial_rays = start;
  RayBatch parallel_rays = start;
  RayPaths serial;
  RayPaths parallel;
  oneapi::tbb::task_arena one(1);
  one.execute(
      [&] { static_cast<void>(SequentialTracer().trace(cs, PathId{0}, serial_rays, serial)); });
  static_cast<void>(SequentialTracer().trace(cs, PathId{0}, parallel_rays, parallel));
  require_same_batch(serial_rays, parallel_rays);
  REQUIRE(serial.ray_count() == parallel.ray_count());
  for (std::size_t r = 0; r < serial.ray_count(); ++r) {
    require_same_record(serial, r, parallel, r);
  }
}
