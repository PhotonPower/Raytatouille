// Virtual reference rays of centrally obscured systems (#229; ADR 0032, addendum #229).
//
// RayBatch::set_ignore_obscurations(true) makes the sequential tracer skip the inner rim of
// circular apertures (inner_radius > 0, a central obscuration) and keep every outer rim, so the
// chief ray of an obscured pupil reaches the image as the ray of the same system without the
// obscuration. stop_radius() is |R_s|, the stop radius of pupil coordinates |p| = 1: with real
// aiming a ray of pupil point p meets the stop at |p| |R_s|, so a centred obscuration of radius
// r_i at the stop is the pupil circle |p| = r_i / |R_s|.
//
// System: tests/reference/m2/paraboloid_stop (stop in front of a paraboloid, f = 100 mm, object
// at infinity) with EPD 10 mm and the stop aperture {5 mm, inner 2 mm}: obscuration 0.4 in
// pupil units.

#include <array>
#include <bit>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/math/types.hpp"
#include "rtt/model/model.hpp"
#include "rtt/trace/apply_event.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/ray_paths.hpp"
#include "rtt/trace/sequential.hpp"
#include "rtt/trace/sources.hpp"

using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::model::Element;
using rtt::model::System;
using rtt::trace::PointListPupil;
using rtt::trace::PupilPoint;
using rtt::trace::RayBatch;
using rtt::trace::RayStatus;

namespace {

System load(const std::string& name) {
  return rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/" + name + ".rtt.json");
}

System paraboloid(double inner_radius) {
  System s = load("m2/paraboloid_stop");
  s.aperture.value = rtt::model::Param(10.0);
  std::get<Element>(s.root.children[0].value).surfaces[0].aperture =
      rtt::model::CircularAperture{5.0, inner_radius};
  return s;
}

CompiledSystem compiled(const System& s) {
  return rtt::compile::compile(s, rtt::material::MaterialLibrary{});
}

/// The M2 systems with their glass catalogue (tests/catalogs/m2/schott.agf).
CompiledSystem compiled_m2(const System& s) {
  rtt::material::MaterialLibrary lib;
  lib.add_catalog(std::string(RTT_CATALOG_DIR) + "/m2/schott.agf");
  return rtt::compile::compile(s, lib);
}

/// Two pupil points in the obscuration (|p| < 0.4), two in the annulus, one beyond the rim.
constexpr std::array<PupilPoint, 5> kPoints{
    {{0.0, 0.0}, {0.2, -0.1}, {0.5, 0.0}, {0.0, -0.9}, {1.2, 0.0}}};

RayBatch traced(const CompiledSystem& cs, bool virtual_rays) {
  const std::array<std::uint16_t, 1> fields{0};
  const PointListPupil list{std::vector<PupilPoint>(kPoints.begin(), kPoints.end())};
  RayBatch rays = rtt::trace::make_rays(cs, PathId{0}, fields, 0, list);
  rays.set_ignore_obscurations(virtual_rays);
  [[maybe_unused]] const auto stats = rtt::trace::SequentialTracer().trace(cs, PathId{0}, rays);
  return rays;
}

bool same_bits(double a, double b) {
  return std::bit_cast<std::uint64_t>(a) == std::bit_cast<std::uint64_t>(b);
}

/// Every traced column of ray i is bitwise the same in both batches.
bool same_ray(const RayBatch& a, const RayBatch& b, std::size_t i) {
  const auto eq = [i](std::span<const double> x, std::span<const double> y) {
    return same_bits(x[i], y[i]);
  };
  const rtt::math::CMat3 pa = a.prt_matrix(i);
  const rtt::math::CMat3 pb = b.prt_matrix(i);
  bool same_p = true;
  for (int r = 0; r < 3; ++r) {
    for (int c = 0; c < 3; ++c) {
      same_p = same_p && same_bits(pa(r, c).real(), pb(r, c).real()) &&
               same_bits(pa(r, c).imag(), pb(r, c).imag());
    }
  }
  return a.status()[i] == b.status()[i] && a.last_surface()[i] == b.last_surface()[i] &&
         eq(a.pos_x(), b.pos_x()) && eq(a.pos_y(), b.pos_y()) && eq(a.pos_z(), b.pos_z()) &&
         eq(a.dir_x(), b.dir_x()) && eq(a.dir_y(), b.dir_y()) && eq(a.dir_z(), b.dir_z()) &&
         eq(a.opl(), b.opl()) && eq(a.weight(), b.weight()) && same_p;
}

std::uint32_t stop_surface(const CompiledSystem& cs) {
  for (const auto& e : cs.path(PathId{0}).events) {
    if (cs.surfaces()[e.surface].element_kind == rtt::model::ElementKind::Stop) return e.surface;
  }
  FAIL("no stop on the path");
  return 0;
}

/// Real aiming puts pupil point p on the stop at |p| R_s, within the aim tolerance
/// kAimTolerance = 1e-9 mm (the aim residual is the distance to the target, which bounds the
/// radial deviation): checked on the recorded path at three points of field 0.
void check_stop_hits(const CompiledSystem& cs, double r_s) {
  REQUIRE(r_s > 0.0);
  const std::vector<PupilPoint> points = {{0.4, 0.0}, {0.0, -1.0}, {0.3, 0.6}};
  const std::array<std::uint16_t, 1> fields{0};
  RayBatch rays = rtt::trace::make_rays(cs, PathId{0}, fields, 0, PointListPupil{points});
  rtt::trace::RayPaths paths;
  [[maybe_unused]] const auto stats =
      rtt::trace::SequentialTracer().trace(cs, PathId{0}, rays, paths);
  const std::uint32_t stop = stop_surface(cs);
  std::size_t slot = 0;
  for (std::size_t k = 0; k < paths.event_surfaces.size(); ++k) {
    if (paths.event_surfaces[k] == stop) slot = k + 1;  // slot k + 1 is after event k
  }
  REQUIRE(slot > 0);
  for (std::size_t r = 0; r < points.size(); ++r) {
    const rtt::math::Vec3 local =
        cs.surfaces()[stop].to_local.apply_point(paths.position_at(r, slot));
    const double expected = std::hypot(points[r].px, points[r].py) * r_s;
    INFO("point " << r << ": stop radius " << std::hypot(local.x(), local.y()) << ", expected "
                  << expected);
    CHECK(std::abs(std::hypot(local.x(), local.y()) - expected) <= 1e-9 + 1e-12 * r_s);
  }
}

}  // namespace

TEST_CASE("RayBatch: virtual reference rays only on request", "[obscuration]") {
  RayBatch rays(3);
  CHECK_FALSE(rays.ignores_obscurations());  // a new batch traces every rim
  rays.set_ignore_obscurations(true);
  CHECK(rays.ignores_obscurations());
  rays.resize(5);
  CHECK(rays.ignores_obscurations());  // a property of the batch, not of its rows
  const RayBatch copy = rays;
  CHECK(copy.ignores_obscurations());
  rays.set_ignore_obscurations(false);
  CHECK_FALSE(rays.ignores_obscurations());
}

TEST_CASE("real rays: the central obscuration vignettes at the stop", "[obscuration]") {
  const CompiledSystem cs = compiled(paraboloid(2.0));
  const std::uint32_t stop = stop_surface(cs);
  const std::uint32_t image = cs.path(PathId{0}).events.back().surface;
  const RayBatch rays = traced(cs, false);
  for (const std::size_t i : {std::size_t{0}, std::size_t{1}, std::size_t{4}}) {
    INFO("ray " << i);
    CHECK(rays.status()[i] == RayStatus::Vignetted);
    CHECK(rays.last_surface()[i] == stop);
  }
  for (const std::size_t i : {std::size_t{2}, std::size_t{3}}) {
    INFO("ray " << i);
    CHECK(rays.status()[i] == RayStatus::Alive);
    CHECK(rays.last_surface()[i] == image);
  }
}

TEST_CASE("virtual rays skip the inner rim and keep the outer rim", "[obscuration]") {
  const CompiledSystem cs = compiled(paraboloid(2.0));
  const std::uint32_t stop = stop_surface(cs);
  const std::uint32_t image = cs.path(PathId{0}).events.back().surface;
  const RayBatch real = traced(cs, false);
  const RayBatch virt = traced(cs, true);
  CHECK(virt.ignores_obscurations());  // the trace keeps the request
  // The rays in the obscuration arrive; the ray beyond the outer rim is still vignetted there.
  for (const std::size_t i : {std::size_t{0}, std::size_t{1}}) {
    INFO("ray " << i);
    CHECK(virt.status()[i] == RayStatus::Alive);
    CHECK(virt.last_surface()[i] == image);
  }
  CHECK(virt.status()[4] == RayStatus::Vignetted);
  CHECK(virt.last_surface()[4] == stop);
  // Rays in the annulus never meet the inner rim: bitwise the same as real rays.
  for (const std::size_t i : {std::size_t{2}, std::size_t{3}}) {
    INFO("ray " << i);
    CHECK(same_ray(real, virt, i));
  }
  // Every virtual ray is bitwise the real ray of the same system without the obscuration: only
  // the test r^2 >= r_i^2 is skipped (ADR 0032, addendum #229).
  const RayBatch open = traced(compiled(paraboloid(0.0)), false);
  REQUIRE(open.status()[0] == RayStatus::Alive);  // content guard: the open chief ray arrives
  for (std::size_t i = 0; i < kPoints.size(); ++i) {
    INFO("ray " << i);
    CHECK(same_ray(open, virt, i));
  }
  // The virtual chief ray of the axial field meets the focus on the axis (content guard).
  CHECK(std::abs(virt.pos_x()[0]) <= 1e-12);
  CHECK(std::abs(virt.pos_y()[0]) <= 1e-12);
}

TEST_CASE("stop_radius: the stop radius of pupil coordinates |p| = 1", "[obscuration]") {
  // Paraboloid: the stop is the first surface in object space, the entrance pupil itself.
  const CompiledSystem para = compiled(paraboloid(2.0));
  CHECK(std::abs(rtt::trace::stop_radius(para, PathId{0}, 0) - 5.0) <= 1e-12 * 5.0);
  // stop_size: |p| = 1 is the stop rim (radius 4 mm), whatever the lens before the stop.
  const CompiledSystem two = compiled(load("m1/two_lenses_stop_between"));
  CHECK(std::abs(rtt::trace::stop_radius(two, PathId{0}, 0) - 4.0) <= 1e-12 * 4.0);
  // Cooke triplet (EPD 12.5 mm, stop between the lenses).
  const CompiledSystem cooke = compiled_m2(load("m2/cooke_triplet"));
  const double r_s = rtt::trace::stop_radius(cooke, PathId{0}, 0);
  CHECK(std::abs(r_s - 6.25) > 1e-3);  // content guard: the stop is not the entrance pupil
  check_stop_hits(cooke, r_s);
  // Object-space telecentric singlet (object_na 0.05, stop in the rear focal plane): pupil
  // coordinates are object-space slopes, R_s = u_m |s| with the stop height s per unit slope
  // (review of #238).
  const CompiledSystem tele = compiled(load("m1/telecentric_singlet"));
  check_stop_hits(tele, rtt::trace::stop_radius(tele, PathId{0}, 0));
  CHECK_THROWS_AS(rtt::trace::stop_radius(cooke, PathId{7}, 0), std::invalid_argument);
  CHECK_THROWS_AS(rtt::trace::stop_radius(cooke, PathId{0}, 9), std::invalid_argument);
}
