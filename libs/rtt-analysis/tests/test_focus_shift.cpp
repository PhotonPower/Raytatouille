// Through focus (ADR 0032, point 5; #202): OpdOptions::focus_shift moves the centre of the
// reference sphere to C' = C + Delta d_c along the chief ray of the reference wavelength.
//
// Expected value (derived in ADR 0032, point 5): a perfect spherical wave converging to C. A ray
// with direction u, cos(alpha) = u . d_c, crosses the sphere about C' (radius R' = |C' - E|) at
// s = -Delta cos(alpha) + sqrt(R'^2 - Delta^2 sin^2(alpha)) before C, the chief ray at
// s_c = R' - Delta, so
//   W = n' (s - s_c) = n' [Delta (1 - cos(alpha)) - Delta^2 sin^2(alpha) / (sqrt(R'^2 -
//       Delta^2 sin^2(alpha)) + R')]                                    (mm; waves: / lambda_ref)
// and n' Delta (1 - cos(alpha)) for R' = infinity. The tests evaluate it without cancellation:
// 1 - cos(alpha) = sin^2(alpha) / (1 + cos(alpha)), sin^2(alpha) = |u x d_c|^2.
//
// The real wave of the paraboloid is spherical only up to the rounding of the trace, so the
// tests compare W(Delta) - W(0) with the formula: the OPL error of each ray is the same at both
// Delta and drops out. Tolerance 1e-9 waves (5.9e-13 mm): the OPL of a ray is up to about 250 mm,
// so each of the few roundings in the OPL to the sphere is about u 250 mm = 2.8e-14 mm; the
// formula itself rounds at u times |W| <= 0.25 mm. 1e-9 waves is about 20 times that.

#include <Eigen/Geometry>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#include "rtt/analysis/opd.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
#include "rtt/paraxial/paraxial.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/sequential.hpp"
#include "rtt/trace/sources.hpp"

using rtt::analysis::OpdOptions;
using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::math::Vec3;
using rtt::model::Element;
using rtt::model::Param;
using rtt::model::System;

namespace {

constexpr PathId kPath{0};

System load(const std::string& name) {
  return rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/" + name);
}

CompiledSystem compiled(const System& s) {
  return rtt::compile::compile(s, rtt::material::MaterialLibrary{});
}

bool same_bits(double a, double b) {
  return std::memcmp(&a, &b, sizeof(double)) == 0;
}

/// The rays of opd_points: make_rays and the sequential tracer, as in the analyses.
rtt::trace::RayBatch traced(const CompiledSystem& cs, const rtt::trace::PupilSampling& sampling) {
  const std::uint16_t field = 0;
  rtt::trace::RayBatch rays =
      rtt::trace::make_rays(cs, kPath, std::span(&field, 1), cs.reference_wavelength(), sampling);
  static_cast<void>(rtt::trace::SequentialTracer().trace(cs, kPath, rays));
  return rays;
}

Vec3 pos(const rtt::trace::RayBatch& r, std::size_t i) {
  return {r.pos_x()[i], r.pos_y()[i], r.pos_z()[i]};
}
Vec3 dir(const rtt::trace::RayBatch& r, std::size_t i) {
  return {r.dir_x()[i], r.dir_y()[i], r.dir_z()[i]};
}

/// Image-space index n' = |Re n| after the last event at the reference wavelength.
double image_index(const CompiledSystem& cs) {
  const auto& last = cs.path(kPath).events.back();
  return std::abs(cs.media()[last.medium_after].index[cs.reference_wavelength()].real());
}

double lambda_mm(const CompiledSystem& cs) {
  return cs.wavelengths_um()[cs.reference_wavelength()] * 1e-3;
}

/// W of the formula above in mm for direction u, chief direction dc, shift delta, radius r.
double expected_w_mm(const Vec3& u, const Vec3& dc, double delta, double r, double n) {
  const double cos_a = u.dot(dc);
  const double sin2 = u.cross(dc).squaredNorm();
  const double one_minus_cos = sin2 / (1.0 + cos_a);
  if (std::isinf(r)) return n * delta * one_minus_cos;
  const double d2s2 = delta * delta * sin2;
  return n * (delta * one_minus_cos - d2s2 / (std::sqrt(r * r - d2s2) + r));
}

/// Checks W(Delta) - W(0) against the formula on the hexapolar points of `cs` (field 0).
void check_against_formula(const CompiledSystem& cs, double delta) {
  const rtt::trace::HexapolarPupil sampling{6};
  OpdOptions shifted;
  shifted.focus_shift = delta;
  const auto w0 = rtt::analysis::opd_points(cs, kPath, 0, cs.reference_wavelength(), sampling);
  const auto wd =
      rtt::analysis::opd_points(cs, kPath, 0, cs.reference_wavelength(), sampling, shifted);
  const rtt::trace::RayBatch rays = traced(cs, sampling);
  const rtt::trace::RayBatch chief = traced(cs, rtt::trace::SinglePupilPoint{0.0, 0.0});
  const Vec3 c = pos(chief, 0);
  const Vec3 dc = dir(chief, 0);
  // The sphere: C' = C + Delta d_c (the same expression, so bit for bit), R' = |C' - E|.
  const Vec3 c_shift = c + delta * dc;
  CHECK(wd.sphere.centre == c_shift);
  const auto fo = rtt::paraxial::first_order(cs, kPath, cs.reference_wavelength());
  REQUIRE(fo.exit_pupil);
  const double r = fo.exit_pupil->z ? (c_shift - Vec3(0.0, 0.0, *fo.exit_pupil->z)).norm()
                                    : std::numeric_limits<double>::infinity();
  CHECK(same_bits(wd.sphere.radius, r));
  REQUIRE(wd.points.size() == rays.size());
  const double n = image_index(cs);
  std::size_t arrived = 0;
  double largest = 0.0;
  for (std::size_t i = 0; i < rays.size(); ++i) {
    INFO("point " << i << " (" << wd.points[i].px << ", " << wd.points[i].py << ")");
    REQUIRE(w0.points[i].status == rtt::trace::RayStatus::Alive);
    REQUIRE(wd.points[i].status == rtt::trace::RayStatus::Alive);
    const double expected = expected_w_mm(dir(rays, i), dc, delta, r, n) / lambda_mm(cs);
    CHECK(std::abs(wd.points[i].w - w0.points[i].w - expected) <= 1e-9);
    largest = std::max(largest, std::abs(expected));
    ++arrived;
  }
  // Content guard: all points arrive and the defocus is many waves at the rim.
  CHECK(arrived == rays.size());
  CHECK(largest > 1.0);
}

}  // namespace

TEST_CASE("focus_shift: W of a paraboloid against the through-focus formula, finite exit pupil",
          "[opd][focus]") {
  // m2/paraboloid_stop: f = 100 mm, EPD 60, stop at z = -50 between mirror and focus, so the
  // exit pupil is the virtual image of the stop behind the mirror (finite R).
  const CompiledSystem cs = compiled(load("m2/paraboloid_stop.rtt.json"));
  REQUIRE(rtt::paraxial::first_order(cs, kPath, 0).exit_pupil->z);
  for (const double delta : {-0.7, 0.3, 5.0}) {
    INFO("Delta = " << delta);
    check_against_formula(cs, delta);
  }
}

TEST_CASE("focus_shift: W of a paraboloid against the formula, exit pupil at infinity",
          "[opd][focus]") {
  // The stop moved into the focal plane z = -100: its image by the mirror lies at infinity
  // (image-space telecentric, R = infinity, #102).
  System s = load("m2/paraboloid_stop.rtt.json");
  std::get<Element>(s.root.children[0].value).pose.position[2] = Param(-100.0);
  const CompiledSystem cs = compiled(s);
  REQUIRE_FALSE(rtt::paraxial::first_order(cs, kPath, 0).exit_pupil->z);
  for (const double delta : {-0.7, 0.3, 5.0}) {
    INFO("Delta = " << delta);
    check_against_formula(cs, delta);
  }
}

TEST_CASE("focus_shift: the same W as an image surface moved by Delta and traced again",
          "[opd][focus]") {
  // The counterpart of #202: on the axis of the singlet the chief ray runs along +z, so moving
  // the detector by Delta puts its chief-ray point exactly on C + Delta d_c with the same
  // radius; the OPL up to the sphere is the same in the homogeneous image space. Only the
  // rounding differs: OPL up to about 120 mm, a few roundings of u 120 mm = 1.3e-14 mm each;
  // tolerance 1e-9 waves = 5.9e-13 mm.
  const double delta = 0.5;
  const System base = load("m1/singlet_const.rtt.json");
  System moved = base;
  bool found = false;
  for (auto& node : moved.root.children) {
    auto* e = std::get_if<Element>(&node.value);
    if (e != nullptr && e->kind == rtt::model::ElementKind::Detector) {
      e->pose.position[2] = Param(e->pose.position[2].value + delta);
      found = true;
    }
  }
  REQUIRE(found);
  const CompiledSystem cs = compiled(base);
  const CompiledSystem cs_moved = compiled(moved);
  const rtt::trace::HexapolarPupil sampling{6};
  OpdOptions shifted;
  shifted.focus_shift = delta;
  const auto a = rtt::analysis::opd_points(cs, kPath, 0, 0, sampling, shifted);
  const auto b = rtt::analysis::opd_points(cs_moved, kPath, 0, 0, sampling);
  CHECK((a.sphere.centre - b.sphere.centre).norm() <= 1e-12);
  CHECK(std::abs(a.sphere.radius - b.sphere.radius) <= 1e-12);
  REQUIRE(a.points.size() == b.points.size());
  double largest = 0.0;
  for (std::size_t i = 0; i < a.points.size(); ++i) {
    INFO("point " << i);
    REQUIRE(a.points[i].status == rtt::trace::RayStatus::Alive);
    REQUIRE(b.points[i].status == rtt::trace::RayStatus::Alive);
    CHECK(std::abs(a.points[i].w - b.points[i].w) <= 1e-9);
    largest = std::max(largest, std::abs(a.points[i].w));
  }
  CHECK(largest > 0.1);  // content guard: an aberrated, defocused wavefront
}

TEST_CASE("focus_shift = 0 gives the sphere and W of before bit for bit", "[opd][focus]") {
  // The default and an explicit 0 (also -0.0) take the same path: no shift is applied.
  const CompiledSystem cs = compiled(load("m1/singlet_const.rtt.json"));
  OpdOptions zero;
  zero.focus_shift = 0.0;
  OpdOptions minus_zero;
  minus_zero.focus_shift = -0.0;
  for (const std::uint16_t field : {std::uint16_t{0}, std::uint16_t{1}}) {
    const auto a = rtt::analysis::opd_map(cs, kPath, field, 0);
    for (const OpdOptions& o : {zero, minus_zero}) {
      const auto b = rtt::analysis::opd_map(cs, kPath, field, 0, o);
      CHECK(a.sphere.centre == b.sphere.centre);
      CHECK(same_bits(a.sphere.radius, b.sphere.radius));
      REQUIRE(a.points.size() == b.points.size());
      for (std::size_t i = 0; i < a.points.size(); ++i) {
        CHECK(same_bits(a.points[i].w, b.points[i].w));
        CHECK(a.points[i].status == b.points[i].status);
      }
    }
  }
}

TEST_CASE("focus_shift: not finite is invalid_argument, C' on the exit pupil an AnalysisError",
          "[opd][focus]") {
  const CompiledSystem cs = compiled(load("m2/paraboloid_stop.rtt.json"));
  for (const double bad :
       {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(),
        -std::numeric_limits<double>::infinity()}) {
    OpdOptions o;
    o.focus_shift = bad;
    CHECK_THROWS_AS(rtt::analysis::opd_map(cs, kPath, 0, 0, o), std::invalid_argument);
    CHECK_THROWS_AS(rtt::analysis::opd_fan(cs, kPath, 0, 0, o), std::invalid_argument);
    CHECK_THROWS_AS(rtt::analysis::opd_points(cs, kPath, 0, 0, rtt::trace::HexapolarPupil{1}, o),
                    std::invalid_argument);
  }
  // Delta = d_c . (E - C) puts C' exactly on the exit-pupil centre E on the axis.
  const rtt::trace::RayBatch chief = traced(cs, rtt::trace::SinglePupilPoint{0.0, 0.0});
  const Vec3 c = pos(chief, 0);
  const Vec3 dc = dir(chief, 0);
  const Vec3 e(0.0, 0.0, *rtt::paraxial::first_order(cs, kPath, 0).exit_pupil->z);
  OpdOptions o;
  o.focus_shift = dc.dot(e - c);
  REQUIRE(c + o.focus_shift * dc == e);  // precondition: exact in double here
  CHECK_THROWS_AS(rtt::analysis::opd_map(cs, kPath, 0, 0, o), rtt::analysis::AnalysisError);
}
