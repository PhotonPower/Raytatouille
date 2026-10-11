// Pupil function (ADR 0032, points 3, 4, 5, 8, 9, 10; #202 PR B).
//
// Main system: m2/paraboloid_stop, a paraboloid f = 100 mm, EPD 60 (r = 30 mm), stop at z = -50,
// focus (image surface) at z = -100, object at infinity. On the axis its wave is spherical up to
// the rounding of the trace (W < 1e-6 waves as in the mandatory OPD case of #29), the mirror is
// an ideal conductor (weight 1).
//
// Mapping (ADR 0032, point 4): a ray at height h = rho r reflects towards the focus at the
// distance L = f + h^2 / (4 f) (paraboloid), so |s| = n' h / L = n' (h / f) / (1 + g) with
// g = (h / 2f)^2, while M p = n' (h / f) p / |p| (M is the derivative at p = 0). Hence
// |s - M p| / ||M||_2 = rho g / (1 + g) for every sample point. M itself comes from central
// differences at p = +-1e-3: its relative error is about (1e-3 r / 2f)^2 = 2.3e-8; the deltas are
// compared to 1e-5 relative, |M| against n' r / f' (paraxial) to 1e-6 relative.

#include <tbb/global_control.h>

#include <Eigen/Geometry>
#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#include "rtt/analysis/opd.hpp"
#include "rtt/analysis/pupil.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
#include "rtt/paraxial/paraxial.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/sequential.hpp"
#include "rtt/trace/sources.hpp"

using rtt::analysis::AnalysisError;
using rtt::analysis::PupilFunction;
using rtt::analysis::PupilFunctionOptions;
using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::math::Vec3;
using rtt::model::Element;
using rtt::model::Param;
using rtt::model::System;
using rtt::trace::RayStatus;

namespace {

constexpr PathId kPath{0};

System load(const std::string& name) {
  return rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/" + name);
}

CompiledSystem compiled(const System& s) {
  rtt::material::MaterialLibrary lib;
  lib.add_catalog(std::string(RTT_CATALOG_DIR) + "/schott.agf");
  return rtt::compile::compile(s, lib);
}

bool same_bits(double a, double b) {
  return std::memcmp(&a, &b, sizeof(double)) == 0;
}

bool has(const std::vector<rtt::model::Diagnostic>& w, const std::string& code) {
  return std::any_of(w.begin(), w.end(), [&](const auto& d) { return d.code == code; });
}

/// All columns, the sphere, the mapping, the losses and the warnings bit for bit.
void check_same(const PupilFunction& a, const PupilFunction& b) {
  REQUIRE(a.fraction.size() == b.fraction.size());
  CHECK(a.sphere.centre == b.sphere.centre);
  CHECK(same_bits(a.sphere.radius, b.sphere.radius));
  CHECK(same_bits(a.focus_shift, b.focus_shift));
  for (std::size_t j = 0; j < a.fraction.size(); ++j) {
    INFO("cell " << j);
    CHECK(same_bits(a.fraction[j], b.fraction[j]));
    CHECK(same_bits(a.px[j], b.px[j]));
    CHECK(same_bits(a.py[j], b.py[j]));
    CHECK(same_bits(a.amplitude[j], b.amplitude[j]));
    CHECK(same_bits(a.w[j], b.w[j]));
    CHECK(a.status[j] == b.status[j]);
  }
  REQUIRE(a.prt.size() == b.prt.size());
  for (std::size_t j = 0; j < a.prt.size(); ++j) {
    CHECK(same_bits(a.prt[j].real(), b.prt[j].real()));
    CHECK(same_bits(a.prt[j].imag(), b.prt[j].imag()));
  }
  REQUIRE(a.mapping.has_value() == b.mapping.has_value());
  if (a.mapping) {
    for (std::size_t k = 0; k < 4; ++k) CHECK(same_bits(a.mapping->m[k], b.mapping->m[k]));
    CHECK(same_bits(a.mapping->s0[0], b.mapping->s0[0]));
    CHECK(same_bits(a.mapping->s0[1], b.mapping->s0[1]));
    CHECK(a.mapping->centre == b.mapping->centre);
    CHECK(a.mapping->e_x == b.mapping->e_x);
    CHECK(a.mapping->e_y == b.mapping->e_y);
    CHECK(a.mapping->d_c == b.mapping->d_c);
    CHECK(same_bits(a.mapping->delta_max, b.mapping->delta_max));
    CHECK(same_bits(a.mapping->delta_rms, b.mapping->delta_rms));
  }
  CHECK(same_bits(a.arrived_fraction, b.arrived_fraction));
  CHECK(a.losses.launched == b.losses.launched);
  CHECK(a.losses.by_status == b.losses.by_status);
  CHECK(a.losses.worst_surface == b.losses.worst_surface);
  REQUIRE(a.warnings.size() == b.warnings.size());
  for (std::size_t k = 0; k < a.warnings.size(); ++k) {
    CHECK(a.warnings[k].code == b.warnings[k].code);
    CHECK(a.warnings[k].location == b.warnings[k].location);
    CHECK(a.warnings[k].message == b.warnings[k].message);
  }
}

/// The cells with f > 0 as a point list, in cell order.
rtt::trace::PointListPupil traced_cells(const PupilFunction& pf) {
  rtt::trace::PointListPupil list;
  for (std::size_t j = 0; j < pf.fraction.size(); ++j) {
    if (pf.fraction[j] > 0.0) list.points.push_back({pf.px[j], pf.py[j]});
  }
  return list;
}

}  // namespace

TEST_CASE("pupil function of a paraboloid on the axis: W = 0, A = 1, every cell arrives",
          "[pupil][function]") {
  const CompiledSystem cs = compiled(load("m2/paraboloid_stop.rtt.json"));
  const PupilFunction pf = rtt::analysis::pupil_function(cs, kPath, 0, 0);
  const int n = 64;
  REQUIRE(pf.grid == n);
  CHECK(pf.path == kPath);
  CHECK(pf.field == 0);
  CHECK(pf.wavelength == 0);
  CHECK(pf.wavelength_um == cs.wavelengths_um()[0]);
  CHECK(pf.reference_wavelength_um == cs.wavelengths_um()[cs.reference_wavelength()]);
  CHECK(pf.focus_shift == 0.0);
  CHECK(pf.prt.empty());
  // The cells are those of pupil_cells.
  const rtt::analysis::PupilCells cells = rtt::analysis::pupil_cells(n);
  REQUIRE(pf.fraction.size() == cells.fraction.size());
  std::size_t traced = 0;
  for (std::size_t j = 0; j < pf.fraction.size(); ++j) {
    INFO("cell " << j);
    CHECK(same_bits(pf.fraction[j], cells.fraction[j]));
    CHECK(same_bits(pf.px[j], cells.px[j]));
    CHECK(same_bits(pf.py[j], cells.py[j]));
    if (pf.fraction[j] == 0.0) {
      CHECK(pf.status[j] == RayStatus::Vignetted);
      CHECK(pf.amplitude[j] == 0.0);
      CHECK(pf.w[j] == 0.0);
      CHECK(std::isnan(pf.end.pos_x[j]));
      continue;
    }
    ++traced;
    REQUIRE(pf.status[j] == RayStatus::Alive);
    CHECK(std::abs(pf.w[j]) < 1e-6);
    CHECK(std::abs(pf.amplitude[j] - 1.0) <= 1e-15);
  }
  CHECK(pf.arrived_fraction == 1.0);
  CHECK(pf.losses.launched == traced);
  CHECK(pf.losses.count(RayStatus::Alive) == traced);
  CHECK(traced > 3000);  // content guard: about pi/4 of 64^2
}

TEST_CASE("pupil function of a paraboloid: mapping, frame and the nonlinearity delta",
          "[pupil][function][mapping]") {
  const CompiledSystem cs = compiled(load("m2/paraboloid_stop.rtt.json"));
  const PupilFunction pf = rtt::analysis::pupil_function(cs, kPath, 0, 0);
  REQUIRE(pf.mapping);
  const rtt::analysis::PupilMapping& m = *pf.mapping;
  // Frame: d_c = -z after the mirror; the image surface is not rotated, so e_x = +x and
  // e_y = d_c x e_x = -y (ADR 0032, point 4: concave mirror on the axis, M = n'|u'| diag(-1, 1)).
  CHECK((m.d_c - Vec3(0.0, 0.0, -1.0)).norm() <= 1e-15);
  CHECK((m.e_x - Vec3(1.0, 0.0, 0.0)).norm() <= 1e-15);
  CHECK((m.e_y - Vec3(0.0, -1.0, 0.0)).norm() <= 1e-15);
  CHECK(m.centre == pf.sphere.centre);
  CHECK(m.s0[0] == 0.0);  // reference wavelength: exactly 0 by construction
  CHECK(m.s0[1] == 0.0);
  const auto fo = rtt::paraxial::first_order(cs, kPath, 0);
  const auto& last = cs.path(kPath).events.back();
  const double n_image = std::abs(cs.media()[last.medium_after].index[0].real());
  const double r = *fo.entrance_pupil->diameter / 2.0;
  const double nu = n_image * r / *fo.rear_focal_length;  // paraxial n'|u'|
  CHECK(std::abs(m.m[0] + nu) <= 1e-6 * nu);
  CHECK(std::abs(m.m[3] - nu) <= 1e-6 * nu);
  CHECK(std::abs(m.m[1]) <= 1e-12);
  CHECK(std::abs(m.m[2]) <= 1e-12);
  // Expected deltas from the sample points: rho g / (1 + g), g = (rho r / 2f)^2.
  const double f = *fo.rear_focal_length / n_image;
  double expected_max = 0.0;
  double sum_f = 0.0;
  double sum_f_d2 = 0.0;
  for (std::size_t j = 0; j < pf.fraction.size(); ++j) {
    if (pf.fraction[j] == 0.0) continue;
    const double rho = std::hypot(pf.px[j], pf.py[j]);
    const double g = (rho * r / (2.0 * f)) * (rho * r / (2.0 * f));
    const double d = rho * g / (1.0 + g);
    expected_max = std::max(expected_max, d);
    sum_f += pf.fraction[j];
    sum_f_d2 += pf.fraction[j] * d * d;
  }
  const double expected_rms = std::sqrt(sum_f_d2 / sum_f);
  CHECK(std::abs(m.delta_max - expected_max) <= 1e-5 * expected_max);
  CHECK(std::abs(m.delta_rms - expected_rms) <= 1e-5 * expected_rms);
  // About 0.0225 at the rim (F/1.7): above the default threshold 1e-2, so the warning is there.
  CHECK(expected_max > 1e-2);
  CHECK(has(pf.warnings, "pupil.mapping_nonlinear"));
  for (const auto& d : pf.warnings) {
    if (d.code == "pupil.mapping_nonlinear") CHECK(d.location == "/paths/0");
  }
  PupilFunctionOptions loose;
  loose.mapping_warning = 0.05;
  CHECK_FALSE(has(rtt::analysis::pupil_function(cs, kPath, 0, 0, loose).warnings,
                  "pupil.mapping_nonlinear"));
}

TEST_CASE("pupil function: W and status bit for bit those of opd_points on the same points",
          "[pupil][function]") {
  // ADR 0032, point 2; also with a focus shift (point 5) and off the axis.
  const CompiledSystem cs = compiled(load("m2/cooke_triplet.rtt.json"));
  for (const double delta : {0.0, 0.3}) {
    for (const std::uint16_t field : {0, 2}) {
      INFO("Delta = " << delta << ", field " << field);
      PupilFunctionOptions o;
      o.grid = 24;
      o.focus_shift = delta;
      const PupilFunction pf = rtt::analysis::pupil_function(cs, kPath, field, 0, o);
      rtt::analysis::OpdOptions oo;
      oo.focus_shift = delta;
      const auto opd = rtt::analysis::opd_points(cs, kPath, field, 0, traced_cells(pf), oo);
      CHECK(opd.sphere.centre == pf.sphere.centre);
      CHECK(same_bits(opd.sphere.radius, pf.sphere.radius));
      std::size_t k = 0;
      std::size_t arrived = 0;
      for (std::size_t j = 0; j < pf.fraction.size(); ++j) {
        if (pf.fraction[j] == 0.0) continue;
        INFO("cell " << j);
        REQUIRE(k < opd.points.size());
        CHECK(same_bits(pf.w[j], opd.points[k].w));
        CHECK(pf.status[j] == opd.points[k].status);
        if (pf.status[j] == RayStatus::Alive) ++arrived;
        ++k;
      }
      CHECK(k == opd.points.size());
      CHECK(pf.losses.by_status == opd.losses.by_status);
      CHECK(arrived > 300);  // content guard
    }
  }
}

TEST_CASE("refocus gives pupil_function with that focus_shift bit for bit", "[pupil][refocus]") {
  // Delta = -190 moves C' to 10 mm from the exit pupil (z = +100): the shifted sphere has
  // R' = 10 mm and the outer rays miss it (their distance from C' exceeds R'), so status,
  // amplitude, W, arrived_fraction and the losses change with Delta.
  const CompiledSystem cs = compiled(load("m2/paraboloid_stop.rtt.json"));
  PupilFunctionOptions o;
  o.grid = 32;
  o.with_prt = true;
  o.apodization = 0.4;
  const PupilFunction base = rtt::analysis::pupil_function(cs, kPath, 0, 0, o);
  for (const double delta : {0.3, -2.5, -190.0}) {
    INFO("Delta = " << delta);
    PupilFunctionOptions od = o;
    od.focus_shift = delta;
    const PupilFunction direct = rtt::analysis::pupil_function(cs, kPath, 0, 0, od);
    const PupilFunction moved = rtt::analysis::refocus(cs, base, delta);
    check_same(moved, direct);
    if (delta == -190.0) {  // content guard: cells drop out by missing the shifted sphere
      CHECK(direct.arrived_fraction < 0.9);
      CHECK(direct.arrived_fraction > 0.0);
    }
  }
  // And back: refocus to 0 gives the base.
  check_same(rtt::analysis::refocus(cs, rtt::analysis::refocus(cs, base, -190.0), 0.0), base);
  CHECK_THROWS_AS(rtt::analysis::refocus(cs, base, std::numeric_limits<double>::quiet_NaN()),
                  std::invalid_argument);
}

TEST_CASE("pupil function: Gaussian apodisation A = exp(-G rho^2)", "[pupil][function]") {
  // Ansys Zemax OpticStudio User Guide, "Apodization Type": A(rho) = exp(-G rho^2). The weight of
  // every ray is 1 (ideal mirror), so A is the apodisation alone, up to the rounding of sqrt(1).
  const CompiledSystem cs = compiled(load("m2/paraboloid_stop.rtt.json"));
  for (const double g : {0.5, 1.0}) {
    INFO("G = " << g);
    PupilFunctionOptions o;
    o.grid = 16;
    o.apodization = g;
    const PupilFunction pf = rtt::analysis::pupil_function(cs, kPath, 0, 0, o);
    double smallest = 1.0;
    for (std::size_t j = 0; j < pf.fraction.size(); ++j) {
      if (pf.fraction[j] == 0.0) continue;
      const double expected = std::exp(-g * (pf.px[j] * pf.px[j] + pf.py[j] * pf.py[j]));
      CHECK(std::abs(pf.amplitude[j] - expected) <= 4e-16);  // sqrt(weight) = 1 up to u, times A
      smallest = std::min(smallest, pf.amplitude[j]);
    }
    CHECK(smallest < std::exp(-0.8 * g));  // content guard: cells near the rim
  }
}

TEST_CASE("pupil function: with_prt copies the accumulated P of each traced cell",
          "[pupil][function]") {
  const CompiledSystem cs = compiled(load("m2/paraboloid_stop.rtt.json"));
  PupilFunctionOptions o;
  o.grid = 12;
  o.with_prt = true;
  const PupilFunction pf = rtt::analysis::pupil_function(cs, kPath, 0, 0, o);
  REQUIRE(pf.prt.size() == 9 * pf.fraction.size());
  // The same rays traced independently.
  const std::uint16_t field = 0;
  rtt::trace::RayBatch rays =
      rtt::trace::make_rays(cs, kPath, std::span(&field, 1), 0, traced_cells(pf));
  static_cast<void>(rtt::trace::SequentialTracer().trace(cs, kPath, rays));
  std::size_t k = 0;
  bool not_identity = false;
  for (std::size_t j = 0; j < pf.fraction.size(); ++j) {
    for (std::size_t e = 0; e < 9; ++e) {
      const std::complex<double> v = pf.prt[9 * j + e];
      if (pf.fraction[j] == 0.0) {
        CHECK(v == std::complex<double>(0.0, 0.0));
        continue;
      }
      CHECK(same_bits(v.real(),
                      rays.prt(static_cast<int>(e / 3), static_cast<int>(e % 3))[k].real()));
      CHECK(same_bits(v.imag(),
                      rays.prt(static_cast<int>(e / 3), static_cast<int>(e % 3))[k].imag()));
      if (e == 0 && std::abs(v - 1.0) > 1e-3) not_identity = true;
    }
    if (pf.fraction[j] > 0.0) ++k;
  }
  CHECK(k == rays.size());
  CHECK(not_identity);  // content guard: the mirror changes P
}

TEST_CASE("pupil function: s0 is the chief ray of the evaluated wavelength", "[pupil][mapping]") {
  // Cooke triplet, field 2 (20 deg): lateral colour moves the chief ray of wavelength 0 against
  // that of the reference wavelength 1. s0 = n' (d . e_x, d . e_y) with the frame of the
  // reference chief ray, recomputed here from independently traced rays.
  const CompiledSystem cs = compiled(load("m2/cooke_triplet.rtt.json"));
  PupilFunctionOptions o;
  o.grid = 8;
  const PupilFunction at_ref = rtt::analysis::pupil_function(cs, kPath, 2, 1, o);
  REQUIRE(at_ref.mapping);
  CHECK(at_ref.mapping->s0[0] == 0.0);
  CHECK(at_ref.mapping->s0[1] == 0.0);
  const PupilFunction pf = rtt::analysis::pupil_function(cs, kPath, 2, 0, o);
  REQUIRE(pf.mapping);
  const std::uint16_t field = 2;
  auto chief = [&](std::uint16_t wl) {
    rtt::trace::RayBatch rays = rtt::trace::make_rays(cs, kPath, std::span(&field, 1), wl,
                                                      rtt::trace::SinglePupilPoint{0.0, 0.0});
    static_cast<void>(rtt::trace::SequentialTracer().trace(cs, kPath, rays));
    return Vec3(rays.dir_x()[0], rays.dir_y()[0], rays.dir_z()[0]);
  };
  const Vec3 dc = chief(1);
  const std::uint32_t image = cs.path(kPath).events.back().surface;
  const Vec3 xb = cs.surfaces()[image].to_global.apply_vector(Vec3::UnitX());
  const Vec3 ex = (xb - xb.dot(dc) * dc).normalized();
  const Vec3 ey = dc.cross(ex);
  const auto& last = cs.path(kPath).events.back();
  const double n = std::abs(cs.media()[last.medium_after].index[0].real());
  const Vec3 d = chief(0);
  CHECK(std::abs(pf.mapping->s0[0] - n * d.dot(ex)) <= 1e-15);
  CHECK(std::abs(pf.mapping->s0[1] - n * d.dot(ey)) <= 1e-15);
  CHECK(std::abs(pf.mapping->s0[1]) > 1e-7);      // content guard: lateral colour
  CHECK(pf.mapping->d_c == at_ref.mapping->d_c);  // the same frame for every wavelength
}

TEST_CASE("pupil function: an afocal image space has no mapping", "[pupil][mapping]") {
  // The singlet with the object in its front focal plane: the rays of an object point leave
  // parallel, s is the same for all of them, M = 0 (ADR 0032, point 4: |det M| <= 1e-12 n'^2).
  System s = load("m1/singlet_const.rtt.json");
  const auto fo = rtt::paraxial::first_order(compiled(s), kPath, 0);
  s.object.at_infinity = false;
  s.object.distance = Param(-*fo.front_focal_z);
  const CompiledSystem cs = compiled(s);
  PupilFunctionOptions o;
  o.grid = 8;
  const PupilFunction pf = rtt::analysis::pupil_function(cs, kPath, 0, 0, o);
  CHECK_FALSE(pf.mapping.has_value());
  CHECK_FALSE(has(pf.warnings, "pupil.mapping_nonlinear"));
  CHECK(pf.arrived_fraction > 0.9);  // content guard: the cells themselves are there
}

TEST_CASE("pupil function: errors with location", "[pupil][function]") {
  SECTION("a chief ray grazing the image surface leaves the PSF plane without a frame") {
    // The detector turned by 90 deg - 1e-5 deg about y: its local x axis lies within 1.7e-7 rad
    // of d_c = -z, below the limit 1e-6 of ADR 0032, point 4.
    System s = load("m2/paraboloid_stop.rtt.json");
    std::get<Element>(s.root.children[2].value).pose.rotation_deg[1] = Param(90.0 - 1e-5);
    const CompiledSystem cs = compiled(s);
    PupilFunctionOptions o;
    o.grid = 4;
    try {
      static_cast<void>(rtt::analysis::pupil_function(cs, kPath, 0, 0, o));
      FAIL("no AnalysisError");
    } catch (const AnalysisError& e) {
      const std::uint32_t image = cs.path(kPath).events.back().surface;
      CHECK(e.location() == std::optional<std::string>(cs.surfaces()[image].location));
      CHECK_FALSE(e.surface().has_value());
      CHECK_FALSE(e.ray_status().has_value());
    }
  }
  SECTION("a path with crystal modes is rejected before any ray is traced") {
    // The calcite plate of the OPD guard in test_opd.cpp, with an ordinary event.
    System s;
    s.name = "calcite plate behind a stop";
    s.environment.medium = "VACUUM";
    s.wavelengths = {{0.5876, 1.0, true}};
    s.aperture = {rtt::model::SystemApertureType::EntrancePupilDiameter, Param(2.0)};
    s.fields = {rtt::model::FieldType::AngleDeg, {{0.0, 0.0, 1.0}}};
    s.root.name = "root";
    rtt::model::Surface sto;
    sto.id = rtt::model::SurfaceId("STO");
    sto.aperture = rtt::model::CircularAperture{1.0, 0.0};
    rtt::model::Surface s1;
    s1.id = rtt::model::SurfaceId("P.S1");
    rtt::model::Surface s2;
    s2.id = rtt::model::SurfaceId("P.S2");
    s2.pose = rtt::model::Pose::along_z(2.0);
    rtt::model::Surface img;
    img.id = rtt::model::SurfaceId("IMG");
    Element plate{"P",
                  rtt::model::ElementKind::Plate,
                  rtt::model::Pose::along_z(10.0),
                  std::nullopt,
                  {s1, s2}};
    plate.crystal = rtt::model::CrystalMaterial{"CONST:1.6584", "CONST:1.4864"};
    plate.optic_axis = std::array{1.0, 0.0, 1.0};
    s.root.children = {{Element{"stop",
                                rtt::model::ElementKind::Stop,
                                rtt::model::Pose::along_z(0.0),
                                std::nullopt,
                                {sto}}},
                       {plate},
                       {Element{"image",
                                rtt::model::ElementKind::Detector,
                                rtt::model::Pose::along_z(20.0),
                                std::nullopt,
                                {img}}}};
    using rtt::model::EventKind;
    using rtt::model::SurfaceId;
    s.paths = {{"o",
                false,
                {{SurfaceId("STO"), EventKind::Transmit, 0},
                 {SurfaceId("P.S1"), EventKind::Ordinary, 0},
                 {SurfaceId("P.S2"), EventKind::Refract, 0},
                 {SurfaceId("IMG"), EventKind::Transmit, 0}}}};
    const CompiledSystem cs = compiled(s);
    try {
      static_cast<void>(rtt::analysis::pupil_function(cs, kPath, 0, 0));
      FAIL("no AnalysisError");
    } catch (const AnalysisError& e) {
      CHECK(e.location() == std::optional<std::string>("/paths/0"));
      CHECK(std::string(e.what()).find("crystal modes") != std::string::npos);
    }
  }
}

TEST_CASE("pupil function: options out of range are invalid_argument", "[pupil][function]") {
  const CompiledSystem cs = compiled(load("m2/paraboloid_stop.rtt.json"));
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();
  const auto bad = [&](auto change) {
    PupilFunctionOptions o;
    o.grid = 4;
    change(o);
    CHECK_THROWS_AS(rtt::analysis::pupil_function(cs, kPath, 0, 0, o), std::invalid_argument);
  };
  bad([](PupilFunctionOptions& o) { o.grid = 1; });
  bad([](PupilFunctionOptions& o) { o.grid = rtt::analysis::kMaxPupilGrid + 1; });
  bad([](PupilFunctionOptions& o) { o.apodization = -0.1; });
  bad([&](PupilFunctionOptions& o) { o.apodization = nan; });
  bad([&](PupilFunctionOptions& o) { o.apodization = inf; });
  bad([&](PupilFunctionOptions& o) { o.focus_shift = nan; });
  bad([](PupilFunctionOptions& o) { o.lost_warning_fraction = 1.5; });
  bad([](PupilFunctionOptions& o) { o.mapping_warning = -0.1; });
  bad([&](PupilFunctionOptions& o) { o.mapping_warning = nan; });
  bad([](PupilFunctionOptions& o) { o.mapping_warning = 1.5; });
}

TEST_CASE("pupil function: a stop that clips the beam gives the loss warnings",
          "[pupil][function]") {
  // The stop aperture shrunk to 10 mm while the system aperture (EPD 60) asks for 30 mm: about
  // (10/30)^2 = 1/9 of the pupil area arrives.
  System s = load("m2/paraboloid_stop.rtt.json");
  std::get<Element>(s.root.children[0].value).surfaces[0].aperture =
      rtt::model::CircularAperture{10.0, 0.0};
  const CompiledSystem cs = compiled(s);
  PupilFunctionOptions o;
  o.grid = 32;
  const PupilFunction pf = rtt::analysis::pupil_function(cs, kPath, 0, 0, o);
  CHECK(has(pf.warnings, "rays.lost"));
  CHECK(has(pf.warnings, "stop.clips_beam"));
  CHECK(std::abs(pf.arrived_fraction - 1.0 / 9.0) < 0.03);
  for (std::size_t j = 0; j < pf.fraction.size(); ++j) {
    if (pf.fraction[j] > 0.0 && pf.status[j] != RayStatus::Alive) {
      CHECK(pf.amplitude[j] == 0.0);
      CHECK(pf.w[j] == 0.0);
    }
  }
}

TEST_CASE("pupil function: bit for bit at 1, 4 and all threads", "[pupil][function]") {
  const CompiledSystem cs = compiled(load("m2/cooke_triplet.rtt.json"));
  PupilFunctionOptions o;
  o.grid = 32;
  o.with_prt = true;
  o.focus_shift = 0.2;
  const PupilFunction all = rtt::analysis::pupil_function(cs, kPath, 2, 0, o);
  for (const int threads : {1, 4}) {
    INFO(threads << " threads");
    const tbb::global_control limit(tbb::global_control::max_allowed_parallelism,
                                    static_cast<std::size_t>(threads));
    check_same(rtt::analysis::pupil_function(cs, kPath, 2, 0, o), all);
  }
}

TEST_CASE("AnalysisError: a location without a lost ray", "[pupil][errors]") {
  // ADR 0032, point 9 and ADR 0034: the two constructors exclude each other.
  const AnalysisError at_path("no frame", std::string("/paths/3"));
  CHECK(at_path.location() == std::optional<std::string>("/paths/3"));
  CHECK_FALSE(at_path.surface().has_value());
  CHECK_FALSE(at_path.ray_status().has_value());
  CHECK_FALSE(at_path.field().has_value());
  CHECK_FALSE(at_path.wavelength().has_value());
  AnalysisError::LostRay ray;
  ray.location = "/root/children/1/surfaces/0";
  ray.ray_status = RayStatus::Vignetted;
  const AnalysisError lost("lost", ray);
  CHECK(lost.location() == std::optional<std::string>("/root/children/1/surfaces/0"));
  CHECK(lost.ray_status() == std::optional<RayStatus>(RayStatus::Vignetted));
  const AnalysisError plain("plain");
  CHECK_FALSE(plain.location().has_value());
}
