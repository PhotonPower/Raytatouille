// Irradiance map on a surface from traced rays (#379, ADR 0048). Reference values are exact
// (rays at cell centres, the ideal beam splitter R = 0.5) or analytic bounds (binning moves a
// ray by at most half a cell).

#include <tbb/global_control.h>

#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numbers>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "rtt/analysis/irradiance.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/path.hpp"
#include "rtt/trace/ray_paths.hpp"
#include "rtt/trace/run_control.hpp"
#include "rtt/trace/sequential.hpp"
#include "rtt/trace/sources.hpp"

using rtt::analysis::irradiance_map;
using rtt::analysis::IrradianceExtent;
using rtt::analysis::IrradianceMap;
using rtt::analysis::IrradianceOptions;
using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::material::MaterialLibrary;
using rtt::model::SurfaceId;
using rtt::trace::RayBatch;
using rtt::trace::RayPaths;
using rtt::trace::RayStatus;

namespace {

constexpr double kEps = std::numeric_limits<double>::epsilon();

CompiledSystem compiled(const std::string& relative) {
  const MaterialLibrary lib;
  return rtt::compile::compile(
      rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/" + relative), lib);
}

std::size_t path_index(const CompiledSystem& cs, const std::string& name) {
  for (std::size_t i = 0; i < cs.paths().size(); ++i) {
    if (cs.paths()[i].name == name) return i;
  }
  FAIL("no path " << name);
  return 0;
}

/// Collimated rays along +z from z = -1, at the centres of 2 x 2 sub-cells of every cell of the
/// grid [-1, 1]^2 with 4 x 4 cells (pitch 0.25 mm, offset 0.125 mm): 64 rays, 4 per cell.
RayBatch cell_centred_bundle() {
  RayBatch rays(64);
  for (std::size_t i = 0; i < 8; ++i) {
    for (std::size_t j = 0; j < 8; ++j) {
      const std::size_t k = 8 * j + i;
      rays.pos_x()[k] = -1.0 + 0.125 + 0.25 * static_cast<double>(i);
      rays.pos_y()[k] = -1.0 + 0.125 + 0.25 * static_cast<double>(j);
      rays.pos_z()[k] = -1.0;
    }
  }
  return rays;
}

/// Collimated rays along +z from z = 0 on a square grid of pitch 0.5 mm inside r <= 2.25 mm
/// (within the EPD of 5 mm of m0/michelson, which has no stop to aim at), weight 1.
RayBatch michelson_bundle() {
  std::vector<std::pair<double, double>> points;
  for (int i = -4; i <= 4; ++i) {
    for (int j = -4; j <= 4; ++j) {
      const double x = 0.5 * i;
      const double y = 0.5 * j;
      if (x * x + y * y <= 2.25 * 2.25) points.emplace_back(x, y);
    }
  }
  RayBatch rays(points.size());
  for (std::size_t k = 0; k < points.size(); ++k) {
    rays.pos_x()[k] = points[k].first;
    rays.pos_y()[k] = points[k].second;
  }
  return rays;
}

IrradianceOptions grid(std::size_t n, IrradianceExtent extent) {
  IrradianceOptions o;
  o.nx = n;
  o.ny = n;
  o.extent = extent;
  return o;
}

/// Sum over the cells times the cell area.
double integral(const IrradianceMap& m) {
  const double dx = m.x_edges[1] - m.x_edges[0];
  const double dy = m.y_edges[1] - m.y_edges[0];
  double sum = 0.0;
  for (const double v : m.values) sum += v;
  return sum * dx * dy;
}

}  // namespace

TEST_CASE("irradiance: a uniform collimated bundle gives a constant map (#379)",
          "[irradiance][379]") {
  // m4/ghost_plates: two plane-parallel plates in front of the detector IMG. At normal incidence
  // the rays keep x and y; every ray loses the same Fresnel share, so all reach IMG with one
  // weight w. Four rays per cell of 0.5 x 0.5 mm: value = 4 w / 0.25 = 16 w, exactly.
  const CompiledSystem cs = compiled("m4/ghost_plates.rtt.json");
  RayBatch rays = cell_centred_bundle();
  RayPaths paths;
  const rtt::trace::SequentialTracer tracer;
  static_cast<void>(tracer.trace(cs, PathId{0}, rays, paths));
  const double w = rays.weight()[0];
  REQUIRE(w > 0.5);
  REQUIRE(w < 1.0);  // content guard: the plates take a Fresnel share
  for (std::size_t k = 0; k < rays.size(); ++k) {
    REQUIRE(rays.status()[k] == RayStatus::Alive);
    REQUIRE(rays.weight()[k] == w);
  }
  const IrradianceOptions o = grid(4, {-1.0, 1.0, -1.0, 1.0});
  for (const IrradianceMap& m :
       {irradiance_map(cs, 0, paths, SurfaceId{"IMG"}, o), irradiance_map(cs, 0, rays, o)}) {
    REQUIRE(m.values.size() == 16);
    for (const double v : m.values) CHECK(v == 16.0 * w);
    CHECK(m.rays_counted == 64);
    // 64 additions of w: rounding of the running sum, at most 64 ulp.
    CHECK(std::abs(m.weight_counted - 64.0 * w) <= 64.0 * kEps * 64.0 * w);
    CHECK(m.x_edges == std::vector<double>{-1.0, -0.5, 0.0, 0.5, 1.0});
    CHECK(m.y_edges == m.x_edges);
    CHECK(m.surface_id.str() == "IMG");
    CHECK(m.vignetted.count == 0);
    CHECK(m.not_reached.count == 0);
    CHECK(m.outside_grid.count == 0);
    // The mean is the power over the area: 64 w / 4 mm^2.
    CHECK(std::abs(integral(m) / 4.0 - 16.0 * w) <= 16.0 * kEps * 16.0 * w);
  }
}

TEST_CASE("irradiance: energy is conserved and the map is centred on the spot (#379)",
          "[irradiance][379]") {
  // m1/singlet_const at its image IMG, on axis and at the field 3.5 degree, 21 x 21 grid pupil.
  const CompiledSystem cs = compiled("m1/singlet_const.rtt.json");
  for (const std::uint16_t field : {std::uint16_t{0}, std::uint16_t{1}}) {
    INFO("field " << field);
    const std::vector<std::uint16_t> fields{field};
    RayBatch rays = rtt::trace::make_rays(cs, PathId{0}, fields, 0, rtt::trace::GridPupil{21});
    const rtt::trace::SequentialTracer tracer;
    static_cast<void>(tracer.trace(cs, PathId{0}, rays));
    // Ray centroid on IMG in local coordinates (IMG lies at z = const, axes parallel to global).
    const auto& img = cs.surfaces()[*cs.find_surface(SurfaceId{"IMG"})];
    double sw = 0.0;
    double sx = 0.0;
    double sy = 0.0;
    std::size_t alive = 0;
    for (std::size_t k = 0; k < rays.size(); ++k) {
      if (rays.status()[k] != RayStatus::Alive) continue;
      const rtt::math::Vec3 p =
          img.to_local.apply_point({rays.pos_x()[k], rays.pos_y()[k], rays.pos_z()[k]});
      sw += rays.weight()[k];
      sx += rays.weight()[k] * p.x();
      sy += rays.weight()[k] * p.y();
      ++alive;
    }
    REQUIRE(alive > 300);  // content guard
    const double cx = sx / sw;
    const double cy = sy / sw;
    // A grid of 0.2 mm around the centroid that holds the whole spot.
    const IrradianceMap m =
        irradiance_map(cs, 0, rays, grid(40, {cx - 4.0, cx + 4.0, cy - 4.0, cy + 4.0}));
    CHECK(m.rays_counted == alive);
    CHECK(m.outside_grid.count == 0);
    CHECK(m.rays_counted + m.vignetted.count + m.not_reached.count == rays.size());
    // Energy: the integral of the map is the counted weight, up to rounding of the sums.
    CHECK(std::abs(integral(m) - m.weight_counted) <= 1e-12 * m.weight_counted);
    CHECK(std::abs(m.weight_counted - sw) <= 1e-12 * sw);
    // Centroid of the map (cell centres): binning moves each ray by at most half a cell, so the
    // centroid moves by at most half a cell in each axis.
    const double d = m.x_edges[1] - m.x_edges[0];
    double mw = 0.0;
    double mx = 0.0;
    double my = 0.0;
    for (std::size_t iy = 0; iy < m.ny; ++iy) {
      for (std::size_t ix = 0; ix < m.nx; ++ix) {
        const double v = m.values[iy * m.nx + ix];
        mw += v;
        mx += v * 0.5 * (m.x_edges[ix] + m.x_edges[ix + 1]);
        my += v * 0.5 * (m.y_edges[iy] + m.y_edges[iy + 1]);
      }
    }
    CHECK(std::abs(mx / mw - cx) <= 0.5 * d + 1e-12);
    CHECK(std::abs(my / mw - cy) <= 0.5 * d + 1e-12);
    if (field == 0) {
      // On axis the chief ray meets IMG at the origin; by symmetry so does the centroid.
      CHECK(std::abs(cx) <= 1e-12);
      CHECK(std::abs(cy) <= 1e-12);
    }
  }
}

TEST_CASE("irradiance: independent of threads and block size; both overloads agree (#379)",
          "[irradiance][379]") {
  const CompiledSystem cs = compiled("m1/singlet_const.rtt.json");
  const std::vector<std::uint16_t> fields{1};
  const RayBatch start = rtt::trace::make_rays(cs, PathId{0}, fields, 0, rtt::trace::GridPupil{15});
  const IrradianceOptions o = grid(16, {-3.0, 3.0, 4.0, 10.0});
  const rtt::trace::SequentialTracer tracer;
  RayBatch plain = start;
  RayPaths recorded;
  static_cast<void>(tracer.trace(cs, PathId{0}, plain, recorded));
  const IrradianceMap reference = irradiance_map(cs, 0, plain, o);
  REQUIRE(reference.rays_counted > 100);  // content guard: the grid holds the spot
  const IrradianceMap from_paths = irradiance_map(cs, 0, recorded, SurfaceId{"IMG"}, o);
  CHECK(from_paths.values == reference.values);  // bitwise
  CHECK(from_paths.weight_counted == reference.weight_counted);
  for (const std::size_t block : {std::size_t{1}, std::size_t{7}, std::size_t{1000}}) {
    INFO("block " << block);
    rtt::trace::RunControl control;
    control.cancel = rtt::trace::CancelToken{};
    control.block_size = block;
    RayBatch rays = start;
    static_cast<void>(tracer.trace(cs, PathId{0}, rays, control));
    const IrradianceMap m = irradiance_map(cs, 0, rays, o);
    CHECK(m.values == reference.values);
    CHECK(m.weight_counted == reference.weight_counted);
  }
}

TEST_CASE("irradiance: the incident weight counts, and the event can be chosen (#379)",
          "[irradiance][379]") {
  // m0/michelson, path "reference arm": BS (transmit), M1, BS (reflect), DET. The ideal beam
  // splitter has R = T = 0.5 and the mirror is lossless, so the rays meet BS first with weight 1
  // and the second time with 0.5. The incident weight counts, not the weight after the event.
  const CompiledSystem cs = compiled("m0/michelson.rtt.json");
  const std::size_t arm = path_index(cs, "reference arm");
  RayBatch rays = michelson_bundle();
  const std::size_t n = rays.size();
  REQUIRE(n > 50);  // content guard
  RayPaths paths;
  const rtt::trace::SequentialTracer tracer;
  static_cast<void>(tracer.trace(cs, PathId{static_cast<std::uint32_t>(arm)}, rays, paths));
  for (std::size_t k = 0; k < n; ++k) REQUIRE(rays.status()[k] == RayStatus::Alive);
  // The weights are products of sqrt(R) amplitudes squared, 0.5 only up to rounding (as in
  // test_paths.cpp): 1e-12 relative.
  const double nd = static_cast<double>(n);
  // A grid that holds the beam on BS (EPD of the file, in the local frame of the tilted BS).
  const IrradianceOptions all = grid(8, {-20.0, 20.0, -20.0, 20.0});
  IrradianceOptions first = all;
  first.event = 0;
  IrradianceOptions second = all;
  second.event = 2;
  const IrradianceMap m0 = irradiance_map(cs, arm, paths, SurfaceId{"BS"}, first);
  const IrradianceMap m2 = irradiance_map(cs, arm, paths, SurfaceId{"BS"}, second);
  const IrradianceMap sum = irradiance_map(cs, arm, paths, SurfaceId{"BS"}, all);
  CHECK(m0.events == std::vector<std::size_t>{0});
  CHECK(m2.events == std::vector<std::size_t>{2});
  CHECK(sum.events == std::vector<std::size_t>{0, 2});
  CHECK(m0.rays_counted == n);
  CHECK(m0.weight_counted == nd);  // incident weight 1 exactly, not T = 0.5 after BS
  CHECK(std::abs(m2.weight_counted - 0.5 * nd) <= 1e-12 * nd);
  CHECK(sum.rays_counted == 2 * n);
  CHECK(std::abs(sum.weight_counted - 1.5 * nd) <= 1e-12 * nd);
  for (std::size_t i = 0; i < sum.values.size(); ++i) {
    CHECK(std::abs(sum.values[i] - (m0.values[i] + m2.values[i])) <= 1e-12 * sum.values[i]);
  }
  // Event 1 is M1, not BS.
  IrradianceOptions wrong = all;
  wrong.event = 1;
  CHECK_THROWS_AS(irradiance_map(cs, arm, paths, SurfaceId{"BS"}, wrong), std::invalid_argument);
}

TEST_CASE("irradiance: input errors (#379)", "[irradiance][379]") {
  const CompiledSystem cs = compiled("m4/ghost_plates.rtt.json");
  RayBatch rays = cell_centred_bundle();
  RayPaths paths;
  const rtt::trace::SequentialTracer tracer;
  static_cast<void>(tracer.trace(cs, PathId{0}, rays, paths));
  const IrradianceOptions ok = grid(4, {-1.0, 1.0, -1.0, 1.0});
  SECTION("no extent and no aperture on the surface") {
    REQUIRE_FALSE(cs.surfaces()[*cs.find_surface(SurfaceId{"IMG"})].aperture.has_value());
    CHECK_THROWS_AS(irradiance_map(cs, 0, paths, SurfaceId{"IMG"}), std::invalid_argument);
    CHECK_THROWS_AS(irradiance_map(cs, 0, rays), std::invalid_argument);
  }
  SECTION("grid") {
    IrradianceOptions zero = ok;
    zero.nx = 0;
    CHECK_THROWS_AS(irradiance_map(cs, 0, rays, zero), std::invalid_argument);
    CHECK_THROWS_AS(irradiance_map(cs, 0, rays, grid(4, {1.0, 1.0, -1.0, 1.0})),
                    std::invalid_argument);
    const double nan = std::numeric_limits<double>::quiet_NaN();
    CHECK_THROWS_AS(irradiance_map(cs, 0, rays, grid(4, {-1.0, 1.0, nan, 1.0})),
                    std::invalid_argument);
  }
  SECTION("surface") {
    CHECK_THROWS_AS(irradiance_map(cs, 0, paths, SurfaceId{"NOPE"}, ok), std::invalid_argument);
    CHECK_THROWS_AS(irradiance_map(cs, 7, rays, ok), std::invalid_argument);
  }
  SECTION("the default extent is the bounding rectangle of the aperture") {
    // The stop STO has a circular aperture: [-r, r]^2.
    const auto& a = cs.surfaces()[*cs.find_surface(SurfaceId{"STO"})];
    REQUIRE(a.aperture.has_value());
    const double r = std::get<rtt::model::CircularAperture>(*a.aperture).radius;
    const IrradianceMap m = irradiance_map(cs, 0, paths, SurfaceId{"STO"});
    CHECK(m.x_edges.front() == -r);
    CHECK(m.x_edges.back() == r);
    CHECK(m.y_edges.front() == -r);
    CHECK(m.y_edges.back() == r);
    CHECK(m.rays_counted == 64);
    CHECK(m.weight_counted == 64.0);  // incident on the first surface: weight 1
  }
}

TEST_CASE("irradiance: the batch overload refuses a last event that changes the weight (#379)",
          "[irradiance][379]") {
  // m0/michelson with an extra path that ends at the beam splitter's Transmit (T = 0.5): the
  // batch then holds the weight after the event, not the incident one. The two arms end at the
  // detector DET (Transmit, Fresnel, the same medium on both sides) and are accepted.
  rtt::model::System s =
      rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m0/michelson.rtt.json");
  rtt::model::Path to_bs;
  to_bs.name = "to BS";
  to_bs.events = {{rtt::model::SurfaceId{"BS"}, rtt::model::EventKind::Transmit}};
  s.paths.push_back(to_bs);
  const MaterialLibrary lib;
  const CompiledSystem cs = rtt::compile::compile(s, lib);
  const IrradianceOptions o = grid(4, {-20.0, 20.0, -20.0, 20.0});
  for (std::size_t p = 0; p < cs.paths().size(); ++p) {
    const auto& surface = cs.surfaces()[cs.paths()[p].events.back().surface];
    INFO(cs.paths()[p].name << " ends at " << surface.id.str());
    RayBatch rays = michelson_bundle();
    const rtt::trace::SequentialTracer tracer;
    static_cast<void>(tracer.trace(cs, PathId{static_cast<std::uint32_t>(p)}, rays));
    if (surface.id.str() == "DET") {
      CHECK_NOTHROW(irradiance_map(cs, p, rays, o));
    } else {
      CHECK(cs.paths()[p].name == "to BS");  // content guard: the refused path is ours
      CHECK_THROWS_AS(irradiance_map(cs, p, rays, o), std::invalid_argument);
    }
  }
}

TEST_CASE("irradiance: independent of the number of threads (#379)", "[irradiance][379]") {
  const CompiledSystem cs = compiled("m1/singlet_const.rtt.json");
  const std::vector<std::uint16_t> fields{1};
  const IrradianceOptions o = grid(16, {-3.0, 3.0, 4.0, 10.0});
  const auto run = [&](std::size_t threads) {
    const tbb::global_control limit(tbb::global_control::max_allowed_parallelism, threads);
    RayBatch rays = rtt::trace::make_rays(cs, PathId{0}, fields, 0, rtt::trace::GridPupil{15});
    RayPaths paths;
    const rtt::trace::SequentialTracer tracer;
    static_cast<void>(tracer.trace(cs, PathId{0}, rays, paths));
    return std::pair{irradiance_map(cs, 0, rays, o),
                     irradiance_map(cs, 0, paths, SurfaceId{"L1.S2"}, o)};
  };
  const auto one = run(1);
  REQUIRE(one.first.rays_counted > 100);  // content guard
  for (const std::size_t threads : {std::size_t{2}, std::size_t{4}, std::size_t{8}}) {
    INFO("threads " << threads);
    const auto other = run(threads);
    CHECK(other.first.values == one.first.values);  // bitwise
    CHECK(other.second.values == one.second.values);
  }
}

TEST_CASE("irradiance: the incident weight includes the absorption on the last path (#379)",
          "[irradiance][379]") {
  // m4/ghost_plates with plate A (5 mm, n = 1.5) absorbing: kappa = 1e-5. At normal incidence
  // the rays cross A from A.S1 to A.S2 over t = 5 mm, so the power that reaches A.S2 is the
  // weight after A.S1 times exp(-4 pi kappa t / lambda) (Byrnes, arXiv:1603.02720v5, Eqs. (1),
  // (2)), lambda = 0.5876 um. The RayPaths slot before A.S2 holds the weight at A.S1, without it.
  rtt::model::System s =
      rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m4/ghost_plates.rtt.json");
  std::get<rtt::model::Element>(s.root.children[1].value).material = "CONST:1.5,1e-5";
  const MaterialLibrary lib;
  const CompiledSystem cs = rtt::compile::compile(s, lib);
  RayBatch rays = cell_centred_bundle();
  RayPaths paths;
  const rtt::trace::SequentialTracer tracer;
  static_cast<void>(tracer.trace(cs, PathId{0}, rays, paths));
  REQUIRE(paths.event_surfaces.size() >= 3);
  REQUIRE(cs.surfaces()[paths.event_surfaces[2]].id.str() == "A.S2");
  const double lambda_mm = 0.5876e-3;
  const double factor = std::exp(-4.0 * std::numbers::pi * 1e-5 * 5.0 / lambda_mm);
  REQUIRE(factor < 0.4);  // content guard: a strong, visible absorption
  double expected = 0.0;
  for (std::size_t r = 0; r < paths.ray_count(); ++r) {
    expected += paths.weight[r * paths.slots + 2] * factor;  // after A.S1, then the 5 mm
  }
  IrradianceOptions o = grid(4, {-1.0, 1.0, -1.0, 1.0});
  SECTION("without the wavelength: an input error") {
    CHECK_THROWS_AS(irradiance_map(cs, 0, paths, SurfaceId{"A.S2"}, o), std::invalid_argument);
  }
  SECTION("with the wavelength") {
    o.wavelength = 0;
    const IrradianceMap m = irradiance_map(cs, 0, paths, SurfaceId{"A.S2"}, o);
    CHECK(m.rays_counted == 64);
    CHECK(std::abs(m.weight_counted - expected) <= 1e-12 * expected);
    // The medium before A.S1 (air) does not absorb: no wavelength needed there.
    o.wavelength.reset();
    CHECK(irradiance_map(cs, 0, paths, SurfaceId{"A.S1"}, o).weight_counted == 64.0);
  }
}

TEST_CASE("irradiance: a point on a cell edge falls into the cell to its right (#379)",
          "[irradiance][379]") {
  // Edges x_i = x_min + i (x_max - x_min) / n over [-4.5, -1.25] with n = 23: the index from the
  // formula floor((x_i - x_min) / (x_max - x_min) n) is i - 1 for 17 of the 22 inner edges
  // (rounding; found by a search in #379). Rays exactly on every inner edge must still land in
  // the cell that starts at their edge, as documented.
  const CompiledSystem cs = compiled("m4/ghost_plates.rtt.json");
  IrradianceOptions o;
  o.nx = 23;
  o.ny = 1;
  o.extent = IrradianceExtent{-4.5, -1.25, -1.0, 1.0};
  const IrradianceMap empty = irradiance_map(cs, 0, cell_centred_bundle(), o);
  constexpr std::size_t kInner = 22;
  std::size_t formula_misses = 0;
  for (std::size_t i = 1; i <= kInner; ++i) {
    const double t = (empty.x_edges[i] + 4.5) / 3.25 * 23.0;
    if (std::floor(t) < static_cast<double>(i)) ++formula_misses;
  }
  REQUIRE(formula_misses >= 10);  // content guard: the plain formula fails at many edges
  RayBatch rays(kInner);
  for (std::size_t i = 0; i < kInner; ++i) {
    rays.pos_x()[i] = empty.x_edges[i + 1];
    rays.pos_z()[i] = -1.0;
  }
  const rtt::trace::SequentialTracer tracer;
  static_cast<void>(tracer.trace(cs, PathId{0}, rays));
  for (std::size_t i = 0; i < kInner; ++i) REQUIRE(rays.pos_x()[i] == empty.x_edges[i + 1]);
  const IrradianceMap m = irradiance_map(cs, 0, rays, o);
  REQUIRE(m.rays_counted == kInner);
  for (std::size_t ix = 0; ix < 23; ++ix) {
    INFO("cell " << ix);
    CHECK((m.values[ix] > 0.0) == (ix >= 1));  // cells 1 .. 22 hold one ray each, cell 0 none
  }
}

TEST_CASE("irradiance: an extent too wide for finite cells is refused (#379)",
          "[irradiance][379]") {
  const CompiledSystem cs = compiled("m4/ghost_plates.rtt.json");
  CHECK_THROWS_AS(irradiance_map(cs, 0, cell_centred_bundle(), grid(4, {-1e308, 1e308, -1.0, 1.0})),
                  std::invalid_argument);
}
