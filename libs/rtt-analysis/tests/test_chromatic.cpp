#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#include "rtt/analysis/chromatic.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
#include "rtt/paraxial/paraxial.hpp"
#include "rtt/trace/sequential.hpp"

using rtt::compile::compile;
using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::material::MaterialLibrary;
using rtt::model::System;

namespace {

System load(const std::string& relative) {
  return rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/" + relative);
}

/// Adds the SCHOTT test catalogue (N-BK7, F2); MaterialLibrary is neither copyable nor movable.
void add_schott(MaterialLibrary& lib) {
  lib.add_catalog(std::string(RTT_CATALOG_DIR) + "/schott.agf");
}

double medium_index(const CompiledSystem& cs, const std::string& reference, std::uint16_t wl) {
  const auto& m = cs.media();
  const auto it =
      std::find_if(m.begin(), m.end(), [&](const auto& x) { return x.reference == reference; });
  REQUIRE(it != m.end());
  return it->index[wl].real();
}

}  // namespace

TEST_CASE("singlet: paraxial longitudinal colour equals the BFL difference", "[colour]") {
  // Issue #31: Delta = BFL(first) - BFL(second) from rtt-paraxial, relative 1e-12.
  MaterialLibrary lib;
  add_schott(lib);
  const CompiledSystem cs = compile(load("m0/singlet.rtt.json"), lib);
  const auto lc = rtt::analysis::longitudinal_colour(cs, PathId{0});
  REQUIRE(lc.pair.first == 0);
  REQUIRE(lc.pair.second == 2);
  const auto f0 = rtt::paraxial::first_order(cs, PathId{0}, 0);
  const auto f2 = rtt::paraxial::first_order(cs, PathId{0}, 2);
  const double expected = *f0.bfl - *f2.bfl;
  REQUIRE(expected != 0.0);
  REQUIRE(std::abs(lc.paraxial - expected) <= 1e-12 * std::abs(expected));
  REQUIRE(lc.foci.size() == 3);
  REQUIRE(lc.foci[0].paraxial_z == *f0.rear_focal_z);
  // Positive lens with normal dispersion: blue focuses first, so F - C < 0.
  REQUIRE(lc.paraxial < 0.0);
}

TEST_CASE("singlet: the real zone focus approaches the paraxial focus", "[colour]") {
  // For a small zone the real axis crossing converges to the paraxial focus (spherical
  // aberration ~ zone^2), so the real colour difference converges to the paraxial one.
  MaterialLibrary lib;
  add_schott(lib);
  const CompiledSystem cs = compile(load("m0/singlet.rtt.json"), lib);
  rtt::analysis::ChromaticOptions options;
  options.zone = 0.01;
  const auto lc = rtt::analysis::longitudinal_colour(cs, PathId{0}, options);
  for (const auto& f : lc.foci) {
    REQUIRE(std::abs(f.real_z - f.paraxial_z) <= 1e-3);
  }
  REQUIRE(std::abs(lc.real - lc.paraxial) <= 1e-3 * std::abs(lc.paraxial));
}

TEST_CASE("achromat N-BK7/F2: longitudinal colour against an independent y-nu trace", "[colour]") {
  // Independent paraxial trace of the cemented doublet in m2/achromat.rtt.json (vertices at
  // z = 5, 11, 14 mm; R = 61.5, -35.2, -136.9 mm), written here without rtt-paraxial:
  // J. E. Greivenkamp, OPTI-201/202 lecture notes, Sec. 9, p. 9-2 (see docs/quellen.md):
  // transfer y' = y + t (n u) / n, refraction n'u' = n u - y phi, phi = c (n' - n). A ray
  // parallel to the axis (y = 1, n u = 0) leaves the last vertex with (y, n'u'); it meets the
  // axis at z_F' = 14 - y n' / (n'u'). Indices from the compiled system, so the test stays valid
  // when the environment medium changes (#25).
  MaterialLibrary lib;
  add_schott(lib);
  const CompiledSystem cs = compile(load("m2/achromat.rtt.json"), lib);
  const auto focus = [&](std::uint16_t wl) {
    const double n0 = cs.media()[cs.environment_medium()].index[wl].real();
    const double n1 = medium_index(cs, "SCHOTT:N-BK7", wl);
    const double n2 = medium_index(cs, "SCHOTT:F2", wl);
    const double z[] = {5.0, 11.0, 14.0};
    const double c[] = {1.0 / 61.5, -1.0 / 35.2, -1.0 / 136.9};
    const double n_before[] = {n0, n1, n2};
    const double n_after[] = {n1, n2, n0};
    double y = 1.0;
    double nu = 0.0;
    for (int i = 0; i < 3; ++i) {
      if (i > 0) y += (z[i] - z[i - 1]) * nu / n_before[i];
      nu -= y * c[i] * (n_after[i] - n_before[i]);
    }
    return 14.0 - y * n0 / nu;
  };
  const auto lc = rtt::analysis::longitudinal_colour(cs, PathId{0});
  for (const auto& f : lc.foci) {
    INFO("wavelength " << f.wavelength);
    REQUIRE(std::abs(f.paraxial_z - focus(f.wavelength)) <= 1e-9);
  }
  REQUIRE(std::abs(lc.paraxial - (focus(0) - focus(2))) <= 1e-9);
  // The achromat corrects F against C far better than the N-BK7 singlet of similar power.
  MaterialLibrary lib2;
  add_schott(lib2);
  const auto singlet =
      rtt::analysis::longitudinal_colour(compile(load("m0/singlet.rtt.json"), lib2), PathId{0});
  REQUIRE(std::abs(lc.paraxial) < 0.1 * std::abs(singlet.paraxial));
}

TEST_CASE("lateral colour: chief ray per wavelength relative to the reference", "[colour]") {
  MaterialLibrary lib;
  add_schott(lib);
  const CompiledSystem cs = compile(load("m0/singlet.rtt.json"), lib);
  const auto lat = rtt::analysis::lateral_colour(cs, PathId{0}, 1);
  REQUIRE(lat.chief.size() == 3);
  const std::uint16_t ref = cs.reference_wavelength();
  REQUIRE(lat.offset[ref].x == 0.0);
  REQUIRE(lat.offset[ref].y == 0.0);
  // Independent: aim and trace each chief ray directly.
  const std::uint32_t image = cs.path(PathId{0}).events.back().surface;
  for (std::uint16_t wl = 0; wl < 3; ++wl) {
    const auto aimed = rtt::trace::aim_ray(cs, PathId{0}, 1, wl, 0.0, 0.0);
    rtt::trace::RayBatch rays(1);
    rays.pos_x()[0] = aimed.ray.pos.x();
    rays.pos_y()[0] = aimed.ray.pos.y();
    rays.pos_z()[0] = aimed.ray.pos.z();
    rays.dir_x()[0] = aimed.ray.dir.x();
    rays.dir_y()[0] = aimed.ray.dir.y();
    rays.dir_z()[0] = aimed.ray.dir.z();
    rays.wl()[0] = wl;
    [[maybe_unused]] const auto stats = rtt::trace::SequentialTracer().trace(cs, PathId{0}, rays);
    const auto p = cs.surfaces()[image].to_local.apply_point(
        rtt::math::Vec3(rays.pos_x()[0], rays.pos_y()[0], rays.pos_z()[0]));
    REQUIRE(lat.chief[wl].y == p.y());
    REQUIRE(lat.offset[wl].y == p.y() - lat.chief[ref].y);
  }
  // Dispersion separates the colours.
  REQUIRE(lat.offset[0].y != 0.0);
  // Constant index in vacuum: no lateral colour at all.
  System s = load("m1/singlet_const.rtt.json");
  s.environment.medium = "VACUUM";
  const MaterialLibrary plain;
  const auto none = rtt::analysis::lateral_colour(compile(s, plain), PathId{0}, 2);
  for (const auto& o : none.offset) {
    REQUIRE(o.x == 0.0);
    REQUIRE(o.y == 0.0);
  }
}

TEST_CASE("invalid colour input", "[colour]") {
  MaterialLibrary lib;
  add_schott(lib);
  const CompiledSystem cs = compile(load("m0/singlet.rtt.json"), lib);
  rtt::analysis::ChromaticOptions options;
  options.pair = rtt::paraxial::ChromaticPair{0, 5};
  REQUIRE_THROWS_AS(rtt::analysis::longitudinal_colour(cs, PathId{0}, options),
                    std::invalid_argument);
  options = {};
  options.zone = 0.0;
  REQUIRE_THROWS_AS(rtt::analysis::longitudinal_colour(cs, PathId{0}, options),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(rtt::analysis::lateral_colour(cs, PathId{0}, 9), std::invalid_argument);
}
