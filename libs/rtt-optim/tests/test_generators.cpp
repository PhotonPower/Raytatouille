// Merit generators rms_spot and rms_wavefront (ADR 0030, point 4 and addendum of #168; tests
// G1-G7) and their place in the merit function and in optimize() (H1-H5).

#include <tbb/global_control.h>

#include <algorithm>
#include <bit>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "rtt/analysis/opd.hpp"
#include "rtt/analysis/spot.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/optimization.hpp"
#include "rtt/model/parameters.hpp"
#include "rtt/optim/generators.hpp"
#include "rtt/optim/merit.hpp"
#include "rtt/optim/optimize.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/sequential.hpp"
#include "rtt/trace/sources.hpp"

using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::material::MaterialLibrary;
using rtt::model::Generator;
using rtt::model::SpotGenerator;
using rtt::model::SpotReference;
using rtt::model::System;
using rtt::model::WavefrontGenerator;
using rtt::optim::GeneratorStats;
using rtt::optim::MeritEvaluation;
using rtt::optim::MeritFunction;
using rtt::optim::OptimError;
using rtt::optim::OptimResult;

namespace {

System load(const std::string& relative) {
  return rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/" + relative);
}

rtt::model::Element& element(System& s, std::size_t child) {
  return std::get<rtt::model::Element>(s.root.children[child].value);
}

/// m1/singlet_const: stop, plano-convex L1 (L1.S1 at z = 5, aperture 12.7 mm), detector; object
/// at infinity, EPD 20 (collimated beam of radius 10 mm up to L1.S1), fields 0, 3.5 and 5 deg,
/// wavelengths 0.4861, 0.5876 (reference), 0.6563 um. `lens` replaces the aperture of L1.S1.
System singlet(std::optional<rtt::model::CircularAperture> lens = std::nullopt) {
  System s = load("m1/singlet_const.rtt.json");
  if (lens) element(s, 1).surfaces[0].aperture = *lens;
  return s;
}

/// m2/paraboloid_stop: paraboloid R = -200, focus at z = -100; one field on axis, one
/// wavelength. `detector_z` moves the detector (defocus).
System paraboloid(double detector_z = -100.0) {
  System s = load("m2/paraboloid_stop.rtt.json");
  s.root.children[2].pose.position[2] = detector_z;
  return s;
}

CompiledSystem compiled(const System& s, std::size_t configuration = 0) {
  return rtt::compile::compile(s, MaterialLibrary{}, configuration);
}

SpotGenerator spot_generator(SpotReference reference = SpotReference::Centroid) {
  SpotGenerator g;
  g.path = "main";
  g.reference = reference;
  return g;
}

WavefrontGenerator wavefront_generator() {
  WavefrontGenerator g;
  g.path = "main";
  return g;
}

std::vector<double> residuals(const Generator& g, const CompiledSystem& cs, GeneratorStats& st) {
  std::vector<double> out(rtt::optim::generator_size(g, cs));
  rtt::optim::generator_residuals(g, cs, out, st);
  return out;
}

bool same(double a, double b) {
  return std::bit_cast<std::uint64_t>(a) == std::bit_cast<std::uint64_t>(b);
}

double largest(const std::vector<double>& r) {
  double m = 0.0;
  for (const double v : r) m = std::max(m, std::abs(v));
  return m;
}

double sum_of_squares(const std::vector<double>& r) {
  double s = 0.0;
  for (const double v : r) s += v * v;
  return s;
}

/// Hit point of one ray on the image surface in its local coordinates; none if it does not
/// arrive Alive there. Traced alone (batch of one), independent of the generator's batches.
std::optional<std::pair<double, double>> hit(
    const CompiledSystem& cs, std::uint16_t field, std::uint16_t wavelength, double px, double py) {
  const PathId path{0};
  rtt::trace::RayBatch rays =
      rtt::trace::make_rays(cs, path, std::span<const std::uint16_t>(&field, 1), wavelength,
                            rtt::trace::SinglePupilPoint{px, py});
  static_cast<void>(rtt::trace::SequentialTracer().trace(cs, path, rays));
  const std::uint32_t image = cs.path(path).events.back().surface;
  if (rays.status()[0] != rtt::trace::RayStatus::Alive || rays.last_surface()[0] != image) {
    return std::nullopt;
  }
  const rtt::math::Vec3 p = cs.surfaces()[image].to_local.apply_point(
      rtt::math::Vec3(rays.pos_x()[0], rays.pos_y()[0], rays.pos_z()[0]));
  return std::pair{p.x(), p.y()};
}

}  // namespace

TEST_CASE("generators G1: sizes", "[optim][generators]") {
  const CompiledSystem cs = compiled(singlet());
  SpotGenerator spot = spot_generator();
  WavefrontGenerator wave = wavefront_generator();
  // None: all 3 fields and all 3 wavelengths; rings 3, arms 6 by default.
  CHECK(rtt::optim::generator_size(spot, cs) == 2 * 3 * 3 * 18);
  CHECK(rtt::optim::generator_size(wave, cs) == 3 * 3 * 18);
  spot.fields = std::vector<std::uint16_t>{0, 2};
  spot.wavelengths = std::vector<std::uint16_t>{1};
  spot.rings = 2;
  spot.arms = 5;
  CHECK(rtt::optim::generator_size(spot, cs) == 2 * 2 * 1 * 10);
  wave.fields = spot.fields;
  wave.rings = 4;
  wave.arms = 1;
  CHECK(rtt::optim::generator_size(wave, cs) == 2 * 3 * 4);

  GeneratorStats st;
  std::vector<double> wrong(rtt::optim::generator_size(spot, cs) + 1);
  CHECK_THROWS_AS(rtt::optim::generator_residuals(spot, cs, wrong, st), std::invalid_argument);
  spot.rings = 0;
  CHECK_THROWS_AS(rtt::optim::generator_size(spot, cs), std::invalid_argument);
}

TEST_CASE("generators G2: a paraboloid at its focus gives residuals 0", "[optim][generators]") {
  // A paraboloid images the axial point at infinity perfectly (stigmatic): spot and wavefront
  // vanish up to rounding. Defocused by 0.5 mm, both are clearly not 0 (content check).
  for (const double z : {-100.0, -100.5}) {
    INFO("detector at z = " << z);
    const CompiledSystem cs = compiled(paraboloid(z));
    for (const Generator& g : {Generator(spot_generator()), Generator(wavefront_generator())}) {
      GeneratorStats st;
      const std::vector<double> r = residuals(g, cs, st);
      CHECK(st.rays_launched == 18);
      CHECK(st.rays_lost == 0);
      CHECK(st.undefined.empty());
      if (z == -100.0) {
        CHECK(largest(r) < 1e-9);
      } else {
        CHECK(largest(r) > 1e-3);
      }
    }
  }
}

TEST_CASE("generators G3: rms_spot against independently traced rays", "[optim][generators]") {
  // Fields {2, 1} with model weights 3 and 1 (W_f = 3/4, 1/4), wavelengths {0, 2} with weights
  // 2 and 1 (W_l = 2/3, 1/3), generator weight 4, GaussPupil{2, 5}. Every ray traced alone.
  System s = singlet();
  s.fields.points[1].weight = 1.0;
  s.fields.points[2].weight = 3.0;
  s.wavelengths[0].weight = 2.0;
  s.wavelengths[2].weight = 1.0;
  const CompiledSystem cs = compiled(s);
  const std::vector<std::uint16_t> fields = {2, 1};
  const std::vector<std::uint16_t> wavelengths = {0, 2};
  const double wf[] = {0.75, 0.25};
  const double wl[] = {2.0 / 3.0, 1.0 / 3.0};
  const double w = 4.0;
  const rtt::trace::GaussPupil gauss{2, 5};
  const std::vector<rtt::trace::PupilPoint> pupil = rtt::trace::pupil_points(gauss);
  const std::vector<double> q = rtt::trace::gauss_pupil_weights(gauss);
  const std::size_t n = pupil.size();

  for (const SpotReference reference : {SpotReference::Centroid, SpotReference::Chief}) {
    INFO("reference " << (reference == SpotReference::Centroid ? "centroid" : "chief"));
    SpotGenerator g = spot_generator(reference);
    g.fields = fields;
    g.wavelengths = wavelengths;
    g.rings = gauss.rings;
    g.arms = gauss.arms;
    g.weight = w;
    GeneratorStats st;
    const std::vector<double> r = residuals(g, cs, st);
    REQUIRE(r.size() == 2 * 2 * 2 * n);
    CHECK(st.rays_launched == 2 * 2 * n);
    CHECK(st.rays_lost == 0);

    double expected_ms = 0.0;
    for (std::size_t f = 0; f < 2; ++f) {
      std::vector<std::pair<double, double>> p;  // wavelength-major, then pupil order
      for (std::size_t l = 0; l < 2; ++l) {
        for (std::size_t k = 0; k < n; ++k) {
          const auto h = hit(cs, fields[f], wavelengths[l], pupil[k].px, pupil[k].py);
          REQUIRE(h);
          p.push_back(*h);
        }
      }
      double rx = 0.0;
      double ry = 0.0;
      if (reference == SpotReference::Centroid) {
        double sw = 0.0;
        for (std::size_t l = 0; l < 2; ++l) {
          for (std::size_t k = 0; k < n; ++k) {
            sw += wl[l] * q[k];
            rx += wl[l] * q[k] * p[l * n + k].first;
            ry += wl[l] * q[k] * p[l * n + k].second;
          }
        }
        rx /= sw;
        ry /= sw;
      } else {
        const auto chief = hit(cs, fields[f], cs.reference_wavelength(), 0.0, 0.0);
        REQUIRE(chief);
        rx = chief->first;
        ry = chief->second;
      }
      for (std::size_t l = 0; l < 2; ++l) {
        for (std::size_t k = 0; k < n; ++k) {
          const double c = std::sqrt(w * wf[f] * wl[l] * q[k]);
          const double dx = p[l * n + k].first - rx;
          const double dy = p[l * n + k].second - ry;
          const std::size_t i = 2 * ((f * 2 + l) * n + k);
          INFO("field " << f << ", wavelength " << l << ", point " << k);
          CHECK(std::abs(r[i] - c * dx) <= 1e-12);
          CHECK(std::abs(r[i + 1] - c * dy) <= 1e-12);
          expected_ms += wf[f] * wl[l] * q[k] * (dx * dx + dy * dy);
        }
      }
    }
    CHECK(std::abs(st.mean_square - expected_ms) <= 1e-12 * expected_ms);
    CHECK(std::abs(sum_of_squares(r) - w * st.mean_square) <= 1e-12 * w * st.mean_square);
    CHECK(largest(r) > 1e-3);  // content: off-axis fields with colour
  }
}

TEST_CASE("generators G3: rms_wavefront against opd_points of single rays", "[optim][generators]") {
  // Same choice and weights as the spot case; W of every point from analysis::opd_points with
  // a single pupil point, W_mean per field and wavelength with the quadrature weights.
  System s = singlet();
  s.fields.points[1].weight = 1.0;
  s.fields.points[2].weight = 3.0;
  s.wavelengths[0].weight = 2.0;
  s.wavelengths[2].weight = 1.0;
  const CompiledSystem cs = compiled(s);
  const std::vector<std::uint16_t> fields = {2, 1};
  const std::vector<std::uint16_t> wavelengths = {0, 2};
  const double wf[] = {0.75, 0.25};
  const double wl[] = {2.0 / 3.0, 1.0 / 3.0};
  const double w = 4.0;
  const rtt::trace::GaussPupil gauss{2, 5};
  const std::vector<rtt::trace::PupilPoint> pupil = rtt::trace::pupil_points(gauss);
  const std::vector<double> q = rtt::trace::gauss_pupil_weights(gauss);
  const std::size_t n = pupil.size();

  WavefrontGenerator g = wavefront_generator();
  g.fields = fields;
  g.wavelengths = wavelengths;
  g.rings = gauss.rings;
  g.arms = gauss.arms;
  g.weight = w;
  GeneratorStats st;
  const std::vector<double> r = residuals(g, cs, st);
  REQUIRE(r.size() == 2 * 2 * n);
  CHECK(st.rays_lost == 0);

  double expected_ms = 0.0;
  for (std::size_t f = 0; f < 2; ++f) {
    for (std::size_t l = 0; l < 2; ++l) {
      std::vector<double> wk;
      double sq = 0.0;
      double mean = 0.0;
      for (std::size_t k = 0; k < n; ++k) {
        const rtt::analysis::OpdPupilPoints one =
            rtt::analysis::opd_points(cs, PathId{0}, fields[f], wavelengths[l],
                                      rtt::trace::SinglePupilPoint{pupil[k].px, pupil[k].py});
        REQUIRE(one.points.size() == 1);
        REQUIRE(one.points[0].status == rtt::trace::RayStatus::Alive);
        wk.push_back(one.points[0].w);
        sq += q[k];
        mean += q[k] * one.points[0].w;
      }
      mean /= sq;
      for (std::size_t k = 0; k < n; ++k) {
        const double c = std::sqrt(w * wf[f] * wl[l] * q[k]);
        const std::size_t i = (f * 2 + l) * n + k;
        INFO("field " << f << ", wavelength " << l << ", point " << k);
        CHECK(std::abs(r[i] - c * (wk[k] - mean)) <= 1e-9);
        expected_ms += wf[f] * wl[l] * q[k] * (wk[k] - mean) * (wk[k] - mean);
      }
    }
  }
  CHECK(std::abs(st.mean_square - expected_ms) <= 1e-9 * expected_ms);
  CHECK(largest(r) > 1e-2);  // content: off axis the wavefront is not flat
}

TEST_CASE("generators G4: mean squares against a reference integration", "[optim][generators]") {
  // Defocused paraboloid (0.5 mm), on axis: the Gaussian quadrature with 8 rings is exact for
  // polynomials in rho^2 up to degree 15, so its mean square is the integral over the pupil up to
  // the series remainder. References: the variance of opd_map on a 401 grid (uniform points in
  // the unit circle) and the RMS of analysis::spot on 200 hexapolar rings (bias about 1 / rings
  // from the outer ring, derived for a spot linear in rho: sum i^3 / N^2 / sum i = (1 + 1/N) / 2
  // against 1/2), each within 1e-2.
  const CompiledSystem cs = compiled(paraboloid(-100.5));
  WavefrontGenerator wave = wavefront_generator();
  wave.rings = 8;
  wave.arms = 8;
  GeneratorStats st;
  static_cast<void>(residuals(wave, cs, st));
  rtt::analysis::OpdOptions opd;
  opd.grid = 401;
  const rtt::analysis::OpdMap map = rtt::analysis::opd_map(cs, PathId{0}, 0, 0, opd);
  INFO("wavefront: generator " << st.mean_square << ", grid " << map.rms * map.rms);
  CHECK(std::abs(st.mean_square - map.rms * map.rms) <= 1e-2 * map.rms * map.rms);
  CHECK(map.rms > 0.1);  // content: clearly defocused

  SpotGenerator spot = spot_generator();
  spot.rings = 8;
  spot.arms = 8;
  static_cast<void>(residuals(spot, cs, st));
  rtt::analysis::SpotOptions so;
  so.sampling = rtt::trace::HexapolarPupil{200};
  const rtt::analysis::SpotDiagram d = rtt::analysis::spot(cs, PathId{0}, 0, 0, so);
  const double rms2 = d.stats.rms_centroid * d.stats.rms_centroid;
  INFO("spot: generator " << st.mean_square << ", hexapolar " << rms2);
  CHECK(std::abs(st.mean_square - rms2) <= 1e-2 * rms2);
  CHECK(d.stats.rms_centroid > 0.05);
}

TEST_CASE("generators G5: lost rays", "[optim][generators]") {
  // L1.S1 with radius 5 mm: of GaussPupil{3, 6} (rho = 0.34, 0.71, 0.94 at heights 3.4, 7.1,
  // 9.4 mm on L1.S1, field 0) only the inner ring arrives; the outer 12 rays give 0.
  SECTION("residuals 0 for the lost rays, the rest from the inner ring") {
    const CompiledSystem cs = compiled(singlet(rtt::model::CircularAperture{5.0, 0.0}));
    SpotGenerator spot = spot_generator();
    spot.fields = std::vector<std::uint16_t>{0};
    spot.wavelengths = std::vector<std::uint16_t>{1};
    WavefrontGenerator wave = wavefront_generator();
    wave.fields = spot.fields;
    wave.wavelengths = spot.wavelengths;
    GeneratorStats st;
    const std::vector<double> rs = residuals(spot, cs, st);
    CHECK(st.rays_launched == 18);
    CHECK(st.rays_lost == 12);
    CHECK(st.undefined.empty());
    for (std::size_t i = 12; i < 36; ++i) CHECK(rs[i] == 0.0);
    CHECK(std::isfinite(st.mean_square));
    const std::vector<double> rw = residuals(wave, cs, st);
    CHECK(st.rays_lost == 12);
    for (std::size_t i = 6; i < 18; ++i) CHECK(rw[i] == 0.0);
    // Content: the inner ring alone still has a spot about its centroid (spherical aberration).
    // Its wavefront on axis is constant over the ring, so W - W_mean is about 0: only finite.
    CHECK(largest(rs) > 1e-4);
    for (const double v : rw) CHECK(std::isfinite(v));
  }
  SECTION("no ray of the field arrives: NaN with a reason") {
    // Radius 1 mm: the chief ray passes, every gauss ring is stopped.
    const CompiledSystem cs = compiled(singlet(rtt::model::CircularAperture{1.0, 0.0}));
    for (const Generator& g :
         {Generator(spot_generator(SpotReference::Centroid)),
          Generator(spot_generator(SpotReference::Chief)), Generator(wavefront_generator())}) {
      GeneratorStats st;
      const std::vector<double> r = residuals(g, cs, st);
      CHECK(std::isnan(r[0]));
      CHECK(std::isnan(st.mean_square));
      CHECK(!st.undefined.empty());
      CHECK(st.rays_lost == st.rays_launched);
    }
  }
  SECTION("the chief ray is lost") {
    // An annulus with inner radius 1 mm stops the chief ray only.
    const CompiledSystem cs = compiled(singlet(rtt::model::CircularAperture{12.7, 1.0}));
    GeneratorStats st;
    const std::vector<double> chief = residuals(spot_generator(SpotReference::Chief), cs, st);
    CHECK(std::isnan(chief[0]));
    CHECK(!st.undefined.empty());
    const std::vector<double> centroid = residuals(spot_generator(SpotReference::Centroid), cs, st);
    CHECK(st.undefined.empty());
    CHECK(st.rays_lost == 0);
    for (const double v : centroid) CHECK(std::isfinite(v));
    // The wavefront needs the chief ray as reference: the analysis error of opd_points.
    CHECK_THROWS_AS(residuals(wavefront_generator(), cs, st), rtt::analysis::AnalysisError);
  }
}

TEST_CASE("generators G6: deterministic, also across thread counts", "[optim][generators]") {
  const CompiledSystem cs = compiled(singlet());
  for (const Generator& g : {Generator(spot_generator()), Generator(wavefront_generator())}) {
    GeneratorStats st;
    std::vector<double> one;
    {
      const tbb::global_control limit(tbb::global_control::max_allowed_parallelism, 1);
      one = residuals(g, cs, st);
    }
    std::vector<double> four;
    {
      const tbb::global_control limit(tbb::global_control::max_allowed_parallelism, 4);
      four = residuals(g, cs, st);
    }
    const std::vector<double> again = residuals(g, cs, st);
    REQUIRE(one.size() == four.size());
    for (std::size_t i = 0; i < one.size(); ++i) {
      CHECK(same(one[i], four[i]));
      CHECK(same(one[i], again[i]));
    }
  }
}

TEST_CASE("generators G7: start checks of the merit function", "[optim][generators]") {
  const MaterialLibrary lib;
  const auto rejected = [&lib](const System& s) {
    try {
      const MeritFunction merit(s, lib, nullptr);
      FAIL("no OptimError");
    } catch (const OptimError& e) {
      REQUIRE(e.diagnostics().size() == 1);
      CHECK(e.diagnostics()[0].code == "merit.operand_unsupported");
      CHECK(e.diagnostics()[0].location == "/optimization/generators/0");
    }
  };
  SECTION("the chosen field weights sum to 0") {
    System s = singlet();
    s.fields.points[1].weight = 0.0;
    SpotGenerator g = spot_generator();
    g.fields = std::vector<std::uint16_t>{1};
    s.optimization.generators = {g};
    rejected(s);
  }
  SECTION("the chosen wavelength weights sum to 0") {
    System s = singlet();
    s.wavelengths[0].weight = 0.0;
    s.wavelengths[2].weight = 0.0;
    WavefrontGenerator g = wavefront_generator();
    g.wavelengths = std::vector<std::uint16_t>{2, 0};
    s.optimization.generators = {g};
    rejected(s);
  }
  SECTION("a path that first_order rejects (crystal)") {
    System s = load("m4/calcite_walkoff.rtt.json");
    SpotGenerator g = spot_generator();
    g.path = "o";
    s.optimization.generators = {g};
    rejected(s);
  }
}

TEST_CASE("merit H1: generator residuals follow the operands", "[optim][generators]") {
  // m5/zoom, configuration "tele" for the generator, an EFL operand in "wide": the residuals of
  // the generator are generator_residuals() on the system compiled in its configuration.
  const MaterialLibrary lib;
  System s = load("m5/zoom.rtt.json");
  rtt::model::FirstOrderOperand efl;
  efl.common = {50.0, 1.0, std::string("wide")};
  efl.quantity = rtt::model::FirstOrderQuantity::Efl;
  efl.path = "main";
  WavefrontGenerator wave = wavefront_generator();
  wave.configuration = "tele";
  wave.weight = 2.0;
  SpotGenerator spot = spot_generator();
  spot.fields = std::vector<std::uint16_t>{1};
  s.optimization.operands = {efl};
  s.optimization.generators = {spot, wave};
  const MeritFunction merit(s, lib, nullptr);

  const CompiledSystem wide = compiled(s, 0);
  const CompiledSystem tele = compiled(s, 1);
  const std::size_t ns = rtt::optim::generator_size(spot, wide);
  const std::size_t nw = rtt::optim::generator_size(wave, tele);
  CHECK(merit.size() == 1 + ns + nw);
  const MeritEvaluation e = merit.evaluate(merit.start());
  REQUIRE(e.residuals.size() == merit.size());
  CHECK(e.values.size() == 1);
  CHECK(e.undefined.size() == 1);
  REQUIRE(e.generators.size() == 2);

  GeneratorStats st;
  const std::vector<double> rs = residuals(spot, wide, st);
  for (std::size_t i = 0; i < ns; ++i) CHECK(same(e.residuals[1 + i], rs[i]));
  CHECK(same(e.generators[0].mean_square, st.mean_square));
  CHECK(e.generators[0].rays_launched == st.rays_launched);
  const std::vector<double> rw = residuals(wave, tele, st);
  for (std::size_t i = 0; i < nw; ++i) CHECK(same(e.residuals[1 + ns + i], rw[i]));
  CHECK(same(e.generators[1].mean_square, st.mean_square));
  CHECK(largest(rw) > 0.0);
  CHECK(e.valid());
}

TEST_CASE("optimize H2: generator table, merit phi and the start check", "[optim][generators]") {
  const MaterialLibrary lib;
  System s = singlet();
  std::get<rtt::model::Conic>(element(s, 1).surfaces[0].shape.base).radius.variable = true;
  rtt::model::FirstOrderOperand efl;
  efl.common = {100.0, 1.0, std::nullopt};
  efl.quantity = rtt::model::FirstOrderQuantity::Efl;
  efl.path = "main";
  SpotGenerator spot = spot_generator();
  spot.fields = std::vector<std::uint16_t>{0};
  spot.wavelengths = std::vector<std::uint16_t>{1};
  spot.weight = 3.0;
  s.optimization.operands = {efl};
  s.optimization.generators = {spot};

  SECTION("table and phi in the final state") {
    const OptimResult r = rtt::optim::optimize(s, lib);
    REQUIRE(r.generators.size() == 1);
    const auto& g = r.generators[0];
    CHECK(g.pointer == "/optimization/generators/0");
    CHECK(g.weight == 3.0);
    CHECK(g.rays_launched == 18);
    CHECK(g.rays_lost == 0);
    GeneratorStats st;
    const std::vector<double> end = residuals(spot, compiled(r.system), st);
    CHECK(std::abs(g.rms - std::sqrt(st.mean_square)) <= 1e-12 * g.rms);
    REQUIRE(r.operands.size() == 1);
    CHECK(std::abs(r.operands[0].contribution + g.contribution - 100.0) <= 1e-9);
    CHECK(g.contribution > 0.0);
    // phi = 2 F / (sum of the operand and generator weights).
    const double efl_end = r.operands[0].value - 100.0;
    const double phi = (efl_end * efl_end + sum_of_squares(end)) / (1.0 + 3.0);
    REQUIRE(!r.history.empty());
    CHECK(std::abs(r.history.back().phi - phi) <= 1e-9 * phi);
    for (const auto& d : r.diagnostics) CHECK(d.code != "optim.rays_lost");
  }
  SECTION("weight 0: the RMS is still reported") {
    std::get<SpotGenerator>(s.optimization.generators[0]).weight = 0.0;
    const OptimResult r = rtt::optim::optimize(s, lib);
    REQUIRE(r.generators.size() == 1);
    CHECK(r.generators[0].rms > 0.0);
    CHECK(r.generators[0].contribution == 0.0);
  }
  SECTION("no ray arrives at the start: the error names the generator") {
    element(s, 1).surfaces[0].aperture = rtt::model::CircularAperture{1.0, 0.0};
    try {
      static_cast<void>(rtt::optim::optimize(s, lib));
      FAIL("no exception");
    } catch (const std::invalid_argument& e) {
      CHECK(std::string(e.what()).find("/optimization/generators/0") != std::string::npos);
    }
  }
}

TEST_CASE("optimize H3: lost rays in the final state give optim.rays_lost", "[optim][generators]") {
  // L1.S1 with radius 5 mm stops 12 of 18 rays for every R1 (collimated beam before L1.S1).
  const MaterialLibrary lib;
  System s = singlet(rtt::model::CircularAperture{5.0, 0.0});
  std::get<rtt::model::Conic>(element(s, 1).surfaces[0].shape.base).radius.variable = true;
  rtt::model::FirstOrderOperand efl;
  efl.common = {100.0, 1.0, std::nullopt};
  efl.quantity = rtt::model::FirstOrderQuantity::Efl;
  efl.path = "main";
  WavefrontGenerator wave = wavefront_generator();
  wave.fields = std::vector<std::uint16_t>{0};
  wave.wavelengths = std::vector<std::uint16_t>{1};
  s.optimization.operands = {efl};
  s.optimization.generators = {wave};
  const OptimResult r = rtt::optim::optimize(s, lib);
  REQUIRE(r.generators.size() == 1);
  CHECK(r.generators[0].rays_lost == 12);
  const auto found = std::find_if(r.diagnostics.begin(), r.diagnostics.end(),
                                  [](const auto& d) { return d.code == "optim.rays_lost"; });
  REQUIRE(found != r.diagnostics.end());
  CHECK(found->location == "/optimization/generators/0");
  CHECK(found->message.find("12") != std::string::npos);
}

TEST_CASE("optimize H4: the spot generator lowers the RMS spot", "[optim][generators]") {
  // Smoke test of the whole chain (the acceptance with analytic values is #170): R1 variable,
  // only a spot generator on axis at the reference wavelength.
  const MaterialLibrary lib;
  System s = singlet();
  std::get<rtt::model::Conic>(element(s, 1).surfaces[0].shape.base).radius.variable = true;
  SpotGenerator spot = spot_generator();
  spot.fields = std::vector<std::uint16_t>{0};
  spot.wavelengths = std::vector<std::uint16_t>{1};
  s.optimization.generators = {spot};
  GeneratorStats start;
  static_cast<void>(residuals(spot, compiled(s), start));
  const OptimResult r = rtt::optim::optimize(s, lib);
  CHECK(r.status != rtt::optim::LmStatus::Failed);
  REQUIRE(r.generators.size() == 1);
  CHECK(r.generators[0].rms < std::sqrt(start.mean_square));
}
