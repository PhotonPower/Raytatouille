#include <oneapi/tbb/task_arena.h>

#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#include "rtt/analysis/spot.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"

using rtt::analysis::AnalysisError;
using rtt::analysis::Point2;
using rtt::analysis::SpotPoint;
using rtt::compile::compile;
using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::material::MaterialLibrary;
using rtt::model::Element;
using rtt::model::Param;
using rtt::model::System;

// Spot statistics are definitions (decided for #28): weighted centroid c = sum w r / sum w,
// RMS = sqrt(sum w |r - ref|^2 / sum w), GEO = max |r - ref|. Reference cases of issue #28.

namespace {

System load(const std::string& relative) {
  return rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/" + relative);
}

Element& element(System& s, std::size_t i) {
  return std::get<Element>(s.root.children[i].value);
}

}  // namespace

TEST_CASE("spot statistics: weighted centroid, RMS and GEO by hand", "[spot]") {
  // A = (1, 0) with weight 1, B = (-1, 2) with weight 3, chief at the origin:
  // c = ((1 - 3) / 4, (0 + 6) / 4) = (-0.5, 1.5)
  // |A - c|^2 = 4.5, |B - c|^2 = 0.5 -> RMS_c = sqrt((4.5 + 3 * 0.5) / 4) = sqrt(1.5)
  // |A|^2 = 1, |B|^2 = 5 -> RMS_chief = sqrt((1 + 15) / 4) = 2
  // GEO_c = sqrt(4.5), GEO_chief = sqrt(5)
  const std::vector<SpotPoint> points{{1.0, 0.0, 0, 1.0}, {-1.0, 2.0, 1, 3.0}};
  const auto s = rtt::analysis::spot_statistics(points, Point2{0.0, 0.0});
  REQUIRE(std::abs(s.centroid.x + 0.5) <= 1e-15);
  REQUIRE(std::abs(s.centroid.y - 1.5) <= 1e-15);
  REQUIRE(std::abs(s.rms_centroid - std::sqrt(1.5)) <= 1e-15);
  REQUIRE(std::abs(s.rms_chief - 2.0) <= 1e-15);
  REQUIRE(std::abs(s.geo_centroid - std::sqrt(4.5)) <= 1e-15);
  REQUIRE(std::abs(s.geo_chief - std::sqrt(5.0)) <= 1e-15);
  // Scaling all weights changes nothing.
  const std::vector<SpotPoint> scaled{{1.0, 0.0, 0, 0.25}, {-1.0, 2.0, 1, 0.75}};
  const auto t = rtt::analysis::spot_statistics(scaled, Point2{0.0, 0.0});
  REQUIRE(std::abs(t.rms_centroid - s.rms_centroid) <= 1e-15);
  REQUIRE_THROWS_AS(rtt::analysis::spot_statistics({}, Point2{}), AnalysisError);
  // A point with weight 0 does not count for GEO either.
  const std::vector<SpotPoint> with_zero{
      {1.0, 0.0, 0, 1.0}, {-1.0, 2.0, 1, 3.0}, {50.0, 50.0, 2, 0.0}};
  const auto u = rtt::analysis::spot_statistics(with_zero, Point2{0.0, 0.0});
  REQUIRE(u.geo_chief == s.geo_chief);
  REQUIRE(u.rms_chief == s.rms_chief);
  const std::vector<SpotPoint> zero{{1.0, 0.0, 0, 0.0}};
  REQUIRE_THROWS_AS(rtt::analysis::spot_statistics(zero, Point2{}), AnalysisError);
}

TEST_CASE("paraboloid on axis: RMS spot below 1e-9 mm", "[spot]") {
  // Architecture, Validierung: paraboloid, object at infinity, on axis: RMS spot < 1e-9 mm.
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(load("m2/paraboloid_stop.rtt.json"), lib);
  const auto spot = rtt::analysis::spot(cs, PathId{0}, 0, 0);
  REQUIRE(spot.rays_launched == 127);  // hexapolar, 6 rings
  REQUIRE(spot.rays_arrived == 127);
  REQUIRE(spot.vignetted_fraction == 0.0);
  REQUIRE(spot.image_surface == 2);
  REQUIRE(spot.stats.rms_centroid < 1e-9);
  REQUIRE(spot.stats.rms_chief < 1e-9);
  REQUIRE(spot.stats.geo_chief < 1e-9);
  REQUIRE(std::hypot(spot.chief.x, spot.chief.y) < 1e-9);
}

TEST_CASE("vignetted fraction counts the rays that do not reach the image", "[spot]") {
  // Mirror aperture 21 mm with EPD 60: hexapolar rings k = 5, 6 (radius 25, 30 mm) are
  // vignetted, i.e. 30 + 36 = 66 of 127 rays.
  System s = load("m2/paraboloid_stop.rtt.json");
  element(s, 1).surfaces[0].aperture = rtt::model::CircularAperture{21.0, 0.0};
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(s, lib);
  const auto spot = rtt::analysis::spot(cs, PathId{0}, 0, 0);
  REQUIRE(spot.rays_launched == 127);
  REQUIRE(spot.rays_arrived == 61);
  REQUIRE(spot.points.size() == 61);
  REQUIRE(std::abs(spot.vignetted_fraction - 66.0 / 127.0) <= 1e-15);
}

TEST_CASE("polychromatic spot uses the normalised wavelength weights", "[spot]") {
  // Constant index: all wavelengths give the same spot, so the polychromatic statistics equal
  // the monochromatic ones, while every point carries w_lambda / sum(w) = 1/4 or 3/4. The
  // weighted centroid itself is checked by hand in "spot statistics".
  // Constant index also for the surroundings: the file's AIR (Ciddor since #25) is dispersive,
  // so VACUUM is used here.
  System s = load("m1/singlet_const.rtt.json");
  s.environment.medium = "VACUUM";
  s.wavelengths = {{0.4861, 1.0, false}, {0.5876, 3.0, true}};
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(s, lib);
  const auto mono = rtt::analysis::spot(cs, PathId{0}, 1, 1);
  const auto poly = rtt::analysis::spot(cs, PathId{0}, 1, std::nullopt);
  REQUIRE(!poly.wavelength.has_value());
  REQUIRE(poly.rays_launched == 2 * mono.rays_launched);
  REQUIRE(poly.points.size() == 2 * mono.points.size());
  double sum = 0.0;
  for (const auto& p : poly.points) {
    REQUIRE(p.weight == (p.wavelength == 0 ? 0.25 : 0.75));
    sum += p.weight;
  }
  REQUIRE(std::abs(sum - static_cast<double>(mono.points.size())) <= 1e-12);
  REQUIRE(std::abs(poly.stats.centroid.y - mono.stats.centroid.y) <= 1e-12);
  REQUIRE(std::abs(poly.stats.rms_centroid - mono.stats.rms_centroid) <= 1e-12);
  REQUIRE(poly.chief.y == mono.chief.y);
}

TEST_CASE("spot coordinates are local to the image surface", "[spot]") {
  // Rotating the image surface by +90 deg about z (still rotationally symmetric for the
  // paraxial pupil) maps local x to global y and local y to global -x.
  System s = load("m1/singlet_const.rtt.json");
  const MaterialLibrary lib;
  const CompiledSystem plain = compile(s, lib);
  element(s, 2).pose.rotation_deg = {Param(0.0), Param(0.0), Param(90.0)};
  const CompiledSystem rotated = compile(s, lib);
  const auto a = rtt::analysis::spot(plain, PathId{0}, 2, 1);
  const auto b = rtt::analysis::spot(rotated, PathId{0}, 2, 1);
  REQUIRE(std::abs(b.chief.x - a.chief.y) <= 1e-12);
  REQUIRE(std::abs(b.chief.y + a.chief.x) <= 1e-12);
  REQUIRE(std::abs(b.stats.centroid.x - a.stats.centroid.y) <= 1e-12);
  REQUIRE(std::abs(b.stats.rms_centroid - a.stats.rms_centroid) <= 1e-12);
}

TEST_CASE("spot is bitwise identical with one thread and with all threads", "[spot]") {
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(load("m1/singlet_const.rtt.json"), lib);
  oneapi::tbb::task_arena one(1);
  const auto serial =
      one.execute([&] { return rtt::analysis::spot(cs, PathId{0}, 2, std::nullopt); });
  const auto parallel = rtt::analysis::spot(cs, PathId{0}, 2, std::nullopt);
  REQUIRE(serial.points.size() == parallel.points.size());
  for (std::size_t i = 0; i < serial.points.size(); ++i) {
    REQUIRE(serial.points[i].x == parallel.points[i].x);
    REQUIRE(serial.points[i].y == parallel.points[i].y);
  }
  REQUIRE(serial.stats.rms_centroid == parallel.stats.rms_centroid);
  REQUIRE(serial.stats.rms_chief == parallel.stats.rms_chief);
}

TEST_CASE("chief ray that does not reach the image is an AnalysisError", "[spot]") {
  // 80 deg on the singlet: the chief ray passes the stop centre but misses the lens aperture.
  System s = load("m1/singlet_const.rtt.json");
  s.fields.points.push_back({0.0, 80.0, 1.0});
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(s, lib);
  REQUIRE_THROWS_AS(rtt::analysis::spot(cs, PathId{0}, 3, 1), AnalysisError);
}

TEST_CASE("invalid spot input throws std::invalid_argument", "[spot]") {
  const MaterialLibrary lib;
  System s = load("m1/singlet_const.rtt.json");
  const CompiledSystem cs = compile(s, lib);
  REQUIRE_THROWS_AS(rtt::analysis::spot(cs, PathId{0}, 7, 1), std::invalid_argument);
  REQUIRE_THROWS_AS(rtt::analysis::spot(cs, PathId{0}, 0, 9), std::invalid_argument);
  REQUIRE_THROWS_AS(rtt::analysis::spot(cs, PathId{3}, 0, 1), std::invalid_argument);
  s.wavelengths = {{0.4861, 0.0, false}, {0.5876, 0.0, true}};
  const CompiledSystem zero = compile(s, lib);
  REQUIRE_THROWS_AS(rtt::analysis::spot(zero, PathId{0}, 0, std::nullopt), std::invalid_argument);
  // Negative or non-finite weights never reach the analysis: compile() rejects them
  // (model::validate), so a CompiledSystem always has weights >= 0; spot() checks again.
  for (const double bad :
       {-1.0, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity()}) {
    s.wavelengths = {{0.4861, bad, false}, {0.5876, 1.0, true}};
    REQUIRE_THROWS_AS(compile(s, lib), rtt::compile::CompileError);
  }
}
