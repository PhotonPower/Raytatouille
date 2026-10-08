#include <catch2/catch_test_macros.hpp>
#include <chrono>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include "rtt/analysis/chromatic.hpp"
#include "rtt/analysis/field.hpp"
#include "rtt/analysis/opd.hpp"
#include "rtt/analysis/spot.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/trace/run_control.hpp"

// Cancellation and progress of the analyses (#83): an active control never changes a result;
// a cancelled analysis throws trace::Cancelled; the callback sees the documented stages.

using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::trace::Cancelled;
using rtt::trace::CancelToken;
using rtt::trace::Progress;
using rtt::trace::RunControl;

namespace {

CompiledSystem singlet() {
  const rtt::material::MaterialLibrary lib;
  return rtt::compile::compile(
      rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m1/singlet_const.rtt.json"), lib);
}

/// Active control that never cancels: a token and a callback, small blocks.
RunControl active_control() {
  RunControl control;
  control.cancel = CancelToken();
  control.progress = [](const Progress&) {};
  control.min_interval = std::chrono::milliseconds{0};
  control.block_size = 7;
  return control;
}

RunControl cancelled_control() {
  RunControl control;
  control.cancel = CancelToken();
  control.cancel->request_cancel();
  return control;
}

}  // namespace

TEST_CASE("analyses with an active control give the same results, bit for bit (#83)",
          "[run_control]") {
  const CompiledSystem cs = singlet();
  const PathId path{0};
  const RunControl control = active_control();

  SECTION("spot (polychromatic) and ray fan") {
    const auto a = rtt::analysis::spot(cs, path, 2, std::nullopt);
    const auto b = rtt::analysis::spot(cs, path, 2, std::nullopt, {}, control);
    REQUIRE(a.points.size() == b.points.size());
    for (std::size_t i = 0; i < a.points.size(); ++i) {
      REQUIRE(a.points[i].x == b.points[i].x);
      REQUIRE(a.points[i].y == b.points[i].y);
      REQUIRE(a.points[i].weight == b.points[i].weight);
      REQUIRE(a.points[i].wavelength == b.points[i].wavelength);
    }
    REQUIRE(a.stats.rms_centroid == b.stats.rms_centroid);
    const auto fa = rtt::analysis::ray_fan(cs, path, 1, 1);
    const auto fb = rtt::analysis::ray_fan(cs, path, 1, 1, {}, control);
    REQUIRE(fa.tangential.size() == fb.tangential.size());
    for (std::size_t i = 0; i < fa.tangential.size(); ++i) {
      REQUIRE(fa.tangential[i].ey == fb.tangential[i].ey);
      REQUIRE(fa.sagittal[i].ex == fb.sagittal[i].ex);
    }
  }
  SECTION("OPD map and fan") {
    const auto a = rtt::analysis::opd_map(cs, path, 2, 1);
    const auto b = rtt::analysis::opd_map(cs, path, 2, 1, {}, control);
    REQUIRE(a.points.size() == b.points.size());
    for (std::size_t i = 0; i < a.points.size(); ++i) REQUIRE(a.points[i].w == b.points[i].w);
    REQUIRE(a.rms == b.rms);
    const auto fa = rtt::analysis::opd_fan(cs, path, 1, 1);
    const auto fb = rtt::analysis::opd_fan(cs, path, 1, 1, {}, control);
    for (std::size_t i = 0; i < fa.tangential.size(); ++i) {
      REQUIRE(fa.tangential[i].w == fb.tangential[i].w);
      REQUIRE(fa.sagittal[i].w == fb.sagittal[i].w);
    }
  }
  SECTION("field sweeps and colour") {
    const auto da = rtt::analysis::distortion(cs, path, 1);
    const auto db = rtt::analysis::distortion(cs, path, 1, {}, control);
    REQUIRE(da.size() == db.size());
    for (std::size_t i = 0; i < da.size(); ++i) {
      REQUIRE(da[i].real_height == db[i].real_height);
      REQUIRE(da[i].percent == db[i].percent);
    }
    const auto ca = rtt::analysis::field_curvature(cs, path, 1);
    const auto cb = rtt::analysis::field_curvature(cs, path, 1, {}, control);
    for (std::size_t i = 0; i < ca.size(); ++i) {
      REQUIRE(ca[i].tangential == cb[i].tangential);
      REQUIRE(ca[i].sagittal == cb[i].sagittal);
    }
    const auto la = rtt::analysis::longitudinal_colour(cs, path);
    const auto lb = rtt::analysis::longitudinal_colour(cs, path, {}, control);
    REQUIRE(la.real == lb.real);
    REQUIRE(la.paraxial == lb.paraxial);
    for (std::size_t i = 0; i < la.foci.size(); ++i)
      REQUIRE(la.foci[i].real_z == lb.foci[i].real_z);
    const auto ta = rtt::analysis::lateral_colour(cs, path, 2);
    const auto tb = rtt::analysis::lateral_colour(cs, path, 2, rtt::trace::Aiming::Real, control);
    for (std::size_t i = 0; i < ta.offset.size(); ++i) {
      REQUIRE(ta.chief[i].x == tb.chief[i].x);
      REQUIRE(ta.chief[i].y == tb.chief[i].y);
    }
  }
}

TEST_CASE("a cancelled analysis throws Cancelled (#83)", "[run_control]") {
  const CompiledSystem cs = singlet();
  const PathId path{0};
  const RunControl control = cancelled_control();
  REQUIRE_THROWS_AS(rtt::analysis::spot(cs, path, 2, std::nullopt, {}, control), Cancelled);
  REQUIRE_THROWS_AS(rtt::analysis::ray_fan(cs, path, 1, 1, {}, control), Cancelled);
  REQUIRE_THROWS_AS(rtt::analysis::opd_map(cs, path, 1, 1, {}, control), Cancelled);
  REQUIRE_THROWS_AS(rtt::analysis::opd_fan(cs, path, 1, 1, {}, control), Cancelled);
  REQUIRE_THROWS_AS(rtt::analysis::distortion(cs, path, 1, {}, control), Cancelled);
  REQUIRE_THROWS_AS(rtt::analysis::field_curvature(cs, path, 1, {}, control), Cancelled);
  REQUIRE_THROWS_AS(rtt::analysis::longitudinal_colour(cs, path, {}, control), Cancelled);
  REQUIRE_THROWS_AS(rtt::analysis::lateral_colour(cs, path, 2, rtt::trace::Aiming::Real, control),
                    Cancelled);
}

TEST_CASE("analyses report their stages, each ending with done == total (#83)", "[run_control]") {
  const CompiledSystem cs = singlet();
  const PathId path{0};
  std::vector<std::string> stages;
  std::vector<Progress> calls;
  RunControl control;
  control.min_interval = std::chrono::milliseconds{0};
  control.progress = [&](const Progress& p) {
    stages.emplace_back(p.stage);
    calls.push_back(p);
  };
  const auto final_of = [&](const std::string& stage) {
    for (std::size_t k = calls.size(); k-- > 0;) {
      if (stages[k] == stage) return calls[k];
    }
    FAIL("stage " << stage << " not reported");
    return Progress{};
  };

  SECTION("spot: aim and trace") {
    const auto d = rtt::analysis::spot(cs, path, 1, std::uint16_t{1}, {}, control);
    REQUIRE(final_of("aim").done == d.rays_launched);
    REQUIRE(final_of("trace").done == d.rays_launched);
    REQUIRE(final_of("trace").total == d.rays_launched);
  }
  SECTION("distortion: field") {
    const auto d = rtt::analysis::distortion(cs, path, 1, {5, rtt::trace::Aiming::Real}, control);
    REQUIRE(final_of("field").done == 5);
    REQUIRE(final_of("field").total == 5);
    REQUIRE(d.size() == 5);
  }
  SECTION("longitudinal colour: wavelength") {
    [[maybe_unused]] const auto lc = rtt::analysis::longitudinal_colour(cs, path, {}, control);
    REQUIRE(final_of("wavelength").done == cs.wavelengths_um().size());
  }
}

TEST_CASE("an exception of the progress callback leaves the analysis (#83)", "[run_control]") {
  const CompiledSystem cs = singlet();
  RunControl control;
  control.progress = [](const Progress&) { throw std::domain_error("stop"); };
  REQUIRE_THROWS_AS(rtt::analysis::spot(cs, PathId{0}, 1, std::uint16_t{1}, {}, control),
                    std::domain_error);
  REQUIRE_THROWS_AS(rtt::analysis::distortion(cs, PathId{0}, 1, {}, control), std::domain_error);
}
