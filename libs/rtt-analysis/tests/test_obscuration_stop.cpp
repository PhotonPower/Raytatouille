// Obscurations in the analyses (#320, ADR 0050, point 3 and point 4, PR 4).
//
// Equality: a centred CircularObscuration of radius r_i at the stop is treated exactly as the
// inner radius r_i of the stop's circular aperture (fan statuses, spot, OPD, PSF, epsilon of the
// pupil function, count of stop.clips_beam). Two levels: reference rays ignore every obscuration
// (the chief ray with a spider is bitwise the one without), the inner-rim cells of the pupil
// function ignore only the centred circles (a spider still blocks them).
//
// System: tests/reference/m2/paraboloid_stop as in test_obscured.cpp (stop 50 mm in front of a
// paraboloid, object at infinity), stop aperture radius 5 mm.

#include <algorithm>
#include <bit>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <optional>
#include <string>
#include <variant>
#include <vector>

#include "rtt/analysis/opd.hpp"
#include "rtt/analysis/psf.hpp"
#include "rtt/analysis/pupil.hpp"
#include "rtt/analysis/spot.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/sources.hpp"

using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::model::CircularAperture;
using rtt::model::CircularObscuration;
using rtt::model::Element;
using rtt::model::Obscuration;
using rtt::model::Param;
using rtt::model::RectangularObscuration;
using rtt::model::System;
using rtt::trace::RayStatus;

namespace {

/// paraboloid_stop with EPD `epd`, the stop aperture {5 mm, inner_radius} and `obscurations`
/// at the stop.
System paraboloid(double epd, double inner_radius, std::vector<Obscuration> obscurations) {
  System s = rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m2/paraboloid_stop.rtt.json");
  s.aperture.value = Param(epd);
  rtt::model::Surface& stop = std::get<Element>(s.root.children[0].value).surfaces[0];
  stop.aperture = CircularAperture{Param(5.0), Param(inner_radius)};
  stop.obscurations = std::move(obscurations);
  return s;
}

CompiledSystem compiled(const System& s) {
  return rtt::compile::compile(s, rtt::material::MaterialLibrary{});
}

bool same_bits(double a, double b) {
  return std::bit_cast<std::uint64_t>(a) == std::bit_cast<std::uint64_t>(b);
}

std::size_t clipped(const std::vector<rtt::model::Diagnostic>& warnings) {
  const auto it = std::find_if(warnings.begin(), warnings.end(),
                               [](const auto& d) { return d.code == "stop.clips_beam"; });
  if (it == warnings.end()) return 0;
  return std::stoul(it->message);  // "<n> rays ended Vignetted at the stop ..."
}

CircularObscuration centred(double r) {
  return CircularObscuration{Param(r), Param(0.0), Param(0.0)};
}

/// A bar through the centre along x: |y| < h over the whole stop.
RectangularObscuration spider(double h) {
  return RectangularObscuration{Param(6.0), Param(h), Param(0.0), Param(0.0), 0.0};
}

}  // namespace

TEST_CASE("a centred obscuration at the stop equals the inner radius (ADR 0050, point 4)",
          "[obscuration]") {
  // EPD 11: the outer ring of the spot sampling is clipped by the outer rim, so stop.clips_beam
  // has a count to compare; the rays in the obscuration do not count for it.
  const CompiledSystem ring = compiled(paraboloid(11.0, 2.0, {}));
  const CompiledSystem disc = compiled(paraboloid(11.0, 0.0, {centred(2.0)}));

  // Fans: statuses and aberrations.
  const rtt::analysis::RayFan fr = rtt::analysis::ray_fan(ring, PathId{0}, 0, 0);
  const rtt::analysis::RayFan fd = rtt::analysis::ray_fan(disc, PathId{0}, 0, 0);
  REQUIRE(fr.tangential.size() == fd.tangential.size());
  std::size_t fan_lost = 0;
  for (std::size_t i = 0; i < fr.tangential.size(); ++i) {
    CHECK(fr.tangential[i].status == fd.tangential[i].status);
    CHECK(same_bits(fr.tangential[i].ey, fd.tangential[i].ey));
    CHECK(fr.sagittal[i].status == fd.sagittal[i].status);
    CHECK(same_bits(fr.sagittal[i].ex, fd.sagittal[i].ex));
    if (fr.tangential[i].status != RayStatus::Alive) ++fan_lost;
  }
  CHECK(fan_lost > 0);  // content guard: the obscuration takes rays of the fan
  CHECK(same_bits(fr.chief.y, fd.chief.y));

  // Spot: positions, losses and stop.clips_beam.
  const rtt::analysis::SpotDiagram sr = rtt::analysis::spot(ring, PathId{0}, 0, std::nullopt);
  const rtt::analysis::SpotDiagram sd = rtt::analysis::spot(disc, PathId{0}, 0, std::nullopt);
  REQUIRE(sr.points.size() == sd.points.size());
  for (std::size_t i = 0; i < sr.points.size(); ++i) {
    CHECK(same_bits(sr.points[i].x, sd.points[i].x));
    CHECK(same_bits(sr.points[i].y, sd.points[i].y));
  }
  CHECK(sr.arrived == sd.arrived);
  CHECK(sr.losses.by_status == sd.losses.by_status);
  CHECK(sr.losses.count(RayStatus::Vignetted) == 19 + 36);  // as in test_obscured.cpp
  CHECK(clipped(sr.warnings) == 36);
  CHECK(clipped(sd.warnings) == 36);

  // OPD points.
  const auto pr = rtt::analysis::opd_points(ring, PathId{0}, 0, 0, rtt::trace::HexapolarPupil{6});
  const auto pd = rtt::analysis::opd_points(disc, PathId{0}, 0, 0, rtt::trace::HexapolarPupil{6});
  REQUIRE(pr.points.size() == pd.points.size());
  for (std::size_t i = 0; i < pr.points.size(); ++i) {
    CHECK(pr.points[i].status == pd.points[i].status);
    CHECK(same_bits(pr.points[i].w, pd.points[i].w));
  }
}

TEST_CASE("pupil function and PSF: a centred obscuration gives the exact ring cells",
          "[obscuration]") {
  const CompiledSystem ring = compiled(paraboloid(10.0, 2.0, {}));
  const CompiledSystem disc = compiled(paraboloid(10.0, 0.0, {centred(2.0)}));
  // The largest centred circle wins over a smaller inner radius (effective r_i, ADR 0050).
  const CompiledSystem both = compiled(paraboloid(10.0, 1.0, {centred(2.0), centred(1.5)}));
  rtt::analysis::PupilFunctionOptions o;
  o.grid = 16;
  const auto fr = rtt::analysis::pupil_function(ring, PathId{0}, 0, 0, o);
  const auto fd = rtt::analysis::pupil_function(disc, PathId{0}, 0, 0, o);
  const auto fb = rtt::analysis::pupil_function(both, PathId{0}, 0, 0, o);
  CHECK(fr.obscuration == 0.4);  // 2 mm / 5 mm
  CHECK(same_bits(fd.obscuration, fr.obscuration));
  CHECK(same_bits(fb.obscuration, fr.obscuration));
  REQUIRE(fr.fraction.size() == fd.fraction.size());
  std::size_t arrived = 0;
  for (std::size_t j = 0; j < fr.fraction.size(); ++j) {
    INFO("cell " << j);
    CHECK(same_bits(fr.fraction[j], fd.fraction[j]));
    CHECK(fr.status[j] == fd.status[j]);
    CHECK(same_bits(fr.w[j], fd.w[j]));
    CHECK(same_bits(fr.amplitude[j], fd.amplitude[j]));
    CHECK(fb.status[j] == fr.status[j]);
    if (fr.status[j] == RayStatus::Alive) ++arrived;
  }
  CHECK(arrived > 0);
  CHECK(arrived < fr.fraction.size());
  // Content guard: ring cells with a partial area at the inner rim exist in both.
  const rtt::analysis::PupilCells cells = rtt::analysis::pupil_cells(16, 0.4);
  std::size_t partial_rim = 0;
  for (std::size_t j = 0; j < cells.fraction.size(); ++j) {
    if (cells.inner_rim[j] != 0 && fd.fraction[j] > 0.0 && fd.fraction[j] < 1.0) ++partial_rim;
  }
  CHECK(partial_rim > 0);

  rtt::analysis::PsfOptions p;
  p.pupil.grid = 16;
  const rtt::analysis::Psf psf_r = rtt::analysis::psf(ring, PathId{0}, 0, p);
  const rtt::analysis::Psf psf_d = rtt::analysis::psf(disc, PathId{0}, 0, p);
  REQUIRE(psf_r.intensity.size() == psf_d.intensity.size());
  for (std::size_t k = 0; k < psf_r.intensity.size(); ++k) {
    if (!same_bits(psf_r.intensity[k], psf_d.intensity[k])) {
      FAIL("PSF differs at " << k);
    }
  }
  CHECK(same_bits(psf_r.strehl, psf_d.strehl));
  CHECK(psf_r.peak.value > 0.0);  // content guard
}

TEST_CASE("two levels: reference rays ignore a spider, inner-rim cells do not (ADR 0050)",
          "[obscuration]") {
  // A spider bar |y| < 1 mm through the centre of the stop.
  const CompiledSystem open = compiled(paraboloid(10.0, 0.0, {}));
  const CompiledSystem with_spider = compiled(paraboloid(10.0, 0.0, {spider(1.0)}));
  // The chief ray runs through the centre of the bar; as a reference ray it is the one of the
  // system without the spider.
  const rtt::analysis::SpotDiagram so = rtt::analysis::spot(open, PathId{0}, 0, std::nullopt);
  const rtt::analysis::SpotDiagram ss =
      rtt::analysis::spot(with_spider, PathId{0}, 0, std::nullopt);
  REQUIRE(std::isfinite(so.chief.x));  // content guard: a real chief ray, not NaN
  REQUIRE(std::isfinite(so.chief.y));
  CHECK(same_bits(so.chief.x, ss.chief.x));
  CHECK(same_bits(so.chief.y, ss.chief.y));
  // Sampled rays on the bar are lost, which is no clipping of the beam.
  CHECK(ss.losses.count(RayStatus::Vignetted) > 0);
  CHECK(so.losses.count(RayStatus::Vignetted) == 0);
  CHECK(clipped(ss.warnings) == 0);

  // Ring (inner radius 2 mm) and spider: a cell on the inner rim whose ray hits the bar ends
  // Vignetted; the other rim cells arrive (they ignore the centred circle).
  const CompiledSystem ring_spider = compiled(paraboloid(10.0, 2.0, {spider(1.0)}));
  rtt::analysis::PupilFunctionOptions o;
  o.grid = 16;
  const auto f = rtt::analysis::pupil_function(ring_spider, PathId{0}, 0, 0, o);
  const rtt::analysis::PupilCells cells = rtt::analysis::pupil_cells(16, 0.4);
  std::size_t rim_on_bar = 0;
  std::size_t rim_free = 0;
  for (std::size_t j = 0; j < f.fraction.size(); ++j) {
    if (!(f.fraction[j] > 0.0) || cells.inner_rim[j] == 0) continue;
    INFO("cell " << j << " at (" << f.px[j] << ", " << f.py[j] << ")");
    // The stop height of pupil point p is |p| R_s with R_s = 5 mm (ADR 0032, addendum #229).
    const double y = 5.0 * f.py[j];
    if (std::abs(y) < 0.9) {
      CHECK(f.end.ray_status[j] == RayStatus::Vignetted);
      ++rim_on_bar;
    } else if (std::abs(y) > 1.1) {
      CHECK(f.end.ray_status[j] == RayStatus::Alive);
      ++rim_free;
    }
  }
  CHECK(rim_on_bar > 0);  // content guards
  CHECK(rim_free > 0);
}

TEST_CASE("stop.clips_beam: a spider arm beyond the outer rim does not hide clipping (ADR 0050)",
          "[obscuration]") {
  // EPD 11: ring 6 of HexapolarPupil{6} (36 rays at 5.5 mm, every 10 deg) lies beyond the outer
  // rim (5 mm). A bar |y| < 1 mm reaching to |x| = 6 mm covers 6 of them (0, 10, 170, 180, 190
  // and 350 deg: |5.5 sin| < 1); they are still clipped, not obscured.
  const CompiledSystem cs = compiled(paraboloid(11.0, 0.0, {spider(1.0)}));
  const rtt::analysis::SpotDiagram s = rtt::analysis::spot(cs, PathId{0}, 0, std::nullopt);
  CHECK(clipped(s.warnings) == 36);
  CHECK(5.5 * std::sin(10.0 * std::numbers::pi / 180.0) < 1.0);  // the hand numbers
  CHECK(5.5 * std::sin(20.0 * std::numbers::pi / 180.0) > 1.0);
}
