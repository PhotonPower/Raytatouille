// Analyses of centrally obscured systems with virtual reference rays (#229; ADR 0032, addendum
// #229). Reference rays (chief, parabasal, +-delta and zone rays) skip the inner radius of
// circular apertures, so every result of an obscured system that depends only on reference rays
// and on rays in the annulus is bitwise that of the same system without the obscuration. Sampled
// rays in the obscuration are lost as before; stop.clips_beam reports the outer rim only.
//
// System: tests/reference/m2/paraboloid_stop (stop 50 mm in front of a paraboloid, f = 100 mm,
// object at infinity) with EPD 10 mm and the stop aperture {5 mm, inner 2 mm}: obscuration 0.4.
// HexapolarPupil{6} has rings at rho = k / 6 with 6 k rays: rings 1 and 2 (rho <= 1/3, heights
// up to 1.67 mm) and the centre lie in the obscuration, 19 of 127 rays.

#include <algorithm>
#include <bit>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <variant>
#include <vector>

#include "rtt/analysis/chromatic.hpp"
#include "rtt/analysis/field.hpp"
#include "rtt/analysis/jones.hpp"
#include "rtt/analysis/opd.hpp"
#include "rtt/analysis/pupil.hpp"
#include "rtt/analysis/spot.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/sources.hpp"

using rtt::analysis::AnalysisError;
using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::model::Element;
using rtt::model::System;
using rtt::trace::RayStatus;

namespace {

constexpr std::size_t kObscuredRays = 19;  // centre and rings 1, 2 of HexapolarPupil{6}

System paraboloid(double epd, double inner_radius) {
  System s = rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m2/paraboloid_stop.rtt.json");
  s.aperture.value = rtt::model::Param(epd);
  std::get<Element>(s.root.children[0].value).surfaces[0].aperture =
      rtt::model::CircularAperture{5.0, inner_radius};
  return s;
}

CompiledSystem compiled(const System& s) {
  return rtt::compile::compile(s, rtt::material::MaterialLibrary{});
}

bool same_bits(double a, double b) {
  return std::bit_cast<std::uint64_t>(a) == std::bit_cast<std::uint64_t>(b);
}

bool has(const std::vector<rtt::model::Diagnostic>& warnings, const std::string& code) {
  return std::any_of(warnings.begin(), warnings.end(),
                     [&](const rtt::model::Diagnostic& d) { return d.code == code; });
}

std::uint32_t stop_surface(const CompiledSystem& cs) {
  return cs.path(PathId{0}).events.front().surface;  // the stop is the first event
}

}  // namespace

TEST_CASE("spot, fans and OPD of an obscured mirror use the virtual chief ray", "[obscured]") {
  const CompiledSystem obscured = compiled(paraboloid(10.0, 2.0));
  const CompiledSystem open = compiled(paraboloid(10.0, 0.0));
  REQUIRE(obscured.surfaces()[stop_surface(obscured)].element_kind ==
          rtt::model::ElementKind::Stop);
  // Spot: the chief ray is the one of the open system; the obscured rays are lost at the stop,
  // which is not a clipped beam.
  const rtt::analysis::SpotDiagram so = rtt::analysis::spot(obscured, PathId{0}, 0, std::nullopt);
  const rtt::analysis::SpotDiagram sp = rtt::analysis::spot(open, PathId{0}, 0, std::nullopt);
  CHECK(same_bits(so.chief.x, sp.chief.x));
  CHECK(same_bits(so.chief.y, sp.chief.y));
  CHECK(so.losses.launched == 127);
  CHECK(so.losses.count(RayStatus::Vignetted) == kObscuredRays);
  CHECK(so.losses.worst_surface == std::optional<std::uint32_t>(stop_surface(obscured)));
  CHECK(sp.losses.count(RayStatus::Vignetted) == 0);  // content guard
  CHECK_FALSE(has(so.warnings, "stop.clips_beam"));
  // Ray fans run with the same chief ray.
  CHECK_NOTHROW(rtt::analysis::ray_fan(obscured, PathId{0}, 0, 0));
  // OPD: the points in the annulus are bitwise those of the open system (same reference sphere,
  // same rays); the obscured points are lost.
  const auto po =
      rtt::analysis::opd_points(obscured, PathId{0}, 0, 0, rtt::trace::HexapolarPupil{6});
  const auto pp = rtt::analysis::opd_points(open, PathId{0}, 0, 0, rtt::trace::HexapolarPupil{6});
  REQUIRE(po.points.size() == pp.points.size());
  std::size_t lost = 0;
  for (std::size_t i = 0; i < po.points.size(); ++i) {
    INFO("point " << i);
    if (po.points[i].status != RayStatus::Alive) {
      ++lost;
      continue;
    }
    CHECK(same_bits(po.points[i].w, pp.points[i].w));
  }
  CHECK(lost == kObscuredRays);
  CHECK_FALSE(has(po.warnings, "stop.clips_beam"));
}

TEST_CASE("pupil function of an obscured mirror: the reference and M of the open system",
          "[obscured]") {
  const CompiledSystem obscured = compiled(paraboloid(10.0, 2.0));
  const CompiledSystem open = compiled(paraboloid(10.0, 0.0));
  rtt::analysis::PupilFunctionOptions o;
  o.grid = 16;
  const rtt::analysis::PupilFunction fo =
      rtt::analysis::pupil_function(obscured, PathId{0}, 0, 0, o);
  const rtt::analysis::PupilFunction fp = rtt::analysis::pupil_function(open, PathId{0}, 0, 0, o);
  // The mapping M comes from virtual parabasal rays at |p| = 1e-3, deep in the obscuration.
  if (!fo.mapping.has_value() || !fp.mapping.has_value()) {
    FAIL("no pupil mapping");
    return;
  }
  const rtt::analysis::PupilMapping& mo = *fo.mapping;
  const rtt::analysis::PupilMapping& mp = *fp.mapping;
  for (std::size_t k = 0; k < 4; ++k) CHECK(same_bits(mo.m[k], mp.m[k]));
  for (int k = 0; k < 3; ++k) {
    CHECK(same_bits(mo.centre[k], mp.centre[k]));
    CHECK(same_bits(mo.d_c[k], mp.d_c[k]));
  }
  // Cells (annular cells, ADR 0032, addendum #229, PR C): f = 0 in the obscuration; a cell on the
  // inner rim arrives as a virtual ray at the centroid of its annular part from pupil_cells(N,
  // epsilon); every other arrived cell is bitwise that of the open system.
  const rtt::analysis::PupilCells ring = rtt::analysis::pupil_cells(16, 0.4);
  std::size_t obscured_cells = 0;
  std::size_t rim_cells = 0;
  for (std::size_t j = 0; j < fo.fraction.size(); ++j) {
    INFO("cell " << j);
    if (fo.fraction[j] == 0.0) {
      if (fp.fraction[j] > 0.0) ++obscured_cells;
      continue;
    }
    CHECK(fo.status[j] == RayStatus::Alive);
    if (ring.inner_rim[j] != 0) {
      ++rim_cells;
      CHECK(same_bits(fo.px[j], ring.px[j]));
      CHECK(same_bits(fo.py[j], ring.py[j]));
    } else {
      CHECK(same_bits(fo.w[j], fp.w[j]));
      CHECK(same_bits(fo.amplitude[j], fp.amplitude[j]));
    }
  }
  CHECK(obscured_cells > 0);  // content guards
  CHECK(rim_cells > 0);
  CHECK_FALSE(has(fo.warnings, "stop.clips_beam"));
}

TEST_CASE("Jones pupil of an obscured mirror: the chief ray keeps its J", "[obscured]") {
  // Convenience form with field amplitudes and the global phase: the chief ray of the global
  // phase is traced again for its power factors (#208 B1). As a virtual reference ray it reaches
  // the image, so there is no pupil.chief_phase_unavailable and psi is that of the open system.
  // m1/singlet_const (EPD 20 mm, stop radius 10 mm) with an obscuration of 3 mm at the stop: its
  // two refractions carry power factors, so a real, lost chief ray would give NaN there.
  const auto singlet = [](double inner_radius) {
    System s = rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m1/singlet_const.rtt.json");
    std::get<Element>(s.root.children[0].value).surfaces[0].aperture =
        rtt::model::CircularAperture{10.0, inner_radius};
    return compiled(s);
  };
  const CompiledSystem obscured = singlet(3.0);
  const CompiledSystem open = singlet(0.0);
  rtt::analysis::JonesPupilOptions o;
  o.sampling = rtt::trace::GridPupil{9};
  o.field_amplitude = true;
  o.remove_global_phase = true;
  const auto jo = rtt::analysis::jones_pupil(obscured, PathId{0}, 0, 0, o);
  const auto jp = rtt::analysis::jones_pupil(open, PathId{0}, 0, 0, o);
  CHECK_FALSE(has(jo.warnings, "pupil.chief_phase_unavailable"));
  CHECK(same_bits(jo.global_phase, jp.global_phase));
  REQUIRE(jo.j.size() == jp.j.size());
  std::size_t compared = 0;
  for (std::size_t i = 0; i < jo.status.size(); ++i) {
    if (jo.status[i] != RayStatus::Alive) continue;
    INFO("point " << i);
    for (std::size_t e = 0; e < 4; ++e) {
      CHECK(same_bits(jo.j[4 * i + e].real(), jp.j[4 * i + e].real()));
      CHECK(same_bits(jo.j[4 * i + e].imag(), jp.j[4 * i + e].imag()));
    }
    ++compared;
  }
  CHECK(compared > 0);                 // content guard: points in the annulus
  CHECK(compared < jo.status.size());  // and points in the obscuration
}

TEST_CASE("distortion, field curvature and colour use virtual reference rays", "[obscured]") {
  // An off-axis field of 0.5 deg and three wavelengths; the longitudinal colour with the zone
  // ray at 0.3, inside the obscuration.
  const auto with_field_and_colours = [](System s) {
    s.fields.points = {rtt::model::Field{}, rtt::model::Field{0.0, 0.5, 1.0}};
    s.wavelengths = {{0.5876, 1.0, true}, {0.4861, 1.0, false}, {0.6563, 1.0, false}};
    return compiled(s);
  };
  const CompiledSystem obscured = with_field_and_colours(paraboloid(10.0, 2.0));
  const CompiledSystem open = with_field_and_colours(paraboloid(10.0, 0.0));
  const rtt::model::Field off_axis{0.0, 0.5, 1.0};
  const auto d_o = rtt::analysis::distortion_at(obscured, PathId{0}, off_axis, 0);
  const auto d_p = rtt::analysis::distortion_at(open, PathId{0}, off_axis, 0);
  CHECK(same_bits(d_o.real_height, d_p.real_height));
  const auto c_o = rtt::analysis::field_curvature_at(obscured, PathId{0}, off_axis, 0);
  const auto c_p = rtt::analysis::field_curvature_at(open, PathId{0}, off_axis, 0);
  CHECK(same_bits(c_o.tangential, c_p.tangential));
  CHECK(same_bits(c_o.sagittal, c_p.sagittal));
  const auto l_o = rtt::analysis::lateral_colour(obscured, PathId{0}, 1);
  const auto l_p = rtt::analysis::lateral_colour(open, PathId{0}, 1);
  REQUIRE(l_o.chief.size() == 3);
  REQUIRE(l_p.chief.size() == 3);
  for (std::size_t k = 0; k < 3; ++k) {
    CHECK(same_bits(l_o.chief[k].x, l_p.chief[k].x));
    CHECK(same_bits(l_o.chief[k].y, l_p.chief[k].y));
  }
  rtt::analysis::ChromaticOptions zone;
  zone.zone = 0.3;
  const auto z_o = rtt::analysis::longitudinal_colour(obscured, PathId{0}, zone);
  const auto z_p = rtt::analysis::longitudinal_colour(open, PathId{0}, zone);
  REQUIRE(z_o.foci.size() == 3);
  REQUIRE(z_p.foci.size() == 3);
  for (std::size_t k = 0; k < 3; ++k) CHECK(same_bits(z_o.foci[k].real_z, z_p.foci[k].real_z));
}

TEST_CASE("stop.clips_beam counts the rays clipped by the outer rim only", "[obscured]") {
  // EPD 11: ring 6 of HexapolarPupil{6} (36 rays, height 5.5 mm) passes the outer rim of the
  // stop (5 mm), the 19 rays up to 1.83 mm lie in the obscuration (2 mm).
  const CompiledSystem cs = compiled(paraboloid(11.0, 2.0));
  const rtt::analysis::SpotDiagram s = rtt::analysis::spot(cs, PathId{0}, 0, std::nullopt);
  CHECK(s.losses.count(RayStatus::Vignetted) == kObscuredRays + 36);
  const auto clips = std::find_if(s.warnings.begin(), s.warnings.end(),
                                  [](const auto& d) { return d.code == "stop.clips_beam"; });
  REQUIRE(clips != s.warnings.end());
  CHECK(clips->message.rfind("36 rays ended Vignetted at the stop", 0) == 0);
}

TEST_CASE("a chief ray lost at an outer rim is still an AnalysisError", "[obscured]") {
  // Field 35 deg: the chief ray through the stop centre meets the mirror plane about
  // 50 tan(35 deg) = 35 mm off axis, beyond its aperture of 31 mm (ADR 0022, point 5).
  System s = paraboloid(10.0, 2.0);
  s.fields.points = {rtt::model::Field{0.0, 35.0, 1.0}};
  const CompiledSystem cs = compiled(s);
  try {
    static_cast<void>(rtt::analysis::spot(cs, PathId{0}, 0, std::nullopt));
    FAIL("no AnalysisError");
  } catch (const AnalysisError& e) {
    CHECK(std::string(e.what()).find("chief ray") != std::string::npos);
    CHECK(e.ray_status() == std::optional<RayStatus>(RayStatus::Vignetted));
    REQUIRE(e.surface().has_value());
    CHECK(e.surface()->str() == "M");
  }
}
