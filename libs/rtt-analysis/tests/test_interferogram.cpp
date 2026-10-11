// Interferogram of two paths (#210, M6, roadmap W1). Reference: the Michelson of
// tests/reference/m4/michelson_offset.rtt.json (VACUUM, ideal splitter R = 0.5, ideal mirrors,
// test mirror Delta = 7.5 mm farther). Sources: Konijnenberg, Adam, Urbach, BSc Optics, 2nd ed.,
// TU Delft Open 2024: two-beam interference Sec. 5.2, Eq. (5.12), equal amplitudes Eqs. (5.10),
// (5.11); plane wave Sec. 1.5.1, Eqs. (1.37)-(1.39) (docs/quellen.md).

#include <tbb/global_control.h>

#include <algorithm>
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "rtt/analysis/paths.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/math/types.hpp"
#include "rtt/model/model.hpp"
#include "rtt/paraxial/paraxial.hpp"
#include "rtt/trace/run_control.hpp"
#include "rtt/trace/sources.hpp"

using rtt::analysis::Interferogram;
using rtt::analysis::InterferogramArrival;
using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::material::MaterialLibrary;
using rtt::model::Element;
using rtt::model::EventKind;
using rtt::model::SurfaceId;
using rtt::model::System;
using rtt::trace::RayBatch;
using rtt::trace::RayStatus;

namespace {

constexpr double kLambdaMm = 0.6328e-3;  // the system wavelength of the Michelson, mm

System load(const std::string& relative) {
  return rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/" + relative);
}

/// tests/reference/m4/michelson_offset.rtt.json: stop at z = 10, splitter under 45 deg at
/// z = 50, reference mirror at z = 100, test mirror at y = 57.5 (rotation -90 deg about x),
/// camera at y = -50. Paths "reference arm" (a) and "test arm" (b) to the camera.
System michelson() {
  return load("m4/michelson_offset.rtt.json");
}

Element& element(System& s, const std::string& name) {
  for (auto& node : s.root.children) {
    if (auto* e = std::get_if<Element>(&node.value); e && e->name == name) return *e;
  }
  throw std::logic_error("no element " + name);
}

rtt::model::Path& model_path(System& s, const std::string& name) {
  for (auto& p : s.paths) {
    if (p.name == name) return p;
  }
  throw std::logic_error("no path " + name);
}

PathId path(const CompiledSystem& cs, const std::string& name) {
  const std::optional<PathId> id = cs.find_path(name);
  REQUIRE(id.has_value());
  return *id;
}

std::uint32_t surface(const CompiledSystem& cs, const char* id) {
  const std::optional<std::uint32_t> s = cs.find_surface(SurfaceId(id));
  REQUIRE(s.has_value());
  return *s;
}

/// Collimated start rays along +z at z = 0 at the given points; pupil labels = r / 2.5.
RayBatch bundle(const std::vector<std::pair<double, double>>& points) {
  RayBatch rays(points.size());
  for (std::size_t k = 0; k < points.size(); ++k) {
    rays.pos_x()[k] = points[k].first;
    rays.pos_y()[k] = points[k].second;
    rays.pupil_x()[k] = points[k].first / 2.5;
    rays.pupil_y()[k] = points[k].second / 2.5;
  }
  return rays;
}

/// The axis ray, then a square grid with a pitch of 0.5 mm inside r <= 2.25 mm.
RayBatch grid_bundle() {
  std::vector<std::pair<double, double>> points{{0.0, 0.0}};
  for (int i = -4; i <= 4; ++i) {
    for (int j = -4; j <= 4; ++j) {
      const double x = 0.5 * i;
      const double y = 0.5 * j;
      if ((i != 0 || j != 0) && x * x + y * y <= 2.25 * 2.25) points.emplace_back(x, y);
    }
  }
  return bundle(points);
}

/// One row along the start y axis (x = 0), from -2.25 to 2.25 mm with a pitch of 0.005 mm.
RayBatch row_bundle() {
  std::vector<std::pair<double, double>> points;
  for (int j = -450; j <= 450; ++j) points.emplace_back(0.0, 0.005 * j);
  return bundle(points);
}

/// `start` traced on `p` with the sequential tracer (positions, OPL in the batch).
RayBatch traced(const CompiledSystem& cs, PathId p, const RayBatch& start) {
  RayBatch rays = start;
  [[maybe_unused]] const auto stats = rtt::trace::SequentialTracer().trace(cs, p, rays);
  return rays;
}

/// Local x, y of ray i on surface `s`, mm.
std::pair<double, double> local(const CompiledSystem& cs,
                                const RayBatch& r,
                                std::size_t i,
                                std::uint32_t s) {
  const rtt::math::Vec3 p = cs.surfaces()[s].to_local.apply_point(
      rtt::math::Vec3(r.pos_x()[i], r.pos_y()[i], r.pos_z()[i]));
  return {p.x(), p.y()};
}

/// The Michelson with the test mirror tilted by alpha (rad) about x.
System tilted_michelson(double alpha) {
  System s = michelson();
  element(s, "test mirror").pose.rotation_deg[0].value += alpha * 180.0 / std::numbers::pi;
  return s;
}

bool same_bits(const Interferogram& u, const Interferogram& v) {
  if (u.points.size() != v.points.size() || u.max_shear != v.max_shear) return false;
  for (std::size_t i = 0; i < u.points.size(); ++i) {
    const auto& p = u.points[i];
    const auto& q = v.points[i];
    if (p.x != q.x || p.y != q.y || p.delta != q.delta || p.waves != q.waves ||
        p.phase != q.phase || p.fringe != q.fringe || p.intensity != q.intensity ||
        p.shear != q.shear || p.arrival != q.arrival || p.status_a != q.status_a ||
        p.status_b != q.status_b) {
      return false;
    }
  }
  return true;
}

}  // namespace

TEST_CASE("interferogram: tilted Michelson arm, the wave of b carried to the hits of a (#210)",
          "[interferogram]") {
  // Test mirror tilted by alpha = 1 mrad about x. Collimated input, plane splitter and plane
  // mirrors: both waves at the camera are exactly plane (Sec. 1.5.1, Eq. (1.37)), wave a along
  // the camera normal, wave b at the angle 2 alpha to it in the y-z plane, i.e. along the
  // camera's local y. So delta(x_a) is affine in the local y with the slope n' sin(2 alpha)
  // (VACUUM, n' = 1) and independent of the local x: fringe period lambda / sin(2 alpha).
  // The rays of b land 2 alpha times their path behind the mirror (about 107 mm) away from those
  // of a, so the shear is about 0.2 mm and the transfer term n' d_b . (x_a - x_b) about 4e-4 mm
  // = 0.7 waves: without the transfer the absolute check below fails.
  // Tolerances: OPL sums of about 200 mm, rounding about 1e-13 per sum (as #122): 1e-10 mm.
  const double alpha = 1e-3;
  const double slope = std::sin(2.0 * alpha);
  const MaterialLibrary lib;
  const CompiledSystem cs = rtt::compile::compile(tilted_michelson(alpha), lib);
  const PathId a = path(cs, "reference arm");
  const PathId b = path(cs, "test arm");
  const std::uint32_t det = surface(cs, "DET");
  const RayBatch start = grid_bundle();
  const Interferogram g = rtt::analysis::interferogram(cs, a, b, start);
  REQUIRE(g.points.size() == start.size());
  REQUIRE(g.image_surface == det);
  REQUIRE(g.path_a == a);
  REQUIRE(g.path_b == b);

  // Hits of a, independently traced: the points lie there.
  const RayBatch on_a = traced(cs, a, start);
  for (std::size_t i = 0; i < start.size(); ++i) {
    INFO("ray " << i);
    const auto& p = g.points[i];
    REQUIRE(p.arrival == InterferogramArrival::Both);
    REQUIRE(p.status_a == RayStatus::Alive);
    REQUIRE(p.status_b == RayStatus::Alive);
    REQUIRE(p.px == start.pupil_x()[i]);
    REQUIRE(p.py == start.pupil_y()[i]);
    const auto [x, y] = local(cs, on_a, i, det);
    REQUIRE(std::abs(p.x - x) <= 1e-12);
    REQUIRE(std::abs(p.y - y) <= 1e-12);
  }

  // Affine in the local y with |slope| = sin(2 alpha), none in x. The sign depends on the
  // orientation of the camera's local y and is taken from the data once.
  const auto& p0 = g.points[0];  // the axis ray
  double sign = 0.0;
  for (const auto& p : g.points) {
    if (std::abs(p.y - p0.y) > 1.0) {
      sign = (p.delta - p0.delta) / (p.y - p0.y) > 0.0 ? 1.0 : -1.0;
      break;
    }
  }
  REQUIRE(sign != 0.0);
  for (std::size_t i = 0; i < start.size(); ++i) {
    INFO("ray " << i);
    const auto& p = g.points[i];
    REQUIRE(std::abs(p.delta - p0.delta - sign * slope * (p.y - p0.y)) <= 1e-10);
  }

  // Absolute: the OPL of b at x_a from the ray of b that really lands at x_a. The start plane
  // z = 0 is a wave front of both waves (all OPL 0 there), and the map from the start point to
  // the hit of b is affine (plane surfaces, parallel rays): three rays fix it, a fourth lands at
  // x_a. Independent of the transfer formula of the analysis.
  const RayBatch basis = traced(cs, b, bundle({{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}}));
  const auto h0 = local(cs, basis, 0, det);
  const auto h1 = local(cs, basis, 1, det);
  const auto h2 = local(cs, basis, 2, det);
  const double m00 = h1.first - h0.first;
  const double m10 = h1.second - h0.second;
  const double m01 = h2.first - h0.first;
  const double m11 = h2.second - h0.second;
  const double det_m = m00 * m11 - m01 * m10;
  REQUIRE(std::abs(det_m) > 0.5);
  std::vector<std::pair<double, double>> aimed;
  for (const auto& p : g.points) {
    const double u = p.x - h0.first;
    const double v = p.y - h0.second;
    aimed.emplace_back((m11 * u - m01 * v) / det_m, (-m10 * u + m00 * v) / det_m);
  }
  const RayBatch on_b = traced(cs, b, bundle(aimed));
  double shear_max = 0.0;
  for (std::size_t i = 0; i < start.size(); ++i) {
    INFO("ray " << i);
    const auto& p = g.points[i];
    REQUIRE(on_b.status()[i] == RayStatus::Alive);
    const auto [x, y] = local(cs, on_b, i, det);
    REQUIRE(std::abs(x - p.x) <= 1e-9);  // the aimed ray of b lands at x_a
    REQUIRE(std::abs(y - p.y) <= 1e-9);
    REQUIRE(std::abs(p.delta - (on_b.opl()[i] - on_a.opl()[i])) <= 1e-10);
    // Content guard: the transfer is not trivial here.
    REQUIRE(p.shear > 0.1);
    REQUIRE(p.shear < 0.4);
    shear_max = std::max(shear_max, p.shear);
    // The derived values.
    REQUIRE(std::abs(p.waves - p.delta / kLambdaMm) <= 1e-9);
    REQUIRE(std::abs(p.phase - 2.0 * std::numbers::pi * p.delta / kLambdaMm) <= 1e-8);
    REQUIRE(std::abs(p.fringe - (1.0 + std::cos(p.phase))) <= 1e-12);
    REQUIRE(std::abs(p.intensity - (0.5 + 0.5 * std::cos(p.phase))) <= 1e-12);  // w = 0.25
  }
  REQUIRE(g.max_shear == shear_max);
}

TEST_CASE("interferogram: fringe count lambda / sin(2 alpha) along a row (#210)",
          "[interferogram]") {
  // Acceptance criterion of #210: the fringe period lambda / (2 alpha) for small alpha. Along a
  // row of length L in the local y, the phase runs over Phi = 2 pi L sin(2 alpha) / lambda, so
  // cos(phase) changes its sign Phi / pi times, up to 1 at the ends. Pitch 0.005 mm, far below
  // half the period (0.16 mm): no crossing is missed.
  const double alpha = 1e-3;
  const MaterialLibrary lib;
  const CompiledSystem cs = rtt::compile::compile(tilted_michelson(alpha), lib);
  const Interferogram g = rtt::analysis::interferogram(cs, path(cs, "reference arm"),
                                                       path(cs, "test arm"), row_bundle());
  REQUIRE(g.points.size() == 901);
  int crossings = 0;
  for (std::size_t i = 1; i < g.points.size(); ++i) {
    REQUIRE(g.points[i].arrival == InterferogramArrival::Both);
    const double c0 = g.points[i - 1].fringe - 1.0;
    const double c1 = g.points[i].fringe - 1.0;
    if ((c0 < 0.0) != (c1 < 0.0)) ++crossings;
  }
  const double length = std::abs(g.points.back().y - g.points.front().y);
  REQUIRE(std::abs(length - 4.5) <= 1e-9);  // collimated 1:1 onto the camera
  const double expected = 2.0 * length * std::sin(2.0 * alpha) / kLambdaMm;
  INFO("crossings " << crossings << ", expected " << expected);
  REQUIRE(expected > 20.0);  // content guard: many fringes
  REQUIRE(std::abs(static_cast<double>(crossings) - expected) <= 1.0);
}

TEST_CASE("interferogram: aligned Michelson is the OPL difference, two wavelengths (#210)",
          "[interferogram]") {
  // Without tilt both rays of a start ray land at the same point: no shear, delta =
  // opl_difference() = 2 Delta = 15 mm for every ray (#122). A second wavelength 0.5 um: waves
  // and phase use the ray's own wavelength. Equal weights 0.25: intensity = 0.5 (1 + cos phase).
  System s = michelson();
  s.wavelengths.push_back({0.5, 1.0, false});
  const MaterialLibrary lib;
  const CompiledSystem cs = rtt::compile::compile(s, lib);
  const PathId a = path(cs, "reference arm");
  const PathId b = path(cs, "test arm");
  RayBatch start = grid_bundle();
  for (std::size_t i = 0; i < start.size(); i += 2) start.wl()[i] = 1;
  const Interferogram g = rtt::analysis::interferogram(cs, a, b, start);
  const auto d = rtt::analysis::opl_difference(cs, a, b, start);
  REQUIRE(g.points.size() == start.size());
  REQUIRE(g.max_shear <= 1e-12);
  for (std::size_t i = 0; i < start.size(); ++i) {
    INFO("ray " << i);
    const auto& p = g.points[i];
    const double lambda = (start.wl()[i] == 1 ? 0.5 : 0.6328) * 1e-3;
    REQUIRE(p.arrival == InterferogramArrival::Both);
    REQUIRE(p.shear <= 1e-12);
    REQUIRE(std::abs(p.delta - d.points[i].delta) <= 1e-12);
    REQUIRE(std::abs(p.delta - 15.0) <= 1e-10);
    REQUIRE(std::abs(p.waves - 15.0 / lambda) <= 1e-6);
    REQUIRE(std::abs(p.phase - 2.0 * std::numbers::pi * p.delta / lambda) <= 1e-8);
    REQUIRE(std::abs(p.fringe - (1.0 + std::cos(p.phase))) <= 1e-12);
    REQUIRE(std::abs(p.intensity - 0.5 * (1.0 + std::cos(p.phase))) <= 1e-12);
  }
  // Content guard: the two wavelengths give different phases.
  REQUIRE(std::abs(g.points[0].waves - g.points[1].waves) > 1000.0);
}

TEST_CASE("interferogram: unequal arms, w_a + w_b + 2 sqrt(w_a w_b) cos phase (#210)",
          "[interferogram]") {
  // An ideal splitter R = 0.5 at normal incidence in the test arm, passed twice in
  // transmission: amplitude sqrt(1 - R) each time, no phase, no change of the polarisation
  // (ADR 0021). So w_b = R T T T = 0.0625 and w_a = 0.25, both from path_transmission(); the
  // OPL is unchanged (thin element in VACUUM), delta = 15 mm. Eq. (5.12) with these weights.
  System s = michelson();
  rtt::model::Surface att;
  att.id = SurfaceId("ATT");
  att.aperture = rtt::model::CircularAperture{12.5, 0.0};
  att.interaction = rtt::model::IdealBeamSplitter{0.5, 0.5};
  rtt::model::Pose pose;  // normal to the test arm, between the splitter and the test mirror
  pose.position[1].value = 30.0;
  pose.position[2].value = 50.0;
  pose.rotation_deg[0].value = -90.0;
  s.root.children.push_back(
      {Element{"attenuator", rtt::model::ElementKind::ThinElement, pose, std::nullopt, {att}}});
  model_path(s, "test arm").events = {
      {SurfaceId("STO"), EventKind::Transmit, 0}, {SurfaceId("BS"), EventKind::Reflect, 0},
      {SurfaceId("ATT"), EventKind::Transmit, 0}, {SurfaceId("M2"), EventKind::Reflect, 0},
      {SurfaceId("ATT"), EventKind::Transmit, 0}, {SurfaceId("BS"), EventKind::Transmit, 0},
      {SurfaceId("DET"), EventKind::Transmit, 0}};
  const MaterialLibrary lib;
  const CompiledSystem cs = rtt::compile::compile(s, lib);
  const PathId a = path(cs, "reference arm");
  const PathId b = path(cs, "test arm");
  const RayBatch start = grid_bundle();
  const auto ta = rtt::analysis::path_transmission(cs, a, start);
  const auto tb = rtt::analysis::path_transmission(cs, b, start);
  const Interferogram g = rtt::analysis::interferogram(cs, a, b, start);
  REQUIRE(g.points.size() == start.size());
  for (std::size_t i = 0; i < start.size(); ++i) {
    INFO("ray " << i);
    const double wa = ta.rays[i].weight;
    const double wb = tb.rays[i].weight;
    REQUIRE(std::abs(wa - 0.25) <= 1e-12);
    REQUIRE(std::abs(wb - 0.0625) <= 1e-12);  // content guard: w_a != w_b
    const auto& p = g.points[i];
    REQUIRE(p.arrival == InterferogramArrival::Both);
    REQUIRE(std::abs(p.delta - 15.0) <= 1e-10);
    REQUIRE(std::abs(p.intensity - (wa + wb + 2.0 * std::sqrt(wa * wb) * std::cos(p.phase))) <=
            1e-12);
    REQUIRE(std::abs(p.fringe - (1.0 + std::cos(p.phase))) <= 1e-12);
  }
}

TEST_CASE("interferogram: rays on one path only, statuses and physical intensity (#210)",
          "[interferogram]") {
  // Reference mirror of radius 1 mm: a loses the start rays with r > 1 (Vignetted). Test mirror
  // 1.2 mm wide in x (half width 0.6; the splitter keeps x): b loses |x| > 0.6. On the grid of
  // pitch 0.5 all four cases occur, e.g. (0, 1) Both, (1, 0) OnlyA, (0.5, 1) OnlyB, (1, 0.5)
  // None. Not Both: delta, waves, phase, fringe, shear 0; intensity w_a, w_b or 0; x, y of the
  // arrived path (of b for OnlyB), 0 for None.
  System s = michelson();
  element(s, "reference mirror").surfaces[0].aperture = rtt::model::CircularAperture{1.0, 0.0};
  element(s, "test mirror").surfaces[0].aperture = rtt::model::RectangularAperture{0.6, 12.5};
  const MaterialLibrary lib;
  const CompiledSystem cs = rtt::compile::compile(s, lib);
  const PathId a = path(cs, "reference arm");
  const PathId b = path(cs, "test arm");
  const std::uint32_t det = surface(cs, "DET");
  const RayBatch start = grid_bundle();
  const Interferogram g = rtt::analysis::interferogram(cs, a, b, start);
  const RayBatch on_a = traced(cs, a, start);
  const RayBatch on_b = traced(cs, b, start);
  REQUIRE(g.points.size() == start.size());
  std::array<std::size_t, 4> count{};  // Both, OnlyA, OnlyB, None
  for (std::size_t i = 0; i < start.size(); ++i) {
    INFO("ray " << i);
    const double x0 = start.pos_x()[i];
    const double y0 = start.pos_y()[i];
    const bool in_a = x0 * x0 + y0 * y0 <= 1.0;
    const bool in_b = std::abs(x0) <= 0.6;
    const auto& p = g.points[i];
    REQUIRE(p.status_a == (in_a ? RayStatus::Alive : RayStatus::Vignetted));
    REQUIRE(p.status_b == (in_b ? RayStatus::Alive : RayStatus::Vignetted));
    if (in_a && in_b) {
      ++count[0];
      REQUIRE(p.arrival == InterferogramArrival::Both);
      REQUIRE(std::abs(p.delta - 15.0) <= 1e-10);
      REQUIRE(std::abs(p.intensity - 0.5 * (1.0 + std::cos(p.phase))) <= 1e-12);
      continue;
    }
    REQUIRE(p.delta == 0.0);
    REQUIRE(p.waves == 0.0);
    REQUIRE(p.phase == 0.0);
    REQUIRE(p.fringe == 0.0);
    REQUIRE(p.shear == 0.0);
    if (in_a) {
      ++count[1];
      REQUIRE(p.arrival == InterferogramArrival::OnlyA);
      REQUIRE(p.intensity == on_a.weight()[i]);
      REQUIRE(std::abs(p.intensity - 0.25) <= 1e-12);
      const auto [x, y] = local(cs, on_a, i, det);
      REQUIRE(p.x == x);
      REQUIRE(p.y == y);
    } else if (in_b) {
      ++count[2];
      REQUIRE(p.arrival == InterferogramArrival::OnlyB);
      REQUIRE(p.intensity == on_b.weight()[i]);
      REQUIRE(std::abs(p.intensity - 0.25) <= 1e-12);
      const auto [x, y] = local(cs, on_b, i, det);
      REQUIRE(p.x == x);
      REQUIRE(p.y == y);
      REQUIRE(std::abs(p.y) > 0.9);  // content guard: not at the origin
    } else {
      ++count[3];
      REQUIRE(p.arrival == InterferogramArrival::None);
      REQUIRE(p.intensity == 0.0);
      REQUIRE(p.x == 0.0);
      REQUIRE(p.y == 0.0);
    }
  }
  for (const std::size_t c : count) REQUIRE(c > 0);
  REQUIRE(g.losses_a.count(RayStatus::Vignetted) == count[2] + count[3]);
  REQUIRE(g.losses_b.count(RayStatus::Vignetted) == count[1] + count[3]);
  bool lost_warning = false;
  for (const auto& w : g.warnings) lost_warning = lost_warning || w.code == "rays.lost";
  REQUIRE(lost_warning);
}

TEST_CASE("interferogram: convenience form, threads, run control (#210)", "[interferogram]") {
  const MaterialLibrary lib;
  // Convenience form on a symmetric path, a path with itself: delta 0, fringe 2, intensity
  // 4 w; bitwise the main form with the rays of make_rays.
  System singlet = load("m1/singlet_const.rtt.json");
  singlet.environment.medium = "VACUUM";
  const CompiledSystem cs1 = rtt::compile::compile(singlet, lib);
  const rtt::analysis::PathOptions options;
  const Interferogram convenient = rtt::analysis::interferogram(cs1, PathId{0}, PathId{0}, 1, 0);
  const std::vector<std::uint16_t> field{1};
  const RayBatch rays =
      rtt::trace::make_rays(cs1, PathId{0}, field, 0, options.sampling, options.aiming);
  const Interferogram main = rtt::analysis::interferogram(cs1, PathId{0}, PathId{0}, rays);
  REQUIRE(convenient.points.size() == rays.size());
  REQUIRE(same_bits(convenient, main));
  const auto t = rtt::analysis::path_transmission(cs1, PathId{0}, rays);
  std::size_t both = 0;
  for (std::size_t i = 0; i < rays.size(); ++i) {
    const auto& p = main.points[i];
    if (p.arrival != InterferogramArrival::Both) continue;
    ++both;
    REQUIRE(p.delta == 0.0);
    REQUIRE(p.fringe == 2.0);
    REQUIRE(std::abs(p.intensity - 4.0 * t.rays[i].weight) <= 1e-12);
  }
  REQUIRE(both > 10);
  const CompiledSystem folded = rtt::compile::compile(michelson(), lib);
  REQUIRE_THROWS_AS(rtt::analysis::interferogram(folded, PathId{0}, PathId{1}, 0, 0),
                    rtt::paraxial::ParaxialError);

  // Threads: bitwise at 1, 4 and all threads (the tracer runs in parallel, the analysis is
  // serial in batch order).
  const CompiledSystem cs = rtt::compile::compile(tilted_michelson(1e-3), lib);
  const PathId a = path(cs, "reference arm");
  const PathId b = path(cs, "test arm");
  const RayBatch start = grid_bundle();
  const Interferogram all = rtt::analysis::interferogram(cs, a, b, start);
  for (const int threads : {1, 4}) {
    INFO(threads << " threads");
    const tbb::global_control limit(tbb::global_control::max_allowed_parallelism,
                                    static_cast<std::size_t>(threads));
    CHECK(same_bits(rtt::analysis::interferogram(cs, a, b, start), all));
  }

  // Run control: the same result, and a cancellation.
  rtt::trace::RunControl control;
  control.cancel = rtt::trace::CancelToken();
  control.block_size = 7;
  std::vector<std::string> stages;
  control.progress = [&](const rtt::trace::Progress& p) { stages.emplace_back(p.stage); };
  REQUIRE(same_bits(rtt::analysis::interferogram(cs, a, b, start, options, control), all));
  REQUIRE(!stages.empty());
  for (const auto& st : stages) REQUIRE(st == "trace");
  control.cancel->request_cancel();
  REQUIRE_THROWS_AS(rtt::analysis::interferogram(cs, a, b, start, options, control),
                    rtt::trace::Cancelled);
}

TEST_CASE("interferogram: input errors (#210)", "[interferogram]") {
  const MaterialLibrary lib;
  const CompiledSystem cs = rtt::compile::compile(michelson(), lib);
  const PathId a = path(cs, "reference arm");
  const PathId b = path(cs, "test arm");
  const RayBatch start = grid_bundle();
  // Different image surfaces (camera and stop).
  REQUIRE_THROWS_AS(rtt::analysis::interferogram(cs, a, path(cs, "reference arm, return"), start),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(rtt::analysis::interferogram(cs, a, b, RayBatch(0)), std::invalid_argument);
  REQUIRE_THROWS_AS(rtt::analysis::interferogram(cs, a, PathId{99}, start), std::invalid_argument);
  rtt::analysis::PathOptions bad;
  bad.lost_warning_fraction = 2.0;
  REQUIRE_THROWS_AS(rtt::analysis::interferogram(cs, a, b, start, bad), std::invalid_argument);

  // A crystal as the image medium: the o and e modes of the calcite plate of
  // m4/calcite_walkoff, both ending at the entrance face P.S1.
  System calcite = load("m4/calcite_walkoff.rtt.json");
  calcite.paths.push_back({"o in", false, {{SurfaceId("P.S1"), EventKind::Ordinary, 0}}});
  calcite.paths.push_back({"e in", false, {{SurfaceId("P.S1"), EventKind::Extraordinary, 0}}});
  const CompiledSystem cc = rtt::compile::compile(calcite, lib);
  RayBatch one(1);
  one.pos_y()[0] = 0.3;
  REQUIRE_THROWS_AS(rtt::analysis::interferogram(cc, path(cc, "o in"), path(cc, "e in"), one),
                    std::invalid_argument);
  // The same image surface, but different image media (review of #241): refracted into the lens
  // of m1/singlet_const and reflected back into the VACUUM at L1.S1. The waves would not meet
  // in one medium, so it is an error that names both paths.
  System singlet = load("m1/singlet_const.rtt.json");
  singlet.environment.medium = "VACUUM";
  singlet.paths.push_back(
      {"into the lens",
       false,
       {{SurfaceId("STO"), EventKind::Transmit, 0}, {SurfaceId("L1.S1"), EventKind::Refract, 0}}});
  singlet.paths.push_back(
      {"reflected",
       false,
       {{SurfaceId("STO"), EventKind::Transmit, 0}, {SurfaceId("L1.S1"), EventKind::Reflect, 0}}});
  const CompiledSystem cl = rtt::compile::compile(singlet, lib);
  const PathId into = path(cl, "into the lens");
  const PathId back = path(cl, "reflected");
  RayBatch ahead(1);
  ahead.pos_z()[0] = -1.0;  // in front of the stop (AGENTS.md)
  std::string message;
  try {
    [[maybe_unused]] const Interferogram bad_media =
        rtt::analysis::interferogram(cl, into, back, ahead);
  } catch (const std::invalid_argument& e) {
    message = e.what();
  }
  INFO(message);
  REQUIRE(message.find("/paths/" + std::to_string(into.index)) != std::string::npos);
  REQUIRE(message.find("/paths/" + std::to_string(back.index)) != std::string::npos);
  REQUIRE(message.find("image medium") != std::string::npos);
  // Through the plate into the vacuum behind it: allowed (the image medium is isotropic).
  const Interferogram oe = rtt::analysis::interferogram(cc, path(cc, "o"), path(cc, "e"), one);
  REQUIRE(oe.points.size() == 1);
  REQUIRE(oe.points[0].arrival == InterferogramArrival::Both);
}
