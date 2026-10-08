#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
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
#include "rtt/model/model.hpp"
#include "rtt/paraxial/paraxial.hpp"
#include "rtt/trace/run_control.hpp"
#include "rtt/trace/sources.hpp"

// Path evaluation (#122): transmission per path and OPL difference of two paths. Reference
// cases are analytic; the derivations are in the test comments.

using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::material::MaterialLibrary;
using rtt::model::Element;
using rtt::model::System;
using rtt::trace::RayBatch;
using rtt::trace::RayStatus;

namespace {

System load(const std::string& relative) {
  return rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/" + relative);
}

/// tests/reference/m4/michelson_offset.rtt.json, VACUUM: stop at z = 10, beam splitter
/// (ideal, R = 0.5 for s and p) at z = 50 under 45 deg, reference mirror 50 mm behind it at
/// z = 100, test mirror at y = 57.5, i.e. Delta = 7.5 mm farther; camera at y = -50. Paths
/// "reference arm" (T, then R) and "test arm" (R, then T) to the camera, "reference arm,
/// return" (T, T) and "test arm, return" (R, R) back through the stop.
System michelson() {
  return load("m4/michelson_offset.rtt.json");
}

Element& element(System& s, const std::string& name) {
  for (auto& node : s.root.children) {
    if (auto* e = std::get_if<Element>(&node.value); e && e->name == name) return *e;
  }
  throw std::logic_error("no element " + name);
}

PathId path(const CompiledSystem& cs, const std::string& name) {
  const std::optional<PathId> id = cs.find_path(name);
  REQUIRE(id.has_value());
  return *id;
}

/// Collimated start bundle along +z at z = 0: the axis ray first, then a square grid with a
/// pitch of 0.5 mm inside r <= 2.25 mm (within the EPD of 5 mm). Pupil coordinates = r / 2.5.
RayBatch collimated_bundle() {
  std::vector<std::pair<double, double>> points{{0.0, 0.0}};
  for (int i = -4; i <= 4; ++i) {
    for (int j = -4; j <= 4; ++j) {
      const double x = 0.5 * i;
      const double y = 0.5 * j;
      if ((i != 0 || j != 0) && x * x + y * y <= 2.25 * 2.25) points.emplace_back(x, y);
    }
  }
  RayBatch rays(points.size());
  for (std::size_t k = 0; k < points.size(); ++k) {
    rays.pos_x()[k] = points[k].first;
    rays.pos_y()[k] = points[k].second;
    rays.pupil_x()[k] = points[k].first / 2.5;
    rays.pupil_y()[k] = points[k].second / 2.5;
  }
  return rays;
}

}  // namespace

TEST_CASE("Michelson: arm weights R T and the balance over both output ports (#122)", "[paths]") {
  // Ideal beam splitter with R = 0.5 for s and p: power-normalized amplitudes sqrt(R) on
  // Reflect and sqrt(1 - R) on Transmit (ADR 0021), ideal mirrors lossless. weight is the
  // power for an unpolarized source, so each camera arm carries T R = 0.25, the return paths
  // T T and R R = 0.25. Per ray the four paths add up to (R + T)^2 = 1: all power leaves
  // through the two output ports. Products of exact values: 1e-12.
  const MaterialLibrary lib;
  const CompiledSystem cs = rtt::compile::compile(michelson(), lib);
  const RayBatch start = collimated_bundle();
  std::vector<double> sum(start.size(), 0.0);
  for (const std::string name :
       {"reference arm", "test arm", "reference arm, return", "test arm, return"}) {
    INFO(name);
    const auto t = rtt::analysis::path_transmission(cs, path(cs, name), start);
    REQUIRE(t.rays_launched == start.size());
    REQUIRE(t.rays_arrived == start.size());
    REQUIRE(t.rays.size() == start.size());
    REQUIRE(std::abs(t.mean - 0.25) <= 1e-12);
    REQUIRE(std::abs(t.min - 0.25) <= 1e-12);
    REQUIRE(std::abs(t.max - 0.25) <= 1e-12);
    REQUIRE(t.losses.launched == start.size());
    for (std::size_t i = 0; i < start.size(); ++i) {
      REQUIRE(t.rays[i].status == RayStatus::Alive);
      REQUIRE(t.rays[i].px == start.pupil_x()[i]);
      REQUIRE(t.rays[i].py == start.pupil_y()[i]);
      sum[i] += t.rays[i].weight;
    }
  }
  for (std::size_t i = 0; i < start.size(); ++i) REQUIRE(std::abs(sum[i] - 1.0) <= 1e-12);
}

TEST_CASE("Michelson: polarizing splitter, the mean of the products (#122)", "[paths]") {
  // R_s = 0.3, R_p = 0.6. The axis ray meets the splitter twice in the same plane of incidence
  // (the splitter is rotated about x, the arms lie in the y-z plane) and the mirrors at normal
  // incidence, so s and p do not mix: the reference arm carries 1/2 (T_s R_s + T_p R_p) =
  // 1/2 (0.21 + 0.24) = 0.225, not the product of the means 0.45 * 0.55 = 0.2475. Return
  // paths: 1/2 (T_s^2 + T_p^2) = 0.325 and 1/2 (R_s^2 + R_p^2) = 0.225; the sum of all four is
  // 1/2 ((R_s + T_s)^2 + (R_p + T_p)^2) = 1.
  System s = michelson();
  element(s, "beam splitter").surfaces[0].interaction = rtt::model::IdealBeamSplitter{0.3, 0.6};
  const MaterialLibrary lib;
  const CompiledSystem cs = rtt::compile::compile(s, lib);
  const RayBatch axis(1);  // at the origin along +z, wavelength 0, Alive
  const auto weight = [&](const std::string& name) {
    const auto t = rtt::analysis::path_transmission(cs, path(cs, name), axis);
    REQUIRE(t.rays.size() == 1);
    REQUIRE(t.rays_arrived == 1);
    return t.rays[0].weight;
  };
  const double ref = weight("reference arm");
  const double test = weight("test arm");
  const double ref_return = weight("reference arm, return");
  const double test_return = weight("test arm, return");
  REQUIRE(std::abs(ref - 0.225) <= 1e-12);
  REQUIRE(std::abs(test - 0.225) <= 1e-12);
  REQUIRE(std::abs(ref_return - 0.325) <= 1e-12);
  REQUIRE(std::abs(test_return - 0.225) <= 1e-12);
  REQUIRE(std::abs(ref + test + ref_return + test_return - 1.0) <= 1e-12);
}

TEST_CASE("Michelson: the OPL difference of the arms is twice the arm difference (#122)",
          "[paths]") {
  // Both arms share the way to the splitter and from the splitter to the camera; every ray
  // returns from its plane mirror at normal incidence to the point where it left the splitter.
  // The arms differ only in splitter -> mirror -> splitter, and the test mirror is Delta =
  // 7.5 mm farther: OPL(test) - OPL(reference) = 2 Delta = 15 mm for every ray (VACUUM).
  // Rounding of OPL sums of about 200 mm: 1e-10 mm.
  const MaterialLibrary lib;
  const CompiledSystem cs = rtt::compile::compile(michelson(), lib);
  const RayBatch start = collimated_bundle();
  const auto d =
      rtt::analysis::opl_difference(cs, path(cs, "reference arm"), path(cs, "test arm"), start);
  REQUIRE(d.points.size() == start.size());
  REQUIRE(d.chief.has_value());
  REQUIRE(std::abs(*d.chief - 15.0) <= 1e-10);
  for (const auto& p : d.points) {
    REQUIRE(p.status == RayStatus::Alive);
    REQUIRE(std::abs(p.delta - 15.0) <= 1e-10);
  }
  REQUIRE(d.losses_a.launched == start.size());
  REQUIRE(d.losses_b.launched == start.size());
}

TEST_CASE("path transmission: lost rays count as launched with weight 0 (#122)", "[paths]") {
  // Test mirror with radius 1 mm: the splitter maps the start ray at (x, y) to (x, -y) on the
  // mirror's local plane (45 deg fold, mirror rotated by -90 deg about x), so the rays with
  // x^2 + y^2 > 1 end Vignetted at M2.
  System s = michelson();
  element(s, "test mirror").surfaces[0].aperture = rtt::model::CircularAperture{1.0, 0.0};
  const MaterialLibrary lib;
  const CompiledSystem cs = rtt::compile::compile(s, lib);
  const RayBatch start = collimated_bundle();
  std::size_t inside = 0;
  for (std::size_t i = 0; i < start.size(); ++i) {
    const double r2 = start.pos_x()[i] * start.pos_x()[i] + start.pos_y()[i] * start.pos_y()[i];
    if (r2 <= 1.0) ++inside;
  }
  REQUIRE(2 * inside < start.size());  // more than half lost: the warning must come
  const auto t = rtt::analysis::path_transmission(cs, path(cs, "test arm"), start);
  REQUIRE(t.rays_arrived == inside);
  REQUIRE(t.losses.count(RayStatus::Vignetted) == start.size() - inside);
  REQUIRE(t.losses.worst_surface == cs.find_surface(rtt::model::SurfaceId("M2")));
  REQUIRE(std::abs(t.mean - 0.25 * static_cast<double>(inside) /
                                static_cast<double>(start.size())) <= 1e-12);
  for (const auto& r : t.rays) {
    if (r.status != RayStatus::Alive) REQUIRE(r.weight == 0.0);
  }
  bool lost_warning = false;
  for (const auto& w : t.warnings) lost_warning = lost_warning || w.code == "rays.lost";
  REQUIRE(lost_warning);
}

TEST_CASE("path transmission: start rays that are not Alive, start weight, wavelength (#122)",
          "[paths]") {
  const MaterialLibrary lib;
  const CompiledSystem cs = rtt::compile::compile(michelson(), lib);
  RayBatch start = collimated_bundle();
  start.status()[1] = RayStatus::NoConvergence;  // e.g. from make_rays: launched and lost
  start.weight()[2] = 0.5;                       // a source apodization is carried along
  const auto t = rtt::analysis::path_transmission(cs, path(cs, "reference arm"), start);
  REQUIRE(t.rays_launched == start.size());
  REQUIRE(t.rays.size() == start.size());
  REQUIRE(t.rays_arrived == start.size() - 1);
  REQUIRE(t.rays[1].status == RayStatus::NoConvergence);
  REQUIRE(t.rays[1].weight == 0.0);
  REQUIRE(t.losses.count(RayStatus::NoConvergence) == 1);
  REQUIRE(std::abs(t.rays[2].weight - 0.125) <= 1e-12);

  RayBatch wrong = collimated_bundle();
  wrong.wl()[0] = 5;  // not a system wavelength
  REQUIRE_THROWS_AS(rtt::analysis::path_transmission(cs, PathId{0}, wrong), std::invalid_argument);
  REQUIRE_THROWS_AS(rtt::analysis::path_transmission(cs, PathId{0}, RayBatch{}),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(rtt::analysis::path_transmission(cs, PathId{9}, start), std::invalid_argument);
  // The OPL difference needs one image surface: camera against the return port.
  REQUIRE_THROWS_AS(rtt::analysis::opl_difference(cs, path(cs, "reference arm"),
                                                  path(cs, "test arm, return"), start),
                    std::invalid_argument);
  rtt::analysis::PathOptions bad;
  bad.lost_warning_fraction = 1.5;
  REQUIRE_THROWS_AS(rtt::analysis::path_transmission(cs, PathId{0}, start, bad),
                    std::invalid_argument);
}

TEST_CASE("path transmission, convenience form: make_rays on a symmetric path (#122)", "[paths]") {
  // Singlet CONST:1.5168 in VACUUM: the axis ray meets both surfaces at normal incidence, so
  // its weight is T^2 with T = 1 - ((n - 1) / (n + 1))^2: Byrnes, arXiv:1603.02720v5, Eq. (6)
  // at normal incidence, T = 1 - R for a lossless interface by Eqs. (21)-(23) (docs/quellen.md).
  // 1e-12.
  System s = load("m1/singlet_const.rtt.json");
  s.environment.medium = "VACUUM";
  const MaterialLibrary lib;
  const CompiledSystem cs = rtt::compile::compile(s, lib);
  rtt::analysis::PathOptions axis;
  axis.sampling = rtt::trace::SinglePupilPoint{0.0, 0.0};
  const auto t = rtt::analysis::path_transmission(cs, PathId{0}, 0, 0, axis);
  REQUIRE(t.rays.size() == 1);
  REQUIRE(t.rays_arrived == 1);
  const double n = 1.5168;
  const double r = (n - 1.0) / (n + 1.0);
  REQUIRE(std::abs(t.rays[0].weight - (1.0 - r * r) * (1.0 - r * r)) <= 1e-12);

  // Bitwise the main form with the rays of make_rays.
  const rtt::analysis::PathOptions options;
  const auto convenient = rtt::analysis::path_transmission(cs, PathId{0}, 1, 0, options);
  const std::vector<std::uint16_t> field{1};
  const RayBatch start =
      rtt::trace::make_rays(cs, PathId{0}, field, 0, options.sampling, options.aiming);
  const auto main = rtt::analysis::path_transmission(cs, PathId{0}, start, options);
  REQUIRE(convenient.rays.size() == main.rays.size());
  for (std::size_t i = 0; i < main.rays.size(); ++i) {
    REQUIRE(convenient.rays[i].weight == main.rays[i].weight);
    REQUIRE(convenient.rays[i].status == main.rays[i].status);
  }
  REQUIRE(convenient.mean == main.mean);

  // The paraxial aiming does not accept the folded Michelson: use the main form there.
  const CompiledSystem folded = rtt::compile::compile(michelson(), lib);
  REQUIRE_THROWS_AS(rtt::analysis::path_transmission(folded, PathId{0}, 0, 0),
                    rtt::paraxial::ParaxialError);
}

TEST_CASE("path evaluation with a run control (#122)", "[paths][run_control]") {
  const MaterialLibrary lib;
  const CompiledSystem cs = rtt::compile::compile(michelson(), lib);
  const RayBatch start = collimated_bundle();
  const PathId ref = path(cs, "reference arm");
  const PathId test = path(cs, "test arm");
  rtt::trace::RunControl control;
  control.cancel = rtt::trace::CancelToken();
  control.block_size = 7;
  std::vector<std::string> stages;
  control.progress = [&](const rtt::trace::Progress& p) { stages.emplace_back(p.stage); };
  const rtt::analysis::PathOptions options;
  const auto plain = rtt::analysis::path_transmission(cs, ref, start);
  const auto controlled = rtt::analysis::path_transmission(cs, ref, start, options, control);
  REQUIRE(plain.rays.size() == start.size());
  REQUIRE(controlled.rays.size() == start.size());
  for (std::size_t i = 0; i < start.size(); ++i) {
    REQUIRE(plain.rays[i].weight == controlled.rays[i].weight);
  }
  REQUIRE(!stages.empty());
  for (const auto& s : stages) REQUIRE(s == "trace");
  const auto d_plain = rtt::analysis::opl_difference(cs, ref, test, start);
  const auto d_controlled = rtt::analysis::opl_difference(cs, ref, test, start, options, control);
  REQUIRE(d_plain.points.size() == start.size());
  REQUIRE(d_controlled.points.size() == start.size());
  for (std::size_t i = 0; i < start.size(); ++i) {
    REQUIRE(d_plain.points[i].delta == d_controlled.points[i].delta);
  }
  control.cancel->request_cancel();
  REQUIRE_THROWS_AS(rtt::analysis::path_transmission(cs, ref, start, options, control),
                    rtt::trace::Cancelled);
  REQUIRE_THROWS_AS(rtt::analysis::opl_difference(cs, ref, test, start, options, control),
                    rtt::trace::Cancelled);
}
