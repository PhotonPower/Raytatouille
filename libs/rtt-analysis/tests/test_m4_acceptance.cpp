// M4 acceptance (#135): the criteria of the roadmap row M4 Multi-Path, each through the tracer
// with a reference system under tests/reference/m4/ (and the feature tour, m0), analytic
// expectations and tolerances fixed beforehand. Sources (docs/quellen.md):
// - S. J. Byrnes, Multilayer optical calculations, arXiv:1603.02720v5, Eq. (6) and
//   Eqs. (21)-(23): R = ((n - 1) / (n + 1))^2, T = 1 - R at normal incidence;
// - C. Palmer, Diffraction Grating Handbook, 7th ed., Newport (2014), Eq. (2-1) and Fig. 2-1:
//   m lambda = d (sin alpha + sin beta), t_x = -sin alpha, t'_x = sin beta;
// - M. Mansuripur, Proc. SPIE 6620, 66200N (2007), Eq. (7b): n2 sigma' = n1 sigma +
//   m lambda0 grad F, F = phi / (2 pi);
// - J. E. Greivenkamp, OPTI-201/202, Sec. 9, p. 9-2: y-nu trace n'u' = nu - y phi.
// - W.-S. T. Lam, Anisotropic Ray Trace, dissertation, Eq. (2.39) for n_e(theta); the walk-off
//   at 45 degree, tan(rho) = (n_O^2 - n_E^2) / (n_O^2 + n_E^2), is derived in #130 from Lam,
//   Eq. (2.41) and p. 107 (docs/quellen.md).

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <numbers>
#include <optional>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "rtt/analysis/ghosts.hpp"
#include "rtt/analysis/paths.hpp"
#include "rtt/coating/catalog.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/compile/ghosts.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/math/types.hpp"
#include "rtt/model/model.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/ray_paths.hpp"
#include "rtt/trace/sequential.hpp"

using rtt::compile::CompiledSystem;
using rtt::compile::GhostSystem;
using rtt::compile::PathId;
using rtt::math::Vec3;
using rtt::trace::RayBatch;
using rtt::trace::RayStatus;

namespace {

constexpr double kPi = std::numbers::pi;

rtt::model::System load(const std::string& relative) {
  return rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/" + relative);
}

PathId path(const CompiledSystem& cs, const std::string& name) {
  const auto id = cs.find_path(name);
  REQUIRE(id.has_value());
  return *id;
}

bool close(double a, double b, double rel) {
  return std::abs(a - b) <= rel * std::max(std::abs(a), std::abs(b));
}

/// Collimated bundle along +z at z = 0: the axis ray first, then a square grid with a pitch of
/// 0.5 mm inside r <= 2.25 mm; pupil labels r / 2.5 (as in libs/rtt-py/tests/paths_cases.cpp).
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

TEST_CASE("M4 Michelson: R + T = 1, arm weights R T, OPL difference 2 Delta", "[m4]") {
  // tests/reference/m4/michelson_offset.rtt.json, VACUUM: ideal beam splitter R = 0.5 for s and
  // p (lossless, ADR 0021), test mirror Delta = 7.5 mm farther than the reference mirror. Each of
  // the four paths to the camera and back carries R T, T R, T T or R R = 0.25 per ray, so the
  // four add up to (R + T)^2 = 1 per start ray; the test arm is 2 Delta = 15 mm longer.
  // Tolerances: weights are a few products of amplitudes, 1e-12; OPL sums of about 200 mm,
  // rounding about 1e-13 per sum, 1e-10.
  const CompiledSystem cs =
      rtt::compile::compile(load("m4/michelson_offset.rtt.json"), rtt::material::MaterialLibrary{});
  const RayBatch start = collimated_bundle();
  std::vector<double> total(start.size(), 0.0);
  for (const char* name :
       {"reference arm", "test arm", "reference arm, return", "test arm, return"}) {
    INFO(name);
    const auto t = rtt::analysis::path_transmission(cs, path(cs, name), start);
    REQUIRE(t.rays_arrived == start.size());
    REQUIRE(std::abs(t.mean - 0.25) <= 1e-12);
    for (std::size_t i = 0; i < start.size(); ++i) {
      REQUIRE(std::abs(t.rays[i].weight - 0.25) <= 1e-12);
      total[i] += t.rays[i].weight;
    }
  }
  for (const double w : total) REQUIRE(std::abs(w - 1.0) <= 1e-12);
  const auto d =
      rtt::analysis::opl_difference(cs, path(cs, "reference arm"), path(cs, "test arm"), start);
  REQUIRE(d.chief.has_value());
  REQUIRE(std::abs(*d.chief - 15.0) <= 1e-10);
  for (const auto& p : d.points) {
    REQUIRE(p.status == RayStatus::Alive);
    REQUIRE(std::abs(p.delta - 15.0) <= 1e-10);
  }
}

TEST_CASE("M4 calcite walk-off: offset t tan(rho) of the e-ray, none for the o-ray", "[m4]") {
  // tests/reference/m4/calcite_walkoff.rtt.json: calcite plate (CONST n_O = 1.6584,
  // n_E = 1.4864, ADR 0026 acceptance), t = 2 mm from z = 10 to 12, optic axis (1, 0, 1) at
  // 45 degree in the x-z plane, detector at z = 20, vacuum, normal incidence. The o-ray goes
  // straight through; the e-ray leaves the plate parallel to the incidence, shifted by
  // t tan(rho) away from the axis direction +x, with tan(rho) = (n_O^2 - n_E^2) /
  // (n_O^2 + n_E^2) (Lam, p. 107). OPL: 18 mm in vacuum plus n t with n_O or
  // n_e(45 deg) = sqrt(2 / (1 / n_O^2 + 1 / n_E^2)) (Lam, Eq. (2.39); in the plate
  // l = t / cos(rho) along S and k . S = cos(rho)). Rounding of a few operations on mm values:
  // 1e-10 (acceptance criterion of #135).
  const double n_o = 1.6584;
  const double n_e = 1.4864;
  const double thickness = 2.0;
  const double tan_rho = (n_o * n_o - n_e * n_e) / (n_o * n_o + n_e * n_e);
  const double n_e45 = std::sqrt(2.0 / (1.0 / (n_o * n_o) + 1.0 / (n_e * n_e)));
  const CompiledSystem cs =
      rtt::compile::compile(load("m4/calcite_walkoff.rtt.json"), rtt::material::MaterialLibrary{});
  for (const char* name : {"o", "e"}) {
    INFO(name);
    const bool e = std::string(name) == "e";
    RayBatch rays(1);
    rays.pos_y()[0] = 0.3;
    [[maybe_unused]] const auto stats =
        rtt::trace::SequentialTracer().trace(cs, path(cs, name), rays);
    REQUIRE(rays.status()[0] == RayStatus::Alive);
    REQUIRE(std::abs(rays.dir_x()[0]) <= 1e-12);
    REQUIRE(std::abs(rays.dir_y()[0]) <= 1e-12);
    REQUIRE(std::abs(rays.dir_z()[0] - 1.0) <= 1e-12);
    REQUIRE(std::abs(rays.pos_x()[0] - (e ? -thickness * tan_rho : 0.0)) <= 1e-10);
    REQUIRE(std::abs(rays.pos_y()[0] - 0.3) <= 1e-10);
    REQUIRE(std::abs(rays.pos_z()[0] - 20.0) <= 1e-10);
    REQUIRE(std::abs(rays.opl()[0] - (18.0 + (e ? n_e45 : n_o) * thickness)) <= 1e-10);
  }
  // The offset is not trivially small: about 0.218 mm.
  REQUIRE(thickness * tan_rho > 0.2);
}

TEST_CASE("M4 grating equation: orders -1, 0, +1 by Palmer (2-1), order +6 evanescent", "[m4]") {
  // tests/reference/m4/grating_transmission.rtt.json: thin grating 300 lines/mm, grooves along
  // y, vacuum, lambda = 0.5876 um. A ray at alpha = 10 deg in the x-z plane (t_x = -sin alpha,
  // Palmer, Fig. 2-1b): sin beta = m lambda G - sin alpha. Order +6 at normal incidence has
  // tau = 6 lambda G = 1.058 > n' = 1 and ends at the grating with Evanescent (ADR 0025,
  // point 7). A few operations on numbers of order 1: 1e-12.
  const CompiledSystem cs = rtt::compile::compile(load("m4/grating_transmission.rtt.json"),
                                                  rtt::material::MaterialLibrary{});
  const double lambda_mm = 0.5876e-3;
  const double alpha = 10.0 * kPi / 180.0;
  for (const auto& [name, m] :
       {std::pair{"order -1", -1}, std::pair{"order 0", 0}, std::pair{"order +1", 1}}) {
    INFO(name);
    RayBatch rays(1);
    rays.pos_z()[0] = -1.0;
    rays.dir_x()[0] = -std::sin(alpha);
    rays.dir_z()[0] = std::cos(alpha);
    [[maybe_unused]] const auto stats =
        rtt::trace::SequentialTracer().trace(cs, path(cs, name), rays);
    REQUIRE(rays.status()[0] == RayStatus::Alive);
    const double sin_beta = m * lambda_mm * 300.0 - std::sin(alpha);
    REQUIRE(std::abs(rays.dir_x()[0] - sin_beta) <= 1e-12);
    REQUIRE(rays.dir_y()[0] == 0.0);
    REQUIRE(std::abs(rays.dir_z()[0] - std::sqrt(1.0 - sin_beta * sin_beta)) <= 1e-12);
  }
  RayBatch normal(1);
  normal.pos_z()[0] = -1.0;
  [[maybe_unused]] const auto stats =
      rtt::trace::SequentialTracer().trace(cs, path(cs, "order +6"), normal);
  REQUIRE(normal.status()[0] == RayStatus::Evanescent);
  REQUIRE(cs.surfaces()[normal.last_surface()[0]].id.str() == "G");
  REQUIRE(std::abs(normal.pos_z()[0] - 10.0) <= 1e-12);
}

TEST_CASE("M4 ghost ranking of two plates: rho is the Fresnel product, in its order", "[m4]") {
  // tests/reference/m4/ghost_plates.rtt.json: plates A (n = 1.5) and B (n = 2.0) at normal
  // incidence in a collimated beam. Every ghost leaves collimated on the rays of the useful
  // image, so r_g = r_b and rho = P_g / P_b = R_i R_j prod_{i<k<j} T_k^2 (the surfaces between
  // i and j are passed three times instead of once), with R = ((n - 1) / (n + 1))^2,
  // T = 1 - R (Byrnes): R_A = 0.04, R_B = 1/9. The hexapolar bundle of 6 rings (127 rays,
  // radius 2k/6 mm) gives r_b^2 = 294/127 (equal weights, centroid 0). Relative 1e-10.
  const rtt::material::MaterialLibrary lib;
  const rtt::coating::CoatingLibrary coatings;
  const GhostSystem g =
      rtt::compile::compile_with_ghosts(load("m4/ghost_plates.rtt.json"), "main", lib, coatings);
  REQUIRE(g.ghosts.size() == 6);
  const double ra = 0.04;
  const double rb = 1.0 / 9.0;
  const double ta = 1.0 - ra;
  const double tb = 1.0 - rb;
  const std::map<std::string, double> expected{
      {"main ghost A.S2/A.S1", ra * ra},
      {"main ghost B.S1/A.S1", ra * rb * ta * ta},
      {"main ghost B.S1/A.S2", ra * rb},
      {"main ghost B.S2/A.S1", ra * rb * ta * ta * tb * tb},
      {"main ghost B.S2/A.S2", ra * rb * tb * tb},
      {"main ghost B.S2/B.S1", rb * rb},
  };
  const auto ranking = rtt::analysis::ghost_ranking(g, 0, 0);
  REQUIRE(ranking.entries.size() == 6);
  REQUIRE(close(ranking.base_rms_radius, std::sqrt(294.0 / 127.0), 1e-10));
  for (const auto& e : ranking.entries) {
    const std::string name = g.system.path(e.path).name;
    INFO(name);
    REQUIRE(expected.count(name) == 1);
    REQUIRE(close(e.relative_irradiance, expected.at(name), 1e-10));
    REQUIRE_FALSE(e.focus_offset.has_value());  // collimated
  }
  const std::vector<std::string> order{"main ghost B.S2/B.S1", "main ghost B.S1/A.S2",
                                       "main ghost B.S1/A.S1", "main ghost B.S2/A.S2",
                                       "main ghost B.S2/A.S1", "main ghost A.S2/A.S1"};
  for (std::size_t k = 0; k < order.size(); ++k) {
    REQUIRE(g.system.path(ranking.entries[k].path).name == order[k]);
  }
}

TEST_CASE("M4 ghost of a singlet: paraxial focus and blur from y-nu", "[m4]") {
  // tests/reference/m4/ghost_singlet.rtt.json: plano-convex n = 1.5, R1 = 64 mm at z = 10,
  // plane S2 at z = 13, image at z = 140; collimated marginal ray y = 2 (EPD 4). The ghost
  // reflects at S2 back and at S1 forward. Independent y-nu (Greivenkamp, p. 9-2; signed index,
  // reflection n' = -n, transfer with global dz): the focus after S2 where y = 0 and the blur
  // |y| at the image. Relative 1e-10.
  const rtt::material::MaterialLibrary lib;
  const rtt::coating::CoatingLibrary coatings;
  const GhostSystem g =
      rtt::compile::compile_with_ghosts(load("m4/ghost_singlet.rtt.json"), "main", lib, coatings);
  REQUIRE(g.ghosts.size() == 1);
  const double c1 = 1.0 / 64.0;
  double y = 2.0;
  double n = 1.0;
  double u = 0.0;
  const auto bend = [&](double c, double n2) {
    u = (n * u - y * c * (n2 - n)) / n2;
    n = n2;
  };
  const auto move = [&](double dz) { y += u * dz; };
  move(10.0);
  bend(c1, 1.5);  // refraction into the glass at S1
  move(3.0);
  bend(0.0, -1.5);  // reflection at the plane S2
  move(-3.0);
  bend(c1, 1.5);  // reflection at S1 from inside (n: -1.5 -> +1.5)
  move(3.0);
  bend(0.0, 1.0);  // refraction out of the glass at S2
  const double focus = 13.0 - y / u;
  move(140.0 - 13.0);
  const auto ranking = rtt::analysis::ghost_ranking(g, 0, 0);
  REQUIRE(ranking.entries.size() == 1);
  const auto& e = ranking.entries[0];
  REQUIRE(e.focus_offset.has_value());
  REQUIRE(close(*e.focus_offset, focus - 140.0, 1e-10));
  REQUIRE(e.paraxial_blur_radius.has_value());
  REQUIRE(close(*e.paraxial_blur_radius, std::abs(y), 1e-10));
}

TEST_CASE("M4 feature tour: the first order at the grating plate follows Mansuripur (7b)", "[m4]") {
  // tests/reference/m0/feature_tour.rtt.json, path "first order": refract with order 1 into the
  // tilted N-BK7 grating plate G (300 lines/mm, orientation 90 deg, so the grating vector is the
  // local y axis; ADR 0025, point 4, refract with order since #135). The recorded ray states
  // before and after G.S1 (slots 5 and 6) in the local frame of G.S1 (plane, N = z) must obey
  // (7b): n2 t'_par = n1 t_par + m lambda0 G y_hat, with n1, n2 the compiled indices at the
  // ray's wavelength, and t'_z > 0 from |t'| = 1. A few operations: 1e-12.
  // The tour shows every element of the format and is not meant to compile as it is: A.S1 has a
  // Zernike sag term (compile supports it from M8 on) and the coating "AR_VIS" without a
  // catalog. The test removes the term (A.S1 stays an even asphere) and makes A.S1 a bare
  // Fresnel surface; the grating plate G and the path are those of the file.
  rtt::model::System tour = load("m0/feature_tour.rtt.json");
  auto& a_s1 = std::get<rtt::model::Element>(tour.root.children[2].value).surfaces[0];
  a_s1.shape.terms.clear();
  a_s1.interaction = rtt::model::Fresnel{};
  rtt::material::MaterialLibrary lib;
  lib.add_catalog(std::string(RTT_CATALOG_DIR) + "/schott.agf");
  const CompiledSystem cs = rtt::compile::compile(tour, lib);
  const PathId first = path(cs, "first order");
  const auto& events = cs.path(first).events;
  REQUIRE(cs.surfaces()[events[5].surface].id.str() == "G.S1");
  REQUIRE(events[5].kind == rtt::model::EventKind::Refract);
  REQUIRE(events[5].order == 1);
  RayBatch rays(1);
  rays.pos_y()[0] = 2.0;  // inside the annular stop (1 mm < r < 6 mm)
  rtt::trace::RayPaths recorded;
  [[maybe_unused]] const auto stats =
      rtt::trace::SequentialTracer().trace(cs, first, rays, recorded);
  REQUIRE(rays.status()[0] == RayStatus::Alive);
  const auto& g1 = cs.surfaces()[events[5].surface];
  const Vec3 t = g1.to_local.apply_vector(recorded.direction_at(0, 5));
  const Vec3 t2 = g1.to_local.apply_vector(recorded.direction_at(0, 6));
  const std::uint16_t wl = rays.wl()[0];
  const double n1 = cs.media()[events[5].medium_before].index[wl].real();
  const double n2 = cs.media()[events[5].medium_after].index[wl].real();
  REQUIRE(n2 > 1.4);  // into the N-BK7 plate
  const double lambda_mm = cs.wavelengths_um()[wl] / 1000.0;
  REQUIRE(std::abs(n2 * t2.x() - n1 * t.x()) <= 1e-12);
  REQUIRE(std::abs(n2 * t2.y() - (n1 * t.y() + 1 * lambda_mm * 300.0)) <= 1e-12);
  REQUIRE(t2.z() > 0.0);
  REQUIRE(std::abs(t2.norm() - 1.0) <= 1e-12);
}
