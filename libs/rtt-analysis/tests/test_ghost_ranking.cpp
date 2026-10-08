#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#include "rtt/analysis/ghosts.hpp"
#include "rtt/coating/catalog.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/compile/ghosts.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
#include "rtt/trace/run_control.hpp"
#include "rtt/trace/sources.hpp"

// Ghost ranking (#124, ADR 0027 addendum): rank value rho = (P_g / P_b) (r_b^2 + r0^2) /
// (r_g^2 + r0^2), paraxial focus and blur of each ghost as diagnostics.

using rtt::compile::CompiledSystem;
using rtt::compile::GhostSystem;
using rtt::material::MaterialLibrary;
using rtt::model::Element;
using rtt::model::ElementKind;
using rtt::model::Param;
using rtt::model::Pose;
using rtt::model::Surface;
using rtt::model::SurfaceId;
using rtt::model::System;

namespace {

Surface plane(const std::string& id, double z) {
  Surface s;
  s.id = SurfaceId(id);
  s.pose = Pose::along_z(z);
  return s;
}

Surface sphere(const std::string& id, double z, double radius) {
  Surface s = plane(id, z);
  s.shape.base = rtt::model::Conic{Param(radius), Param(0.0)};
  return s;
}

/// VACUUM, 0.5876 um, object at infinity on axis, EPD 4 mm, stop (r = 5) at z = 0, detector
/// at z = `image_z`, automatic path "main".
System base_system(double image_z) {
  System s;
  s.name = "ghost ranking";
  s.environment.medium = "VACUUM";
  s.wavelengths = {{0.5876, 1.0, true}};
  s.aperture = {rtt::model::SystemApertureType::EntrancePupilDiameter, Param(4.0)};
  s.fields = {rtt::model::FieldType::AngleDeg, {{0.0, 0.0, 1.0}}};
  s.root.name = "root";
  Surface stop = plane("STO", 0.0);
  stop.aperture = rtt::model::CircularAperture{5.0, 0.0};
  s.root.children = {{Element{"stop", ElementKind::Stop, Pose::along_z(0.0), std::nullopt, {stop}}},
                     {Element{"image",
                              ElementKind::Detector,
                              Pose::along_z(image_z),
                              std::nullopt,
                              {plane("IMG", 0.0)}}}};
  s.paths = {{"main", true, {}}};
  return s;
}

/// Two plane plates: A (CONST:1.5) from z = 10 to 15, B (CONST:2.0) from z = 20 to 25.
System two_plates() {
  System s = base_system(40.0);
  s.root.children.insert(s.root.children.begin() + 1,
                         {{Element{"A",
                                   ElementKind::Plate,
                                   Pose::along_z(10.0),
                                   "CONST:1.5",
                                   {plane("A.S1", 0.0), plane("A.S2", 5.0)}}},
                          {Element{"B",
                                   ElementKind::Plate,
                                   Pose::along_z(20.0),
                                   "CONST:2.0",
                                   {plane("B.S1", 0.0), plane("B.S2", 5.0)}}}});
  return s;
}

/// Plano-convex lens CONST:1.5, R1 = 64 mm at z = 10, plane S2 at z = 13, image at z = 140.
System singlet() {
  System s = base_system(140.0);
  s.root.children.insert(s.root.children.begin() + 1,
                         {Element{"L",
                                  ElementKind::Lens,
                                  Pose::along_z(10.0),
                                  "CONST:1.5",
                                  {sphere("L.S1", 0.0, 64.0), plane("L.S2", 3.0)}}});
  return s;
}

/// Sets the aperture of the detector surface (the last child of the root).
void detector_aperture(System& s, const rtt::model::CircularAperture& aperture) {
  std::get<Element>(s.root.children.back().value).surfaces[0].aperture = aperture;
}

GhostSystem with_ghosts(const System& s) {
  const MaterialLibrary lib;
  const rtt::coating::CoatingLibrary coatings;
  return rtt::compile::compile_with_ghosts(s, "main", lib, coatings);
}

std::string ghost_name(const CompiledSystem& cs, rtt::compile::PathId id) {
  return cs.path(id).name;
}

bool close(double a, double b, double rel) {
  return std::abs(a - b) <= rel * std::max(std::abs(a), std::abs(b));
}

}  // namespace

TEST_CASE("ghost ranking of two plates: rho is the Fresnel product, in its order (#124)",
          "[ghosts][ranking]") {
  // Collimated bundle at normal incidence: every ghost leaves collimated on the rays of the
  // useful image (normal reflection and refraction keep x and y), so r_g = r_b and
  // rho = P_g / P_b. A ghost (j, i) passes the surfaces between i and j three times instead of
  // once: P_g / P_b = R_i R_j prod_{i<k<j} T_k^2, with R = ((n - 1) / (n + 1))^2 and T = 1 - R
  // at normal incidence (Byrnes, arXiv:1603.02720v5, Eq. (6), T = 1 - R by Eqs. (21)-(23);
  // docs/quellen.md). R_A = 0.04 (n = 1.5), R_B = 1/9 (n = 2.0). Relative 1e-10.
  const GhostSystem g = with_ghosts(two_plates());
  REQUIRE(g.ghosts.size() == 6);
  const double ra = (0.5 / 2.5) * (0.5 / 2.5);
  const double rb = (1.0 / 3.0) * (1.0 / 3.0);
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
  REQUIRE(ranking.base_power > 0.0);
  // r_b: the collimated hexapolar bundle (6 rings, 127 rays, ring k of 6k rays on radius
  // 2k/6 mm at EPD 4) keeps its positions through plates at normal incidence, with equal
  // weights and centroid 0: r_b^2 = sum_k 6k (2k/6)^2 / 127 = 294/127.
  REQUIRE(close(ranking.base_rms_radius, std::sqrt(294.0 / 127.0), 1e-10));
  for (const auto& e : ranking.entries) {
    const std::string name = ghost_name(g.system, e.path);
    INFO(name);
    REQUIRE(expected.count(name) == 1);
    REQUIRE(close(e.relative_power, expected.at(name), 1e-10));
    REQUIRE(close(e.relative_irradiance, expected.at(name), 1e-10));
    REQUIRE(close(e.rms_radius, ranking.base_rms_radius, 1e-10));
    REQUIRE_FALSE(e.focus_offset.has_value());  // collimated
  }
  // Order: B2/B1 (1/81), B1/A2, B1/A1, B2/A2, B2/A1, A2/A1 (R_A^2).
  const std::vector<std::string> order{"main ghost B.S2/B.S1", "main ghost B.S1/A.S2",
                                       "main ghost B.S1/A.S1", "main ghost B.S2/A.S2",
                                       "main ghost B.S2/A.S1", "main ghost A.S2/A.S1"};
  for (std::size_t k = 0; k < order.size(); ++k) {
    REQUIRE(ghost_name(g.system, ranking.entries[k].path) == order[k]);
  }
}

TEST_CASE("ghost ranking: paraxial focus and blur of a singlet ghost (#124)", "[ghosts][ranking]") {
  // singlet(): collimated marginal ray y = 2 (EPD 4, the stop is the entrance pupil). The ghost
  // reflects at S2 back and at S1 forward. Independent y-nu (Greivenkamp, OPTI-201/202, Sec. 9,
  // p. 9-2, docs/quellen.md, as docs/architecture.md: signed index, n'u' = n u - y phi,
  // phi = c (n' - n), reflection n' = -n, transfer with global dz): focus after S2 where y = 0,
  // blur = |y| at the image.
  const GhostSystem g = with_ghosts(singlet());
  REQUIRE(g.ghosts.size() == 1);
  const double c1 = 1.0 / 64.0;
  // Each step: refraction or reflection at curvature c from index n to n2, then transfer dz.
  double y = 2.0;
  double n = 1.0;
  double u = 0.0;
  const auto bend = [&](double c, double n2) {
    u = (n * u - y * c * (n2 - n)) / n2;
    n = n2;
  };
  const auto move = [&](double dz) { y += u * dz; };
  move(10.0);                         // stop -> S1
  bend(c1, 1.5);                      // refraction into the glass at S1
  move(3.0);                          // S1 -> S2
  bend(0.0, -1.5);                    // reflection at the plane S2
  move(-3.0);                         // S2 -> S1
  bend(c1, 1.5);                      // reflection at S1 from inside (n: -1.5 -> +1.5)
  move(3.0);                          // S1 -> S2
  bend(0.0, 1.0);                     // refraction out of the glass at S2
  const double focus = 13.0 - y / u;  // where the ray meets the axis after S2
  move(140.0 - 13.0);                 // S2 -> image
  const auto ranking = rtt::analysis::ghost_ranking(g, 0, 0);
  REQUIRE(ranking.entries.size() == 1);
  const auto& e = ranking.entries[0];
  REQUIRE(e.focus_offset.has_value());
  REQUIRE(close(*e.focus_offset, focus - 140.0, 1e-10));
  REQUIRE(e.paraxial_blur_radius.has_value());
  REQUIRE(close(*e.paraxial_blur_radius, std::abs(y), 1e-10));
}

TEST_CASE("ghost ranking: rho follows its definition, also for another r0 (#124)",
          "[ghosts][ranking]") {
  // Plates: r_g = r_b, rho does not depend on r0. Singlet: the ghost spot is far wider than
  // the useful image (paraxial blur about 12 mm), so rho grows with r0 as the formula says.
  const GhostSystem plates = with_ghosts(two_plates());
  const GhostSystem lens = with_ghosts(singlet());
  for (const GhostSystem* g : {&plates, &lens}) {
    std::vector<double> rho;
    for (const double r0 : {0.005, 0.05}) {
      INFO("r0 = " << r0);
      rtt::analysis::GhostRankingOptions options;
      options.resolution_radius = r0;
      const auto ranking = rtt::analysis::ghost_ranking(*g, 0, 0, options);
      REQUIRE(ranking.resolution_radius == r0);
      REQUIRE(ranking.entries.size() == g->ghosts.size());
      const double rb2 = ranking.base_rms_radius * ranking.base_rms_radius;
      for (const auto& e : ranking.entries) {
        REQUIRE(close(e.relative_power, e.power / ranking.base_power, 1e-12));
        const double expected =
            e.relative_power * (rb2 + r0 * r0) / (e.rms_radius * e.rms_radius + r0 * r0);
        REQUIRE(close(e.relative_irradiance, expected, 1e-12));
      }
      for (std::size_t k = 1; k < ranking.entries.size(); ++k) {
        REQUIRE(ranking.entries[k - 1].relative_irradiance >=
                ranking.entries[k].relative_irradiance);
      }
      if (g == &lens) {
        REQUIRE(ranking.entries[0].rms_radius > 10.0 * ranking.base_rms_radius);
        rho.push_back(ranking.entries[0].relative_irradiance);
      }
    }
    if (g == &lens) {
      REQUIRE(rho.size() == 2);
      REQUIRE(rho[1] > 2.0 * rho[0]);  // a larger r0 dims the sharp useful image more
    }
  }
}

TEST_CASE("ghost ranking of the reference plate: rho = R^2 (#124)", "[ghosts][ranking]") {
  // tests/reference/m3/fresnel_bk7.rtt.json (plate CONST:1.5168, VACUUM, object at infinity,
  // field 0 on axis), path "main": one ghost, on the rays of the useful image at normal
  // incidence, so rho = P_g / P_b = R^2 with R = ((n - 1) / (n + 1))^2 (Byrnes,
  // arXiv:1603.02720v5, Eq. (6); docs/quellen.md). Relative 1e-10.
  const MaterialLibrary lib;
  const rtt::coating::CoatingLibrary coatings;
  const GhostSystem g = rtt::compile::compile_with_ghosts(
      rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m3/fresnel_bk7.rtt.json"), "main",
      lib, coatings);
  const auto ranking = rtt::analysis::ghost_ranking(g, 0, 0);
  REQUIRE(ranking.entries.size() == 1);
  const double n = 1.5168;
  const double r = ((n - 1.0) / (n + 1.0)) * ((n - 1.0) / (n + 1.0));
  REQUIRE(close(ranking.entries[0].relative_irradiance, r * r, 1e-10));
  REQUIRE_FALSE(ranking.entries[0].focus_offset.has_value());
}

TEST_CASE("ghost ranking: lost ghost rays do not warn (#124)", "[ghosts][ranking]") {
  // Detector of radius 5 mm: the useful image (near focus) arrives completely, the ghost spot
  // (paraxial blur about 12 mm) loses most rays. Ghosts do not warn about lost rays (ADR 0027,
  // addendum #124); their losses are in the entry.
  System s = singlet();
  detector_aperture(s, rtt::model::CircularAperture{5.0, 0.0});
  const GhostSystem g = with_ghosts(s);
  const auto ranking = rtt::analysis::ghost_ranking(g, 0, 0);
  REQUIRE(ranking.warnings.empty());
  REQUIRE(ranking.entries.size() == 1);
  const auto& e = ranking.entries[0];
  REQUIRE(e.losses.launched == 127);
  REQUIRE(e.rays_arrived > 0);
  REQUIRE(2 * e.rays_arrived < e.losses.launched);  // more than the default warning fraction
}

TEST_CASE("ghost ranking: errors and run control (#124)", "[ghosts][ranking]") {
  const GhostSystem g = with_ghosts(two_plates());
  rtt::analysis::GhostRankingOptions bad;
  bad.resolution_radius = -1.0;
  REQUIRE_THROWS_AS(rtt::analysis::ghost_ranking(g, 0, 0, bad), std::invalid_argument);
  REQUIRE_THROWS_AS(rtt::analysis::ghost_ranking(g, 5, 0), std::invalid_argument);
  REQUIRE_THROWS_AS(rtt::analysis::ghost_ranking(g, 0, 3), std::invalid_argument);
  rtt::analysis::GhostRankingOptions zero;
  zero.resolution_radius = 0.0;
  REQUIRE_THROWS_AS(rtt::analysis::ghost_ranking(g, 0, 0, zero), std::invalid_argument);
  rtt::analysis::GhostRankingOptions fraction;
  fraction.lost_warning_fraction = 1.5;
  REQUIRE_THROWS_AS(rtt::analysis::ghost_ranking(g, 0, 0, fraction), std::invalid_argument);
  rtt::analysis::GhostRankingOptions no_rays;
  no_rays.sampling = rtt::trace::RandomPupil{0, 0};
  REQUIRE_THROWS_AS(rtt::analysis::ghost_ranking(g, 0, 0, no_rays), std::invalid_argument);
  // No ghosts: stop and detector only.
  REQUIRE_THROWS_AS(rtt::analysis::ghost_ranking(with_ghosts(base_system(40.0)), 0, 0),
                    std::invalid_argument);
  // No useful image: an annular detector outside the focused bundle.
  System annulus = singlet();
  detector_aperture(annulus, rtt::model::CircularAperture{50.0, 5.0});
  REQUIRE_THROWS_AS(rtt::analysis::ghost_ranking(with_ghosts(annulus), 0, 0),
                    rtt::analysis::AnalysisError);

  rtt::trace::RunControl control;
  control.cancel = rtt::trace::CancelToken();
  control.progress = [](const rtt::trace::Progress&) {};
  const auto plain = rtt::analysis::ghost_ranking(g, 0, 0);
  const auto controlled = rtt::analysis::ghost_ranking(g, 0, 0, {}, control);
  REQUIRE(plain.entries.size() == controlled.entries.size());
  for (std::size_t k = 0; k < plain.entries.size(); ++k) {
    REQUIRE(plain.entries[k].path == controlled.entries[k].path);
    REQUIRE(plain.entries[k].relative_irradiance == controlled.entries[k].relative_irradiance);
  }
  control.cancel->request_cancel();
  REQUIRE_THROWS_AS(rtt::analysis::ghost_ranking(g, 0, 0, {}, control), rtt::trace::Cancelled);
}
