// Real rays with the stop on a surface (#53, ADR 0039, point 5 and the acceptance reading of the
// maintainer on #53): the chief rays agree with the workaround of #34 within a derived bound; all
// other rays are checked against the definition, each aimed ray meets the stop surface L2.S1
// locally at sigma p R_s. The workaround is a different stop for them (a plane stop element
// 1e-8 mm behind the vertex), so they differ physically by about u * sag.

#include <algorithm>
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <optional>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "rtt/coating/catalog.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/compile/perturbation.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/math/types.hpp"
#include "rtt/model/model.hpp"
#include "rtt/trace/ray_paths.hpp"
#include "rtt/trace/sequential.hpp"
#include "rtt/trace/sources.hpp"

using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::material::MaterialLibrary;
using rtt::math::Vec3;
using rtt::model::SurfaceId;
using rtt::model::System;
using rtt::trace::kAimTolerance;
using rtt::trace::PupilPoint;

namespace {

System load(const std::string& relative) {
  return rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/" + relative);
}

const MaterialLibrary& schott_m2() {
  struct Lib {  // MaterialLibrary is not movable
    MaterialLibrary lib;
    Lib() { lib.add_catalog(std::string(RTT_CATALOG_DIR) + "/m2/schott.agf"); }
  };
  static const Lib l;
  return l.lib;
}

const CompiledSystem& surface_stop() {
  static const CompiledSystem cs =
      rtt::compile::compile(load("m2/cooke_triplet_surface_stop.rtt.json"), schott_m2());
  return cs;
}

const CompiledSystem& workaround() {
  static const CompiledSystem cs =
      rtt::compile::compile(load("m2/cooke_triplet.rtt.json"), schott_m2());
  return cs;
}

/// One traced ray: the local hit on L2.S1 (first event there) and the hit on the image surface.
struct Traced {
  bool alive = false;
  Vec3 l2_local;
  Vec3 image;
};

/// Real-aimed rays of `field` at the reference wavelength for `points`, traced with recording.
std::vector<Traced> trace(const CompiledSystem& cs,
                          std::uint16_t field,
                          const std::vector<PupilPoint>& points) {
  const PathId path{0};
  const std::array<std::uint16_t, 1> fields{field};
  rtt::trace::RayBatch rays = rtt::trace::make_rays(cs, path, fields, cs.reference_wavelength(),
                                                    rtt::trace::PointListPupil{points});
  rtt::trace::RayPaths paths;
  [[maybe_unused]] const auto stats = rtt::trace::SequentialTracer().trace(cs, path, rays, paths);
  const std::uint32_t l2 = cs.find_surface(SurfaceId("L2.S1")).value();
  std::size_t l2_slot = 0;
  for (std::size_t k = 0; k < paths.event_surfaces.size(); ++k) {
    if (paths.event_surfaces[k] == l2) {
      l2_slot = k + 1;  // slot k + 1 is the state after event k
      break;
    }
  }
  REQUIRE(l2_slot > 0);
  std::vector<Traced> out;
  for (std::size_t r = 0; r < paths.ray_count(); ++r) {
    Traced t;
    t.alive = paths.lost_at[r] == -1 && paths.count[r] == paths.slots;
    t.l2_local = cs.surfaces()[l2].to_local.apply_point(paths.position_at(r, l2_slot));
    t.image = paths.position_at(r, paths.slots - 1);
    out.push_back(t);
  }
  return out;
}

}  // namespace

TEST_CASE("surface stop: the chief rays agree with the workaround", "[trace][aperture-stop]") {
  // Bound per field (derived on #53, ADR 0039, "Schranke für die Hauptstrahlen"): the workaround
  // aims the chief ray at the centre of STO, 1e-8 mm behind the vertex of L2.S1 in the glass;
  // the surface stop aims it at the centre of L2.S1. On L2.S1 the two targets are eps_f apart,
  // measured at run time as the local |y| of the workaround's chief ray on L2.S1; each aiming
  // adds up to kAimTolerance. The image point depends on the hit on L2.S1 through k_f, measured
  // as the difference quotient dy_image / dy_L2 of the workaround's rays at p = (0, +-1e-3)
  // (both measured, linear to O(eps^2) over 1e-8 mm). B_f = 2 (eps_f + 2 kAimTolerance) |k_f| +
  // 1e-12 mm; the factor 2 covers the measurement of k_f, 1e-12 mm the rounding of the image
  // coordinate (|y| <= 20 mm, eps_M * 20 ~ 4e-15 mm). Expected below 2e-10 mm.
  const std::vector<PupilPoint> chief{{0.0, 0.0}};
  const std::vector<PupilPoint> pair{{0.0, 1e-3}, {0.0, -1e-3}};
  for (const std::uint16_t field : {std::uint16_t{0}, std::uint16_t{1}, std::uint16_t{2}}) {
    INFO("field " << field);
    const Traced a = trace(surface_stop(), field, chief).at(0);
    const Traced b = trace(workaround(), field, chief).at(0);
    const std::vector<Traced> d = trace(workaround(), field, pair);
    REQUIRE(a.alive);
    REQUIRE(b.alive);
    REQUIRE(d.at(0).alive);
    REQUIRE(d.at(1).alive);
    const double eps = std::hypot(b.l2_local.x(), b.l2_local.y());
    const double k = (d[0].image.y() - d[1].image.y()) / (d[0].l2_local.y() - d[1].l2_local.y());
    const double bound = 2.0 * (eps + 2.0 * kAimTolerance) * std::abs(k) + 1e-12;
    INFO(std::setprecision(17) << "eps_f " << eps << ", k_f " << k << ", B_f " << bound
                               << ", image y " << a.image.y() << " / " << b.image.y());
    CHECK(std::abs(a.image.x() - b.image.x()) <= bound);
    CHECK(std::abs(a.image.y() - b.image.y()) <= bound);
    CHECK(bound < 2e-10);
    // Content guard: the surface stop's chief ray meets L2.S1 at its centre, the workaround's
    // off it on the oblique fields (|u_bar| delta, about 2e-9 and 3e-9 mm at 14 and 20 deg).
    CHECK(std::hypot(a.l2_local.x(), a.l2_local.y()) <= kAimTolerance);
    if (field > 0) CHECK(eps > kAimTolerance);
  }
}

TEST_CASE("surface stop: every aimed ray meets L2.S1 at sigma p R_s", "[trace][aperture-stop]") {
  // The aim target of ADR 0007 on the stop surface (ADR 0039, point 5): local (x, y) of the hit
  // = sigma (px, py) R_s within kAimTolerance (plus rounding, 1e-12 mm). sigma R_s is taken from
  // the sign of the ray at p = (0, 1). Content guard: the workaround, a plane stop, meets L2.S1
  // at the same pupil point measurably elsewhere (about u * sag, > 1e-3 mm at |p| = 1).
  const CompiledSystem& cs = surface_stop();
  const double r_s = rtt::trace::stop_radius(cs, PathId{0}, cs.reference_wavelength());
  std::vector<PupilPoint> points{{0.0, 1.0}};
  for (const double p : {-1.0, -0.5, 0.5, 1.0}) {
    points.push_back({0.0, p});
    points.push_back({p, 0.0});
  }
  points.push_back({0.6, -0.6});
  for (const std::uint16_t field : {std::uint16_t{0}, std::uint16_t{1}, std::uint16_t{2}}) {
    INFO("field " << field);
    const std::vector<Traced> t = trace(cs, field, points);
    REQUIRE(t.at(0).alive);
    const double sigma = t[0].l2_local.y() > 0.0 ? 1.0 : -1.0;
    double max_dev = 0.0;
    for (std::size_t i = 0; i < points.size(); ++i) {
      INFO("p = (" << points[i].px << ", " << points[i].py << ")");
      REQUIRE(t[i].alive);
      const double dx = t[i].l2_local.x() - sigma * points[i].px * r_s;
      const double dy = t[i].l2_local.y() - sigma * points[i].py * r_s;
      CHECK(std::hypot(dx, dy) <= kAimTolerance + 1e-12);
    }
    const std::vector<Traced> w = trace(workaround(), field, points);
    for (std::size_t i = 0; i < points.size(); ++i) {
      if (std::hypot(points[i].px, points[i].py) < 1.0) continue;
      REQUIRE(w[i].alive);
      const double dx = w[i].l2_local.x() - sigma * points[i].px * r_s;
      const double dy = w[i].l2_local.y() - sigma * points[i].py * r_s;
      max_dev = std::max(max_dev, std::hypot(dx, dy));
    }
    INFO(std::setprecision(6) << "largest deviation of the workaround at |p| = 1: " << max_dev);
    CHECK(max_dev > 1e-3);
  }
}

TEST_CASE("surface stop: a mirror as the stop", "[trace][aperture-stop]") {
  // m2/paraboloid_stop without its stop element, the stop on the mirror M: the chief ray of the
  // axial field meets the mirror at its centre (local 0, 0) within kAimTolerance.
  System s = load("m2/paraboloid_stop.rtt.json");
  std::erase_if(s.root.children, [](const rtt::model::Node& node) {
    const auto* e = std::get_if<rtt::model::Element>(&node.value);
    return e != nullptr && e->kind == rtt::model::ElementKind::Stop;
  });
  s.aperture.stop = SurfaceId("M");
  const CompiledSystem cs = rtt::compile::compile(s, MaterialLibrary{});
  const rtt::trace::AimedRay chief =
      rtt::trace::aim_ray(cs, PathId{0}, 0, cs.reference_wavelength(), 0.0, 0.0);
  REQUIRE(chief.ray.status == rtt::trace::RayStatus::Alive);
  const rtt::trace::AimedRay rim =
      rtt::trace::aim_ray(cs, PathId{0}, 0, cs.reference_wavelength(), 0.0, 1.0);
  REQUIRE(rim.ray.status == rtt::trace::RayStatus::Alive);
  // The rim ray of the 60 mm entrance pupil starts 30 mm off the axis (parallel to it).
  CHECK(std::abs(std::abs(rim.ray.pos.y()) - 30.0) <= 1e-9);
  CHECK(std::abs(chief.ray.pos.y()) <= 1e-9);
}

TEST_CASE("surface stop: the stop follows the perturbed lens",
          "[trace][aperture-stop][tolerance]") {
  // ADR 0039, point 7: a trial aims with the nominal aiming data (ADR 0036, point 2), the target
  // lies in the local coordinates of the stop surface of the perturbed system. Decentring L2 by
  // 0.05 mm moves L2.S1 and with it the stop: the chief ray still meets L2.S1 at its centre, which
  // now lies 0.05 mm off the axis. Decentring L1 leaves the stop where it is: the chief ray meets
  // L2.S1 at its centre on the axis. Content guard: the global hits differ by the decenter.
  const System s = load("m2/cooke_triplet_surface_stop.rtt.json");
  const CompiledSystem& nominal = surface_stop();
  const std::uint16_t wl = nominal.reference_wavelength();
  const rtt::trace::AimReference reference = rtt::trace::aim_reference(nominal, PathId{0}, wl);
  const rtt::coating::CoatingLibrary coatings;
  const std::uint32_t l2 = nominal.find_surface(SurfaceId("L2.S1")).value();
  auto chief_on_l2 = [&](const std::string& node, std::uint16_t field) {
    rtt::compile::Perturbation p;
    p.nodes[node].translation_mm = Vec3(0.0, 0.05, 0.0);
    const CompiledSystem cs = rtt::compile::compile(s, schott_m2(), coatings, 0, p);
    const std::array<std::uint16_t, 1> fields{field};
    rtt::trace::RayBatch batch =
        rtt::trace::make_rays(cs, PathId{0}, fields, wl, rtt::trace::SinglePupilPoint{},
                              rtt::trace::Aiming::Real, reference);
    REQUIRE(batch.status()[0] == rtt::trace::RayStatus::Alive);
    rtt::trace::RayPaths paths;
    [[maybe_unused]] const auto stats =
        rtt::trace::SequentialTracer().trace(cs, PathId{0}, batch, paths);
    REQUIRE(paths.lost_at[0] == -1);
    const Vec3 global = paths.position_at(0, 3);  // after event 2, L2.S1
    return std::pair<Vec3, Vec3>{global, cs.surfaces()[l2].to_local.apply_point(global)};
  };
  for (const std::uint16_t field : {std::uint16_t{0}, std::uint16_t{2}}) {
    INFO("field " << field);
    const auto [g2, l2_local] = chief_on_l2("L2", field);
    CHECK(std::hypot(l2_local.x(), l2_local.y()) <= kAimTolerance);
    CHECK(std::abs(g2.y() - 0.05) <= kAimTolerance + 1e-12);
    const auto [g1, l1_local] = chief_on_l2("L1", field);
    CHECK(std::hypot(l1_local.x(), l1_local.y()) <= kAimTolerance);
    CHECK(std::abs(g1.y()) <= kAimTolerance + 1e-12);
  }
}
