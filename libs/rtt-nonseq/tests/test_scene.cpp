// Scene of the non-sequential tracer (ADR 0057, points 2 and 3; PR 2): bodies with media on
// their sides, generated mantles and shoulders, thin faces, the error codes nonseq.*, and the
// nearest hit through the BVH at analytic positions and against a brute-force search.
//
// Tolerances: the poses are shifts along z (exact), so a hit parameter differs from its analytic
// value only by the rounding of a few operations (a quadratic root, one sqrt, the transforms):
// a few ulp, i.e. 1e-15 relative for plane and shoulder hits and 1e-14 for walls and spheres,
// whose roots pass through sqrt and a division. The Newton case states its own bound.

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <variant>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
#include "rtt/nonseq/scene.hpp"
#include "rtt/trace/apply_event.hpp"
#include "test_support.hpp"

using Catch::Approx;
using rtt::compile::CompiledSystem;
using rtt::geom::HitStatus;
using rtt::material::MaterialLibrary;
using rtt::math::Vec3;
using rtt::model::CircularAperture;
using rtt::model::Element;
using rtt::model::ElementKind;
using rtt::model::Node;
using rtt::model::Pose;
using rtt::model::RectangularAperture;
using rtt::model::Surface;
using rtt::model::SurfaceId;
using rtt::model::System;
using rtt::nonseq::FaceKind;
using rtt::nonseq::HitRegion;
using rtt::nonseq::kNoFace;
using rtt::nonseq::Scene;
using rtt::nonseq::SceneError;
using rtt::nonseq::SceneHit;
using rtt::nonseq::SceneOptions;
using rtt::nonseq::WallShape;
using rtt::trace::RayState;

namespace {

constexpr double kR = 51.68;  // radius of L1.S1 in make_singlet, mm
constexpr double kRim = 12.7;

/// The reference singlet in CONST:1.5168 (make_singlet): stop r = 10 at z = 0, L1 with a convex
/// sphere at z = 5 and a plane at z = 9 (both r = 12.7), detector without aperture.
System singlet() {
  System s = rtt::model::test::make_singlet();
  rtt::model::test::element(s, 1).material = "CONST:1.5168";
  return s;
}

CompiledSystem compiled(const System& s) {
  const MaterialLibrary lib;
  return rtt::compile::compile(s, lib);
}

RayState ray(const Vec3& pos, const Vec3& dir) {
  RayState r;
  r.pos = pos;
  r.dir = dir.normalized();
  return r;
}

/// Sag of a sphere of radius R at height r, z = R - sqrt(R^2 - r^2).
double sphere_sag(double radius, double r) {
  return radius - std::sqrt(radius * radius - r * r);
}

/// Index of the first face of `kind` of element `element`.
std::uint32_t face_of(const Scene& scene,
                      FaceKind kind,
                      std::uint32_t element,
                      std::uint32_t nth = 0) {
  for (std::uint32_t f = 0; f < scene.faces().size(); ++f) {
    if (scene.faces()[f].kind == kind && scene.faces()[f].element == element && nth-- == 0) {
      return f;
    }
  }
  return kNoFace;
}

/// The scene error with exactly one diagnostic, or an empty diagnostic if there is none.
rtt::model::Diagnostic scene_error(const CompiledSystem& system, const SceneOptions& o = {}) {
  try {
    (void)Scene::build(system, o);
  } catch (const SceneError& e) {
    REQUIRE(e.diagnostics().size() == 1);
    return e.diagnostics().front();
  }
  return {};
}

}  // namespace

TEST_CASE("scene of the singlet: faces, media on their sides, mantle", "[nonseq][scene]") {
  const CompiledSystem system = compiled(singlet());
  const Scene scene = Scene::build(system);
  const std::uint32_t glass = system.elements()[1].media[0];
  const std::uint32_t env = system.environment_medium();
  REQUIRE(glass != env);
  // Stop, L1.S1, L1.S2, the mantle of segment 0, detector.
  REQUIRE(scene.faces().size() == 5);
  const auto& f = scene.faces();
  CHECK(f[0].kind == FaceKind::Stop);
  CHECK(f[1].kind == FaceKind::Body);
  CHECK(f[2].kind == FaceKind::Body);
  CHECK(f[3].kind == FaceKind::Mantle);
  CHECK(f[4].kind == FaceKind::Detector);
  // Segment 0 lies on the +z side of S1 and on the -z side of S2 (ADR 0057, point 2).
  CHECK(f[1].medium_positive == glass);
  CHECK(f[1].medium_negative == env);
  CHECK(f[2].medium_positive == env);
  CHECK(f[2].medium_negative == glass);
  CHECK(f[3].generated.wall == WallShape::Circle);
  CHECK(f[3].generated.a == kRim);
  CHECK_FALSE(f[3].generated.constant_z);
  CHECK(f[3].medium_positive == env);
  CHECK(f[3].medium_negative == glass);
  CHECK(std::holds_alternative<rtt::model::Fresnel>(scene.surfaces()[f[3].surface].interaction));
  // The detector has no aperture: it stays outside the BVH (#450, decision 3 of PR 2).
  REQUIRE(scene.unbounded().size() == 1);
  CHECK(scene.unbounded()[0] == 4);
  CHECK_FALSE(scene.bvh().empty());
}

TEST_CASE("nearest hits of the singlet at analytic positions", "[nonseq][scene]") {
  const CompiledSystem system = compiled(singlet());
  const Scene scene = Scene::build(system);
  const std::uint32_t glass = system.elements()[1].media[0];
  const std::uint32_t env = system.environment_medium();

  SECTION("axial ray through the stop opening onto the vertex of S1") {
    const SceneHit h = scene.closest_hit(ray(Vec3(0.0, 0.0, -1.0), Vec3::UnitZ()));
    REQUIRE(h.status == HitStatus::Hit);
    CHECK(h.face == 1);
    CHECK(h.region == HitRegion::Clear);
    CHECK(h.local.t == Approx(6.0).epsilon(1e-15));
    CHECK_FALSE(h.from_positive);  // d . n > 0: from the -z side, the environment
    CHECK(h.medium_before == env);
    CHECK(h.medium_beyond == glass);
  }
  SECTION("from S1 inside the glass: the self-hit is excluded, S2 is next") {
    const SceneHit h = scene.closest_hit(ray(Vec3(0.0, 0.0, 5.0), Vec3::UnitZ()), 1);
    REQUIRE(h.status == HitStatus::Hit);
    CHECK(h.face == 2);
    CHECK(h.local.t == Approx(4.0).epsilon(1e-15));
    CHECK(h.medium_before == glass);
    CHECK(h.medium_beyond == env);
  }
  SECTION("inside the glass towards the rim: the mantle at r = 12.7") {
    const SceneHit h = scene.closest_hit(ray(Vec3(0.0, 0.0, 7.5), Vec3::UnitX()));
    REQUIRE(h.status == HitStatus::Hit);
    CHECK(h.face == 3);
    CHECK(h.local.t == Approx(kRim).epsilon(1e-14));
    CHECK(h.point.isApprox(Vec3(kRim, 0.0, 7.5), 1e-14));
    CHECK(h.normal.isApprox(Vec3::UnitX(), 1e-15));
    CHECK(h.medium_before == glass);
    CHECK(h.medium_beyond == env);
  }
  SECTION("from outside in front of the S1 rim: past the mantle onto S1") {
    // At z = 6 (local z = 1 of S1) the mantle is not there (the glass starts at the rim height
    // 1.585 mm); the ray meets S1 where its sag is 1, r^2 = 2 R - 1.
    const SceneHit h = scene.closest_hit(ray(Vec3(20.0, 0.0, 6.0), -Vec3::UnitX()));
    REQUIRE(h.status == HitStatus::Hit);
    CHECK(h.face == 1);
    CHECK(h.local.t == Approx(20.0 - std::sqrt(2.0 * kR - 1.0)).epsilon(1e-13));
    CHECK(h.medium_before == env);
    CHECK(h.medium_beyond == glass);
  }
  SECTION("on the stop screen: absorbed") {
    const SceneHit h = scene.closest_hit(ray(Vec3(0.0, 13.0, -1.0), Vec3::UnitZ()));
    REQUIRE(h.status == HitStatus::Hit);
    CHECK(h.face == 0);
    CHECK(h.region == HitRegion::Absorbed);
    CHECK(h.local.t == Approx(1.0).epsilon(1e-15));
  }
  SECTION("beyond the screen (2 x 10 mm) and the lens: onto the detector") {
    const SceneHit h = scene.closest_hit(ray(Vec3(0.0, 25.0, -1.0), Vec3::UnitZ()));
    REQUIRE(h.status == HitStatus::Hit);
    CHECK(h.face == 4);
    CHECK(h.local.t == Approx(107.363).epsilon(1e-14));
  }
  SECTION("away from everything: Missed") {
    const SceneHit h = scene.closest_hit(ray(Vec3(0.0, 0.0, -1.0), -Vec3::UnitZ()));
    CHECK(h.status == HitStatus::Missed);
    CHECK(h.face == kNoFace);
  }
}

TEST_CASE("nearest hit through the BVH equals the brute-force search", "[nonseq][scene]") {
  const CompiledSystem system = compiled(singlet());
  const Scene scene = Scene::build(system);
  std::size_t hits = 0;
  for (int i = 0; i < 12; ++i) {
    for (int k = 0; k < 12; ++k) {
      // Rays from points on and around the lens into a fan of directions.
      const Vec3 o(-15.0 + 2.5 * i, 9.0 - 1.5 * k, -2.0 + 1.0 * ((i + k) % 11));
      const double phi = 0.52 * i + 0.31 * k;
      const double theta = 0.27 * (i + 2 * k);
      const Vec3 d(std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi),
                   std::cos(theta));
      const RayState r = ray(o, d);
      const SceneHit h = scene.closest_hit(r);
      SceneHit best;
      double best_t = std::numeric_limits<double>::infinity();
      for (std::uint32_t f = 0; f < scene.faces().size(); ++f) {
        const SceneHit c = scene.intersect_face(r, f, 0.0, best_t);
        if (c.status == HitStatus::Hit) {
          best = c;
          best_t = c.local.t;
        }
      }
      REQUIRE(h.status == best.status);
      CHECK(h.face == best.face);
      if (h.status == HitStatus::Hit) {
        ++hits;
        CHECK(h.local.t == best.local.t);
      }
    }
  }
  CHECK(hits > 40);  // guard: many rays do hit something
}

TEST_CASE("step between unequal circles: wall at the larger radius and a shoulder",
          "[nonseq][scene]") {
  System s = singlet();
  rtt::model::test::element(s, 1).surfaces[1].aperture = CircularAperture{10.0, 0.0};
  const CompiledSystem system = compiled(s);
  const Scene scene = Scene::build(system);
  const std::uint32_t glass = system.elements()[1].media[0];
  const std::uint32_t env = system.environment_medium();
  const std::uint32_t wall = face_of(scene, FaceKind::Mantle, 1);
  const std::uint32_t shoulder = face_of(scene, FaceKind::Shoulder, 1);
  REQUIRE(wall != kNoFace);
  REQUIRE(shoulder != kNoFace);
  const auto& w = scene.faces()[wall].generated;
  CHECK(w.a == kRim);
  CHECK(w.constant_z);
  CHECK(w.z_lo == Approx(sphere_sag(kR, kRim)).epsilon(1e-14));  // rim of S1, frame of S1
  CHECK(w.z_hi == Approx(4.0).epsilon(1e-15));                   // rim of S2
  const auto& sh = scene.faces()[shoulder].generated;
  CHECK(sh.r_in == 10.0);
  CHECK(sh.r_out == kRim);
  CHECK(sh.z_lo == Approx(4.0).epsilon(1e-15));
  CHECK(sh.orientation == 1.0);  // the glass lies below the ring (towards S1)

  SECTION("from behind onto the shoulder at z = 9") {
    const SceneHit h = scene.closest_hit(ray(Vec3(11.0, 0.0, 20.0), -Vec3::UnitZ()));
    REQUIRE(h.status == HitStatus::Hit);
    CHECK(h.face == shoulder);
    CHECK(h.local.t == Approx(11.0).epsilon(1e-15));
    CHECK(h.from_positive);
    CHECK(h.medium_before == env);
    CHECK(h.medium_beyond == glass);
  }
  SECTION("inside the glass onto the wall at r = 12.7") {
    const SceneHit h = scene.closest_hit(ray(Vec3(0.0, 0.0, 7.5), Vec3::UnitX()));
    REQUIRE(h.status == HitStatus::Hit);
    CHECK(h.face == wall);
    CHECK(h.local.t == Approx(kRim).epsilon(1e-14));
    CHECK(h.medium_before == glass);
  }
  SECTION("inside the glass under the shoulder towards +z: the shoulder from below") {
    const SceneHit h = scene.closest_hit(ray(Vec3(11.0, 0.0, 8.0), Vec3::UnitZ()));
    REQUIRE(h.status == HitStatus::Hit);
    CHECK(h.face == shoulder);
    CHECK(h.local.t == Approx(1.0).epsilon(1e-14));
    CHECK_FALSE(h.from_positive);
    CHECK(h.medium_before == glass);
  }
}

TEST_CASE("rectangular plate: four edge walls with outward normals", "[nonseq][scene]") {
  System s = singlet();
  Element& l1 = rtt::model::test::element(s, 1);
  l1.kind = ElementKind::Plate;
  l1.surfaces[0].shape.base = rtt::model::Plane{};
  l1.surfaces[0].aperture = RectangularAperture{4.0, 3.0};
  l1.surfaces[1].aperture = RectangularAperture{4.0, 3.0};
  const Scene scene = Scene::build(compiled(s));
  std::size_t edges = 0;
  for (const auto& f : scene.faces()) {
    if (f.kind == FaceKind::Mantle && f.generated.wall == WallShape::Edge) ++edges;
  }
  CHECK(edges == 4);
  const SceneHit hx = scene.closest_hit(ray(Vec3(0.0, 0.0, 7.0), Vec3::UnitX()));
  REQUIRE(hx.status == HitStatus::Hit);
  CHECK(scene.faces()[hx.face].kind == FaceKind::Mantle);
  CHECK(hx.local.t == Approx(4.0).epsilon(1e-15));
  CHECK(hx.normal.isApprox(Vec3::UnitX(), 1e-15));
  const SceneHit hy = scene.closest_hit(ray(Vec3(0.0, 0.0, 7.0), -Vec3::UnitY()));
  REQUIRE(hy.status == HitStatus::Hit);
  CHECK(hy.local.t == Approx(3.0).epsilon(1e-15));
  CHECK(hy.normal.isApprox(-Vec3::UnitY(), 1e-15));
  // Outside the plate in front of it: the wall is not there, the ray passes.
  const SceneHit miss = scene.closest_hit(ray(Vec3(-10.0, 0.0, 4.5), Vec3::UnitX()));
  CHECK(miss.status == HitStatus::Missed);
}

namespace {

/// A triangular prism of CONST:1.5 (plate of one material with three plane surfaces) in the
/// y-z plane: base z = 0 for |y| <= 5, faces through (5, 0)-(0, 5) and (-5, 0)-(0, 5), x in
/// [-5, 5]; with an explicit path through stop and detector only (compile rejects a prism on
/// an automatic path). `offset` moves the second face along z (non-convex for offset < 0).
System prism(double offset = 0.0) {
  System s = singlet();
  Surface base;
  base.id = SurfaceId("P.S1");
  base.aperture = RectangularAperture{5.0, 5.0};
  Surface right;
  right.id = SurfaceId("P.S2");
  right.pose.position = {0.0, 2.5, 2.5 + offset};
  right.pose.rotation_deg = {-45.0, 0.0, 0.0};
  right.aperture = RectangularAperture{5.0, 2.5 * std::sqrt(2.0)};
  Surface left;
  left.id = SurfaceId("P.S3");
  left.pose.position = {0.0, -2.5, 2.5};
  left.pose.rotation_deg = {45.0, 0.0, 0.0};
  left.aperture = RectangularAperture{5.0, 2.5 * std::sqrt(2.0)};
  s.root.children[1] =
      Node{Element{"P", ElementKind::Plate, Pose::along_z(5.0), "CONST:1.5", {base, right, left}}};
  s.paths = {{"p",
              false,
              {{SurfaceId("STO"), rtt::model::EventKind::Transmit},
               {SurfaceId("IMG"), rtt::model::EventKind::Transmit}}}};
  return s;
}

}  // namespace

TEST_CASE("prism: the inside is the side of the corner centroid", "[nonseq][scene]") {
  const CompiledSystem system = compiled(prism());
  const Scene scene = Scene::build(system);
  const std::uint32_t glass = system.elements()[1].media[0];
  const std::uint32_t env = system.environment_medium();
  CHECK(face_of(scene, FaceKind::Mantle, 1) == kNoFace);  // closed by its faces
  const SceneHit in = scene.closest_hit(ray(Vec3(0.0, 1.0, 4.0), Vec3::UnitZ()));
  REQUIRE(in.status == HitStatus::Hit);
  CHECK(in.face == 1);  // the base at z = 5
  CHECK(in.local.t == Approx(1.0).epsilon(1e-15));
  CHECK(in.medium_before == env);
  CHECK(in.medium_beyond == glass);
  // Inside, upwards at y = 1: the face y + z = 5 (local) at z = 4, t = 3.
  const SceneHit out = scene.closest_hit(ray(Vec3(0.0, 1.0, 6.0), Vec3::UnitZ()));
  REQUIRE(out.status == HitStatus::Hit);
  CHECK(out.face == 2);
  CHECK(out.local.t == Approx(3.0).epsilon(1e-14));
  CHECK(out.medium_before == glass);
  CHECK(out.medium_beyond == env);
}

TEST_CASE("thin faces: mirror substrate, transparent detector, absorbing rims, inactive",
          "[nonseq][scene]") {
  System s = singlet();
  Surface m;
  m.id = SurfaceId("M");
  m.aperture = CircularAperture{5.0, 0.0};
  s.root.children.push_back(
      Node{Element{"mirror", ElementKind::Mirror, Pose::along_z(200.0), "CONST:1.5", {m}}});
  const CompiledSystem system = compiled(s);
  SceneOptions o;
  o.absorbing_rims = {1};
  o.transparent_detectors = {2};
  const Scene scene = Scene::build(system, o);
  const std::uint32_t mirror = face_of(scene, FaceKind::Mirror, 3);
  REQUIRE(mirror != kNoFace);
  CHECK(scene.faces()[mirror].substrate.has_value());
  CHECK(scene.faces()[face_of(scene, FaceKind::Detector, 2)].transparent);
  const std::uint32_t mantle = face_of(scene, FaceKind::Mantle, 1);
  CHECK(std::holds_alternative<rtt::model::Absorber>(
      scene.surfaces()[scene.faces()[mantle].surface].interaction));

  System off = singlet();
  rtt::model::test::element(off, 1).active = false;
  const Scene without = Scene::build(compiled(off));
  CHECK(without.faces().size() == 2);  // stop and detector
  CHECK(face_of(without, FaceKind::Body, 1) == kNoFace);

  SceneOptions bad;
  bad.screen_scale = 1.0;
  CHECK_THROWS_AS(Scene::build(system, bad), std::invalid_argument);
  SceneOptions wrong;
  wrong.absorbing_rims = {0};  // the stop
  CHECK_THROWS_AS(Scene::build(system, wrong), std::invalid_argument);
}

TEST_CASE("scene errors nonseq.* with their locations", "[nonseq][scene]") {
  SECTION("nonseq.aperture_missing") {
    System s = singlet();
    rtt::model::test::element(s, 1).surfaces[1].aperture.reset();
    const auto d = scene_error(compiled(s));
    CHECK(d.code == "nonseq.aperture_missing");
    CHECK(d.location == "/root/children/1/surfaces/1");
  }
  SECTION("nonseq.segment_thickness") {
    System s = singlet();
    rtt::model::test::element(s, 1).surfaces[1].pose = Pose::along_z(0.0);
    const auto d = scene_error(compiled(s));
    CHECK(d.code == "nonseq.segment_thickness");
    CHECK(d.location == "/root/children/1");
  }
  SECTION("nonseq.rim_unsupported: different aperture kinds") {
    System s = singlet();
    rtt::model::test::element(s, 1).surfaces[1].aperture = RectangularAperture{12.7, 12.7};
    const auto d = scene_error(compiled(s));
    CHECK(d.code == "nonseq.rim_unsupported");
    CHECK(d.location == "/root/children/1");
  }
  SECTION("nonseq.rim_unsupported: tilted axis") {
    System s = singlet();
    rtt::model::test::element(s, 1).surfaces[1].pose.rotation_deg = {1.0, 0.0, 0.0};
    CHECK(scene_error(compiled(s)).code == "nonseq.rim_unsupported");
  }
  SECTION("nonseq.rim_unsupported: inner radius on one surface") {
    System s = singlet();
    rtt::model::test::element(s, 1).surfaces[1].aperture = CircularAperture{12.7, 2.0};
    CHECK(scene_error(compiled(s)).code == "nonseq.rim_unsupported");
  }
  SECTION("nonseq.rim_unsupported: step at a surface that is not rotationally symmetric") {
    System s = singlet();
    Surface& s2 = rtt::model::test::element(s, 1).surfaces[1];
    s2.shape.base = rtt::model::Biconic{-200.0, -300.0, 0.0, 0.0};
    s2.aperture = CircularAperture{10.0, 0.0};
    CHECK(scene_error(compiled(s)).code == "nonseq.rim_unsupported");
  }
  SECTION("nonseq.plate_not_convex") {
    const auto d = scene_error(compiled(prism(-1.0)));
    CHECK(d.code == "nonseq.plate_not_convex");
    CHECK(d.location == "/root/children/1/surfaces/0");
  }
}

namespace {

/// The singlet with S1 an even asphere of the same radius and A4 r^4.
System asphere_singlet(double a4) {
  System s = singlet();
  rtt::model::test::element(s, 1).surfaces[0].shape.base =
      rtt::model::EvenAsphere{kR, 0.0, {rtt::model::Param(a4)}};
  return s;
}

/// Sag of the asphere of asphere_singlet() at height r (ADR 0043 form with k = 0).
double asphere_sag(double a4, double r) {
  return sphere_sag(kR, r) + a4 * r * r * r * r;
}

}  // namespace

TEST_CASE("Newton surface: the next crossing after one outside the aperture is bracketed",
          "[nonseq][scene][newton]") {
  // ADR 0057, Nachtrag PR 2, i (#450): before, intersect_surface(t_min = t1) restarted Newton at
  // the base conic, fell back onto the same crossing and gave NoConvergence for both rays below.
  SECTION("A4 > 0: inside the glass the mantle at r = 12.7 is found, not NoConvergence") {
    const Scene scene = Scene::build(compiled(asphere_singlet(1e-7)));
    const SceneHit h = scene.closest_hit(ray(Vec3(0.0, 0.0, 7.5), Vec3::UnitX()));
    REQUIRE(h.status == HitStatus::Hit);
    CHECK(scene.faces()[h.face].kind == FaceKind::Mantle);
    CHECK(h.local.t == Approx(kRim).epsilon(1e-14));
  }
  SECTION("A4 < 0: beside the lens in the corner of its box the ray reaches the detector") {
    const Scene scene = Scene::build(compiled(asphere_singlet(-1e-7)));
    const SceneHit h = scene.closest_hit(ray(Vec3(12.0, 12.0, 1.0), Vec3::UnitZ()));
    REQUIRE(h.status == HitStatus::Hit);
    CHECK(scene.faces()[h.face].kind == FaceKind::Detector);
    CHECK(h.local.t == Approx(105.363).epsilon(1e-14));
  }
  SECTION("a second crossing inside the aperture is found at the analytic point") {
    // A thin asphere (no mantle) at z = 5, aperture 12.7; the chord through the surface points
    // at r = -14 (outside the aperture) and r = 8 (inside) crosses the surface exactly there
    // (the sag is convex, so a chord has no third crossing).
    constexpr double a4 = 1e-7;
    System s = singlet();
    Surface t;
    t.id = SurfaceId("T");
    t.shape.base = rtt::model::EvenAsphere{kR, 0.0, {rtt::model::Param(a4)}};
    t.aperture = CircularAperture{kRim, 0.0};
    s.root.children[1] = Node{Element{"T", ElementKind::ThinElement, Pose::along_z(5.0), {}, {t}}};
    const Scene scene = Scene::build(compiled(s));
    const Vec3 p1(-14.0, 0.0, 5.0 + asphere_sag(a4, 14.0));
    const Vec3 p2(8.0, 0.0, 5.0 + asphere_sag(a4, 8.0));
    const Vec3 dir = (p2 - p1).normalized();
    const SceneHit h = scene.closest_hit(ray(p1 - 2.0 * dir, dir));
    REQUIRE(h.status == HitStatus::Hit);
    CHECK(scene.faces()[h.face].kind == FaceKind::Thin);
    // Bisection stops at |F| <= 1e-12 mm; with |dF/dt| ~ 0.2 along this chord the t error stays
    // below 1e-11 mm, well inside 1e-9 mm.
    CHECK(std::abs(h.local.t - (2.0 + (p2 - p1).norm())) < 1e-9);
    CHECK((h.point - p2).norm() < 1e-9);
  }
}

TEST_CASE("Newton surface: a bracketed second crossing with its normal",
          "[nonseq][scene][newton]") {
  // As the chord section above, with A4 < 0: the root of the base conic then lies behind the
  // crossing at r = 8, so restarting Newton there (the path before Nachtrag PR 2, i) does not
  // find it; the bracketing does. Normal from the gradient, (-dz/dx, 0, 1) / |...| at r = 8 with
  // dz/dr = r / sqrt(R^2 - r^2) + 4 A4 r^3.
  constexpr double a4 = -1e-7;
  System s = singlet();
  Surface t;
  t.id = SurfaceId("T");
  t.shape.base = rtt::model::EvenAsphere{kR, 0.0, {rtt::model::Param(a4)}};
  t.aperture = CircularAperture{kRim, 0.0};
  s.root.children[1] = Node{Element{"T", ElementKind::ThinElement, Pose::along_z(5.0), {}, {t}}};
  const Scene scene = Scene::build(compiled(s));
  const Vec3 p1(-14.0, 0.0, 5.0 + asphere_sag(a4, 14.0));
  const Vec3 p2(8.0, 0.0, 5.0 + asphere_sag(a4, 8.0));
  const Vec3 dir = (p2 - p1).normalized();
  const SceneHit h = scene.closest_hit(ray(p1 - 2.0 * dir, dir));
  REQUIRE(h.status == HitStatus::Hit);
  CHECK(scene.faces()[h.face].kind == FaceKind::Thin);
  CHECK((h.point - p2).norm() < 1e-9);  // bound as in the section above
  const double slope = 8.0 / std::sqrt(kR * kR - 64.0) + 4.0 * a4 * 512.0;
  CHECK(h.normal.isApprox(Vec3(-slope, 0.0, 1.0).normalized(), 1e-12));
  CHECK(h.from_positive);  // the chord runs down onto the surface from its +z side
}

TEST_CASE("cemented doublet: glass on both sides of the cemented surface", "[nonseq][scene]") {
  System s = singlet();
  Element& l1 = rtt::model::test::element(s, 1);
  Surface s3 = l1.surfaces[1];
  s3.id = SurfaceId("L1.S3");
  s3.pose = Pose::along_z(7.0);
  l1.surfaces.push_back(s3);
  l1.material.reset();
  l1.segment_materials = {"CONST:1.5", "CONST:1.7"};
  const CompiledSystem system = compiled(s);
  const Scene scene = Scene::build(system);
  const std::uint32_t g1 = system.elements()[1].media[0];
  const std::uint32_t g2 = system.elements()[1].media[1];
  REQUIRE(g1 != g2);
  const auto& cemented = scene.faces()[2];
  CHECK(cemented.kind == FaceKind::Body);
  CHECK(cemented.medium_negative == g1);
  CHECK(cemented.medium_positive == g2);
  // Two mantles, one per segment, each with its own glass.
  CHECK(scene.faces()[face_of(scene, FaceKind::Mantle, 1, 0)].medium_negative == g1);
  CHECK(scene.faces()[face_of(scene, FaceKind::Mantle, 1, 1)].medium_negative == g2);

  // Segments on the same side of the inner surface overlap: S2 behind S3.
  l1.surfaces[1].pose = Pose::along_z(8.0);
  const auto d = scene_error(compiled(s));
  CHECK(d.code == "nonseq.segment_thickness");
  CHECK(d.location == "/root/children/1/surfaces/1");
}

TEST_CASE("bore: an inner radius on both surfaces gives an inner mantle", "[nonseq][scene]") {
  SECTION("equal inner radii: a wall at r = 2 with its normal towards the axis") {
    System s = singlet();
    for (auto& surf : rtt::model::test::element(s, 1).surfaces) {
      surf.aperture = CircularAperture{kRim, 2.0};
    }
    const CompiledSystem system = compiled(s);
    const Scene scene = Scene::build(system);
    const std::uint32_t bore = face_of(scene, FaceKind::Mantle, 1, 1);
    REQUIRE(bore != kNoFace);
    CHECK(scene.faces()[bore].generated.a == 2.0);
    CHECK(scene.faces()[bore].generated.orientation == -1.0);
    const SceneHit h = scene.closest_hit(ray(Vec3(5.0, 0.0, 7.5), -Vec3::UnitX()));
    REQUIRE(h.status == HitStatus::Hit);
    CHECK(h.face == bore);
    CHECK(h.local.t == Approx(3.0).epsilon(1e-15));
    CHECK(h.normal.isApprox(-Vec3::UnitX(), 1e-15));
    CHECK(h.medium_before == system.elements()[1].media[0]);
    CHECK(h.medium_beyond == system.environment_medium());
  }
  SECTION("unequal inner radii: the wall at the smaller one and a shoulder at the other rim") {
    System s = singlet();
    auto& surfs = rtt::model::test::element(s, 1).surfaces;
    surfs[0].aperture = CircularAperture{kRim, 2.0};
    surfs[1].aperture = CircularAperture{kRim, 3.0};
    const CompiledSystem system = compiled(s);
    const Scene scene = Scene::build(system);
    const std::uint32_t shoulder = face_of(scene, FaceKind::Shoulder, 1);
    REQUIRE(shoulder != kNoFace);
    const auto& g = scene.faces()[shoulder].generated;
    CHECK(g.r_in == 2.0);
    CHECK(g.r_out == 3.0);
    CHECK(g.z_lo == Approx(4.0).epsilon(1e-15));  // inner rim of S2, frame of S1
    CHECK(g.orientation == 1.0);                  // the glass lies towards S1
    const std::uint32_t wall = face_of(scene, FaceKind::Mantle, 1, 1);
    CHECK(scene.faces()[wall].generated.a == 2.0);
    CHECK(scene.faces()[wall].generated.constant_z);
    // From behind through the hole of S2 onto the shoulder at z = 9.
    const SceneHit h = scene.closest_hit(ray(Vec3(2.5, 0.0, 20.0), -Vec3::UnitZ()));
    REQUIRE(h.status == HitStatus::Hit);
    CHECK(h.face == shoulder);
    CHECK(h.local.t == Approx(11.0).epsilon(1e-15));
    CHECK(h.from_positive);
    CHECK(h.medium_beyond == system.elements()[1].media[0]);
  }
}

TEST_CASE("elliptical plate: the wall (x / a)^2 + (y / b)^2 = 1 and its normal",
          "[nonseq][scene]") {
  System s = singlet();
  Element& l1 = rtt::model::test::element(s, 1);
  l1.kind = ElementKind::Plate;
  l1.surfaces[0].shape.base = rtt::model::Plane{};
  l1.surfaces[0].aperture = rtt::model::EllipticalAperture{4.0, 3.0};
  l1.surfaces[1].aperture = rtt::model::EllipticalAperture{4.0, 3.0};
  const Scene scene = Scene::build(compiled(s));
  const std::uint32_t wall = face_of(scene, FaceKind::Mantle, 1);
  REQUIRE(wall != kNoFace);
  CHECK(scene.faces()[wall].generated.wall == WallShape::Ellipse);
  // Along the diagonal: (t / sqrt(2))^2 (1/16 + 1/9) = 1, t = 12 sqrt(2) / 5, x = y = 12/5.
  const SceneHit h = scene.closest_hit(ray(Vec3(0.0, 0.0, 7.0), Vec3(1.0, 1.0, 0.0)));
  REQUIRE(h.status == HitStatus::Hit);
  CHECK(h.face == wall);
  CHECK(h.local.t == Approx(12.0 * std::sqrt(2.0) / 5.0).epsilon(1e-14));
  // Normal along (x / a^2, y / b^2) = (2.4 / 16, 2.4 / 9).
  CHECK(h.normal.isApprox(Vec3(2.4 / 16.0, 2.4 / 9.0, 0.0).normalized(), 1e-14));
}

TEST_CASE("stop: inner radius and obscurations absorb on the screen, beyond it nothing",
          "[nonseq][scene]") {
  System s = singlet();
  Surface& stop = rtt::model::test::element(s, 0).surfaces[0];
  stop.aperture = CircularAperture{10.0, 2.0};
  stop.obscurations = {rtt::model::CircularObscuration{1.0, 5.0, 0.0},
                       rtt::model::CircularObscuration{1.0, 30.0, 0.0}};
  const Scene scene = Scene::build(compiled(s));
  for (const double x : {0.0, 5.0}) {
    const SceneHit h = scene.closest_hit(ray(Vec3(x, 0.0, -1.0), Vec3::UnitZ()));
    REQUIRE(h.status == HitStatus::Hit);
    CHECK(h.face == 0);
    CHECK(h.region == HitRegion::Absorbed);
  }
  // The obscuration at x = 30 lies beyond the screen (2 x 10 mm): the stop is not there
  // (ADR 0057, Nachtrag PR 2, e), the ray reaches the detector.
  const SceneHit beyond = scene.closest_hit(ray(Vec3(30.0, 0.0, -1.0), Vec3::UnitZ()));
  REQUIRE(beyond.status == HitStatus::Hit);
  CHECK(scene.faces()[beyond.face].kind == FaceKind::Detector);
  // Open between inner radius and obscuration: through to S1.
  const SceneHit open = scene.closest_hit(ray(Vec3(3.0, 0.0, -1.0), Vec3::UnitZ()));
  REQUIRE(open.status == HitStatus::Hit);
  CHECK(open.face == 1);
}

TEST_CASE("equal rims outside the domain of a surface: nonseq.rim_unsupported", "[nonseq][scene]") {
  // S1 a sphere of radius R = -10 (centre on the -z side) with a rim of 12.7 > |R|: the rim lies
  // beyond the hemisphere, the body would be open there.
  System s = singlet();
  rtt::model::test::element(s, 1).surfaces[0].shape.base =
      rtt::model::Conic{rtt::model::Param(-10.0), rtt::model::Param(0.0)};
  CHECK(scene_error(compiled(s)).code == "nonseq.rim_unsupported");
}

TEST_CASE("NoConvergence counts only if the hidden crossing could be the nearest",
          "[nonseq][scene][newton]") {
  // A thin asphere with a strong A4 = 1e-4 (sag 4.2 mm at the rim instead of 1.6 mm for the
  // base sphere): Newton from the base conic fails for the rays below (status NoConvergence on
  // the face itself). ADR 0057, Nachtrag PR 2, g: the ray ends with NoConvergence only if the
  // face's box starts before the best hit.
  System s = singlet();
  Surface t;
  t.id = SurfaceId("T");
  t.shape.base = rtt::model::EvenAsphere{kR, 0.0, {rtt::model::Param(1e-4)}};
  t.aperture = CircularAperture{kRim, 0.0};
  s.root.children[1] = Node{Element{"T", ElementKind::ThinElement, Pose::along_z(5.0), {}, {t}}};
  const Scene scene = Scene::build(compiled(s));
  const std::uint32_t thin = face_of(scene, FaceKind::Thin, 1);
  const auto from = [](double x0, double theta) {
    return ray(Vec3(x0, 0.0, -1.0), Vec3(std::sin(theta), 0.0, std::cos(theta)));
  };
  constexpr double kInfT = std::numeric_limits<double>::infinity();
  SECTION("nothing nearer: the ray ends with NoConvergence") {
    const RayState r = from(2.0, 0.84);
    REQUIRE(scene.intersect_face(r, thin, 0.0, kInfT).status == HitStatus::NoConvergence);
    const SceneHit h = scene.closest_hit(r);
    CHECK(h.status == HitStatus::NoConvergence);
    CHECK(h.face == thin);
  }
  SECTION("the stop screen is hit before the face's box: the hit counts") {
    const RayState r = from(10.0, 0.96);
    REQUIRE(scene.intersect_face(r, thin, 0.0, kInfT).status == HitStatus::NoConvergence);
    const SceneHit h = scene.closest_hit(r);
    REQUIRE(h.status == HitStatus::Hit);
    CHECK(h.face == 0);
    CHECK(h.region == HitRegion::Absorbed);
    CHECK(h.local.t == Approx(1.0 / std::cos(0.96)).epsilon(1e-14));
  }
  SECTION("the ray never enters the face's box: it reaches the detector") {
    const RayState r = from(8.0, 0.72);
    REQUIRE(scene.intersect_face(r, thin, 0.0, kInfT).status == HitStatus::NoConvergence);
    const SceneHit h = scene.closest_hit(r);
    REQUIRE(h.status == HitStatus::Hit);
    CHECK(scene.faces()[h.face].kind == FaceKind::Detector);
    CHECK(h.local.t == Approx(107.363 / std::cos(0.72)).epsilon(1e-14));
  }
}

TEST_CASE("known limit: the first Newton crossing at grazing incidence",
          "[nonseq][scene][newton]") {
  // ADR 0057, point 3, "Bekannte Grenze": the Newton intersection starts at the base conic and
  // may find the farther of two crossings inside the aperture. A chord of a thin asphere that
  // grazes it between r = 5 and r = 5.5 crosses it at both points; stage 1 accepts either.
  constexpr double a4 = 1e-6;
  System s = singlet();
  Surface t;
  t.id = SurfaceId("T");
  t.shape.base = rtt::model::EvenAsphere{kR, 0.0, {rtt::model::Param(a4)}};
  t.aperture = CircularAperture{kRim, 0.0};
  s.root.children[1] = Node{Element{"T", ElementKind::ThinElement, Pose::along_z(5.0), {}, {t}}};
  const Scene scene = Scene::build(compiled(s));
  const Vec3 p1(5.0, 0.0, 5.0 + asphere_sag(a4, 5.0));
  const Vec3 p2(5.5, 0.0, 5.0 + asphere_sag(a4, 5.5));
  const Vec3 dir = (p2 - p1).normalized();
  const SceneHit h = scene.closest_hit(ray(p1 - 2.0 * dir, dir));
  REQUIRE(h.status == HitStatus::Hit);
  CHECK(scene.faces()[h.face].kind == FaceKind::Thin);
  CHECK(std::min((h.point - p1).norm(), (h.point - p2).norm()) < 1e-9);
}
