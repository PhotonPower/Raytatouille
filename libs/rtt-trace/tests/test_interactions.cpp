// Interactions in the tracer: P and weight (#61, ADR 0021). Sources: S. J. Byrnes,
// arXiv:1603.02720v5 (Fresnel, power factors, absorption), W.-S. T. Lam (PRT), see
// docs/quellen.md.

#include <Eigen/Geometry>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "rtt/coating/catalog.hpp"
#include "rtt/coating/transfer_matrix.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
#include "rtt/polar/fresnel.hpp"
#include "rtt/polar/ideal.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/sequential.hpp"

using Cx = std::complex<double>;
using rtt::compile::compile;
using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::material::MaterialLibrary;
using rtt::math::CMat3;
using rtt::math::Mat3;
using rtt::math::Vec3;
using rtt::model::Element;
using rtt::model::ElementKind;
using rtt::model::EventKind;
using rtt::model::Param;
using rtt::model::Pose;
using rtt::model::Surface;
using rtt::model::SurfaceId;
using rtt::model::System;
using rtt::trace::RayBatch;
using rtt::trace::RayStatus;
using rtt::trace::SequentialTracer;
using CVec3 = Eigen::Vector3cd;

namespace {

constexpr double kPi = std::numbers::pi;
constexpr double kGlass = 1.5168;  // N-BK7 at d (architecture, acceptance table)

/// System in vacuum (n = 1 exactly) with one wavelength and an automatic path.
System base_system(double wavelength_um = 0.55) {
  System s;
  s.name = "interactions";
  s.environment.medium = "VACUUM";
  s.wavelengths = {{wavelength_um, 1.0, true}};
  s.aperture = {rtt::model::SystemApertureType::EntrancePupilDiameter, Param(10.0)};
  s.fields = {rtt::model::FieldType::AngleDeg, {{0.0, 0.0, 1.0}}};
  s.root.name = "root";
  s.paths = {{"main", true, {}}};
  return s;
}

Surface plane(const std::string& id, double z_mm = 0.0) {
  Surface s;
  s.id = SurfaceId(id);
  s.pose = Pose::along_z(z_mm);
  return s;
}

/// Plane-parallel plate from z to z + thickness.
Element plate(const std::string& name, double z, double thickness, const std::string& material) {
  return Element{name,
                 ElementKind::Plate,
                 Pose::along_z(z),
                 material,
                 {plane(name + ".S1"), plane(name + ".S2", thickness)}};
}

Element thin(const std::string& name,
             double z,
             rtt::model::Interaction interaction,
             double rotation_z_deg = 0.0) {
  Surface s = plane(name);
  s.interaction = std::move(interaction);
  Pose pose = Pose::along_z(z);
  pose.rotation_deg[2] = Param(rotation_z_deg);
  return Element{name, ElementKind::ThinElement, pose, std::nullopt, {s}};
}

Element detector(double z) {
  return Element{"D", ElementKind::Detector, Pose::along_z(z), std::nullopt, {plane("IMG")}};
}

/// Traces one ray from (0, 0, z0) along dir on path `path`.
RayBatch trace_one(const CompiledSystem& cs,
                   const Vec3& start,
                   const Vec3& dir,
                   PathId path = PathId{0}) {
  RayBatch rays(1);
  rays.pos_x()[0] = start.x();
  rays.pos_y()[0] = start.y();
  rays.pos_z()[0] = start.z();
  rays.dir_x()[0] = dir.x();
  rays.dir_y()[0] = dir.y();
  rays.dir_z()[0] = dir.z();
  [[maybe_unused]] const auto stats = SequentialTracer().trace(cs, path, rays);
  return rays;
}

/// |P E|^2 for a unit transverse field E: power fraction for that polarization (ADR 0021).
double polarized_power(const RayBatch& rays, const Vec3& e) {
  return (rays.prt_matrix(0) * e.cast<Cx>()).squaredNorm();
}

}  // namespace

TEST_CASE("N-BK7 plate at normal incidence: weight = (1 - R)^2", "[interactions]") {
  // Issue #61: uncoated plate, no multiple reflections; R = ((n - 1)/(n + 1))^2 (Byrnes, Eq. (6)
  // at normal incidence, #56). Unpolarized and every polarization give (1 - R)^2. Tolerance
  // 1e-12.
  System s = base_system();
  s.root.children = {{plate("P", 0.0, 5.0, "CONST:1.5168")}, {detector(10.0)}};
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(s, lib);
  const RayBatch rays = trace_one(cs, Vec3(0.0, 0.0, -1.0), Vec3(0.0, 0.0, 1.0));
  const double r = (kGlass - 1.0) / (kGlass + 1.0);
  const double expected = (1.0 - r * r) * (1.0 - r * r);
  REQUIRE(rays.status()[0] == RayStatus::Alive);
  REQUIRE(std::abs(rays.weight()[0] - expected) <= 1e-12);
  REQUIRE(std::abs(polarized_power(rays, Vec3(1.0, 0.0, 0.0)) - expected) <= 1e-12);
  REQUIRE(std::abs(polarized_power(rays, Vec3(0.0, 1.0, 0.0)) - expected) <= 1e-12);
}

TEST_CASE("oblique plate: s and p transmittances, Brewster p without loss", "[interactions]") {
  // A plane-parallel plate keeps the plane of incidence, so s and p stay decoupled:
  // power(s) = T_s1 T_s2, power(p) = T_p1 T_p2 with the Fresnel transmittances of #56 at the
  // same tangential invariant; unpolarized weight = (power(s) + power(p)) / 2. At Brewster's
  // angle arctan(n) both p transmittances are 1. Tolerance 1e-12.
  System s = base_system();
  s.root.children = {{plate("P", 0.0, 5.0, "CONST:1.5168")}, {detector(20.0)}};
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(s, lib);
  for (const double theta : {0.7, std::atan(kGlass)}) {
    const Vec3 dir(0.0, std::sin(theta), std::cos(theta));
    const RayBatch rays = trace_one(cs, Vec3(0.0, 0.0, -1.0), dir);
    const double xi = std::sin(theta);
    const auto in = rtt::polar::fresnel_power(Cx(1.0), Cx(kGlass), xi);
    const auto out = rtt::polar::fresnel_power(Cx(kGlass), Cx(1.0), xi);
    const double ps = in.transmittance_s * out.transmittance_s;
    const double pp = in.transmittance_p * out.transmittance_p;
    // s is perpendicular to the plane of incidence (y-z plane): x; p = k x s.
    const Vec3 e_s(1.0, 0.0, 0.0);
    const Vec3 e_p = dir.cross(e_s);
    INFO("theta " << theta);
    REQUIRE(std::abs(polarized_power(rays, e_s) - ps) <= 1e-12);
    REQUIRE(std::abs(polarized_power(rays, e_p) - pp) <= 1e-12);
    REQUIRE(std::abs(rays.weight()[0] - 0.5 * (ps + pp)) <= 1e-12);
    if (theta != 0.7) REQUIRE(std::abs(pp - 1.0) <= 1e-12);
  }
}

TEST_CASE("Malus through the tracer: two ideal polarizers", "[interactions]") {
  // Issue #61: thin elements with IdealPolarizer, axis (1, 0, 0) in element coordinates; the
  // second element is rotated by theta about z, so its axis is rotated by theta globally
  // (compile, ADR 0021). Unpolarized source: 1/2 after the first polarizer, cos^2 theta of that
  // after the second: weight = cos^2(theta) / 2. Crossed polarizers give 0, not 1/4 (the product
  // of per-surface averages would). Tolerance 1e-12.
  const MaterialLibrary lib;
  for (const double deg : {0.0, 30.0, 60.0, 90.0}) {
    System s = base_system();
    s.root.children = {{thin("POL1", 0.0, rtt::model::IdealPolarizer{{1.0, 0.0, 0.0}, 0.0})},
                       {thin("POL2", 5.0, rtt::model::IdealPolarizer{{1.0, 0.0, 0.0}, 0.0}, deg)},
                       {detector(10.0)}};
    const CompiledSystem cs = compile(s, lib);
    const RayBatch rays = trace_one(cs, Vec3(0.0, 0.0, -1.0), Vec3(0.0, 0.0, 1.0));
    const double c = std::cos(deg * kPi / 180.0);
    INFO("theta " << deg);
    REQUIRE(std::abs(rays.weight()[0] - 0.5 * c * c) <= 1e-12);
    REQUIRE(std::abs(polarized_power(rays, Vec3(1.0, 0.0, 0.0)) - c * c) <= 1e-12);
  }
}

TEST_CASE("absorbing plate: Beer-Lambert exp(-4 pi kappa L / lambda)", "[interactions]") {
  // Issue #61: plate with kappa > 0 between ideal anti-reflection surfaces (T = 1), so only the
  // volume absorption remains: weight = exp(-4 pi kappa d / lambda_vac) (Byrnes, Eqs. (1), (2),
  // (17), (18)); also at oblique incidence with the geometric path d / cos(theta_t), theta_t
  // from Snell with Re(n). Tolerance 1e-12 relative.
  const double kappa = 1e-6;
  const double d = 5.0;  // mm
  const double lambda_mm = 0.55e-3;
  System s = base_system(0.55);
  Element p = plate("P", 0.0, d, "CONST:1.5,1e-6");
  for (auto& surf : p.surfaces) surf.interaction = rtt::model::IdealAntiReflection{};
  s.root.children = {{p}, {detector(20.0)}};
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(s, lib);
  for (const double theta : {0.0, 0.5}) {
    const RayBatch rays =
        trace_one(cs, Vec3(0.0, 0.0, -1.0), Vec3(0.0, std::sin(theta), std::cos(theta)));
    const double cos_t = std::sqrt(1.0 - std::pow(std::sin(theta) / 1.5, 2));
    const double expected = std::exp(-4.0 * kPi * kappa * (d / cos_t) / lambda_mm);
    INFO("theta " << theta);
    REQUIRE(std::abs(rays.weight()[0] - expected) <= 1e-12 * expected);
  }
}

TEST_CASE("ideal mirror and a mirror without material: no loss, E_r = -(I - 2 N N^T) E",
          "[interactions]") {
  // ADR 0021: IdealMirror and Fresnel on a Mirror without material reflect with r_s = -1,
  // r_p = +1, the field of an ideal conductor (docs/architecture.md). Tolerance 1e-14.
  const MaterialLibrary lib;
  for (const bool ideal : {true, false}) {
    System s = base_system();
    Surface m = plane("M");
    if (ideal) m.interaction = rtt::model::IdealMirror{};
    Pose pose = Pose::along_z(5.0);
    pose.rotation_deg[0] = Param(30.0);  // tilted about x
    Element mirror{"M", ElementKind::Mirror, pose, std::nullopt, {m}};
    s.root.children = {{mirror}};
    const CompiledSystem cs = compile(s, lib);
    const RayBatch rays = trace_one(cs, Vec3(0.0, 0.0, 0.0), Vec3(0.0, 0.0, 1.0));
    const Vec3 n = cs.surfaces()[0].to_global.apply_vector(Vec3(0.0, 0.0, 1.0));
    const Mat3 conductor = -(Mat3::Identity() - 2.0 * n * n.transpose());
    INFO("ideal " << ideal);
    REQUIRE(std::abs(rays.weight()[0] - 1.0) <= 1e-14);
    for (const Vec3 e : {Vec3(1.0, 0.0, 0.0), Vec3(0.0, 1.0, 0.0)}) {
      const CVec3 out = rays.prt_matrix(0) * e.cast<Cx>();
      REQUIRE((out - (conductor * e).cast<Cx>()).norm() <= 1e-14);
    }
  }
}

TEST_CASE("ideal beam splitter: reflected and transmitted powers add up to 1", "[interactions]") {
  // IdealBeamSplitter (R_s, R_p) with the path choosing Reflect or Transmit (ADR 0021):
  // unpolarized weights (R_s + R_p)/2 and 1 - (R_s + R_p)/2. Tolerance 1e-12.
  System s = base_system();
  Surface bs = plane("BS");
  bs.interaction = rtt::model::IdealBeamSplitter{0.3, 0.6};
  Pose pose = Pose::along_z(5.0);
  pose.rotation_deg[0] = Param(45.0);
  s.root.children = {{Element{"BS", ElementKind::ThinElement, pose, std::nullopt, {bs}}}};
  s.paths = {{"transmit", false, {{SurfaceId("BS"), EventKind::Transmit, 0}}},
             {"reflect", false, {{SurfaceId("BS"), EventKind::Reflect, 0}}}};
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(s, lib);
  const RayBatch t = trace_one(cs, Vec3(0.0, 0.0, 0.0), Vec3(0.0, 0.0, 1.0), PathId{0});
  const RayBatch r = trace_one(cs, Vec3(0.0, 0.0, 0.0), Vec3(0.0, 0.0, 1.0), PathId{1});
  REQUIRE(std::abs(r.weight()[0] - 0.45) <= 1e-12);
  REQUIRE(std::abs(t.weight()[0] - 0.55) <= 1e-12);
}

TEST_CASE("coating: forward stack on the way in, reversed stack from the substrate",
          "[interactions]") {
  // ADR 0019: the substrate is the inside of the element; light from inside sees the reversed
  // stack. Plate of n = 1.52 with V_AR (demo catalogue, two different layers) on its exit
  // surface and IdealAntiReflection on its entrance surface; the ray goes in, is reflected at
  // the coated surface from the substrate side and leaves again. At normal incidence the
  // transverse part of P is then r_s of the reversed stack (rtt-coating, Byrnes, #58), which
  // differs from r_s of the forward stack. Tolerance 1e-12.
  rtt::coating::CoatingLibrary coatings;
  coatings.add_catalog(std::string(RTT_CATALOG_DIR) + "/coatings/demo.json");
  System s = base_system(0.55);
  Element p = plate("P", 0.0, 5.0, "CONST:1.52");
  p.surfaces[0].interaction = rtt::model::IdealAntiReflection{};
  p.surfaces[1].interaction = rtt::model::CoatingRef{"DEMO:V_AR"};
  s.root.children = {{p}, {detector(-2.0)}};
  s.paths = {{"ghost",
              false,
              {{SurfaceId("P.S1"), EventKind::Refract, 0},
               {SurfaceId("P.S2"), EventKind::Reflect, 0},
               {SurfaceId("P.S1"), EventKind::Refract, 0},
               {SurfaceId("IMG"), EventKind::Transmit, 0}}}};
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(s, lib, coatings);
  const RayBatch rays = trace_one(cs, Vec3(0.0, 0.0, -1.0), Vec3(0.0, 0.0, 1.0));
  REQUIRE(rays.status()[0] == RayStatus::Alive);
  const auto& layers = cs.coatings()[0].layers[0];
  const std::vector<rtt::coating::Layer<double>> reversed(layers.rbegin(), layers.rend());
  const auto back = rtt::coating::stack_amplitudes<double>(Cx(1.52), reversed, Cx(1.0), 0.0, 0.55);
  const auto front = rtt::coating::stack_amplitudes<double>(Cx(1.0), layers, Cx(1.52), 0.0, 0.55);
  const Cx p00 = rays.prt_matrix(0)(0, 0);
  INFO("P(0,0) " << p00 << ", r_s reversed " << back.rs << ", forward " << front.rs);
  REQUIRE(std::abs(p00 - back.rs) <= 1e-12);
  REQUIRE(std::abs(p00 - front.rs) > 1e-3);
  REQUIRE(std::abs(rays.weight()[0] - std::norm(back.rs)) <= 1e-12);
}

TEST_CASE("absorber and impossible interaction events", "[interactions]") {
  // ADR 0021: Absorber -> Absorbed with weight 0; IdealMirror cannot refract and an ideal
  // polarizer acts only on Transmit -> EventImpossible.
  const MaterialLibrary lib;
  {
    System s = base_system();
    s.root.children = {{thin("A", 0.0, rtt::model::Absorber{})}};
    const RayBatch rays = trace_one(compile(s, lib), Vec3(0.0, 0.0, -1.0), Vec3(0.0, 0.0, 1.0));
    REQUIRE(rays.status()[0] == RayStatus::Absorbed);
    REQUIRE(rays.weight()[0] == 0.0);
  }
  {
    System s = base_system();
    Element p = plate("P", 0.0, 5.0, "CONST:1.5");
    p.surfaces[0].interaction = rtt::model::IdealMirror{};
    s.root.children = {{p}};
    const RayBatch rays = trace_one(compile(s, lib), Vec3(0.0, 0.0, -1.0), Vec3(0.0, 0.0, 1.0));
    REQUIRE(rays.status()[0] == RayStatus::EventImpossible);
  }
  {
    System s = base_system();
    s.root.children = {{thin("POL", 0.0, rtt::model::IdealPolarizer{{1.0, 0.0, 0.0}, 0.0})}};
    s.paths = {{"bounce", false, {{SurfaceId("POL"), EventKind::Reflect, 0}}}};
    const RayBatch rays = trace_one(compile(s, lib), Vec3(0.0, 0.0, -1.0), Vec3(0.0, 0.0, 1.0));
    REQUIRE(rays.status()[0] == RayStatus::EventImpossible);
  }
}

TEST_CASE("Fresnel reflection without a change of medium vanishes; Transmit is a dummy passage",
          "[interactions]") {
  // ADR 0021 (decided for #61): Reflect with Fresnel at a surface that is not a Mirror and has
  // the same medium on both sides gives r_s = r_p = 0 (Byrnes, Eq. (6) with n_i = n_t), i.e.
  // weight 0 with status Alive; Transmit at a Fresnel or coated surface leaves P and weight.
  const MaterialLibrary lib;
  {
    System s = base_system();
    s.root.children = {{thin("T", 0.0, rtt::model::Fresnel{})}};
    s.paths = {{"bounce", false, {{SurfaceId("T"), EventKind::Reflect, 0}}}};
    const RayBatch rays = trace_one(compile(s, lib), Vec3(0.0, 0.0, -1.0), Vec3(0.0, 0.0, 1.0));
    REQUIRE(rays.status()[0] == RayStatus::Alive);
    REQUIRE(rays.weight()[0] == 0.0);
  }
  {
    System s = base_system();
    s.root.children = {{plate("P", 0.0, 5.0, "CONST:1.5")}};
    s.paths = {{"dummy", false, {{SurfaceId("P.S1"), EventKind::Transmit, 0}}}};
    const RayBatch rays = trace_one(compile(s, lib), Vec3(0.0, 0.0, -1.0), Vec3(0.0, 0.0, 1.0));
    REQUIRE(rays.status()[0] == RayStatus::Alive);
    REQUIRE(rays.weight()[0] == 1.0);
    REQUIRE(rays.prt_matrix(0) == CMat3::Identity());
  }
}

TEST_CASE("medium_beyond: the other side of a surface, also for Reflect", "[interactions]") {
  // ADR 0021: CompiledEvent::medium_beyond is the medium a Refract would enter; for a ghost
  // reflection inside a plate it is the environment, for the entrance reflection the plate.
  System s = base_system();
  s.root.children = {{plate("P", 0.0, 5.0, "CONST:1.5")}};
  s.paths = {{"ghost",
              false,
              {{SurfaceId("P.S1"), EventKind::Reflect, 0},
               {SurfaceId("P.S1"), EventKind::Refract, 0},
               {SurfaceId("P.S2"), EventKind::Reflect, 0}}}};
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(s, lib);
  const auto& events = cs.path(PathId{0}).events;
  const std::uint32_t env = cs.environment_medium();
  REQUIRE(events[0].medium_before == env);
  REQUIRE(events[0].medium_beyond != env);  // the plate
  REQUIRE(events[1].medium_beyond == events[1].medium_after);
  REQUIRE(events[2].medium_before != env);
  REQUIRE(events[2].medium_beyond == env);
  REQUIRE(!events[0].from_inside);
  REQUIRE(!events[1].from_inside);
  REQUIRE(events[2].from_inside);
}

TEST_CASE("Fresnel reflection against the medium beyond: entrance, ghost and TIR",
          "[interactions]") {
  // N-BK7 (n = 1.5168) at normal incidence, R = ((n - 1)/(n + 1))^2 for s and p (Byrnes,
  // Eq. (6), #56); every polarization gives the same power.
  // - Reflect at the entrance surface from outside: weight = R.
  // - Ghost: Refract S1, Reflect S2 from inside (against the environment), Refract S1 back:
  //   weight = (1 - R) R (1 - R).
  // - Prism: S2 tilted by 45 deg about x; inside the glass the ray meets S2 at 45 deg, beyond
  //   the critical angle asin(1/n) = 41.25 deg, so the Fresnel reflection against the
  //   environment is total (|r_s| = |r_p| = 1): weight = 1 - R from the entrance only.
  // Tolerance 1e-12.
  const double r = (kGlass - 1.0) / (kGlass + 1.0);
  const double big_r = r * r;
  const MaterialLibrary lib;
  {
    System s = base_system();
    s.root.children = {{plate("P", 0.0, 5.0, "CONST:1.5168")}};
    s.paths = {{"entrance", false, {{SurfaceId("P.S1"), EventKind::Reflect, 0}}}};
    const RayBatch rays = trace_one(compile(s, lib), Vec3(0.0, 0.0, -1.0), Vec3(0.0, 0.0, 1.0));
    REQUIRE(rays.status()[0] == RayStatus::Alive);
    REQUIRE(std::abs(rays.weight()[0] - big_r) <= 1e-12);
    REQUIRE(std::abs(polarized_power(rays, Vec3(1.0, 0.0, 0.0)) - big_r) <= 1e-12);
  }
  {
    System s = base_system();
    s.root.children = {{plate("P", 0.0, 5.0, "CONST:1.5168")}, {detector(-2.0)}};
    s.paths = {{"ghost",
                false,
                {{SurfaceId("P.S1"), EventKind::Refract, 0},
                 {SurfaceId("P.S2"), EventKind::Reflect, 0},
                 {SurfaceId("P.S1"), EventKind::Refract, 0},
                 {SurfaceId("IMG"), EventKind::Transmit, 0}}}};
    const RayBatch rays = trace_one(compile(s, lib), Vec3(0.0, 0.0, -1.0), Vec3(0.0, 0.0, 1.0));
    const double expected = (1.0 - big_r) * big_r * (1.0 - big_r);
    REQUIRE(rays.status()[0] == RayStatus::Alive);
    REQUIRE(std::abs(rays.weight()[0] - expected) <= 1e-12);
    REQUIRE(std::abs(polarized_power(rays, Vec3(0.0, 1.0, 0.0)) - expected) <= 1e-12);
  }
  {
    System s = base_system();
    Element prism = plate("P", 0.0, 5.0, "CONST:1.5168");
    prism.surfaces[1].pose.rotation_deg[0] = Param(45.0);
    s.root.children = {{prism}};
    s.paths = {
        {"tir",
         false,
         {{SurfaceId("P.S1"), EventKind::Refract, 0}, {SurfaceId("P.S2"), EventKind::Reflect, 0}}}};
    const RayBatch rays = trace_one(compile(s, lib), Vec3(0.0, 0.0, -1.0), Vec3(0.0, 0.0, 1.0));
    REQUIRE(rays.status()[0] == RayStatus::Alive);
    REQUIRE(std::abs(rays.weight()[0] - (1.0 - big_r)) <= 1e-12);
    REQUIRE(std::abs(polarized_power(rays, Vec3(1.0, 0.0, 0.0)) - (1.0 - big_r)) <= 1e-12);
    REQUIRE(std::abs(polarized_power(rays, Vec3(0.0, 1.0, 0.0)) - (1.0 - big_r)) <= 1e-12);
  }
}

TEST_CASE("coating on refraction: AR_MGF2 on both sides of a plate", "[interactions]") {
  // Single quarter-wave layer n_c = 1.38 on n_s = 1.52 at its design wavelength 0.55 um, normal
  // incidence, lossless: r = (n_0 n_s - n_c^2)/(n_0 n_s + n_c^2) from Byrnes, Eqs. (6), (8),
  // (11), (13), (15) (derived in rtt-coating's test "quarter-wave MgF2 anti-reflection", #58;
  // R = 1.26 % in the acceptance table of docs/architecture.md), T = 1 - R. The exit surface is
  // crossed from the substrate (reversed stack, power factors from glass to vacuum); r is
  // symmetric in n_0 and n_s, so the plate transmits (1 - R)^2. The power factors cancel over
  // the plate (c = n_s/n_0 in, n_0/n_s out), so the path "in" ends inside the glass: weight =
  // 1 - R there checks them for a single surface. Tolerance 1e-12.
  rtt::coating::CoatingLibrary coatings;
  coatings.add_catalog(std::string(RTT_CATALOG_DIR) + "/coatings/demo.json");
  System s = base_system(0.55);
  Element p = plate("P", 0.0, 5.0, "CONST:1.52");
  for (auto& surf : p.surfaces) surf.interaction = rtt::model::CoatingRef{"DEMO:AR_MGF2"};
  s.root.children = {{p}, {detector(10.0)}};
  s.paths = {{"main", true, {}}, {"in", false, {{SurfaceId("P.S1"), EventKind::Refract, 0}}}};
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(s, lib, coatings);
  const RayBatch rays = trace_one(cs, Vec3(0.0, 0.0, -1.0), Vec3(0.0, 0.0, 1.0));
  const double nc2 = 1.38 * 1.38;
  const double big_r = std::pow((1.52 - nc2) / (1.52 + nc2), 2);
  REQUIRE(rays.status()[0] == RayStatus::Alive);
  REQUIRE(std::abs(rays.weight()[0] - (1.0 - big_r) * (1.0 - big_r)) <= 1e-12);
  REQUIRE(std::abs(polarized_power(rays, Vec3(1.0, 0.0, 0.0)) - (1.0 - big_r) * (1.0 - big_r)) <=
          1e-12);
  const RayBatch in = trace_one(cs, Vec3(0.0, 0.0, -1.0), Vec3(0.0, 0.0, 1.0), PathId{1});
  REQUIRE(in.status()[0] == RayStatus::Alive);
  REQUIRE(std::abs(in.weight()[0] - (1.0 - big_r)) <= 1e-12);
}

TEST_CASE("reference example m3/absorbing_ar_plate: AR coating and volume absorption",
          "[interactions]") {
  // tests/reference/m3/absorbing_ar_plate.rtt.json: plate n = 1.52 + 1e-6 i, d = 10 mm,
  // DEMO:AR_MGF2 on both sides, lambda = 0.55 um (design wavelength), vacuum, axial ray.
  // Expected: weight = (1 - R)^2 exp(-4 pi kappa d / lambda) with R of the quarter-wave layer on
  // the real substrate, r = (n_0 n_s - n_c^2)/(n_0 n_s + n_c^2) (Byrnes, Eqs. (6)-(15), derived
  // in rtt-coating's test, #58), and the absorption of Byrnes, Eqs. (1), (2), (17), (18).
  // kappa changes r, t and the power factors only in second order: at normal incidence they are
  // real functions of n_s, a purely imaginary change i kappa changes them by i kappa f' (f'
  // real), so |r|^2, |t|^2 and Re(n) change by O(kappa^2) = 1e-12. Tolerance 1e-10 relative
  // (second-order terms with coefficients up to ~10, plus rounding).
  rtt::coating::CoatingLibrary coatings;
  coatings.add_catalog(std::string(RTT_CATALOG_DIR) + "/coatings/demo.json");
  const MaterialLibrary lib;
  const System s =
      rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m3/absorbing_ar_plate.rtt.json");
  const CompiledSystem cs = compile(s, lib, coatings);
  const RayBatch rays = trace_one(cs, Vec3(0.0, 0.0, 0.0), Vec3(0.0, 0.0, 1.0));
  const double nc2 = 1.38 * 1.38;
  const double big_r = std::pow((1.52 - nc2) / (1.52 + nc2), 2);
  const double absorption = std::exp(-4.0 * kPi * 1e-6 * 10.0 / 0.55e-3);
  const double expected = (1.0 - big_r) * (1.0 - big_r) * absorption;
  REQUIRE(rays.status()[0] == RayStatus::Alive);
  REQUIRE(rays.last_surface()[0] == 2);
  REQUIRE(std::abs(rays.weight()[0] - expected) <= 1e-10 * expected);
}

TEST_CASE("reflection at an ideal anti-reflection surface vanishes", "[interactions]") {
  // Decided for #61 (ADR 0021): like a Fresnel surface without a change of medium, Reflect at
  // IdealAntiReflection is a vanishing, not an impossible reflection: amplitude 0, weight 0,
  // status Alive (a ghost path over an ideal AR surface simply has weight 0).
  System s = base_system();
  Element p = plate("P", 0.0, 5.0, "CONST:1.5");
  p.surfaces[0].interaction = rtt::model::IdealAntiReflection{};
  s.root.children = {{p}};
  s.paths = {{"ghost", false, {{SurfaceId("P.S1"), EventKind::Reflect, 0}}}};
  const MaterialLibrary lib;
  const RayBatch rays = trace_one(compile(s, lib), Vec3(0.0, 0.0, -1.0), Vec3(0.0, 0.0, 1.0));
  REQUIRE(rays.status()[0] == RayStatus::Alive);
  REQUIRE(rays.weight()[0] == 0.0);
  REQUIRE(std::abs(rays.dir_z()[0] + 1.0) <= 1e-15);  // reflected
}

TEST_CASE("ideal retarder through the tracer: quarter-wave plate at 45 deg gives circular light",
          "[interactions]") {
  // x polarizer, then IdealRetarder with retardance 1/4 wave and fast axis (1, 0, 0) in element
  // coordinates, element rotated by 45 deg about z (axis rotated by compile, ADR 0021). Linear
  // light at 45 deg to the fast axis leaves circularly polarized: |S3| = S0 (docs/architecture.md,
  // acceptance table), S0 = 1 for the incident x state; unpolarized weight 1/2. Tolerance 1e-12.
  System s = base_system();
  s.root.children = {{thin("POL", 0.0, rtt::model::IdealPolarizer{{1.0, 0.0, 0.0}, 0.0})},
                     {thin("QWP", 5.0, rtt::model::IdealRetarder{{1.0, 0.0, 0.0}, 0.25}, 45.0)},
                     {detector(10.0)}};
  const MaterialLibrary lib;
  const RayBatch rays = trace_one(compile(s, lib), Vec3(0.0, 0.0, -1.0), Vec3(0.0, 0.0, 1.0));
  const Vec3 k(0.0, 0.0, 1.0);
  const CVec3 e = rays.prt_matrix(0) * Vec3(1.0, 0.0, 0.0).cast<Cx>();
  const auto st = rtt::polar::stokes<double>(e, Vec3(1.0, 0.0, 0.0), k);
  REQUIRE(std::abs(st.s0 - 1.0) <= 1e-12);
  REQUIRE(std::abs(std::abs(st.s3) - st.s0) <= 1e-12);
  REQUIRE(std::abs(rays.weight()[0] - 0.5) <= 1e-12);
}
