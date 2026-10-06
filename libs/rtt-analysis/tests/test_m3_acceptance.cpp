// M3 acceptance (#63): the M3 rows of the table "Validierung und Tests" in docs/architecture.md
// end to end through the tracer, with the systems of tests/reference/m3 (read with rtt-io). All
// media are vacuum or constant indices, so neither dispersion nor Ciddor air enters the expected
// values.
//
// Observables are convention-free (docs/architecture.md, "Retardance einzelner Reflexionen"):
// powers |P E|^2, Stokes quantities of the output state, Malus transmission, circularity after a
// quarter-wave plate, the rotation angle after a half-wave plate and, for metal mirrors, the
// ellipsometric ratio against the ideal conductor and the total of a mirror sequence. P is
// power-normalised (ADR 0021), so |P E|^2 is the reflected or transmitted power fraction of the
// unit input state E directly, including the n cos factors of a transmission.
//
// Expected values are closed forms evaluated here, independently of the code under test:
// Fresnel amplitudes after S. J. Byrnes, arXiv:1603.02720v5, Eq. (6) (Convention A, as in
// docs/architecture.md) with the root branch of Appendix D, powers after Eqs. (21)-(23), the
// ellipsometric ratio after Eq. (16); Jones matrices of ideal elements after W.-S. T. Lam,
// dissertation, Eqs. (3.4)/(3.6), Table 2.1 (handedness), Fig. 1.1 (quarter wave) and p. 123
// (half wave); see docs/quellen.md.
//
// Tolerances (a priori): the expected values are double closed forms (a few ulp, ~1e-16); the
// tracer multiplies at most five complex 3x3 matrices, normalises directions and intersects plane
// surfaces, each ~1e-15 relative. 1e-12 therefore leaves about three orders of magnitude; the
// anti-reflection coating and the metal mirror use the 1e-10 of the acceptance table.

#include <Eigen/Geometry>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <complex>
#include <cstdint>
#include <memory>
#include <numbers>
#include <optional>
#include <string>
#include <variant>
#include <vector>

#include "rtt/coating/catalog.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
#include "rtt/polar/ideal.hpp"
#include "rtt/polar/prt.hpp"
#include "rtt/polar/prt_analysis.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/sequential.hpp"
#include "rtt/trace/sources.hpp"

using Cx = std::complex<double>;
using CVec3 = Eigen::Vector3cd;
using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::math::Vec3;
using rtt::model::Element;
using rtt::model::Param;
using rtt::model::System;
using rtt::trace::Aiming;
using rtt::trace::RayBatch;
using rtt::trace::RayStatus;
using rtt::trace::SequentialTracer;

namespace {

constexpr double kPi = std::numbers::pi;
constexpr double kDeg = kPi / 180.0;
constexpr double kGlass = 1.5168;  // N-BK7 at d (acceptance table)

System load(const std::string& name) {
  return rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m3/" + name);
}

/// Compiles with the demo and the M3 test coating catalogues (vacuum and CONST media only).
CompiledSystem build(const System& s) {
  const rtt::material::MaterialLibrary materials;
  const auto coatings = std::make_unique<rtt::coating::CoatingLibrary>();
  coatings->add_catalog(std::string(RTT_CATALOG_DIR) + "/coatings/demo.json");
  coatings->add_catalog(std::string(RTT_CATALOG_DIR) + "/coatings/m3.json");
  return rtt::compile::compile(s, materials, *coatings);
}

PathId path_named(const CompiledSystem& cs, const std::string& name) {
  for (std::size_t i = 0; i < cs.paths().size(); ++i) {
    if (cs.paths()[i].name == name) return PathId{static_cast<std::uint32_t>(i)};
  }
  FAIL("no path " << name);
  return PathId{0};
}

Element& element_named(System& s, const std::string& name) {
  for (auto& child : s.root.children) {
    if (auto* e = std::get_if<Element>(&child.value); e && e->name == name) return *e;
  }
  FAIL("no element " << name);
  return std::get<Element>(s.root.children.front().value);
}

/// Fresnel amplitudes of an interface n_i -> n_t for the tangential invariant xi = n_i sin(theta)
/// (Byrnes, Eq. (6), Convention A; q = n cos(theta) = sqrt(n^2 - xi^2) with Im q >= 0, and
/// Re q >= 0 if Im q = 0, Appendix D).
struct Amplitudes {
  Cx rs;
  Cx rp;
  Cx ts;
  Cx tp;
};
Cx normal_component(Cx n, double xi) {
  Cx v = std::sqrt(n * n - xi * xi);
  if (v.imag() < 0.0 || (v.imag() == 0.0 && v.real() < 0.0)) v = -v;
  return v;
}
Amplitudes fresnel_reference(Cx n_i, Cx n_t, double xi) {
  const Cx q_i = normal_component(n_i, xi);
  const Cx q_t = normal_component(n_t, xi);
  const Cx cos_i = q_i / n_i;
  const Cx cos_t = q_t / n_t;
  const Cx dp = n_t * cos_i + n_i * cos_t;
  return {(q_i - q_t) / (q_i + q_t), (n_t * cos_i - n_i * cos_t) / dp, 2.0 * q_i / (q_i + q_t),
          2.0 * n_i * cos_i / dp};
}

/// Layer of a stack: index and physical thickness in um.
struct Layer {
  Cx index;
  double thickness_um;
};

/// Reflectances R_s, R_p of a layer stack from the ambient n_0 to the substrate n_sub, after the
/// transfer-matrix method of Byrnes (interface form): M = (1/t_01) [[1, r_01], [r_01, 1]] times,
/// for each layer n, diag(e^(-i delta_n), e^(i delta_n)) (1/t_(n,n+1)) [[1, r_(n,n+1)],
/// [r_(n,n+1), 1]], with delta_n = 2 pi q_n d_n / lambda and q_n = n_n cos(theta_n)
/// (Eqs. (8), (11), (13)); r = M_10 / M_00 (Eq. (15)). Interface amplitudes from Eq. (6). Own
/// implementation, independent of rtt-coating.
std::pair<double, double> stack_reflectance(
    Cx n_0, const std::vector<Layer>& layers, Cx n_sub, double xi, double wavelength_um) {
  const auto interface = [](Cx r, Cx t) {
    Eigen::Matrix2cd m;
    m << 1.0, r, r, 1.0;
    return Eigen::Matrix2cd(m / t);
  };
  std::vector<Cx> n{n_0};
  for (const Layer& l : layers) n.push_back(l.index);
  n.push_back(n_sub);
  const Amplitudes first = fresnel_reference(n[0], n[1], xi);
  Eigen::Matrix2cd ms = interface(first.rs, first.ts);
  Eigen::Matrix2cd mp = interface(first.rp, first.tp);
  const Cx i(0.0, 1.0);
  for (std::size_t j = 0; j < layers.size(); ++j) {
    const Cx delta =
        2.0 * kPi * normal_component(layers[j].index, xi) * layers[j].thickness_um / wavelength_um;
    Eigen::Matrix2cd phase = Eigen::Matrix2cd::Zero();
    phase(0, 0) = std::exp(-i * delta);
    phase(1, 1) = std::exp(i * delta);
    const Amplitudes next = fresnel_reference(n[j + 1], n[j + 2], xi);
    ms = ms * phase * interface(next.rs, next.ts);
    mp = mp * phase * interface(next.rp, next.tp);
  }
  return {std::norm(ms(1, 0) / ms(0, 0)), std::norm(mp(1, 0) / mp(0, 0))};
}

/// Rays of one field (hexapolar, one ring: 7 rays) with real aiming, traced on `path`. The start
/// directions are returned in `k_in`.
RayBatch trace_field(const CompiledSystem& cs,
                     PathId path,
                     std::uint16_t field,
                     std::vector<Vec3>* k_in = nullptr) {
  const std::vector<std::uint16_t> fields{field};
  RayBatch rays =
      rtt::trace::make_rays(cs, path, fields, 0, rtt::trace::HexapolarPupil{1}, Aiming::Real);
  if (k_in != nullptr) {
    k_in->clear();
    for (std::size_t i = 0; i < rays.size(); ++i) {
      k_in->emplace_back(rays.dir_x()[i], rays.dir_y()[i], rays.dir_z()[i]);
    }
  }
  [[maybe_unused]] const auto stats = SequentialTracer().trace(cs, path, rays);
  return rays;
}

/// One ray from `start` along `dir` on `path` (tilted systems, which the paraxial aiming of
/// rtt::trace::make_rays does not accept).
RayBatch trace_one(const CompiledSystem& cs, PathId path, const Vec3& start, const Vec3& dir) {
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

CVec3 out(const RayBatch& rays, std::size_t i, const CVec3& e) {
  return rays.prt_matrix(i) * e;
}
CVec3 out(const RayBatch& rays, std::size_t i, const Vec3& e) {
  return out(rays, i, CVec3(e.cast<Cx>()));
}
double power(const RayBatch& rays, std::size_t i, const Vec3& e) {
  return out(rays, i, e).squaredNorm();
}
Vec3 direction(const RayBatch& rays, std::size_t i) {
  return {rays.dir_x()[i], rays.dir_y()[i], rays.dir_z()[i]};
}

}  // namespace

TEST_CASE("M3 Fresnel at N-BK7: R = ((n - 1)/(n + 1))^2 at normal incidence", "[m3]") {
  // tests/reference/m3/fresnel_bk7.rtt.json, field 0. Byrnes, Eq. (6) at normal incidence:
  // r_s = -r_p = (1 - n)/(1 + n); R = 4.22 % for every polarization (acceptance table).
  const CompiledSystem cs = build(load("fresnel_bk7.rtt.json"));
  const double r = (kGlass - 1.0) / (kGlass + 1.0);
  const double big_r = r * r;
  REQUIRE(std::abs(big_r - 0.0422) < 5e-5);
  const RayBatch reflected = trace_field(cs, path_named(cs, "reflect S1"), 0);
  const RayBatch transmitted = trace_field(cs, path_named(cs, "main"), 0);
  REQUIRE(reflected.size() == 7);
  for (std::size_t i = 0; i < reflected.size(); ++i) {
    INFO("ray " << i);
    REQUIRE(reflected.status()[i] == RayStatus::Alive);
    REQUIRE(transmitted.status()[i] == RayStatus::Alive);
    for (const Vec3& e : {Vec3(1.0, 0.0, 0.0), Vec3(0.0, 1.0, 0.0)}) {
      REQUIRE(std::abs(power(reflected, i, e) - big_r) <= 1e-12);
      REQUIRE(std::abs(power(transmitted, i, e) - (1.0 - big_r) * (1.0 - big_r)) <= 1e-12);
    }
    REQUIRE(std::abs(reflected.weight()[i] - big_r) <= 1e-12);
  }
}

TEST_CASE("M3 Brewster angle atan(n): R_p = 0, R_s from Fresnel", "[m3]") {
  // Field 1 of fresnel_bk7: theta_B = atan(n) = 56.6038... deg. In Byrnes, Eq. (6),
  // r_p = (n cos(theta_i) - cos(theta_t)) / (n cos(theta_i) + cos(theta_t)) vanishes for
  // n cos(theta_i) = cos(theta_t); with Snell (sin(theta_i) = n sin(theta_t)) that is
  // tan(theta_i) = n. The field lies in the y-z plane, so s = x and p_in = k_in x s.
  const CompiledSystem cs = build(load("fresnel_bk7.rtt.json"));
  const double theta_b = std::atan(kGlass);
  const Amplitudes a = fresnel_reference(1.0, kGlass, std::sin(theta_b));
  REQUIRE(std::abs(a.rp) <= 1e-15);  // the reference itself
  std::vector<Vec3> k_in;
  const RayBatch reflected = trace_field(cs, path_named(cs, "reflect S1"), 1, &k_in);
  for (std::size_t i = 0; i < reflected.size(); ++i) {
    INFO("ray " << i);
    REQUIRE(reflected.status()[i] == RayStatus::Alive);
    REQUIRE(std::abs(k_in[i].y() / k_in[i].z() - kGlass) <= 1e-12);  // incidence at theta_B
    const Vec3 e_s(1.0, 0.0, 0.0);
    const Vec3 e_p = k_in[i].cross(e_s).normalized();
    REQUIRE(power(reflected, i, e_p) <= 1e-12);
    REQUIRE(std::abs(power(reflected, i, e_s) - std::norm(a.rs)) <= 1e-12);
    REQUIRE(std::norm(a.rs) > 0.1);  // R_s is far from zero there
  }
}

TEST_CASE("M3 total internal reflection: |r| = 1 and the s/p phase jump", "[m3]") {
  // tests/reference/m3/tir_prism.rtt.json: normal entrance (t_s = t_p, power 1 - R each), then
  // the face P.S2 tilted by 45 deg about x, met inside at 45 deg > asin(1/n) = 41.25 deg. With
  // xi = n sin(45 deg) > 1, Byrnes, Eq. (6) with q_t = i sqrt(xi^2 - 1) (Appendix D) gives
  // |r_s| = |r_p| = 1 and a phase difference Delta = arg(r_p / r_s). The input (x + y)/sqrt(2)
  // (s = x, p = y) leaves with S1 = 0, |S2| = S0 |cos(Delta)| and |S3| = S0 |sin(Delta)|; the
  // magnitudes do not depend on the sign convention of r_p.
  const CompiledSystem cs = build(load("tir_prism.rtt.json"));
  const double r = (kGlass - 1.0) / (kGlass + 1.0);
  const double entrance = 1.0 - r * r;
  const Amplitudes a = fresnel_reference(kGlass, 1.0, kGlass * std::sin(45.0 * kDeg));
  REQUIRE(std::abs(std::abs(a.rs) - 1.0) <= 1e-15);
  REQUIRE(std::abs(std::abs(a.rp) - 1.0) <= 1e-15);
  const double delta = std::arg(a.rp / a.rs);
  REQUIRE(std::abs(std::sin(delta)) > 0.5);  // a real phase jump, not 0 or pi
  const RayBatch rays =
      trace_one(cs, path_named(cs, "tir"), Vec3(0.0, 0.0, -1.0), Vec3(0.0, 0.0, 1.0));
  REQUIRE(rays.status()[0] == RayStatus::Alive);
  REQUIRE(std::abs(power(rays, 0, Vec3(1.0, 0.0, 0.0)) - entrance) <= 1e-12);
  REQUIRE(std::abs(power(rays, 0, Vec3(0.0, 1.0, 0.0)) - entrance) <= 1e-12);
  REQUIRE(std::abs(rays.weight()[0] - entrance) <= 1e-12);
  const CVec3 e = out(rays, 0, Vec3(1.0, 1.0, 0.0).normalized());
  const auto st = rtt::polar::stokes<double>(e, Vec3(1.0, 0.0, 0.0), direction(rays, 0));
  REQUIRE(std::abs(st.s0 - entrance) <= 1e-12);
  REQUIRE(std::abs(st.s1) <= 1e-12);
  REQUIRE(std::abs(std::abs(st.s2) - st.s0 * std::abs(std::cos(delta))) <= 1e-12);
  REQUIRE(std::abs(std::abs(st.s3) - st.s0 * std::abs(std::sin(delta))) <= 1e-12);
  // Convention check (Convention A and the Stokes definition of docs/architecture.md): with
  // e1 = x = s and e2 = k_out x e1 = p_out = k_out x s, E1 ~ r_s and E2 ~ r_p with the same
  // positive factor, so S2 = S0 cos(Delta) and S3 = -S0 sin(Delta).
  REQUIRE(std::abs(st.s2 - st.s0 * std::cos(delta)) <= 1e-12);
  REQUIRE(std::abs(st.s3 + st.s0 * std::sin(delta)) <= 1e-12);
}

TEST_CASE("M3 total internal reflection: the critical angle asin(1/n)", "[m3]") {
  // P.S2 of tir_prism tilted to the angle beta with sin(beta) = (1/n)(1 +- 1e-12): a Refract at
  // P.S2 must stop the ray as Tir above the critical angle and pass below it. The ray enters
  // P.S1 at normal incidence along +z, so inside it meets P.S2 at beta; direction and rotation
  // carry rounding errors of ~1e-16 relative (normalisation, deg -> rad, rotation matrix), so
  // xi = n sin(beta) is off by ~1e-16 relative, far below the margin 1e-12.
  const double sin_c = 1.0 / kGlass;
  for (const double side : {+1.0, -1.0}) {
    INFO("sin(beta) = (1/n)(1 " << (side > 0 ? "+" : "-") << " 1e-12)");
    System s = load("tir_prism.rtt.json");
    Element& prism = element_named(s, "P");
    const double beta = std::asin(sin_c * (1.0 + side * 1e-12));
    prism.surfaces[1].pose.rotation_deg[0] = Param(beta / kDeg);
    s.paths = {{"exit",
                false,
                {{rtt::model::SurfaceId("P.S1"), rtt::model::EventKind::Refract, 0},
                 {rtt::model::SurfaceId("P.S2"), rtt::model::EventKind::Refract, 0}}}};
    const CompiledSystem cs = build(s);
    const RayBatch rays = trace_one(cs, PathId{0}, Vec3(0.0, 0.0, -1.0), Vec3(0.0, 0.0, 1.0));
    REQUIRE(rays.status()[0] == (side > 0 ? RayStatus::Tir : RayStatus::Alive));
  }
}

TEST_CASE("M3 quarter-wave MgF2 anti-reflection: R = 1.26 % at the design wavelength", "[m3]") {
  // tests/reference/m3/ar_mgf2.rtt.json: DEMO:AR_MGF2 (n_c = 1.38, one quarter wave at 0.55 um)
  // on n_s = 1.52 in vacuum (n_0 = 1), normal incidence: r = (n_0 n_s - n_c^2)/(n_0 n_s + n_c^2)
  // from Byrnes, Eqs. (6), (8), (11), (13), (15) with the layer phase pi/2 (derived in
  // rtt-coating's test "quarter-wave MgF2 anti-reflection", #58). Tolerance 1e-10 (table).
  const CompiledSystem cs = build(load("ar_mgf2.rtt.json"));
  const double nc2 = 1.38 * 1.38;
  const double big_r = std::pow((1.52 - nc2) / (1.52 + nc2), 2);
  REQUIRE(std::abs(big_r - 0.0126) < 5e-5);
  const RayBatch reflected = trace_field(cs, path_named(cs, "reflect S1"), 0);
  for (std::size_t i = 0; i < reflected.size(); ++i) {
    INFO("ray " << i);
    REQUIRE(reflected.status()[i] == RayStatus::Alive);
    REQUIRE(std::abs(power(reflected, i, Vec3(1.0, 0.0, 0.0)) - big_r) <= 1e-10);
    REQUIRE(std::abs(power(reflected, i, Vec3(0.0, 1.0, 0.0)) - big_r) <= 1e-10);
  }
}

TEST_CASE("M3 lossless stack: R + T = 1 for s and p, R from the transfer matrix", "[m3]") {
  // tests/reference/m3/lossless_stack.rtt.json: M3:HL5 (five lossless quarter-wave layers H L H
  // L H at 0.55 um, H = 2.35, L = 1.38, tests/catalogs/coatings/m3.json) on n = 1.52; the exit
  // face is an ideal anti-reflection surface (amplitude 1), so the path "main" transmits exactly
  // the stack's T. Byrnes, Eqs. (21)-(23): R + T = 1 for a lossless incidence medium and lossless
  // layers, for s and p, at 0, 30 and 60 deg. Since that holds for every lossless stack, R_s and
  // R_p are also compared with the transfer matrix of stack_reflectance (own implementation of
  // Byrnes, Eqs. (6), (8), (11), (13), (15); second review of #87): it fixes the layer phase
  // 2 pi q d / lambda with q = n cos(theta), not n. Tolerance 1e-12 for both: about 15 complex
  // 2x2 products, exponentials and square roots on numbers of order 1 (|r| <= 1, |1/t| < 10 for
  // these indices), each ~1e-16 relative, in the tracer and in the reference alike.
  const CompiledSystem cs = build(load("lossless_stack.rtt.json"));
  const std::vector<Layer> hl5{{2.35, 0.55 / (4.0 * 2.35)},
                               {1.38, 0.55 / (4.0 * 1.38)},
                               {2.35, 0.55 / (4.0 * 2.35)},
                               {1.38, 0.55 / (4.0 * 1.38)},
                               {2.35, 0.55 / (4.0 * 2.35)}};
  // Plausibility of the reference itself against the values of the second review (rounded).
  struct Check {
    double theta_deg;
    double rs;
    double rp;
    double tol;
  };
  for (const Check c : {Check{0.0, 0.877245, 0.877245, 5e-7}, Check{30.0, 0.9072, 0.8240, 5e-5},
                        Check{60.0, 0.9510, 0.4358, 5e-5}}) {
    const auto [rs, rp] = stack_reflectance(1.0, hl5, 1.52, std::sin(c.theta_deg * kDeg), 0.55);
    INFO("theta " << c.theta_deg << ": R_s " << rs << ", R_p " << rp);
    REQUIRE(std::abs(rs - c.rs) <= c.tol);
    REQUIRE(std::abs(rp - c.rp) <= c.tol);
  }
  for (std::uint16_t field = 0; field < 3; ++field) {
    std::vector<Vec3> k_in;
    const RayBatch reflected = trace_field(cs, path_named(cs, "reflect S1"), field, &k_in);
    const RayBatch transmitted = trace_field(cs, path_named(cs, "main"), field);
    for (std::size_t i = 0; i < reflected.size(); ++i) {
      INFO("field " << field << ", ray " << i);
      REQUIRE(reflected.status()[i] == RayStatus::Alive);
      REQUIRE(transmitted.status()[i] == RayStatus::Alive);
      const Vec3 e_s(1.0, 0.0, 0.0);
      const Vec3 e_p = k_in[i].cross(e_s).normalized();
      const double r_s = power(reflected, i, e_s);
      const double r_p = power(reflected, i, e_p);
      REQUIRE(std::abs(r_s + power(transmitted, i, e_s) - 1.0) <= 1e-12);
      REQUIRE(std::abs(r_p + power(transmitted, i, e_p) - 1.0) <= 1e-12);
      const double xi = std::hypot(k_in[i].x(), k_in[i].y());  // vacuum, normal along z
      const auto [ref_s, ref_p] = stack_reflectance(1.0, hl5, 1.52, xi, 0.55);
      REQUIRE(std::abs(r_s - ref_s) <= 1e-12);
      REQUIRE(std::abs(r_p - ref_p) <= 1e-12);
      REQUIRE(r_s > 0.5);                                   // a real reflector, not a trivial case
      if (field == 2) REQUIRE(std::abs(r_s - r_p) > 1e-2);  // s and p differ at 60 deg
    }
  }
}

TEST_CASE("M3 two polarizers: Malus T = cos^2(theta)", "[m3]") {
  // tests/reference/m3/malus.rtt.json (second polarizer at 30 deg), and the second polarizer
  // turned to 0, 45, 60 and 90 deg. An ideal polarizer with axis a is the projection a a^T on
  // the transverse plane (Lam, Eqs. (3.4)/(3.6), Jones matrix diag(1, 0)); the x state after
  // both: |a . x|^2 = cos^2(theta). Unpolarized: weight = cos^2(theta)/2. Tolerance 1e-12.
  for (const double theta : {30.0, 0.0, 45.0, 60.0, 90.0}) {
    INFO("theta = " << theta << " deg");
    System s = load("malus.rtt.json");
    element_named(s, "polarizer 2").pose.rotation_deg[2] = Param(theta);
    const CompiledSystem cs = build(s);
    const RayBatch rays = trace_field(cs, PathId{0}, 0);
    const double c2 = std::pow(std::cos(theta * kDeg), 2);
    for (std::size_t i = 0; i < rays.size(); ++i) {
      REQUIRE(rays.status()[i] == RayStatus::Alive);
      REQUIRE(std::abs(power(rays, i, Vec3(1.0, 0.0, 0.0)) - c2) <= 1e-12);
      REQUIRE(std::abs(rays.weight()[i] - 0.5 * c2) <= 1e-12);
    }
  }
}

TEST_CASE("M3 quarter-wave plate at 45 deg: linear -> circular", "[m3]") {
  // tests/reference/m3/quarter_wave.rtt.json: x polarizer, then a retarder of 1/4 wave with its
  // fast axis at 45 deg (Lam, Fig. 1.1, p. 48). The x state leaves circular: S1 = S2 = 0,
  // |S3| = S0 = 1 (convention-free). Separately, with the documented convention (fields
  // exp(i(k.r - omega t)), S3 > 0 right circular after Lam, Table 2.1; docs/architecture.md,
  // "Händigkeit und Stokes"): retardation of the slow axis at -45 deg gives
  // E ~ (1, -i)/sqrt(2) in (x, y), i.e. S3 = +S0. Tolerance 1e-12.
  const CompiledSystem cs = build(load("quarter_wave.rtt.json"));
  const RayBatch rays = trace_field(cs, PathId{0}, 0);
  for (std::size_t i = 0; i < rays.size(); ++i) {
    INFO("ray " << i);
    REQUIRE(rays.status()[i] == RayStatus::Alive);
    const auto st = rtt::polar::stokes<double>(out(rays, i, Vec3(1.0, 0.0, 0.0)),
                                               Vec3(1.0, 0.0, 0.0), direction(rays, i));
    REQUIRE(std::abs(st.s0 - 1.0) <= 1e-12);
    REQUIRE(std::abs(st.s1) <= 1e-12);
    REQUIRE(std::abs(st.s2) <= 1e-12);
    REQUIRE(std::abs(std::abs(st.s3) - st.s0) <= 1e-12);
    REQUIRE(std::abs(st.s3 - st.s0) <= 1e-12);  // convention check: right circular
  }
}

TEST_CASE("M3 half-wave plate at theta: rotation by 2 theta", "[m3]") {
  // tests/reference/m3/half_wave.rtt.json (fast axis at 22.5 deg) and the plate turned to 10 and
  // 37 deg: a half-wave retarder with fast axis at theta reflects linear polarization at its
  // axis (Lam, p. 123), so the x state leaves linear at 2 theta: S3 = 0, S0 = 1 and
  // atan2(S2, S1)/2 = 2 theta (mod pi) against the x axis. Tolerance 1e-12.
  for (const double theta : {22.5, 10.0, 37.0}) {
    INFO("theta = " << theta << " deg");
    System s = load("half_wave.rtt.json");
    element_named(s, "half wave plate").pose.rotation_deg[2] = Param(theta);
    const CompiledSystem cs = build(s);
    const RayBatch rays = trace_field(cs, PathId{0}, 0);
    for (std::size_t i = 0; i < rays.size(); ++i) {
      REQUIRE(rays.status()[i] == RayStatus::Alive);
      const auto st = rtt::polar::stokes<double>(out(rays, i, Vec3(1.0, 0.0, 0.0)),
                                                 Vec3(1.0, 0.0, 0.0), direction(rays, i));
      REQUIRE(std::abs(st.s0 - 1.0) <= 1e-12);
      REQUIRE(std::abs(st.s3) <= 1e-12);
      const double angle = 0.5 * std::atan2(st.s2, st.s1);
      REQUIRE(std::abs(std::remainder(angle - 2.0 * theta * kDeg, kPi)) <= 1e-12);
    }
  }
}

TEST_CASE("M3 metal mirror at 45 deg: ellipsometric ratio and diattenuation from Fresnel", "[m3]") {
  // tests/reference/m3/metal_periscope.rtt.json, path "single": mirror with the invented test
  // metal n = 1.2 + 7.26 i (as in the #56 tests; no material data), tilted 45 deg about x, ray
  // along +z, so s = x and p_in = y. Convention-free reference: the same mirror without material
  // (ideal conductor, docs/architecture.md: E_r = (-I + 2 N N^T) E, no retardance). Relative to
  // it the metal multiplies the s component by rho_s = r_s / r_s,ideal = -r_s and the p
  // component by rho_p = r_p / r_p,ideal = r_p (Convention A; the ratios do not depend on the
  // convention): the ellipsometric ratio tan(psi) e^(i Delta) = rho_p / rho_s = -r_p / r_s
  // (Byrnes, Eq. (16)). Checked: the whole complex output field, R_s = |r_s|^2, R_p = |r_p|^2
  // and the diattenuation (R_s - R_p)/(R_s + R_p) from P (equal to the power diattenuation for a
  // reflection, rtt/polar/prt_analysis.hpp). Tolerance 1e-10 (table).
  const Cx metal(1.2, 7.26);
  const Amplitudes a = fresnel_reference(1.0, metal, std::sin(45.0 * kDeg));
  const Cx rho_s = -a.rs;
  const Cx rho_p = a.rp;
  REQUIRE(std::abs(std::norm(a.rs) - std::norm(a.rp)) > 1e-3);  // a real diattenuation
  REQUIRE(std::abs(std::sin(std::arg(rho_p / rho_s))) > 1e-2);  // and a real retardance

  System s = load("metal_periscope.rtt.json");
  const CompiledSystem metal_cs = build(s);
  for (auto& child : s.root.children) {
    if (auto* e = std::get_if<Element>(&child.value);
        e && e->kind == rtt::model::ElementKind::Mirror) {
      e->material = std::nullopt;  // ideal conductor
    }
  }
  const CompiledSystem ideal_cs = build(s);
  const Vec3 start(0.0, 0.0, 0.0);
  const Vec3 k_in(0.0, 0.0, 1.0);
  const RayBatch m = trace_one(metal_cs, path_named(metal_cs, "single"), start, k_in);
  const RayBatch ideal = trace_one(ideal_cs, path_named(ideal_cs, "single"), start, k_in);
  REQUIRE(m.status()[0] == RayStatus::Alive);
  REQUIRE(ideal.status()[0] == RayStatus::Alive);
  const Vec3 k_out = direction(m, 0);
  REQUIRE((k_out - direction(ideal, 0)).norm() <= 1e-15);
  REQUIRE((k_out - Vec3(0.0, 1.0, 0.0)).norm() <= 1e-12);

  const Vec3 s_hat(1.0, 0.0, 0.0);
  const std::vector<CVec3> inputs{CVec3(1.0, 0.0, 0.0), CVec3(0.0, 1.0, 0.0),
                                  CVec3(CVec3(Cx(1.0, 0.0), Cx(0.0, 1.0), 0.0) / std::sqrt(2.0)),
                                  CVec3(Cx(0.6, 0.0), Cx(0.0, -0.8), 0.0)};
  for (const CVec3& e_in : inputs) {
    const CVec3 e_ideal = out(ideal, 0, e_in);
    const CVec3 s_part = s_hat.cast<Cx>() * s_hat.cast<Cx>().dot(e_ideal);
    const CVec3 expected = rho_s * s_part + rho_p * (e_ideal - s_part);
    REQUIRE((out(m, 0, e_in) - expected).norm() <= 1e-10);
  }
  // Ellipsometric ratio measured from the output states for the input (x + y)/sqrt(2).
  const CVec3 e45 = CVec3(1.0, 1.0, 0.0) / std::sqrt(2.0);
  const CVec3 em = out(m, 0, e45);
  const CVec3 ei = out(ideal, 0, e45);
  const Vec3 p_hat = k_out.cross(s_hat);
  const Cx measured = (p_hat.cast<Cx>().dot(em) / s_hat.cast<Cx>().dot(em)) /
                      (p_hat.cast<Cx>().dot(ei) / s_hat.cast<Cx>().dot(ei));
  REQUIRE(std::abs(measured - rho_p / rho_s) <= 1e-10);

  const double big_rs = std::norm(a.rs);
  const double big_rp = std::norm(a.rp);
  REQUIRE(std::abs(power(m, 0, Vec3(1.0, 0.0, 0.0)) - big_rs) <= 1e-10);
  REQUIRE(std::abs(power(m, 0, Vec3(0.0, 1.0, 0.0)) - big_rp) <= 1e-10);
  const auto d = rtt::polar::diattenuation(m.prt_matrix(0), k_in, k_out);
  REQUIRE(std::abs(d.value - std::abs(big_rs - big_rp) / (big_rs + big_rp)) <= 1e-10);
  // With the documented convention (docs/architecture.md, "Retardance einzelner Reflexionen":
  // the physical retardance of a single reflection, from Q^-1 P, is |arg(-r_p / r_s)|), also the
  // retardance of this one reflection: the ellipsometric Delta, about 10.8 deg here.
  const Vec3 normal(0.0, -std::sin(45.0 * kDeg), std::cos(45.0 * kDeg));  // R_x(45 deg) e_z
  const auto q = rtt::polar::geometric_transform<double>(k_in, k_out, normal, true);
  const auto phys = rtt::polar::physical_retardance(m.prt_matrix(0), q, k_in);
  REQUIRE(std::abs(phys.value - std::abs(std::arg(-a.rp / a.rs))) <= 1e-10);
}

TEST_CASE("M3 metal periscope: two 45 deg reflections in total", "[m3]") {
  // Path "periscope" of metal_periscope: two parallel metal mirrors in the same plane of
  // incidence (s = x for both), the ray leaves along +z again. Against the ideal periscope the
  // s and p components get rho_s^2 and rho_p^2 (see the single mirror). The path maps k onto
  // itself, so its retardance from P is convention-free (the pi of the reflection reference
  // cancels over two reflections, docs/architecture.md): delta = |arg(rho_p^2 / rho_s^2)|
  // reduced to [0, pi]. Tolerance 1e-10.
  const Amplitudes a = fresnel_reference(1.0, Cx(1.2, 7.26), std::sin(45.0 * kDeg));
  const Cx rho_s = -a.rs;
  const Cx rho_p = a.rp;
  System s = load("metal_periscope.rtt.json");
  const CompiledSystem metal_cs = build(s);
  for (auto& child : s.root.children) {
    if (auto* e = std::get_if<Element>(&child.value);
        e && e->kind == rtt::model::ElementKind::Mirror) {
      e->material = std::nullopt;
    }
  }
  const CompiledSystem ideal_cs = build(s);
  const Vec3 start(0.0, 0.0, 0.0);
  const Vec3 k(0.0, 0.0, 1.0);
  const RayBatch m = trace_one(metal_cs, path_named(metal_cs, "periscope"), start, k);
  const RayBatch ideal = trace_one(ideal_cs, path_named(ideal_cs, "periscope"), start, k);
  REQUIRE(m.status()[0] == RayStatus::Alive);
  REQUIRE((direction(m, 0) - k).norm() <= 1e-12);
  const Vec3 s_hat(1.0, 0.0, 0.0);
  const std::vector<CVec3> inputs{CVec3(1.0, 0.0, 0.0), CVec3(0.0, 1.0, 0.0),
                                  CVec3(CVec3(Cx(1.0, 0.0), Cx(0.0, 1.0), 0.0) / std::sqrt(2.0))};
  for (const CVec3& e_in : inputs) {
    const CVec3 e_ideal = out(ideal, 0, e_in);
    const CVec3 s_part = s_hat.cast<Cx>() * s_hat.cast<Cx>().dot(e_ideal);
    const CVec3 expected = rho_s * rho_s * s_part + rho_p * rho_p * (e_ideal - s_part);
    REQUIRE((out(m, 0, e_in) - expected).norm() <= 1e-10);
  }
  const double expected_delta = std::abs(std::arg((rho_p * rho_p) / (rho_s * rho_s)));
  const auto ret = rtt::polar::retardance(m.prt_matrix(0), k);
  REQUIRE(std::abs(ret.value - expected_delta) <= 1e-10);
  REQUIRE(std::abs(power(m, 0, Vec3(1.0, 0.0, 0.0)) - std::norm(a.rs * a.rs)) <= 1e-10);
  REQUIRE(std::abs(power(m, 0, Vec3(0.0, 1.0, 0.0)) - std::norm(a.rp * a.rp)) <= 1e-10);
}

TEST_CASE("M3 Fresnel rhomb: two total internal reflections act as a quarter-wave plate", "[m3]") {
  // tests/reference/m3/fresnel_rhomb.rtt.json: plate of N-BK7 with four faces (toggle rule on an
  // explicit path: every face separates glass and vacuum), normal entrance and exit, two total
  // internal reflections on parallel faces at theta_R = 55.218468500096165 deg. theta_R solves
  // cos(2 arg(r_p / r_s)) = 0 with Byrnes, Eq. (6) and Appendix D for n = 1.5168 (offline
  // bisection, script in the PR), i.e. 45 deg (mod 90 deg) per reflection; the first check
  // below verifies theta_R itself (cos(delta_total) is linear in an error of theta_R). Total:
  // quarter wave, so the input (x + y)/sqrt(2) leaves circular with S0 = (1 - R)^2 (two normal
  // faces), and the path, which maps +z onto +z, has the retardance pi/2 (convention-free over
  // two reflections). The expected Stokes values and the retardance are computed from theta_R of
  // the file (S2 = S0 cos(delta_total), retardance |delta_total| reduced to [0, pi]), not from
  // the ideal 45 deg. Tolerance 1e-12.
  const double theta_r = 55.218468500096165 * kDeg;
  const Amplitudes a = fresnel_reference(kGlass, 1.0, kGlass * std::sin(theta_r));
  const double delta_total = 2.0 * std::arg(a.rp / a.rs);
  REQUIRE(std::abs(std::cos(delta_total)) <= 1e-14);  // theta_R gives a quarter wave
  const CompiledSystem cs = build(load("fresnel_rhomb.rtt.json"));
  const Vec3 k(0.0, 0.0, 1.0);
  const RayBatch rays = trace_one(cs, path_named(cs, "rhomb"), Vec3(0.0, 0.0, -1.0), k);
  REQUIRE(rays.status()[0] == RayStatus::Alive);
  REQUIRE((direction(rays, 0) - k).norm() <= 1e-12);
  const double r = (kGlass - 1.0) / (kGlass + 1.0);
  const double through = (1.0 - r * r) * (1.0 - r * r);
  REQUIRE(std::abs(rays.weight()[0] - through) <= 1e-12);
  const auto st = rtt::polar::stokes<double>(out(rays, 0, Vec3(1.0, 1.0, 0.0).normalized()),
                                             Vec3(1.0, 0.0, 0.0), direction(rays, 0));
  REQUIRE(std::abs(st.s0 - through) <= 1e-12);
  REQUIRE(std::abs(st.s1) <= 1e-12);
  REQUIRE(std::abs(st.s2 - st.s0 * std::cos(delta_total)) <= 1e-12);
  REQUIRE(std::abs(std::abs(st.s3) - st.s0) <= 1e-12);
  // Convention check (as for the single total internal reflection): s = x at both faces and the
  // output basis (x, y), so S3 = -S0 sin(delta_total) = +S0 here.
  REQUIRE(std::abs(st.s3 + st.s0 * std::sin(delta_total)) <= 1e-12);
  const auto ret = rtt::polar::retardance(rays.prt_matrix(0), k);
  REQUIRE(std::abs(ret.value - std::abs(std::remainder(delta_total, 2.0 * kPi))) <= 1e-12);
}
