#include <Eigen/Geometry>
#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numbers>
#include <optional>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/geom/phase.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/math/types.hpp"
#include "rtt/model/model.hpp"
#include "rtt/trace/apply_event.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/sequential.hpp"
#include "rtt/trace/sources.hpp"

// Diffraction orders in the tracer (#127, ADR 0025). Reference values come from the grating
// equation in the sources, not from the implementation:
// - C. Palmer, Diffraction Grating Handbook, 7th ed., Newport (2014): Eq. (2-1)
//   m lambda = d (sin alpha + sin beta), Fig. 2-1 for the signs (with x towards the "+" side
//   t_x = -sin alpha for the incident and t'_x = sin beta for the diffracted ray, reflection and
//   transmission), Eq. (2-3) conical diffraction G m lambda = cos(eps) (sin alpha + sin beta),
//   Eq. (2-5) Littrow m lambda = 2 d sin alpha;
// - M. Mansuripur, Proc. SPIE 6620, 66200N (2007): Eq. (7b) (n2 sigma'_x = n1 sigma_x +
//   m lambda0 dF/dx, transmission), Eq. (8) (reflection, n1 on both sides), Eq. (5) (phase
//   2 pi m F of the order);
// - J. Diebel, Representing Attitude, Stanford (2006), Eqs. (183)-(187): rotation about an axis.
// docs/quellen.md. G in lines/mm, lambda0 in mm; F = phi / (2 pi) = G x for a linear grating.

using rtt::compile::CompiledSurface;
using rtt::geom::LinearGratingPhase;
using rtt::geom::RadialPhasePolynomial;
using rtt::math::Mat3;
using rtt::math::Vec3;
using rtt::model::EventKind;
using rtt::trace::EventMedia;
using rtt::trace::RayState;
using rtt::trace::RayStatus;

namespace {

constexpr double kPi = std::numbers::pi;
constexpr double kLambdaUm = 0.5876;              // um
constexpr double kLambdaMm = kLambdaUm / 1000.0;  // mm
/// 2^-11 mm = 0.48828125 um: with 1024 or 2048 lines/mm, m lambda0 G is exactly 0.5 or 1.
constexpr double kExactLambdaUm = 0.48828125;

bool near(const Vec3& a, const Vec3& b, double tol) {
  return (a - b).cwiseAbs().maxCoeff() <= tol;
}

/// max |(R^T R - I)_ij|, in scalar code (an Eigen expression of this form triggers a false
/// -Wnull-dereference in GCC at -O2).
double orthonormality_error(const Mat3& r) {
  double e = 0.0;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      double dot = 0.0;
      for (int k = 0; k < 3; ++k) dot += r(k, i) * r(k, j);
      e = std::max(e, std::abs(dot - (i == j ? 1.0 : 0.0)));
    }
  }
  return e;
}

/// Plane surface at the global origin (local = global) with a linear grating of G lines/mm,
/// grating vector (cos psi, sin psi) (ADR 0025, point 1).
CompiledSurface grating(double lines_per_mm, double psi = 0.0) {
  CompiledSurface s;
  s.shape = rtt::geom::Plane<double>{};
  s.phase_functions = {LinearGratingPhase<double>(lines_per_mm, psi)};
  return s;
}

/// Media of one event: real indices, wavelength in um.
EventMedia media(double n_before, double n_after, double wavelength_um = kLambdaUm) {
  EventMedia m;
  m.before = n_before;
  m.after = n_after;
  m.beyond = n_after;
  m.wavelength_um = wavelength_um;
  return m;
}

/// Ray with unit direction `d` aimed at the global origin from z = -10 mm.
RayState ray_along(const Vec3& d) {
  RayState r;
  r.dir = d;
  r.pos = -10.0 / d.z() * d;
  return r;
}

RayState step(
    const RayState& ray, const CompiledSurface& s, EventKind kind, int order, const EventMedia& m) {
  return rtt::trace::sequential_step(ray, s, 3, kind, order, m);
}

double deg(double d) {
  return d * kPi / 180.0;
}

}  // namespace

TEST_CASE("transmission grating in vacuum: grating equation, Palmer (2-1) (#127)",
          "[diffraction]") {
  // Grating vector +x (psi = 0), G = 300/mm, vacuum, rays in the x-z plane with t_x = -sin alpha
  // (Palmer, Fig. 2-1b): sin beta = m lambda G - sin alpha, t' = (sin beta, 0, cos beta). Every
  // order up to |m| = 2 propagates (|sin beta| <= 0.78). The arithmetic is a few operations on
  // numbers of order 1: 1e-12.
  const CompiledSurface g = grating(300.0);
  for (const double alpha_deg : {0.0, 10.0, -25.0}) {
    const double alpha = deg(alpha_deg);
    const Vec3 d(-std::sin(alpha), 0.0, std::cos(alpha));
    for (int m = -2; m <= 2; ++m) {
      INFO("alpha = " << alpha_deg << " deg, m = " << m);
      const RayState out = step(ray_along(d), g, EventKind::Transmit, m, media(1.0, 1.0));
      REQUIRE(out.status == RayStatus::Alive);
      REQUIRE(out.last_surface == 3);
      const double sin_beta = m * kLambdaMm * 300.0 - std::sin(alpha);
      REQUIRE(near(out.dir, Vec3(sin_beta, 0.0, std::sqrt(1.0 - sin_beta * sin_beta)), 1e-12));
    }
  }
}

TEST_CASE("grating on a refracting surface: n' sin beta = n sin alpha + m lambda G, (7b) (#127)",
          "[diffraction]") {
  // Mansuripur, Eq. (7b), with n1 = 1, n2 = 1.5 and dF/dx = G: n2 t'_x = n1 t_x + m lambda0 G,
  // t'_y = t_y n1 / n2 = 0, t'_z > 0 from |t'| = 1. Fresnel refraction at the interface (ADR 0025,
  // point 6: interface as for order 0, then the layer). 1e-12.
  const CompiledSurface g = grating(600.0);
  for (const double alpha_deg : {0.0, 20.0, -35.0}) {
    const Vec3 d(-std::sin(deg(alpha_deg)), 0.0, std::cos(deg(alpha_deg)));
    for (const int m : {-1, 1, 2}) {
      INFO("alpha = " << alpha_deg << " deg, m = " << m);
      const RayState out = step(ray_along(d), g, EventKind::Refract, m, media(1.0, 1.5));
      REQUIRE(out.status == RayStatus::Alive);
      const double tx = (d.x() + m * kLambdaMm * 600.0) / 1.5;
      REQUIRE(near(out.dir, Vec3(tx, 0.0, std::sqrt(1.0 - tx * tx)), 1e-12));
    }
  }
}

TEST_CASE("reflection grating, Mansuripur (8), and Littrow, Palmer (2-5) (#127)", "[diffraction]") {
  // Mansuripur, Eq. (8): sigma'_x = sigma_x + (m lambda0 / n1) dF/dx on the incidence side,
  // here n1 = 1, ideal mirror; t'_z < 0. Littrow (Palmer, Eq. (2-5)): m lambda = 2 d sin alpha
  // sends the order back along the incident ray, k_out = -k_in. 1e-12.
  CompiledSurface g = grating(1200.0);
  g.interaction = rtt::model::IdealMirror{};
  for (const int m : {-1, 1}) {
    const Vec3 d(-std::sin(deg(15.0)), 0.0, std::cos(deg(15.0)));
    const RayState out = step(ray_along(d), g, EventKind::Reflect, m, media(1.0, 1.0));
    REQUIRE(out.status == RayStatus::Alive);
    const double tx = d.x() + m * kLambdaMm * 1200.0;
    REQUIRE(near(out.dir, Vec3(tx, 0.0, -std::sqrt(1.0 - tx * tx)), 1e-12));
  }
  // Littrow for m = 1: sin alpha = lambda G / 2.
  const double sin_alpha = kLambdaMm * 1200.0 / 2.0;
  const Vec3 d(-sin_alpha, 0.0, std::sqrt(1.0 - sin_alpha * sin_alpha));
  const RayState littrow = step(ray_along(d), g, EventKind::Reflect, 1, media(1.0, 1.0));
  REQUIRE(littrow.status == RayStatus::Alive);
  REQUIRE(near(littrow.dir, -d, 1e-12));
}

TEST_CASE("conical diffraction: the order stays on the cone, Palmer (2-3) (#127)",
          "[diffraction]") {
  // Grating vector u = (cos psi, sin psi, 0) with psi = 30 deg, grooves along w = (-sin psi,
  // cos psi, 0); incident ray with a component along the grooves. In the frame (u, w, z) Palmer's
  // eps is the angle to the plane perpendicular to the grooves, sin eps = t . w, and the angles
  // in that plane are sin alpha = -(t . u) / cos eps, sin beta = (t' . u) / cos eps (Fig. 2-1).
  // Eq. (2-3): G m lambda = cos eps (sin alpha + sin beta); the order keeps t' . w = t . w (the
  // spectra lie on a cone). 1e-12.
  const double psi = deg(30.0);
  const Vec3 u(std::cos(psi), std::sin(psi), 0.0);
  const Vec3 w(-std::sin(psi), std::cos(psi), 0.0);
  const CompiledSurface g = grating(500.0, psi);
  const Vec3 d = Vec3(0.15, -0.3, 0.9).normalized();
  for (const int m : {-1, 1, 2}) {
    INFO("m = " << m);
    const RayState out = step(ray_along(d), g, EventKind::Transmit, m, media(1.0, 1.0));
    REQUIRE(out.status == RayStatus::Alive);
    REQUIRE(std::abs(out.dir.norm() - 1.0) <= 1e-14);
    REQUIRE(std::abs(out.dir.dot(w) - d.dot(w)) <= 1e-12);
    const double cos_eps = std::sqrt(1.0 - d.dot(w) * d.dot(w));
    const double sin_alpha = -d.dot(u) / cos_eps;
    const double sin_beta = out.dir.dot(u) / cos_eps;
    REQUIRE(std::abs(cos_eps * (sin_alpha + sin_beta) - 500.0 * m * kLambdaMm) <= 1e-12);
  }
}

TEST_CASE("evanescent orders: below, above and exactly on the limit |tau| = n' (#127)",
          "[diffraction]") {
  // ADR 0025, point 7: the order is evanescent for |n t_par + m lambda0 g_par / (2 pi)| >= n'.
  // lambda0 = 2^-11 mm and G = 1024/mm give m lambda0 G = 0.5 exactly per order (no rounding),
  // so the limit can be hit exactly at oblique incidence t_x = 0.5: tau = 0.5 + 0.5 = 1 = n'.
  const CompiledSurface g = grating(1024.0);
  const auto ray_tx = [](double tx) { return ray_along(Vec3(tx, 0.0, std::sqrt(1.0 - tx * tx))); };
  const EventMedia vac = media(1.0, 1.0, kExactLambdaUm);
  // Exactly on the limit (grazing): evanescent.
  const RayState in = ray_tx(0.5);
  const RayState on = step(in, g, EventKind::Transmit, 1, vac);
  REQUIRE(on.status == RayStatus::Evanescent);
  REQUIRE(on.last_surface == 3);
  REQUIRE(near(on.pos, Vec3::Zero(), 1e-12));  // at the hit point, like Tir
  REQUIRE(on.dir == in.dir);
  REQUIRE(on.weight == in.weight);
  REQUIRE(on.prt == in.prt);
  // Just inside: propagates, nearly grazing; just outside: evanescent.
  const RayState below = step(ray_tx(0.5 - 1e-6), g, EventKind::Transmit, 1, vac);
  REQUIRE(below.status == RayStatus::Alive);
  REQUIRE(below.dir.z() > 0.0);
  REQUIRE(std::abs(below.dir.x() - (1.0 - 1e-6)) <= 1e-12);
  REQUIRE(step(ray_tx(0.5 + 1e-6), g, EventKind::Transmit, 1, vac).status == RayStatus::Evanescent);
  // The opposite order is fine: tau = 0.
  const RayState back = step(in, g, EventKind::Transmit, -1, vac);
  REQUIRE(back.status == RayStatus::Alive);
  REQUIRE(near(back.dir, Vec3::UnitZ(), 1e-15));
  // Normal incidence with G = 2048/mm: m lambda0 G = 1 = n' for m = +-1, evanescent; m = 0
  // passes.
  const CompiledSurface g2 = grating(2048.0);
  const RayState normal = ray_along(Vec3::UnitZ());
  REQUIRE(step(normal, g2, EventKind::Transmit, 1, vac).status == RayStatus::Evanescent);
  REQUIRE(step(normal, g2, EventKind::Transmit, -1, vac).status == RayStatus::Evanescent);
  REQUIRE(step(normal, g2, EventKind::Transmit, 0, vac).status == RayStatus::Alive);
}

TEST_CASE("Tir of order 0 comes before an order that would propagate (#127)", "[diffraction]") {
  // ADR 0025, points 6 and 7: without a transmitted order 0 there is no field to diffract. Glass
  // (n = 1.5) to vacuum with t_x = 0.8: n t_x = 1.2 > 1 is Tir; the order m = -1 with
  // m lambda0 G = -0.5 would have tau = 0.7 < 1, real, but the ray ends with Tir.
  const CompiledSurface g = grating(1024.0);
  const RayState in = ray_along(Vec3(0.8, 0.0, 0.6));
  const RayState out = step(in, g, EventKind::Refract, -1, media(1.5, 1.0, kExactLambdaUm));
  REQUIRE(out.status == RayStatus::Tir);
}

TEST_CASE("order 0 at a phase surface is bitwise the surface without phase layer (#127)",
          "[diffraction]") {
  // ADR 0025, point 2: for m = 0 the tracer does not evaluate the phase.
  CompiledSurface plain;
  plain.shape = rtt::geom::Conic<double>(1.0 / 40.0, 0.0);
  CompiledSurface phased = plain;
  phased.phase_functions = {LinearGratingPhase<double>(300.0, 0.2),
                            RadialPhasePolynomial<double>(5.0, {12.0, -3.0})};
  const RayState in = ray_along(Vec3(0.1, -0.2, 1.0).normalized());
  for (const EventKind kind : {EventKind::Refract, EventKind::Reflect, EventKind::Transmit}) {
    const RayState a = step(in, plain, kind, 0, media(1.0, 1.5));
    const RayState b = step(in, phased, kind, 0, media(1.0, 1.5));
    REQUIRE(a.status == b.status);
    REQUIRE(a.pos == b.pos);
    REQUIRE(a.dir == b.dir);
    REQUIRE(a.opl == b.opl);
    REQUIRE(a.prt == b.prt);
    REQUIRE(a.weight == b.weight);
    REQUIRE(a.last_surface == b.last_surface);
  }
}

TEST_CASE("diffraction efficiency multiplies the weight only (#127)", "[diffraction]") {
  // ADR 0025, point 5: without the list every order has efficiency 1; with it, listed orders
  // their value and the others 0 (also order 0), status Alive; P is unchanged.
  CompiledSurface with = grating(300.0);
  with.diffraction_efficiency = std::vector<rtt::model::DiffractionEfficiency>{{1, 0.8}, {0, 0.25}};
  const CompiledSurface without = grating(300.0);
  const RayState in = ray_along(Vec3(0.1, 0.0, 1.0).normalized());
  for (const auto& [m, eta] : {std::pair{1, 0.8}, std::pair{0, 0.25}, std::pair{-1, 0.0}}) {
    INFO("m = " << m);
    const RayState a = step(in, with, EventKind::Refract, m, media(1.0, 1.5));
    const RayState b = step(in, without, EventKind::Refract, m, media(1.0, 1.5));
    REQUIRE(a.status == RayStatus::Alive);
    REQUIRE(b.status == RayStatus::Alive);
    REQUIRE(b.weight > 0.0);
    REQUIRE(a.weight == b.weight * eta);
    REQUIRE(a.prt == b.prt);
    REQUIRE(a.dir == b.dir);
  }
}

TEST_CASE("optical path of an order: m phi lambda0 / (2 pi) and the eikonal condition (#127)",
          "[diffraction]") {
  // ADR 0025, point 3: Delta OPL = m phi lambda0 / (2 pi) (Mansuripur, Eq. (5): phase 2 pi m F).
  // Eikonal condition: d(Delta OPL)/dx equals the jump n' t'_x - n t_x of the tangential
  // momentum (point 2), checked by a central difference with h = 1e-3 mm. Delta OPL is a
  // polynomial in x; its third derivative m lambda0 / (2 pi) * 24 c_2 x / R^4 is below 5e-5 mm^-2
  // for x <= 4 mm, so the truncation h^2 / 6 * 5e-5 < 1e-11; the rounding of the OPL (about
  // 10 mm, ~2e-15 each) over 2 h gives ~2e-12. Tolerance 1e-10.
  CompiledSurface s;
  s.shape = rtt::geom::Plane<double>{};
  s.phase_functions = {RadialPhasePolynomial<double>(10.0, {-500.0, 50.0})};
  const int m = 1;
  const auto delta_opl = [&](double x) {
    const RayState in = ray_along(Vec3::UnitZ());
    RayState shifted = in;
    shifted.pos.x() = x;
    const RayState zero = step(shifted, s, EventKind::Transmit, 0, media(1.0, 1.0));
    const RayState first = step(shifted, s, EventKind::Transmit, m, media(1.0, 1.0));
    return std::pair{first.opl - zero.opl, first.dir.x()};
  };
  for (const double x : {0.0, 1.5, 4.0}) {
    INFO("x = " << x);
    const double u = x * x / 100.0;
    const double phi = -500.0 * u + 50.0 * u * u;
    const auto [dopl, tx] = delta_opl(x);
    // Delta OPL is the difference of two OPLs of about 10 mm: rounding ~4e-15, tolerance 1e-14.
    REQUIRE(std::abs(dopl - m * phi * kLambdaMm / (2.0 * kPi)) <= 1e-14 + 1e-12 * std::abs(dopl));
    constexpr double h = 1e-3;
    const double slope = (delta_opl(x + h).first - delta_opl(x - h).first) / (2.0 * h);
    REQUIRE(std::abs(slope - tx) <= 1e-10);  // n = n' = 1, t_x = 0
  }
}

TEST_CASE("radial phase as a thin lens: f = -pi R^2 / (m lambda0 c1) (#127)", "[diffraction]") {
  // ADR 0025, point 3 (from point 2): phi = c1 rho^2 gives g_y = 2 c1 y / R^2, so a ray parallel
  // to the axis at height y leaves with t'_y = m lambda0 c1 y / (pi R^2) = -y / f exactly (the
  // relation is linear in y; only t'_z carries the nonparaxial part). c1 < 0 with m = +1:
  // converging; c1 > 0 diverging (t'_y has the sign of y). 1e-12 relative.
  const double radius = 10.0;
  const double f = 100.0;
  const double c1 = -kPi * radius * radius / (kLambdaMm * f);
  for (const double sign : {1.0, -1.0}) {
    CompiledSurface s;
    s.shape = rtt::geom::Plane<double>{};
    s.phase_functions = {RadialPhasePolynomial<double>(radius, {sign * c1})};
    for (const double y : {0.5, -2.0, 4.0}) {
      INFO("sign = " << sign << ", y = " << y);
      RayState in = ray_along(Vec3::UnitZ());
      in.pos.y() = y;
      const RayState out = step(in, s, EventKind::Transmit, 1, media(1.0, 1.0));
      REQUIRE(out.status == RayStatus::Alive);
      REQUIRE(std::abs(out.dir.y() - (-sign * y / f)) <= 1e-12 * std::abs(y / f));
      REQUIRE(out.dir.x() == 0.0);
      REQUIRE(out.dir.z() > 0.0);
    }
  }
}

TEST_CASE("polarization of an order: P = R(k_0 -> k_m) P_0 without rotation about k (#127)",
          "[diffraction]") {
  // ADR 0025, point 6. Diffractive lens at normal incidence, thin (Transmit at a Fresnel surface
  // is a dummy passage, P_0 = I): P = R(z -> k_m), orthonormal, P z = k_m; for k_m in the y-z
  // plane the axis is x, so x stays x (no rotation about k); on the axis k_m = z and P = I
  // exactly.
  CompiledSurface lens;
  lens.shape = rtt::geom::Plane<double>{};
  lens.phase_functions = {RadialPhasePolynomial<double>(10.0, {-5000.0})};
  for (const double y : {0.0, 1.0, -3.0}) {
    INFO("y = " << y);
    RayState in = ray_along(Vec3::UnitZ());
    in.pos.y() = y;
    const RayState out = step(in, lens, EventKind::Transmit, 1, media(1.0, 1.0));
    REQUIRE(out.status == RayStatus::Alive);
    const rtt::math::CMat3& p = out.prt;
    REQUIRE(p.imag().cwiseAbs().maxCoeff() == 0.0);
    const Mat3 r = p.real();
    REQUIRE(orthonormality_error(r) <= 1e-14);
    REQUIRE(near(r * Vec3::UnitZ(), out.dir, 1e-14));
    REQUIRE(near(r * Vec3::UnitX(), Vec3::UnitX(), 1e-14));
    if (y == 0.0) REQUIRE(p == rtt::math::CMat3::Identity());
  }
  // Refraction with Fresnel amplitudes: P_m = R(k_0 -> k_m) P_0 with P_0 that of order 0, and the
  // power fraction (weight) is that of order 0 (R keeps the transverse norm, ADR 0021). 1e-14.
  const CompiledSurface g = grating(300.0);
  const RayState in = ray_along(Vec3(0.2, 0.1, 1.0).normalized());
  const RayState zero = step(in, g, EventKind::Refract, 0, media(1.0, 1.5));
  const RayState first = step(in, g, EventKind::Refract, 1, media(1.0, 1.5));
  const rtt::math::CMat3 expected =
      rtt::trace::rotation_between(zero.dir, first.dir).cast<std::complex<double>>() * zero.prt;
  REQUIRE((first.prt - expected).cwiseAbs().maxCoeff() <= 1e-14);
  REQUIRE(std::abs(first.weight - zero.weight) <= 1e-14 * zero.weight);
}

TEST_CASE("rotation_between: R a = b, orthonormal, det 1, also close to antiparallel (#127)",
          "[diffraction]") {
  // Diebel (2006), Eqs. (183)-(187), transposed (active). Errors: the axis n = v / |v| carries
  // about eps / s (v = a x b is accurate to ~eps absolutely, s = |v| = sin alpha), and R is
  // linear in n with factors <= 2, so |R^T R - I| and |R a - b| stay below about 8 eps / s; the
  // tolerance is 64 eps / s, plus 64 eps for s near 1.
  const double eps = std::numeric_limits<double>::epsilon();
  REQUIRE(rtt::trace::rotation_between(Vec3::UnitZ(), Vec3::UnitZ()) == Mat3::Identity());
  std::vector<std::pair<Vec3, Vec3>> pairs = {
      {Vec3::UnitZ(), Vec3::UnitX()},
      {Vec3(0.3, -0.4, 0.8).normalized(), Vec3(-0.1, 0.9, 0.2).normalized()},
      {Vec3::UnitZ(), Vec3(1e-9, 0.0, 1.0).normalized()},  // nearly parallel
  };
  for (const double one_plus_c : {1e-3, 1e-6, 1e-9}) {
    // b with a . b = c = -1 + one_plus_c: b = (s, 0, c) with s = sqrt(1 - c^2).
    const double c = -1.0 + one_plus_c;
    pairs.emplace_back(Vec3::UnitZ(), Vec3(std::sqrt(1.0 - c * c), 0.0, c));
  }
  for (const auto& [a, b] : pairs) {
    INFO("a = " << a.transpose() << ", b = " << b.transpose());
    const Mat3 r = rtt::trace::rotation_between(a, b);
    const double s = a.cross(b).norm();
    const double tol = 64.0 * eps / s + 64.0 * eps;
    REQUIRE((r * a - b).cwiseAbs().maxCoeff() <= tol);
    REQUIRE(orthonormality_error(r) <= tol);
    REQUIRE(std::abs(r.determinant() - 1.0) <= tol);
    // The axis a x b stays fixed: no rotation about it.
    const Vec3 n = a.cross(b) / s;
    REQUIRE((r * n - n).cwiseAbs().maxCoeff() <= tol);
  }
}

TEST_CASE("an order without wavelength or phase layer gives EventImpossible (#127)",
          "[diffraction]") {
  // Guards of apply_event: lambda0 is needed for m != 0; validate() allows orders only at
  // surfaces with phase layers, a direct call without them is a modelling error.
  const CompiledSurface g = grating(300.0);
  const RayState in = ray_along(Vec3::UnitZ());
  REQUIRE(step(in, g, EventKind::Transmit, 1, media(1.0, 1.0, 0.0)).status ==
          RayStatus::EventImpossible);
  CompiledSurface bare;
  bare.shape = rtt::geom::Plane<double>{};
  REQUIRE(step(in, bare, EventKind::Transmit, 1, media(1.0, 1.0)).status ==
          RayStatus::EventImpossible);
  REQUIRE(step(in, bare, EventKind::Transmit, 0, media(1.0, 1.0)).status == RayStatus::Alive);
}

TEST_CASE("compile resolves the phase layers for the tracer (#127)", "[diffraction]") {
  // ADR 0025, point 1: the same numbers, parameters as values, orientation_deg in rad.
  rtt::model::System s;
  s.name = "phases";
  s.environment.medium = "VACUUM";
  s.wavelengths = {{kLambdaUm, 1.0, true}};
  s.aperture = {rtt::model::SystemApertureType::EntrancePupilDiameter, rtt::model::Param(4.0)};
  s.fields = {rtt::model::FieldType::AngleDeg, {{0.0, 0.0, 1.0}}};
  s.root.name = "root";
  rtt::model::Surface surface;
  surface.id = rtt::model::SurfaceId("G");
  surface.phases = {rtt::model::LinearGrating{rtt::model::Param(300.0), 30.0},
                    rtt::model::RadialPhase{rtt::model::Param(8.0),
                                            {rtt::model::Param(2.0), rtt::model::Param(-1.0)}}};
  s.root.children = {{rtt::model::Element{"G",
                                          rtt::model::ElementKind::ThinElement,
                                          rtt::model::Pose::along_z(0.0),
                                          std::nullopt,
                                          {surface}}}};
  s.paths = {{"main", true, {}}};
  const rtt::material::MaterialLibrary lib;
  const auto cs = rtt::compile::compile(s, lib);
  const auto& f = cs.surfaces()[0].phase_functions;
  REQUIRE(f.size() == 2);
  const auto* g = std::get_if<LinearGratingPhase<double>>(&f[0]);
  REQUIRE(g != nullptr);
  REQUIRE(g->lines_per_mm() == 300.0);
  REQUIRE(std::abs(g->orientation() - kPi / 6.0) <= 1e-15);
  const auto* r = std::get_if<RadialPhasePolynomial<double>>(&f[1]);
  REQUIRE(r != nullptr);
  REQUIRE(r->normalization_radius() == 8.0);
  REQUIRE(r->coefficients() == std::vector<double>{2.0, -1.0});
}

TEST_CASE("reference file m4/grating_transmission: orders -1, 0, +1 by Palmer (2-1) (#127)",
          "[diffraction]") {
  // tests/reference/m4/grating_transmission.rtt.json: thin grating, 300 lines/mm, grooves along
  // y, vacuum, lambda = 0.5876 um, screen 50 mm behind it. A ray at alpha = 10 deg in the x-z
  // plane (t_x = -sin alpha): sin beta = m lambda G - sin alpha (Palmer, Eq. (2-1), Fig. 2-1b).
  const rtt::material::MaterialLibrary lib;
  const auto cs = rtt::compile::compile(
      rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m4/grating_transmission.rtt.json"),
      lib);
  const double alpha = deg(10.0);
  const Vec3 d(-std::sin(alpha), 0.0, std::cos(alpha));
  for (const auto& [name, m] :
       {std::pair{"order -1", -1}, std::pair{"order 0", 0}, std::pair{"order +1", 1}}) {
    INFO(name);
    const auto path = cs.find_path(name);
    REQUIRE(path.has_value());
    rtt::trace::RayBatch rays(1);
    rays.pos_x()[0] = 0.0;
    rays.pos_y()[0] = 0.0;
    rays.pos_z()[0] = -1.0;
    rays.dir_x()[0] = d.x();
    rays.dir_y()[0] = d.y();
    rays.dir_z()[0] = d.z();
    [[maybe_unused]] const auto stats = rtt::trace::SequentialTracer().trace(cs, *path, rays);
    REQUIRE(rays.status()[0] == RayStatus::Alive);
    const double sin_beta = m * kLambdaMm * 300.0 - std::sin(alpha);
    REQUIRE(std::abs(rays.dir_x()[0] - sin_beta) <= 1e-12);
    REQUIRE(std::abs(rays.dir_z()[0] - std::sqrt(1.0 - sin_beta * sin_beta)) <= 1e-12);
  }
}

TEST_CASE("ray aiming through a phase surface of order 0 before the stop (#127)", "[diffraction]") {
  // ADR 0025, point 4 (start rays): first_order accepts phase surfaces at order 0, and stop_hit
  // traces them like the tracer, also with efficiencies (which only scale the weight). The
  // aimed rays are bitwise those of the same system without phase layer and efficiency (order 0
  // takes the code path of a plain surface, ADR 0025, point 2).
  const auto system = [](bool phased) {
    rtt::model::System s;
    s.name = "aim";
    s.environment.medium = "VACUUM";
    s.wavelengths = {{kLambdaUm, 1.0, true}};
    s.aperture = {rtt::model::SystemApertureType::EntrancePupilDiameter, rtt::model::Param(4.0)};
    s.fields = {rtt::model::FieldType::AngleDeg, {{0.0, 0.0, 1.0}, {0.0, 3.0, 1.0}}};
    s.root.name = "root";
    rtt::model::Surface g;
    g.id = rtt::model::SurfaceId("G");
    if (phased) {
      g.phases = {rtt::model::LinearGrating{rtt::model::Param(300.0), 90.0}};
      g.diffraction_efficiency = std::vector<rtt::model::DiffractionEfficiency>{{0, 0.5}};
    }
    rtt::model::Surface stop;
    stop.id = rtt::model::SurfaceId("STO");
    stop.aperture = rtt::model::CircularAperture{3.0, 0.0};
    rtt::model::Surface image;
    image.id = rtt::model::SurfaceId("IMG");
    s.root.children = {{rtt::model::Element{"G",
                                            rtt::model::ElementKind::ThinElement,
                                            rtt::model::Pose::along_z(-20.0),
                                            std::nullopt,
                                            {g}}},
                       {rtt::model::Element{"stop",
                                            rtt::model::ElementKind::Stop,
                                            rtt::model::Pose::along_z(0.0),
                                            std::nullopt,
                                            {stop}}},
                       {rtt::model::Element{"image",
                                            rtt::model::ElementKind::Detector,
                                            rtt::model::Pose::along_z(30.0),
                                            std::nullopt,
                                            {image}}}};
    s.paths = {{"main", true, {}}};
    const rtt::material::MaterialLibrary lib;
    return rtt::compile::compile(s, lib);
  };
  const auto a = system(true);
  const auto b = system(false);
  const std::vector<std::uint16_t> fields{0, 1};
  const auto ra = rtt::trace::make_rays(a, rtt::compile::PathId{0}, fields, 0,
                                        rtt::trace::HexapolarPupil{2}, rtt::trace::Aiming::Real);
  const auto rb = rtt::trace::make_rays(b, rtt::compile::PathId{0}, fields, 0,
                                        rtt::trace::HexapolarPupil{2}, rtt::trace::Aiming::Real);
  REQUIRE(ra.size() == rb.size());
  for (std::size_t i = 0; i < ra.size(); ++i) {
    REQUIRE(ra.pos_x()[i] == rb.pos_x()[i]);
    REQUIRE(ra.pos_y()[i] == rb.pos_y()[i]);
    REQUIRE(ra.pos_z()[i] == rb.pos_z()[i]);
    REQUIRE(ra.dir_x()[i] == rb.dir_x()[i]);
    REQUIRE(ra.dir_y()[i] == rb.dir_y()[i]);
    REQUIRE(ra.dir_z()[i] == rb.dir_z()[i]);
  }
}

TEST_CASE("ray aiming through an absorbing medium: directions from Re(n) only (#127)",
          "[diffraction]") {
  // stop_hit uses the complex indices of the tracer (event_media.hpp) since #127. Direction and
  // hit point depend on Re(n) only, so the aimed rays through an absorbing glass before the stop
  // are bitwise those through the same glass without absorption.
  const auto system = [](const std::string& glass) {
    rtt::model::System s;
    s.name = "aim";
    s.environment.medium = "VACUUM";
    s.wavelengths = {{kLambdaUm, 1.0, true}};
    s.aperture = {rtt::model::SystemApertureType::EntrancePupilDiameter, rtt::model::Param(4.0)};
    s.fields = {rtt::model::FieldType::AngleDeg, {{0.0, 0.0, 1.0}, {0.0, 3.0, 1.0}}};
    s.root.name = "root";
    rtt::model::Surface s1;
    s1.id = rtt::model::SurfaceId("L.S1");
    s1.shape.base = rtt::model::Conic{rtt::model::Param(40.0), rtt::model::Param(0.0)};
    rtt::model::Surface s2;
    s2.id = rtt::model::SurfaceId("L.S2");
    s2.pose = rtt::model::Pose::along_z(4.0);
    rtt::model::Surface stop;
    stop.id = rtt::model::SurfaceId("STO");
    stop.aperture = rtt::model::CircularAperture{3.0, 0.0};
    rtt::model::Surface image;
    image.id = rtt::model::SurfaceId("IMG");
    s.root.children = {
        {rtt::model::Element{
            "L", rtt::model::ElementKind::Lens, rtt::model::Pose::along_z(-20.0), glass, {s1, s2}}},
        {rtt::model::Element{"stop",
                             rtt::model::ElementKind::Stop,
                             rtt::model::Pose::along_z(0.0),
                             std::nullopt,
                             {stop}}},
        {rtt::model::Element{"image",
                             rtt::model::ElementKind::Detector,
                             rtt::model::Pose::along_z(60.0),
                             std::nullopt,
                             {image}}}};
    s.paths = {{"main", true, {}}};
    const rtt::material::MaterialLibrary lib;
    return rtt::compile::compile(s, lib);
  };
  const auto a = system("CONST:1.5,1e-3");
  const auto b = system("CONST:1.5");
  const std::vector<std::uint16_t> fields{0, 1};
  const auto ra = rtt::trace::make_rays(a, rtt::compile::PathId{0}, fields, 0,
                                        rtt::trace::HexapolarPupil{2}, rtt::trace::Aiming::Real);
  const auto rb = rtt::trace::make_rays(b, rtt::compile::PathId{0}, fields, 0,
                                        rtt::trace::HexapolarPupil{2}, rtt::trace::Aiming::Real);
  REQUIRE(ra.size() == rb.size());
  for (std::size_t i = 0; i < ra.size(); ++i) {
    REQUIRE(ra.pos_x()[i] == rb.pos_x()[i]);
    REQUIRE(ra.pos_y()[i] == rb.pos_y()[i]);
    REQUIRE(ra.dir_x()[i] == rb.dir_x()[i]);
    REQUIRE(ra.dir_y()[i] == rb.dir_y()[i]);
    REQUIRE(ra.dir_z()[i] == rb.dir_z()[i]);
  }
}
