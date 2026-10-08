// Refraction into and out of a uniaxial medium (#130, ADR 0026, points 4 and 5). Source:
// W.-S. T. Lam, "Anisotropic Ray Trace", dissertation, University of Arizona (docs/quellen.md):
// Eq. (2.11) (tangential phase matching), Eq. (2.39) (n_e(theta)), Eq. (2.41) (K-surface),
// p. 107 (S normal to the K-surface, E_e perpendicular to S_e), p. 104 (D_o, D_e), Eqs. (3.23),
// (3.24), (3.27) (P at isotropic/anisotropic intercepts), Section 4.2 (KTP example).

#include <Eigen/Geometry>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <numbers>
#include <optional>

#include "rtt/math/types.hpp"
#include "rtt/polar/birefringence.hpp"

using rtt::math::CMat3;
using rtt::math::Complex;
using rtt::math::Vec3;
using rtt::polar::eigen_polarization;
using rtt::polar::isotropic_from_tangential;
using rtt::polar::Mode;
using rtt::polar::ModeSolution;
using rtt::polar::prt_crystal_entry;
using rtt::polar::prt_crystal_exit;
using rtt::polar::uniaxial_mode;

namespace {

constexpr double kCalciteO = 1.6584;  // CONST indices of the ADR 0026 acceptance case
constexpr double kCalciteE = 1.4864;
const Vec3 kNormal = Vec3::UnitZ();  // surface normal into the crystal

/// Tangential part n t_par of a ray from a medium of index n at angle alpha in the y-z plane.
Vec3 tangential(double n, double alpha) {
  return {0.0, n * std::sin(alpha), 0.0};
}

double distance(const Vec3& a, const Vec3& b) {
  return (a - b).norm();
}

/// Lam, Eq. (2.41), e-mode factor, written in the principal frame of the crystal (axis = third
/// coordinate), as a function of k in global coordinates: F(k) = k1^2/nE^2 + k2^2/nE^2 +
/// k3^2/nO^2 - 1. Independent of the vector form used by uniaxial_mode().
double k_surface(const Vec3& k, const Vec3& axis, double n_o, double n_e) {
  const Vec3 a = axis.normalized();
  const Vec3 u = a.unitOrthogonal();
  const Vec3 v = a.cross(u);
  const double k1 = k.dot(u);
  const double k2 = k.dot(v);
  const double k3 = k.dot(a);
  return k1 * k1 / (n_e * n_e) + k2 * k2 / (n_e * n_e) + k3 * k3 / (n_o * n_o) - 1.0;
}

}  // namespace

TEST_CASE("ordinary mode: Snell with n_O, S = k, E perpendicular to axis and k", "[biref]") {
  const Vec3 axis = Vec3(0.3, 0.4, 0.866).normalized();
  for (const double alpha : {0.0, 0.3, 0.6, 1.2}) {
    INFO("alpha = " << alpha);
    const std::optional<ModeSolution<double>> m =
        uniaxial_mode(Mode::Ordinary, tangential(1.0, alpha), kNormal, kCalciteO, kCalciteE, axis);
    REQUIRE(m.has_value());
    const double s = std::sin(alpha) / kCalciteO;
    const Vec3 expected(0.0, s, std::sqrt(1.0 - s * s));  // Snell (de Greve, Eq. (22))
    REQUIRE(m->n == kCalciteO);
    REQUIRE(distance(m->k, expected) <= 1e-12);
    REQUIRE(distance(m->s, m->k) <= 1e-15);
    REQUIRE(std::abs(m->e.norm() - 1.0) <= 1e-15);
    REQUIRE(std::abs(m->e.dot(axis)) <= 1e-15);
    REQUIRE(std::abs(m->e.dot(m->k)) <= 1e-15);
  }
}

TEST_CASE("extraordinary mode, axis along the normal: closed form from (2.41)", "[biref]") {
  // k = (0, sin a, b) with sin^2 a / nE^2 + b^2 / nO^2 = 1 (axis z); S normal to the K-surface.
  for (const double alpha : {0.0, 0.4, 0.9}) {
    INFO("alpha = " << alpha);
    const auto m = uniaxial_mode(Mode::Extraordinary, tangential(1.0, alpha), kNormal, kCalciteO,
                                 kCalciteE, kNormal);
    REQUIRE(m.has_value());
    const double t = std::sin(alpha);
    const double b = kCalciteO * std::sqrt(1.0 - t * t / (kCalciteE * kCalciteE));
    const double n = std::sqrt(t * t + b * b);
    REQUIRE(std::abs(m->n - n) <= 1e-12);
    REQUIRE(distance(m->k, Vec3(0.0, t, b) / n) <= 1e-12);
    const Vec3 s = Vec3(0.0, t / (kCalciteE * kCalciteE), b / (kCalciteO * kCalciteO)).normalized();
    REQUIRE(distance(m->s, s) <= 1e-12);
  }
  // Normal incidence along the axis: n = n_O, no walk-off.
  const auto m = uniaxial_mode(Mode::Extraordinary, Vec3(Vec3::Zero()), kNormal, kCalciteO,
                               kCalciteE, kNormal);
  REQUIRE(m.has_value());
  REQUIRE(std::abs(m->n - kCalciteO) <= 1e-15);
  REQUIRE(distance(m->s, kNormal) <= 1e-15);
}

TEST_CASE("extraordinary mode, axis in the surface perpendicular to the plane of incidence",
          "[biref]") {
  // k is always perpendicular to the axis (x): n_e = n_E, Snell with n_E, no walk-off.
  for (const double alpha : {0.0, 0.5, 1.0}) {
    INFO("alpha = " << alpha);
    const auto m = uniaxial_mode(Mode::Extraordinary, tangential(1.0, alpha), kNormal, kCalciteO,
                                 kCalciteE, Vec3(Vec3::UnitX()));
    REQUIRE(m.has_value());
    const double s = std::sin(alpha) / kCalciteE;
    REQUIRE(std::abs(m->n - kCalciteE) <= 1e-12);
    REQUIRE(distance(m->k, Vec3(0.0, s, std::sqrt(1.0 - s * s))) <= 1e-12);
    REQUIRE(distance(m->s, m->k) <= 1e-12);
  }
}

TEST_CASE("calcite, axis at 45 degree, normal incidence: n_e(45) and walk-off", "[biref]") {
  const Vec3 axis = Vec3(0.0, 1.0, 1.0).normalized();
  const auto m =
      uniaxial_mode(Mode::Extraordinary, Vec3(Vec3::Zero()), kNormal, kCalciteO, kCalciteE, axis);
  REQUIRE(m.has_value());
  // Lam, Eq. (2.39) at theta = 45 degree.
  const double n45 =
      std::sqrt(2.0 / (1.0 / (kCalciteO * kCalciteO) + 1.0 / (kCalciteE * kCalciteE)));
  REQUIRE(std::abs(m->n - n45) <= 1e-12);
  REQUIRE(distance(m->k, kNormal) <= 1e-15);
  // Walk-off from S normal to the K-surface (p. 107): tan(theta_S) = (nO^2 / nE^2) tan(theta_k),
  // so at 45 degree tan(rho) = (nO^2 - nE^2) / (nO^2 + nE^2); for calcite (nE < nO) S turns
  // away from the axis (S_y < 0 for an axis with a_y > 0).
  const double o2 = kCalciteO * kCalciteO;
  const double e2 = kCalciteE * kCalciteE;
  const double rho = std::atan((o2 - e2) / (o2 + e2));
  REQUIRE(std::abs(rho * 180.0 / std::numbers::pi - 6.2241) <= 1e-4);
  REQUIRE(distance(m->s, Vec3(0.0, -std::sin(rho), std::cos(rho))) <= 1e-12);
  // Cross-check by geometry: the K-surface is the polar curve r(theta) = n_e(theta) about the
  // axis; its normal deviates from the radius by rho with tan(rho) = -(1/r) dr/dtheta. Here n_e
  // from (2.39) differentiated numerically, independent of the closed form above.
  const auto ne = [&](double theta) {
    const double c = std::cos(theta);
    const double s = std::sin(theta);
    return 1.0 / std::sqrt(c * c / o2 + s * s / e2);
  };
  const double h = 1e-6;
  const double theta = std::numbers::pi / 4.0;
  const double tan_rho = -(ne(theta + h) - ne(theta - h)) / (2.0 * h) / ne(theta);
  REQUIRE(std::abs(tan_rho - std::tan(rho)) <= 1e-8);
}

TEST_CASE("extraordinary S is the normal of the K-surface (2.41), numerical gradient", "[biref]") {
  // Independent of the closed form: S from the numerical gradient of Lam's Eq. (2.41) in the
  // principal frame, for an oblique ray and a tilted axis.
  const Vec3 axis = Vec3(0.2, -0.5, 0.8).normalized();
  const Vec3 t = Vec3(0.1, 0.45, 0.0);  // n t_par from air, |t| < 1
  const auto m = uniaxial_mode(Mode::Extraordinary, t, kNormal, kCalciteO, kCalciteE, axis);
  REQUIRE(m.has_value());
  const Vec3 k = m->n * m->k;
  REQUIRE(std::abs(k_surface(k, axis, kCalciteO, kCalciteE)) <= 1e-14);  // k on the surface
  REQUIRE(distance(Vec3(k.x(), k.y(), 0.0), t) <= 1e-14);                // phase matching
  const double h = 1e-6;
  Vec3 gradient;
  for (int i = 0; i < 3; ++i) {
    const Vec3 d = h * Vec3::Unit(i);
    gradient(i) = (k_surface(k + d, axis, kCalciteO, kCalciteE) -
                   k_surface(k - d, axis, kCalciteO, kCalciteE)) /
                  (2.0 * h);
  }
  REQUIRE(distance(m->s, gradient.normalized()) <= 1e-9);
  REQUIRE(m->s.dot(kNormal) > 0.0);  // energy flows into the crystal
  REQUIRE(m->s.dot(m->k) > 0.0);
}

TEST_CASE("Lam Section 4.2 (KTP): slow mode as the e-mode of a uniaxial medium", "[biref]") {
  // KTP with principal axes along x, y, z; a 500 nm ray from air at 35 degree in the y-z plane.
  // For k in the y-z plane the slow mode sees the ellipse ky^2/nS^2 + kz^2/nM^2 = 1: the e-mode
  // of a uniaxial medium with axis z, n_O = n_M, n_E = n_S. Printed values (3 digits): n = 1.807
  // (p. 144), E_ts = (0, -0.958, 0.286) (p. 146), S_ts = (0, 0.286, 0.958) from the column of
  // P_ts (p. 150); p. 146 prints S_ts with the value of the fast mode by mistake.
  const double n_m = 1.797;
  const double n_s = 1.902;
  const double alpha = 35.0 * std::numbers::pi / 180.0;
  const auto m =
      uniaxial_mode(Mode::Extraordinary, tangential(1.0, alpha), kNormal, n_m, n_s, kNormal);
  REQUIRE(m.has_value());
  // Lam p. 144: n' = sqrt(nM^2 + (1 - nM^2/nS^2) sin^2(theta)), to rounding.
  const double sin2 = std::sin(alpha) * std::sin(alpha);
  REQUIRE(std::abs(m->n - std::sqrt(n_m * n_m + (1.0 - n_m * n_m / (n_s * n_s)) * sin2)) <= 1e-12);
  REQUIRE(std::abs(m->n - 1.807) <= 5e-4);
  REQUIRE(distance(m->s, Vec3(0.0, 0.286, 0.958)) <= 1e-3);
  REQUIRE(std::abs(std::abs(m->e.dot(Vec3(0.0, -0.958, 0.286))) - 1.0) <= 1e-3);
  // The fast mode is an o-mode with n = n_F = 1.786: S_tf = (0, 0.321, 0.947) (p. 145).
  const auto f =
      uniaxial_mode(Mode::Ordinary, tangential(1.0, alpha), kNormal, 1.786, n_s, kNormal);
  REQUIRE(f.has_value());
  REQUIRE(distance(f->s, Vec3(0.0, 0.321, 0.947)) <= 1e-3);
  REQUIRE(std::abs(std::abs(f->e.dot(Vec3::UnitX())) - 1.0) <= 1e-15);  // E_tf = (1, 0, 0)
}

TEST_CASE("entry and exit through parallel faces restore the direction", "[biref]") {
  // Phase matching at both faces (Lam, Eq. (2.11)): the exit uses n k of the mode.
  const Vec3 axis = Vec3(0.2, -0.5, 0.8).normalized();
  const Vec3 in = Vec3(0.1, 0.45, std::sqrt(1.0 - 0.01 - 0.2025));
  for (const Mode mode : {Mode::Ordinary, Mode::Extraordinary}) {
    const Vec3 t(in.x(), in.y(), 0.0);
    const auto m = uniaxial_mode(mode, t, kNormal, kCalciteO, kCalciteE, axis);
    REQUIRE(m.has_value());
    const Vec3 nk = m->n * m->k;
    const auto out = isotropic_from_tangential(Vec3(nk.x(), nk.y(), 0.0), kNormal, 1.0);
    REQUIRE(out.has_value());
    REQUIRE(distance(*out, in) <= 1e-12);
    // With a diffraction term g at the entry and -g at the exit (ADR 0025, point 2).
    const Vec3 g(0.05, -0.02, 0.0);
    const auto mg = uniaxial_mode(mode, Vec3(t + g), kNormal, kCalciteO, kCalciteE, axis);
    REQUIRE(mg.has_value());
    const Vec3 nkg = mg->n * mg->k;
    const auto outg =
        isotropic_from_tangential(Vec3(Vec3(nkg.x(), nkg.y(), 0.0) - g), kNormal, 1.0);
    REQUIRE(outg.has_value());
    REQUIRE(distance(*outg, in) <= 1e-12);
  }
}

TEST_CASE("no real solution: total reflection and grazing count as none", "[biref]") {
  const Vec3 axis = Vec3(0.0, 1.0, 1.0).normalized();
  // From a dense isotropic medium (n = 2) at a steep angle into calcite.
  const Vec3 t = tangential(2.0, 1.2);  // |t| = 1.864 > nO
  REQUIRE_FALSE(uniaxial_mode(Mode::Ordinary, t, kNormal, kCalciteO, kCalciteE, axis).has_value());
  REQUIRE_FALSE(
      uniaxial_mode(Mode::Extraordinary, t, kNormal, kCalciteO, kCalciteE, axis).has_value());
  // Exit into air beyond the critical angle.
  REQUIRE_FALSE(isotropic_from_tangential(Vec3(0.0, 1.2, 0.0), kNormal, 1.0).has_value());
  // Just inside the limits a solution exists (the counter-check of the cases above).
  REQUIRE(isotropic_from_tangential(Vec3(0.0, 0.999, 0.0), kNormal, 1.0).has_value());
  REQUIRE(uniaxial_mode(Mode::Ordinary, Vec3(0.0, 0.999 * kCalciteO, 0.0), kNormal, kCalciteO,
                        kCalciteE, axis)
              .has_value());
  REQUIRE(
      uniaxial_mode(Mode::Extraordinary, tangential(2.0, 0.5), kNormal, kCalciteO, kCalciteE, axis)
          .has_value());
  // Exactly grazing: no propagating direction (ADR 0025, point 7).
  REQUIRE_FALSE(isotropic_from_tangential(Vec3(0.0, 1.0, 0.0), kNormal, 1.0).has_value());
  REQUIRE_FALSE(
      uniaxial_mode(Mode::Ordinary, Vec3(0.0, kCalciteO, 0.0), kNormal, kCalciteO, kCalciteE, axis)
          .has_value());
}

TEST_CASE("eigenpolarizations: orthogonality, plane of k and axis, degenerate case", "[biref]") {
  const Vec3 axis = Vec3(0.2, -0.5, 0.8).normalized();
  const Vec3 t(0.1, 0.45, 0.0);
  const auto o = uniaxial_mode(Mode::Ordinary, t, kNormal, kCalciteO, kCalciteE, axis);
  const auto e = uniaxial_mode(Mode::Extraordinary, t, kNormal, kCalciteO, kCalciteE, axis);
  REQUIRE(o.has_value());
  REQUIRE(e.has_value());
  // E_e is perpendicular to S_e (Lam p. 107) and lies in the plane of k_e and the axis (p. 104).
  REQUIRE(std::abs(e->e.dot(e->s)) <= 1e-15);
  REQUIRE(std::abs(e->e.dot(e->k.cross(axis).normalized())) <= 1e-15);
  REQUIRE(std::abs(e->e.norm() - 1.0) <= 1e-15);
  // The sign of the axis does not change the solution (ADR 0026, point 2).
  const auto e_minus =
      uniaxial_mode(Mode::Extraordinary, t, kNormal, kCalciteO, kCalciteE, Vec3(-axis));
  REQUIRE(e_minus.has_value());
  REQUIRE(distance(e_minus->s, e->s) <= 1e-15);
  REQUIRE(distance(e_minus->e, e->e) <= 1e-15);
  // k along the axis: both modes coincide; E_o by the axis rule of prt.hpp (smallest |k . e|,
  // ties x, y, z): k = z gives E_o = z x x = y and E_e = S x E_o = z x y = -x.
  const auto d = uniaxial_mode(Mode::Extraordinary, Vec3(Vec3::Zero()), kNormal, kCalciteO,
                               kCalciteE, kNormal);
  REQUIRE(d.has_value());
  REQUIRE(distance(eigen_polarization(Mode::Ordinary, d->k, d->s, kNormal), Vec3::UnitY()) <=
          1e-15);
  REQUIRE(distance(d->e, -Vec3::UnitX()) <= 1e-15);
}

TEST_CASE("P of entry and exit: directions, power, sign of the axis", "[biref]") {
  const Vec3 axis = Vec3(0.2, -0.5, 0.8).normalized();
  const Vec3 in = Vec3(0.1, 0.45, std::sqrt(1.0 - 0.01 - 0.2025));
  const Vec3 t(in.x(), in.y(), 0.0);
  for (const Mode mode : {Mode::Ordinary, Mode::Extraordinary}) {
    const auto m = uniaxial_mode(mode, t, kNormal, kCalciteO, kCalciteE, axis);
    REQUIRE(m.has_value());
    const CMat3 p = prt_crystal_entry(in, *m);
    // Lam, Eq. (3.23): P S_in = S'.
    REQUIRE((p * in.cast<Complex>() - m->s.cast<Complex>()).norm() <= 1e-15);
    // Lossless up to the projection: |P E| <= 1 for every unit E perpendicular to S_in.
    for (const double phi : {0.0, 0.7, 1.4, 2.1, 2.8}) {
      const Vec3 u = in.cross(Vec3::UnitX()).normalized();
      const Vec3 v = in.cross(u);
      const Vec3 field = std::cos(phi) * u + std::sin(phi) * v;
      REQUIRE((p * field.cast<Complex>()).norm() <= 1.0 + 1e-15);
    }
    // The sign of the axis does not change P (ADR 0026, point 5).
    const auto m_minus = uniaxial_mode(mode, t, kNormal, kCalciteO, kCalciteE, Vec3(-axis));
    REQUIRE((prt_crystal_entry(in, *m_minus) - p).norm() <= 1e-15);
    // Exit into air (Lam, Eqs. (3.24), (3.27)): P S_m = S_out, the absent mode carries nothing.
    const Vec3 nk = m->n * m->k;
    const auto out = isotropic_from_tangential(Vec3(nk.x(), nk.y(), 0.0), kNormal, 1.0);
    REQUIRE(out.has_value());
    const CMat3 q = prt_crystal_exit(m->s, m->e, *out);
    REQUIRE((q * m->s.cast<Complex>() - out->cast<Complex>()).norm() <= 1e-15);
    const Vec3 absent = m->s.cross(m->e).normalized();
    REQUIRE((q * absent.cast<Complex>()).norm() <= 1e-15);
  }
}

TEST_CASE("power of both modes: exactly 1 at normal incidence, 1 +- |e_o . e_e| oblique",
          "[biref]") {
  const auto coupling = [](const Vec3& in, const Vec3& t, const Vec3& axis) {
    const auto o = uniaxial_mode(Mode::Ordinary, t, kNormal, kCalciteO, kCalciteE, axis);
    const auto e = uniaxial_mode(Mode::Extraordinary, t, kNormal, kCalciteO, kCalciteE, axis);
    REQUIRE(o.has_value());
    REQUIRE(e.has_value());
    // Power of a unit field E perpendicular to S_in in both modes: |P_o E|^2 + |P_e E|^2 =
    // E^T (u u^T + v v^T) E with u, v the unit projections of E_o, E_e perpendicular to S_in
    // (P = E_v e_v^T + S S_in^T, ADR 0026, point 5).
    const CMat3 g = prt_crystal_entry(in, *o).adjoint() * prt_crystal_entry(in, *o) +
                    prt_crystal_entry(in, *e).adjoint() * prt_crystal_entry(in, *e);
    return g;
  };
  // Normal incidence, axis at 45 degree: the projections of E_o and E_e are x and y, the sum is
  // 1 for every polarization.
  {
    const Vec3 axis = Vec3(0.0, 1.0, 1.0).normalized();
    const CMat3 g = coupling(kNormal, Vec3(Vec3::Zero()), axis);
    for (const Vec3& field :
         {Vec3(1.0, 0.0, 0.0), Vec3(0.0, 1.0, 0.0), Vec3(Vec3(1.0, 1.0, 0.0).normalized())}) {
      // E is perpendicular to S_in, so the S S_in^T terms do not contribute.
      const double power = (field.cast<Complex>().adjoint() * g * field.cast<Complex>())(0).real();
      REQUIRE(std::abs(power - 1.0) <= 1e-14);
    }
    // Each mode on its own (the M4 acceptance of ADR 0026): E_o = z x a ~ x, E_e ~ y.
    const auto o =
        uniaxial_mode(Mode::Ordinary, Vec3(Vec3::Zero()), kNormal, kCalciteO, kCalciteE, axis);
    const auto e =
        uniaxial_mode(Mode::Extraordinary, Vec3(Vec3::Zero()), kNormal, kCalciteO, kCalciteE, axis);
    REQUIRE(o.has_value());
    REQUIRE(e.has_value());
    const CMat3 po = prt_crystal_entry(kNormal, *o);
    const CMat3 pe = prt_crystal_entry(kNormal, *e);
    const Vec3 x = Vec3::UnitX();
    const Vec3 y = Vec3::UnitY();
    REQUIRE(std::abs((po * x.cast<Complex>()).norm() - 1.0) <= 1e-15);
    REQUIRE((po * y.cast<Complex>()).norm() <= 1e-15);
    REQUIRE(std::abs((pe * y.cast<Complex>()).norm() - 1.0) <= 1e-15);
    REQUIRE((pe * x.cast<Complex>()).norm() <= 1e-15);
  }
  // Oblique incidence (H2 of ADR 0026, point 5): u and v are unit vectors in the plane
  // perpendicular to S_in, but not orthogonal, so the transverse part u u^T + v v^T has the
  // eigenvalues 1 + |u . v| and 1 - |u . v|: some polarization gets more than all of its power.
  // M4 accepts this; the anisotropic Fresnel coefficients (follow-up) remove it. Measured for
  // this case (calcite, alpha = 1.2 rad, azimuth (0.6, 0.8), axis (1, 0, 0.3)):
  // |u . v| = 0.20584, so the summed power lies between 0.79416 and 1.20584.
  {
    const Vec3 axis = Vec3(1.0, 0.0, 0.3).normalized();
    const double alpha = 1.2;
    const Vec3 in(0.6 * std::sin(alpha), 0.8 * std::sin(alpha), std::cos(alpha));
    const Vec3 t(in.x(), in.y(), 0.0);
    const auto o = uniaxial_mode(Mode::Ordinary, t, kNormal, kCalciteO, kCalciteE, axis);
    const auto e = uniaxial_mode(Mode::Extraordinary, t, kNormal, kCalciteO, kCalciteE, axis);
    REQUIRE(o.has_value());
    REQUIRE(e.has_value());
    const auto project = [&](const Vec3& x) { return Vec3((x - x.dot(in) * in).normalized()); };
    const double overlap = std::abs(project(o->e).dot(project(e->e)));
    INFO("|u . v| = " << overlap);
    REQUIRE(std::abs(overlap - 0.20584) <= 1e-5);  // the documented value above
    const CMat3 g = coupling(in, t, axis);
    const Vec3 a = in.cross(Vec3::UnitX()).normalized();
    const Vec3 b = in.cross(a);
    Eigen::Matrix2d transverse;
    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < 2; ++j) {
        const Vec3& x = i == 0 ? a : b;
        const Vec3& y = j == 0 ? a : b;
        transverse(i, j) = (x.cast<Complex>().adjoint() * g * y.cast<Complex>())(0).real();
      }
    }
    // Eigenvalues of the symmetric 2x2 matrix: half the trace -+ root, with the root written
    // without cancellation: sqrt(((t00 - t11) / 2)^2 + t01 t10) (= sqrt(half^2 - det)).
    const double half = 0.5 * (transverse(0, 0) + transverse(1, 1));
    const double diff = 0.5 * (transverse(0, 0) - transverse(1, 1));
    const double root = std::sqrt(diff * diff + transverse(0, 1) * transverse(1, 0));
    REQUIRE(std::abs((half - root) - (1.0 - overlap)) <= 1e-14);
    REQUIRE(std::abs((half + root) - (1.0 + overlap)) <= 1e-14);
  }
}
