// GCC (13 in CI, 16 locally) reports a false -Wnull-dereference inside Eigen's SSE code
// (emmintrin.h, BinaryFunctors.h) for the fixed-size SVD and eigen solvers at -O2/-O3 in Release.
// The warning is attributed to the Eigen headers, so the suppression has to cover the includes;
// it is limited to this translation unit and to GCC (cf. the GCC 13 workaround in agf.cpp, #24).
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wnull-dereference"
#endif

// 3D polarization ray tracing matrix (#57). Source: W.-S. T. Lam, "Anisotropic ray trace",
// dissertation, University of Arizona, Eqs. (3.1)-(3.9), (4.3)-(4.6), Sec. 4.5.2
// (docs/quellen.md); conventions in rtt/polar/prt.hpp and docs/architecture.md.

#include <Eigen/Geometry>
#include <algorithm>
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <random>
#include <utility>

#include "rtt/polar/fresnel.hpp"
#include "rtt/polar/prt.hpp"
#include "rtt/polar/prt_analysis.hpp"

using Cx = std::complex<double>;
using rtt::math::CMat3;
using rtt::math::Mat3;
using rtt::math::Vec3;
using CVec3 = Eigen::Vector3cd;

namespace {

constexpr double kPi = std::numbers::pi;

/// Deterministic uniform numbers in [0, 1) from a fixed seed: (x >> 11) * 2^-53 as for the
/// random pupil of rtt-trace.
class Uniform {
 public:
  explicit Uniform(std::uint64_t seed) : engine_(seed) {}
  double operator()() { return static_cast<double>(engine_() >> 11) * 0x1.0p-53; }

 private:
  std::mt19937_64 engine_;
};

Vec3 random_unit(Uniform& u) {
  const double z = 2.0 * u() - 1.0;
  const double phi = 2.0 * kPi * u();
  const double r = std::sqrt(1.0 - z * z);
  return {r * std::cos(phi), r * std::sin(phi), z};
}

/// Random complex field transverse to k, |E| = 1.
CVec3 random_transverse(Uniform& u, const Vec3& k) {
  CVec3 e(Cx(u() - 0.5, u() - 0.5), Cx(u() - 0.5, u() - 0.5), Cx(u() - 0.5, u() - 0.5));
  const CVec3 kc = k.cast<Cx>();
  e -= kc * kc.dot(e);
  return e.normalized();
}

/// Refracted direction, B. de Greve, "Reflections and Refractions in Ray Tracing" (2006),
/// Eqs. (22), (23), (28) (docs/quellen.md): t = mu d + (mu cos_i - cos_t) N with N oriented into
/// the incident medium. Written here independently of rtt-trace.
Vec3 refract(const Vec3& d, Vec3 n, double n1, double n2) {
  if (n.dot(d) > 0.0) n = -n;
  const double mu = n1 / n2;
  const double cos_i = -n.dot(d);
  const double cos_t = std::sqrt(1.0 - mu * mu * (1.0 - cos_i * cos_i));
  return mu * d + (mu * cos_i - cos_t) * n;
}

Vec3 reflect(const Vec3& d, const Vec3& n) {
  return d - 2.0 * d.dot(n) * n;
}

/// Rotation matrix from a random unit quaternion.
Mat3 random_rotation(Uniform& u) {
  Eigen::Quaterniond q(u() - 0.5, u() - 0.5, u() - 0.5, u() - 0.5);
  return q.normalized().toRotationMatrix();
}

double max_abs(const CMat3& m) {
  return m.cwiseAbs().maxCoeff();
}

}  // namespace

TEST_CASE("PRT basis: right-handed s, p, k and the normal-incidence axis rule", "[prt]") {
  // Lam, Eq. (3.3): s = k x N / |k x N|, p = k x s. Normal incidence (decided for #57): s = k x e
  // with the global axis e of smallest |k . e|, ties in the order x, y, z.
  Uniform u(1);
  for (int i = 0; i < 100; ++i) {
    const Vec3 k = random_unit(u);
    const Vec3 n = random_unit(u);
    const auto b = rtt::polar::prt_basis(k, n);
    REQUIRE(std::abs(b.s.norm() - 1.0) <= 1e-15);
    REQUIRE(std::abs(b.s.dot(k)) <= 1e-15);
    REQUIRE(std::abs(b.s.dot(n)) <= 1e-14);
    REQUIRE((b.s.cross(b.p) - k).norm() <= 1e-15);
  }
  const auto z = rtt::polar::prt_basis(Vec3(0.0, 0.0, 1.0), Vec3(0.0, 0.0, -1.0));
  // Tolerance 1e-15 (unit vectors from cross products of coordinate axes).
  REQUIRE((z.s - Vec3(0.0, 1.0, 0.0)).norm() <= 1e-15);  // e = x (tie x, y): z x x = y
  REQUIRE((z.p - Vec3(-1.0, 0.0, 0.0)).norm() <= 1e-15);
  const auto x = rtt::polar::prt_basis(Vec3(1.0, 0.0, 0.0), Vec3(1.0, 0.0, 0.0));
  REQUIRE((x.s - Vec3(0.0, 0.0, 1.0)).norm() <= 1e-15);  // e = y (tie y, z): x x y = z
  REQUIRE((x.p - Vec3(0.0, -1.0, 0.0)).norm() <= 1e-15);
}

TEST_CASE("P maps k_in to k_out and transverse fields to transverse fields", "[prt]") {
  // Lam, Eqs. (3.1), (3.2): E' = P E, k' = P k. Arbitrary complex a_s, a_p; refraction and
  // reflection. Tolerance 1e-14 (rounding of unit vectors).
  Uniform u(2);
  for (int i = 0; i < 200; ++i) {
    const Vec3 k = random_unit(u);
    const Vec3 n = random_unit(u);
    if (std::abs(k.dot(n)) < 0.05) continue;  // avoid grazing incidence
    const Cx a_s(u() - 0.5, u() - 0.5);
    const Cx a_p(u() - 0.5, u() - 0.5);
    for (const bool reflection : {false, true}) {
      const Vec3 k_out = reflection ? reflect(k, n) : refract(k, n, 1.0, 1.6);
      const CMat3 p = rtt::polar::prt_matrix(k, k_out, n, a_s, a_p);
      REQUIRE((p * k.cast<Cx>() - k_out.cast<Cx>()).norm() <= 1e-14);
      const CVec3 e = random_transverse(u, k);
      REQUIRE(std::abs(k_out.cast<Cx>().dot(p * e)) <= 1e-14);
      const Mat3 q = rtt::polar::geometric_transform(k, k_out, n, reflection);
      REQUIRE((q * k - k_out).norm() <= 1e-14);
      REQUIRE((q.transpose() * q - Mat3::Identity()).cwiseAbs().maxCoeff() <= 1e-14);
      REQUIRE(std::abs(q.determinant() - (reflection ? -1.0 : 1.0)) <= 1e-14);
    }
  }
}

TEST_CASE("energy at a lossless interface: |P_r E|^2 + c |P_t E|^2 = |E|^2", "[prt]") {
  // Fresnel amplitudes of #56 (Byrnes, Eqs. (6), (21)-(23)): R + T = 1 for s and p with the same
  // factor c = n_t cos_t / (n_i cos_i) = Re(q_t) / Re(q_i) for lossless media. The out bases are
  // orthonormal, so |P E|^2 = |a_s E_s|^2 + |a_p E_p|^2 for E transverse to k_in. Both
  // directions (also total internal reflection, c = 0). Tolerance 1e-12.
  Uniform u(3);
  for (const auto& [n_i, n_t] : {std::pair{1.0, 1.5168}, std::pair{1.5168, 1.0}}) {
    for (int i = 0; i < 200; ++i) {
      const Vec3 k = random_unit(u);
      const Vec3 n = random_unit(u);
      const double cos_i = std::abs(k.dot(n));
      if (cos_i < 0.05) continue;
      const double xi = n_i * std::sqrt(1.0 - cos_i * cos_i);
      const auto a = rtt::polar::fresnel(Cx(n_i), Cx(n_t), xi);
      const Cx q_i = rtt::polar::normal_component(Cx(n_i), xi);
      const Cx q_t = rtt::polar::normal_component(Cx(n_t), xi);
      const double c = q_t.real() / q_i.real();
      const CMat3 p_r = rtt::polar::prt_matrix(k, reflect(k, n), n, a.rs, a.rp);
      const CVec3 e = random_transverse(u, k);
      double power = (p_r * e).squaredNorm();
      if (c > 0.0) {
        const CMat3 p_t = rtt::polar::prt_matrix(k, refract(k, n, n_i, n_t), n, a.ts, a.tp);
        power += c * (p_t * e).squaredNorm();
      }
      INFO("n_i " << n_i << ", cos_i " << cos_i);
      REQUIRE(std::abs(power - 1.0) <= 1e-12);
    }
  }
}

TEST_CASE("P does not depend on the orientation of N", "[prt]") {
  // Reversing N reverses s and both p together (prt.hpp); P and Q are unchanged.
  Uniform u(4);
  for (int i = 0; i < 100; ++i) {
    const Vec3 k = random_unit(u);
    const Vec3 n = random_unit(u);
    const Vec3 k_out = reflect(k, n);
    const Cx a_s(u(), u());
    const Cx a_p(u(), -u());
    REQUIRE(max_abs(rtt::polar::prt_matrix(k, k_out, n, a_s, a_p) -
                    rtt::polar::prt_matrix(k, k_out, Vec3(-n), a_s, a_p)) <= 1e-15);
    REQUIRE((rtt::polar::geometric_transform(k, k_out, n, true) -
             rtt::polar::geometric_transform(k, k_out, Vec3(-n), true))
                .cwiseAbs()
                .maxCoeff() <= 1e-15);
  }
}

TEST_CASE("normal incidence: P is independent of the choice of s", "[prt]") {
  // Transverse isotropic J (prt.hpp): reflection with a_p = -a_s gives P = a_s (I - k k^T) - k k^T,
  // transmission with a_p = a_s gives P = a_s (I - k k^T) + k k^T, whatever s is. Also just below
  // and above kNormalIncidence (tilts 1e-13 and 1e-11 rad), where s is ill-conditioned, with the
  // Fresnel amplitudes of N-BK7. There the exact P differs from the normal-incidence form by
  // first order in the tilt: k_out turns by 2 tilt (reflection; less for refraction), which
  // moves the k_out and p_out columns by <= 2 tilt each, and the amplitudes differ by O(tilt^2);
  // so |P - P_normal| <= 4 tilt. Tolerance 4 tilt + 1e-14. Without making s exactly transverse
  // (prt_basis) the error at tilt 1e-11 would be ~1e-16 / 1e-11 = 1e-5.
  Uniform u(5);
  for (int i = 0; i < 50; ++i) {
    const Vec3 k = random_unit(u);
    const Mat3 kk = k * k.transpose();
    const CMat3 transverse = (Mat3::Identity() - kk).cast<Cx>();
    for (const double tilt : {0.0, 1e-13, 1e-11}) {
      const Vec3 side = k.unitOrthogonal();
      const Vec3 n = (k * std::cos(tilt) + side * std::sin(tilt)).normalized();
      const auto a = rtt::polar::fresnel(Cx(1.0), Cx(1.5168), std::sin(tilt));
      const CMat3 p_r = rtt::polar::prt_matrix(k, reflect(k, n), n, a.rs, a.rp);
      const CMat3 p_t = rtt::polar::prt_matrix(k, refract(k, n, 1.0, 1.5168), n, a.ts, a.tp);
      INFO("tilt " << tilt);
      const double tolerance = 4.0 * tilt + 1e-14;
      REQUIRE(max_abs(p_r - (a.rs * transverse - kk.cast<Cx>())) <= tolerance);
      REQUIRE(max_abs(p_t - (a.ts * transverse + kk.cast<Cx>())) <= tolerance);
    }
  }
}

TEST_CASE("rigid rotation: P(R k_in, R k_out, R N) = R P R^T", "[prt]") {
  // Invariance under rigid motion (architecture, test levels). Tolerance 1e-14.
  Uniform u(6);
  for (int i = 0; i < 100; ++i) {
    const Vec3 k = random_unit(u);
    const Vec3 n = random_unit(u);
    const Vec3 k_out = refract(k, n, 1.0, 1.5);
    const Mat3 r = random_rotation(u);
    const Cx a_s(u(), u());
    const Cx a_p(u(), u());
    const CMat3 p = rtt::polar::prt_matrix(k, k_out, n, a_s, a_p);
    const CMat3 rotated =
        rtt::polar::prt_matrix(Vec3(r * k), Vec3(r * k_out), Vec3(r * n), a_s, a_p);
    REQUIRE(max_abs(rotated - r.cast<Cx>() * p * r.transpose().cast<Cx>()) <= 1e-14);
  }
}

TEST_CASE("ideal mirror: no diattenuation, no retardance, E_r = -(I - 2 N N^T) E", "[prt]") {
  // An ideal conductor reflects E_r = -E_tan + E_norm = -(I - 2 N N^T) E, which in the basis
  // p = k x s is a_s = -1, a_p = +1 (docs/architecture.md, Polarisation und Fresnel). With
  // Lam's Q for reflection, Q^-1 P = O_in diag(-1, -1, 1) O_in^T: retardance 0. Tolerance 1e-14
  // for the matrices, 1e-12 for the decompositions.
  Uniform u(7);
  for (int i = 0; i < 100; ++i) {
    const Vec3 k = random_unit(u);
    const Vec3 n = random_unit(u);
    if (std::abs(k.dot(n)) < 0.05) continue;
    const Vec3 k_out = reflect(k, n);
    const CMat3 p = rtt::polar::prt_matrix(k, k_out, n, Cx(-1.0), Cx(1.0));
    const Mat3 mirror = -(Mat3::Identity() - 2.0 * n * n.transpose());
    const CVec3 e = random_transverse(u, k);
    REQUIRE((p * e - mirror.cast<Cx>() * e).norm() <= 1e-14);
    REQUIRE(rtt::polar::diattenuation(p, k, k_out).value <= 1e-12);
    const Mat3 q = rtt::polar::geometric_transform(k, k_out, n, true);
    REQUIRE(rtt::polar::physical_retardance(p, q, k).value <= 1e-12);
  }
}

TEST_CASE("consistency with #56: metal and glass at oblique incidence", "[prt]") {
  // Single Fresnel reflection: Q^-1 P = O_in diag(r_s, -r_p, 1) O_in^T, so the physical
  // retardance is |arg(-r_p / r_s)| = |Delta| of #56 and the fast axis is s if Delta > 0
  // (p delayed); the diattenuation from the singular values |r_s|, |r_p| equals the reflection
  // diattenuation of #56. Transmission: retardance |arg(t_p / t_s)|. Tolerance 1e-10
  // (architecture, acceptance table).
  const Vec3 k(0.0, std::sqrt(0.5), std::sqrt(0.5));  // 45 deg onto the plane z = 0
  const Vec3 n(0.0, 0.0, 1.0);
  const Vec3 k_r = reflect(k, n);
  const auto basis = rtt::polar::prt_basis(k, n);
  for (const Cx metal : {Cx(1.2, 7.26), Cx(0.2, 3.4)}) {
    const auto a = rtt::polar::fresnel(Cx(1.0), metal, std::sqrt(0.5));
    const auto power = rtt::polar::fresnel_power(Cx(1.0), metal, std::sqrt(0.5));
    const CMat3 p = rtt::polar::prt_matrix(k, k_r, n, a.rs, a.rp);
    const Mat3 q = rtt::polar::geometric_transform(k, k_r, n, true);
    const auto ret = rtt::polar::physical_retardance(p, q, k);
    const double delta = rtt::polar::reflection_phase_difference(a);
    INFO("metal " << metal << ": retardance " << ret.value << ", Delta " << delta);
    REQUIRE(std::abs(ret.value - std::abs(delta)) <= 1e-10);
    REQUIRE(delta > 0.0);
    REQUIRE(std::abs(std::abs(ret.fast_axis.dot(basis.s.cast<Cx>())) - 1.0) <= 1e-10);
    const auto d = rtt::polar::diattenuation(p, k, k_r);
    REQUIRE(std::abs(d.value - rtt::polar::reflection_diattenuation(power)) <= 1e-10);
    // R_s > R_p for a metal at 45 deg: the axis of maximum transmission is s.
    REQUIRE(power.reflectance_s > power.reflectance_p);
    REQUIRE(std::abs(std::abs(d.axis.dot(basis.s.cast<Cx>())) - 1.0) <= 1e-10);
  }
  const auto t = rtt::polar::fresnel(Cx(1.0), Cx(1.5168), std::sqrt(0.5));
  const Vec3 k_t = refract(k, n, 1.0, 1.5168);
  const CMat3 p_t = rtt::polar::prt_matrix(k, k_t, n, t.ts, t.tp);
  const Mat3 q_t = rtt::polar::geometric_transform(k, k_t, n, false);
  // t_s, t_p real > 0; Q^-1 P is real and its transverse eigenvalue may coincide with that of k:
  // the fast axis must still be transverse.
  const auto ret_t = rtt::polar::physical_retardance(p_t, q_t, k);
  REQUIRE(ret_t.value <= 1e-10);
  REQUIRE(std::abs(ret_t.fast_axis.dot(k.cast<Cx>())) <= 1e-14);
  const double ds = std::norm(t.ts);
  const double dp = std::norm(t.tp);
  REQUIRE(std::abs(rtt::polar::diattenuation(p_t, k, k_t).value - std::abs(ds - dp) / (ds + dp)) <=
          1e-10);
}

TEST_CASE("retardance with a degenerate eigenvalue: retarder diag(1, i) about a random k",
          "[prt]") {
  // m = O diag(1, i, 1) O^T with k = O e_z: a quarter-wave retarder whose fast axis O e_x has the
  // same eigenvalue 1 as k. Lam, Eq. (4.4): delta = arg(i) - arg(1) = pi/2, fast axis = the
  // eigenstate with the smaller phase, O e_x. The decomposition must stay in the transverse plane
  // (physik-reviewer of #57). Tolerance 1e-12.
  Uniform u(8);
  for (int i = 0; i < 50; ++i) {
    const Mat3 o = random_rotation(u);
    const Vec3 k = o.col(2);
    CMat3 diag = CMat3::Identity();
    diag(1, 1) = Cx(0.0, 1.0);
    const CMat3 m = o.cast<Cx>() * diag * o.transpose().cast<Cx>();
    const auto r = rtt::polar::retardance(m, k);
    REQUIRE(std::abs(r.value - kPi / 2.0) <= 1e-12);
    REQUIRE(std::abs(std::abs(r.fast_axis.dot(o.col(0).cast<Cx>())) - 1.0) <= 1e-12);
    REQUIRE(std::abs(r.fast_axis.dot(k.cast<Cx>())) <= 1e-12);
  }
}

TEST_CASE("golden: three prisms of Lam, Fig. 4.9 and Table 4.3", "[prt]") {
  // Lam, Sec. 4.5.1 and Fig. 4.9: a ray k = (0, 0, 1) through three glass prisms (n = 1.515) with
  // the six surface normals below exits in its incident direction; Lam gives
  //   Q = [[0.996, -0.092, 0], [0.092, 0.996, 0], [0, 0, 1]] (rotation by +5.27 deg about z),
  //   circular retardance of Q 10.546 deg (Table 4.3), and for Q^-1 P a linear diattenuator with
  //   maximum/minimum amplitude transmission 0.845/0.792 (uncoated, Fresnel).
  // Q is purely geometric (no sign convention involved). Tolerances a priori from the printed
  // digits, propagated through the nominal geometry (script lam_bound.py in the PR of #57,
  // independent of rtt): the normals have 3 decimals (component error <= 5e-4, direction error
  // dN <= sqrt(3) 5e-4 = 8.7e-4 rad). At each refraction d theta_t / d theta_i = A =
  // n_i cos_i / (n_t cos_t) (Snell) and out of the plane the factor is n_i / n_t, so the
  // direction error grows as delta_j = a_j delta_(j-1) + b_j dN with a = max(A, n_i/n_t),
  // b = max(|1 - A|, |1 - n_i/n_t|); here A = 0.62, 2.34, 0.50, 2.06, 0.42, 1.61 (incidence
  // 25.0, 33.8, 49.4, 30.7, 58.8, 16.4 deg), giving an exit direction within 6.05e-3 rad of z.
  // Each Q_j changes by <= 2 (delta_(j-1) + delta_j) + 2 dN, so ||dQ|| <= 0.0759 (alignment
  // included): entries of Q within 0.076 + 5e-4, rotation within 0.076 rad, retardance within
  // 0.152 rad = 8.7 deg. Fresnel amplitudes: relative change <= sum g_j (delta_(j-1) + dN) =
  // 0.013 with g = max |dt/dtheta| / t (1.63 at most), plus basis mixing <= ||dQ||: within
  // 0.845 * (0.013 + 0.076) = 0.075. These bounds are loose, so this golden test is qualitative
  // (fast axis, sign of the rotation) and only roughly quantitative; the retardance is secured
  // quantitatively by the property tests and the consistency with #56. The sign of the rotation
  // (Q(0,1) = -0.092, a sign error would be off by 0.18) and the fast axis of Q are checked
  // sharply. Lam, Table 4.3: the eigenvalues of Q are e^(-0.092 i) with eigenvector (1, i)/sqrt(2)
  // (fast) and e^(+0.092 i); after the alignment Q is a rotation about z, whose eigenvectors are
  // exactly (1, +-i)/sqrt(2) whatever the rounding of the normals: fast axis checked to 1e-12.
  const std::array<Vec3, 6> normals{Vec3(0.0, -0.423, 0.906),   Vec3(0.0, 0.423, 0.906),
                                    Vec3(0.366, 0.211, 0.906),  Vec3(-0.366, -0.211, 0.906),
                                    Vec3(-0.380, 0.198, 0.903), Vec3(0.378, -0.201, 0.903)};
  const double glass = 1.515;
  Vec3 k(0.0, 0.0, 1.0);
  const Vec3 k0 = k;
  CMat3 p = CMat3::Identity();
  Mat3 q = Mat3::Identity();
  for (std::size_t j = 0; j < 6; ++j) {
    const Vec3 n = normals[j].normalized();
    const bool entering = j % 2 == 0;
    const double n_i = entering ? 1.0 : glass;
    const double n_t = entering ? glass : 1.0;
    const Vec3 k_out = refract(k, n, n_i, n_t);
    const auto a = rtt::polar::fresnel_at_angle(Cx(n_i), Cx(n_t), std::abs(k.dot(n)));
    p = rtt::polar::prt_matrix(k, k_out, n, a.ts, a.tp) * p;  // Lam, Eq. (3.9)
    q = rtt::polar::geometric_transform(k, k_out, n, false) * q;
    k = k_out;
  }
  // Order of the chain (Lam, Eq. (3.9)): P_6 ... P_1 maps k0 onto the exit direction.
  REQUIRE((p * k0.cast<Cx>() - k.cast<Cx>()).norm() <= 1e-14);
  REQUIRE((q * k0 - k).norm() <= 1e-14);
  INFO("exit direction " << k.transpose());
  REQUIRE(std::acos(std::min(1.0, k.dot(k0))) <= 6.05e-3);
  // Rotate the small exit offset back onto z (minimal rotation), so that both matrices map z
  // onto itself as in Lam's figure.
  const Mat3 align = Eigen::Quaterniond::FromTwoVectors(k, k0).toRotationMatrix();
  const Mat3 qa = align * q;
  const CMat3 pa = align.cast<Cx>() * p;
  INFO("Q = \n" << qa);
  const double entry_tolerance = 0.076 + 5e-4;
  REQUIRE(std::abs(qa(0, 0) - 0.996) <= entry_tolerance);
  REQUIRE(std::abs(qa(0, 1) + 0.092) <= entry_tolerance);
  REQUIRE(std::abs(qa(1, 0) - 0.092) <= entry_tolerance);
  REQUIRE(std::abs(qa(1, 1) - 0.996) <= entry_tolerance);
  const auto ret_q = rtt::polar::retardance(qa.cast<Cx>(), k0);
  INFO("retardance of Q " << ret_q.value * 180.0 / kPi << " deg");
  REQUIRE(std::abs(ret_q.value * 180.0 / kPi - 10.546) <= 8.7);
  const CVec3 lam_fast = CVec3(Cx(1.0), Cx(0.0, 1.0), Cx(0.0)) / std::sqrt(2.0);
  REQUIRE(std::abs(std::abs(ret_q.fast_axis.dot(lam_fast)) - 1.0) <= 1e-12);
  const auto d = rtt::polar::diattenuation(Mat3(qa.transpose()).cast<Cx>() * pa, k0, k0);
  INFO("amplitude transmission " << d.maximum << " / " << d.minimum);
  REQUIRE(std::abs(d.maximum - 0.845) <= 0.075);
  REQUIRE(std::abs(d.minimum - 0.792) <= 0.075);
}

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif
