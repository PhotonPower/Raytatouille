// Fresnel amplitudes with phase (#56). Source: S. J. Byrnes, "Multilayer optical
// calculations", arXiv:1603.02720v5 (docs/quellen.md); conventions in rtt/polar/fresnel.hpp and
// docs/architecture.md ("Polarisation und Fresnel").

#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <complex>
#include <limits>
#include <numbers>
#include <stdexcept>

#include "rtt/polar/fresnel.hpp"

using Cx = std::complex<double>;
using rtt::polar::fresnel;
using rtt::polar::fresnel_power;

namespace {

constexpr double kPi = std::numbers::pi;
constexpr double kNbk7 = 1.5168;  // N-BK7 at the d line (architecture, acceptance table)

double xi_of(double n_i, double theta) {
  return n_i * std::sin(theta);
}

}  // namespace

TEST_CASE("normal incidence on N-BK7: R = ((n - 1) / (n + 1))^2, r_p = -r_s", "[fresnel]") {
  // Byrnes, Eq. (6) at theta = 0: r_s = (n_i - n_t) / (n_i + n_t), r_p = -r_s (Convention A,
  // App. A); R = |r|^2 (Eq. (23)) ~ 4.22 % (architecture, acceptance table). Tolerance 1e-12.
  const auto a = fresnel(Cx(1.0), Cx(kNbk7), 0.0);
  const auto p = fresnel_power(Cx(1.0), Cx(kNbk7), 0.0);
  const double r = (1.0 - kNbk7) / (1.0 + kNbk7);
  REQUIRE(std::abs(a.rs - r) <= 1e-12);
  REQUIRE(std::abs(a.rp + r) <= 1e-12);
  REQUIRE(std::abs(p.reflectance_s - r * r) <= 1e-12);
  REQUIRE(std::abs(p.reflectance_p - r * r) <= 1e-12);
  REQUIRE(std::abs(p.reflectance_s - 0.0422) <= 5e-5);
  REQUIRE(std::abs(a.ts - 2.0 / (1.0 + kNbk7)) <= 1e-12);
  REQUIRE(std::abs(a.tp - 2.0 / (1.0 + kNbk7)) <= 1e-12);
  // Ellipsometric phase difference arg(-r_p / r_s) = 0 at normal incidence (Eq. (16)).
  REQUIRE(std::abs(rtt::polar::reflection_phase_difference(a)) <= 1e-15);
  REQUIRE(std::abs(rtt::polar::transmission_phase_difference(a)) <= 1e-15);
}

TEST_CASE("Brewster angle arctan(n): R_p = 0", "[fresnel]") {
  // Brewster angle theta_B = arctan(n_t / n_i) (architecture, acceptance table); there
  // n_t cos(theta_i) = n_i cos(theta_t) and r_p = 0 by Byrnes, Eq. (6). Tolerance 1e-12, from
  // both sides.
  for (const double n_i : {1.0, kNbk7}) {
    const double n_t = n_i == 1.0 ? kNbk7 : 1.0;
    const double theta_b = std::atan(n_t / n_i);
    const auto p = fresnel_power(Cx(n_i), Cx(n_t), xi_of(n_i, theta_b));
    INFO("n_i " << n_i << ": R_p " << p.reflectance_p);
    REQUIRE(p.reflectance_p <= 1e-12);
    REQUIRE(p.reflectance_s > 1e-3);
    REQUIRE(std::abs(rtt::polar::reflection_diattenuation(p) - 1.0) <= 1e-10);
  }
}

TEST_CASE("total internal reflection: critical angle and phase jumps", "[fresnel]") {
  // Glass -> air. Critical angle arcsin(n_t / n_i) (architecture, acceptance table): there
  // xi = n_i sin(theta_c) = n_t and q_t = 0, so r_s = r_p = 1 and T = 0 by Eq. (6). The test uses
  // the canonical input xi = n_t exactly (an angle computed with asin/sin lands a few ulp off,
  // where q_t ~ 1.5e-8 sqrt(ulps) makes |r| - 1 ~ 6e-8 sqrt(ulps)), and checks separately that
  // n_i sin(arcsin(n_t / n_i)) equals n_t to rounding (1e-15).
  // Beyond it q_t = i b with b = sqrt(xi^2 - n_t^2) (branch Im q >= 0, Byrnes D.3.1), and Eq. (6)
  // in q form gives
  //   r_s = (q_i - i b) / (q_i + i b), arg r_s = -2 atan(b / q_i),
  //   r_p = (n_t^2 q_i - n_i^2 i b) / (n_t^2 q_i + n_i^2 i b),
  //   arg r_p = -2 atan(n_i^2 b / (n_t^2 q_i)),
  // with |r_s| = |r_p| = 1 and T = 0 (Re q_t = 0). Tolerance 1e-12.
  const double n_i = kNbk7;
  const double n_t = 1.0;
  const double theta_c = std::asin(n_t / n_i);
  REQUIRE(std::abs(xi_of(n_i, theta_c) - n_t) <= 1e-15);
  {
    const auto a = fresnel(Cx(n_i), Cx(n_t), n_t);
    const auto p = fresnel_power(Cx(n_i), Cx(n_t), n_t);
    REQUIRE(std::abs(a.rs - 1.0) <= 1e-12);
    REQUIRE(std::abs(a.rp - 1.0) <= 1e-12);
    REQUIRE(p.transmittance_s <= 1e-12);
    REQUIRE(p.transmittance_p <= 1e-12);
    // Just below the critical angle light is transmitted.
    const auto below = fresnel_power(Cx(n_i), Cx(n_t), xi_of(n_i, theta_c - 1e-3));
    REQUIRE(below.transmittance_s > 1e-2);
  }
  for (const double theta_deg : {45.0, 50.0, 60.0, 75.0, 89.0}) {
    const double theta = theta_deg * kPi / 180.0;
    const double xi = xi_of(n_i, theta);
    const double q_i = n_i * std::cos(theta);
    const double b = std::sqrt(xi * xi - n_t * n_t);
    const auto a = fresnel(Cx(n_i), Cx(n_t), xi);
    const auto p = fresnel_power(Cx(n_i), Cx(n_t), xi);
    INFO("theta " << theta_deg << " deg");
    REQUIRE(std::abs(std::abs(a.rs) - 1.0) <= 1e-12);
    REQUIRE(std::abs(std::abs(a.rp) - 1.0) <= 1e-12);
    REQUIRE(std::abs(std::arg(a.rs) - (-2.0 * std::atan(b / q_i))) <= 1e-12);
    REQUIRE(std::abs(std::arg(a.rp) - (-2.0 * std::atan(n_i * n_i * b / (n_t * n_t * q_i)))) <=
            1e-12);
    REQUIRE(p.transmittance_s == 0.0);
    REQUIRE(p.transmittance_p == 0.0);
    // Evanescent wave: q_t = i b, decaying (Im q > 0).
    const Cx q_t = rtt::polar::normal_component(Cx(n_t), xi);
    REQUIRE(std::abs(q_t - Cx(0.0, b)) <= 1e-12);
  }
}

TEST_CASE("lossless media: R + T = 1 for s and p", "[fresnel]") {
  // Byrnes, Eqs. (21)-(23) with lossless incident medium. Both directions, also through the
  // critical angle; and an absorbing second medium (power crossing the interface is conserved
  // as long as the incident medium is lossless; this needs the conjugate cos* in T_p, Eq. (22)).
  // Tolerance 1e-12.
  const Cx metal(0.2, 3.4);
  for (const auto& [n_i, n_t] : {std::pair{Cx(1.0), Cx(kNbk7)}, std::pair{Cx(kNbk7), Cx(1.0)},
                                 std::pair{Cx(1.0), metal}, std::pair{Cx(1.0), Cx(1.5, 0.3)}}) {
    for (int deg = 0; deg < 90; deg += 3) {
      const double xi = n_i.real() * std::sin(deg * kPi / 180.0);
      const auto p = fresnel_power(n_i, n_t, xi);
      INFO("n_i " << n_i << ", n_t " << n_t << ", theta " << deg);
      REQUIRE(std::abs(p.reflectance_s + p.transmittance_s - 1.0) <= 1e-12);
      REQUIRE(std::abs(p.reflectance_p + p.transmittance_p - 1.0) <= 1e-12);
      REQUIRE(p.transmittance_s >= 0.0);
      REQUIRE(p.transmittance_p >= 0.0);
    }
  }
}

TEST_CASE("reciprocity: r_ab = -r_ba and t_ab t_ba - r_ab r_ba = 1", "[fresnel]") {
  // Identities that follow from Byrnes, Eq. (6) (stated after his Eq. (9)), for s and p, for
  // complex indices and the same tangential invariant on both sides. Tolerance 1e-12.
  const Cx a_index(1.3, 0.0);
  const Cx b_index(1.9, 0.4);
  for (int deg = 0; deg < 90; deg += 5) {
    const double xi = a_index.real() * std::sin(deg * kPi / 180.0);
    const auto ab = fresnel(a_index, b_index, xi);
    const auto ba = fresnel(b_index, a_index, xi);
    INFO("theta " << deg);
    REQUIRE(std::abs(ab.rs + ba.rs) <= 1e-12);
    REQUIRE(std::abs(ab.rp + ba.rp) <= 1e-12);
    REQUIRE(std::abs(ab.ts * ba.ts - ab.rs * ba.rs - 1.0) <= 1e-12);
    REQUIRE(std::abs(ab.tp * ba.tp - ab.rp * ba.rp - 1.0) <= 1e-12);
  }
}

TEST_CASE("metal at 45 deg: Abeles r_p = r_s^2, diattenuation and phase difference", "[fresnel]") {
  // From air at 45 deg (cos = sin = c, 2 c^2 = 1, xi^2 = 1/2), with q = sqrt(n^2 - 1/2), Eq. (6):
  //   r_s^2 = (c - q)^2 / (c + q)^2 = (n^2 - 2 c q) / (n^2 + 2 c q)   (c^2 + q^2 = n^2),
  //   r_p = (n^2 c - q) / (n^2 c + q) = (n^2 - 2 c q) / (n^2 + 2 c q)  (times 2c, 2 c^2 = 1),
  // so r_p = r_s^2 for every complex n (Abeles relation, holds in Convention A only). With it the
  // diattenuation and the ellipsometric phase difference follow from r_s alone:
  //   D = (|r_s|^2 - |r_s|^4) / (|r_s|^2 + |r_s|^4) = (1 - |r_s|^2) / (1 + |r_s|^2),
  //   Delta = arg(-r_p / r_s) = arg(-r_s).
  // r_s is evaluated here directly with std::sqrt (principal branch, Im >= 0 for these n).
  // Tolerance 1e-10 (architecture, acceptance table).
  const double c = std::sqrt(0.5);
  for (const Cx n : {Cx(1.2, 7.26), Cx(0.2, 3.4), Cx(kNbk7, 0.0)}) {
    const Cx q = std::sqrt(n * n - 0.5);
    const Cx rs = (c - q) / (c + q);
    const auto a = fresnel(Cx(1.0), n, c);
    const auto p = fresnel_power(Cx(1.0), n, c);
    INFO("n " << n);
    REQUIRE(std::abs(a.rs - rs) <= 1e-12);
    REQUIRE(std::abs(a.rp - rs * rs) <= 1e-12);
    const double r2 = std::norm(rs);
    REQUIRE(std::abs(rtt::polar::reflection_diattenuation(p) - (1.0 - r2) / (1.0 + r2)) <= 1e-10);
    REQUIRE(std::abs(rtt::polar::reflection_phase_difference(a) - std::arg(-rs)) <= 1e-10);
  }
}

TEST_CASE("metal with large kappa is an ideal mirror: r_s = -1, r_p = +1", "[fresnel]") {
  // For n_t = 1 + i kappa, kappa -> infinity: q_t ~ n_t, r_s = (c - q_t)/(c + q_t) = -1 + O(c /
  // kappa) and r_p = (n_t^2 c - q_t)/(n_t^2 c + q_t) = 1 + O(1 / (kappa c)) with c = cos(theta).
  // Independently, an ideal conductor reflects E_r = -E_tan + E_norm, which in the basis p = k x s
  // is a_s = -1, a_p = +1 at every angle (no diattenuation, no retardance; architecture.md). A
  // sign error in r_p (Convention B) would give r_p -> -1 here. Tolerance |r_s + 1| <= 3 / kappa,
  // |r_p - 1| <= 3 / (kappa c).
  const double kappa = 1e6;
  for (const double deg : {0.0, 30.0, 60.0, 85.0}) {
    const double cos_i = std::cos(deg * kPi / 180.0);
    const auto a = rtt::polar::fresnel_at_angle(Cx(1.0), Cx(1.0, kappa), cos_i);
    const auto p = fresnel_power(Cx(1.0), Cx(1.0, kappa), std::sin(deg * kPi / 180.0));
    INFO("theta " << deg);
    REQUIRE(std::abs(a.rs + 1.0) <= 3.0 / kappa);
    REQUIRE(std::abs(a.rp - 1.0) <= 3.0 / (kappa * cos_i));
    REQUIRE(std::abs(rtt::polar::reflection_phase_difference(a)) <= 6.0 / (kappa * cos_i));
    REQUIRE(std::abs(rtt::polar::reflection_diattenuation(p)) <= 6.0 / (kappa * cos_i));
  }
}

TEST_CASE("branch of the normal component: decaying or forward", "[fresnel]") {
  // Byrnes, App. D: Im q >= 0, and Re q >= 0 if Im q = 0.
  REQUIRE(rtt::polar::normal_component(Cx(1.5), 0.0) == Cx(1.5));
  const Cx absorbing = rtt::polar::normal_component(Cx(1.5, 0.2), 0.9);
  REQUIRE(absorbing.imag() > 0.0);
  REQUIRE(absorbing.real() > 0.0);
  REQUIRE(std::abs(absorbing * absorbing - (Cx(1.5, 0.2) * Cx(1.5, 0.2) - 0.81)) <= 1e-14);
  // A negative zero imaginary part must not flip the branch.
  const Cx tir = rtt::polar::normal_component(Cx(1.0, -0.0), 1.2);
  REQUIRE(tir.imag() > 0.0);
  REQUIRE(tir.real() == 0.0);
}

TEST_CASE("checked variants reject invalid input", "[fresnel]") {
  using rtt::polar::fresnel_checked;
  using rtt::polar::fresnel_power_checked;
  const double nan = std::numeric_limits<double>::quiet_NaN();
  REQUIRE_NOTHROW(fresnel_checked(Cx(1.0), Cx(1.5, 0.1), 0.5));
  REQUIRE_THROWS_AS(fresnel_checked(Cx(1.0), Cx(1.5, -0.1), 0.5), std::invalid_argument);
  REQUIRE_THROWS_AS(fresnel_checked(Cx(1.0, -1e-3), Cx(1.5), 0.5), std::invalid_argument);
  REQUIRE_THROWS_AS(fresnel_checked(Cx(0.0), Cx(1.5), 0.0), std::invalid_argument);
  REQUIRE_THROWS_AS(fresnel_checked(Cx(1.0), Cx(-1.5), 0.0), std::invalid_argument);
  REQUIRE_THROWS_AS(fresnel_checked(Cx(1.0), Cx(1.5), -0.1), std::invalid_argument);
  REQUIRE_THROWS_AS(fresnel_checked(Cx(1.0), Cx(1.5), 1.1), std::invalid_argument);
  REQUIRE_THROWS_AS(fresnel_checked(Cx(nan), Cx(1.5), 0.1), std::invalid_argument);
  REQUIRE_THROWS_AS(fresnel_checked(Cx(1.0), Cx(1.5), nan), std::invalid_argument);
  REQUIRE_THROWS_AS(fresnel_checked(Cx(1.2), Cx(1.2), 1.2), std::invalid_argument);
  REQUIRE_NOTHROW(fresnel_checked(Cx(1.0), Cx(1.5), 1.0));  // grazing, r = -1
  REQUIRE_THROWS_AS(fresnel_power_checked(Cx(1.0), Cx(1.5), 1.0), std::invalid_argument);
  REQUIRE_NOTHROW(fresnel_power_checked(Cx(1.0), Cx(1.5), 0.99));
}
