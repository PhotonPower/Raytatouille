// Power-normalised interface PRT matrices for the tracer (#61, ADR 0021). Source: S. J. Byrnes,
// arXiv:1603.02720v5, Eqs. (6), (21)-(23) (docs/quellen.md).

#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <complex>
#include <cstdint>
#include <numbers>
#include <random>
#include <utility>

#include "rtt/polar/fresnel.hpp"
#include "rtt/polar/interface.hpp"
#include "rtt/polar/prt_analysis.hpp"

using Cx = std::complex<double>;
using rtt::math::CMat3;
using rtt::math::Vec3;
using CVec3 = Eigen::Vector3cd;

namespace {

constexpr double kPi = std::numbers::pi;

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

/// Snell in vector form with Re(n) (B. de Greve, Eqs. (22), (23), (28); docs/quellen.md);
/// none for total internal reflection.
std::pair<bool, Vec3> refract(const Vec3& d, Vec3 n, double n1, double n2) {
  if (n.dot(d) > 0.0) n = -n;
  const double mu = n1 / n2;
  const double cos_i = -n.dot(d);
  const double sin_t2 = mu * mu * (1.0 - cos_i * cos_i);
  if (sin_t2 > 1.0) return {false, d};
  return {true, mu * d + (mu * cos_i - std::sqrt(1.0 - sin_t2)) * n};
}

Vec3 reflect(const Vec3& d, const Vec3& n) {
  return d - 2.0 * d.dot(n) * n;
}

CVec3 random_transverse(Uniform& u, const Vec3& k) {
  CVec3 e(Cx(u() - 0.5, u() - 0.5), Cx(u() - 0.5, u() - 0.5), Cx(u() - 0.5, u() - 0.5));
  const CVec3 kc = k.cast<Cx>();
  e -= kc * kc.dot(e);
  return e.normalized();
}

}  // namespace

TEST_CASE("tangential invariant and power factors", "[interface]") {
  // xi = Re(n_i) sin(theta_i) with sin = |k x N|; the power factors are T / |t|^2 of Byrnes,
  // Eqs. (21), (22) as implemented in fresnel_power (#56). Tolerance 1e-14 relative.
  const Vec3 n(0.0, 0.0, 1.0);
  for (const double deg : {0.0, 20.0, 55.0, 80.0}) {
    const double t = deg * kPi / 180.0;
    const Vec3 k(0.0, std::sin(t), std::cos(t));
    REQUIRE(std::abs(rtt::polar::tangential_invariant(k, n, Cx(1.5, 0.01)) - 1.5 * std::sin(t)) <=
            1e-15);
    for (const auto& [n_i, n_t] : {std::pair{Cx(1.0), Cx(1.5168)}, std::pair{Cx(1.0), Cx(1.5, 0.3)},
                                   std::pair{Cx(1.3, 0.02), Cx(1.7, 0.1)}}) {
      const double xi = n_i.real() * std::sin(t);
      const auto a = rtt::polar::fresnel(n_i, n_t, xi);
      const auto p = rtt::polar::fresnel_power(n_i, n_t, xi);
      const auto c = rtt::polar::transmission_power_factors(n_i, n_t, xi);
      INFO("theta " << deg << ", n_i " << n_i << ", n_t " << n_t);
      REQUIRE(std::abs(c.s * std::norm(a.ts) - p.transmittance_s) <= 1e-14 * p.transmittance_s);
      REQUIRE(std::abs(c.p * std::norm(a.tp) - p.transmittance_p) <= 1e-14 * p.transmittance_p);
      REQUIRE(c.s > 0.0);
      REQUIRE(c.p > 0.0);
    }
  }
}

TEST_CASE("power-normalised Fresnel PRT: |P_r E|^2 + |P_t E|^2 = |E|^2", "[interface]") {
  // With a' = a sqrt(c), |P E|^2 is the power fraction for every incident state; for a lossless
  // incident medium R + T = 1 for s and p (Byrnes, Eqs. (21)-(23)), also into an absorbing
  // medium, and s and p are orthonormal: |P_r E|^2 + |P_t E|^2 = 1 for |E| = 1. Total internal
  // reflection: |P_r E|^2 = 1. Random geometry with a fixed seed; tolerance 1e-12.
  Uniform u(11);
  for (const auto& [n_i, n_t] : {std::pair{Cx(1.0), Cx(1.5168)}, std::pair{Cx(1.5168), Cx(1.0)},
                                 std::pair{Cx(1.0), Cx(1.5, 0.3)}}) {
    for (int i = 0; i < 200; ++i) {
      const Vec3 k = random_unit(u);
      const Vec3 n = random_unit(u);
      if (std::abs(k.dot(n)) < 0.05) continue;
      const CVec3 e = random_transverse(u, k);
      const CMat3 p_r = rtt::polar::fresnel_prt(k, reflect(k, n), n, n_i, n_t, true);
      double power = (p_r * e).squaredNorm();
      const auto [transmitted, k_t] = refract(k, n, n_i.real(), n_t.real());
      if (transmitted) {
        const CMat3 p_t = rtt::polar::fresnel_prt(k, k_t, n, n_i, n_t, false);
        power += (p_t * e).squaredNorm();
        REQUIRE((p_t * k.cast<Cx>() - k_t.cast<Cx>()).norm() <= 1e-14);
      }
      INFO("n_i " << n_i << ", n_t " << n_t);
      REQUIRE(std::abs(power - 1.0) <= 1e-12);
    }
  }
}

TEST_CASE("power-normalised transmission: diattenuation from P equals the power diattenuation",
          "[interface]") {
  // With the power factors in P the singular values of P_T are sqrt(T_s), sqrt(T_p), so the
  // diattenuation of #57 equals the transmission diattenuation of #56 also for an absorbing
  // medium (where the field-amplitude diattenuation differs). diattenuation() gives |D| (Lam,
  // Eq. (4.3), singular values L1 >= L2), transmission_diattenuation() the signed
  // (T_s - T_p) / (T_s + T_p), negative here. Tolerance 1e-12.
  const Vec3 n(0.0, 0.0, 1.0);
  const double t = 50.0 * kPi / 180.0;
  const Vec3 k(0.0, std::sin(t), std::cos(t));
  for (const Cx n_t : {Cx(1.5168), Cx(1.5, 0.3)}) {
    const auto [ok, k_t] = refract(k, n, 1.0, n_t.real());
    REQUIRE(ok);
    const CMat3 p = rtt::polar::fresnel_prt(k, k_t, n, Cx(1.0), n_t, false);
    const auto power = rtt::polar::fresnel_power(Cx(1.0), n_t, std::sin(t));
    REQUIRE(std::abs(rtt::polar::diattenuation(p, k, k_t).value -
                     std::abs(rtt::polar::transmission_diattenuation(power))) <= 1e-12);
  }
}

TEST_CASE("interface_prt with unit power factors is prt_matrix", "[interface]") {
  Uniform u(12);
  for (int i = 0; i < 20; ++i) {
    const Vec3 k = random_unit(u);
    const Vec3 n = random_unit(u);
    const Vec3 k_r = reflect(k, n);
    const Cx a_s(u(), u());
    const Cx a_p(u(), -u());
    const CMat3 diff = rtt::polar::interface_prt(k, k_r, n, a_s, a_p, {1.0, 1.0}) -
                       rtt::polar::prt_matrix(k, k_r, n, a_s, a_p);
    REQUIRE(diff.cwiseAbs().maxCoeff() <= 1e-15);
  }
}
