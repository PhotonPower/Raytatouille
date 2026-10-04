#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <complex>
#include <numbers>
#include <span>
#include <stdexcept>
#include <vector>

#include "rtt/coating/thickness.hpp"
#include "rtt/coating/transfer_matrix.hpp"

using Catch::Matchers::WithinAbs;
using Complex = std::complex<double>;
using rtt::coating::Amplitudes;
using rtt::coating::Layer;
using rtt::coating::Powers;

namespace {

constexpr double kPi = std::numbers::pi;

double deg(double d) {
  return d * kPi / 180.0;
}

Amplitudes<double> stack(Complex n_in,
                         const std::vector<Layer<double>>& layers,
                         Complex n_out,
                         double xi,
                         double wavelength_um) {
  rtt::coating::check_stack(n_in, layers, n_out, xi, wavelength_um);
  return rtt::coating::stack_amplitudes<double>(n_in, layers, n_out, xi, wavelength_um);
}

Powers<double> powers(Complex n_in,
                      const std::vector<Layer<double>>& layers,
                      Complex n_out,
                      double xi,
                      double wavelength_um) {
  return rtt::coating::stack_powers(stack(n_in, layers, n_out, xi, wavelength_um), n_in, n_out, xi);
}

/// The same layers seen from the exit side: reversed order (light from the substrate).
std::vector<Layer<double>> reversed(std::vector<Layer<double>> layers) {
  return {layers.rbegin(), layers.rend()};
}

void require_close(Complex actual, Complex expected, double tol) {
  REQUIRE_THAT(actual.real(), WithinAbs(expected.real(), tol));
  REQUIRE_THAT(actual.imag(), WithinAbs(expected.imag(), tol));
}

/// Quarter-wave stack (HL)^pairs H at lambda0 = 0.55 um: H = 2.35 (TiO2-like), L = 1.38
/// (MgF2-like), lossless.
std::vector<Layer<double>> hl_stack(int pairs) {
  const double l0 = 0.55;
  std::vector<Layer<double>> layers;
  for (int k = 0; k < pairs; ++k) {
    layers.push_back({2.35, rtt::coating::quarter_wave_thickness_um(1.0, l0, 2.35)});
    layers.push_back({1.38, rtt::coating::quarter_wave_thickness_um(1.0, l0, 1.38)});
  }
  layers.push_back({2.35, rtt::coating::quarter_wave_thickness_um(1.0, l0, 2.35)});
  return layers;
}

}  // namespace

TEST_CASE("kernels are noexcept", "[coating]") {
  // Rule 3: no exceptions in tracing loops. (Lvalue arguments: the std::complex constructors of
  // the standard library are not noexcept themselves.)
  const Complex n(1.5);
  const std::span<const Layer<double>> layers;
  const Amplitudes<double> a;
  const double xi = 0.0;
  static_assert(noexcept(rtt::coating::stack_amplitudes<double>(n, layers, n, xi, xi)));
  static_assert(noexcept(rtt::coating::stack_powers<double>(a, n, n, xi)));
  static_assert(noexcept(rtt::coating::normal_component<double>(n, xi)));
  static_assert(noexcept(rtt::coating::interface_amplitudes<double>(n, n, xi)));
  SUCCEED();
}

TEST_CASE("branch of q = n cos(theta) (Byrnes App. D)", "[coating]") {
  // Real index, no TIR: q real and positive.
  const Complex q_glass = rtt::coating::normal_component<double>(1.5, 0.5);
  REQUIRE(q_glass.imag() == 0.0);
  REQUIRE_THAT(q_glass.real(), WithinAbs(std::sqrt(1.5 * 1.5 - 0.25), 1e-14));
  // TIR (xi > n): q imaginary with Im q > 0, the evanescent wave decays.
  const Complex q_tir = rtt::coating::normal_component<double>(1.0, 1.2);
  REQUIRE(q_tir.real() == 0.0);
  REQUIRE_THAT(q_tir.imag(), WithinAbs(std::sqrt(1.44 - 1.0), 1e-14));
  // Absorbing medium: Im q > 0 and (theorem of App. D.2) Re q > 0.
  const Complex q_metal = rtt::coating::normal_component<double>({0.2, 3.0}, 0.9);
  REQUIRE(q_metal.imag() > 0.0);
  REQUIRE(q_metal.real() > 0.0);
}

TEST_CASE("quarter-wave MgF2 anti-reflection on n = 1.52: R = 1.26 %", "[coating]") {
  // Architecture, reference case (1e-10). Derivation from Byrnes Eq. (6), (8), (11), (13), (15)
  // for N = 3 media (air n0 = 1, MgF2 n1 = 1.38, glass n2 = 1.52), normal incidence:
  // a quarter-wave layer has delta = 2 pi n1 d / lambda0 = pi / 2 (Eq. (8) with d = lambda0 /
  // (4 n1)), so diag(e^{-i delta}, e^{i delta}) = diag(-i, i) and
  //   M~ = 1 / (t01 t12) [[1, r01], [r01, 1]] [[-i, -i r12], [i r12, i]]
  //      = 1 / (t01 t12) [[-i + i r01 r12, ...], [-i r01 + i r12, ...]],
  //   r = M~10 / M~00 = (r01 - r12) / (1 - r01 r12).
  // With r01 = (n0 - n1) / (n0 + n1), r12 = (n1 - n2) / (n1 + n2) (Eq. (6), cos = 1):
  //   r01 - r12 = 2 (n0 n2 - n1^2) / D, 1 - r01 r12 = 2 (n0 n2 + n1^2) / D,
  //   D = (n0 + n1)(n1 + n2), so r = (n0 n2 - n1^2) / (n0 n2 + n1^2) and R = r^2.
  const double n0 = 1.0;
  const double n1 = 1.38;
  const double n2 = 1.52;
  const double l0 = 0.55;
  const double d = rtt::coating::physical_thickness_um(rtt::coating::QuarterWaves{1.0, l0}, n1);
  REQUIRE_THAT(d, WithinAbs(l0 / (4.0 * n1), 1e-16));
  const std::vector<Layer<double>> ar = {{n1, d}};
  const double r_expected = (n0 * n2 - n1 * n1) / (n0 * n2 + n1 * n1);
  const double r2 = r_expected * r_expected;
  REQUIRE_THAT(r2, WithinAbs(0.0126, 5e-5));  // "R = 1.26 %" (architecture), rounded
  const Powers<double> p = powers(n0, ar, n2, 0.0, l0);
  REQUIRE_THAT(p.reflectance_s, WithinAbs(r2, 1e-10));
  REQUIRE_THAT(p.reflectance_p, WithinAbs(r2, 1e-10));
  const Amplitudes<double> a = stack(n0, ar, n2, 0.0, l0);
  require_close(a.rs, r_expected, 1e-12);
  require_close(a.rp, -r_expected, 1e-12);  // Convention A: r_p = -r_s at normal incidence
  REQUIRE_THAT(p.reflectance_s + p.transmittance_s, WithinAbs(1.0, 1e-12));
}

TEST_CASE("lossless stack: R + T = 1 for s and p", "[coating]") {
  // Architecture, reference case (1e-12); Byrnes Eq. (21)-(23) without absorption.
  const std::vector<Layer<double>> layers = hl_stack(5);
  for (const double theta : {0.0, 15.0, 30.0, 45.0, 60.0, 70.0}) {
    for (const double wl : {0.40, 0.48, 0.55, 0.63, 0.80}) {
      const double xi = rtt::coating::tangential_invariant(1.0, deg(theta));
      const Powers<double> p = powers(1.0, layers, 1.52, xi, wl);
      INFO("theta " << theta << " deg, wavelength " << wl << " um");
      REQUIRE_THAT(p.reflectance_s + p.transmittance_s, WithinAbs(1.0, 1e-12));
      REQUIRE_THAT(p.reflectance_p + p.transmittance_p, WithinAbs(1.0, 1e-12));
      REQUIRE(p.reflectance_s > 0.0);
    }
  }
  // At the design wavelength and normal incidence the quarter-wave stack is a good mirror.
  REQUIRE(powers(1.0, layers, 1.52, 0.0, 0.55).reflectance_s > 0.9);
}

TEST_CASE("no layers: Fresnel amplitudes of Eq. (6)", "[coating]") {
  // Architecture: zero layers equal Fresnel (#56; until its merge against Eq. (6) written out
  // here). Air -> glass n = 1.5168 at several angles.
  const double n1 = 1.0;
  const double n2 = 1.5168;
  for (const double theta : {0.0, 20.0, 45.0, 70.0, 85.0}) {
    const double c1 = std::cos(deg(theta));
    const double s1 = std::sin(deg(theta));
    const double c2 = std::sqrt(1.0 - (n1 * s1 / n2) * (n1 * s1 / n2));
    const Amplitudes<double> a = stack(n1, {}, n2, n1 * s1, 0.55);
    INFO("theta " << theta);
    require_close(a.rs, (n1 * c1 - n2 * c2) / (n1 * c1 + n2 * c2), 1e-14);
    require_close(a.rp, (n2 * c1 - n1 * c2) / (n2 * c1 + n1 * c2), 1e-14);
    require_close(a.ts, 2.0 * n1 * c1 / (n1 * c1 + n2 * c2), 1e-14);
    require_close(a.tp, 2.0 * n1 * c1 / (n2 * c1 + n1 * c2), 1e-14);
  }
  // Normal incidence on N-BK7: R = 4.22 % (architecture), r_p = -r_s.
  const Amplitudes<double> normal = stack(1.0, {}, n2, 0.0, 0.55);
  REQUIRE_THAT(std::norm(normal.rs), WithinAbs(std::pow((n2 - 1.0) / (n2 + 1.0), 2), 1e-14));
  REQUIRE(std::abs(normal.rp + normal.rs) < 1e-14);
  // Brewster angle arctan(n2 / n1): r_p = 0.
  const double brewster = std::atan(n2 / n1);
  REQUIRE(std::abs(stack(n1, {}, n2, n1 * std::sin(brewster), 0.55).rp) < 1e-14);
}

TEST_CASE("total internal reflection through a lossless stack", "[coating]") {
  // Glass (1.52) -> stack -> air at 60 deg: beyond the critical angle of air the exit medium
  // carries no power (q_out imaginary, Eq. (21)/(22) give T = 0) and |r| = 1 without loss.
  const std::vector<Layer<double>> layers = reversed(hl_stack(2));
  const double xi = rtt::coating::tangential_invariant(1.52, deg(60.0));
  const Amplitudes<double> a = stack(1.52, layers, 1.0, xi, 0.55);
  const Powers<double> p = rtt::coating::stack_powers(a, Complex(1.52), Complex(1.0), xi);
  REQUIRE_THAT(std::abs(a.rs), WithinAbs(1.0, 1e-12));
  REQUIRE_THAT(std::abs(a.rp), WithinAbs(1.0, 1e-12));
  REQUIRE(p.transmittance_s == 0.0);
  REQUIRE(p.transmittance_p == 0.0);
}

TEST_CASE("half-wave layer (delta = pi) is absent, also at oblique incidence", "[coating]") {
  // Derivation from Eq. (8), (11), (13), (15): delta = pi gives diag(-1, -1), so
  //   M~ = -1 / (t01 t12) [[1, r01], [r01, 1]] [[1, r12], [r12, 1]],
  //   r = (r01 + r12) / (1 + r01 r12), t = -t01 t12 / (1 + r01 r12).
  // Eq. (6) with cos = q / n reads r_s,ij = (q_i - q_j) / (q_i + q_j) and r_p,ij = (a_i - a_j) /
  // (a_i + a_j) with a = q / n^2; t_s,ij = 1 + r_s,ij and t_p,ij = (n_i / n_j)(1 + r_p,ij). With
  // x = q (s) or a (p): r01 + r12 = 2 x1 (x0 - x2) / D, 1 + r01 r12 = 2 x1 (x0 + x2) / D,
  // D = (x0 + x1)(x1 + x2), so r = (x0 - x2) / (x0 + x2) = r02, and (1 + r01)(1 + r12) /
  // (1 + r01 r12) = 1 + r02 gives t = -t02 for s and p (the factors n0/n1 n1/n2 = n0/n2 of
  // t_p match). The layer adds the phase pi and nothing else, at every angle, if its
  // thickness is d = lambda / (2 q1) (delta = pi by Eq. (8) with (2)).
  const double n0 = 1.0;
  const double n1 = 2.0;
  const double n2 = 1.52;
  const double l0 = 0.6;
  for (const double theta : {0.0, 30.0, 50.0, 75.0}) {
    const double xi = std::sin(deg(theta));
    const double q1 = std::sqrt(n1 * n1 - xi * xi);
    const std::vector<Layer<double>> half = {{n1, l0 / (2.0 * q1)}};
    const Amplitudes<double> with = stack(n0, half, n2, xi, l0);
    const Amplitudes<double> bare = stack(n0, {}, n2, xi, l0);
    INFO("theta " << theta);
    require_close(with.rs, bare.rs, 1e-14);
    require_close(with.rp, bare.rp, 1e-14);
    require_close(with.ts, -bare.ts, 1e-14);
    require_close(with.tp, -bare.tp, 1e-14);
  }
  // At normal incidence this is 2 QWOT.
  REQUIRE_THAT(rtt::coating::quarter_wave_thickness_um(2.0, l0, n1),
               WithinAbs(l0 / (2.0 * n1), 1e-16));
}

TEST_CASE("single layer: r = (r01 + r12 e^{2i delta}) / (1 + r01 r12 e^{2i delta})", "[coating]") {
  // Eq. (11), (13), (15) for one layer: M~ = 1 / (t01 t12) [[1, r01], [r01, 1]]
  // [[e^{-i delta}, r12 e^{-i delta}], [r12 e^{i delta}, e^{i delta}]], so
  //   M~00 = (e^{-i delta} + r01 r12 e^{i delta}) / (t01 t12),
  //   M~10 = (r01 e^{-i delta} + r12 e^{i delta}) / (t01 t12),
  //   r = (r01 + r12 e^{2i delta}) / (1 + r01 r12 e^{2i delta}),
  //   t = t01 t12 e^{i delta} / (1 + r01 r12 e^{2i delta}),
  // with delta = 2 pi q1 d / lambda (Eq. (8) with (2)). This fixes the phase sign (e^{+i delta}
  // for the forward wave, exp(i(k.r - omega t))) and the cos(theta) in delta at oblique
  // incidence, also for an absorbing layer.
  const std::complex<double> i(0.0, 1.0);
  for (const Complex n1 : {Complex(2.1, 0.0), Complex(1.8, 0.04)}) {
    for (const double theta : {0.0, 35.0, 65.0}) {
      for (const double d : {0.037, 0.137, 0.61}) {
        const double wl = 0.55;
        const double xi = std::sin(deg(theta));
        const Complex delta = 2.0 * kPi * rtt::coating::normal_component<double>(n1, xi) * d / wl;
        const Amplitudes<double> i01 = rtt::coating::interface_amplitudes<double>(1.0, n1, xi);
        const Amplitudes<double> i12 = rtt::coating::interface_amplitudes<double>(n1, 1.52, xi);
        const Complex e1 = std::exp(i * delta);
        const Complex e2 = std::exp(2.0 * i * delta);
        const Amplitudes<double> a = stack(1.0, {{n1, d}}, 1.52, xi, wl);
        INFO("n1 " << n1 << ", theta " << theta << ", d " << d);
        require_close(a.rs, (i01.rs + i12.rs * e2) / (1.0 + i01.rs * i12.rs * e2), 1e-14);
        require_close(a.rp, (i01.rp + i12.rp * e2) / (1.0 + i01.rp * i12.rp * e2), 1e-14);
        require_close(a.ts, i01.ts * i12.ts * e1 / (1.0 + i01.rs * i12.rs * e2), 1e-14);
        require_close(a.tp, i01.tp * i12.tp * e1 / (1.0 + i01.rp * i12.rp * e2), 1e-14);
      }
    }
  }
  // Index-matched layer: r = 0 and t = e^{i delta}, the bare propagation phase.
  const double xi = std::sin(deg(40.0));
  const double q = std::sqrt(1.5 * 1.5 - xi * xi);
  const Amplitudes<double> matched = stack(1.5, {{1.5, 0.2}}, 1.5, xi, 0.55);
  require_close(matched.rs, 0.0, 1e-14);
  require_close(matched.ts, std::exp(i * (2.0 * kPi * q * 0.2 / 0.55)), 1e-14);
  require_close(matched.tp, std::exp(i * (2.0 * kPi * q * 0.2 / 0.55)), 1e-14);
}

TEST_CASE("absorbing substrate: R + T = 1 for a non-absorbing incident medium and stack",
          "[coating]") {
  // Byrnes App. B: for light from one side the power entering the last interface is T, and
  // with a real incident index the power entering the stack is 1 - R; without absorption in the
  // layers both are the same power flow (Eq. (18)/(20) are continuous), so R + T = 1. This
  // checks the conjugation in Eq. (22) (without it R_p + T_p = 0.98 at 60 deg here).
  const Complex metal(0.2, 3.0);
  for (const auto& layers : {std::vector<Layer<double>>{}, hl_stack(1)}) {
    for (const double theta : {0.0, 30.0, 60.0, 80.0}) {
      const Powers<double> p = powers(1.0, layers, metal, std::sin(deg(theta)), 0.55);
      INFO("layers " << layers.size() << ", theta " << theta);
      REQUIRE_THAT(p.reflectance_s + p.transmittance_s, WithinAbs(1.0, 1e-12));
      REQUIRE_THAT(p.reflectance_p + p.transmittance_p, WithinAbs(1.0, 1e-12));
      REQUIRE(p.transmittance_p > 0.0);
    }
  }
}

TEST_CASE("thick absorbing layer reflects like its first interface", "[coating]") {
  // Layer n1 = 1.5 + 0.1i, 50 um at 0.5 um: Im delta = 2 pi 0.1 50 / 0.5 = 62.8, so in M~ the
  // backward term e^{i delta} is e^{-2 Im delta} ~ 1e-55 smaller than e^{-i delta}: r = r01
  // (Eq. (13), (15)) up to that factor, and nothing is transmitted.
  const Complex n1(1.5, 0.1);
  for (const double theta : {0.0, 40.0}) {
    const double xi = std::sin(deg(theta));
    const Amplitudes<double> a = stack(1.0, {{n1, 50.0}}, 1.52, xi, 0.5);
    const Amplitudes<double> first = rtt::coating::interface_amplitudes<double>(1.0, n1, xi);
    require_close(a.rs, first.rs, 1e-14);
    require_close(a.rp, first.rp, 1e-14);
    REQUIRE(std::abs(a.ts) < 1e-25);
    REQUIRE(std::abs(a.tp) < 1e-25);
  }
}

TEST_CASE("thin absorbing layer: 0 < A < 1", "[coating]") {
  // 20 nm of a silver-like metal (n = 0.05 + 3.0i at 0.55 um) on glass.
  const std::vector<Layer<double>> metal = {{{0.05, 3.0}, 0.020}};
  for (const double theta : {0.0, 30.0, 60.0}) {
    const Powers<double> p = powers(1.0, metal, 1.52, std::sin(deg(theta)), 0.55);
    INFO("theta " << theta);
    REQUIRE(p.absorptance_s() > 0.0);
    REQUIRE(p.absorptance_p() > 0.0);
    REQUIRE(p.reflectance_s + p.transmittance_s < 1.0);
    REQUIRE(p.transmittance_s > 0.0);
  }
}

TEST_CASE("single interface: r_ab = -r_ba and t_ab t_ba - r_ab r_ba = 1", "[coating]") {
  // The identities Byrnes uses after Eq. (9) (they follow from Eq. (6)); they hold for one
  // interface only, not for a stack (a film in air has r_ab = r_ba).
  const Complex na = 1.0;
  const Complex nb(1.7, 0.02);
  for (const double theta : {0.0, 30.0, 75.0}) {
    const double xi = std::sin(deg(theta));
    const Amplitudes<double> ab = stack(na, {}, nb, xi, 0.55);
    const Amplitudes<double> ba = rtt::coating::stack_amplitudes<double>(nb, {}, na, xi, 0.55);
    INFO("theta " << theta);
    require_close(ab.rs, -ba.rs, 1e-14);
    require_close(ab.rp, -ba.rp, 1e-14);
    require_close(ab.ts * ba.ts - ab.rs * ba.rs, 1.0, 1e-14);
    require_close(ab.tp * ba.tp - ab.rp * ba.rp, 1.0, 1e-14);
  }
}

TEST_CASE("stack reciprocity", "[coating]") {
  // Derivation from Byrnes Eq. (6), (11), (13): as a two-port, M~ = (1/t) [[1, -r'], [r, t t' -
  // r r']] with t', r' for light from the other side, so det M~ = t'/t. On the other hand
  // det M~ = prod_k t_{k+1,k} / t_{k,k+1} (phase matrices have det 1, det [[1, r], [r, 1]] / t^2
  // = (1 - r^2) / t^2 = t_ba t_ab / t_ab^2 with the single-interface identities). By Eq. (6)
  // each factor is q_{k+1} / q_k for s (t_s ratio n2 cos2 / n1 cos1) and for p (t_p ratio
  // (2 n2 cos2) / (2 n1 cos1)), so the product telescopes:
  //   (1) t_ba / t_ab = q_out / q_in for s and p, also with absorbing layers.
  // With real outer media Eq. (21)/(22) then give
  //   (2) T_ab = T_ba (also with absorption in the stack), and for lossless stacks
  //   (3) |r_ab| = |r_ba| because R = 1 - T on both sides.
  // The single-interface relations r_ab = -r_ba, t t' - r r' = 1 do NOT hold for stacks.
  const double n_in = 1.0;
  const double n_out = 1.52;
  std::vector<Layer<double>> absorbing = hl_stack(2);
  absorbing.insert(absorbing.begin() + 2, Layer<double>{{1.9, 0.05}, 0.07});
  for (const double theta : {0.0, 25.0, 50.0}) {
    for (const double wl : {0.45, 0.6}) {
      const double xi = n_in * std::sin(deg(theta));
      INFO("theta " << theta << ", wavelength " << wl);
      for (const auto& layers : {hl_stack(3), absorbing}) {
        const Amplitudes<double> ab = stack(n_in, layers, n_out, xi, wl);
        const Amplitudes<double> ba = stack(n_out, reversed(layers), n_in, xi, wl);
        const Complex q_in = rtt::coating::normal_component<double>(n_in, xi);
        const Complex q_out = rtt::coating::normal_component<double>(n_out, xi);
        require_close(ba.ts / ab.ts, q_out / q_in, 1e-12);  // (1)
        require_close(ba.tp / ab.tp, q_out / q_in, 1e-12);
        const Powers<double> pab =
            rtt::coating::stack_powers(ab, Complex(n_in), Complex(n_out), xi);
        const Powers<double> pba =
            rtt::coating::stack_powers(ba, Complex(n_out), Complex(n_in), xi);
        REQUIRE_THAT(pab.transmittance_s, WithinAbs(pba.transmittance_s, 1e-12));  // (2)
        REQUIRE_THAT(pab.transmittance_p, WithinAbs(pba.transmittance_p, 1e-12));
      }
      const Amplitudes<double> ab = stack(n_in, hl_stack(3), n_out, xi, wl);
      const Amplitudes<double> ba = stack(n_out, reversed(hl_stack(3)), n_in, xi, wl);
      REQUIRE_THAT(std::abs(ab.rs), WithinAbs(std::abs(ba.rs), 1e-12));  // (3)
      REQUIRE_THAT(std::abs(ab.rp), WithinAbs(std::abs(ba.rp), 1e-12));
      // A film in air looks the same from both sides: r_ab = r_ba != 0, so r_ab != -r_ba.
      const std::vector<Layer<double>> film = {{1.6, 0.3}};
      const Amplitudes<double> front = stack(1.0, film, 1.0, xi, wl);
      const Amplitudes<double> back = stack(1.0, reversed(film), 1.0, xi, wl);
      require_close(front.rs, back.rs, 1e-14);
      REQUIRE(std::abs(front.rs) > 0.01);
    }
  }
}

TEST_CASE("metal with large kappa approaches the ideal mirror r_s = -1, r_p = +1", "[coating]") {
  // Eq. (6) with |n2| -> infinity: cos2 -> 1, r_s -> -1 and r_p -> +1 (Convention A); the
  // deviation is O(1/|n2|).
  for (const double theta : {0.0, 30.0, 60.0}) {
    const double xi = std::sin(deg(theta));
    INFO("theta " << theta);
    const Amplitudes<double> metal = stack(1.0, {}, Complex(1.0, 1e5), xi, 0.55);
    require_close(metal.rs, -1.0, 1e-4);
    require_close(metal.rp, 1.0, 1e-4);
    const Amplitudes<double> conductor = stack(1.0, {}, Complex(1e8, 0.0), xi, 0.55);
    require_close(conductor.rs, -1.0, 1e-7);
    require_close(conductor.rp, 1.0, 1e-7);
  }
}

TEST_CASE("check_stack rejects invalid input at the API boundary", "[coating]") {
  const std::vector<Layer<double>> ok = {{1.38, 0.1}};
  REQUIRE_NOTHROW(rtt::coating::check_stack(1.0, ok, 1.52, 0.5, 0.55));
  REQUIRE_THROWS_AS(rtt::coating::check_stack(1.0, ok, 1.52, 0.5, 0.0), std::invalid_argument);
  const std::vector<Layer<double>> negative_thickness = {{{1.38, 0.1}, -0.1}};
  REQUIRE_THROWS_AS(rtt::coating::check_stack(1.0, negative_thickness, 1.52, 0.5, 0.55),
                    std::invalid_argument);
  const std::vector<Layer<double>> gain = {{{1.38, -0.1}, 0.1}};
  REQUIRE_THROWS_AS(rtt::coating::check_stack(1.0, gain, 1.52, 0.5, 0.55), std::invalid_argument);
  REQUIRE_THROWS_AS(rtt::coating::check_stack(1.0, ok, 0.0, 0.5, 0.55), std::invalid_argument);
  REQUIRE_THROWS_AS(rtt::coating::check_stack(1.0, ok, 1.52, 1.0, 0.55),
                    std::invalid_argument);  // grazing
  REQUIRE_THROWS_AS(rtt::coating::check_stack(1.0, ok, 1.52, -0.1, 0.55), std::invalid_argument);
  REQUIRE_THROWS_AS(rtt::coating::check_stack(1.0, ok, 1.52, std::nan(""), 0.55),
                    std::invalid_argument);
  // Re n <= 0 (negative-index media, Byrnes App. D.4/D.5) is not supported.
  REQUIRE_THROWS_AS(rtt::coating::check_stack(1.0, ok, Complex(-1.5, 0.1), 0.5, 0.55),
                    std::invalid_argument);
  // A non-absorbing layer exactly at its critical angle (q = 0) would divide by t = 0.
  const std::vector<Layer<double>> critical = {{1.2, 0.1}};
  REQUIRE_THROWS_AS(rtt::coating::check_stack(1.5, critical, 1.52, 1.2, 0.55),
                    std::invalid_argument);
  REQUIRE_NOTHROW(rtt::coating::check_stack(1.5, critical, 1.52, 1.3, 0.55));  // evanescent
  // An absorbing incident medium allows any real xi (Byrnes Sec. 2: n sin(theta) real).
  REQUIRE_NOTHROW(rtt::coating::check_stack(Complex(1.0, 0.1), ok, 1.52, 1.0, 0.55));
}
