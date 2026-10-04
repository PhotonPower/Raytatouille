#pragma once

/// @file fresnel.hpp
/// Fresnel amplitude coefficients with phase at a planar interface between two isotropic,
/// non-magnetic media with complex refractive indices n = n' + i kappa (kappa >= 0), including
/// total internal reflection (docs/architecture.md, Konventionen, "Polarisation und Fresnel").
///
/// Source: S. J. Byrnes, "Multilayer optical calculations", arXiv:1603.02720v5 (2020); see
/// docs/quellen.md. Conventions (decided for #56 and #58):
/// - Fields ~ exp(i(k.r - omega t)); kappa > 0 is absorption; a larger phase arg means a delay.
/// - Tangential invariant xi = n sin(theta), real and equal in both media (Byrnes, Eq. (3)).
///   For an absorbing incident medium xi = Re(n_i) sin(theta_i), as in the refraction of the
///   tracer.
/// - Normal component q_j = n_j cos(theta_j) = sqrt(n_j^2 - xi^2) with Im q >= 0, and Re q >= 0
///   if Im q = 0: the transmitted wave decays, or carries power forwards (Byrnes, App. D).
/// - "Convention A" for p (Byrnes, Eq. (6) and App. A; App. A, footnote 6 lists textbooks that
///   use it):
///   r_s = (n_i cos_i - n_t cos_t) / (n_i cos_i + n_t cos_t),
///   r_p = (n_t cos_i - n_i cos_t) / (n_t cos_i + n_i cos_t),
///   t_s = 2 n_i cos_i / (n_i cos_i + n_t cos_t),  t_p = 2 n_i cos_i / (n_t cos_i + n_i cos_t).
///   With the interface normal N pointing into the medium t, the p unit vectors of Byrnes,
///   Eq. (5) are p = k x s for the incident, reflected and transmitted wave with the same
///   s = k_in x N / |k_in x N| (with N reversed, s and all p change sign together and r, t stay
///   the same), so r and t are directly the
///   amplitudes a_s, a_p of the polarization ray tracing matrix J = diag(a_s, a_p, 1)
///   (architecture.md, rtt-polar). r_p = -r_s at normal incidence; an ideal mirror without
///   retardance has r_s = -1, r_p = +1.

#include <cassert>
#include <cmath>
#include <complex>
#include <numbers>

#include "rtt/math/real.hpp"

namespace rtt::polar {

/// Fresnel amplitude coefficients (dimensionless, complex) of one interface, from the incident
/// medium i into the medium t, in Convention A (file comment).
template <math::Real T>
struct FresnelAmplitudes {
  std::complex<T> rs;  ///< reflection, s (E perpendicular to the plane of incidence)
  std::complex<T> rp;  ///< reflection, p (amplitude along p = k x s)
  std::complex<T> ts;  ///< transmission, s
  std::complex<T> tp;  ///< transmission, p
};

/// Power fractions of one interface, relative to the incident power (dimensionless).
template <math::Real T>
struct FresnelPower {
  T reflectance_s = T(0);    ///< R_s = |r_s|^2
  T reflectance_p = T(0);    ///< R_p = |r_p|^2
  T transmittance_s = T(0);  ///< T_s = |t_s|^2 Re(q_t) / Re(q_i)
  T transmittance_p = T(0);  ///< T_p = |t_p|^2 Re(n_t cos*_t) / Re(n_i cos*_i)
};

/// Normal component q = n cos(theta) = sqrt(n^2 - xi^2) of the wave vector in units of the
/// vacuum wave number, on the branch Im q >= 0 and, if Im q = 0, Re q >= 0 (Byrnes, App. D.2,
/// D.3). Total internal reflection gives q = i sqrt(xi^2 - n^2).
/// @pre Re(n) > 0, Im(n) >= 0, xi >= 0, all finite
template <math::Real T>
[[nodiscard]] std::complex<T> normal_component(std::complex<T> n, T xi) noexcept {
  assert(n.real() > T(0) && n.imag() >= T(0) && xi >= T(0));
  using std::sqrt;
  const std::complex<T> z = n * n - xi * xi;
  // A negative zero imaginary part would put std::sqrt on the lower side of its branch cut
  // (q = -i b for total internal reflection); Im(n^2 - xi^2) = 2 n' kappa >= 0 here.
  std::complex<T> q = sqrt(std::complex<T>(z.real(), z.imag() == T(0) ? T(0) : z.imag()));
  if (q.imag() < T(0) || (q.imag() == T(0) && q.real() < T(0))) q = -q;
  return q;
}

/// Fresnel amplitudes for the tangential invariant xi = Re(n_i) sin(theta_i) (canonical form,
/// valid for absorbing media on both sides; Byrnes, Eq. (6) in the q form of the file comment).
/// Never throws (Rule 3).
/// @param n_i complex index of the incident medium, kappa >= 0
/// @param n_t complex index of the medium behind the interface, kappa >= 0
/// @param xi  tangential invariant, 0 <= xi <= Re(n_i)
/// @pre Re(n_i) > 0, Re(n_t) > 0, Im >= 0, 0 <= xi <= Re(n_i), finite, and not xi = n_i = n_t
///      real (grazing incidence between equal media). Checked by fresnel_checked().
template <math::Real T>
[[nodiscard]] FresnelAmplitudes<T> fresnel(std::complex<T> n_i,
                                           std::complex<T> n_t,
                                           T xi) noexcept {
  assert(xi <= n_i.real());
  const std::complex<T> q_i = normal_component(n_i, xi);
  const std::complex<T> q_t = normal_component(n_t, xi);
  const std::complex<T> ni2 = n_i * n_i;
  const std::complex<T> nt2 = n_t * n_t;
  // Byrnes, Eq. (6) with cos_j = q_j / n_j; the p terms multiplied by n_i n_t.
  const std::complex<T> ds = q_i + q_t;
  const std::complex<T> dp = nt2 * q_i + ni2 * q_t;
  assert(ds != std::complex<T>(0) && dp != std::complex<T>(0));
  return {(q_i - q_t) / ds, (nt2 * q_i - ni2 * q_t) / dp, T(2) * q_i / ds,
          T(2) * n_i * n_t * q_i / dp};
}

/// Convenience form of fresnel() for a NON-absorbing incident medium: xi = n_i sin(theta_i)
/// from cos(theta_i).
/// @param n_i   real index of the incident medium (Im(n_i) = 0)
/// @param n_t   complex index of the medium behind the interface, kappa >= 0
/// @param cos_i cosine of the angle of incidence, in [0, 1]
/// @pre Im(n_i) = 0, 0 <= cos_i <= 1, and the preconditions of fresnel()
template <math::Real T>
[[nodiscard]] FresnelAmplitudes<T> fresnel_at_angle(std::complex<T> n_i,
                                                    std::complex<T> n_t,
                                                    T cos_i) noexcept {
  using std::sqrt;
  assert(n_i.imag() == T(0) && cos_i >= T(0) && cos_i <= T(1));
  return fresnel(n_i, n_t, n_i.real() * sqrt(T(1) - cos_i * cos_i));
}

/// Reflectances and transmittances (Byrnes, Eqs. (21)-(23)). For an absorbing incident medium
/// R + T can differ from 1 without gain (Byrnes, Sec. 4.4.1, App. B); for lossless incident
/// media R + T = 1 for s and p, also if the second medium absorbs.
/// @pre as fresnel(), and Re(q_i) > 0 (fails only for xi = n_i real, grazing incidence)
template <math::Real T>
[[nodiscard]] FresnelPower<T> fresnel_power(std::complex<T> n_i,
                                            std::complex<T> n_t,
                                            T xi) noexcept {
  const FresnelAmplitudes<T> a = fresnel(n_i, n_t, xi);
  const std::complex<T> q_i = normal_component(n_i, xi);
  const std::complex<T> q_t = normal_component(n_t, xi);
  assert(q_i.real() > T(0));
  // Byrnes, Eq. (21): Re(n cos) = Re(q); Eq. (22): Re(n cos*) = Re(n q* / n*) = Re(n^2 q*) / |n|^2.
  const auto p_factor = [](std::complex<T> n, std::complex<T> q) {
    return std::real(n * n * std::conj(q)) / std::norm(n);
  };
  FresnelPower<T> p;
  p.reflectance_s = std::norm(a.rs);
  p.reflectance_p = std::norm(a.rp);
  p.transmittance_s = std::norm(a.ts) * q_t.real() / q_i.real();
  p.transmittance_p = std::norm(a.tp) * p_factor(n_t, q_t) / p_factor(n_i, q_i);
  return p;
}

/// Phase difference of reflection Delta = arg(-r_p / r_s) in rad, in (-pi, pi] (ellipsometric
/// Delta, Byrnes, Eq. (16)); 0 at normal incidence. Delta > 0 means p is delayed against s. Returns
/// 0 if r_s or r_p is 0. The retardance with its slow axis follows from the PRT matrix (#57).
template <math::Real T>
[[nodiscard]] T reflection_phase_difference(const FresnelAmplitudes<T>& a) noexcept {
  if (a.rs == std::complex<T>(0) || a.rp == std::complex<T>(0)) return T(0);
  const T phase = std::arg(-a.rp / a.rs);
  return phase <= -std::numbers::pi_v<T> ? phase + T(2) * std::numbers::pi_v<T> : phase;
}

/// Phase difference of transmission Delta = arg(t_p / t_s) in rad, in (-pi, pi]; Delta > 0 means p
/// is delayed. Returns 0 if t_s or t_p is 0.
template <math::Real T>
[[nodiscard]] T transmission_phase_difference(const FresnelAmplitudes<T>& a) noexcept {
  if (a.ts == std::complex<T>(0) || a.tp == std::complex<T>(0)) return T(0);
  const T phase = std::arg(a.tp / a.ts);
  return phase <= -std::numbers::pi_v<T> ? phase + T(2) * std::numbers::pi_v<T> : phase;
}

/// Diattenuation of reflection D = (R_s - R_p) / (R_s + R_p) in [-1, 1]; 0 if R_s + R_p = 0.
template <math::Real T>
[[nodiscard]] T reflection_diattenuation(const FresnelPower<T>& p) noexcept {
  const T sum = p.reflectance_s + p.reflectance_p;
  return sum > T(0) ? (p.reflectance_s - p.reflectance_p) / sum : T(0);
}

/// Diattenuation of transmission D = (T_s - T_p) / (T_s + T_p) in [-1, 1]; 0 if T_s + T_p = 0.
template <math::Real T>
[[nodiscard]] T transmission_diattenuation(const FresnelPower<T>& p) noexcept {
  const T sum = p.transmittance_s + p.transmittance_p;
  return sum > T(0) ? (p.transmittance_s - p.transmittance_p) / sum : T(0);
}

/// fresnel() with the preconditions checked, for calls at the API boundary.
/// @throws std::invalid_argument if an index or xi is not finite, Re(n) <= 0 or Im(n) < 0, xi < 0
///         or xi > Re(n_i), or a Fresnel denominator vanishes (xi = n_i = n_t real: grazing
///         incidence between equal media)
[[nodiscard]] FresnelAmplitudes<double> fresnel_checked(std::complex<double> n_i,
                                                        std::complex<double> n_t,
                                                        double xi);

/// fresnel_power() with the preconditions checked.
/// @throws std::invalid_argument as fresnel_checked(), and if Re(q_i) = 0, i.e. xi = n_i real
///         (grazing incidence, no incident power through the interface)
[[nodiscard]] FresnelPower<double> fresnel_power_checked(std::complex<double> n_i,
                                                         std::complex<double> n_t,
                                                         double xi);

}  // namespace rtt::polar
