#pragma once

/// @file interface.hpp
/// PRT matrix of an interface for the tracer (ADR 0021): Fresnel or any other amplitudes a_s,
/// a_p (e.g. of a coating stack) put into the basis of rtt/polar/prt.hpp and POWER-NORMALISED,
/// a' = a sqrt(c) with the power factors c of S. J. Byrnes, arXiv:1603.02720v5, Eqs. (21), (22)
/// (docs/quellen.md). Then |P E|^2 is the power fraction for every incident state E (|E| = 1);
/// phases, hence retardance, are unchanged. P is therefore not the field-amplitude PRT matrix of
/// Chipman/Lam: the field amplitude differs by the factors sqrt(c) of the interfaces passed
/// (docs/architecture.md, Strahl-Batch; relevant for the Jones pupil, M6).
///
/// Conventions as in rtt/polar/fresnel.hpp: tangential invariant xi = Re(n_i) sin(theta_i) with
/// sin(theta_i) = |k_in x N|; Convention A; fields ~ exp(i(k.r - omega t)).

#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>

#include "rtt/math/real.hpp"
#include "rtt/math/types.hpp"
#include "rtt/polar/fresnel.hpp"
#include "rtt/polar/prt.hpp"

namespace rtt::polar {

/// Power factors c_s = Re(q_t) / Re(q_i) and c_p = Re(n_t cos*_t) / Re(n_i cos*_i) of a
/// transmission from n_i into n_t (Byrnes, Eqs. (21), (22)); both > 0 for a transmitted wave
/// (Byrnes, App. D.2). Reflection has c = 1.
template <math::Real T>
struct PowerFactors {
  T s = T(1);
  T p = T(1);
};

/// Power factors of a transmission at the tangential invariant xi.
/// @pre preconditions of fresnel_power(), Re(q_t) > 0 (no total internal reflection)
template <math::Real T>
[[nodiscard]] PowerFactors<T> transmission_power_factors(std::complex<T> n_i,
                                                         std::complex<T> n_t,
                                                         T xi) noexcept {
  const std::complex<T> q_i = normal_component(n_i, xi);
  const std::complex<T> q_t = normal_component(n_t, xi);
  assert(q_i.real() > T(0));
  // Byrnes, Eq. (21): Re(n cos) = Re(q); Eq. (22): Re(n cos*) = Re(n^2 q*) / |n|^2, as in
  // fresnel_power().
  const auto p_factor = [](std::complex<T> n, std::complex<T> q) {
    return std::real(n * n * std::conj(q)) / std::norm(n);
  };
  return {q_t.real() / q_i.real(), p_factor(n_t, q_t) / p_factor(n_i, q_i)};
}

/// Tangential invariant xi = Re(n_i) sin(theta_i) with sin(theta_i) = |k_in x N| (accurate near
/// normal incidence).
/// @param k_in   unit incident direction, global
/// @param normal unit surface normal, global, either orientation
/// @param n_i    complex index of the incident medium
template <math::Real T>
[[nodiscard]] T tangential_invariant(const math::Vec3T<T>& k_in,
                                     const math::Vec3T<T>& normal,
                                     std::complex<T> n_i) noexcept {
  using std::min;
  // |k x N| <= 1 up to rounding; the clamp keeps xi <= Re(n_i) (precondition of fresnel()).
  return n_i.real() * min(k_in.cross(normal).norm(), T(1));
}

/// Power-normalised PRT matrix of an interface: prt_matrix(k_in, k_out, N, a_s sqrt(c_s),
/// a_p sqrt(c_p)).
/// @param a_s, a_p amplitudes in Convention A (basis p = k x s)
/// @param c        power factors (c = 1 for reflection)
template <math::Real T>
[[nodiscard]] CMat3T<T> interface_prt(const math::Vec3T<T>& k_in,
                                      const math::Vec3T<T>& k_out,
                                      const math::Vec3T<T>& normal,
                                      std::complex<T> a_s,
                                      std::complex<T> a_p,
                                      const PowerFactors<T>& c) noexcept {
  using std::sqrt;
  return prt_matrix(k_in, k_out, normal, a_s * sqrt(c.s), a_p * sqrt(c.p));
}

/// Power-normalised PRT matrix of a Fresnel interface between n_i and n_t (Byrnes, Eq. (6) via
/// rtt/polar/fresnel.hpp): reflection with r_s, r_p (c = 1), transmission with t_s, t_p and the
/// power factors. Snell's law with Re(n) determines k_out (the caller's refraction).
/// @param k_in       unit incident direction, global
/// @param k_out      unit reflected or refracted direction, global
/// @param normal     unit surface normal, global, either orientation
/// @param n_i        complex index of the incident medium
/// @param n_t        complex index on the other side of the surface
/// @param reflection true: reflected wave, false: transmitted wave
/// @pre transmission only without total internal reflection
template <math::Real T>
[[nodiscard]] CMat3T<T> fresnel_prt(const math::Vec3T<T>& k_in,
                                    const math::Vec3T<T>& k_out,
                                    const math::Vec3T<T>& normal,
                                    std::complex<T> n_i,
                                    std::complex<T> n_t,
                                    bool reflection) noexcept {
  const T xi = tangential_invariant(k_in, normal, n_i);
  const FresnelAmplitudes<T> a = fresnel(n_i, n_t, xi);
  if (reflection) return interface_prt(k_in, k_out, normal, a.rs, a.rp, PowerFactors<T>{});
  return interface_prt(k_in, k_out, normal, a.ts, a.tp, transmission_power_factors(n_i, n_t, xi));
}

}  // namespace rtt::polar
