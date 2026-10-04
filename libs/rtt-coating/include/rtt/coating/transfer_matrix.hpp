#pragma once

/// @file transfer_matrix.hpp
/// Reflection and transmission amplitudes of planar thin-film stacks with complex indices, by
/// the transfer-matrix method in the interface formulation of S. J. Byrnes, "Multilayer optical
/// calculations", arXiv:1603.02720v5 (2020) (docs/quellen.md). Equation numbers below refer to
/// that paper.
///
/// Conventions (docs/architecture.md, "Konventionen"; decided for #56/#58):
/// - Fields ~ exp(i(k·r - omega t)); complex index n + i kappa with kappa >= 0 for absorption
///   (Byrnes Sec. 2, 2.1).
/// - Amplitudes are the Jones-matrix entries a_s, a_p of rtt-polar in the basis [s, p, k] with
///   s = k_in x n normalised and p = k x s. This is Byrnes' "Convention A" (Eq. (6), App. A,
///   as Jackson (7.41), Hecht (4.38)): r_p = (n2 cos1 - n1 cos2) / (n2 cos1 + n1 cos2), so
///   r_p = -r_s at normal incidence; an ideal conductor gives r_s = -1, r_p = +1.
/// - Tangential invariant xi = n sin(theta), real and the same in every medium (Eq. (3)); it is
///   the canonical angle input. For a non-absorbing incident medium xi = n_in sin(theta_in).
/// - Normal component q = n cos(theta) = sqrt(n^2 - xi^2) with Im q >= 0, and Re q >= 0 if
///   Im q = 0 (App. D: the forward wave decays or carries energy forwards; under total internal
///   reflection q is imaginary with Im q > 0).
/// - Layers are listed from the incident medium to the exit medium (substrate); thicknesses
///   and the vacuum wavelength in um.
///
/// Error handling (rule 3): the kernels are noexcept and never throw; their preconditions are
/// documented and asserted in debug builds. check_stack() validates inputs at the API boundary.

#include <cassert>
#include <cmath>
#include <complex>
#include <cstddef>
#include <numbers>
#include <span>

#include "rtt/math/real.hpp"

namespace rtt::coating {

/// One layer at one wavelength.
template <math::Real T>
struct Layer {
  std::complex<T> index;  ///< complex index n + i kappa (kappa >= 0 absorbs), dimensionless
  T thickness_um{0};      ///< physical thickness, um, >= 0
};

/// Complex amplitudes in the basis [s, p, k] of rtt-polar (Convention A, see file comment):
/// reflected / incident and transmitted / incident field amplitude.
template <math::Real T>
struct Amplitudes {
  std::complex<T> rs{0};
  std::complex<T> rp{0};
  std::complex<T> ts{0};
  std::complex<T> tp{0};
};

/// Reflected and transmitted power as fractions of the incident power (Eq. (21)-(23)).
template <math::Real T>
struct Powers {
  T reflectance_s{0};
  T reflectance_p{0};
  T transmittance_s{0};
  T transmittance_p{0};
  /// A = 1 - R - T for s. Absorption in the stack only if the incident medium does not absorb
  /// (App. B: with an absorbing incident medium R + T + A = 1 does not describe absorption).
  [[nodiscard]] T absorptance_s() const noexcept { return T(1) - reflectance_s - transmittance_s; }
  /// A = 1 - R - T for p, see absorptance_s().
  [[nodiscard]] T absorptance_p() const noexcept { return T(1) - reflectance_p - transmittance_p; }
};

/// Tangential invariant xi = n sin(theta) for an incident medium WITHOUT absorption.
/// @param n_in      real index of the incident medium
/// @param theta_rad angle of incidence from the surface normal, rad, in [0, pi/2)
template <math::Real T>
[[nodiscard]] T tangential_invariant(T n_in, T theta_rad) noexcept {
  using std::sin;
  return n_in * sin(theta_rad);
}

/// Normal component q = n cos(theta) = sqrt(n^2 - xi^2) in a medium of index n, with the
/// branch of App. D: Im q >= 0, and Re q >= 0 if Im q = 0.
/// @param n  complex index
/// @param xi tangential invariant (real)
template <math::Real T>
[[nodiscard]] std::complex<T> normal_component(std::complex<T> n, T xi) noexcept {
  using std::sqrt;
  std::complex<T> q = sqrt(n * n - std::complex<T>(xi * xi));
  if (q.imag() < T(0) || (q.imag() == T(0) && q.real() < T(0))) q = -q;
  return q;
}

/// Fresnel amplitudes of one interface from medium 1 into medium 2 (Eq. (6)), with
/// cos(theta_j) = q_j / n_j (q from normal_component()):
/// r_s = (n1 cos1 - n2 cos2) / (n1 cos1 + n2 cos2), r_p = (n2 cos1 - n1 cos2) / (n2 cos1 + n1
/// cos2), t_s = 2 n1 cos1 / (n1 cos1 + n2 cos2),           t_p = 2 n1 cos1 / (n2 cos1 + n1 cos2).
/// @pre n1, n2 != 0; the denominators are not 0 (not at grazing incidence)
template <math::Real T>
[[nodiscard]] Amplitudes<T> interface_amplitudes(std::complex<T> n1,
                                                 std::complex<T> n2,
                                                 T xi) noexcept {
  const std::complex<T> q1 = normal_component(n1, xi);
  const std::complex<T> q2 = normal_component(n2, xi);
  const std::complex<T> cos1 = q1 / n1;
  const std::complex<T> cos2 = q2 / n2;
  const std::complex<T> ds = n1 * cos1 + n2 * cos2;
  const std::complex<T> dp = n2 * cos1 + n1 * cos2;
  return {(n1 * cos1 - n2 * cos2) / ds, (n2 * cos1 - n1 * cos2) / dp, T(2) * n1 * cos1 / ds,
          T(2) * n1 * cos1 / dp};
}

namespace detail {

/// 2x2 complex matrix [[a, b], [c, d]].
template <math::Real T>
struct Matrix2 {
  std::complex<T> a, b, c, d;

  [[nodiscard]] Matrix2 operator*(const Matrix2& o) const noexcept {
    return {a * o.a + b * o.c, a * o.b + b * o.d, c * o.a + d * o.c, c * o.b + d * o.d};
  }
};

/// (1/t) [[1, r], [r, 1]] of one interface (factors of Eq. (11) and (13)).
template <math::Real T>
[[nodiscard]] Matrix2<T> interface_matrix(std::complex<T> r, std::complex<T> t) noexcept {
  const std::complex<T> inv = T(1) / t;
  return {inv, r * inv, r * inv, inv};
}

}  // namespace detail

/// Amplitudes r_s, r_p, t_s, t_p of a stack (Eq. (8), (11), (13), (15)):
/// delta_n = 2 pi q_n d_n / lambda (q = n cos theta), M_n = diag(e^{-i delta_n}, e^{i delta_n})
/// (1/t_{n,n+1}) [[1, r_{n,n+1}], [r_{n,n+1}, 1]], M~ = (1/t_01) [[1, r_01], [r_01, 1]]
/// M_1 ... M_{N-2}, t = 1 / M~_00, r = M~_10 / M~_00, separately for s and p with the interface
/// amplitudes of Eq. (6). Without layers the result is interface_amplitudes(n_in, n_out, xi).
/// @param n_in          complex index of the incident medium (semi-infinite)
/// @param layers        layers from the incident medium to the exit medium
/// @param n_out         complex index of the exit medium (substrate, semi-infinite)
/// @param xi            tangential invariant n sin(theta), real
/// @param wavelength_um vacuum wavelength, um
/// @pre wavelength_um > 0; thicknesses >= 0; indices != 0; no medium with q = 0 exactly
///      (grazing in the incident medium or at the critical angle of a layer: the interface
///      formulation divides by t = 0 there and the result is not finite)
template <math::Real T>
[[nodiscard]] Amplitudes<T> stack_amplitudes(std::complex<T> n_in,
                                             std::span<const Layer<T>> layers,
                                             std::complex<T> n_out,
                                             T xi,
                                             T wavelength_um) noexcept {
  assert(wavelength_um > T(0));
  using std::exp;
  const auto next_index = [&](std::size_t j) {
    return j + 1 < layers.size() ? layers[j + 1].index : n_out;
  };
  const Amplitudes<T> first =
      interface_amplitudes(n_in, layers.empty() ? n_out : layers[0].index, xi);
  detail::Matrix2<T> ms = detail::interface_matrix(first.rs, first.ts);
  detail::Matrix2<T> mp = detail::interface_matrix(first.rp, first.tp);
  const std::complex<T> i(T(0), T(1));
  for (std::size_t j = 0; j < layers.size(); ++j) {
    assert(layers[j].thickness_um >= T(0));
    const std::complex<T> q = normal_component(layers[j].index, xi);
    const std::complex<T> delta =
        T(2) * std::numbers::pi_v<T> * q * layers[j].thickness_um / wavelength_um;  // Eq. (8)
    const std::complex<T> forward = exp(-i * delta);
    const std::complex<T> backward = exp(i * delta);
    const detail::Matrix2<T> phase{forward, std::complex<T>(0), std::complex<T>(0), backward};
    const Amplitudes<T> next = interface_amplitudes(layers[j].index, next_index(j), xi);
    ms = ms * (phase * detail::interface_matrix(next.rs, next.ts));  // Eq. (11), (13)
    mp = mp * (phase * detail::interface_matrix(next.rp, next.tp));
  }
  return {ms.c / ms.a, mp.c / mp.a, T(1) / ms.a, T(1) / mp.a};  // Eq. (15)
}

/// Reflected and transmitted power for amplitudes of stack_amplitudes() (Eq. (21)-(23)):
/// R = |r|^2, T_s = |t_s|^2 Re[n_out cos_out] / Re[n_in cos_in],
/// T_p = |t_p|^2 Re[n_out conj(cos_out)] / Re[n_in conj(cos_in)].
/// Without absorption in the stack R + T = 1. With an absorbing incident medium the values
/// follow the definitions of App. B (R + T may exceed 1).
/// @pre not at grazing incidence (Re[n_in cos_in] != 0)
template <math::Real T>
[[nodiscard]] Powers<T> stack_powers(const Amplitudes<T>& a,
                                     std::complex<T> n_in,
                                     std::complex<T> n_out,
                                     T xi) noexcept {
  using std::conj;
  using std::norm;
  const std::complex<T> q_in = normal_component(n_in, xi);
  const std::complex<T> q_out = normal_component(n_out, xi);
  const std::complex<T> cos_in = q_in / n_in;
  const std::complex<T> cos_out = q_out / n_out;
  return {norm(a.rs), norm(a.rp), norm(a.ts) * q_out.real() / q_in.real(),
          norm(a.tp) * (n_out * conj(cos_out)).real() / (n_in * conj(cos_in)).real()};
}

/// Validates the inputs of stack_amplitudes() at the API boundary.
/// @throws std::invalid_argument if the wavelength is not finite and > 0, a thickness is not
///         finite and >= 0, an index is not finite or 0, an index has kappa < 0 (gain), xi is
///         not finite or < 0, or xi >= n_in for a non-absorbing incident medium (grazing or
///         beyond)
void check_stack(std::complex<double> n_in,
                 std::span<const Layer<double>> layers,
                 std::complex<double> n_out,
                 double xi,
                 double wavelength_um);

}  // namespace rtt::coating
