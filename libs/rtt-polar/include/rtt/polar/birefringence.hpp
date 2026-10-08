#pragma once

/// @file birefringence.hpp
/// Refraction into and out of a uniaxial medium and the PRT matrices of these intercepts in the
/// projection model of M4 (ADR 0026, points 4 and 5).
///
/// Source: W.-S. T. Lam, "Anisotropic Ray Trace", dissertation, University of Arizona
/// (docs/quellen.md):
/// - Eq. (2.11): tangential phase matching; n k_inc,par is kept, the normal part follows.
/// - Eq. (2.39): 1/n_e(theta)^2 = cos^2(theta)/n_O^2 + sin^2(theta)/n_E^2, theta = angle(k, axis).
/// - Eq. (2.41): K-surface of the e-mode, k1^2/n_E^2 + k2^2/n_E^2 + k3^2/n_O^2 = 1 (k in units
///   of k_0, axis along the third coordinate); here in vector form
///   |k - (k.a) a|^2 / n_E^2 + (k.a)^2 / n_O^2 = 1.
/// - p. 107: S is normal to the K-surface at k; k_o is parallel to S_o; E_e is perpendicular to
///   S_e. p. 104: D_e lies in the plane of k and the axis, D_o is perpendicular to D_e and k.
/// - Eqs. (3.23), (3.24), (3.27): P at isotropic/anisotropic intercepts; here with the coupling
///   factors replaced by projections (ADR 0026, point 5).
/// Conventions: all vectors are global unit vectors (right-handed, optical axis +z), except the
/// tangential parts n t_par, which are dimensionless with length up to n; indices are
/// real (absorbing crystals are rejected by rtt-compile); the sign of the optic axis has no
/// meaning. Functions are noexcept; "no real solution" is std::nullopt, never an exception.

#include <Eigen/Geometry>
#include <cmath>
#include <complex>
#include <cstdint>
#include <optional>

#include "rtt/math/real.hpp"
#include "rtt/math/types.hpp"
#include "rtt/polar/prt.hpp"

namespace rtt::polar {

/// Eigenmode of a uniaxial medium.
enum class Mode : std::uint8_t {
  Ordinary,       ///< polarized perpendicular to the optic axis, index n_O, S parallel to k
  Extraordinary,  ///< polarized in the plane of k and the axis, index n_e(theta), walk-off
};

/// A propagating eigenmode after refraction into a uniaxial medium.
template <math::Real T>
struct ModeSolution {
  math::Vec3T<T> k;  ///< unit wave normal (phase direction), global
  math::Vec3T<T> s;  ///< unit Poynting vector (energy direction, the ray), global
  T n{};             ///< refractive index of the mode along k, dimensionless
  math::Vec3T<T> e;  ///< unit eigenpolarization E, global, perpendicular to s
};

namespace detail {

/// `tangential` projected into the plane perpendicular to the unit `normal`.
template <math::Real T>
[[nodiscard]] math::Vec3T<T> tangential_part(const math::Vec3T<T>& tangential,
                                             const math::Vec3T<T>& normal) noexcept {
  return tangential - tangential.dot(normal) * normal;
}

}  // namespace detail

/// Unit eigenpolarization E of a mode with wave normal `k` and Poynting vector `s` (Lam, p. 104
/// and p. 107):
/// - E_o = k x a / |k x a| with the axis a oriented so that a . k >= 0 (the sign of the axis has
///   no meaning; this only fixes the intermediate sign, ADR 0026, point 5);
/// - E_e = s x E_o, exact: E_e lies in the plane of k and a and is perpendicular to s.
/// If k is parallel to the axis (|k x a| < kNormalIncidence) both modes coincide; E_o then uses
/// the axis rule of prt.hpp (k x e with e the global axis of smallest |k . e|, ties x, y, z).
/// @param mode which eigenmode
/// @param k    unit wave normal of the mode, global
/// @param s    unit Poynting vector of the mode, global (= k for the o-mode)
/// @param axis optic axis, global, non-zero (normalised here)
/// @return unit electric field direction of the mode, global, perpendicular to s
template <math::Real T>
[[nodiscard]] math::Vec3T<T> eigen_polarization(Mode mode,
                                                const math::Vec3T<T>& k,
                                                const math::Vec3T<T>& s,
                                                const math::Vec3T<T>& axis) noexcept {
  math::Vec3T<T> a = axis.normalized();
  if (a.dot(k) < T(0)) a = -a;
  math::Vec3T<T> c = k.cross(a);
  if (!(c.norm() >= T(kNormalIncidence))) c = prt_basis(k, k).s;  // axis rule of prt_basis
  c -= c.dot(k) * k;
  const math::Vec3T<T> e_o = c.normalized();
  if (mode == Mode::Ordinary) return e_o;
  return s.cross(e_o).normalized();
}

/// Refraction into a uniaxial medium (Lam, Eqs. (2.11), (2.41) and p. 107): the wave vector
/// k = kappa + beta N keeps the tangential part kappa of `tangential`; for the o-mode |k| = n_O,
/// for the e-mode k lies on the K-surface |k - (k.a) a|^2 / n_E^2 + (k.a)^2 / n_O^2 = 1, a
/// quadratic equation in beta. S is the normal of the K-surface, S ~ (k - (k.a) a)/n_E^2 +
/// (k.a) a / n_O^2; of the two roots the one with S . N > 0 (energy into the medium) is taken.
/// @param tangential n t_par of the incident ray (n its index, t its unit direction), plus a
///                   diffraction term m lambda_0/(2 pi) g_par if any (ADR 0025); only its part
///                   perpendicular to `normal` is used. Dimensionless.
/// @param normal     unit surface normal pointing into the uniaxial medium, global
/// @param n_o        ordinary principal index n_O > 0
/// @param n_e        extraordinary principal index n_E > 0
/// @param axis       optic axis, global, non-zero; its sign has no meaning
/// @return the mode, or std::nullopt if it does not propagate (no real solution, or grazing with
///         S . N = 0; ADR 0025, point 7)
template <math::Real T>
[[nodiscard]] std::optional<ModeSolution<T>> uniaxial_mode(Mode mode,
                                                           const math::Vec3T<T>& tangential,
                                                           const math::Vec3T<T>& normal,
                                                           T n_o,
                                                           T n_e,
                                                           const math::Vec3T<T>& axis) noexcept {
  using std::sqrt;
  const math::Vec3T<T> kappa = detail::tangential_part(tangential, normal);
  const math::Vec3T<T> a = axis.normalized();
  ModeSolution<T> m;
  if (mode == Mode::Ordinary) {
    const T beta2 = n_o * n_o - kappa.squaredNorm();
    if (!(beta2 > T(0))) return std::nullopt;
    const math::Vec3T<T> k = kappa + sqrt(beta2) * normal;
    m.n = n_o;
    m.k = k / n_o;
    m.s = m.k;
  } else {
    // A |k|^2 + B (k.a)^2 = 1 with A = 1/n_E^2, B = 1/n_O^2 - 1/n_E^2 and k = kappa + beta N:
    // a2 beta^2 + 2 a1 beta + a0 = 0. For the root beta, S . N ~ a2 beta + a1 = +-sqrt(disc), so
    // the + root is the one with energy into the medium.
    const T big_a = T(1) / (n_e * n_e);
    const T big_b = T(1) / (n_o * n_o) - big_a;
    const T c = normal.dot(a);
    const T d = kappa.dot(a);
    const T a2 = big_a + big_b * c * c;
    const T a1 = big_b * d * c;
    const T a0 = big_a * kappa.squaredNorm() + big_b * d * d - T(1);
    const T disc = a1 * a1 - a2 * a0;
    if (!(disc > T(0))) return std::nullopt;
    const T beta = (-a1 + sqrt(disc)) / a2;
    const math::Vec3T<T> k = kappa + beta * normal;
    m.n = k.norm();
    m.k = k / m.n;
    m.s = (big_a * k + big_b * k.dot(a) * a).normalized();
  }
  m.e = eigen_polarization(mode, m.k, m.s, axis);
  return m;
}

/// Refraction out of a uniaxial medium into an isotropic one (Lam, Eq. (2.11)): the unit
/// direction t' with n t'_par = the tangential part of `tangential` and t' . N > 0.
/// @param tangential n k_par of the mode (mode index times wave normal), plus a diffraction term
///                   if any; only its part perpendicular to `normal` is used. Dimensionless.
/// @param normal     unit surface normal pointing into the isotropic medium, global
/// @param n          real index of the isotropic medium, > 0
/// @return the unit direction, or std::nullopt for total internal reflection or grazing exit
template <math::Real T>
[[nodiscard]] std::optional<math::Vec3T<T>> isotropic_from_tangential(
    const math::Vec3T<T>& tangential, const math::Vec3T<T>& normal, T n) noexcept {
  using std::sqrt;
  const math::Vec3T<T> kappa = detail::tangential_part(tangential, normal);
  const T beta2 = n * n - kappa.squaredNorm();
  if (!(beta2 > T(0))) return std::nullopt;
  return math::Vec3T<T>((kappa + sqrt(beta2) * normal) / n);
}

/// PRT matrix of the entry into mode `m` (ADR 0026, point 5; structure of Lam, Eq. (3.23)):
/// P = E_v e_v^T + S_v S_in^T, with e_v the unit projection of E_v perpendicular to S_in, so that
/// |P E| = |e_v . E| for a field E perpendicular to S_in and P S_in = S_v. If the projection is
/// shorter than kNormalIncidence, the mode does not couple (P = S_v S_in^T).
/// @param s_in unit direction of the incident ray (isotropic medium: S = k), global
/// @param m    the mode in the crystal (uniaxial_mode())
/// @return the power-normalised PRT matrix of the entry, dimensionless
template <math::Real T>
[[nodiscard]] CMat3T<T> prt_crystal_entry(const math::Vec3T<T>& s_in,
                                          const ModeSolution<T>& m) noexcept {
  CMat3T<T> p = detail::outer(m.s, s_in).template cast<std::complex<T>>();
  const math::Vec3T<T> projected = m.e - m.e.dot(s_in) * s_in;
  const T length = projected.norm();
  if (length >= T(kNormalIncidence)) {
    p += detail::outer(m.e, math::Vec3T<T>(projected / length)).template cast<std::complex<T>>();
  }
  return p;
}

/// PRT matrix of the exit from a mode into an isotropic medium (ADR 0026, point 5; structure of
/// Lam, Eqs. (3.24), (3.27)): P = e_out E_m^T + S_out S_m^T, with e_out the unit projection of
/// E_m perpendicular to S_out. The absent mode S_m x E_m carries no power, P S_m = S_out.
/// @param s_m   unit Poynting vector of the mode in the crystal, global
/// @param e_m   unit eigenpolarization of the mode (perpendicular to s_m), global
/// @param s_out unit direction in the isotropic medium (isotropic_from_tangential()), global
/// @return the power-normalised PRT matrix of the exit, dimensionless
template <math::Real T>
[[nodiscard]] CMat3T<T> prt_crystal_exit(const math::Vec3T<T>& s_m,
                                         const math::Vec3T<T>& e_m,
                                         const math::Vec3T<T>& s_out) noexcept {
  CMat3T<T> p = detail::outer(s_out, s_m).template cast<std::complex<T>>();
  const math::Vec3T<T> projected = e_m - e_m.dot(s_out) * s_out;
  const T length = projected.norm();
  if (length >= T(kNormalIncidence)) {
    p += detail::outer(math::Vec3T<T>(projected / length), e_m).template cast<std::complex<T>>();
  }
  return p;
}

}  // namespace rtt::polar
