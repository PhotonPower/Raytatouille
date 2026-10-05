#pragma once

/// @file ideal.hpp
/// PRT matrices of ideal thin polarizing elements: a general Jones matrix embedded in 3D, the
/// linear polarizer with extinction ratio, the linear and the circular retarder (#60;
/// docs/architecture.md, rtt-polar and Konventionen "Polarisation und Fresnel").
///
/// Source: W.-S. T. Lam, "Anisotropic ray trace", dissertation, University of Arizona
/// (docs/quellen.md): the matrices are the spectral forms of the definitions of diattenuation
/// (Eq. (4.3)) and retardance (Eq. (4.4), the eigenstate with the smaller phase is the fast
/// axis); handedness after Lam, Table 2.1 (p. 69) and Fig. 1.1 (p. 48); the half-wave plate of
/// Lam, p. 123 is reproduced up to a global phase; the embedding of a Jones matrix is Lam,
/// - Fields ~ exp(i(k.r - omega t)), so a larger phase is a delay; positive retardance delays the
///   slow axis.
/// - Thin elements do not deflect the ray: P k = k. The transverse basis is e1 (an axis projected
///   perpendicular to k and normalised) and e2 = k x e1, right-handed like [s, p, k].
/// - Right circular polarization is E = (i e1 + e2)/sqrt(2) (Lam, Table 2.1: (i, 1, 0) "evolving
///   clockwise in time as one looks into the beam"); left circular is (e1 + i e2)/sqrt(2).
/// - Retarders split the phase symmetrically, e^(-i delta/2) and e^(+i delta/2): an ideal retarder
///   then adds no piston to the wavefront, so OPD maps through ideal retarders are not shifted
///   by delta/2 (decided for #60).
///
/// Templates over rtt::math::Real (ADR 0006, 0015); the kernels never throw (Rule 3), the
/// *_checked variants check the preconditions at the API boundary.

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cassert>
#include <cmath>
#include <complex>
#include <numbers>

#include "rtt/math/real.hpp"
#include "rtt/math/types.hpp"
#include "rtt/polar/prt.hpp"

namespace rtt::polar {

/// 2x2 complex Jones matrix over a generic scalar, in the transverse basis (e1, e2).
template <math::Real T>
using Jones2T = Eigen::Matrix<std::complex<T>, 2, 2>;

/// Smallest admissible length of an element axis projected perpendicular to k (unit vectors);
/// below it the axis is (nearly) parallel to the ray and the element is undefined.
inline constexpr double kAxisAlongK = 1e-12;

/// Unit vector of `axis` projected perpendicular to the unit direction k.
/// @pre axis != 0 and |axis - (axis . k) k| >= kAxisAlongK |axis|
template <math::Real T>
[[nodiscard]] math::Vec3T<T> transverse_axis(const math::Vec3T<T>& axis,
                                             const math::Vec3T<T>& k) noexcept {
  const math::Vec3T<T> a = axis - axis.dot(k) * k;
  assert(axis.norm() > T(0) && a.norm() >= T(kAxisAlongK) * axis.norm());
  return a / a.norm();
}

/// PRT matrix of a thin element with Jones matrix J in the transverse basis e1 = axis projected
/// perpendicular to k, e2 = k x e1: P = sum_ij J_ij e_i e_j^T + k k^T (P k = k); this is Lam,
/// Eq. (3.4), P = O_out J_3D O_in^-1 with O = [e1 e2 k] and J_3D of Eq. (3.6).
/// @param j    2x2 Jones matrix, J(i, j) maps component j onto component i
/// @param axis direction of e1, global, any length, not parallel to k
/// @param k    unit propagation direction, global
template <math::Real T>
[[nodiscard]] CMat3T<T> embed_jones(const Jones2T<T>& j,
                                    const math::Vec3T<T>& axis,
                                    const math::Vec3T<T>& k) noexcept {
  using C = std::complex<T>;
  const math::Vec3T<T> e1 = transverse_axis(axis, k);
  const math::Vec3T<T> e2 = k.cross(e1);
  const auto outer = [](const math::Vec3T<T>& a, const math::Vec3T<T>& b) {
    return CMat3T<T>((a * b.transpose()).template cast<C>());
  };
  return j(0, 0) * outer(e1, e1) + j(0, 1) * outer(e1, e2) + j(1, 0) * outer(e2, e1) +
         j(1, 1) * outer(e2, e2) + outer(k, k);
}

/// Ideal linear polarizer: transmission axis a projected perpendicular to k with amplitude 1 (no
/// loss, no phase) and the crossed axis k x a with amplitude sqrt(epsilon), where epsilon =
/// T_min / T_max is the extinction ratio in POWER (model::IdealPolarizer). Diattenuation
/// D = (1 - epsilon) / (1 + epsilon).
/// @param axis             transmission axis, global, not parallel to k
/// @param extinction_ratio epsilon in [0, 1], 0 = perfect polarizer
/// @param k                unit propagation direction, global
template <math::Real T>
[[nodiscard]] CMat3T<T> ideal_polarizer(const math::Vec3T<T>& axis,
                                        T extinction_ratio,
                                        const math::Vec3T<T>& k) noexcept {
  using std::sqrt;
  assert(extinction_ratio >= T(0) && extinction_ratio <= T(1));
  Jones2T<T> j = Jones2T<T>::Zero();
  j(0, 0) = std::complex<T>(T(1));
  j(1, 1) = std::complex<T>(sqrt(extinction_ratio));
  return embed_jones(j, axis, k);
}

/// Ideal linear retarder with fast axis f (projected perpendicular to k) and slow axis k x f:
/// Jones matrix diag(e^(-i delta/2), e^(+i delta/2)) with delta = 2 pi retardance_waves (symmetric
/// phase, file comment). Positive retardance delays the slow axis (model::IdealRetarder); a
/// negative retardance exchanges the roles of fast and slow axis.
/// @param fast_axis        fast axis, global, not parallel to k
/// @param retardance_waves retardance in waves (delta / (2 pi))
/// @param k                unit propagation direction, global
template <math::Real T>
[[nodiscard]] CMat3T<T> linear_retarder(const math::Vec3T<T>& fast_axis,
                                        T retardance_waves,
                                        const math::Vec3T<T>& k) noexcept {
  using std::cos;
  using std::sin;
  const T half = std::numbers::pi_v<T> * retardance_waves;  // delta / 2
  Jones2T<T> j = Jones2T<T>::Zero();
  j(0, 0) = std::complex<T>(cos(half), -sin(half));  // fast: e^(-i delta/2)
  j(1, 1) = std::complex<T>(cos(half), sin(half));   // slow: e^(+i delta/2), delayed
  return embed_jones(j, fast_axis, k);
}

/// Ideal circular retarder: right circular polarization (file comment) is delayed by delta =
/// 2 pi retardance_waves against left circular, with the symmetric phase e^(+-i delta/2). This is
/// a rotation of linear polarization by +delta/2 from e1 towards e2 (counter-clockwise about k):
/// J = [[cos(delta/2), -sin(delta/2)], [sin(delta/2), cos(delta/2)]], independent of e1.
/// @param retardance_waves circular retardance in waves
/// @param k                unit propagation direction, global
template <math::Real T>
[[nodiscard]] CMat3T<T> circular_retarder(T retardance_waves, const math::Vec3T<T>& k) noexcept {
  using std::cos;
  using std::sin;
  const T half = std::numbers::pi_v<T> * retardance_waves;  // delta / 2
  Jones2T<T> j;
  j << std::complex<T>(cos(half)), std::complex<T>(-sin(half)), std::complex<T>(sin(half)),
      std::complex<T>(cos(half));
  // Any e1 transverse to k gives the same matrix (a rotation about k).
  const math::Vec3T<T> e1 = k.unitOrthogonal();
  return embed_jones(j, e1, k);
}

/// Stokes parameters of a field in the transverse basis (e1, e2), in units of |E|^2. Project
template <math::Real T>
struct Stokes {
  T s0 = T(0);  ///< |E1|^2 + |E2|^2
  T s1 = T(0);  ///< |E1|^2 - |E2|^2
  T s2 = T(0);  ///< 2 Re(E1* E2)
  T s3 = T(0);  ///< 2 Im(E1 E2*); > 0 for right circular (Lam, Table 2.1)
};

/// Stokes parameters of the complex field e (transverse to k) in the basis e1 = axis projected
/// perpendicular to k, e2 = k x e1. The linear polarization angle from e1 towards e2 is
/// atan2(s2, s1) / 2; s3 > 0 is right circular, "clockwise as one looks into the beam" (Lam,
/// Table 2.1, p. 69; docs/architecture.md).
template <math::Real T>
[[nodiscard]] Stokes<T> stokes(const Eigen::Matrix<std::complex<T>, 3, 1>& e,
                               const math::Vec3T<T>& axis,
                               const math::Vec3T<T>& k) noexcept {
  using C = std::complex<T>;
  const math::Vec3T<T> u1 = transverse_axis(axis, k);
  const math::Vec3T<T> u2 = k.cross(u1);
  const C e1 = u1.template cast<C>().dot(e);  // real basis: no conjugation effect
  const C e2 = u2.template cast<C>().dot(e);
  Stokes<T> s;
  s.s0 = std::norm(e1) + std::norm(e2);
  s.s1 = std::norm(e1) - std::norm(e2);
  s.s2 = T(2) * std::real(std::conj(e1) * e2);
  s.s3 = T(2) * std::imag(e1 * std::conj(e2));
  return s;
}

/// ideal_polarizer() with the preconditions checked, for calls at the API boundary.
/// @throws std::invalid_argument if a vector or the ratio is not finite, |k| is not 1 (to 1e-9),
///         the axis is parallel to k, or the extinction ratio is outside [0, 1]
[[nodiscard]] math::CMat3 ideal_polarizer_checked(const math::Vec3& axis,
                                                  double extinction_ratio,
                                                  const math::Vec3& k);

/// linear_retarder() with the preconditions checked.
/// @throws std::invalid_argument if a vector or the retardance is not finite, |k| is not 1 (to
///         1e-9), or the fast axis is parallel to k
[[nodiscard]] math::CMat3 linear_retarder_checked(const math::Vec3& fast_axis,
                                                  double retardance_waves,
                                                  const math::Vec3& k);

/// circular_retarder() with the preconditions checked.
/// @throws std::invalid_argument if the retardance or k is not finite or |k| is not 1 (to 1e-9)
[[nodiscard]] math::CMat3 circular_retarder_checked(double retardance_waves, const math::Vec3& k);

}  // namespace rtt::polar
