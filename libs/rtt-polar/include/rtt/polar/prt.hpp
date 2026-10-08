#pragma once

/// @file prt.hpp
/// Three-dimensional polarization ray tracing (PRT) matrix of one ray intercept and its
/// geometric transformation (docs/architecture.md, rtt-polar and Konventionen "Polarisation
/// und Fresnel").
///
/// Source: W.-S. T. Lam, "Anisotropic ray trace", dissertation, University of Arizona (Chipman
/// group), Eqs. (3.1)-(3.9) and (4.5)-(4.6); see docs/quellen.md. The calculus goes back to Yun,
/// McClain, Chipman, "Three-dimensional polarization ray-tracing calculus I/II", Applied Optics
/// 2011 (not freely available; its content is checked via Lam). Conventions:
/// - All vectors are unit vectors in global coordinates (right-handed, optical axis +z); fields
///   ~ exp(i(k.r - omega t)).
/// - Basis of an isotropic intercept (Lam, Eq. (3.3)): s = k_in x N / |k_in x N|, p_in =
///   k_in x s, p_out = k_out x s, the SAME s for the incident and the exiting wave (s' = s for
///   an isotropic interface). Reversing N reverses s and both p together; P does not change.
/// - P = O_out J O_in^T with O_in = [s p_in k_in], O_out = [s p_out k_out] (columns) and
///   J = diag(a_s, a_p, 1) (Lam, Eqs. (3.4)-(3.8)), i.e.
///   P = a_s s s^T + a_p p_out p_in^T + k_out k_in^T. Then P k_in = k_out and E' = P E. With
///   Fresnel (rtt/polar/fresnel.hpp) a_s, a_p are r_s, r_p or t_s, t_p in Convention A.
/// - Normal incidence (|k_in x N| < kNormalIncidence): s = k_in x e / |k_in x e| with e the
///   global axis x, y or z with the smallest |k_in . e| (ties in the order x, y, z). This is
///   deterministic and stateless. P does not depend on this choice whenever J is isotropic in
///   the transverse plane, which holds for Fresnel and coatings at normal incidence: reflection
///   has a_p = -a_s and p_out = -p_in, so P = a_s (I - k k^T) - k k^T; transmission has
///   a_p = a_s, so P = a_s (I - k k^T) + k k^T.
/// - Accumulation along a ray path: P_total = P_N ... P_2 P_1 (Lam, Eq. (3.9)).
/// - Geometric transformation (Lam, Eqs. (4.5), (4.6)): Q = O_out diag(1, +-1, 1) O_in^T with
///   +1 for refraction and -1 for reflection. Lam prints the (3,3) element of diag(...) as 0, which
///   contradicts P k = k' (Eq. (3.2)) and the use of Q^-1; it is 1 here (docs/quellen.md).
///
/// The constructions are templates over rtt::math::Real (ADR 0006, 0015). The evaluations by
/// singular value and eigen decomposition (prt_analysis.hpp) are double only: they are analysis
/// outputs, not part of the derivative path of the optimisation.

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cmath>
#include <complex>

#include "rtt/math/real.hpp"
#include "rtt/math/types.hpp"

namespace rtt::polar {

/// 3x3 complex matrix over a generic scalar (PRT matrices).
template <math::Real T>
using CMat3T = Eigen::Matrix<std::complex<T>, 3, 3>;

/// Below this value of |k_in x N| (both unit vectors) the incidence is treated as normal and s
/// follows the axis rule of the file comment (decided for #57).
inline constexpr double kNormalIncidence = 1e-12;

namespace detail {

/// Outer product u v^T, evaluated element by element (m_ij = u_i v_j). Eigen's lazy outer product
/// gives the same values, but GCC -O2 reports a false -Wnull-dereference inside Eigen for every
/// translation unit that instantiates it (#35).
template <math::Real T>
[[nodiscard]] math::Mat3T<T> outer(const math::Vec3T<T>& u, const math::Vec3T<T>& v) noexcept {
  math::Mat3T<T> m;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) m(i, j) = u(i) * v(j);
  }
  return m;
}

}  // namespace detail

/// Right-handed basis (s, p, k) of a ray at an intercept, unit vectors in global coordinates.
template <math::Real T>
struct PrtBasis {
  math::Vec3T<T> s;  ///< perpendicular to the plane of incidence
  math::Vec3T<T> p;  ///< p = k x s
  math::Vec3T<T> k;  ///< propagation direction
};

/// Basis (s, p = k x s, k) of a ray with direction k at a surface with normal N (Lam, Eq. (3.3)),
/// with the normal-incidence rule of the file comment.
/// @param k      unit propagation direction, global
/// @param normal unit surface normal, global, either orientation
template <math::Real T>
[[nodiscard]] PrtBasis<T> prt_basis(const math::Vec3T<T>& k,
                                    const math::Vec3T<T>& normal) noexcept {
  using std::abs;
  math::Vec3T<T> c = k.cross(normal);
  T length = c.norm();
  if (!(length >= T(kNormalIncidence))) {
    int axis = 0;
    T best = abs(k.x());
    if (abs(k.y()) < best) {
      axis = 1;
      best = abs(k.y());
    }
    if (abs(k.z()) < best) axis = 2;
    c = k.cross(math::Vec3T<T>::Unit(axis));
  }
  // Just above kNormalIncidence the rounding of k x N (~1e-16) is large relative to |k x N|, so
  // c is made exactly transverse to k before normalising; its direction in the transverse plane
  // may then be inaccurate, which does not matter there (file comment).
  c -= c.dot(k) * k;
  length = c.norm();
  const math::Vec3T<T> s = c / length;
  return {s, k.cross(s), k};
}

/// PRT matrix P = a_s s s^T + a_p p_out p_in^T + k_out k_in^T of one intercept (Lam, Eq. (3.4)).
/// @param k_in  unit incident direction, global
/// @param k_out unit exiting direction (refracted or reflected), global
/// @param normal unit surface normal, global, either orientation
/// @param a_s   complex amplitude coefficient for s (e.g. r_s or t_s), dimensionless
/// @param a_p   complex amplitude coefficient for p in the basis p = k x s (Convention A)
/// @pre unit vectors; k_out in the plane of k_in and N (isotropic interface)
template <math::Real T>
[[nodiscard]] CMat3T<T> prt_matrix(const math::Vec3T<T>& k_in,
                                   const math::Vec3T<T>& k_out,
                                   const math::Vec3T<T>& normal,
                                   std::complex<T> a_s,
                                   std::complex<T> a_p) noexcept {
  const PrtBasis<T> in = prt_basis(k_in, normal);
  const math::Vec3T<T> p_out = k_out.cross(in.s);
  const CMat3T<T> ss = detail::outer(in.s, in.s).template cast<std::complex<T>>();
  const CMat3T<T> pp = detail::outer(p_out, in.p).template cast<std::complex<T>>();
  const CMat3T<T> kk = detail::outer(k_out, k_in).template cast<std::complex<T>>();
  return a_s * ss + a_p * pp + kk;
}

/// Geometric transformation Q (non-polarizing P matrix) of one intercept (Lam, Eqs. (4.5),
/// (4.6), with the (3,3) element 1): Q = s s^T +- p_out p_in^T + k_out k_in^T, real orthogonal,
/// det Q = +1 for refraction and -1 for reflection. Q^-1 P carries the physical retardance.
/// @param k_in       unit incident direction, global
/// @param k_out      unit exiting direction, global
/// @param normal     unit surface normal, global, either orientation
/// @param reflection true for a reflected k_out, false for a refracted one
template <math::Real T>
[[nodiscard]] math::Mat3T<T> geometric_transform(const math::Vec3T<T>& k_in,
                                                 const math::Vec3T<T>& k_out,
                                                 const math::Vec3T<T>& normal,
                                                 bool reflection) noexcept {
  const PrtBasis<T> in = prt_basis(k_in, normal);
  const math::Vec3T<T> p_out = k_out.cross(in.s);
  const T sign = reflection ? T(-1) : T(1);
  const math::Vec3T<T> signed_p_out = sign * p_out;
  return detail::outer(in.s, in.s) + detail::outer(signed_p_out, in.p) + detail::outer(k_out, k_in);
}

}  // namespace rtt::polar
