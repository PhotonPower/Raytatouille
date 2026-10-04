#pragma once

/// @file prt_analysis.hpp
/// Diattenuation and retardance of a PRT matrix (Lam, "Anisotropic ray trace", Eqs. (4.3),
/// (4.4) and Sec. 4.5.2; docs/quellen.md). Double only: these evaluations use singular value and
/// eigen decompositions and are analysis outputs, not part of the derivative path of the
/// optimisation (decided for #57). Conventions as in rtt/polar/prt.hpp; fields ~ exp(i(k.r -
/// omega t)), so a larger phase is a delay and positive retardance means the slow axis is
/// delayed (docs/architecture.md).

#include "rtt/math/types.hpp"

namespace rtt::polar {

/// Diattenuation of a PRT matrix.
struct Diattenuation {
  double value = 0.0;    ///< D = (L1^2 - L2^2) / (L1^2 + L2^2) in [0, 1]
  double maximum = 0.0;  ///< L1: largest transverse amplitude transmission (singular value)
  double minimum = 0.0;  ///< L2: smallest transverse amplitude transmission
  /// Incident polarization (unit, transverse to k_in, complex) with the largest transmission:
  /// right singular vector of L1. Arbitrary within the transverse plane if D = 0.
  Eigen::Vector3cd axis = Eigen::Vector3cd::Zero();
};

/// Retardance of a PRT matrix that maps a direction onto itself.
struct Retardance {
  double value = 0.0;  ///< delta in [0, pi], rad: phase difference of the eigenpolarizations
  /// Eigenpolarization with the smaller phase (fast axis), unit, complex, transverse to k. For
  /// delta = 0 or delta = pi it is not unique; the solver order then decides deterministically.
  Eigen::Vector3cd fast_axis = Eigen::Vector3cd::Zero();
};

/// Diattenuation from the singular values L1 >= L2 of the transverse part
/// P_T = P - k_out k_in^T (rank 2) (Lam, Eq. (4.3)). The transverse part is used because
/// transmission amplitudes can exceed 1 (glass -> air), so the singular value 1 of k could
/// otherwise be confused with a transverse one. D is based on FIELD amplitudes: it equals the
/// power diattenuation of rtt/polar/fresnel.hpp for every reflection (R = |r|^2) and for
/// transmission between lossless media, and differs only for transmission with absorption
/// (different power factors for s and p). Never throws.
/// @param p     PRT matrix with p k_in = k_out
/// @param k_in  unit incident direction, global
/// @param k_out unit exiting direction, global
[[nodiscard]] Diattenuation diattenuation(const math::CMat3& p,
                                          const math::Vec3& k_in,
                                          const math::Vec3& k_out) noexcept;

/// Retardance of a matrix m that maps k onto itself (m k = k), e.g. Q^-1 P (physical
/// retardance) or P of a path that leaves in its incident direction (Lam, Sec. 4.5.2): polar
/// decomposition m = M_R M_D from the SVD m = U D V^H with M_R = U V^H and eigen decomposition of
/// the unitary M_R, done on the 2x2 block of m in an orthonormal basis of the plane
/// perpendicular to k (the eigenvalue 1 of k is left out, and the fast axis is always transverse
/// even if a transverse eigenvalue equals 1); of the two eigenvalues, delta = |arg(lambda_a /
/// lambda_b)| in [0, pi] (Lam, Eq. (4.4)) and the fast axis is the eigenvector with the smaller
/// phase (Lam, Sec. 4.5). Unlike Lam (e.g. delta = 3.519 rad on p. 165), delta is reduced to
/// [0, pi], which may swap fast and slow axis. Never throws.
/// @pre m k = k with |k| = 1, m maps the plane perpendicular to k onto itself, and m is
///      invertible there (no ideal polarizer); otherwise M_R is not unique
[[nodiscard]] Retardance retardance(const math::CMat3& m, const math::Vec3& k) noexcept;

/// Physical retardance of an intercept or path: retardance of Q^-1 P (Lam, Sec. 4.5.1), which
/// maps k_in onto itself. Q is real orthogonal, so Q^-1 = Q^T.
/// @param p    PRT matrix (P k_in = k_out)
/// @param q    geometric transformation of the same intercept or path (Q k_in = k_out)
/// @param k_in unit incident direction, global
[[nodiscard]] Retardance physical_retardance(const math::CMat3& p,
                                             const math::Mat3& q,
                                             const math::Vec3& k_in) noexcept;

}  // namespace rtt::polar
