// GCC (13 in CI, 16 locally) reports a false -Wnull-dereference inside Eigen's SSE code
// (emmintrin.h, BinaryFunctors.h) for the fixed-size SVD and eigen solvers at -O2/-O3 in Release.
// The warning is attributed to the Eigen headers, so the suppression has to cover the includes;
// it is limited to this translation unit and to GCC (cf. the GCC 13 workaround in agf.cpp, #24).
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wnull-dereference"
#endif

#include "rtt/polar/prt_analysis.hpp"

#include <Eigen/Eigenvalues>
#include <Eigen/SVD>
#include <cmath>
#include <complex>

namespace rtt::polar {

Diattenuation diattenuation(const math::CMat3& p,
                            const math::Vec3& k_in,
                            const math::Vec3& k_out) noexcept {
  const math::CMat3 transverse = p - (k_out * k_in.transpose()).cast<math::Complex>();
  const Eigen::JacobiSVD<math::CMat3> svd(transverse, Eigen::ComputeFullU | Eigen::ComputeFullV);
  const Eigen::Vector3d sigma = svd.singularValues();  // descending
  Diattenuation d;
  d.maximum = sigma(0);
  d.minimum = sigma(1);
  const double a = sigma(0) * sigma(0);
  const double b = sigma(1) * sigma(1);
  d.value = a + b > 0.0 ? (a - b) / (a + b) : 0.0;
  d.axis = svd.matrixV().col(0);
  return d;
}

Retardance retardance(const math::CMat3& m, const math::Vec3& k) noexcept {
  // m maps k onto k and the transverse plane onto itself (precondition), so its polar
  // decomposition is block diagonal and the transverse 2x2 block carries the retardance. Working
  // in that block keeps the eigenvectors transverse even when a transverse eigenvalue equals the
  // eigenvalue 1 of k (degenerate eigenspace, e.g. lossless refraction or a retarder with phase 0
  // on one axis).
  const math::Vec3 u1 = k.unitOrthogonal();
  const math::Vec3 u2 = k.cross(u1);
  Eigen::Matrix<math::Complex, 3, 2> basis;
  basis.col(0) = u1.cast<math::Complex>();
  basis.col(1) = u2.cast<math::Complex>();
  const Eigen::Matrix2cd block = basis.adjoint() * m * basis;
  const Eigen::JacobiSVD<Eigen::Matrix2cd> svd(block, Eigen::ComputeFullU | Eigen::ComputeFullV);
  const Eigen::Matrix2cd unitary = svd.matrixU() * svd.matrixV().adjoint();  // M_R
  const Eigen::ComplexEigenSolver<Eigen::Matrix2cd> eigen(unitary);
  const auto& values = eigen.eigenvalues();
  const auto& vectors = eigen.eigenvectors();
  // arg(lambda_0 / lambda_1) in [-pi, pi]; a larger phase is a delay (slow axis).
  const double difference = std::arg(values(0) / values(1));
  Retardance r;
  r.value = std::abs(difference);
  r.fast_axis = (basis * (difference > 0.0 ? vectors.col(1) : vectors.col(0))).normalized();
  return r;
}

Retardance physical_retardance(const math::CMat3& p,
                               const math::Mat3& q,
                               const math::Vec3& k_in) noexcept {
  return retardance(q.transpose().cast<math::Complex>() * p, k_in);
}

}  // namespace rtt::polar

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif
