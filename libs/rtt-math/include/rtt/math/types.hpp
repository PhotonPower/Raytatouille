#pragma once

/// @file types.hpp
/// Fundamental linear-algebra and complex types (ADR 03). All core computations use double.

#include <Eigen/Core>
#include <complex>

#include "rtt/math/real.hpp"

namespace rtt::math {

using Vec3 = Eigen::Vector3d;          ///< 3-vector, global or local coordinates in mm.
using Mat3 = Eigen::Matrix3d;          ///< 3x3 real matrix.
using Complex = std::complex<double>;  ///< Complex amplitude or complex refractive index.
using CMat3 = Eigen::Matrix3cd;        ///< 3x3 complex matrix (polarization ray tracing).

template <Real T>
using Vec3T = Eigen::Matrix<T, 3, 1>;  ///< 3-vector over a generic scalar.

template <Real T>
using Mat3T = Eigen::Matrix<T, 3, 3>;  ///< 3x3 matrix over a generic scalar.

}  // namespace rtt::math
