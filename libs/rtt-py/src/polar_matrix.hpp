#pragma once

/// @file polar_matrix.hpp
/// The rtt-polar templates prt_matrix() and geometric_transform() for double, instantiated in
/// their own translation unit (polar_matrix.cpp) because GCC -O2 reports -Wnull-dereference
/// inside Eigen for them (#62; the cause in rtt-polar is tracked in #35). Same parameters,
/// units and preconditions as rtt/polar/prt.hpp.

#include <complex>

#include "rtt/math/types.hpp"

namespace rtt::py::polar_matrix {

/// rtt::polar::prt_matrix<double>: P = a_s s s^T + a_p p_out p_in^T + k_out k_in^T (global
/// unit vectors; Lam, Eqs. (3.1)-(3.9)).
[[nodiscard]] math::CMat3 prt_matrix(const math::Vec3& k_in,
                                     const math::Vec3& k_out,
                                     const math::Vec3& normal,
                                     std::complex<double> a_s,
                                     std::complex<double> a_p);

/// rtt::polar::geometric_transform<double>: Q = s s^T +/- p_out p_in^T + k_out k_in^T.
[[nodiscard]] math::Mat3 geometric_transform(const math::Vec3& k_in,
                                             const math::Vec3& k_out,
                                             const math::Vec3& normal,
                                             bool reflection);

}  // namespace rtt::py::polar_matrix
