// GCC -O2 reports -Wnull-dereference inside Eigen's lazy products when rtt::polar::prt_matrix()
// and geometric_transform() (templates in rtt/polar/prt.hpp) are instantiated (CI of #62, as in
// #57 for prt_analysis.cpp and test_prt.cpp). The warning is attributed to the Eigen headers, so
// the suppression has to cover the includes; it is limited to this small translation unit and to
// GCC. The cause in rtt-polar is tracked in #35; remove this pragma when it is fixed there.
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wnull-dereference"
#endif

#include "polar_matrix.hpp"

#include "rtt/polar/prt.hpp"

namespace rtt::py::polar_matrix {

math::CMat3 prt_matrix(const math::Vec3& k_in,
                       const math::Vec3& k_out,
                       const math::Vec3& normal,
                       std::complex<double> a_s,
                       std::complex<double> a_p) {
  return rtt::polar::prt_matrix(k_in, k_out, normal, a_s, a_p);
}

math::Mat3 geometric_transform(const math::Vec3& k_in,
                               const math::Vec3& k_out,
                               const math::Vec3& normal,
                               bool reflection) {
  return rtt::polar::geometric_transform(k_in, k_out, normal, reflection);
}

}  // namespace rtt::py::polar_matrix

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif
