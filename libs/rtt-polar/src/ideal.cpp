#include "rtt/polar/ideal.hpp"

#include <cmath>
#include <stdexcept>

namespace rtt::polar {

namespace {

void check_direction(const math::Vec3& k) {
  if (!k.allFinite()) throw std::invalid_argument("ideal element: k is not finite");
  if (std::abs(k.norm() - 1.0) > 1e-9) {
    throw std::invalid_argument("ideal element: k must be a unit vector");
  }
}

void check_axis(const math::Vec3& axis, const math::Vec3& k) {
  if (!axis.allFinite()) throw std::invalid_argument("ideal element: axis is not finite");
  const math::Vec3 transverse = axis - axis.dot(k) * k;
  if (!(transverse.norm() >= kAxisAlongK * axis.norm()) || axis.norm() == 0.0) {
    throw std::invalid_argument("ideal element: axis is parallel to the ray");
  }
}

}  // namespace

math::CMat3 ideal_polarizer_checked(const math::Vec3& axis,
                                    double extinction_ratio,
                                    const math::Vec3& k) {
  check_direction(k);
  check_axis(axis, k);
  if (!std::isfinite(extinction_ratio) || extinction_ratio < 0.0 || extinction_ratio > 1.0) {
    throw std::invalid_argument("ideal polarizer: extinction ratio must lie in [0, 1]");
  }
  return ideal_polarizer(axis, extinction_ratio, k);
}

math::CMat3 linear_retarder_checked(const math::Vec3& fast_axis,
                                    double retardance_waves,
                                    const math::Vec3& k) {
  check_direction(k);
  check_axis(fast_axis, k);
  if (!std::isfinite(retardance_waves)) {
    throw std::invalid_argument("linear retarder: retardance is not finite");
  }
  return linear_retarder(fast_axis, retardance_waves, k);
}

math::CMat3 circular_retarder_checked(double retardance_waves, const math::Vec3& k) {
  check_direction(k);
  if (!std::isfinite(retardance_waves)) {
    throw std::invalid_argument("circular retarder: retardance is not finite");
  }
  return circular_retarder(retardance_waves, k);
}

}  // namespace rtt::polar
