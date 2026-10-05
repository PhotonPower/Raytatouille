#include "polar_batch.hpp"

#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>

#include "rtt/polar/ideal.hpp"
#include "rtt/polar/prt_analysis.hpp"

namespace rtt::py::polar {
namespace {

using math::CMat3;
using math::Complex;
using math::Vec3;
using trace::RayBatch;

Vec3 direction(const RayBatch& rays, std::size_t i) {
  return {rays.dir_x()[i], rays.dir_y()[i], rays.dir_z()[i]};
}

/// k0 = Re(P^T k) (see initial_directions()).
Vec3 initial_direction(const CMat3& p, const Vec3& k) {
  return (p.transpose() * k.cast<Complex>()).real();
}

bool finite(const CVec3& v) {
  return v.allFinite();
}

/// The state for ray i: the only one, or the i-th.
const CVec3& state(std::span<const CVec3> states, std::size_t i) {
  return states.size() == 1 ? states[0] : states[i];
}

/// Checks states against the initial directions (API boundary, ADR 0009).
void check_states(std::span<const CVec3> states, std::span<const Vec3> k0) {
  if (states.size() != 1 && states.size() != k0.size()) {
    throw std::invalid_argument("polarization must be one state (3,) or one per ray (N, 3) for " +
                                std::to_string(k0.size()) + " rays, got " +
                                std::to_string(states.size()) + " states");
  }
  for (std::size_t i = 0; i < k0.size(); ++i) {
    const CVec3& e = state(states, i);
    const std::string where =
        states.size() == 1 ? "polarization" : "polarization of ray " + std::to_string(i);
    if (!finite(e)) throw std::invalid_argument(where + " is not finite");
    if (!(std::abs(e.norm() - 1.0) <= kStateTolerance)) {
      throw std::invalid_argument(where + " is not a unit vector (|E| = 1 within 1e-12)");
    }
    // Complex E . k0 without conjugation: the component of E along the real k0.
    if (!(std::abs(k0[i].cast<Complex>().dot(e)) <= kStateTolerance)) {
      throw std::invalid_argument(where + " is not transverse to the initial direction of ray " +
                                  std::to_string(i) +
                                  " (use transverse_polarization to project it)");
    }
  }
}

/// ||P_T||_F^2 with P_T = (I - k k^T) P (ADR 0021; as in rtt-trace apply_event.cpp).
double transverse_norm2(const CMat3& p, const Vec3& k) {
  const math::Mat3 projector = math::Mat3::Identity() - k * k.transpose();
  return (projector.cast<Complex>() * p).squaredNorm();
}

}  // namespace

std::vector<Vec3> initial_directions(const RayBatch& rays) {
  std::vector<Vec3> k0(rays.size());
  for (std::size_t i = 0; i < rays.size(); ++i) {
    k0[i] = initial_direction(rays.prt_matrix(i), direction(rays, i));
  }
  return k0;
}

std::vector<CVec3> transverse_polarization(const RayBatch& rays, const CVec3& e) {
  if (!finite(e)) throw std::invalid_argument("polarization is not finite");
  const std::vector<Vec3> k0 = initial_directions(rays);
  std::vector<CVec3> states(rays.size());
  for (std::size_t i = 0; i < rays.size(); ++i) {
    const CVec3 k = k0[i].cast<Complex>();
    const CVec3 projected = e - k * k.dot(e);
    if (!(projected.norm() >= kStateTolerance * e.norm()) || e.norm() == 0.0) {
      throw std::invalid_argument("polarization is parallel to the initial direction of ray " +
                                  std::to_string(i));
    }
    states[i] = projected / projected.norm();
  }
  return states;
}

std::vector<double> transmission(const RayBatch& rays, std::span<const CVec3> states) {
  const auto weight = rays.weight();
  if (states.empty()) return {weight.begin(), weight.end()};
  const std::vector<Vec3> k0 = initial_directions(rays);
  check_states(states, k0);
  std::vector<double> power(rays.size());
  for (std::size_t i = 0; i < rays.size(); ++i) {
    const CMat3 p = rays.prt_matrix(i);
    const double half_norm = 0.5 * transverse_norm2(p, direction(rays, i));
    power[i] = half_norm > 0.0 ? weight[i] * (p * state(states, i)).squaredNorm() / half_norm : 0.0;
  }
  return power;
}

Diattenuations diattenuation(const RayBatch& rays) {
  Diattenuations d;
  d.value.resize(rays.size());
  d.maximum.resize(rays.size());
  d.minimum.resize(rays.size());
  d.axis.resize(rays.size());
  for (std::size_t i = 0; i < rays.size(); ++i) {
    const CMat3 p = rays.prt_matrix(i);
    const Vec3 k = direction(rays, i);
    const rtt::polar::Diattenuation di = rtt::polar::diattenuation(p, initial_direction(p, k), k);
    d.value[i] = di.value;
    d.maximum[i] = di.maximum;
    d.minimum[i] = di.minimum;
    d.axis[i] = di.axis;
  }
  return d;
}

Retardances retardance(const RayBatch& rays) {
  constexpr double nan = std::numeric_limits<double>::quiet_NaN();
  Retardances r;
  r.value.resize(rays.size());
  r.fast_axis.resize(rays.size());
  for (std::size_t i = 0; i < rays.size(); ++i) {
    const CMat3 p = rays.prt_matrix(i);
    const Vec3 k = direction(rays, i);
    if (k.dot(initial_direction(p, k)) >= 1.0 - kSameDirection) {
      const rtt::polar::Retardance ri = rtt::polar::retardance(p, k);
      r.value[i] = ri.value;
      r.fast_axis[i] = ri.fast_axis;
    } else {
      r.value[i] = nan;
      r.fast_axis[i] = CVec3::Constant(Complex(nan, nan));
    }
  }
  return r;
}

std::vector<std::array<double, 4>> stokes(const RayBatch& rays,
                                          std::span<const CVec3> states,
                                          const Vec3& axis) {
  if (!axis.allFinite()) throw std::invalid_argument("axis is not finite");
  const std::vector<Vec3> k0 = initial_directions(rays);
  check_states(states, k0);
  std::vector<std::array<double, 4>> s(rays.size());
  for (std::size_t i = 0; i < rays.size(); ++i) {
    const Vec3 k = direction(rays, i);
    if (!(axis.cross(k).norm() >= rtt::polar::kAxisAlongK * axis.norm()) || axis.norm() == 0.0) {
      throw std::invalid_argument("axis is parallel to the direction of ray " + std::to_string(i));
    }
    const CVec3 e = rays.prt_matrix(i) * state(states, i);
    const rtt::polar::Stokes<double> si = rtt::polar::stokes<double>(e, axis, k);
    s[i] = {si.s0, si.s1, si.s2, si.s3};
  }
  return s;
}

}  // namespace rtt::py::polar
