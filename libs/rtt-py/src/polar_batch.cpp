#include "polar_batch.hpp"

#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>

#include "rtt/polar/ideal.hpp"
#include "rtt/polar/prt_analysis.hpp"

namespace rtt::py::polar_batch {
namespace {

using math::CMat3;
using math::Complex;
using math::Vec3;
using trace::RayBatch;

Vec3 direction(const RayBatch& rays, std::size_t i) {
  return {rays.dir_x()[i], rays.dir_y()[i], rays.dir_z()[i]};
}

// Products and dot products are written as explicit loops over the three components: GCC -O2
// reports -Wnull-dereference inside Eigen's lazy product evaluators for them (CI of #62).

/// Sum_i a_i b_i of a real and a complex vector (no conjugation).
Complex dot(const Vec3& a, const CVec3& b) {
  Complex sum(0.0);
  for (int i = 0; i < 3; ++i) sum += a[i] * b[i];
  return sum;
}

/// Sum_i a_i b_i of two real vectors.
double dot(const Vec3& a, const Vec3& b) {
  double sum = 0.0;
  for (int i = 0; i < 3; ++i) sum += a[i] * b[i];
  return sum;
}

/// P E.
CVec3 times(const CMat3& p, const CVec3& e) {
  CVec3 out;
  for (int r = 0; r < 3; ++r) {
    Complex sum(0.0);
    for (int c = 0; c < 3; ++c) sum += p(r, c) * e[c];
    out[r] = sum;
  }
  return out;
}

/// k0 = Re(P^T k) (see initial_directions()).
Vec3 initial_direction(const CMat3& p, const Vec3& k) {
  Vec3 k0;
  for (int c = 0; c < 3; ++c) {
    double sum = 0.0;
    for (int r = 0; r < 3; ++r) sum += p(r, c).real() * k[r];
    k0[c] = sum;
  }
  return k0;
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
  if (states.empty() || (states.size() != 1 && states.size() != k0.size())) {
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
    if (!(std::abs(dot(k0[i], e)) <= kStateTolerance)) {
      throw std::invalid_argument(where + " is not transverse to the initial direction of ray " +
                                  std::to_string(i) +
                                  " (use transverse_polarization to project it)");
    }
  }
}

/// ||P_T||_F^2 with P_T = (I - k k^T) P (ADR 0021; as in rtt-trace apply_event.cpp): every
/// column of P minus its component along k.
double transverse_norm2(const CMat3& p, const Vec3& k) {
  double sum = 0.0;
  for (int c = 0; c < 3; ++c) {
    const CVec3 column(p(0, c), p(1, c), p(2, c));
    const Complex along = dot(k, column);
    for (int r = 0; r < 3; ++r) sum += std::norm(column[r] - k[r] * along);
  }
  return sum;
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
    // Normalised first: |k0| = 1 only up to rounding, and for a projection just above the
    // threshold the residual k0 . E would otherwise be (|k0|^2 - 1) |k0 . e| / |projection|.
    const Vec3 k = k0[i].normalized();
    const Complex along = dot(k, e);
    CVec3 projected;
    for (int c = 0; c < 3; ++c) projected[c] = e[c] - k[c] * along;
    if (!(projected.norm() >= kStateTolerance * e.norm()) || e.norm() == 0.0) {
      throw std::invalid_argument("polarization is parallel to the initial direction of ray " +
                                  std::to_string(i));
    }
    states[i] = projected / projected.norm();
  }
  return states;
}

std::vector<double> transmission(const RayBatch& rays) {
  const auto weight = rays.weight();
  return {weight.begin(), weight.end()};
}

std::vector<double> transmission(const RayBatch& rays, std::span<const CVec3> states) {
  const auto weight = rays.weight();
  const std::vector<Vec3> k0 = initial_directions(rays);
  check_states(states, k0);
  std::vector<double> power(rays.size());
  for (std::size_t i = 0; i < rays.size(); ++i) {
    const CMat3 p = rays.prt_matrix(i);
    const double half_norm = 0.5 * transverse_norm2(p, direction(rays, i));
    power[i] =
        half_norm > 0.0 ? weight[i] * times(p, state(states, i)).squaredNorm() / half_norm : 0.0;
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
    if (dot(k, initial_direction(p, k)) >= 1.0 - kSameDirection) {
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
    // The precondition of rtt::polar::transverse_axis(), with the same expression.
    const Vec3 projected = axis - dot(axis, k) * k;
    if (!(projected.norm() >= rtt::polar::kAxisAlongK * axis.norm()) || axis.norm() == 0.0) {
      throw std::invalid_argument("axis is parallel to the direction of ray " + std::to_string(i));
    }
    const CVec3 e = times(rays.prt_matrix(i), state(states, i));
    const rtt::polar::Stokes<double> si = rtt::polar::stokes<double>(e, axis, k);
    s[i] = {si.s0, si.s1, si.s2, si.s3};
  }
  return s;
}

}  // namespace rtt::py::polar_batch
