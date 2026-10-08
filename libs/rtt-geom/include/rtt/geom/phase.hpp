#pragma once

/// @file phase.hpp
/// Phase functions of the phase layers of a surface and their gradient (ADR 0025, #126).
///
/// Conventions (ADR 0025, point 1; docs/architecture.md):
/// - The phase phi(x, y) is in rad, a function of the local surface coordinates (x, y) in mm,
///   i.e. of the lateral position, not of the arc length on a curved surface. It describes a
///   wavelength-independent groove structure: phi / (2 pi) is the characteristic function F of
///   M. Mansuripur, Proc. SPIE 6620, 66200N (2007), Sec. 2 (F counts periods).
/// - The phase of a surface is the sum of its layers (phase(), phase_grad()).
/// - Sign: phi grows along the gradient; order m adds m phi lambda_0 / (2 pi) to the optical
///   path (a larger phase is a delay, ADR 0025, point 3) and m lambda_0 / (2 pi) times the
///   tangential gradient to n t (ADR 0025, point 2, Mansuripur Eq. (7b)).
/// - The model (rtt::model::LinearGrating, RadialPhase) and the file format store the same
///   numbers: lines per mm, coefficients in rad (no factor 2 pi, not in waves), coefficient k
///   for rho^(2k) starting at k = 1. A conversion only turns parameters into values and
///   orientation_deg into rad, without a factor 2 pi and without a change of sign.
/// - Gradients are in rad/mm. The functions do not throw; input errors are rejected by the
///   constructors (std::invalid_argument).

#include <cmath>
#include <cstddef>
#include <numbers>
#include <span>
#include <stdexcept>
#include <utility>
#include <variant>
#include <vector>

#include "rtt/math/real.hpp"
#include "rtt/math/types.hpp"

namespace rtt::geom {

/// Straight-line grating: phi = 2 pi G (x cos psi + y sin psi) (ADR 0025, point 1). The grating
/// vector (cos psi, sin psi) is normal to the grooves; psi = 0 puts the grooves parallel to the
/// local y axis. phi(0, 0) = 0.
template <rtt::math::Real T>
class LinearGratingPhase {
 public:
  /// @param lines_per_mm G in 1/mm, finite (a negative G reverses the grating vector)
  /// @param orientation  psi in rad, from the local x axis towards the local y axis (right-handed
  ///                     about the local +z axis), finite
  /// @throws std::invalid_argument if a value is not finite
  LinearGratingPhase(T lines_per_mm, T orientation);

  /// G in 1/mm.
  [[nodiscard]] T lines_per_mm() const noexcept { return lines_per_mm_; }
  /// psi in rad.
  [[nodiscard]] T orientation() const noexcept { return orientation_; }

  /// phi(x, y) in rad at the local point (x, y) in mm.
  [[nodiscard]] T phase(T x, T y) const noexcept { return T(0) * x * y; }  // STUB

  /// (dphi/dx, dphi/dy) = 2 pi G (cos psi, sin psi) in rad/mm, independent of (x, y).
  [[nodiscard]] std::pair<T, T> grad(T /*x*/, T /*y*/) const noexcept {
    return {T(0) * gx_, T(0) * gy_};
  }  // STUB

 private:
  T lines_per_mm_;
  T orientation_;
  T gx_;  ///< 2 pi G cos psi, rad/mm
  T gy_;  ///< 2 pi G sin psi, rad/mm
};

/// Rotationally symmetric phase phi = sum_k c_k rho^(2k), k = 1..N, rho = r / R with
/// r^2 = x^2 + y^2 and R the normalization radius (ADR 0025, point 1). With u = rho^2 the
/// gradient is (2 x / R^2, 2 y / R^2) dphi/du, dphi/du = sum_k k c_k u^(k-1); it is regular at
/// r = 0.
///
/// Relation to Mansuripur, Proc. SPIE 6620, 66200N (2007), Eq. (9), F(r) = sum_n a_n r^n:
/// phi = 2 pi F gives a_n = c_(n/2) / (2 pi R^n) for even n and a_n = 0 for odd n.
template <rtt::math::Real T>
class RadialPhasePolynomial {
 public:
  /// @param normalization_radius R in mm, finite and > 0
  /// @param coefficients         c_1, c_2, ... in rad for rho^2, rho^4, ...; finite; may be
  ///                             empty (phi = 0)
  /// @throws std::invalid_argument for R <= 0 or a value that is not finite
  RadialPhasePolynomial(T normalization_radius, std::vector<T> coefficients);

  /// R in mm.
  [[nodiscard]] T normalization_radius() const noexcept { return radius_; }
  /// c_1, c_2, ... in rad.
  [[nodiscard]] const std::vector<T>& coefficients() const noexcept { return coefficients_; }

  /// phi(x, y) in rad at the local point (x, y) in mm.
  [[nodiscard]] T phase(T x, T y) const noexcept;

  /// (dphi/dx, dphi/dy) in rad/mm at (x, y) in mm.
  [[nodiscard]] std::pair<T, T> grad(T x, T y) const noexcept;

 private:
  T radius_;
  std::vector<T> coefficients_;
};

/// One phase layer of a surface.
template <rtt::math::Real T>
using PhaseFunction = std::variant<LinearGratingPhase<T>, RadialPhasePolynomial<T>>;

/// Phase of a surface: the sum of its layers at (x, y) in mm, in rad; 0 without layers.
template <rtt::math::Real T>
[[nodiscard]] T phase(std::span<const PhaseFunction<T>> layers, T x, T y) noexcept;

/// Gradient of the phase of a surface: the sum of the layer gradients at (x, y) in mm, in
/// rad/mm; (0, 0) without layers.
template <rtt::math::Real T>
[[nodiscard]] std::pair<T, T> phase_grad(std::span<const PhaseFunction<T>> layers,
                                         T x,
                                         T y) noexcept;

/// Tangential part g_par = (I - N N^T) g of the lateral phase gradient g = (gx, gy, 0) at a
/// point with unit normal N, all in local surface coordinates (ADR 0025, point 2). This is the
/// surface gradient of phi: g_par . dr = dphi for every tangential displacement dr. On a
/// surface of revolution with rotationally symmetric phi it equals Mansuripur, Proc. SPIE 6620,
/// 66200N (2007), Eq. (11), dF/ds = (dF/dr) / sqrt(1 + (dh/dr)^2), along the meridional tangent.
/// @param g           (dphi/dx, dphi/dy) in rad/mm
/// @param unit_normal N, |N| = 1 (either orientation; not checked)
/// @return g_par in rad/mm
template <rtt::math::Real T>
[[nodiscard]] math::Vec3T<T> tangential_gradient(std::pair<T, T> g,
                                                 const math::Vec3T<T>& unit_normal) noexcept;

// ------------------------------------------------------------------- implementation -----

template <rtt::math::Real T>
LinearGratingPhase<T>::LinearGratingPhase(T lines_per_mm, T orientation)
    : lines_per_mm_(lines_per_mm), orientation_(orientation), gx_(T(0)), gy_(T(0)) {
  using std::cos;
  using std::isfinite;
  using std::sin;
  if (!isfinite(lines_per_mm) || !isfinite(orientation)) {
    throw std::invalid_argument("geom: grating lines_per_mm and orientation must be finite");
  }
  const T k = T(2.0 * std::numbers::pi) * lines_per_mm;
  gx_ = k * cos(orientation);
  gy_ = k * sin(orientation);
}

template <rtt::math::Real T>
RadialPhasePolynomial<T>::RadialPhasePolynomial(T normalization_radius, std::vector<T> coefficients)
    : radius_(normalization_radius), coefficients_(std::move(coefficients)) {
  using std::isfinite;
  if (!isfinite(radius_) || !(radius_ > T(0))) {
    throw std::invalid_argument("geom: radial phase normalization_radius must be finite and > 0");
  }
  for (const T& c : coefficients_) {
    if (!isfinite(c)) {
      throw std::invalid_argument("geom: radial phase coefficients must be finite");
    }
  }
}

// STUB for the red run of the #126 tests: zero phase, zero gradient.
template <rtt::math::Real T>
T RadialPhasePolynomial<T>::phase(T x, T y) const noexcept {
  return T(0) * x * y;
}

template <rtt::math::Real T>
std::pair<T, T> RadialPhasePolynomial<T>::grad(T x, T y) const noexcept {
  return {T(0) * x, T(0) * y};
}

template <rtt::math::Real T>
T phase(std::span<const PhaseFunction<T>> layers, T x, T y) noexcept {
  return T(0) * x * y * static_cast<double>(layers.size());
}

template <rtt::math::Real T>
std::pair<T, T> phase_grad(std::span<const PhaseFunction<T>> layers, T x, T y) noexcept {
  const T zero = T(0) * static_cast<double>(layers.size());
  return {zero * x, zero * y};
}

template <rtt::math::Real T>
math::Vec3T<T> tangential_gradient(std::pair<T, T> g, const math::Vec3T<T>& unit_normal) noexcept {
  return math::Vec3T<T>(T(0) * g.first, T(0) * g.second, T(0) * unit_normal.x());
}

}  // namespace rtt::geom
