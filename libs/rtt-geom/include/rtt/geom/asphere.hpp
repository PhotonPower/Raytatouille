#pragma once

/// @file asphere.hpp
/// Even asphere: conic of revolution plus an even polynomial in r (lengths in mm).

#include <optional>
#include <utility>
#include <vector>

#include "rtt/geom/conic.hpp"
#include "rtt/geom/shape.hpp"

namespace rtt::geom {

/// Even asphere z(r) = c r^2 / (1 + sqrt(1 - (1 + k) c^2 r^2)) + A4 r^4 + A6 r^6 + ...,
/// r^2 = x^2 + y^2 (conic part: W. T. Welford, Aberrations of Optical Systems, Ch. 2; full
/// form: ISO 10110-12, equation of the rotationally symmetric aspheric surface).
/// Sign conventions of c and k as for Conic. The coefficients follow the model
/// (rtt/model/surface.hpp): coefficients[i] = A_(2 i + 4) in mm^(-(2 i + 3)).
template <rtt::math::Real T>
class EvenAsphere final : public Shape<T> {
 public:
  /// @param curvature      c = 1/R of the base conic in 1/mm
  /// @param conic_constant k of the base conic, dimensionless
  /// @param coefficients   A4, A6, A8, ... (coefficients[0] = A4), A_n in mm^(1 - n)
  EvenAsphere(T curvature, T conic_constant, std::vector<T> coefficients)
      : conic_(curvature, conic_constant), coefficients_(std::move(coefficients)) {}

  /// Polynomial coefficients A4, A6, ... in mm^(1 - n).
  [[nodiscard]] const std::vector<T>& coefficients() const noexcept { return coefficients_; }

  /// Sag z in mm at (x, y) in mm; NaN for r > max_radius().
  [[nodiscard]] T sag(T x, T y) const override;

  /// Slopes (dz/dx, dz/dy); NaN for r > max_radius(), +-infinity at r = max_radius().
  [[nodiscard]] std::pair<T, T> grad(T x, T y) const override;

  /// (c, k) of the base conic.
  [[nodiscard]] std::pair<T, T> base_conic() const override { return conic_.base_conic(); }

  /// Domain of the base conic; the polynomial is defined everywhere.
  [[nodiscard]] std::optional<T> max_radius() const override { return conic_.max_radius(); }

 private:
  Conic<T> conic_;
  std::vector<T> coefficients_;
};

template <rtt::math::Real T>
T EvenAsphere<T>::sag(T x, T y) const {
  // Polynomial part sum_i A_(2i+4) r^(2i+4) = u^2 sum_i A_(2i+4) u^i with u = r^2, evaluated
  // with Horner's scheme. Without coefficients it is exactly 0, so the sag equals the conic's.
  const T u = x * x + y * y;
  T acc = T(0);
  for (auto it = coefficients_.rbegin(); it != coefficients_.rend(); ++it) {
    acc = acc * u + *it;
  }
  return conic_.sag(x, y) + acc * u * u;
}

template <rtt::math::Real T>
std::pair<T, T> EvenAsphere<T>::grad(T x, T y) const {
  // d/dx of A_n r^n is n A_n r^(n-2) x, so the polynomial contributes
  // x u sum_i (2i + 4) A_(2i+4) u^i (and the same with y), u = r^2.
  const T u = x * x + y * y;
  T acc = T(0);
  int i = static_cast<int>(coefficients_.size());
  for (auto it = coefficients_.rbegin(); it != coefficients_.rend(); ++it) {
    --i;
    acc = acc * u + T(2 * i + 4) * *it;
  }
  const auto [cx, cy] = conic_.grad(x, y);
  return {cx + x * u * acc, cy + y * u * acc};
}

}  // namespace rtt::geom
