#pragma once

/// @file shape.hpp
/// Abstract surface shape z = sag(x, y) in the local coordinate system of a surface.
///
/// Conventions (docs/architecture.md): lengths in mm, right-handed local coordinates with the
/// surface vertex at the origin and the optical axis along +z. The sag is measured along +z;
/// a positive curvature c = 1/R means that the centre of curvature lies on the +z side.

#include <optional>
#include <utility>

#include "rtt/math/real.hpp"

namespace rtt::geom {

/// Rotationally symmetric or freeform surface shape given as a height function over (x, y).
/// Implementations are immutable and therefore safe to share between threads.
template <rtt::math::Real T>
class Shape {
 public:
  Shape() = default;
  Shape(const Shape&) = default;
  Shape(Shape&&) noexcept = default;
  Shape& operator=(const Shape&) = default;
  Shape& operator=(Shape&&) noexcept = default;
  virtual ~Shape() = default;

  /// Surface height z at the local point (x, y), all in mm. Outside the domain of the shape
  /// (see max_radius()) the result is NaN; the tracer must check the domain or use a status.
  [[nodiscard]] virtual T sag(T x, T y) const = 0;

  /// Surface slopes (dz/dx, dz/dy) at (x, y) in mm, dimensionless. NaN outside the domain.
  [[nodiscard]] virtual std::pair<T, T> grad(T x, T y) const = 0;

  /// Best-fit base conic (curvature c in 1/mm, conic constant k) used as starting value for
  /// the intersection.
  [[nodiscard]] virtual std::pair<T, T> base_conic() const = 0;

  /// Largest radial distance r = sqrt(x^2 + y^2) in mm on which the shape is defined, or
  /// std::nullopt if the shape is defined for all r.
  [[nodiscard]] virtual std::optional<T> max_radius() const = 0;
};

}  // namespace rtt::geom
