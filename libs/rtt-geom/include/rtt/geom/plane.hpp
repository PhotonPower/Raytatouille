#pragma once

/// @file plane.hpp
/// Flat surface z = 0 in local coordinates (lengths in mm).

#include <optional>
#include <utility>

#include "rtt/geom/shape.hpp"

namespace rtt::geom {

/// Plane through the local origin, perpendicular to the local z axis. Flat surfaces are always
/// a Plane, never a conic with radius 0 or infinity (docs/architecture.md, Konventionen).
template <rtt::math::Real T>
class Plane final : public Shape<T> {
 public:
  /// Height is zero everywhere, mm.
  [[nodiscard]] T sag(T /*x*/, T /*y*/) const override { return T(0); }

  /// Slopes are zero everywhere.
  [[nodiscard]] std::pair<T, T> grad(T /*x*/, T /*y*/) const override { return {T(0), T(0)}; }

  /// Curvature 0 and conic constant 0.
  [[nodiscard]] std::pair<T, T> base_conic() const override { return {T(0), T(0)}; }

  /// Unbounded.
  [[nodiscard]] std::optional<T> max_radius() const override { return std::nullopt; }
};

}  // namespace rtt::geom
