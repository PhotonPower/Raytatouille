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
class Plane final {
 public:
  /// Height is zero everywhere, mm.
  [[nodiscard]] T sag(T /*x*/, T /*y*/) const noexcept { return T(0); }

  /// Slopes are zero everywhere.
  [[nodiscard]] std::pair<T, T> grad(T /*x*/, T /*y*/) const noexcept { return {T(0), T(0)}; }

  /// Curvature 0 and conic constant 0.
  [[nodiscard]] std::pair<T, T> base_conic() const { return {T(0), T(0)}; }

  /// Unbounded.
  [[nodiscard]] std::optional<T> max_radius() const { return std::nullopt; }

  /// Flat at the origin: value 0, slope 0, Hessian 0, order 4, rotational, cap K = 0, A4 = 0
  /// (ADR 0041, point 2).
  [[nodiscard]] VertexExpansion<T> vertex_expansion() const noexcept {
    VertexExpansion<T> e;
    e.order = 4;
    e.rotational4 = true;
    return e;
  }
};

}  // namespace rtt::geom
