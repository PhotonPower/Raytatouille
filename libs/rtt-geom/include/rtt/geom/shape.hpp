#pragma once

/// @file shape.hpp
/// Requirements on a surface shape z = sag(x, y) in the local coordinate system of a surface.
///
/// Conventions (docs/architecture.md): lengths in mm, right-handed local coordinates with the
/// optical axis along +z. The sag is measured along +z; a positive curvature c = 1/R means that
/// the centre of curvature lies on the +z side. The vertex of a surface is the local point
/// (0, 0, sag(0, 0)); it lies at the origin unless the shape has a constant part (`value` of its
/// VertexExpansion, ADR 0041, point 1).
///
/// There is no common base class (ADR 0041, point 4): every shape is a final class, the tracer
/// visits the variant compile::CompiledShape statically, and generic code is constrained by the
/// concepts below. Shapes are immutable and therefore safe to share between threads.

#include <algorithm>
#include <array>
#include <concepts>
#include <cstddef>
#include <optional>
#include <span>
#include <utility>

#include "rtt/math/real.hpp"

namespace rtt::geom {

/// A shape the Newton intersection (intersect_newton() in intersect.hpp) can work on:
/// - `sag(x, y)`: surface height z in mm at the local point (x, y) in mm. Outside the domain of
///   the shape (see `max_radius()`) the result is NaN; the tracer checks the domain or uses a
///   status. Must not throw.
/// - `grad(x, y)`: slopes (dz/dx, dz/dy), dimensionless. NaN outside the domain; a component may
///   be +-infinity where the surface is vertical (rim r = max_radius()). Must not throw.
/// - `base_conic()`: best-fit base conic (curvature c in 1/mm, conic constant k), the starting
///   value of the intersection; not the paraxial curvature (ADR 0041, point 1). Must not throw
///   (called from the noexcept intersection).
/// - `max_radius()`: largest radial distance r = sqrt(x^2 + y^2) in mm on which the shape is
///   defined, or std::nullopt if it is defined for all r. Must not throw (as base_conic()).
template <class S, class T>
concept IntersectableShape = rtt::math::Real<T> && requires(const S& shape, T x, T y) {
  { shape.sag(x, y) } noexcept -> std::same_as<T>;
  { shape.grad(x, y) } noexcept -> std::same_as<std::pair<T, T>>;
  { shape.base_conic() } -> std::same_as<std::pair<T, T>>;
  { shape.max_radius() } -> std::same_as<std::optional<T>>;
};

/// Taylor expansion of the sag z(x, y) at the vertex (x, y) = (0, 0) in the local coordinates of
/// the shape (ADR 0041, point 1). The paraxial rule of rtt-paraxial reads it (ADR 0041, point 2);
/// it decides from it, not from base_conic(), the curvature, the vertex position and the Seidel
/// cap of a surface.
template <rtt::math::Real T>
struct VertexExpansion {
  /// z(0, 0), the constant part, mm.
  T value = T(0);
  /// (dz/dx, dz/dy) at the vertex, dimensionless; without meaning for order 0.
  std::pair<T, T> slope{T(0), T(0)};
  /// (d2z/dx2, d2z/dxdy, d2z/dy2) at the vertex, 1/mm; without meaning for order < 2. A
  /// rotationally symmetric shape forms the common curvature once and assigns it to both
  /// diagonal elements, so that they compare equal bit for bit.
  std::array<T, 3> hessian{T(0), T(0), T(0)};
  /// Regularity: the highest order up to which the expansion exists, at most 4 (0 for a cone at
  /// the vertex, 2 with an |r|^3 term, 4 for a smooth shape).
  int order = 0;
  /// The expansion is rotationally symmetric up to fourth order, decided from the input
  /// parameters of the shape, not from computed coefficients.
  bool rotational4 = false;
  /// Fourth-order part as the cap of the Seidel sums, a = (K c^3 + 8 A4) y^4 dn (Sasian,
  /// OPTI 518, L14, pp. 16 and 17; docs/quellen.md): the conic constant K (dimensionless) and
  /// A4 (mm^-3), with c the vertex curvature from `hessian`. Only valid for order 4 and
  /// `rotational4`.
  T conic = T(0);
  T a4 = T(0);
};

/// A shape the paraxial rule can work on: `vertex_expansion()` gives its VertexExpansion<T>.
template <class S, class T>
concept ParaxialShape = rtt::math::Real<T> && requires(const S& shape) {
  { shape.vertex_expansion() } -> std::same_as<VertexExpansion<T>>;
};

/// Expansion of a stack of a base shape and additive terms (ShapeStack, ADR 0041, point 1):
/// `value`, `slope`, `hessian` and the r^4 coefficients b4 add, `order` is the minimum and
/// `rotational4` the logical and. The cap of the sum is K = 0 and A4 = b4 - c^3 / 8 with the sum
/// c of the curvatures and the sum b4; without terms it is the expansion of the base unchanged,
/// (K, A4) included.
///
/// b4 of each part is (1 + K) c^3 / 8 + A4 with its c = d2z/dx2: with eps = (1 + K) c^2 r^2 and
/// sqrt(1 - eps) = 1 - eps/2 - eps^2/8 + ... the r^4 coefficient of a conic with A4 is
/// (1 + K) c^3 / 8 + A4, so K c^3 + 8 A4 = 8 b4 - c^3 for every K (ADR 0041, point 1;
/// docs/quellen.md, row "Asphärenbeitrag"; equivalent to arXiv:2202.11378, Eq. (8)). Only valid
/// where the parts have order 4 and are rotational.
/// @param base  expansion of the base shape
/// @param terms expansions of the terms, in their order in the stack
template <rtt::math::Real T>
[[nodiscard]] VertexExpansion<T> stack_expansion(const VertexExpansion<T>& base,
                                                 std::span<const VertexExpansion<T>> terms) {
  if (terms.empty()) return base;
  const auto b4 = [](const VertexExpansion<T>& e) {
    const T c = e.hessian[0];
    return (T(1) + e.conic) * c * c * c / T(8) + e.a4;
  };
  VertexExpansion<T> sum = base;
  T b4_sum = b4(base);
  for (const VertexExpansion<T>& term : terms) {
    sum.value += term.value;
    sum.slope.first += term.slope.first;
    sum.slope.second += term.slope.second;
    for (std::size_t i = 0; i < sum.hessian.size(); ++i) sum.hessian[i] += term.hessian[i];
    sum.order = std::min(sum.order, term.order);
    sum.rotational4 = sum.rotational4 && term.rotational4;
    b4_sum += b4(term);
  }
  const T c = sum.hessian[0];
  sum.conic = T(0);
  sum.a4 = b4_sum - c * c * c / T(8);
  return sum;
}

}  // namespace rtt::geom
