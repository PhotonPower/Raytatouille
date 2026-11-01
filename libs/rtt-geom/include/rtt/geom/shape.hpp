#pragma once

/// @file shape.hpp
/// Requirements on a surface shape z = sag(x, y) in the local coordinate system of a surface.
///
/// Conventions (docs/architecture.md): lengths in mm, right-handed local coordinates with the
/// surface vertex at the origin and the optical axis along +z. The sag is measured along +z;
/// a positive curvature c = 1/R means that the centre of curvature lies on the +z side.
///
/// There is no common base class (ADR 0041, point 4): every shape is a final class, the tracer
/// visits the variant compile::CompiledShape statically, and generic code is constrained by the
/// concepts below. Shapes are immutable and therefore safe to share between threads.

#include <concepts>
#include <optional>
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

}  // namespace rtt::geom
