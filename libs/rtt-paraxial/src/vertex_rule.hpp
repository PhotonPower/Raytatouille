#pragma once

// Internal to rtt-paraxial: the rule for the rotationally symmetric paraxial model of one
// surface (ADR 0041, point 2), a pure function over the vertex expansion of its shape, so that
// the tests reach it with shapes that exist only there. Not installed, not part of the public API.

#include <optional>
#include <string>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/geom/shape.hpp"

namespace rtt::paraxial::detail {

/// The paraxial model of one surface in its local orientation (lengths in mm).
struct VertexRule {
  /// Why the surface has no first-order model; empty if it has one.
  std::optional<std::string> first_order_error;
  /// Why the surface has no fourth-order cap for the Seidel sums; empty if it has one. Always set
  /// when first_order_error is.
  std::optional<std::string> fourth_order_error;
  double value = 0.0;  ///< z(0, 0): shift of the vertex along the local z axis, mm
  double c = 0.0;      ///< vertex curvature d2z/dx2 = d2z/dy2, 1/mm
  double conic = 0.0;  ///< K of the cap a = (K c^3 + 8 A4) y^4 dn, dimensionless
  double a4 = 0.0;     ///< A4 of the cap, mm^-3
};

/// ADR 0041, point 2, checked in this order on the local expansion:
/// - order < 2: no paraxial model at the vertex;
/// - slope != 0: tilted at the vertex;
/// - d2z/dx2 != d2z/dy2 or d2z/dxdy != 0: not rotationally symmetric in second order;
/// - otherwise c = d2z/dx2 and the vertex shifts by `value`; the Seidel cap (K, A4) needs
///   order 4 and `rotational4` in addition.
/// The comparisons are exact: a rotationally symmetric shape assigns the same value to both
/// diagonal elements (ADR 0041, point 1), and +0 == -0.
[[nodiscard]] VertexRule vertex_rule(const geom::VertexExpansion<double>& e);

/// The rule for the shape of a compiled surface (std::visit over every alternative).
[[nodiscard]] VertexRule vertex_rule(const compile::CompiledSurface& s);

/// Global z of the vertex of a surface on the z axis (ADR 0041, point 2): origin_z +
/// sign(axis_z) value, with the same factor sign(axis_z) that turns c into global orientation.
/// For value == 0 it is origin_z bit for bit, also the sign of a zero.
/// @param origin_z global z of the local origin, mm
/// @param axis_z   global z component of the local z axis (+-1 on a rotationally symmetric path)
/// @param value    VertexRule::value, mm
[[nodiscard]] double vertex_z(double origin_z, double axis_z, double value);

}  // namespace rtt::paraxial::detail
