#pragma once

// Internal to rtt-paraxial: the rule for the rotationally symmetric paraxial model of one
// surface (ADR 0041, point 2) and its extension per section for anamorphic surfaces (ADR 0042,
// point 2), pure functions over the vertex expansion of a shape, so that the tests reach them
// with shapes that exist only there. Not installed, not part of the public API.

#include <optional>
#include <string>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/geom/shape.hpp"
#include "rtt/math/types.hpp"
#include "rtt/paraxial/paraxial.hpp"

namespace rtt::paraxial::detail {

/// Largest off-axis vertex position (mm) and largest transverse component of a surface's local
/// z axis that still count as rotationally symmetric (decided for #7); also the relative limit
/// for principal axes aligned with x and y (ADR 0042, point 2).
inline constexpr double kSymmetryTolerance = 1e-12;

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

/// The vertex expansion of the shape of a compiled surface (std::visit over every alternative).
[[nodiscard]] geom::VertexExpansion<double> vertex_expansion(const compile::CompiledSurface& s);

/// The rule for the shape of a compiled surface (std::visit over every alternative).
[[nodiscard]] VertexRule vertex_rule(const compile::CompiledSurface& s);

/// True if the expansion has order >= 2 and d2z/dx2 != d2z/dy2 or d2z/dxdy != 0 (exact
/// comparisons, as vertex_rule()): the shape is anamorphic in second order (ADR 0042, point 2,
/// "Erkennung"). For order < 2 the Hessian has no meaning and the shape counts as not
/// anamorphic; vertex_rule() and section_rule() reject it alike.
[[nodiscard]] bool is_anamorphic(const geom::VertexExpansion<double>& e);

/// Orientation of a surface on the z axis (ADR 0042, point 2): p_g = o + M p_l with
/// M = diag(Q, sigma); Q is the transverse 2x2 block of the rotation (row-major, q_ij = M_ij for
/// i, j in {x, y}) and sigma = +1 if the local z axis points along global +z, else -1.
struct SectionFrame {
  double q00 = 1.0;
  double q01 = 0.0;
  double q10 = 0.0;
  double q11 = 1.0;
  double sigma = 1.0;
};

/// The frame of `rotation` (the rotation of CompiledSurface::to_global); sigma from the sign of
/// its element (2, 2), as first_order() takes the sign of the local z axis.
[[nodiscard]] SectionFrame section_frame(const math::Mat3& rotation);

/// The paraxial model of one surface in one section, in global orientation (lengths in mm).
struct SectionRule {
  /// Why the surface has no model in this section (ends with the ADR); empty if it has one.
  std::optional<std::string> error;
  double value = 0.0;  ///< z(0, 0) along the local z axis, mm (as VertexRule::value)
  double c = 0.0;      ///< curvature in the section in global orientation, 1/mm
  /// Magnitude of c for the afocal threshold, 1/mm: |c| for a rotationally symmetric surface;
  /// |h_xx| + 2 |h_xy| + |h_yy| of the local Hessian for an anamorphic one, because c is
  /// computed from all of it and rounds relative to it.
  double c_scale = 0.0;
};

/// ADR 0042, point 2, "Krümmung einer Fläche im Schnitt":
/// - rotationally symmetric in second order (is_anamorphic() false): vertex_rule(), and c in
///   global orientation exactly as first_order() takes it (sigma c), without the rounding of
///   Q Q^T;
/// - anamorphic: order >= 2 and slope = 0 as vertex_rule(); the global Hessian
///   H_g = sigma Q H_l Q^T must have aligned principal axes, |H_g,xy| <= kSymmetryTolerance
///   max(|H_g,xx|, |H_g,yy|); c = H_g,xx in section x, H_g,yy in section y.
[[nodiscard]] SectionRule section_rule(const geom::VertexExpansion<double>& e,
                                       const SectionFrame& frame,
                                       Section section);

/// The ideal cylinder lens in one section (ADR 0042, point 2, "Ideale Zylinderlinse").
struct CylinderLensSection {
  /// Why the lens has no model per section (ends with the ADR); empty if it has one.
  std::optional<std::string> error;
  /// b_s^2 with b = Q b_l: the power of the section is phi b_s^2, phi = |n|/f; exactly 0 for
  /// b_s^2 <= kSymmetryTolerance (the section without power, ADR 0042, Nachtrag #315 PR C1).
  double factor = 0.0;
};

/// The power tensor phi b b^T with b = Q b_l (no sigma: phi counts along the propagation and
/// b b^T is even in b); its principal axes must lie along x and y,
/// |phi b_x b_y| <= kSymmetryTolerance |phi|, i.e. |b_x b_y| <= kSymmetryTolerance.
/// @param power_axis b_l = CompiledIdealLens::power_axis, in the local frame of the surface
[[nodiscard]] CylinderLensSection cylinder_lens_section(const math::Vec3& power_axis,
                                                        const SectionFrame& frame,
                                                        Section section);

/// Global z of the vertex of a surface on the z axis (ADR 0041, point 2): origin_z +
/// sign(axis_z) value, with the same factor sign(axis_z) that turns c into global orientation.
/// For value == 0 it is origin_z bit for bit, also the sign of a zero.
/// @param origin_z global z of the local origin, mm
/// @param axis_z   global z component of the local z axis (+-1 on a rotationally symmetric path)
/// @param value    VertexRule::value, mm
[[nodiscard]] double vertex_z(double origin_z, double axis_z, double value);

}  // namespace rtt::paraxial::detail
