#pragma once

/// @file layout.hpp
/// Geometry export for drawing a system (#81, docs/architecture.md, Layout): sag and normals of
/// the compiled surfaces, profile polylines in a section plane and closed outlines of the glass
/// segments of an element. Lengths in mm; "local" is the coordinate system of a surface
/// (CompiledSurface::to_global), "global" the system coordinates (right-handed, optical axis
/// +z).

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/math/types.hpp"

namespace rtt::compile {

/// Coordinate frame of results.
enum class Frame : std::uint8_t {
  Local,   ///< coordinates of the surface
  Global,  ///< system coordinates
};

/// A section plane in global coordinates: the points p with (p - point) . normal = 0.
struct SectionPlane {
  math::Vec3 point = math::Vec3::Zero();    ///< a point of the plane, mm
  math::Vec3 normal = math::Vec3::UnitX();  ///< plane normal, any length > 0
};

/// A polyline of points in global coordinates, mm.
using Polyline = std::vector<math::Vec3>;

/// Sag z = shape(x, y) of `surface` at the local points (x[i], y[i]), mm, written to z[i]; NaN
/// outside the domain of the shape (rtt-geom). The aperture is not applied.
/// @throws std::out_of_range if `surface` is not a surface index
/// @throws std::invalid_argument if x, y and z differ in length
void surface_sag(const CompiledSystem& system,
                 std::uint32_t surface,
                 std::span<const double> x,
                 std::span<const double> y,
                 std::span<double> z);

/// Unit normal of `surface` at the local points (x[i], y[i]): n = (-dz/dx, -dz/dy, 1) / |...|,
/// pointing to +z at the vertex as in rtt-geom; in local coordinates or rotated into global
/// coordinates. Written to normals[3 i + c]; NaN outside the domain of the shape.
/// @throws std::out_of_range if `surface` is not a surface index
/// @throws std::invalid_argument if the lengths do not match (normals: 3 per point)
void surface_normal(const CompiledSystem& system,
                    std::uint32_t surface,
                    std::span<const double> x,
                    std::span<const double> y,
                    Frame frame,
                    std::span<double> normals);

/// Profile of `surface` in a section plane: the curve z = sag along the line in which the plane
/// cuts the local x-y plane, clipped to the aperture (an annulus gives two pieces) and to the
/// domain of the shape, in global coordinates. The plane must be parallel to the local z axis of
/// the surface (|n . z_local| <= 1e-12 for the unit normal n), which covers meridional sections,
/// tilts within the section plane, fold mirrors and planes offset from the axis; other planes
/// need a general surface section (with the 3D meshes, M10). Each piece has `samples` points,
/// equally spaced along the line in the direction t = z_local x n (local), from one end to the
/// other; ends on the aperture edge are exact. Where the domain of the shape (max_radius, e.g.
/// the hemisphere of a sphere) limits the piece, it ends 1e-12 max_radius inside it, since the
/// sag is NaN beyond and has an infinite slope there. No piece if the line misses the aperture.
/// A surface without aperture is cut at the domain of its shape, which is rarely what a drawing
/// wants: give lens surfaces apertures (semi-diameters).
/// @throws std::out_of_range if `surface` is not a surface index
/// @throws std::invalid_argument if samples < 2, the plane normal is zero or not finite, the
///         plane is not parallel to the local z axis, or the surface has neither an aperture nor
///         a bounded shape
[[nodiscard]] std::vector<Polyline> surface_profile(const CompiledSystem& system,
                                                    std::uint32_t surface,
                                                    const SectionPlane& plane,
                                                    std::size_t samples);

/// Closed outlines of the glass segments of `element` in a section plane, global coordinates:
/// for segment j between surfaces j and j + 1, the profile of surface j, the edge to surface
/// j + 1, the profile of surface j + 1 backwards and the edge back; the last point repeats the
/// first. Both profiles are taken along the section direction of surface j (a turned-over
/// surface j + 1 is reversed). Edges, in the section coordinates of surface j: if the two rims
/// end at different positions along the section line, the edge is a step, parallel to the axis
/// of surface j at the outermost of the two rims (as seen from the glass), then along the
/// section line to the other rim. At an outer rim this is the cylindrical edge at the larger
/// radius with a flat shoulder, at a central hole the bore at the smaller radius. An annulus on
/// both surfaces gives one outline per side of the hole.
///
/// Only lenses and segmented plates, whose surfaces bound the segments pairwise, have outlines.
/// All plates of one material (windows, prisms, rhombs), mirrors, thin elements, stops and
/// detectors give an empty list; draw their surface profiles. Pieces that the section plane
/// misses on both surfaces of a segment give no outline.
/// @throws std::out_of_range if `element` is not an element index
/// @throws std::invalid_argument as surface_profile(), or if the two surfaces of a segment give
///         a different number of pieces (annulus on one surface only, or the plane misses one
///         surface)
[[nodiscard]] std::vector<Polyline> element_outlines(const CompiledSystem& system,
                                                     std::uint32_t element,
                                                     const SectionPlane& plane,
                                                     std::size_t samples);

}  // namespace rtt::compile
