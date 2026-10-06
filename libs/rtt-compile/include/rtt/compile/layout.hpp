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
/// the surface (|n . z_local| < 1e-12 |n|), which covers meridional sections, tilts within the
/// section plane, fold mirrors and planes offset from the axis; other planes need a general
/// surface section (with the 3D meshes, M10). Each piece has `samples` points, equally spaced
/// along the line, from one aperture edge to the other in the direction t = z_local x n
/// (local); the end points lie on the aperture edge. No piece if the line misses the aperture.
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
/// first. Edges between rims of different radius are a step: parallel to the axis of surface j
/// at the larger radius up to the axial position of the smaller rim, then along the section
/// line to it (as lens drawings show a cylindrical edge and a flat shoulder). An annulus on both
/// surfaces gives one outline per side of the hole.
///
/// Only elements whose surfaces bound the segments pairwise: lenses and segmented plates. Plates,
/// prisms and rhombs of one material, mirrors, thin elements, stops and detectors give no
/// outline (empty list); draw their surface profiles. Pieces that the section plane misses on
/// both surfaces of a segment give no outline.
/// @throws std::out_of_range if `element` is not an element index
/// @throws std::invalid_argument as surface_profile(), or if the two surfaces of a segment give
///         a different number of pieces (annulus on one surface only, or the plane misses one
///         surface)
[[nodiscard]] std::vector<Polyline> element_outlines(const CompiledSystem& system,
                                                     std::uint32_t element,
                                                     const SectionPlane& plane,
                                                     std::size_t samples);

}  // namespace rtt::compile
