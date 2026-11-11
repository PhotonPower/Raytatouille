#pragma once

/// @file irradiance.hpp
/// Irradiance map on a surface from traced rays (#379, ADR 0048): the hit points of the rays on
/// one surface, binned into an nx x ny grid of the local x-y plane of that surface.
///
/// Conventions (decisions of the maintainer in #379, ADR 0048):
/// - **Incident power:** a ray counts with the power that reaches the surface, at its hit point:
///   its weight after the path to the surface (with the volume absorption of the medium before
///   it, as the trace applies it before the event) and before the event. What happens at the
///   surface (Fresnel, coating, absorber, total internal reflection) does not change its
///   contribution.
/// - **Which rays count:** those that hit the surface inside its aperture: status after the
///   event Alive, Tir, Absorbed, EventImpossible or Evanescent (all carry the hit point,
///   rtt/trace/ray_paths.hpp). Not counted: Vignetted (hit outside the aperture), Missed and
///   NoConvergence (no hit), rays lost before the surface, and hits outside the grid; each group
///   is reported with its number and weight.
/// - **Cell area:** the area of the cell in the local x-y plane (dx dy), i.e. the projection of
///   the surface onto that plane. On a curved surface the true area of a cell is larger by the
///   factor 1/cos of the surface slope; that correction is not applied (a later extension).
/// - **Units:** positions and edges in mm, local surface coordinates; values in weight per mm^2.
///   The weight is the power of a ray relative to its start weight (1 for an unpolarized ray
///   from a source, ADR 0021), so a value is the incident power per mm^2 in the units of the
///   start weights.
/// - **Binning:** cell ix covers [x_ix, x_ix+1) with x_i = x_min + i (x_max - x_min) / nx, the
///   last cell includes x_max; likewise in y. The index is checked against these edges, so a
///   point on an edge falls into the cell to its right. Sums run over the rays in batch order,
///   so the map does not depend on the number of threads or the block size of the trace
///   (ADR 0004).

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/math/isometry.hpp"
#include "rtt/model/surface.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/ray_paths.hpp"

namespace rtt::analysis {

/// Rectangle of the grid in local surface coordinates, mm.
struct IrradianceExtent {
  double x_min = 0.0;
  double x_max = 0.0;
  double y_min = 0.0;
  double y_max = 0.0;
};

/// Options of irradiance_map().
struct IrradianceOptions {
  std::size_t nx = 64;  ///< cells in local x, >= 1
  std::size_t ny = 64;  ///< cells in local y, >= 1
  /// The grid; none: the bounding rectangle of the surface aperture (a circle or annulus of
  /// radius r: [-r, r]^2; a rectangle: [-hx, hx] x [-hy, hy]; an ellipse: [-a, a] x [-b, b]).
  std::optional<IrradianceExtent> extent;
  /// RayPaths only: index of the path event at the surface (0 = first event of the path); none:
  /// every event of the path at the surface, summed. Ignored by the RayBatch overload.
  std::optional<std::size_t> event;
  /// RayPaths only: system wavelength index of the recorded rays. Needed if the medium before a
  /// counted event absorbs (kappa != 0 at any system wavelength): the incident weight then
  /// includes exp(-4 pi kappa t / lambda) over the last path t, since the RayPaths slots hold
  /// the state at the previous surface. Ignored otherwise and by the RayBatch overload.
  std::optional<std::uint16_t> wavelength;
};

/// Rays that did not count, with their summed incident weight (at the surface; for rays lost
/// before it, the weight before the event where they were lost). Without a chosen event a ray
/// counts once per event at the surface, here too.
struct IrradianceLosses {
  std::size_t count = 0;
  double weight = 0.0;
};

/// Result of irradiance_map() (ADR 0048).
struct IrradianceMap {
  std::size_t nx = 0;
  std::size_t ny = 0;
  /// Incident weight per mm^2 of each cell, row-major [iy * nx + ix] (ny rows of nx values).
  std::vector<double> values;
  std::vector<double> x_edges;  ///< nx + 1 cell edges in local x, mm, increasing
  std::vector<double> y_edges;  ///< ny + 1 cell edges in local y, mm, increasing
  std::uint32_t surface = 0;    ///< surface index in the compiled system
  model::SurfaceId surface_id;  ///< stable id, e.g. "IMG"
  std::string pointer;          ///< JSON pointer of the surface in the system file
  /// Local surface coordinates -> global (CompiledSurface::to_global), to place the cells on the
  /// geometry.
  math::Isometry3 to_global;
  std::vector<std::size_t> events;  ///< path events that were counted, increasing
  std::size_t rays_counted = 0;     ///< hits inside the aperture and the grid (over all events)
  double weight_counted = 0.0;      ///< their incident weight; sum(values) dx dy, up to rounding
  IrradianceLosses outside_grid;    ///< hits inside the aperture but outside the grid
  IrradianceLosses vignetted;       ///< hits outside the aperture (status Vignetted)
  IrradianceLosses not_reached;     ///< rays lost before the surface, or Missed/NoConvergence
};

/// Irradiance map on `surface` from recorded ray paths (any surface of the path).
///
/// A ray contributes once per counted event: for event k at the hit point of slot k + 1, with
/// the weight of slot k times the volume absorption over the path from slot k to slot k + 1
/// (rtt/trace/ray_paths.hpp; Byrnes, arXiv:1603.02720v5, Eqs. (1), (2), as move_to_hit).
/// @param system  the compiled system the paths were traced in
/// @param path    index of the compiled path the paths were traced on
/// @param paths   recorded paths (SequentialTracer::trace with RayPaths)
/// @param surface id of the surface, e.g. "IMG"
/// @param options grid, extent, event and wavelength
/// @throws std::invalid_argument for an invalid path or paths that do not belong to it, if the
///         surface is unknown or not an event of the path, if `options.event` is not an event
///         at the surface, if a counted event follows an absorbing medium and
///         `options.wavelength` is missing or invalid, for the grid errors (nx or ny 0; an
///         extent that is empty, not finite or so wide that a cell width overflows), or if there
///         is no extent and the surface has no aperture
[[nodiscard]] IrradianceMap irradiance_map(const compile::CompiledSystem& system,
                                           std::size_t path,
                                           const trace::RayPaths& paths,
                                           const model::SurfaceId& surface,
                                           const IrradianceOptions& options = {});

/// Irradiance map on the last surface of path `path` from the batch after the trace, without a
/// recording, so for any number of rays (decision of #379).
///
/// The batch holds the state after the last event, so its weight is the incident weight only if
/// that event cannot change it. Allowed: a Transmit event of order 0 at a surface with the
/// interaction Fresnel or IdealAntiReflection between media of the same index (a detector or a
/// dummy surface), without an ideal lens and without a diffraction efficiency; the weight then
/// stays the same up to rounding (a few ulp on a tilted surface, from the direction through the
/// local frame). Otherwise use the RayPaths overload. Counted: rays with status Alive; Vignetted
/// with that surface as `last_surface` is reported as `vignetted` (also a vignetting at an
/// earlier event at the same surface), every other status as `not_reached`.
/// @pre the trace ran to its end (not cancelled): rays of blocks that never started are still
///      Alive at their start and would count
/// @param system the compiled system the batch was traced in
/// @param path   index of the compiled path
/// @param rays   the batch after SequentialTracer::trace on that path
/// @param options grid and extent (`event` is ignored)
/// @throws std::invalid_argument for an invalid path, if the last event may change the weight
///         (see above), or for the grid errors of the other overload
[[nodiscard]] IrradianceMap irradiance_map(const compile::CompiledSystem& system,
                                           std::size_t path,
                                           const trace::RayBatch& rays,
                                           const IrradianceOptions& options = {});

}  // namespace rtt::analysis
