#pragma once

/// @file reports.hpp
/// Reports as data (#177): a raytrace report (ray x event), a system data report (paraxial
/// prescription and system settings) and a dimension report (thicknesses and diameters per
/// element). No text tables and no plots in the library; Python writes CSV (#177, part B).

#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <string>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/model/validate.hpp"
#include "rtt/paraxial/prescription.hpp"
#include "rtt/trace/ray_batch.hpp"

namespace rtt::analysis {

// ------------------------------------------------------------------------- raytrace report ----

/// One row of the raytrace report: the state of one ray after one event of the path (slot
/// s >= 1) or at its start (slot 0). Positions in mm, directions as unit vectors; "local" is the
/// coordinate frame of the event's surface (CompiledSurface::to_local), NaN in slot 0 and after
/// the ray stopped (the slots beyond RayPaths::count).
struct RaytraceRow {
  std::size_t ray = 0;                        ///< index of the ray in the start batch
  std::size_t slot = 0;                       ///< 0 = start, s = after event s - 1
  std::uint32_t surface = trace::kNoSurface;  ///< surface of event s - 1; kNoSurface in slot 0
  double x = 0.0, y = 0.0, z = 0.0;           ///< position, global, mm
  double dx = 0.0, dy = 0.0, dz = 0.0;        ///< direction, global, unit
  double local_x = 0.0, local_y = 0.0, local_z = 0.0;     ///< position, surface frame, mm
  double local_dx = 0.0, local_dy = 0.0, local_dz = 0.0;  ///< direction, surface frame, unit
  double opl = 0.0;     ///< optical path length from the start, mm
  double weight = 0.0;  ///< power, unpolarized source = 1
  trace::RayStatus status = trace::RayStatus::Alive;
};

/// Raytrace report of a path: one row per ray and slot (ray-major), up to the slot where each
/// ray stopped (RayPaths::count). The global position, direction, OPL, weight and status are
/// bitwise those of trace::RayPaths for the same rays (#80).
struct RaytraceReport {
  compile::PathId path{0};
  std::vector<RaytraceRow> rows;
  std::size_t rays = 0;   ///< rays of the start batch
  std::size_t slots = 0;  ///< events of the path + 1
};

/// Traces a copy of `start` along `path` with path recording (trace::SequentialTracer, #80) and
/// returns the report; every ray of `start` is recorded (about 65 bytes per ray and slot, so
/// the report is meant for a few rays, e.g. chief and marginal rays).
/// @throws std::invalid_argument for an unknown path or an empty start batch
[[nodiscard]] RaytraceReport raytrace_report(const compile::CompiledSystem& system,
                                             compile::PathId path,
                                             const trace::RayBatch& start);

// --------------------------------------------------------------------- system data report ----

/// System data of one path at one wavelength: the system settings and the paraxial
/// prescription (first-order data, marginal and chief ray per event).
struct SystemReport {
  compile::PathId path{0};
  std::uint16_t wavelength = 0;        ///< index into CompiledSystem::wavelengths_um()
  std::vector<double> wavelengths_um;  ///< all system wavelengths, um (vacuum)
  std::uint16_t reference_wavelength = 0;
  std::size_t field_count = 0;
  std::size_t surface_count = 0;      ///< surfaces of the compiled system
  std::size_t event_count = 0;        ///< events of the path
  std::optional<std::uint32_t> stop;  ///< stop surface on the path, if any
  /// Paraxial data; none if the path has none (not rotationally symmetric, crystal, order != 0),
  /// with the reason in `warnings` (report.paraxial_unavailable).
  std::optional<paraxial::Prescription> prescription;
  std::vector<model::Diagnostic> warnings;
};

/// System data report of `path` at `wavelength`.
/// @throws std::invalid_argument for an unknown path or wavelength index
[[nodiscard]] SystemReport system_report(const compile::CompiledSystem& system,
                                         compile::PathId path,
                                         std::uint16_t wavelength);

// ----------------------------------------------------------------------- dimension report ----

/// Kind of the aperture a semi-diameter comes from.
enum class ApertureKind : std::uint8_t { None, Circular, Rectangular, Elliptical };

/// Dimensions of one segment of an element: the glass between surfaces j and j + 1.
struct SegmentDimensions {
  std::uint32_t element = 0;        ///< index into CompiledSystem::elements()
  std::uint32_t first_surface = 0;  ///< surface j (index into CompiledSystem::surfaces())
  /// True if surface j + 1 shares the axis of surface j (z axes parallel or antiparallel within
  /// 1e-12, vertex on the axis within 1e-9 mm); otherwise the thicknesses are NaN.
  bool coaxial = true;
  /// Centre thickness: distance from the vertex of surface j to that of surface j + 1 along the
  /// z axis of surface j, mm (from the compiled geometry, so relative placement and
  /// configurations are included).
  double centre_thickness = 0.0;
  /// Semi-diameters of surfaces j and j + 1, mm: the radius of a circular aperture (the outer
  /// radius of an annulus), the half diagonal of a rectangle, the major semi-axis of an ellipse
  /// (the circumscribed radius); NaN without aperture.
  double semi_diameter_first = 0.0;
  double semi_diameter_second = 0.0;
  ApertureKind aperture_first = ApertureKind::None;
  ApertureKind aperture_second = ApertureKind::None;
  /// Edge thickness at the larger of the two semi-diameters h (where the body ends; a mechanical
  /// diameter is not part of the model): the distance along the z axis of surface j between the
  /// two surfaces at the meridional height y = h of surface j, mm; NaN if a surface does not
  /// reach h (domain of the shape) or the segment is not coaxial.
  double edge_thickness = 0.0;
  double diameter = 0.0;  ///< 2 h, mm
};

/// Dimension report: one entry per segment of every lens and plate, in element order.
struct DimensionReport {
  std::vector<SegmentDimensions> segments;
};

/// Dimensions of all lenses and plates of `system` (elements with at least two surfaces).
[[nodiscard]] DimensionReport dimension_report(const compile::CompiledSystem& system);

}  // namespace rtt::analysis
