// Reports as data (#177, rtt/analysis/reports.hpp) from Python: the raytrace report, the system
// data report and the dimension report with their column classes. raytatouille.analysis wraps
// the functions (System or CompiledSystem, warnings); raytatouille.reports writes CSV.

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/variant.h>
#include <nanobind/stl/vector.h>

#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <vector>

#include "bindings.hpp"
#include "rtt/analysis/reports.hpp"
#include "rtt/paraxial/prescription.hpp"
#include "rtt/trace/ray_batch.hpp"

namespace nb = nanobind;
using namespace nb::literals;

namespace rtt::py {
namespace {

// Lists inside the reports, bound with columns().
struct RaytraceRows {
  std::vector<analysis::RaytraceRow> points;
};
struct DimensionSegments {
  std::vector<analysis::SegmentDimensions> points;
};

[[nodiscard]] std::uint8_t aperture_value(analysis::ApertureKind kind) {
  return static_cast<std::uint8_t>(kind);
}

/// Read-only NumPy bool copy of the coaxial flags (read_only_array needs vector::data(), which
/// std::vector<bool> lacks).
ReadOnlyArray<bool> coaxial_flags(const DimensionSegments& s) {
  const std::size_t n = s.points.size();
  auto data = std::make_unique<bool[]>(n);
  for (std::size_t i = 0; i < n; ++i) data[i] = s.points[i].coaxial;
  const std::size_t shape[1] = {n};
  const bool* values = data.get();
  // The capsule owns the copy from here on and frees it with the array.
  nb::capsule owner(data.get(), [](void* p) noexcept {
    const std::unique_ptr<bool[]> owned(static_cast<bool*>(p));
  });
  static_cast<void>(data.release());
  return ReadOnlyArray<bool>(values, 1, shape, owner);
}

}  // namespace

void bind_reports(nb::module_& m) {
  using namespace analysis;

  nb::enum_<ApertureKind>(m, "ApertureKind", nb::is_arithmetic(),
                          "Kind of the aperture a semi-diameter comes from; the aperture_* "
                          "columns of DimensionSegments hold these values as uint8.")
      .value("NONE", ApertureKind::None, "no aperture (semi-diameter NaN)")
      .value("CIRCULAR", ApertureKind::Circular, "radius (outer radius of an annulus)")
      .value("RECTANGULAR", ApertureKind::Rectangular, "half diagonal (circumscribed radius)")
      .value("ELLIPTICAL", ApertureKind::Elliptical, "major semi-axis (circumscribed radius)");

  auto rows = columns<RaytraceRows>(
      m, "RaytraceRows",
      "One entry per ray and slot (ray-major), up to the slot where each ray stopped "
      "(RayPaths.count); read-only NumPy copies. Slot 0 is the start, slot s the state after "
      "event s - 1. Positions in mm, directions as unit vectors; 'local' is the frame of the "
      "event's surface (NaN in slot 0).");
  column<std::uint64_t>(
      rows, "ray", [](const RaytraceRow& r) { return r.ray; },
      "Index of the ray in the start batch (copy).");
  column<std::uint64_t>(
      rows, "slot", [](const RaytraceRow& r) { return r.slot; },
      "0 = start, s = after event s - 1 (copy).");
  column<std::uint32_t>(
      rows, "surface", [](const RaytraceRow& r) { return r.surface; },
      "Surface of event s - 1, index into surface_ids; NO_SURFACE in slot 0 (copy).");
  column<double>(
      rows, "x", [](const RaytraceRow& r) { return r.x; }, "Position x, global, mm (copy).");
  column<double>(
      rows, "y", [](const RaytraceRow& r) { return r.y; }, "Position y, global, mm (copy).");
  column<double>(
      rows, "z", [](const RaytraceRow& r) { return r.z; }, "Position z, global, mm (copy).");
  column<double>(
      rows, "dx", [](const RaytraceRow& r) { return r.dx; }, "Direction x, global (copy).");
  column<double>(
      rows, "dy", [](const RaytraceRow& r) { return r.dy; }, "Direction y, global (copy).");
  column<double>(
      rows, "dz", [](const RaytraceRow& r) { return r.dz; }, "Direction z, global (copy).");
  column<double>(
      rows, "local_x", [](const RaytraceRow& r) { return r.local_x; },
      "Position x in the frame of the event's surface, mm; NaN in slot 0 (copy).");
  column<double>(
      rows, "local_y", [](const RaytraceRow& r) { return r.local_y; },
      "Position y in the frame of the event's surface, mm; NaN in slot 0 (copy).");
  column<double>(
      rows, "local_z", [](const RaytraceRow& r) { return r.local_z; },
      "Position z in the frame of the event's surface, mm; NaN in slot 0 (copy).");
  column<double>(
      rows, "local_dx", [](const RaytraceRow& r) { return r.local_dx; },
      "Direction x in the frame of the event's surface; NaN in slot 0 (copy).");
  column<double>(
      rows, "local_dy", [](const RaytraceRow& r) { return r.local_dy; },
      "Direction y in the frame of the event's surface; NaN in slot 0 (copy).");
  column<double>(
      rows, "local_dz", [](const RaytraceRow& r) { return r.local_dz; },
      "Direction z in the frame of the event's surface; NaN in slot 0 (copy).");
  column<double>(
      rows, "opl", [](const RaytraceRow& r) { return r.opl; },
      "Optical path length from the start, mm (copy).");
  column<double>(
      rows, "weight", [](const RaytraceRow& r) { return r.weight; },
      "Power, unpolarized source = 1 (copy).");
  column<std::uint8_t>(
      rows, "status", [](const RaytraceRow& r) { return status_value(r.status); },
      "RayStatus values as uint8 (copy).");

  nb::class_<RaytraceReport>(
      m, "RaytraceReport",
      "Raytrace report of a path (#177, rtt/analysis/reports.hpp): one row per ray and slot. "
      "The global position, direction, OPL, weight and status are bitwise those of RayPaths "
      "for the same rays (#80).")
      .def_prop_ro(
          "path", [](const RaytraceReport& r) { return r.path.index; }, "The traced path.")
      .def_prop_ro(
          "rows", [](const RaytraceReport& r) { return RaytraceRows{r.rows}; },
          "One entry per ray and slot, ray-major.")
      .def_ro("rays", &RaytraceReport::rays, "Rays of the start batch.")
      .def_ro("slots", &RaytraceReport::slots, "Events of the path + 1.");

  nb::class_<SystemReport>(
      m, "SystemReport",
      "System data of one path at one wavelength (#177): the system settings and the paraxial "
      "prescription.")
      .def_prop_ro(
          "path", [](const SystemReport& r) { return r.path.index; }, "The reported path.")
      .def_ro("wavelength", &SystemReport::wavelength,
              "Wavelength index of the prescription (index into wavelengths_um).")
      .def_prop_ro(
          "wavelengths_um",
          [](const SystemReport& r) {
            return read_only_array<double>(r.wavelengths_um, [](double w) { return w; });
          },
          nb::rv_policy::reference, "All system wavelengths, um (vacuum); read-only copy.")
      .def_ro("reference_wavelength", &SystemReport::reference_wavelength,
              "Index of the reference wavelength.")
      .def_ro("field_count", &SystemReport::field_count, "Number of field points.")
      .def_ro("surface_count", &SystemReport::surface_count, "Surfaces of the compiled system.")
      .def_ro("event_count", &SystemReport::event_count, "Events of the path.")
      .def_ro("stop", &SystemReport::stop,
              "Stop surface on the path (index into surface_ids); None without a stop.")
      .def_prop_ro(
          "prescription",
          [](const SystemReport& r) -> std::optional<paraxial::Prescription> {
            return r.prescription;
          },
          "Paraxial prescription (as raytatouille.paraxial.prescription); None if the path has "
          "no paraxial data, with the reason in warnings (report.paraxial_unavailable).")
      .def_ro("warnings", &SystemReport::warnings,
              "Warnings with codes (report.paraxial_unavailable, ADR 0023); "
              "raytatouille.analysis also issues them as RaytatouilleWarning.");

  auto segments = columns<DimensionSegments>(
      m, "DimensionSegments",
      "One entry per segment (the glass between surfaces j and j + 1) of every lens and plate, "
      "in element order; read-only NumPy copies, NaN where a dimension is undefined.");
  column<std::uint32_t>(
      segments, "element", [](const SegmentDimensions& s) { return s.element; },
      "Index into layout.elements() (copy).");
  column<std::uint32_t>(
      segments, "first_surface", [](const SegmentDimensions& s) { return s.first_surface; },
      "Surface j, index into surface_ids (copy).");
  segments.def_prop_ro(
      "coaxial", &coaxial_flags, nb::rv_policy::reference,
      "True if surface j + 1 shares the axis of surface j (z axes parallel or antiparallel "
      "within a tilt of 1e-12 rad, vertex on the axis within 1e-9 mm); otherwise the "
      "thicknesses are NaN (copy).");
  column<double>(
      segments, "centre_thickness", [](const SegmentDimensions& s) { return s.centre_thickness; },
      "Vertex distance along the z axis of surface j, mm (copy).");
  column<double>(
      segments, "semi_diameter_first",
      [](const SegmentDimensions& s) { return s.semi_diameter_first; },
      "Circumscribed radius of the aperture of surface j, mm; NaN without aperture (copy).");
  column<double>(
      segments, "semi_diameter_second",
      [](const SegmentDimensions& s) { return s.semi_diameter_second; },
      "Circumscribed radius of the aperture of surface j + 1, mm; NaN without aperture (copy).");
  column<std::uint8_t>(
      segments, "aperture_first",
      [](const SegmentDimensions& s) { return aperture_value(s.aperture_first); },
      "ApertureKind values of surface j as uint8 (copy).");
  column<std::uint8_t>(
      segments, "aperture_second",
      [](const SegmentDimensions& s) { return aperture_value(s.aperture_second); },
      "ApertureKind values of surface j + 1 as uint8 (copy).");
  column<double>(
      segments, "edge_thickness", [](const SegmentDimensions& s) { return s.edge_thickness; },
      "Distance along the z axis of surface j between the two surfaces at the meridional "
      "height y = h of surface j, h the larger semi-diameter, mm; NaN if a surface does not "
      "reach h or the segment is not coaxial (copy).");
  column<double>(
      segments, "diameter", [](const SegmentDimensions& s) { return s.diameter; },
      "2 h, mm; NaN if a surface has no aperture (copy).");

  nb::class_<DimensionReport>(
      m, "DimensionReport",
      "Dimension report (#177): thicknesses and diameters of every lens and plate, from the "
      "compiled geometry (relative placement and configurations included).")
      .def_prop_ro(
          "segments", [](const DimensionReport& r) { return DimensionSegments{r.segments}; },
          "One entry per segment, in element order.");

  // Functions; raytatouille.analysis wraps them.
  m.def(
      "raytrace_report",
      [](const compile::CompiledSystem& s, const PathArg& path, const trace::RayBatch& start,
         std::optional<int> threads) {
        const compile::PathId id = path_id(s, path);
        return released(threads, [&] { return raytrace_report(s, id, start); });
      },
      "system"_a, "path"_a, "start"_a, "threads"_a.none(),
      "Raytrace report of `path` for the start rays `start` (copied, not changed).");
  m.def(
      "system_report",
      [](const compile::CompiledSystem& s, const PathArg& path,
         std::optional<std::uint16_t> wavelength) {
        return system_report(s, path_id(s, path), wavelength_index(s, wavelength));
      },
      "system"_a, "path"_a, "wavelength"_a.none(),
      "System data report of `path` at `wavelength` (None: reference).");
  m.def(
      "dimension_report", [](const compile::CompiledSystem& s) { return dimension_report(s); },
      "system"_a, "Dimension report of all lenses and plates.");
}

}  // namespace rtt::py
