// Reports as data (#177): raytrace, system data and dimensions.

#include "rtt/analysis/reports.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <variant>
#include <vector>

#include "rtt/compile/errors.hpp"
#include "rtt/compile/layout.hpp"
#include "rtt/diagnostics/codes.hpp"
#include "rtt/paraxial/paraxial.hpp"
#include "rtt/trace/ray_paths.hpp"
#include "rtt/trace/sequential.hpp"

namespace rtt::analysis {
namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

const compile::CompiledPath& checked_path(const compile::CompiledSystem& system,
                                          compile::PathId path,
                                          const char* who) {
  if (path.index >= system.paths().size()) {
    throw std::invalid_argument(std::string(who) + ": path " + std::to_string(path.index) +
                                " does not exist");
  }
  return system.path(path);
}

/// Semi-diameter of an aperture and its kind: the circumscribed radius (#177, decision of the
/// coordinator): circle (outer radius of an annulus), half diagonal of a rectangle, major
/// semi-axis of an ellipse; NaN without aperture.
std::pair<double, ApertureKind> semi_diameter(const std::optional<model::Aperture>& aperture) {
  if (!aperture) return {kNaN, ApertureKind::None};
  if (const auto* c = std::get_if<model::CircularAperture>(&*aperture)) {
    return {c->radius, ApertureKind::Circular};
  }
  if (const auto* r = std::get_if<model::RectangularAperture>(&*aperture)) {
    return {std::hypot(r->half_width_x, r->half_width_y), ApertureKind::Rectangular};
  }
  const auto& e = std::get<model::EllipticalAperture>(*aperture);
  return {std::max(e.semi_axis_x, e.semi_axis_y), ApertureKind::Elliptical};
}

double sag(const compile::CompiledSystem& system, std::uint32_t surface, double x, double y) {
  const double xs[] = {x};
  const double ys[] = {y};
  double z[] = {0.0};
  compile::surface_sag(system, surface, xs, ys, z);
  return z[0];
}

/// Tolerances of the coaxiality test (documented at SegmentDimensions::coaxial).
constexpr double kAxisTolerance = 1e-12;   // 1 - |z_j . z_k|
constexpr double kVertexTolerance = 1e-9;  // mm, lateral offset of vertex k from the axis of j

SegmentDimensions segment(const compile::CompiledSystem& system,
                          std::uint32_t element,
                          std::uint32_t j) {
  const compile::CompiledSurface& first = system.surfaces()[j];
  const compile::CompiledSurface& second = system.surfaces()[j + 1];
  SegmentDimensions d;
  d.element = element;
  d.first_surface = j;
  std::tie(d.semi_diameter_first, d.aperture_first) = semi_diameter(first.aperture);
  std::tie(d.semi_diameter_second, d.aperture_second) = semi_diameter(second.aperture);
  // The larger semi-diameter; undefined (NaN) if a surface has no aperture.
  const double h = std::isnan(d.semi_diameter_first) || std::isnan(d.semi_diameter_second)
                       ? kNaN
                       : std::max(d.semi_diameter_first, d.semi_diameter_second);
  d.diameter = 2.0 * h;

  // Vertex and axis of the second surface in the frame of the first.
  const math::Vec3 vertex =
      first.to_local.apply_point(second.to_global.apply_point(math::Vec3::Zero()));
  const math::Vec3 axis =
      first.to_local.apply_vector(second.to_global.apply_vector(math::Vec3::UnitZ()));
  d.coaxial = 1.0 - std::abs(axis.z()) <= kAxisTolerance &&
              std::hypot(vertex.x(), vertex.y()) <= kVertexTolerance;
  if (!d.coaxial) {
    d.centre_thickness = kNaN;
    d.edge_thickness = kNaN;
    return d;
  }
  d.centre_thickness = vertex.z();
  // Edge at the meridional height y = h of the first surface: the line through (0, h) parallel
  // to the common axis meets the second surface at its local (x', y'); both rim points are
  // compared along the axis of the first surface.
  const double z1 = sag(system, j, 0.0, h);
  const math::Vec3 foot =
      second.to_local.apply_point(first.to_global.apply_point(math::Vec3(0.0, h, 0.0)));
  const double z2_local = sag(system, j + 1, foot.x(), foot.y());
  const math::Vec3 rim2 = first.to_local.apply_point(
      second.to_global.apply_point(math::Vec3(foot.x(), foot.y(), z2_local)));
  d.edge_thickness = rim2.z() - z1;  // NaN if h is NaN or outside a shape's domain
  return d;
}

}  // namespace

RaytraceReport raytrace_report(const compile::CompiledSystem& system,
                               compile::PathId path,
                               const trace::RayBatch& start) {
  static_cast<void>(checked_path(system, path, "raytrace_report"));
  if (start.size() == 0) throw std::invalid_argument("raytrace_report: no start rays");
  trace::RayBatch rays = start;
  trace::RayPaths paths;
  static_cast<void>(
      trace::SequentialTracer{}.trace(system, path, rays, paths, std::nullopt, start.size()));

  RaytraceReport report;
  report.path = path;
  report.rays = start.size();
  report.slots = paths.slots;
  for (std::size_t r = 0; r < paths.ray_count(); ++r) {
    for (std::size_t s = 0; s < paths.count[r]; ++s) {
      const std::size_t i = r * paths.slots + s;
      const math::Vec3 p = paths.position_at(r, s);
      const math::Vec3 d = paths.direction_at(r, s);
      RaytraceRow row;
      row.ray = paths.ray_indices[r];
      row.slot = s;
      row.x = p.x();
      row.y = p.y();
      row.z = p.z();
      row.dx = d.x();
      row.dy = d.y();
      row.dz = d.z();
      row.opl = paths.opl[i];
      row.weight = paths.weight[i];
      row.status = paths.status[i];
      if (s == 0) {
        row.local_x = row.local_y = row.local_z = kNaN;
        row.local_dx = row.local_dy = row.local_dz = kNaN;
      } else {
        row.surface = paths.event_surfaces[s - 1];
        const compile::CompiledSurface& surface = system.surfaces()[row.surface];
        const math::Vec3 lp = surface.to_local.apply_point(p);
        const math::Vec3 ld = surface.to_local.apply_vector(d);
        row.local_x = lp.x();
        row.local_y = lp.y();
        row.local_z = lp.z();
        row.local_dx = ld.x();
        row.local_dy = ld.y();
        row.local_dz = ld.z();
      }
      report.rows.push_back(row);
    }
  }
  return report;
}

SystemReport system_report(const compile::CompiledSystem& system,
                           compile::PathId path,
                           std::uint16_t wavelength) {
  const compile::CompiledPath& p = checked_path(system, path, "system_report");
  if (wavelength >= system.wavelengths_um().size()) {
    throw std::invalid_argument("system_report: wavelength " + std::to_string(wavelength) +
                                " does not exist");
  }
  SystemReport r;
  r.path = path;
  r.wavelength = wavelength;
  r.wavelengths_um = system.wavelengths_um();
  r.reference_wavelength = system.reference_wavelength();
  r.field_count = system.fields().points.size();
  r.surface_count = system.surfaces().size();
  r.event_count = p.events.size();
  for (const compile::CompiledEvent& e : p.events) {
    if (system.surfaces()[e.surface].element_kind == model::ElementKind::Stop) {
      r.stop = e.surface;
      break;
    }
  }
  try {
    r.prescription = paraxial::prescription(system, path, wavelength);
  } catch (const paraxial::ParaxialError& e) {
    constexpr diagnostics::DiagnosticCode kCode = "report.paraxial_unavailable";
    r.warnings.push_back({kCode.severity(), "/paths/" + std::to_string(path.index),
                          std::string("no paraxial data: ") + e.what(), std::string(kCode.str())});
  }
  return r;
}

DimensionReport dimension_report(const compile::CompiledSystem& system) {
  DimensionReport report;
  for (std::uint32_t e = 0; e < system.elements().size(); ++e) {
    const compile::CompiledElement& element = system.elements()[e];
    if (element.kind != model::ElementKind::Lens && element.kind != model::ElementKind::Plate) {
      continue;
    }
    for (std::uint32_t k = 0; k + 1 < element.surface_count; ++k) {
      report.segments.push_back(segment(system, e, element.first_surface + k));
    }
  }
  return report;
}

}  // namespace rtt::analysis
