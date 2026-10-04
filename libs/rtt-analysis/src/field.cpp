#include "rtt/analysis/field.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <string>

#include "common.hpp"
#include "rtt/paraxial/paraxial.hpp"

namespace rtt::analysis {

namespace {

using compile::CompiledSystem;
using compile::PathId;
using detail::arrived;
using detail::check_path;
using detail::check_wavelength;
using detail::image_surface;
using detail::trace_aimed;

/// Largest radial field value of the system's field points (units of its field type).
double max_field(const CompiledSystem& system) {
  double m = 0.0;
  for (const model::Field& f : system.fields().points) m = std::max(m, std::hypot(f.x, f.y));
  return m;
}

/// Field values of the sweep: (0, f * max) for f = i / (samples - 1).
std::vector<model::Field> sweep(const CompiledSystem& system, int samples) {
  if (samples < 2) throw std::invalid_argument("analysis: a field sweep needs samples >= 2");
  const double m = max_field(system);
  if (!(m > 0.0)) throw std::invalid_argument("analysis: the system has no off-axis field point");
  std::vector<model::Field> fields;
  for (int i = 0; i < samples; ++i) {
    const double f = static_cast<double>(i) / (samples - 1);
    fields.push_back({0.0, f * m, 1.0});
  }
  return fields;
}

double fraction(const CompiledSystem& system, const model::Field& field) {
  const double m = max_field(system);
  return m > 0.0 ? std::hypot(field.x, field.y) / m : 0.0;
}

/// Paraxial chief ray of `field` in the plane z (global x, y offsets from the axis, mm).
/// Linear in the field value: a unit ray is traced with rtt-paraxial and scaled, exactly as the
/// field types are converted for ray aiming in rtt-trace (#8): unit slope through the entrance
/// pupil centre for an object at infinity, unit object height otherwise.
Point2 paraxial_chief(const CompiledSystem& system,
                      PathId path,
                      std::uint16_t wavelength,
                      const model::Field& field,
                      double z) {
  const paraxial::FirstOrder fo = paraxial::first_order(system, path, wavelength);
  if (!fo.entrance_pupil || !fo.entrance_pupil->z) {
    throw std::invalid_argument("analysis: the entrance pupil is not defined");
  }
  const double z_ep = *fo.entrance_pupil->z;
  const bool infinite = system.object().at_infinity;
  const double z_obj = infinite ? 0.0 : -system.object().distance.value;
  const std::vector<paraxial::RayAtEvent> unit =
      infinite ? paraxial::trace_ray(system, path, wavelength, z_ep, 0.0, 1.0)
               : paraxial::trace_ray(system, path, wavelength, z_obj, 1.0, -1.0 / (z_ep - z_obj));
  const auto height = [&unit](double plane) {
    const auto& last = unit.back();
    return last.y + (plane - last.z) * last.u;
  };
  double sx = 0.0;
  double sy = 0.0;
  switch (system.fields().type) {
    case model::FieldType::AngleDeg: {
      const double tx = std::tan(field.x * std::numbers::pi / 180.0);
      const double ty = std::tan(field.y * std::numbers::pi / 180.0);
      // Object at infinity: unit slope; finite object: object point on the chief ray through
      // the entrance pupil centre (as in rtt-trace).
      sx = infinite ? tx : (z_obj - z_ep) * tx;
      sy = infinite ? ty : (z_obj - z_ep) * ty;
      break;
    }
    case model::FieldType::ObjectHeight:
      sx = field.x;
      sy = field.y;
      break;
    case model::FieldType::ParaxialImageHeight: {
      if (!fo.image_z) {
        throw std::invalid_argument("analysis: paraxial image height needs a finite image");
      }
      const double unit_image = height(*fo.image_z);
      sx = field.x / unit_image;
      sy = field.y / unit_image;
      break;
    }
  }
  const double h = height(z);
  return {sx * h, sy * h};
}

/// Aims and traces one ray of `field` to the image surface; its global position and direction.
std::pair<math::Vec3, math::Vec3> trace_to_image(const CompiledSystem& system,
                                                 PathId path,
                                                 const model::Field& field,
                                                 std::uint16_t wavelength,
                                                 double px,
                                                 double py,
                                                 trace::Aiming aiming) {
  const trace::AimedRay aimed = trace::aim_ray(system, path, field, wavelength, px, py, aiming);
  const trace::RayBatch rays = trace_aimed(system, path, aimed, wavelength);
  if (!arrived(rays, 0, image_surface(system, path))) {
    throw AnalysisError("analysis: a ray at pupil (" + std::to_string(px) + ", " +
                        std::to_string(py) + ") does not reach the image surface");
  }
  return {math::Vec3(rays.pos_x()[0], rays.pos_y()[0], rays.pos_z()[0]),
          math::Vec3(rays.dir_x()[0], rays.dir_y()[0], rays.dir_z()[0])};
}

/// Midpoint of the shortest connection of the lines p1 + s d1 and p2 + t d2 (unit d).
math::Vec3 closest_point(const math::Vec3& p1,
                         const math::Vec3& d1,
                         const math::Vec3& p2,
                         const math::Vec3& d2) {
  const math::Vec3 w = p1 - p2;
  const double b = d1.dot(d2);
  const double d = d1.dot(w);
  const double e = d2.dot(w);
  const double den = 1.0 - b * b;
  if (!(den > 1e-30)) {
    throw AnalysisError("analysis: neighbour rays are parallel in image space (afocal)");
  }
  const double s = (b * e - d) / den;
  const double t = (e - b * d) / den;
  return 0.5 * (p1 + s * d1 + p2 + t * d2);
}

}  // namespace

DistortionPoint distortion_at(const compile::CompiledSystem& system,
                              compile::PathId path,
                              const model::Field& field,
                              std::uint16_t wavelength,
                              trace::Aiming aiming) {
  check_path(system, path);
  check_wavelength(system, wavelength);
  DistortionPoint p;
  p.field = field;
  p.fraction = fraction(system, field);
  const auto [pos, dir] = trace_to_image(system, path, field, wavelength, 0.0, 0.0, aiming);
  const math::Vec3 vertex = system.surfaces()[image_surface(system, path)].to_global.translation();
  const Point2 par = paraxial_chief(system, path, wavelength, field, vertex.z());
  // Heights along the field direction in the global x-y plane through the image vertex.
  const double r = std::hypot(field.x, field.y);
  const double ex = r > 0.0 ? field.x / r : 0.0;
  const double ey = r > 0.0 ? field.y / r : 1.0;
  p.real_height = (pos.x() - vertex.x()) * ex + (pos.y() - vertex.y()) * ey;
  p.paraxial_height = par.x * ex + par.y * ey;
  if (r == 0.0) {
    p.percent = 0.0;  // on axis by definition
  } else if (p.paraxial_height == 0.0) {
    throw AnalysisError(
        "analysis: the paraxial chief ray height is 0 at the image surface for an off-axis "
        "field (image surface in the exit-pupil plane); distortion is not defined");
  } else {
    p.percent = (p.real_height - p.paraxial_height) / p.paraxial_height * 100.0;
  }
  return p;
}

std::vector<DistortionPoint> distortion(const compile::CompiledSystem& system,
                                        compile::PathId path,
                                        std::uint16_t wavelength,
                                        const FieldSweepOptions& options) {
  check_path(system, path);
  check_wavelength(system, wavelength);
  std::vector<DistortionPoint> points;
  for (const model::Field& f : sweep(system, options.samples)) {
    points.push_back(distortion_at(system, path, f, wavelength, options.aiming));
  }
  return points;
}

FieldCurvaturePoint field_curvature_at(const compile::CompiledSystem& system,
                                       compile::PathId path,
                                       const model::Field& field,
                                       std::uint16_t wavelength,
                                       const FieldCurvatureOptions& options) {
  check_path(system, path);
  check_wavelength(system, wavelength);
  const double delta = options.delta;
  if (!(delta > 0.0 && delta < 1.0)) {
    throw std::invalid_argument("analysis: delta must lie in (0, 1)");
  }
  FieldCurvaturePoint p;
  p.field = field;
  p.fraction = fraction(system, field);
  // Tangential direction in the pupil: towards the field point (+y on axis); sagittal: normal.
  const double r = std::hypot(field.x, field.y);
  const double tx = r > 0.0 ? field.x / r : 0.0;
  const double ty = r > 0.0 ? field.y / r : 1.0;
  const auto [chief_pos, chief_dir] =
      trace_to_image(system, path, field, wavelength, 0.0, 0.0, options.aiming);
  const double sign = chief_dir.z() >= 0.0 ? 1.0 : -1.0;  // image-space propagation along z
  const math::Vec3 vertex = system.surfaces()[image_surface(system, path)].to_global.translation();
  const auto focus = [&](double ux, double uy) {
    const auto [p1, d1] =
        trace_to_image(system, path, field, wavelength, delta * ux, delta * uy, options.aiming);
    const auto [p2, d2] =
        trace_to_image(system, path, field, wavelength, -delta * ux, -delta * uy, options.aiming);
    return (closest_point(p1, d1, p2, d2).z() - vertex.z()) * sign;
  };
  p.tangential = focus(tx, ty);
  p.sagittal = focus(-ty, tx);
  p.astigmatism = p.tangential - p.sagittal;
  return p;
}

std::vector<FieldCurvaturePoint> field_curvature(const compile::CompiledSystem& system,
                                                 compile::PathId path,
                                                 std::uint16_t wavelength,
                                                 const FieldCurvatureOptions& options) {
  check_path(system, path);
  check_wavelength(system, wavelength);
  std::vector<FieldCurvaturePoint> points;
  for (const model::Field& f : sweep(system, options.samples)) {
    points.push_back(field_curvature_at(system, path, f, wavelength, options));
  }
  return points;
}

}  // namespace rtt::analysis
