#include "paraxial_rays.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <string>
#include <string_view>
#include <variant>

namespace rtt::paraxial::detail {
namespace {

using compile::CompiledSystem;
using compile::PathId;

std::string field_error(std::string_view caller, const std::string& what) {
  return "paraxial: " + std::string(caller) + ": " + what;
}

}  // namespace

/// The paraxial trace (paraxial.cpp) uses c in global orientation: 1/R times the sign of the
/// local z axis along global z. The fourth-order sag c^3 r^4 / 8 (K + 1) + A4 r^4 changes sign
/// with the axis as well, so A4 flips with c while K, a ratio, stays.
ShapeData shape_data(const compile::CompiledSurface& s) {
  ShapeData d;
  const auto [c_local, k] =
      std::visit([](const auto& shape) { return shape.base_conic(); }, s.shape);
  double a4_local = 0.0;
  if (const auto* asphere = std::get_if<geom::EvenAsphere<double>>(&s.shape)) {
    if (!asphere->coefficients().empty()) a4_local = asphere->coefficients().front();
  }
  const double sign = s.to_global.apply_vector(math::Vec3::UnitZ()).z() > 0.0 ? 1.0 : -1.0;
  d.c = sign * c_local;
  d.conic = k;
  d.a4 = sign * a4_local;
  return d;
}

RayStart marginal_start(bool object_at_infinity, double z_obj, double z_ep, double r_ep) {
  // Axial object point to the rim of the entrance pupil at +y (as in seidel() before #84).
  return object_at_infinity ? RayStart{z_ep, r_ep, 0.0}
                            : RayStart{z_obj, 0.0, r_ep / (z_ep - z_obj)};
}

RayStart chief_start(const CompiledSystem& system,
                     PathId path,
                     std::uint16_t wavelength,
                     const FirstOrder& fo,
                     double z_ep,
                     double z_obj,
                     std::string_view caller) {
  const bool infinite = system.object().at_infinity;
  const auto& points = system.fields().points;
  // Largest radial field value: tan theta for angles (chief ray d ~ (tan theta_x, tan theta_y,
  // 1), docs/architecture.md "Feldwinkel und Pupille"), mm for heights.
  double value = 0.0;
  for (std::size_t i = 0; i < points.size(); ++i) {
    const model::Field& f = points[i];
    double radial = std::hypot(f.x, f.y);
    if (system.fields().type == model::FieldType::AngleDeg) {
      if (!(std::abs(f.x) < 90.0 && std::abs(f.y) < 90.0)) {
        // Place in the system file (ADR 0022, #98).
        throw ParaxialError(field_error(caller, "field angles must lie in (-90, 90) degree"),
                            "/fields/points/" + std::to_string(i));
      }
      radial = std::hypot(std::tan(f.x * std::numbers::pi / 180.0),
                          std::tan(f.y * std::numbers::pi / 180.0));
    }
    value = std::max(value, radial);
  }

  // Ray through the EP centre with unit field value: unit slope (angle; also for a finite
  // object, whose object point then lies on this ray as in rtt-trace, #8) or unit object height.
  const auto unit_ray = [&](bool angle) -> RayStart {
    if (angle) return {z_ep, 0.0, 1.0};
    return {z_obj, 1.0, -1.0 / (z_ep - z_obj)};
  };
  switch (system.fields().type) {
    case model::FieldType::AngleDeg: {
      const RayStart r = unit_ray(true);
      return {r.z, 0.0, value};
    }
    case model::FieldType::ObjectHeight: {
      if (infinite) {
        throw ParaxialError(field_error(caller, "an object height needs a finite object distance"),
                            "/fields/type");
      }
      const RayStart r = unit_ray(false);
      return {r.z, value, value * r.u};
    }
    case model::FieldType::ParaxialImageHeight: {
      // Linear in the field value: scale the unit chief ray to the requested paraxial image
      // height (as rtt-trace does for #8).
      if (!fo.image_z) {
        throw ParaxialError(
            field_error(caller, "a paraxial image height needs a finite paraxial image"));
      }
      const RayStart r = unit_ray(infinite);
      const auto ray = trace_ray(system, path, wavelength, r.z, r.y, r.u);
      const RayAtEvent& last = ray.back();
      const double unit = last.y + (*fo.image_z - last.z) * last.u;
      if (unit == 0.0)
        throw ParaxialError(field_error(caller, "the chief ray does not reach the image"));
      const double scale = value / unit;
      return {r.z, scale * r.y, scale * r.u};
    }
  }
  return {};
}

}  // namespace rtt::paraxial::detail
