#include "paraxial_rays.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <optional>
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

bool reference_pupil_usable(const CompiledSystem& system,
                            PathId path,
                            std::uint16_t wavelength,
                            double z_obj) {
  if (wavelength == system.reference_wavelength()) return true;
  const model::FieldType type = system.fields().type;
  const bool infinite = system.object().at_infinity;
  const bool image_height = type == model::FieldType::ParaxialImageHeight;
  if (!image_height && !(type == model::FieldType::AngleDeg && !infinite)) return true;
  const FirstOrder fo_ref = first_order(system, path, system.reference_wavelength());
  const auto& ep = fo_ref.entrance_pupil;
  if (!ep) return false;
  // At infinity only a paraxial image height with a finite object works (parallel chief ray).
  if (!ep->z) return image_height && !infinite;
  return infinite || *ep->z != z_obj;
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

  // Field conversion at the reference wavelength, as in rtt-trace (#31, #50; fixed here for
  // #35): a paraxial image height, and a field angle with a finite object, depend on the
  // wavelength through the paraxial chief ray and the entrance pupil. They are converted at the
  // reference wavelength into a slope (object at infinity) or an object height (finite object);
  // the chief ray at `wavelength` then passes through the centre of its own entrance pupil
  // z_ep. At the reference wavelength the expressions are the earlier ones (bitwise the same).
  const std::uint16_t ref = system.reference_wavelength();
  const bool at_ref = wavelength == ref;
  std::optional<FirstOrder> fo_ref_storage;
  const auto reference_first_order = [&]() -> const FirstOrder& {
    if (at_ref) return fo;
    if (!fo_ref_storage) fo_ref_storage.emplace(first_order(system, path, ref));
    return *fo_ref_storage;
  };
  // Ray through the EP centre at z_ep_at with unit field value: unit slope (angle) or unit
  // object height (finite object).
  const auto unit_ray = [&](bool angle, double z_ep_at) -> RayStart {
    if (angle) return {z_ep_at, 0.0, 1.0};
    return {z_obj, 1.0, -1.0 / (z_ep_at - z_obj)};
  };
  switch (system.fields().type) {
    case model::FieldType::AngleDeg: {
      const RayStart r = unit_ray(true, z_ep);
      if (infinite || at_ref) return {r.z, 0.0, value};
      // Finite object: the object point lies on the chief ray through the entrance pupil at the
      // reference wavelength, h = (z_obj - z_EP,ref) tan theta (as rtt-trace, #31); the chief
      // ray at `wavelength` runs from it through z_ep.
      const auto& ep_ref = reference_first_order().entrance_pupil;
      if (!ep_ref || !ep_ref->z) {
        throw ParaxialError(field_error(caller,
                                        "a field angle with a finite object needs a finite "
                                        "entrance pupil at the reference wavelength"),
                            "/fields/type");
      }
      if (*ep_ref->z == z_obj) {
        throw ParaxialError(field_error(
            caller, "the entrance pupil at the reference wavelength lies in the object plane"));
      }
      const double h = (z_obj - *ep_ref->z) * value;
      return {z_ep, 0.0, h / (z_obj - z_ep)};
    }
    case model::FieldType::ObjectHeight: {
      if (infinite) {
        throw ParaxialError(field_error(caller, "an object height needs a finite object distance"),
                            "/fields/type");
      }
      const RayStart r = unit_ray(false, z_ep);
      return {r.z, value, value * r.u};
    }
    case model::FieldType::ParaxialImageHeight: {
      // Linear in the field value: the unit chief ray at the reference wavelength reaches the
      // reference paraxial image at h'_1, so the field is value / h'_1 as a slope or an object
      // height (as rtt-trace, #8, #31); the chief ray at `wavelength` scales its own unit ray.
      const FirstOrder& fo_ref = reference_first_order();
      if (!fo_ref.image_z) {
        throw ParaxialError(
            field_error(caller, "a paraxial image height needs a finite paraxial image"));
      }
      RayStart r_ref;
      if (at_ref) {
        r_ref = unit_ray(infinite, z_ep);
      } else if (fo_ref.entrance_pupil && fo_ref.entrance_pupil->z) {
        if (!infinite && *fo_ref.entrance_pupil->z == z_obj) {
          throw ParaxialError(field_error(
              caller, "the entrance pupil at the reference wavelength lies in the object plane"));
        }
        r_ref = unit_ray(infinite, *fo_ref.entrance_pupil->z);
      } else if (!infinite) {
        // Entrance pupil at infinity at the reference wavelength: the chief ray of unit object
        // height is parallel to the axis (as rtt-trace, #96).
        r_ref = {z_obj, 1.0, 0.0};
      } else {
        throw ParaxialError(
            field_error(caller, "the entrance pupil at the reference wavelength lies at infinity"));
      }
      const auto ray = trace_ray(system, path, ref, r_ref.z, r_ref.y, r_ref.u);
      const RayAtEvent& last = ray.back();
      const double unit = last.y + (*fo_ref.image_z - last.z) * last.u;
      if (unit == 0.0)
        throw ParaxialError(field_error(caller, "the chief ray does not reach the image"));
      const double scale = value / unit;
      const RayStart r = unit_ray(infinite, z_ep);
      return {r.z, scale * r.y, scale * r.u};
    }
  }
  return {};
}

}  // namespace rtt::paraxial::detail
