#include "rtt/trace/apply_event.hpp"

#include <cmath>
#include <type_traits>
#include <variant>

#include "rtt/geom/intersect.hpp"

namespace rtt::trace {

namespace {

/// True if the local point (x, y) lies inside the aperture; the boundary belongs to it.
bool inside(const model::Aperture& aperture, double x, double y) noexcept {
  return std::visit(
      [x, y](const auto& a) {
        using A = std::decay_t<decltype(a)>;
        if constexpr (std::is_same_v<A, model::CircularAperture>) {
          const double r2 = x * x + y * y;
          return r2 <= a.radius * a.radius && r2 >= a.inner_radius * a.inner_radius;
        } else if constexpr (std::is_same_v<A, model::RectangularAperture>) {
          return std::abs(x) <= a.half_width_x && std::abs(y) <= a.half_width_y;
        } else {
          const double u = x / a.semi_axis_x;
          const double v = y / a.semi_axis_y;
          return u * u + v * v <= 1.0;
        }
      },
      aperture);
}

/// Intersection in local coordinates: analytic for plane and conic, Newton for other shapes.
geom::Intersection<double> intersect(const compile::CompiledShape& shape,
                                     const math::Vec3& origin,
                                     const math::Vec3& direction) {
  return std::visit([&](const auto& s) { return geom::intersect(s, origin, direction); }, shape);
}

}  // namespace

std::optional<math::Vec3> refract(const math::Vec3& d,
                                  const math::Vec3& normal,
                                  double n1,
                                  double n2) noexcept {
  // Born & Wolf, Principles of Optics, 7th ed., Sec. 3.2.2: n2 t - n1 d is parallel to the
  // normal and n2 sin(theta_t) = n1 sin(theta_i). With n oriented against d (cos_i = -d . n
  // >= 0) and mu = n1 / n2 this gives t = mu d + (mu cos_i - cos_t) n,
  // cos_t = sqrt(1 - mu^2 (1 - cos_i^2)); a negative radicand is total internal reflection.
  math::Vec3 n = normal;
  double cos_i = -d.dot(n);
  if (cos_i < 0.0) {
    n = -n;
    cos_i = -cos_i;
  }
  const double mu = n1 / n2;
  const double cos_t2 = 1.0 - mu * mu * (1.0 - cos_i * cos_i);
  if (cos_t2 < 0.0) {
    return std::nullopt;
  }
  return math::Vec3(mu * d + (mu * cos_i - std::sqrt(cos_t2)) * n).normalized();
}

math::Vec3 reflect(const math::Vec3& d, const math::Vec3& normal) noexcept {
  // Born & Wolf, Principles of Optics, 7th ed., Sec. 3.2.2 (law of reflection in vector form).
  return d - 2.0 * d.dot(normal) * normal;
}

RayState apply_event(const RayState& ray,
                     const compile::CompiledSurface& surface,
                     std::uint32_t surface_index,
                     model::EventKind kind,
                     double n_before,
                     double n_after) {
  if (ray.status != RayStatus::Alive) {
    return ray;
  }
  // Local coordinates of the surface (docs/architecture.md, Engine 2).
  const math::Vec3 origin = surface.to_local.apply_point(ray.pos);
  const math::Vec3 direction = surface.to_local.apply_vector(ray.dir);
  const geom::Intersection<double> hit = intersect(surface.shape, origin, direction);
  RayState out = ray;
  switch (hit.status) {
    case geom::HitStatus::Hit:
      break;
    case geom::HitStatus::Missed:
      out.status = RayStatus::Missed;
      return out;
    case geom::HitStatus::NoConvergence:
      out.status = RayStatus::NoConvergence;
      return out;
  }

  // The ray reaches the surface: move it to the hit point (t is the geometric path in mm
  // because |direction| = 1) and add the optical path in the medium before the surface.
  out.pos = surface.to_global.apply_point(hit.point);
  out.opl += n_before * hit.t;
  out.last_surface = surface_index;

  if (surface.aperture && !inside(*surface.aperture, hit.point.x(), hit.point.y())) {
    out.status = RayStatus::Vignetted;
    return out;
  }
  if (std::holds_alternative<model::Absorber>(surface.interaction)) {
    out.status = RayStatus::Absorbed;
    return out;
  }

  math::Vec3 local_dir = direction;
  switch (kind) {
    case model::EventKind::Refract: {
      const auto t = refract(direction, hit.normal, n_before, n_after);
      if (!t) {
        out.status = RayStatus::Tir;
        return out;
      }
      local_dir = *t;
      break;
    }
    case model::EventKind::Reflect:
      local_dir = reflect(direction, hit.normal);
      break;
    case model::EventKind::Transmit:
      break;
    case model::EventKind::Diffract:
    case model::EventKind::Ordinary:
    case model::EventKind::Extraordinary:
      out.status = RayStatus::EventImpossible;  // M4
      return out;
  }
  out.dir = surface.to_global.apply_vector(local_dir);
  return out;
}

}  // namespace rtt::trace
