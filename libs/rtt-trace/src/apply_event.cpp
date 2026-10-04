#include "rtt/trace/apply_event.hpp"

#include <cmath>
#include <type_traits>
#include <variant>

namespace rtt::trace {

namespace {

template <class>
inline constexpr bool kUnhandledAperture = false;

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
        } else if constexpr (std::is_same_v<A, model::EllipticalAperture>) {
          const double u = x / a.semi_axis_x;
          const double v = y / a.semi_axis_y;
          return u * u + v * v <= 1.0;
        } else {
          static_assert(kUnhandledAperture<A>, "aperture type not handled");
        }
      },
      aperture);
}

}  // namespace

std::optional<math::Vec3> refract(const math::Vec3& d,
                                  const math::Vec3& normal,
                                  double n1,
                                  double n2) noexcept {
  // Derived from Snell's law in vector form (cf. Born & Wolf, Principles of Optics, 7th ed.,
  // Sec. 3.2.2, equation number: see #19): n2 t - n1 d is parallel to the normal and
  // n2 sin(theta_t) = n1 sin(theta_i). With n oriented against d (cos_i = -d . n
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
  // Law of reflection in vector form, cf. Born & Wolf, 7th ed., Sec. 3.2.2 (see #19).
  return d - 2.0 * d.dot(normal) * normal;
}

SurfaceHit intersect_surface(const RayState& ray,
                             const compile::CompiledSurface& surface) noexcept {
  // Local coordinates of the surface (docs/architecture.md, Engine 2).
  const math::Vec3 origin = surface.to_local.apply_point(ray.pos);
  const math::Vec3 direction = surface.to_local.apply_vector(ray.dir);
  const geom::Intersection<double> hit = std::visit(
      [&](const auto& shape) { return geom::intersect(shape, origin, direction); }, surface.shape);
  return {hit.status, hit.t, hit.point, hit.normal, direction};
}

bool inside_aperture(const compile::CompiledSurface& surface, const SurfaceHit& hit) noexcept {
  return !surface.aperture || inside(*surface.aperture, hit.point.x(), hit.point.y());
}

RayState move_to_hit(const RayState& ray,
                     const compile::CompiledSurface& surface,
                     const SurfaceHit& hit,
                     std::uint32_t surface_index,
                     double n_before) noexcept {
  // t is the geometric path in mm because |dir| = 1; OPL = n * path (M1: real part of n).
  RayState out = ray;
  out.pos = surface.to_global.apply_point(hit.point);
  out.opl += n_before * hit.t;
  out.last_surface = surface_index;
  return out;
}

RayState apply_event(const RayState& ray,
                     const compile::CompiledSurface& surface,
                     const SurfaceHit& hit,
                     std::uint32_t surface_index,
                     model::EventKind kind,
                     double n_before,
                     double n_after) noexcept {
  if (ray.status != RayStatus::Alive) {
    return ray;
  }
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
  out = move_to_hit(ray, surface, hit, surface_index, n_before);
  if (std::holds_alternative<model::Absorber>(surface.interaction)) {
    out.status = RayStatus::Absorbed;
    return out;
  }

  math::Vec3 local_dir = hit.direction;
  switch (kind) {
    case model::EventKind::Refract: {
      const auto t = refract(hit.direction, hit.normal, n_before, n_after);
      if (!t) {
        out.status = RayStatus::Tir;
        return out;
      }
      local_dir = *t;
      break;
    }
    case model::EventKind::Reflect:
      local_dir = reflect(hit.direction, hit.normal);
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

RayState sequential_step(const RayState& ray,
                         const compile::CompiledSurface& surface,
                         std::uint32_t surface_index,
                         model::EventKind kind,
                         double n_before,
                         double n_after) noexcept {
  if (ray.status != RayStatus::Alive) {
    return ray;
  }
  const SurfaceHit hit = intersect_surface(ray, surface);
  if (hit.status == geom::HitStatus::Hit && !inside_aperture(surface, hit)) {
    RayState out = move_to_hit(ray, surface, hit, surface_index, n_before);
    out.status = RayStatus::Vignetted;
    return out;
  }
  return apply_event(ray, surface, hit, surface_index, kind, n_before, n_after);
}

}  // namespace rtt::trace
