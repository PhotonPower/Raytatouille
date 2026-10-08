#include "rtt/trace/apply_event.hpp"

#include <cmath>
#include <numbers>
#include <type_traits>
#include <variant>

#include "rtt/math/units.hpp"
#include "rtt/polar/ideal.hpp"
#include "rtt/polar/interface.hpp"

namespace rtt::trace {

namespace {

template <class>
inline constexpr bool kUnhandledAperture = false;

/// True if the local point (x, y) lies inside the aperture; the boundary belongs to it.
bool inside(const model::Aperture& aperture, double x, double y) noexcept {
  // Rims are inclusive with kApertureTolerance (decided for #50).
  constexpr double tol = kApertureTolerance;
  return std::visit(
      [x, y](const auto& a) {
        using A = std::decay_t<decltype(a)>;
        if constexpr (std::is_same_v<A, model::CircularAperture>) {
          const double r2 = x * x + y * y;
          const double outer = a.radius + tol;
          const double inner = a.inner_radius - tol;
          return r2 <= outer * outer && (inner <= 0.0 || r2 >= inner * inner);
        } else if constexpr (std::is_same_v<A, model::RectangularAperture>) {
          return std::abs(x) <= a.half_width_x + tol && std::abs(y) <= a.half_width_y + tol;
        } else if constexpr (std::is_same_v<A, model::EllipticalAperture>) {
          const double u = x / (a.semi_axis_x + tol);
          const double v = y / (a.semi_axis_y + tol);
          return u * u + v * v <= 1.0;
        } else {
          static_assert(kUnhandledAperture<A>, "aperture type not handled");
        }
      },
      aperture);
}

template <class>
inline constexpr bool kUnhandledInteraction = false;

/// ||P_T||_F^2 of an accumulated PRT matrix P with current direction k (unit). Every factor maps
/// k onto k' (Lam, Eq. (3.2)) and, as O_out J_3D O_in^-1 with orthonormal O and
/// J_3D = diag(a_s, a_p, 1) (Lam, Eq. (3.4); docs/quellen.md), the transverse plane of k onto that
/// of k'; so does the product. Hence P = P_T + k k_0^T with P_T k_0 = 0 and k^T P_T = 0, which
/// gives k_0 = P^T k and P_T = (I - k k^T) P. Computed this way, the error is relative to P_T
/// (not absolute as for ||P||^2 - 1), so small weights (ghosts, crossed polarizers) stay accurate.
double transverse_norm2(const math::CMat3& p, const math::Vec3& k) noexcept {
  const math::Mat3 projector = math::Mat3::Identity() - k * k.transpose();
  return (projector.cast<math::Complex>() * p).squaredNorm();
}

/// PRT matrix with power-normalised amplitudes a_s, a_p given directly (ideal elements).
math::CMat3 ideal_prt(const math::Vec3& k_in,
                      const math::Vec3& k_out,
                      const math::Vec3& normal,
                      math::Complex a_s,
                      math::Complex a_p) noexcept {
  return polar::interface_prt(k_in, k_out, normal, a_s, a_p, polar::PowerFactors<double>{});
}

/// Transverse part of an ideal element's axis; none if it is (nearly) parallel to k.
std::optional<math::Vec3> transverse(const std::optional<math::Vec3>& axis,
                                     const math::Vec3& k) noexcept {
  if (!axis || !(axis->norm() > 0.0)) return std::nullopt;
  const math::Vec3 a = *axis - axis->dot(k) * k;
  if (!(a.norm() >= polar::kAxisAlongK * axis->norm())) return std::nullopt;
  return a;
}

/// Power-normalised PRT matrix of the surface interaction for `kind` (ADR 0021); identity if
/// the event does not change the polarization state, none if the combination is impossible.
/// Directions and normal are global unit vectors; k_out is the direction after the event.
std::optional<math::CMat3> interaction_prt(const compile::CompiledSurface& surface,
                                           model::EventKind kind,
                                           const math::Vec3& k_in,
                                           const math::Vec3& k_out,
                                           const math::Vec3& normal,
                                           const EventMedia& media) noexcept {
  using model::EventKind;
  const bool reflect = kind == EventKind::Reflect;
  const bool refract = kind == EventKind::Refract;
  const bool transmit = kind == EventKind::Transmit;
  return std::visit(
      [&](const auto& interaction) -> std::optional<math::CMat3> {
        using I = std::decay_t<decltype(interaction)>;
        if constexpr (std::is_same_v<I, model::Fresnel>) {
          if (transmit) return math::CMat3::Identity();
          if (refract)
            return polar::fresnel_prt(k_in, k_out, normal, media.before, media.after, false);
          // A mirror without material has no medium beyond the surface: ideal conductor
          // (kappa -> infinity), r_s = -1, r_p = +1 (docs/architecture.md, Polarisation).
          if (surface.element_kind == model::ElementKind::Mirror && media.beyond == media.before) {
            return ideal_prt(k_in, k_out, normal, -1.0, 1.0);
          }
          return polar::fresnel_prt(k_in, k_out, normal, media.before, media.beyond, true);
        } else if constexpr (std::is_same_v<I, model::IdealMirror>) {
          if (!reflect) return std::nullopt;
          return ideal_prt(k_in, k_out, normal, -1.0, 1.0);
        } else if constexpr (std::is_same_v<I, model::IdealAntiReflection>) {
          if (transmit) return math::CMat3::Identity();
          // Reflect: a vanishing, not an impossible reflection (R = 0, decided for #61).
          if (reflect) return ideal_prt(k_in, k_out, normal, 0.0, 0.0);
          return ideal_prt(k_in, k_out, normal, 1.0, 1.0);
        } else if constexpr (std::is_same_v<I, model::IdealBeamSplitter>) {
          if (reflect) {
            return ideal_prt(k_in, k_out, normal, -std::sqrt(interaction.reflectance_s),
                             std::sqrt(interaction.reflectance_p));
          }
          return ideal_prt(k_in, k_out, normal, std::sqrt(1.0 - interaction.reflectance_s),
                           std::sqrt(1.0 - interaction.reflectance_p));
        } else if constexpr (std::is_same_v<I, model::CoatingRef>) {
          if (transmit) return math::CMat3::Identity();
          if (!surface.coating) return std::nullopt;
          // Byrnes, Eqs. (6)-(15) via rtt-coating; the layers come in the ray's order (ADR 0019).
          const double xi = polar::tangential_invariant(k_in, normal, media.before);
          const math::Complex n_out = reflect ? media.beyond : media.after;
          const auto a = coating::stack_amplitudes<double>(media.before, media.layers, n_out, xi,
                                                           media.wavelength_um);
          if (reflect) return ideal_prt(k_in, k_out, normal, a.rs, a.rp);
          return polar::interface_prt(k_in, k_out, normal, a.ts, a.tp,
                                      polar::transmission_power_factors(media.before, n_out, xi));
        } else if constexpr (std::is_same_v<I, model::IdealPolarizer>) {
          const auto axis = transverse(surface.ideal_axis, k_in);
          if (!transmit || !axis) return std::nullopt;
          return polar::ideal_polarizer(*axis, interaction.extinction_ratio, k_in);
        } else if constexpr (std::is_same_v<I, model::IdealRetarder>) {
          const auto axis = transverse(surface.ideal_axis, k_in);
          if (!transmit || !axis) return std::nullopt;
          return polar::linear_retarder(*axis, interaction.retardance_waves, k_in);
        } else if constexpr (std::is_same_v<I, model::Absorber>) {
          return std::nullopt;  // handled before the event
        } else {
          static_assert(kUnhandledInteraction<I>, "interaction not handled");
        }
      },
      surface.interaction);
}

}  // namespace

std::optional<math::Vec3> refract(const math::Vec3& d,
                                  const math::Vec3& normal,
                                  double n1,
                                  double n2) noexcept {
  // B. de Greve, Reflections and Refractions in Ray Tracing (2006), Sec. 6, Eqs. (22), (23),
  // (28): t = mu d + (mu cos_i - cos_t) n, mu = n1 / n2, sin_t^2 = mu^2 (1 - cos_i^2),
  // cos_t = sqrt(1 - sin_t^2); total internal reflection for sin_t^2 > 1 (Eq. (24)).
  // de Greve's n points into the incident medium, so cos_i = -d . n >= 0 (Sec. 6). Note: the
  // text after Eq. (20a) writes cos_i = i . n, which contradicts that orientation; Sec. 6 is
  // correct. rtt-geom returns n with +z at the vertex, so n is flipped when d . n > 0.
  math::Vec3 n = normal;
  double cos_i = -d.dot(n);
  if (cos_i < 0.0) {
    n = -n;
    cos_i = -cos_i;
  }
  const double mu = n1 / n2;
  const double sin_t2 = mu * mu * (1.0 - cos_i * cos_i);
  if (sin_t2 > 1.0) {
    return std::nullopt;
  }
  const double cos_t2 = 1.0 - sin_t2;
  return math::Vec3(mu * d + (mu * cos_i - std::sqrt(cos_t2)) * n).normalized();
}

math::Vec3 reflect(const math::Vec3& d, const math::Vec3& normal) noexcept {
  // de Greve (2006), Eq. (13); valid for either orientation of n.
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
                     const EventMedia& media) noexcept {
  // t is the geometric path in mm because |dir| = 1; OPL = Re(n) * path.
  RayState out = ray;
  out.pos = surface.to_global.apply_point(hit.point);
  out.opl += media.before.real() * hit.t;
  out.last_surface = surface_index;
  // Byrnes, arXiv:1603.02720v5, Eqs. (1), (2): E ~ exp(i 2 pi n z / lambda_vac), so |E|^2, the
  // power by Eqs. (17), (18), decays as exp(-4 pi kappa z / lambda_vac); lambda in mm like t.
  const double kappa = media.before.imag();
  if (kappa != 0.0) {
    out.weight *=
        std::exp(-4.0 * std::numbers::pi * kappa * hit.t / math::um_to_mm(media.wavelength_um));
  }
  return out;
}

RayState move_to_hit(const RayState& ray,
                     const compile::CompiledSurface& surface,
                     const SurfaceHit& hit,
                     std::uint32_t surface_index,
                     double n_before) noexcept {
  EventMedia media;
  media.before = n_before;
  media.after = n_before;
  media.beyond = n_before;
  return move_to_hit(ray, surface, hit, surface_index, media);
}

RayState apply_event(const RayState& ray,
                     const compile::CompiledSurface& surface,
                     const SurfaceHit& hit,
                     std::uint32_t surface_index,
                     model::EventKind kind,
                     const EventMedia& media) noexcept {
  const double n_before = media.before.real();
  const double n_after = media.after.real();
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
  out = move_to_hit(ray, surface, hit, surface_index, media);
  if (std::holds_alternative<model::Absorber>(surface.interaction)) {
    out.status = RayStatus::Absorbed;
    out.weight = 0.0;
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
    case model::EventKind::Ordinary:
    case model::EventKind::Extraordinary:
      out.status = RayStatus::EventImpossible;  // M4
      return out;
  }
  const math::Vec3 k_out = surface.to_global.apply_vector(local_dir);
  const math::Vec3 normal = surface.to_global.apply_vector(hit.normal);
  const std::optional<math::CMat3> p =
      interaction_prt(surface, kind, ray.dir, k_out, normal, media);
  // Non-finite amplitudes (grazing incidence with q_i = 0, a coating layer with q = 0 exactly:
  // preconditions of rtt-polar and rtt-coating) stop the ray instead of spreading NaN.
  if (!p || !p->allFinite()) {
    out.status = RayStatus::EventImpossible;
    return out;
  }
  // weight = s ||P_T||^2 / 2 (ADR 0021): the event multiplies ||P_T||^2 by new / old, s stays.
  const double before = transverse_norm2(ray.prt, ray.dir);
  out.prt = *p * ray.prt;
  if (before > 0.0) {
    out.weight *= transverse_norm2(out.prt, k_out) / before;
  }
  out.dir = k_out;
  return out;
}

RayState apply_event(const RayState& ray,
                     const compile::CompiledSurface& surface,
                     const SurfaceHit& hit,
                     std::uint32_t surface_index,
                     model::EventKind kind,
                     double n_before,
                     double n_after) noexcept {
  EventMedia media;
  media.before = n_before;
  media.after = n_after;
  media.beyond = n_after;
  return apply_event(ray, surface, hit, surface_index, kind, media);
}

RayState sequential_step(const RayState& ray,
                         const compile::CompiledSurface& surface,
                         std::uint32_t surface_index,
                         model::EventKind kind,
                         const EventMedia& media) noexcept {
  if (ray.status != RayStatus::Alive) {
    return ray;
  }
  const SurfaceHit hit = intersect_surface(ray, surface);
  if (hit.status == geom::HitStatus::Hit && !inside_aperture(surface, hit)) {
    RayState out = move_to_hit(ray, surface, hit, surface_index, media);
    out.status = RayStatus::Vignetted;
    return out;
  }
  return apply_event(ray, surface, hit, surface_index, kind, media);
}

RayState sequential_step(const RayState& ray,
                         const compile::CompiledSurface& surface,
                         std::uint32_t surface_index,
                         model::EventKind kind,
                         double n_before,
                         double n_after) noexcept {
  EventMedia media;
  media.before = n_before;
  media.after = n_after;
  media.beyond = n_after;
  return sequential_step(ray, surface, surface_index, kind, media);
}

}  // namespace rtt::trace
