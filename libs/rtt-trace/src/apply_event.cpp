#include "rtt/trace/apply_event.hpp"

#include <cmath>
#include <numbers>
#include <type_traits>
#include <variant>

#include "rtt/geom/phase.hpp"
#include "rtt/math/units.hpp"
#include "rtt/polar/birefringence.hpp"
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
        } else if constexpr (std::is_same_v<I, model::Absorber> ||
                             std::is_same_v<I, model::IdealLens> ||
                             std::is_same_v<I, model::IdealCylinderLens>) {
          // Absorber: handled before the event. Ideal lenses: interim state of #178 until the
          // tracer applies them (ADR 0031, point 8): EventImpossible, never silently undeflected.
          return std::nullopt;
        } else {
          static_assert(kUnhandledInteraction<I>, "interaction not handled");
        }
      },
      surface.interaction);
}

/// 2 pi as in rtt::geom::LinearGratingPhase, so that g / (2 pi) is exact for 2^k lines/mm.
constexpr double kTwoPi = 2.0 * std::numbers::pi;

/// Power fraction of `order` at `surface` (ADR 0025, point 5): 1 without
/// diffraction_efficiency, otherwise the listed value or 0 for an order not listed.
double order_efficiency(const compile::CompiledSurface& surface, int order) noexcept {
  if (!surface.diffraction_efficiency) return 1.0;
  for (const model::DiffractionEfficiency& e : *surface.diffraction_efficiency) {
    if (e.order == order) return e.efficiency;
  }
  return 0.0;
}

/// Scales the weight of `out` after its P became p_event * P (ADR 0021): weight = s ||P_T||^2 / 2,
/// so the event multiplies ||P_T||^2 by new / old and s stays; the efficiency of a diffraction
/// order is a polarization-independent factor of s (ADR 0025, point 5). `out.dir` is the new
/// direction.
void apply_prt(RayState& out,
               const RayState& ray,
               const math::CMat3& p_event,
               const compile::CompiledSurface& surface,
               int order) noexcept {
  const double before = transverse_norm2(ray.prt, ray.dir);
  out.prt = p_event * ray.prt;
  if (before > 0.0) {
    out.weight *= transverse_norm2(out.prt, out.dir) / before;
  }
  if (surface.diffraction_efficiency) {
    out.weight *= order_efficiency(surface, order);
  }
}

/// Entry into, exit from or passage through a surface of a uniaxial crystal (ADR 0026, points 3
/// to 5), for events where `media` has crystal data on either side. `ray` is the incoming ray
/// (reading rule applied), `out` the ray moved to the hit point. The modes are solved in the
/// local coordinates of the surface, where incident_tangential and order_momentum live; P and
/// the eigenpolarizations are built from the global vectors.
RayState crystal_event(const RayState& ray,
                       RayState out,
                       const compile::CompiledSurface& surface,
                       const SurfaceHit& hit,
                       model::EventKind kind,
                       int order,
                       const EventMedia& media) noexcept {
  using model::EventKind;
  const auto stop = [&out](RayStatus status) {
    out.status = status;
    return out;
  };
  // Interactions at crystal surfaces in M4 (point 5): Fresnel and IdealAntiReflection, both as
  // the projection without reflection loss; compile rejects the others
  // (crystal.interaction_unsupported).
  if (!std::holds_alternative<model::Fresnel>(surface.interaction) &&
      !std::holds_alternative<model::IdealAntiReflection>(surface.interaction)) {
    return stop(RayStatus::EventImpossible);
  }
  if (kind == EventKind::Transmit) {
    // Inside the crystal with order 0 (point 4): mode, wave, dir and P stay (dummy passage); the
    // efficiency of order 0 applies as at any surface (ADR 0025, point 5).
    if (!media.crystal_before || !media.crystal_after || order != 0) {
      return stop(RayStatus::EventImpossible);
    }
    if (surface.diffraction_efficiency) out.weight *= order_efficiency(surface, 0);
    return out;
  }
  const bool entry = (kind == EventKind::Ordinary || kind == EventKind::Extraordinary) &&
                     !media.crystal_before && media.crystal_after;
  const bool exit = kind == EventKind::Refract && media.crystal_before && !media.crystal_after &&
                    media.crystal_mode != compile::CrystalMode::None && ray.mode_index > 0.0;
  if (!entry && !exit) return stop(RayStatus::EventImpossible);
  if (order != 0 && (!(media.wavelength_um > 0.0) || surface.phase_functions.empty())) {
    return stop(RayStatus::EventImpossible);
  }

  // Tangential phase matching (Lam, Eq. (2.11)) with the term of the order (ADR 0025, point 2);
  // the normal is oriented into the medium after the event, which the energy direction S enters.
  const math::Vec3 normal_after =
      hit.direction.dot(hit.normal) >= 0.0 ? hit.normal : math::Vec3(-hit.normal);
  const math::Vec3 tau_zero = incident_tangential(ray, surface, hit, media.before.real());
  const math::Vec3 tau =
      order == 0 ? tau_zero : tau_zero + order_momentum(surface, hit, order, media.wavelength_um);

  // No real solution: Tir for order 0, Evanescent for an order whose order 0 exists (ADR 0026,
  // point 4; ADR 0025, point 7).
  math::CMat3 p_event;
  if (entry) {
    const CrystalSide& c = *media.crystal_after;
    const polar::Mode mode =
        kind == EventKind::Ordinary ? polar::Mode::Ordinary : polar::Mode::Extraordinary;
    const math::Vec3 axis = surface.to_local.apply_vector(c.axis);
    const auto solve = [&](const math::Vec3& t) {
      return polar::uniaxial_mode<double>(mode, t, normal_after, c.n_ordinary, c.n_extraordinary,
                                          axis);
    };
    std::optional<polar::ModeSolution<double>> m = solve(tau_zero);
    if (!m) return stop(RayStatus::Tir);
    if (order != 0) {
      m = solve(tau);
      if (!m) return stop(RayStatus::Evanescent);
    }
    // Global vectors; E from the global axis, so that the degenerate case follows the global
    // axis rule of prt.hpp (point 5).
    polar::ModeSolution<double> g;
    g.k = surface.to_global.apply_vector(m->k);
    g.s = surface.to_global.apply_vector(m->s);
    g.n = m->n;
    g.e = polar::eigen_polarization<double>(mode, g.k, g.s, c.axis);
    p_event = polar::prt_crystal_entry<double>(ray.dir, g);
    out.dir = g.s;
    out.wave = g.k;
    out.mode_index = g.n;
  } else {
    const CrystalSide& c = *media.crystal_before;
    const double n_out = media.after.real();
    if (!polar::isotropic_from_tangential<double>(tau_zero, normal_after, n_out)) {
      return stop(RayStatus::Tir);
    }
    const std::optional<math::Vec3> t =
        polar::isotropic_from_tangential<double>(tau, normal_after, n_out);
    if (!t) return stop(RayStatus::Evanescent);
    const polar::Mode mode = media.crystal_mode == compile::CrystalMode::Ordinary
                                 ? polar::Mode::Ordinary
                                 : polar::Mode::Extraordinary;
    const math::Vec3 e_m = polar::eigen_polarization<double>(mode, ray.wave, ray.dir, c.axis);
    out.dir = surface.to_global.apply_vector(*t);
    p_event = polar::prt_crystal_exit<double>(ray.dir, e_m, out.dir);
    out.wave = out.dir;
    out.mode_index = 0.0;
  }
  if (!p_event.allFinite()) return stop(RayStatus::EventImpossible);
  out.opl += order_opl(surface, hit, order, media.wavelength_um);
  apply_prt(out, ray, p_event, surface, order);
  return out;
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
  // t is the geometric path in mm because |dir| = 1; OPL = Re(n) * path. In a crystal (mode
  // index > 0) the path runs along S = dir and the phase along k = wave: OPL = n l (k . S)
  // (Lam, Eq. (2.17), Fig. 2.13; ADR 0026, point 3). Isotropic media keep the first form.
  RayState out = ray;
  out.pos = surface.to_global.apply_point(hit.point);
  if (ray.mode_index > 0.0) {
    out.opl += ray.mode_index * hit.t * ray.wave.dot(ray.dir);
  } else {
    out.opl += media.before.real() * hit.t;
  }
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

math::Vec3 incident_tangential(const RayState& ray,
                               const compile::CompiledSurface& surface,
                               const SurfaceHit& hit,
                               double n_before) noexcept {
  const math::Vec3& n = hit.normal;
  if (ray.mode_index > 0.0) {
    // In a crystal the phase follows k: n k with n the mode index (ADR 0026, point 3).
    const math::Vec3 k = ray.mode_index * surface.to_local.apply_vector(ray.wave);
    return k - n * n.dot(k);
  }
  const math::Vec3& d = hit.direction;
  return n_before * (d - n * n.dot(d));
}

math::Vec3 order_momentum(const compile::CompiledSurface& surface,
                          const SurfaceHit& hit,
                          int order,
                          double wavelength_um) noexcept {
  if (order == 0) return math::Vec3::Zero();
  // ADR 0025, point 2: m lambda0 / (2 pi) g_par (Mansuripur, Eq. (7b), with F = phi / (2 pi)).
  // lambda0 by division: exact for wavelengths 1000 * 2^-k um, unlike um_to_mm (* 1e-3).
  const auto [gx, gy] =
      geom::phase_grad<double>(surface.phase_functions, hit.point.x(), hit.point.y());
  const math::Vec3 periods =
      geom::tangential_gradient<double>({gx / kTwoPi, gy / kTwoPi}, hit.normal);
  return (static_cast<double>(order) * (wavelength_um / 1000.0)) * periods;
}

std::optional<math::Vec3> order_direction(const math::Vec3& tau,
                                          const math::Vec3& unit_normal,
                                          double n_out,
                                          double side) noexcept {
  // |t'| = 1 with the tangential part tau / n_out (ADR 0025, point 2); equality is evanescent
  // (point 7): the order would run along the surface.
  const double tau2 = tau.squaredNorm();
  const double n2 = n_out * n_out;
  if (!(tau2 < n2)) return std::nullopt;
  return (tau + (side * std::sqrt(n2 - tau2)) * unit_normal) / n_out;
}

double order_opl(const compile::CompiledSurface& surface,
                 const SurfaceHit& hit,
                 int order,
                 double wavelength_um) noexcept {
  if (order == 0) return 0.0;
  // ADR 0025, point 3: the order carries the phase m phi (Mansuripur, Eq. (5)), i.e. the path
  // m phi lambda0 / (2 pi) (Mansuripur, Sec. 2, before Eq. (3a): Phi = 2 pi OPD / lambda0).
  const double phi = geom::phase<double>(surface.phase_functions, hit.point.x(), hit.point.y());
  return (static_cast<double>(order) * (wavelength_um / 1000.0)) * (phi / kTwoPi);
}

math::Mat3 rotation_between(const math::Vec3& a, const math::Vec3& b) noexcept {
  // Diebel (2006), Eqs. (183)-(187), transposed (active rotation, see the header): with
  // c = cos(alpha) = a . b, v = a x b = sin(alpha) n, R = c I + [v]_x + (1 - c) n n^T. The unit
  // axis n = v / |v| keeps the last term accurate also close to a . b = -1, where the form
  // v v^T / (1 + c) would divide two small numbers.
  if (a == b) return math::Mat3::Identity();
  const math::Vec3 v = a.cross(b);
  const double s = v.norm();
  const double c = a.dot(b);
  if (!(s > 0.0)) return math::Mat3::Identity();  // parallel within rounding (c > -1 by @pre)
  // Element by element (Diebel, Eqs. (185)-(187) with the half-angle products written as
  // cos(alpha) and 1 - cos(alpha), transposed): R_ij = c delta_ij + k n_i n_j + eps_ikj v_k,
  // k = 1 - c. Scalar code instead of an Eigen expression with an outer product, which GCC
  // flags with -Wnull-dereference at -O2 (a false positive inside Eigen).
  const math::Vec3 n = v / s;
  const double k = 1.0 - c;
  math::Mat3 r;
  r(0, 0) = c + k * n.x() * n.x();
  r(0, 1) = k * n.x() * n.y() - v.z();
  r(0, 2) = k * n.x() * n.z() + v.y();
  r(1, 0) = k * n.y() * n.x() + v.z();
  r(1, 1) = c + k * n.y() * n.y();
  r(1, 2) = k * n.y() * n.z() - v.x();
  r(2, 0) = k * n.z() * n.x() - v.y();
  r(2, 1) = k * n.z() * n.y() + v.x();
  r(2, 2) = c + k * n.z() * n.z();
  return r;
}

RayState apply_event(const RayState& incoming,
                     const compile::CompiledSurface& surface,
                     const SurfaceHit& hit,
                     std::uint32_t surface_index,
                     model::EventKind kind,
                     int order,
                     const EventMedia& media) noexcept {
  const double n_before = media.before.real();
  const double n_after = media.after.real();
  if (incoming.status != RayStatus::Alive) {
    return incoming;
  }
  // Reading rule (ADR 0026, point 3): wave and mode_index count only for mode_index > 0.
  RayState ray = incoming;
  if (!(ray.mode_index > 0.0)) {
    ray.wave = ray.dir;
    ray.mode_index = 0.0;
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
  if (media.crystal_before || media.crystal_after) {
    return crystal_event(ray, out, surface, hit, kind, order, media);
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
      out.status = RayStatus::EventImpossible;  // no crystal after the event
      return out;
  }

  // Diffraction order m != 0 (ADR 0025): order 0 exists (Tir was decided above, point 7); the
  // order needs the wavelength and phase layers (validate allows orders only there).
  math::Vec3 local_order = local_dir;
  if (order != 0) {
    if (!(media.wavelength_um > 0.0) || surface.phase_functions.empty()) {
      out.status = RayStatus::EventImpossible;
      return out;
    }
    const double n_out = kind == model::EventKind::Reflect ? n_before : n_after;
    const math::Vec3 tau = incident_tangential(ray, surface, hit, n_before) +
                           order_momentum(surface, hit, order, media.wavelength_um);
    const double side = local_dir.dot(hit.normal) >= 0.0 ? 1.0 : -1.0;
    const std::optional<math::Vec3> t = order_direction(tau, hit.normal, n_out, side);
    if (!t) {
      out.status = RayStatus::Evanescent;
      return out;
    }
    local_order = *t;
    out.opl += order_opl(surface, hit, order, media.wavelength_um);
  }

  const math::Vec3 k_zero = surface.to_global.apply_vector(local_dir);
  const math::Vec3 normal = surface.to_global.apply_vector(hit.normal);
  const std::optional<math::CMat3> p =
      interaction_prt(surface, kind, ray.dir, k_zero, normal, media);
  // Non-finite amplitudes (grazing incidence with q_i = 0, a coating layer with q = 0 exactly:
  // preconditions of rtt-polar and rtt-coating) stop the ray instead of spreading NaN.
  if (!p || !p->allFinite()) {
    out.status = RayStatus::EventImpossible;
    return out;
  }
  math::CMat3 p_event = *p;
  math::Vec3 k_out = k_zero;
  if (order != 0) {
    // ADR 0025, point 6: interface as for order 0, then the thin layer turns k_0 into k_m
    // without rotating the polarization about the direction.
    k_out = surface.to_global.apply_vector(local_order);
    p_event = rotation_between(k_zero, k_out).cast<math::Complex>() * p_event;
  }
  out.dir = k_out;
  apply_prt(out, ray, p_event, surface, order);
  // Isotropic on both sides: no mode, wave = dir (ADR 0026, point 3).
  out.wave = out.dir;
  out.mode_index = 0.0;
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
  return apply_event(ray, surface, hit, surface_index, kind, 0, media);
}

RayState sequential_step(const RayState& ray,
                         const compile::CompiledSurface& surface,
                         std::uint32_t surface_index,
                         model::EventKind kind,
                         int order,
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
  return apply_event(ray, surface, hit, surface_index, kind, order, media);
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
  return sequential_step(ray, surface, surface_index, kind, 0, media);
}

}  // namespace rtt::trace
