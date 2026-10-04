#pragma once

/// @file apply_event.hpp
/// Building blocks of one surface interaction, independent of paths (docs/architecture.md,
/// Engine 2): intersect_surface() finds the hit, inside_aperture() checks the aperture and
/// apply_event(ray, hit, kind) executes the event. The sequential tracer combines them per path
/// event; the non-sequential tracer (M9) reuses intersect_surface() and apply_event() unchanged
/// but treats a hit outside the aperture as "surface not there" instead of Vignetted.
///
/// Conventions (docs/architecture.md, Konventionen): positions in mm, unit directions, both in
/// global coordinates (right-handed, optical axis +z) unless stated as local; optical path
/// length in mm.

#include <cstdint>
#include <optional>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/geom/intersect.hpp"
#include "rtt/math/types.hpp"
#include "rtt/model/path.hpp"
#include "rtt/trace/ray_batch.hpp"

namespace rtt::trace {

/// State of one ray in global coordinates (one row of a RayBatch, without the fields that
/// M1 does not change).
struct RayState {
  math::Vec3 pos = math::Vec3::Zero();   ///< position, mm, global
  math::Vec3 dir = math::Vec3::UnitZ();  ///< unit direction, global (|dir| = 1 is required)
  double opl = 0.0;                      ///< accumulated optical path length, mm
  RayStatus status = RayStatus::Alive;
  std::uint32_t last_surface = kNoSurface;  ///< index of the last surface hit
};

/// Hit of a ray on a surface, in the local coordinates of the surface.
struct SurfaceHit {
  geom::HitStatus status = geom::HitStatus::Missed;
  double t = 0.0;  ///< geometric path from the ray position to the hit point, mm
  math::Vec3 point = math::Vec3::Zero();      ///< hit point, local, mm
  math::Vec3 normal = math::Vec3::Zero();     ///< unit normal, local, +z at the vertex
  math::Vec3 direction = math::Vec3::Zero();  ///< incident unit direction, local
};

/// Refracted direction (B. de Greve, Reflections and Refractions in Ray Tracing, 2006, Sec. 6,
/// Eqs. (22), (23), (28); see also M. Born, E. Wolf, Principles of Optics, 7th ed., Sec. 3.2.2):
/// t = mu d + (mu cos_i - cos_t) n with mu = n1 / n2, cos_i = -d . n,
/// sin_t^2 = mu^2 (1 - cos_i^2), cos_t = sqrt(1 - sin_t^2). As in de Greve, n is oriented into
/// the incident medium (n . d < 0); the given normal is flipped if necessary.
/// Total internal reflection for sin_t^2 > 1 (de Greve, Eq. (24)).
/// @param d      unit incident direction
/// @param normal unit surface normal, either orientation
/// @param n1     refractive index before the surface (real part)
/// @param n2     refractive index after the surface (real part)
/// @return unit refracted direction, or std::nullopt for total internal reflection.
[[nodiscard]] std::optional<math::Vec3> refract(const math::Vec3& d,
                                                const math::Vec3& normal,
                                                double n1,
                                                double n2) noexcept;

/// Reflected direction r = d - 2 (d . n) n (B. de Greve, Reflections and Refractions in Ray
/// Tracing, 2006, Eq. (13); independent of the orientation of n; see also Born & Wolf,
/// Principles of Optics, 7th ed., Sec. 3.2.2).
/// @param d      unit incident direction
/// @param normal unit surface normal, either orientation
[[nodiscard]] math::Vec3 reflect(const math::Vec3& d, const math::Vec3& normal) noexcept;

/// Intersects the ray with the surface in its local coordinates: analytic for plane and conic,
/// Newton for the even asphere (rtt::geom::intersect). Status Missed or NoConvergence if there
/// is no valid hit; never throws (ADR 0009).
/// @param ray     ray in global coordinates, |dir| = 1
/// @param surface compiled surface (pose and shape are used)
[[nodiscard]] SurfaceHit intersect_surface(const RayState& ray,
                                           const compile::CompiledSurface& surface) noexcept;

/// Tolerance of the aperture check, mm: a hit counts as inside if it lies at most this far
/// outside a rim (outer and inner rims). Equal to the convergence limit of ray aiming
/// (rtt::trace::kAimTolerance), so rays aimed exactly at a stop rim pass (decided for #50).
inline constexpr double kApertureTolerance = 1e-9;

/// True if the local hit point lies inside the surface aperture or the surface has no aperture;
/// rims are inclusive with kApertureTolerance. Circular (r <= R + tol, with an inner radius
/// r >= R_i - tol), rectangular (|x| <= a + tol, |y| <= b + tol) and elliptical apertures
/// (semi-axes enlarged by tol) in local x, y, mm.
[[nodiscard]] bool inside_aperture(const compile::CompiledSurface& surface,
                                   const SurfaceHit& hit) noexcept;

/// Moves the ray to the hit point: position (global), OPL += n_before * t, last_surface.
/// Direction and status are unchanged.
/// @pre hit.status == Hit
[[nodiscard]] RayState move_to_hit(const RayState& ray,
                                   const compile::CompiledSurface& surface,
                                   const SurfaceHit& hit,
                                   std::uint32_t surface_index,
                                   double n_before) noexcept;

/// Executes `kind` at the hit (docs/architecture.md: apply_event(ray, hit, kind)).
///
/// A ray that is not Alive is returned unchanged; a hit with status Missed or NoConvergence
/// gives that status and leaves the ray unchanged. Otherwise the ray moves to the hit point
/// (move_to_hit) and then: Absorber interaction -> Absorbed; Diffract, Ordinary, Extraordinary
/// -> EventImpossible (M4); Refract -> refracted direction or Tir beyond the critical angle;
/// Reflect -> reflected direction; Transmit -> unchanged direction. A stopped ray stays at the
/// hit point with last_surface = surface_index. The aperture is not checked here.
///
/// M1 limits: refraction and OPL use Re(n); absorption (kappa) and Fresnel weights follow in M3.
/// Interactions other than Absorber (Fresnel, coatings, polarizers, retarders, ideal mirror,
/// AR, beam splitter) do not change the ray before M3; the event kind alone decides between
/// refraction and reflection. Phase layers are ignored before M4.
///
/// Never throws; physical problems are status flags (ADR 0009).
/// @param ray           incoming ray in global coordinates, |dir| = 1
/// @param surface       compiled surface (pose and interaction are used)
/// @param hit           result of intersect_surface(ray, surface)
/// @param surface_index index of `surface` in CompiledSystem::surfaces(), stored in last_surface
/// @param kind          event at this surface
/// @param n_before      refractive index (real part) of the medium before the surface
/// @param n_after       refractive index (real part) of the medium after the surface
[[nodiscard]] RayState apply_event(const RayState& ray,
                                   const compile::CompiledSurface& surface,
                                   const SurfaceHit& hit,
                                   std::uint32_t surface_index,
                                   model::EventKind kind,
                                   double n_before,
                                   double n_after) noexcept;

/// One sequential step: intersect_surface(), then Vignetted (ray at the hit point) if the hit
/// lies outside the aperture, otherwise apply_event(). Parameters as for apply_event().
[[nodiscard]] RayState sequential_step(const RayState& ray,
                                       const compile::CompiledSurface& surface,
                                       std::uint32_t surface_index,
                                       model::EventKind kind,
                                       double n_before,
                                       double n_after) noexcept;

}  // namespace rtt::trace
