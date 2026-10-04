#pragma once

/// @file apply_event.hpp
/// Execution of one surface event for one ray, independent of paths (docs/architecture.md,
/// Engine 2). The sequential tracer calls it per path event; the non-sequential tracer (M9)
/// reuses it unchanged.
///
/// Conventions (docs/architecture.md, Konventionen): positions in mm, unit directions, both in
/// global coordinates (right-handed, optical axis +z); optical path length in mm.

#include <cstdint>
#include <optional>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/math/types.hpp"
#include "rtt/model/path.hpp"
#include "rtt/trace/ray_batch.hpp"

namespace rtt::trace {

/// State of one ray in global coordinates (one row of a RayBatch, without the fields that
/// M1 does not change).
struct RayState {
  math::Vec3 pos = math::Vec3::Zero();   ///< position, mm
  math::Vec3 dir = math::Vec3::UnitZ();  ///< unit direction
  double opl = 0.0;                      ///< accumulated optical path length, mm
  RayStatus status = RayStatus::Alive;
  std::uint32_t last_surface = kNoSurface;  ///< index of the last surface hit
};

/// Refracted direction by the law of refraction in vector form
/// (M. Born, E. Wolf, Principles of Optics, 7th ed., Sec. 3.2.2):
/// t = mu d + (mu cos_i - cos_t) n with mu = n1 / n2, cos_i = -d . n, n oriented against d.
/// @param d      unit incident direction
/// @param normal unit surface normal, either orientation
/// @param n1     refractive index before the surface (real part)
/// @param n2     refractive index after the surface (real part)
/// @return unit refracted direction, or std::nullopt for total internal reflection.
[[nodiscard]] std::optional<math::Vec3> refract(const math::Vec3& d,
                                                const math::Vec3& normal,
                                                double n1,
                                                double n2) noexcept;

/// Reflected direction r = d - 2 (d . n) n (Born & Wolf, Principles of Optics, Sec. 3.2.2).
/// @param d      unit incident direction
/// @param normal unit surface normal, either orientation
[[nodiscard]] math::Vec3 reflect(const math::Vec3& d, const math::Vec3& normal) noexcept;

/// Applies one event at `surface` to `ray`.
///
/// Steps (docs/architecture.md, Engine 2): transform into the local coordinates of the surface,
/// intersect (analytic for plane and conic, Newton for the even asphere), check the aperture in
/// local x, y, apply the event, transform back. The OPL grows by n_before times the geometric
/// path to the hit point.
///
/// Outcome: a ray that is not Alive is returned unchanged. Missed and NoConvergence keep the
/// ray as it was (position, direction, OPL, last_surface). Otherwise the ray moves to the hit
/// point and last_surface becomes `surface_index`; then Vignetted (outside the aperture),
/// Absorbed (Absorber interaction), EventImpossible (Diffract, Ordinary, Extraordinary before
/// M4) and Tir (Refract beyond the critical angle) stop it there. Refract, Reflect and Transmit
/// set the new direction.
///
/// M1 limits: refraction and OPL use Re(n); absorption (kappa) and Fresnel weights follow in M3.
/// Interactions other than Absorber (Fresnel, coatings, polarizers, retarders, ideal mirror,
/// AR, beam splitter) do not change the ray before M3; the event kind alone decides between
/// refraction and reflection. Phase layers are ignored before M4.
///
/// Never throws for physical problems; they are status flags (ADR 0009).
/// @param ray           incoming ray in global coordinates
/// @param surface       compiled surface (pose, shape, aperture, interaction)
/// @param surface_index index of `surface` in CompiledSystem::surfaces(), stored in last_surface
/// @param kind          event at this surface
/// @param n_before      refractive index (real part) of the medium before the surface
/// @param n_after       refractive index (real part) of the medium after the surface
[[nodiscard]] RayState apply_event(const RayState& ray,
                                   const compile::CompiledSurface& surface,
                                   std::uint32_t surface_index,
                                   model::EventKind kind,
                                   double n_before,
                                   double n_after);

}  // namespace rtt::trace
