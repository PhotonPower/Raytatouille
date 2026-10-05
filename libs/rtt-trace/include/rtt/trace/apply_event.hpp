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
#include <span>

#include "rtt/coating/transfer_matrix.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/geom/intersect.hpp"
#include "rtt/math/types.hpp"
#include "rtt/model/path.hpp"
#include "rtt/trace/ray_batch.hpp"

namespace rtt::trace {

/// State of one ray in global coordinates (one row of a RayBatch without wavelength, field and
/// pupil columns).
struct RayState {
  math::Vec3 pos = math::Vec3::Zero();   ///< position, mm, global
  math::Vec3 dir = math::Vec3::UnitZ();  ///< unit direction, global (|dir| = 1 is required)
  double opl = 0.0;                      ///< accumulated optical path length, mm
  RayStatus status = RayStatus::Alive;
  std::uint32_t last_surface = kNoSurface;  ///< index of the last surface hit
  /// Accumulated power-normalised PRT matrix (ADR 0021): |P E|^2 is the power fraction for an
  /// incident state E with |E| = 1; P k_0 = k with k_0 the initial direction.
  math::CMat3 prt = math::CMat3::Identity();
  /// Power for an unpolarized source (source normalised to 1): weight = s ||P_T||^2 / 2 with
  /// P_T = P - k k_0^T and s the polarization-independent factors (ADR 0021).
  double weight = 1.0;
};

/// Physical inputs of one event, independent of the path (ADR 0021; the sequential tracer
/// fills them from the compiled path, the non-sequential tracer (M9) from the geometry).
struct EventMedia {
  math::Complex before{1.0};  ///< complex index of the medium the ray is in before the event
  math::Complex after{1.0};   ///< complex index after the event (equal to before unless crossing)
  /// Complex index on the other side of the surface (Fresnel and coating partner for Reflect;
  /// equal to after for Refract).
  math::Complex beyond{1.0};
  /// Vacuum wavelength, um; must be > 0 if `before` absorbs (Im > 0) or `layers` is not empty.
  /// No default wavelength, so that a forgotten value cannot pass unnoticed.
  double wavelength_um = 0.0;
  /// Coating layers at this wavelength in the order seen by the ray: from the incident medium to
  /// the other side (the caller reverses the stack for light from the substrate, ADR 0019).
  /// Empty if the surface has no coating.
  std::span<const coating::Layer<double>> layers;
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

/// Tolerance of the aperture check, mm: a hit counts as inside if it lies about this far outside
/// a rim (outer and inner rims; exact for circles and rectangle edges, see inside_aperture()).
/// Equal to the convergence limit of ray aiming (rtt::trace::kAimTolerance), so rays aimed exactly
/// at a stop rim pass (decided for #50).
inline constexpr double kApertureTolerance = 1e-9;

/// True if the local hit point lies inside the surface aperture or the surface has no aperture;
/// rims are inclusive with kApertureTolerance. Circular (r <= R + tol, with an inner radius
/// r >= R_i - tol), rectangular (|x| <= a + tol, |y| <= b + tol, i.e. up to sqrt(2) tol at the
/// corners) and elliptical apertures (semi-axes a + tol, b + tol: the normal allowance is tol on
/// the axes and between 2 sqrt(ab) / (a + b) tol and tol elsewhere, never more than tol), all in
/// local x, y, mm.
[[nodiscard]] bool inside_aperture(const compile::CompiledSurface& surface,
                                   const SurfaceHit& hit) noexcept;

/// Moves the ray to the hit point: position (global), OPL += Re(n_before) t, last_surface, and
/// volume absorption weight *= exp(-4 pi kappa t / lambda_vac) with kappa = Im(n_before), t the
/// geometric path in mm and lambda_vac in mm (Byrnes, arXiv:1603.02720v5, Eqs. (1), (2):
/// E ~ exp(i k.r) with |k| = 2 pi n / lambda_vac, and power ~ |E|^2 by Eqs. (17), (18); kappa is
/// the dimensionless extinction coefficient, not the absorption coefficient alpha = 4 pi kappa /
/// lambda). Direction, P and status are unchanged.
/// @pre hit.status == Hit
[[nodiscard]] RayState move_to_hit(const RayState& ray,
                                   const compile::CompiledSurface& surface,
                                   const SurfaceHit& hit,
                                   std::uint32_t surface_index,
                                   const EventMedia& media) noexcept;

/// move_to_hit() for a non-absorbing medium of real index n_before (M1 form).
[[nodiscard]] RayState move_to_hit(const RayState& ray,
                                   const compile::CompiledSurface& surface,
                                   const SurfaceHit& hit,
                                   std::uint32_t surface_index,
                                   double n_before) noexcept;

/// Executes `kind` at the hit (docs/architecture.md: apply_event(ray, hit, kind)) and applies the
/// interaction of the surface to P and weight (ADR 0021).
///
/// A ray that is not Alive is returned unchanged; a hit with status Missed or NoConvergence
/// gives that status and leaves the ray unchanged. Otherwise the ray moves to the hit point
/// (move_to_hit, with volume absorption) and then: Absorber interaction -> Absorbed with
/// weight 0; Diffract, Ordinary, Extraordinary -> EventImpossible (M4); Refract -> refracted
/// direction (Snell with Re(n)) or Tir beyond the critical angle; Reflect -> reflected direction;
/// Transmit -> unchanged direction. A stopped ray stays at the hit point with last_surface =
/// surface_index. The aperture is not checked here.
///
/// Interactions (ADR 0021), P_event power-normalised (rtt/polar/interface.hpp), P = P_event P,
/// weight scaled by ||P_T,new||^2 / ||P_T,old||^2 with ||P_T||^2 = ||P||_F^2 - 1:
/// - Fresnel: Refract between before and after; Reflect against beyond. A Mirror without
///   material (beyond = before) reflects as an ideal conductor, r_s = -1, r_p = +1; any other
///   surface with beyond = before gives r = 0 (vanishing reflection, weight 0, status Alive).
/// - IdealMirror: Reflect (-1, +1). IdealAntiReflection: Refract (1, 1), Reflect (0, 0)
///   (vanishing reflection: weight 0, status Alive).
///   IdealBeamSplitter: Reflect (-sqrt(R_s), +sqrt(R_p)), Refract and Transmit (sqrt(1 - R_s),
///   sqrt(1 - R_p)). CoatingRef: Refract or Reflect with the stack of `media.layers`
///   (rtt-coating, Byrnes). IdealPolarizer and IdealRetarder: Transmit only, axis
///   CompiledSurface::ideal_axis projected perpendicular to the ray (rtt/polar/ideal.hpp).
/// - Transmit at a Fresnel, IdealAntiReflection or CoatingRef surface is a dummy passage: P and
///   weight are unchanged.
/// - Any other combination of interaction and event kind, a CoatingRef surface without compiled
///   coating (CompiledSurface::coating empty), an ideal axis parallel to the ray and non-finite
///   amplitudes (grazing incidence, a coating layer exactly at q = 0) give EventImpossible.
/// Phase layers are ignored before M4.
///
/// Never throws; physical problems are status flags (ADR 0009).
/// @param ray           incoming ray in global coordinates, |dir| = 1
/// @param surface       compiled surface (pose and interaction are used)
/// @param hit           result of intersect_surface(ray, surface)
/// @param surface_index index of `surface` in CompiledSystem::surfaces(), stored in last_surface
/// @param kind          event at this surface
/// @param media         complex indices, wavelength and coating layers of this event
[[nodiscard]] RayState apply_event(const RayState& ray,
                                   const compile::CompiledSurface& surface,
                                   const SurfaceHit& hit,
                                   std::uint32_t surface_index,
                                   model::EventKind kind,
                                   const EventMedia& media) noexcept;

/// apply_event() with real indices and no coating layers (M1 form): before = n_before,
/// after = beyond = n_after, no absorption. For Reflect pass n_after = the index beyond the
/// surface; with n_after = n_before a Fresnel reflection at a non-mirror surface has r = 0
/// (weight 0). A CoatingRef surface acts as a bare interface here (no layers).
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
                                       const EventMedia& media) noexcept;

/// sequential_step() with real indices and no coating layers (M1 form, as apply_event()).
[[nodiscard]] RayState sequential_step(const RayState& ray,
                                       const compile::CompiledSurface& surface,
                                       std::uint32_t surface_index,
                                       model::EventKind kind,
                                       double n_before,
                                       double n_after) noexcept;

}  // namespace rtt::trace
