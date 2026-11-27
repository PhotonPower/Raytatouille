#pragma once

/// @file scene.hpp
/// Scene of the non-sequential tracer (ADR 0057, points 2 and 3): the active elements of a
/// compiled system as faces with a medium on each side, closed bodies with generated mantles,
/// and a BVH for the nearest hit.
///
/// Conventions (docs/architecture.md, Konventionen): lengths in mm, global coordinates
/// right-handed with the optical axis +z. Every face has an orientation: its normal points to
/// the positive side (for a surface the side +z at the vertex, as in rtt-geom; for a generated
/// face away from the glass). A ray with d . n < 0 arrives from the positive side. The ray
/// carries no medium; the side it arrives from decides the media of an event (ADR 0057,
/// point 2).

#include <array>
#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/math/types.hpp"
#include "rtt/model/validate.hpp"
#include "rtt/nonseq/aabb.hpp"
#include "rtt/nonseq/bvh.hpp"
#include "rtt/trace/apply_event.hpp"

namespace rtt::nonseq {

/// Index of no face.
inline constexpr std::uint32_t kNoFace = std::numeric_limits<std::uint32_t>::max();

/// What a face is (ADR 0057, point 2).
enum class FaceKind : std::uint8_t {
  Body,      ///< surface of a Lens or Plate; the media of the body on its two sides
  Mantle,    ///< generated mantle of a segment: a wall over the aperture rim
  Shoulder,  ///< generated flat ring of a segment between unequal circular rims
  Thin,      ///< surface of a ThinElement; its interaction on both sides
  Mirror,    ///< surface of a Mirror; the front side reflects, the back side absorbs
  Stop,      ///< surface of a Stop; open inside the aperture, an absorbing screen around it
  Detector,  ///< surface of a Detector
};

/// Cross-section of a generated wall in the frame of its segment's first surface (local x, y).
enum class WallShape : std::uint8_t {
  Circle,   ///< x^2 + y^2 = a^2
  Ellipse,  ///< (x / a)^2 + (y / b)^2 = 1
  Edge,     ///< the strip over the segment p0 -> p1
};

/// Geometry of a generated face (Mantle or Shoulder), in the frame of Scene::surfaces()[surface].
struct GeneratedGeometry {
  WallShape wall = WallShape::Circle;  ///< Mantle only
  double a = 0.0;                      ///< circle radius or semi-axis x, mm
  double b = 0.0;                      ///< ellipse semi-axis y, mm
  std::array<double, 2> p0{};          ///< edge start (x, y), mm
  std::array<double, 2> p1{};          ///< edge end (x, y), mm
  /// Mantle: +1 if the normal points away from the axis (outer mantle, edges of the outer rim),
  /// -1 towards it (bore). Shoulder: +1 if the normal is +z of the frame, -1 if -z.
  double orientation = 1.0;
  /// Mantle: true if a hit counts for z in [z_lo, z_hi] (a step between unequal circles, whose
  /// rim heights are constant); false if it counts between the two surfaces of the segment
  /// (ADR 0057, point 2), pointwise.
  bool constant_z = false;
  double z_lo = 0.0;   ///< mm, frame z
  double z_hi = 0.0;   ///< mm, frame z; Shoulder: z_lo = z_hi = height of the ring
  double r_in = 0.0;   ///< Shoulder: inner radius of the ring, mm
  double r_out = 0.0;  ///< Shoulder: outer radius of the ring, mm
};

/// One face of the scene.
struct Face {
  FaceKind kind = FaceKind::Body;
  /// Index into Scene::surfaces(): frame, shape, aperture and interaction. A generated face has a
  /// surface of its own with the frame of its segment's first surface, no aperture and the rim
  /// interaction (Fresnel or Absorber).
  std::uint32_t surface = 0;
  /// Stop only: index into Scene::surfaces() of the screen, a copy of the stop surface with the
  /// aperture scaled by SceneOptions::screen_scale about the local origin, without inner radius
  /// and obscurations.
  std::uint32_t screen = 0;
  /// Index into CompiledSystem::surfaces(): the surface itself; for a generated face the first
  /// surface j of its segment.
  std::uint32_t system_surface = 0;
  std::uint32_t element = 0;  ///< index into CompiledSystem::elements()
  std::uint32_t segment = 0;  ///< generated face: the segment j; otherwise 0
  /// Generated face: the scene surfaces of the segment's two surfaces j and j + 1.
  std::uint32_t lower = 0;
  std::uint32_t upper = 0;
  /// Generated face: +1 if the segment lies on the +z side of surface j resp. j + 1 (their own
  /// frames), -1 otherwise.
  double lower_side = 1.0;
  double upper_side = 1.0;
  /// Media on the positive and the negative side, indices into CompiledSystem::media().
  std::uint32_t medium_positive = 0;
  std::uint32_t medium_negative = 0;
  /// Mirror with a substrate: its medium (the Fresnel partner of the front reflection); none
  /// otherwise. The substrate is not a body (ADR 0057, point 2).
  std::optional<std::uint32_t> substrate;
  bool transparent = false;     ///< Detector switched to transparent (ADR 0057, point 2)
  GeneratedGeometry generated;  ///< Mantle and Shoulder only
};

/// How a hit acts (ADR 0057, points 2 and 3); regions where a face is not there (outside the
/// aperture, the opening of a stop) are never hits.
enum class HitRegion : std::uint8_t {
  Clear,     ///< the event of the face happens
  Absorbed,  ///< inside an obscuration or on the screen of a stop: the ray is absorbed there
};

/// The nearest hit of a ray in the scene.
struct SceneHit {
  /// Hit; Missed if no face is hit (the ray escapes); NoConvergence if a Newton intersection did
  /// not converge on a face that could be nearer than every hit found (ADR 0057, Nachtrag PR 2, g).
  geom::HitStatus status = geom::HitStatus::Missed;
  std::uint32_t face = kNoFace;
  HitRegion region = HitRegion::Clear;
  /// The hit in the frame of Scene::surfaces()[face.surface]: t (mm), point, unit normal along
  /// the face orientation, incident direction; for apply_event (rtt/trace/apply_event.hpp).
  trace::SurfaceHit local;
  math::Vec3 point = math::Vec3::Zero();   ///< global, mm
  math::Vec3 normal = math::Vec3::Zero();  ///< global unit normal, face orientation
  /// True if the ray arrives from the positive side (d . n < 0).
  bool from_positive = false;
  std::uint32_t medium_before = 0;  ///< medium of the side the ray arrives from
  std::uint32_t medium_beyond = 0;  ///< medium of the other side
};

/// Options of the scene from the section `nonsequential` (ADR 0057, point 9; until the format
/// exists, filled by the caller).
struct SceneOptions {
  /// Outer size of the screen of a Stop as a multiple of its aperture (ADR 0057, point 2), > 1.
  double screen_scale = 2.0;
  /// Elements (indices into CompiledSystem::elements()) whose mantles absorb instead of Fresnel.
  std::vector<std::uint32_t> absorbing_rims;
  /// Detector elements switched to transparent.
  std::vector<std::uint32_t> transparent_detectors;
};

/// Thrown by Scene::build() for a system the scene cannot model; diagnostics() holds the errors
/// with codes nonseq.* and JSON pointers into the system file (ADR 0022).
class SceneError : public std::invalid_argument {
 public:
  explicit SceneError(std::vector<model::Diagnostic> diagnostics);
  [[nodiscard]] const std::vector<model::Diagnostic>& diagnostics() const noexcept {
    return diagnostics_;
  }

 private:
  std::vector<model::Diagnostic> diagnostics_;
};

/// Immutable scene; safe to read from many threads.
class Scene {
 public:
  /// Builds the scene of the active elements of `system` (ADR 0057, point 2; inactive elements,
  /// ADR 0049, are not part of it):
  /// - Lens, segmented Plate and Plate of one material with two surfaces (segment rule, ADR 0017;
  ///   ADR 0057, Nachtrag PR 2, b): the segment j between surfaces j and j + 1 lies on the side of
  ///   surface j given by the sign of n_j . (s_{j+1} - s_j), n_j the axis of surface j and s_j its
  ///   vertex; every segment gets a mantle over the aperture rim, with a step and a flat
  ///   shoulder between unequal circles of rotationally symmetric surfaces, and the same for an
  ///   inner radius.
  /// - Plate of one material with more than two surfaces (prism): convex, the inside is the side
  ///   of the centroid of all aperture corners; no mantle.
  /// - Mirror, ThinElement, Stop and Detector: thin faces with the environment on both sides.
  /// Faces with a finite box go into the BVH, the others are tested for every ray (ADR 0057,
  /// Nachtrag PR 2, c); a Stop without aperture has no face.
  /// @throws SceneError with nonseq.aperture_missing (a body surface without aperture),
  ///         nonseq.segment_thickness (centre thickness 0, or overlapping segments),
  ///         nonseq.rim_unsupported (mantle impossible: different aperture kinds, axes not
  ///         parallel and coaxial within 1e-12, different non-circular apertures, a step at a
  ///         surface that is not rotationally symmetric, a rim outside the domain of the shape)
  ///         and nonseq.plate_not_convex
  /// @throws std::invalid_argument for screen_scale <= 1 or not finite, an element index in the
  ///         options that is out of range or of the wrong kind, and from Bvh (a tree deeper than
  ///         its traversal stack, more than 65535 faces with one box centre)
  [[nodiscard]] static Scene build(const compile::CompiledSystem& system,
                                   const SceneOptions& options = {});

  [[nodiscard]] const std::vector<Face>& faces() const noexcept { return faces_; }
  /// Surfaces of the faces: copies of the compiled surfaces and the frames of generated faces.
  [[nodiscard]] const std::vector<compile::CompiledSurface>& surfaces() const noexcept {
    return surfaces_;
  }
  /// Global box of every face (the box of its local box, margin included); not finite for faces
  /// outside the BVH.
  [[nodiscard]] const std::vector<Aabb>& bounds() const noexcept { return bounds_; }
  [[nodiscard]] const Bvh& bvh() const noexcept { return bvh_; }
  /// Faces outside the BVH (no finite box), tested for every ray.
  [[nodiscard]] const std::vector<std::uint32_t>& unbounded() const noexcept { return unbounded_; }
  /// The environment medium, index into CompiledSystem::media().
  [[nodiscard]] std::uint32_t environment_medium() const noexcept { return environment_; }

  /// Nearest hit of `ray` (global pos and dir, |dir| = 1) over all faces (ADR 0057, point 3).
  /// Per face the nearest crossing with t > t_min, t_min = geom::kDefaultTMin for `last_face`
  /// (self-hit, Nachtrag PR 1, e) and 0 otherwise; a crossing where the face is not there (outside
  /// the aperture, outside the bounds of a mantle, the opening of a stop) is passed and the
  /// search on that face goes on from it. Never throws.
  [[nodiscard]] SceneHit closest_hit(const trace::RayState& ray,
                                     std::uint32_t last_face = kNoFace) const noexcept;

  /// The nearest hit of `ray` on face `face` with t in (t_min, t_max), or status Missed;
  /// NoConvergence as in closest_hit(). Used by closest_hit() and by the tests.
  [[nodiscard]] SceneHit intersect_face(const trace::RayState& ray,
                                        std::uint32_t face,
                                        double t_min,
                                        double t_max) const noexcept;

 private:
  Scene() = default;
  /// The next crossing of the Newton surface of face `face` after t1 and before t_max, by
  /// bracketing F = z - sag inside the face's local box (ADR 0057, Nachtrag PR 2, i).
  [[nodiscard]] trace::SurfaceHit next_crossing(const trace::RayState& ray,
                                                std::uint32_t face,
                                                double t1,
                                                double t_max) const noexcept;
  friend class SceneBuilder;

  std::vector<Face> faces_;
  std::vector<compile::CompiledSurface> surfaces_;
  std::vector<Aabb> bounds_;
  std::vector<Aabb> local_bounds_;  ///< box of every face in the frame of its surface
  Bvh bvh_;
  std::vector<std::uint32_t> unbounded_;
  std::vector<std::uint32_t> bounded_;  ///< face of each BVH primitive
  std::uint32_t environment_ = 0;
};

}  // namespace rtt::nonseq
