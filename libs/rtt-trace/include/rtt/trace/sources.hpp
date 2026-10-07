#pragma once

/// @file sources.hpp
/// Rays from the system definition: field and aperture types, ray aiming and pupil sampling
/// (docs/architecture.md, Engine 2; conventions in "Konventionen", Feldwinkel).
///
/// Conventions:
/// - Field angle (theta_x, theta_y) in degree: chief-ray direction in object space
///   d = (tan theta_x, tan theta_y, 1) normalised; theta_y > 0 means the ray rises towards +y.
/// - Object height (x, y) in mm: object point (x, y, -object.distance) in global coordinates;
///   positive height = +y.
/// - Paraxial image height (x', y') in mm: converted into a field angle (object at infinity) or
///   an object height (finite object) with the paraxial chief ray (rtt-paraxial), then aimed.
/// - Field values are converted into a direction or an object point at the reference wavelength
///   (paraxial data of the reference wavelength, e.g. the entrance pupil for an angle with a
///   finite object), so a field point is the same physical point for every wavelength. Pupil,
///   stop target and start plane belong to the wavelength of each ray (#31, fix of #8).
/// - Normalised pupil coordinates (px, py): the unit circle is the rim of the paraxial entrance
///   pupil, +y is the meridional direction.
/// - Positions in mm and unit directions in global coordinates (right-handed, optical axis +z).

#include <cstddef>
#include <cstdint>
#include <span>
#include <variant>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/trace/apply_event.hpp"
#include "rtt/trace/ray_batch.hpp"

namespace rtt::trace {

/// How a ray is aimed into the stop.
enum class Aiming : std::uint8_t {
  /// Newton iteration on the real ray until it hits the target point on the stop surface
  /// (|residual| < kAimTolerance), starting from the paraxial solution.
  Real,
  /// Straight line through the target point on the paraxial entrance pupil, no iteration.
  Paraxial,
};

/// Convergence limit of real aiming: distance between the hit and the target point on the stop
/// surface, in local stop coordinates, mm (issue #8).
inline constexpr double kAimTolerance = 1e-9;
static_assert(kAimTolerance == kApertureTolerance,
              "rays aimed at a stop rim must pass the aperture check (#50)");

/// Maximum number of Newton steps of real aiming before the ray gets NoConvergence.
inline constexpr int kMaxAimIterations = 20;

/// Step of the central difference for the aiming Jacobian, relative to the paraxial stop radius
/// R_s that belongs to the entrance pupil (ADR 0007): h = kAimStepRelative * R_s.
inline constexpr double kAimStepRelative = 1e-6;

/// Lower bound of that step relative to the pupil radius r_ep in pupil units (mm on the EP
/// plane, or the marginal slope u_m for an object-space telecentric system):
/// h = max(kAimStepRelative |R_s|, kAimStepPupilFloor r_ep) (ADR 0007, addendum #96). The step
/// is sized in stop units but applied to pupil coordinates; for an entrance pupil far away,
/// |R_s / r_ep| -> 0, the first term falls below the resolution of the pupil coordinates. The
/// floor acts only for |R_s / r_ep| < kAimStepPupilFloor / kAimStepRelative = 1e-4. For an
/// object-space telecentric system the step is h = kAimStepRelative u_m in slope units.
inline constexpr double kAimStepPupilFloor = 1e-10;

/// One ray at normalised pupil coordinates (px, py).
struct SinglePupilPoint {
  double px = 0.0;
  double py = 0.0;
};

/// Centre ray plus rings k = 1..rings with 6k rays on radius k / rings, angles 2 pi j / (6k)
/// starting at +y and turning towards +x.
struct HexapolarPupil {
  int rings = 6;
};

/// n x n points evenly spaced on [-1, 1]^2 (rows from py = -1 to 1, px fastest), only points with
/// px^2 + py^2 <= 1 (rim included). n = 1 gives the centre.
struct GridPupil {
  int n = 11;
};

/// n points evenly spaced on px in [-1, 1] at py = 0. n = 1 gives the centre.
struct FanXPupil {
  int n = 11;
};

/// n points evenly spaced on py in [-1, 1] at px = 0. n = 1 gives the centre.
struct FanYPupil {
  int n = 11;
};

/// `count` points uniformly distributed in the unit disk. Generator std::mt19937_64 seeded with
/// `seed`; per ray two draws u1, u2 in [0, 1) as (x >> 11) * 2^-53 (in this order), then
/// r = sqrt(u1), phi = 2 pi u2, (px, py) = (r sin phi, r cos phi). The draws are identical on
/// all platforms; px and py may differ by a few ulp because std::sin and std::cos are not
/// correctly rounded.
struct RandomPupil {
  std::size_t count = 100;
  std::uint64_t seed = 0;
};

/// Pupil sampling of one field.
using PupilSampling =
    std::variant<SinglePupilPoint, HexapolarPupil, GridPupil, FanXPupil, FanYPupil, RandomPupil>;

/// Normalised pupil coordinates.
struct PupilPoint {
  double px = 0.0;
  double py = 0.0;
};

/// Points of a pupil sampling in the documented order.
/// @throws std::invalid_argument for rings < 0 or n < 1
[[nodiscard]] std::vector<PupilPoint> pupil_points(const PupilSampling& sampling);

/// Result of aiming one ray.
struct AimedRay {
  /// Start of the ray in global coordinates (OPL 0). Status Alive, or NoConvergence if real
  /// aiming failed (the ray then holds the last iterate).
  RayState ray;
  int iterations = 0;  ///< Newton steps taken (0 for paraxial aiming)
  /// Distance to the target on the stop surface, mm (real aiming; 0 for paraxial aiming;
  /// +infinity if no iterate reached the stop).
  double residual = 0.0;
};

/// Aims one ray of field `field` at normalised pupil coordinates (px, py).
///
/// Object at infinity: the direction is the field direction; the ray starts on a plane
/// perpendicular to it, placed per field so that the bundle (up to 3 EP radii around the chief
/// ray) starts at least 1 mm before the entrance pupil and before every surface of the path,
/// sampled within its aperture and shape domain (OPL starts on this plane wave). A field angle
/// with a finite object places the object point on the chief ray through the EP centre at the
/// reference wavelength.
/// Finite object: the ray starts in the object point (OPL 0) and runs along the line through
/// it and the EP point into the system (+z), also for a virtual EP on the far side of the
/// object (z_ep < z_obj, #93). Object-space telecentric (EP at infinity, #96): (px, py) are
/// object-space slopes, d ~ (px u_m, py u_m, 1) with the paraxial marginal slope u_m (NA / n
/// for object_na, r_stop / |s| for stop_size); the paraxial chief ray is parallel to the axis.
/// Pupil points are oriented at the reference wavelength: if the EP of `wavelength` lies on the
/// other side of the object than that of the reference wavelength, (px, py) is used as
/// -(px, py), so it reaches the same side of the stop (#96).
/// The target is (px R_s, py R_s) in the local coordinates of the stop surface (first Stop event
/// of the path), R_s = paraxial stop radius belonging to the entrance pupil. Apertures are
/// ignored while aiming; vignetting is left to the tracer.
/// @throws std::invalid_argument if the path, field or wavelength does not exist; px or py is
///         not finite; the field type does not fit the object (object height at infinity); a
///         field angle is not in (-90, 90) degree; no start plane exists before an unbounded
///         surface (no aperture and an unbounded shape domain, e.g. a paraboloid, or an
///         aperture of 1 km or more) that curves back against a steep field (object at
///         infinity); a paraxial image height is requested without a finite paraxial image; or
///         the entrance pupil is not defined at the ray's wavelength (pupil at infinity with the
///         object at infinity, no diameter for this aperture type, pupil in the object plane;
///         object-space telecentric with an entrance pupil diameter or image F-number as system
///         aperture) or, for paraxial image heights and angles with a finite object, at the
///         reference wavelength (pupil in the object plane; a field angle with a finite object
///         and the pupil at infinity)
/// @throws rtt::compile::NoStopError (a std::invalid_argument) if the path has no stop; checked
///         by rtt::compile::require_stop before any ray is traced (ADR 0022)
/// @throws rtt::paraxial::ParaxialError if the path is not rotationally symmetric or the stop
///         aperture is not circular
[[nodiscard]] AimedRay aim_ray(const compile::CompiledSystem& system,
                               compile::PathId path,
                               std::uint16_t field,
                               std::uint16_t wavelength,
                               double px,
                               double py,
                               Aiming aiming = Aiming::Real);

/// Aims one ray of an arbitrary field value, e.g. for field sweeps between the model's field
/// points (#31). `field` is interpreted with the field type of the system
/// (CompiledSystem::fields().type) exactly like a field point of the model; with a value equal
/// to field point k the result is identical to aim_ray(system, path, k, ...).
/// @param field field value: angle in degree, object height in mm or paraxial image height in
///              mm, as given by the system's field type (Field::weight is ignored)
/// @throws as the overload with a field index (except the field-index check), and
///         std::invalid_argument if a field value is not finite
[[nodiscard]] AimedRay aim_ray(const compile::CompiledSystem& system,
                               compile::PathId path,
                               const model::Field& field,
                               std::uint16_t wavelength,
                               double px,
                               double py,
                               Aiming aiming = Aiming::Real);

/// Rays for the given fields and one wavelength: for every field (outer loop) every pupil point
/// (inner loop) of `sampling`, aimed with `aiming`. Sets pos, dir, wl, field, pupil_x, pupil_y,
/// status (Alive or NoConvergence); OPL 0, weight 1, P = identity, last_surface = kNoSurface.
/// @throws as aim_ray()
[[nodiscard]] RayBatch make_rays(const compile::CompiledSystem& system,
                                 compile::PathId path,
                                 std::span<const std::uint16_t> fields,
                                 std::uint16_t wavelength,
                                 const PupilSampling& sampling,
                                 Aiming aiming = Aiming::Real);

}  // namespace rtt::trace
