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
///   pupil, +y is the meridional direction. With a finite object py > 0 is a ray rising in object
///   space against the chief ray (sigma, ADR 0007 addendum #35 B; aim_ray).
/// - Positions in mm and unit directions in global coordinates (right-handed, optical axis +z).

#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <utility>
#include <variant>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/trace/apply_event.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/run_control.hpp"

namespace rtt::trace {

/// How a ray is aimed into the stop.
enum class Aiming : std::uint8_t {
  /// Newton iteration on the real ray until it hits the target point on the stop surface
  /// (|residual| < kAimTolerance), starting from the paraxial solution. If that run fails (strong
  /// pupil aberration, e.g. a stop deep inside the system), a fallback follows the target by
  /// continuation from the real chief ray (#287, ADR 0007 addendum); rays that converge without
  /// it are unchanged by it.
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

/// Gaussian quadrature in the pupil (#168, ADR 0030 point 4): `rings` rings at the Gauss-Legendre
/// nodes in u = rho^2 and `arms` evenly spaced arms per ring. Ring k has radius
/// rho_k = sqrt((1 + x_k) / 2), x_k the k-th node of the `rings`-point Gauss-Legendre rule on
/// [-1, 1] in ascending order (NIST DLMF, Eqs. (3.5.15), (3.5.21)); arm j has the angle
/// theta_j = 2 pi j / arms from +y towards +x, so (px, py) = (rho_k sin theta_j, rho_k cos
/// theta_j). Order: ring-major, rings with ascending rho. The quadrature weights are
/// gauss_pupil_weights(); with them the mean of f over the unit disk is exact for f = u^p e^{i m
/// theta} with p <= 2 rings - 1 and |m| < arms. The analyses that average over their rays without
/// these weights (spot, ghost ranking, path analyses) reject this sampling.
struct GaussPupil {
  int rings = 3;  ///< >= 1
  int arms = 6;   ///< >= 1
};

/// Normalised pupil coordinates.
struct PupilPoint {
  double px = 0.0;
  double py = 0.0;
};

/// An explicit list of pupil points, in this order (ADR 0032, point 2): for samplings computed
/// elsewhere, e.g. the cells of the pupil function. Each point is aimed like a SinglePupilPoint;
/// points outside the unit circle are allowed and aimed like any other point (whether a ray is
/// clipped, e.g. by a stop aperture, is left to the tracer). The analyses that reject GaussPupil
/// accept this list.
struct PointListPupil {
  std::vector<PupilPoint> points;  ///< non-empty, finite px and py
};

/// Pupil sampling of one field. New alternatives are appended (the indices of the others stay).
using PupilSampling = std::variant<SinglePupilPoint,
                                   HexapolarPupil,
                                   GridPupil,
                                   FanXPupil,
                                   FanYPupil,
                                   RandomPupil,
                                   GaussPupil,
                                   PointListPupil>;

/// Points of a pupil sampling in the documented order.
/// @throws std::invalid_argument for rings < 0 (hexapolar), n < 1, rings or arms < 1 (gauss), or
///         an empty point list or a point with a non-finite coordinate (PointListPupil; the
///         message gives the index)
[[nodiscard]] std::vector<PupilPoint> pupil_points(const PupilSampling& sampling);

/// Quadrature weights of the points of `gauss` in the order of pupil_points(): q = w_k / (2 arms)
/// for every point of ring k, w_k the Gauss-Legendre weight of node x_k (sum of w_k = 2), so the
/// weights sum to 1 and sum q f(p) approximates the mean of f over the unit disk. Dimensionless.
/// They are quadrature weights, not ray powers: make_rays() leaves the ray weight at the power
/// (ADR 0021). Nodes and weights from the monic Legendre recurrence (DLMF (3.5.30), (3.5.33_1))
/// and w_k = 1 / sum_{m < n} q_m(x_k)^2 with the orthonormal q_m (derived from DLMF (3.5.30_5)
/// to (3.5.32); docs/quellen.md).
/// @throws std::invalid_argument for rings or arms < 1
[[nodiscard]] std::vector<double> gauss_pupil_weights(const GaussPupil& gauss);

/// Result of aiming one ray.
struct AimedRay {
  /// Start of the ray in global coordinates (OPL 0). Status Alive; NoConvergence if real aiming
  /// failed; Vignetted with last_surface = the stop surface if it failed and the target lies
  /// outside the outer rim of the stop aperture, so that no ray can pass there (ADR 0007,
  /// addendum #310). A failed ray holds the last iterate.
  RayState ray;
  /// Newton steps taken (0 for paraxial aiming); with the continuation fallback of real aiming
  /// all steps of all runs, so it can exceed kMaxAimIterations (#287). For a failed ray
  /// (NoConvergence or Vignetted) the steps of the run from the paraxial solution.
  int iterations = 0;
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
/// Pupil points with a finite object are oriented at the object-space slope (ADR 0007, addendum
/// #35 B; replaces the orientation at the reference wavelength of #96): py > 0 is always a ray
/// rising in object space, for every wavelength and on both sides of an entrance pupil at
/// infinity, so the rays change continuously when the EP passes through infinity. For a virtual
/// EP on the far side of the object (z_ep < z_obj) the EP point of (px, py) is therefore
/// -(px, py) r_ep. With the object at infinity (px, py) is the EP point (#8).
/// The target is sigma (px R_s, py R_s) in the local coordinates of the stop surface (first Stop
/// event of the path), sigma = -1 for a virtual EP beyond the object, else +1, and R_s = paraxial
/// stop radius belonging to the entrance pupil (negative where the stop sees the EP inverted).
/// Apertures are ignored while aiming; vignetting is left to the tracer.
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

/// Radius of the stop belonging to pupil coordinates |(px, py)| = 1, mm, > 0: |R_s| of the aim
/// target sigma (px R_s, py R_s) of aim_ray() at `wavelength` (paraxial: the entrance pupil
/// radius times the paraxial stop height per unit entrance-pupil height; for an object-space
/// telecentric system the marginal slope times the stop height per unit slope). With real
/// aiming a ray of pupil point p meets the stop at |p| |R_s| within kAimTolerance, so a centred
/// obscuration of radius r_i at the stop is the pupil circle |p| = r_i / |R_s| (ADR 0032,
/// addendum #229).
/// @throws as aim_ray() for the path, the wavelength, the stop and the entrance pupil
[[nodiscard]] double stop_radius(const compile::CompiledSystem& system,
                                 compile::PathId path,
                                 std::uint16_t wavelength);

/// Rays for the given fields and one wavelength: for every field (outer loop) every pupil point
/// (inner loop) of `sampling`, aimed with `aiming`. Sets pos, dir, wl, field, pupil_x, pupil_y,
/// status (Alive, NoConvergence or Vignetted, as AimedRay) and last_surface (kNoSurface, the stop
/// surface for a Vignetted ray, #310); OPL 0, weight 1, P = identity.
/// The rays are aimed in parallel with oneTBB (#119), each from its own paraxial start value;
/// the result is bitwise the same for every number of threads (ADR 0004, addendum #119).
/// @throws as aim_ray()
[[nodiscard]] RayBatch make_rays(const compile::CompiledSystem& system,
                                 compile::PathId path,
                                 std::span<const std::uint16_t> fields,
                                 std::uint16_t wavelength,
                                 const PupilSampling& sampling,
                                 Aiming aiming = Aiming::Real);

/// make_rays() with cancellation and progress (run_control.hpp, #83): stage "aim", one work
/// item per ray, in parallel blocks of control.block_size rays as in SequentialTracer::trace()
/// (#119). With an empty `control` it is make_rays() exactly; the rays never depend on the
/// control.
/// @throws as make_rays(), Cancelled, or the exception of the progress callback
[[nodiscard]] RayBatch make_rays(const compile::CompiledSystem& system,
                                 compile::PathId path,
                                 std::span<const std::uint16_t> fields,
                                 std::uint16_t wavelength,
                                 const PupilSampling& sampling,
                                 Aiming aiming,
                                 const RunControl& control);

namespace detail {
struct AimData;  // sources.cpp
}  // namespace detail

/// Paraxial aiming data of one path at one wavelength, computed from a compiled system (ADR 0007,
/// addendum #252; ADR 0036, point 2): entrance pupil (position and radius, or the marginal slope
/// of an object-space telecentric system), stop scale, orientation sigma of the pupil labels, the
/// stop event and the data of the field conversion at the reference wavelength. For the
/// perturbed trials of a tolerance analysis, which are not rotationally symmetric, these data
/// come from the nominal system (rtt::paraxial::first_order rejects the trials). An immutable
/// value: cheap to copy, safe to read from many threads. Created with aim_reference().
class AimReference {
 public:
  /// Path the data belong to.
  [[nodiscard]] compile::PathId path() const noexcept;
  /// Index of the wavelength the data belong to.
  [[nodiscard]] std::uint16_t wavelength() const noexcept;
  /// Internal data (sources.cpp); not part of the API.
  [[nodiscard]] const detail::AimData& data() const noexcept { return *data_; }

 private:
  friend AimReference aim_reference(const compile::CompiledSystem& system,
                                    compile::PathId path,
                                    std::uint16_t wavelength);
  explicit AimReference(std::shared_ptr<const detail::AimData> data) : data_(std::move(data)) {}
  std::shared_ptr<const detail::AimData> data_;
};

/// The aiming data of `path` at `wavelength` of `system` (usually the nominal system), computed
/// once: the same paraxial computation that aim_ray() and make_rays() do without a reference.
/// @throws as aim_ray() for the path, the wavelength, the stop and the entrance pupil (also
///         rtt::paraxial::ParaxialError for a path that is not rotationally symmetric)
[[nodiscard]] AimReference aim_reference(const compile::CompiledSystem& system,
                                         compile::PathId path,
                                         std::uint16_t wavelength);

/// aim_ray() with the aiming data of `reference` instead of the paraxial data of `system`
/// (ADR 0007, addendum #252; ADR 0036, point 2). The entrance pupil, the stop scale, sigma, the
/// Newton start and the field conversion come from the reference, so a field point means the
/// same direction or object point as in the reference system. The ray is traced through
/// `system`, the target sigma (px R_s, py R_s) lies in the local coordinates of the stop surface
/// of `system`, and the start plane (object at infinity) lies before the surfaces of `system`
/// (bundle reserve of 3 entrance-pupil radii around the reference chief ray, and 1 mm; a
/// decenter well below the pupil radius stays inside it). Real aiming thus converges onto a
/// perturbed stop: the chief ray passes through its centre. With the reference of `system`
/// itself the result is bitwise that of aim_ray() without reference.
/// @throws std::invalid_argument if `reference` does not fit: another path or wavelength than
///         `path` and `wavelength`; in `system` another name, number of events or surface id,
///         kind or order of an event of the path; other wavelength values or another reference
///         wavelength; other field type or field values (x, y); another object (at infinity,
///         object_z(), ADR 0038) or system aperture (type, value, stop surface; ADR 0039) (ADR
///         0007, addendum #252, precision #255);
///         and as aim_ray() for px, py, the field and the start plane
[[nodiscard]] AimedRay aim_ray(const compile::CompiledSystem& system,
                               compile::PathId path,
                               std::uint16_t field,
                               std::uint16_t wavelength,
                               double px,
                               double py,
                               Aiming aiming,
                               const AimReference& reference);

/// aim_ray() for an arbitrary field value with the aiming data of `reference` (see above).
/// @throws as the overload above, and std::invalid_argument if a field value is not finite
[[nodiscard]] AimedRay aim_ray(const compile::CompiledSystem& system,
                               compile::PathId path,
                               const model::Field& field,
                               std::uint16_t wavelength,
                               double px,
                               double py,
                               Aiming aiming,
                               const AimReference& reference);

/// |R_s| of the aim target with the data of `reference`, mm (see stop_radius() above); equal to
/// stop_radius(system, path, wavelength) of the system the reference was computed from.
[[nodiscard]] double stop_radius(const AimReference& reference);

/// make_rays() with the aiming data of `reference` (see aim_ray() with a reference); bitwise the
/// same for every number of threads.
/// @throws as aim_ray() with a reference
[[nodiscard]] RayBatch make_rays(const compile::CompiledSystem& system,
                                 compile::PathId path,
                                 std::span<const std::uint16_t> fields,
                                 std::uint16_t wavelength,
                                 const PupilSampling& sampling,
                                 Aiming aiming,
                                 const AimReference& reference);

/// make_rays() with the aiming data of `reference`, cancellation and progress (see above).
/// @throws as make_rays() with a reference, Cancelled, or the exception of the progress callback
[[nodiscard]] RayBatch make_rays(const compile::CompiledSystem& system,
                                 compile::PathId path,
                                 std::span<const std::uint16_t> fields,
                                 std::uint16_t wavelength,
                                 const PupilSampling& sampling,
                                 Aiming aiming,
                                 const AimReference& reference,
                                 const RunControl& control);

}  // namespace rtt::trace
