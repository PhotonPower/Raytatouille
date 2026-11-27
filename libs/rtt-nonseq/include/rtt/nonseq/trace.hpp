#pragma once

/// @file trace.hpp
/// Trace core of the non-sequential tracer (ADR 0057, points 4-8; PR 3 and 4): rays from sources
/// through the Scene, branching at every hit by Monte Carlo or full splitting, with a power
/// balance per source and wavelength, detector fields and a ray path report. Sources and
/// detectors are plain C++ structures until the format carries them (ADR 0057, point 12).
///
/// Conventions (docs/architecture.md, Konventionen): positions in mm and unit directions in
/// global coordinates, powers in W, wavelengths as indices into CompiledSystem::wavelengths_um().
/// The sources are unpolarized: a ray carries power = power of its start ray * weight, with the
/// weight of ADR 0021 starting at 1.

#include <cstdint>
#include <numbers>
#include <optional>
#include <span>
#include <string>
#include <variant>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/math/isometry.hpp"
#include "rtt/math/types.hpp"
#include "rtt/model/path.hpp"
#include "rtt/nonseq/scene.hpp"
#include "rtt/nonseq/source.hpp"
#include "rtt/trace/run_control.hpp"

namespace rtt::nonseq {

/// How a ray continues at a hit with several branches (ADR 0057, point 4).
enum class Splitting : std::uint8_t {
  /// One branch i, chosen with p_i = w_i / W; the ray keeps P_i and gets the weight W
  /// (ADR 0057, Nachtrag Review PR 1, a). Unbiased for every polarization state.
  MonteCarlo,
  /// Every branch with w_i > 0, depth first in the fixed order of the candidates.
  Full,
};

/// One start ray of an explicit list (RaySource).
struct StartRay {
  math::Vec3 position = math::Vec3::Zero();    ///< global, mm
  math::Vec3 direction = math::Vec3::UnitZ();  ///< global, normalised by trace()
  double power = 0.0;                          ///< W, finite and >= 0
  std::uint16_t wavelength = 0;                ///< index into CompiledSystem::wavelengths_um()
};

/// A source as an explicit list of start rays, e.g. the rays of a sequential trace; its index in
/// the list passed to trace() is the source number of the random counter (ADR 0057, point 5) and
/// of the balance, the index of a ray in `rays` its ray index.
struct RaySource {
  std::vector<StartRay> rays;
};

/// Options of a run (ADR 0057, points 4 and 5).
struct TraceOptions {
  Splitting splitting = Splitting::MonteCarlo;
  /// Key of the Philox generator (low and high 32 bits); only Monte Carlo draws.
  std::uint64_t seed = 0;
  /// A branch ends when its number of events exceeds this (ADR 0057, point 4).
  std::uint32_t max_depth = 64;
  /// A branch ends when its weight falls below this times the start weight 1; >= 0.
  double min_relative_power = 1e-6;
  /// Collect the ray path report (ADR 0057, point 8; Nachtrag PR 4, d).
  bool path_report = true;
  /// Paths with a smaller share of the source power are left out of the report; in [0, 1].
  double report_min_share = 1e-4;
  /// Record every hit on an irradiance detector in NonseqResult::hits.
  bool record_hits = false;
};

/// Irradiance detector (ADR 0057, point 7; Nachtrag PR 4, b): nx * ny pixels over the
/// rectangular aperture of a plane Detector face, counted from both sides.
struct IrradianceDetector {
  std::uint32_t face = kNoFace;  ///< index into Scene::faces(), a face of kind Detector
  std::uint32_t nx = 1;          ///< pixels along local x, > 0
  std::uint32_t ny = 1;          ///< pixels along local y, > 0
};

/// Polar intensity detector (far field, ADR 0057, point 7; Nachtrag PR 4, h): counts the rays
/// that hit no face any more (status Missed) by their direction, theta in [0, max_angle] from
/// the local +z axis of `pose` and phi = atan2(y, x) in [0, 2 pi) from its local +x axis.
struct IntensityDetector {
  /// Rotation local -> global, orthonormal (its transpose is the inverse); translation unused.
  math::Isometry3 pose;
  double max_angle = std::numbers::pi / 2.0;  ///< rad, in (0, pi]
  std::uint32_t n_theta = 1;                  ///< cells in theta, > 0
  std::uint32_t n_phi = 1;                    ///< cells in phi, > 0
};

using Detector = std::variant<IrradianceDetector, IntensityDetector>;

/// Field of one detector. Cell (w, r, c) is at (w * rows + r) * cols + c, w the wavelength index
/// of CompiledSystem::wavelengths_um().
struct DetectorField {
  /// Irradiance: ny rows along local y and nx columns along local x; pixel (r, c) covers
  /// x in [-h_x + c Delta x, -h_x + (c + 1) Delta x], likewise y. Intensity: n_theta rows,
  /// n_phi columns; cell (r, c) covers theta in [r Delta theta, (r + 1) Delta theta] and
  /// phi in [c Delta phi, (c + 1) Delta phi].
  std::uint32_t rows = 0;
  std::uint32_t cols = 0;
  /// Irradiance E = sum(Phi_i) / A_pixel in W/mm^2 (PBRT 4ed, §4.1, Eq. (4.1)), or intensity
  /// I = sum(Phi_i) / Delta Omega in W/sr with Delta Omega = 2 Delta phi sin((theta_2 -
  /// theta_1)/2) sin((theta_1 + theta_2)/2) (NBS TN 910-1, Eq. (A2-9)); summed over sources.
  std::vector<double> value;
  std::vector<std::uint64_t> hits;  ///< number of rays (branches) counted per cell
  std::vector<double> weights;      ///< sum of the weights (ADR 0021) per cell
  /// Irradiance: one entry, A_pixel in mm^2. Intensity: Delta Omega in sr per row.
  std::vector<double> cell_size;
};

/// One recorded hit on an irradiance detector (TraceOptions::record_hits).
struct DetectorHit {
  std::uint32_t detector = 0;  ///< index into the detectors passed to trace()
  std::uint32_t source = 0;
  std::uint64_t index = 0;                     ///< ray index within the source
  std::uint16_t wavelength = 0;                ///< index into CompiledSystem::wavelengths_um()
  math::Vec3 point = math::Vec3::Zero();       ///< global, mm
  math::Vec3 direction = math::Vec3::UnitZ();  ///< global, unit
  double power = 0.0;                          ///< W
};

/// How a ray path ends (ADR 0057, point 8; Nachtrag PR 4, d).
enum class PathEnd : std::uint8_t {
  Detected,     ///< absorbed by a detector face
  Absorbed,     ///< absorbed at a face (also a face where every branch has weight 0)
  Escaped,      ///< hit no face any more
  Cut,          ///< ended by max_depth or min_relative_power
  NotModelled,  ///< TIR inside a crystal, NoConvergence, no possible event
};

/// One event of a ray path: face, kind and diffraction order.
struct PathStep {
  std::uint32_t face = kNoFace;  ///< index into Scene::faces()
  model::EventKind kind = model::EventKind::Refract;
  int order = 0;  ///< diffraction order (ADR 0025); 0 without diffraction
  bool operator==(const PathStep&) const = default;
};

/// Branches of one source with the same sequence of events and the same end.
struct RayPath {
  std::vector<PathStep> steps;
  std::uint64_t rays = 0;  ///< branches that took this path
  double power = 0.0;      ///< W at the end, summed over the wavelengths of the source
  double share = 0.0;      ///< power / power of the source
};

/// The paths ending at one end point, by power descending.
struct PathEndpoint {
  PathEnd end = PathEnd::Escaped;
  /// Face of Detected, Absorbed and NotModelled at a face; kNoFace otherwise.
  std::uint32_t face = kNoFace;
  double power = 0.0;      ///< W of all paths to this end point, listed or not
  std::uint64_t rays = 0;  ///< branches of all paths to this end point
  /// Paths with share >= TraceOptions::report_min_share, by power descending.
  std::vector<RayPath> paths;
};

/// Ray path report of one source (Nachtrag PR 4, d): end points by power descending. Volume
/// absorption and event losses on the way belong to no path.
struct PathReport {
  /// W, power of the source over its wavelengths: the power of the start rays traced, so after
  /// a cancellation only that of the finished blocks; RayPath::share refers to it.
  double source = 0.0;
  std::vector<PathEndpoint> endpoints;
};

/// Power balance of one source at one wavelength (ADR 0057, point 8), W:
/// source = sum(detected) + sum(absorbed) + volume + escaped + cut + not_modelled up to rounding.
struct Balance {
  double source = 0.0;  ///< sum of the powers of the start rays
  /// Absorbed by each absorbing detector face (index of Scene::faces(); 0 for other faces).
  std::vector<double> detected;
  /// Absorbed at each face: Absorber interaction, back of a mirror, obscurations, stop screen, and
  /// the loss of an event w_in - W (coating absorption, power a polarizer rejects, diffraction
  /// orders without efficiency); the mantle faces separately.
  std::vector<double> absorbed;
  double volume = 0.0;        ///< absorbed in the media along the paths (Im n > 0)
  double escaped = 0.0;       ///< rays that hit no face any more (status Missed)
  double cut = 0.0;           ///< branches ended by max_depth or min_relative_power
  double not_modelled = 0.0;  ///< TIR inside a crystal, NoConvergence, impossible events
  /// Not part of the sum: power through each transparent detector, counted in both directions.
  std::vector<double> passed;

  /// Sum of all sinks (detected, absorbed, volume, escaped, cut, not_modelled), W.
  [[nodiscard]] double sinks() const noexcept;
};

/// Result of trace().
struct NonseqResult {
  /// Balance per source and wavelength: balance[source * wavelengths + wavelength], with
  /// `wavelengths` = CompiledSystem::wavelengths_um().size().
  std::vector<Balance> balance;
  /// Fields in the order of the detectors passed to trace().
  std::vector<DetectorField> detectors;
  /// Ray path report per source (TraceOptions::path_report; empty otherwise).
  std::vector<PathReport> paths;
  /// Hits on irradiance detectors in the order of the start rays (TraceOptions::record_hits).
  std::vector<DetectorHit> hits;
  std::uint64_t rays = 0;    ///< start rays traced (each wavelength counts)
  std::uint64_t events = 0;  ///< events executed (every branch counts its own)
  /// True if the run stopped after a cancellation request: the result holds the complete blocks
  /// of 4096 start rays that were finished (ADR 0057, point 10).
  bool cancelled = false;
};

/// Start rays per block of the parallel run; independent of the thread count, so that the sums
/// run in a fixed order and the results are bitwise the same for every thread count within a
/// build (ADR 0057, point 5; rule 7).
inline constexpr std::size_t kNonseqBlockSize = 4096;

/// Traces the start rays of `sources` through `scene` and counts them on `detectors`
/// (ADR 0057, points 4-8).
///
/// Start rays: source s with wavelengths L_s (all if empty) and M = source_ray_count() runs the
/// M rays of source_ray(source, s, i, seed) once per wavelength, wavelength-major, each with
/// power Phi_lambda / M (Nachtrag PR 4, f and g).
///
/// At each hit (Scene::closest_hit, self-hit excluded) the ray moves to the hit point (volume
/// absorption into Balance::volume) and the face decides:
/// - outside rays: Missed -> escaped, counted by every intensity detector; NoConvergence ->
///   not_modelled;
/// - region Absorbed (obscuration, stop screen), Absorber interaction, back of a mirror ->
///   absorbed at the face; an absorbing detector -> detected (with or without an irradiance
///   detector on it, Nachtrag PR 4, c); a transparent detector records Balance::passed and lets
///   the ray pass without an event; an irradiance detector counts both;
/// - otherwise the candidates of ADR 0057, point 4, each computed with trace::apply_event and the
///   media of the side the ray arrives from: Refract and Reflect between different media for
///   Fresnel, coatings and ideal_anti_reflection (only Reflect at TIR), a pass without event
///   between equal media; ideal_beam_splitter always Reflect and Transmit or Refract (ADR 0057,
///   Nachtrag PR 3, a); Reflect at the front of a mirror and at ideal_mirror; Transmit at
///   polarizers, retarders and ideal lenses; every listed diffraction order with every kind (order
///   0 only without efficiencies); Ordinary and Extraordinary into a crystal, Refract only inside
///   it (Tir there -> not_modelled).
/// The loss w_in - W of an event is absorbed at the face; if no candidate is possible
/// (EventImpossible), the ray is not_modelled. Monte Carlo draws one number per event from
/// Philox-4x32-10 with key = seed and counter (source 2^16 + wavelength, index low, index high,
/// 3 + event number) (ADR 0057, point 5, Nachtrag Review PR 1, c).
///
/// The start rays of all sources, in order, are split into blocks of kNonseqBlockSize; blocks
/// run in parallel (oneTBB) and are summed in block order. `control` cancels between blocks and
/// reports progress per block; a cancelled run returns with NonseqResult::cancelled.
/// @param scene built from `system` (Scene::build); the media indices of its faces refer to it
/// @throws SceneError with nonseq.detector_aperture (aperture not rectangular) or
///         nonseq.detector_shape (not a plane) for an irradiance detector
/// @throws std::invalid_argument for more than 65535 sources, a source that check_source()
///         rejects, an irradiance detector on a face that is no Detector, a detector with 0
///         cells or an intensity max_angle outside (0, pi], max_depth 0, min_relative_power < 0
///         or not finite, or report_min_share outside [0, 1]
/// @throws whatever the progress callback throws (RunMonitor)
[[nodiscard]] NonseqResult trace(const compile::CompiledSystem& system,
                                 const Scene& scene,
                                 std::span<const Source> sources,
                                 std::span<const Detector> detectors,
                                 const TraceOptions& options = {},
                                 const trace::RunControl& control = {});

/// As above, with explicit start rays: ray i of `sources[s]` has ray index i, power and
/// wavelength as given; its direction is normalised.
/// @throws std::invalid_argument as above, and for a wavelength index out of range, a power that
///         is negative or not finite, a direction that is zero or not finite or a position that
///         is not finite
[[nodiscard]] NonseqResult trace(const compile::CompiledSystem& system,
                                 const Scene& scene,
                                 std::span<const RaySource> sources,
                                 std::span<const Detector> detectors,
                                 const TraceOptions& options,
                                 const trace::RunControl& control = {});

/// A path of the report as an explicit sequential path (ADR 0027; ADR 0057, point 8), e.g. for a
/// stray-light ghost: one event per step at the model surface of its face, with kind and order,
/// and for a path that ends at a detector a Transmit event at the detector surface.
/// @param end  the end point the path belongs to
/// @return none if a step lies on a generated face (mantle or ring shoulder), which has no
///         surface in the model
/// @throws std::out_of_range for a face index outside Scene::faces()
[[nodiscard]] std::optional<model::Path> to_sequential_path(const compile::CompiledSystem& system,
                                                            const Scene& scene,
                                                            const PathEndpoint& end,
                                                            const RayPath& path,
                                                            std::string name);

/// Explicit start rays without detectors (PR 3).
[[nodiscard]] NonseqResult trace(const compile::CompiledSystem& system,
                                 const Scene& scene,
                                 std::span<const RaySource> sources,
                                 const TraceOptions& options = {},
                                 const trace::RunControl& control = {});

}  // namespace rtt::nonseq
