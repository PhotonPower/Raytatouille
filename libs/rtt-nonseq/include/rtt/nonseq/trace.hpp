#pragma once

/// @file trace.hpp
/// Trace core of the non-sequential tracer (ADR 0057, points 4, 5 and 8; PR 3): rays from start
/// rays through the Scene, branching at every hit by Monte Carlo or full splitting, with a power
/// balance per source and wavelength. Sources and detectors are plain C++ structures until the
/// format carries them (ADR 0057, point 12).
///
/// Conventions (docs/architecture.md, Konventionen): positions in mm and unit directions in
/// global coordinates, powers in W, wavelengths as indices into CompiledSystem::wavelengths_um().
/// The sources are unpolarized: a ray carries power = power of its start ray * weight, with the
/// weight of ADR 0021 starting at 1.

#include <cstdint>
#include <span>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/math/types.hpp"
#include "rtt/nonseq/scene.hpp"
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

/// One start ray (until PR 4, sources are given as their rays).
struct StartRay {
  math::Vec3 position = math::Vec3::Zero();    ///< global, mm
  math::Vec3 direction = math::Vec3::UnitZ();  ///< global, normalised by trace()
  double power = 0.0;                          ///< W, finite and >= 0
  std::uint16_t wavelength = 0;                ///< index into CompiledSystem::wavelengths_um()
};

/// A source as a list of start rays; its index in the list passed to trace() is the source
/// number of the random counter (ADR 0057, point 5) and of the balance.
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
  std::uint64_t rays = 0;    ///< start rays traced
  std::uint64_t events = 0;  ///< events executed (every branch counts its own)
  /// True if the run stopped after a cancellation request: the result holds the complete blocks
  /// of 4096 start rays that were finished (ADR 0057, point 10).
  bool cancelled = false;
};

/// Start rays per block of the parallel run; independent of the thread count, so that the sums
/// run in a fixed order and the results are bitwise the same for every thread count within a
/// build (ADR 0057, point 5; rule 7).
inline constexpr std::size_t kNonseqBlockSize = 4096;

/// Traces the start rays of `sources` through `scene` (ADR 0057, points 4, 5 and 8).
///
/// At each hit (Scene::closest_hit, self-hit excluded) the ray moves to the hit point (volume
/// absorption into Balance::volume) and the face decides:
/// - outside rays: Missed -> escaped; NoConvergence -> not_modelled;
/// - region Absorbed (obscuration, stop screen), Absorber interaction, back of a mirror ->
///   absorbed at the face; an absorbing detector -> detected; a transparent detector records
///   Balance::passed and lets the ray pass without an event;
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
/// @throws std::invalid_argument for more than 65535 sources, a wavelength index out of range, a
///         power that is negative or not finite, a direction that is zero or not finite, a
///         position that is not finite, max_depth 0, or min_relative_power < 0 or not finite
/// @throws whatever the progress callback throws (RunMonitor)
[[nodiscard]] NonseqResult trace(const compile::CompiledSystem& system,
                                 const Scene& scene,
                                 std::span<const RaySource> sources,
                                 const TraceOptions& options = {},
                                 const trace::RunControl& control = {});

}  // namespace rtt::nonseq
