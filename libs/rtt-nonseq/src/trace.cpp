#include "rtt/nonseq/trace.hpp"

#include <oneapi/tbb/blocked_range.h>
#include <oneapi/tbb/parallel_for.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <variant>
#include <vector>

#include "rtt/coating/transfer_matrix.hpp"
#include "rtt/diagnostics/codes.hpp"
#include "rtt/math/philox.hpp"
#include "rtt/model/path.hpp"
#include "rtt/model/surface.hpp"
#include "rtt/trace/apply_event.hpp"

namespace rtt::nonseq {

double Balance::sinks() const noexcept {
  double s = volume + escaped + cut + not_modelled;
  for (const double d : detected) s += d;
  for (const double a : absorbed) s += a;
  return s;
}

namespace {

using compile::CompiledSurface;
using compile::CrystalMode;
using model::EventKind;
using trace::RayState;

/// Crystal data of `medium` at wavelength `wl` (ADR 0026, point 3), as rtt-trace fills it for a
/// compiled event; none for an isotropic medium.
std::optional<trace::CrystalSide> crystal_side(const compile::CompiledMedium& medium,
                                               std::uint16_t wl) noexcept {
  if (!medium.is_crystal()) return std::nullopt;
  return trace::CrystalSide{medium.index[wl].real(), medium.index_extraordinary[wl],
                            medium.optic_axis.value_or(math::Vec3::UnitZ())};
}

/// One candidate branch of an event.
struct Candidate {
  EventKind kind = EventKind::Refract;
  int order = 0;
};

/// A ray in flight: its state, the mode it entered a crystal with, its depth and the number of
/// Monte Carlo draws made for it.
struct Branch {
  RayState ray;
  CrystalMode mode = CrystalMode::None;
  std::uint32_t last_face = kNoFace;
  std::uint32_t depth = 0;
  std::uint32_t draws = 0;
  std::vector<PathStep> steps;  ///< events so far, only with TraceOptions::path_report
};

/// Read-only data of a run, shared by all blocks.
struct Run {
  const compile::CompiledSystem& system;
  const Scene& scene;
  const TraceOptions& options;
  std::size_t wavelengths = 0;
  /// Coating stacks seen from the substrate (ADR 0019), per coating and wavelength, prepared
  /// once so that the loop does not allocate.
  std::vector<std::vector<std::vector<coating::Layer<double>>>> reversed;
  /// Per face with ideal_anti_reflection or ideal_beam_splitter: a copy of its surface with
  /// Fresnel interaction, for the reflection at TIR (ADR 0057, Nachtrag PR 3, b); none otherwise.
  std::vector<std::optional<CompiledSurface>> fresnel_twin;
  /// Irradiance detectors (index into the detectors of trace()) per face.
  std::vector<std::vector<std::uint32_t>> irradiance_of_face;
  /// Intensity detectors (index into the detectors of trace()).
  std::vector<std::uint32_t> intensity;
  std::span<const Detector> detectors;
  /// Half widths (x, y) of the rectangle of each irradiance detector, mm (0 for others).
  std::vector<std::array<double, 2>> half_width;
};

/// One count of a detector cell; cell = (wavelength * rows + row) * cols + col.
struct CellCount {
  std::uint32_t detector = 0;
  std::uint64_t cell = 0;
  double power = 0.0;
  double weight = 0.0;
};

/// FNV-1a, 64 bit (Fowler, Noll, Vo; offset basis 14695981039346656037, prime 1099511628211),
/// over the bytes of `value`, least significant first (independent of the platform).
template <typename T>
void fnv1a(std::uint64_t& h, T value) noexcept {
  auto v = static_cast<std::uint64_t>(value);
  for (std::size_t k = 0; k < sizeof(T); ++k) {
    h ^= (v >> (8 * k)) & 0xffU;
    h *= 1099511628211ULL;
  }
}

/// Paths of one source: equal (end, face, steps) are merged; the FNV-1a hash finds candidates,
/// which are compared in full (ADR 0057, point 8).
class PathTable {
 public:
  struct Entry {
    PathEnd end = PathEnd::Escaped;
    std::uint32_t face = kNoFace;
    std::vector<PathStep> steps;
    std::uint64_t rays = 0;
    double power = 0.0;
  };

  void add(PathEnd end,
           std::uint32_t face,
           const std::vector<PathStep>& steps,
           double power,
           std::uint64_t rays) {
    std::uint64_t h = 14695981039346656037ULL;
    fnv1a(h, static_cast<std::uint8_t>(end));
    fnv1a(h, face);
    for (const PathStep& st : steps) {
      fnv1a(h, st.face);
      fnv1a(h, static_cast<std::uint8_t>(st.kind));
      fnv1a(h, static_cast<std::uint32_t>(st.order));
    }
    std::vector<std::uint32_t>& slots = index_[h];
    for (const std::uint32_t k : slots) {
      Entry& e = entries_[k];
      if (e.end == end && e.face == face && e.steps == steps) {
        e.rays += rays;
        e.power += power;
        return;
      }
    }
    slots.push_back(static_cast<std::uint32_t>(entries_.size()));
    entries_.push_back({end, face, steps, rays, power});
  }

  [[nodiscard]] const std::vector<Entry>& entries() const noexcept { return entries_; }

 private:
  std::vector<Entry> entries_;
  std::unordered_map<std::uint64_t, std::vector<std::uint32_t>> index_;
};

/// Sums of one block.
struct Partial {
  std::vector<Balance> balance;
  std::vector<CellCount> cells;
  std::vector<DetectorHit> hits;
  std::vector<PathTable> paths;  ///< per source
  std::uint64_t rays = 0;
  std::uint64_t events = 0;
  bool done = false;
};

/// Where a start ray's results go.
struct Sink {
  Balance& balance;
  Partial& partial;
  std::uint32_t source = 0;
  std::uint64_t index = 0;
};

Balance empty_balance(std::size_t faces) {
  Balance b;
  b.detected.assign(faces, 0.0);
  b.absorbed.assign(faces, 0.0);
  b.passed.assign(faces, 0.0);
  return b;
}

/// The candidates of ADR 0057, point 4, for a hit on `face` that is not absorbing. `inside` is
/// the crystal of the medium before (exit or pass inside), `entering` the crystal beyond.
void candidates(const Face& face,
                const CompiledSurface& s,
                bool same_medium,
                bool inside,
                bool entering,
                std::vector<Candidate>& out) {
  out.clear();
  std::array<EventKind, 2> kinds{};
  std::size_t n = 0;
  if (face.kind == FaceKind::Mirror) {
    kinds[n++] = EventKind::Reflect;
  } else if (inside) {
    // Inside a crystal only the exit (ADR 0057, Nachtrag Review PR 1, b); a cemented surface
    // between segments of the same crystal is passed.
    kinds[n++] = same_medium ? EventKind::Transmit : EventKind::Refract;
  } else if (entering) {
    kinds[n++] = EventKind::Ordinary;
    kinds[n++] = EventKind::Extraordinary;
  } else {
    std::visit(
        [&](const auto& interaction) {
          using I = std::remove_cvref_t<decltype(interaction)>;
          if constexpr (std::is_same_v<I, model::IdealMirror>) {
            kinds[n++] = EventKind::Reflect;
          } else if constexpr (std::is_same_v<I, model::IdealBeamSplitter>) {
            kinds[n++] =
                same_medium ? EventKind::Transmit : EventKind::Refract;  // Nachtrag PR 3, a
            kinds[n++] = EventKind::Reflect;
          } else if constexpr (std::is_same_v<I, model::IdealPolarizer> ||
                               std::is_same_v<I, model::IdealRetarder> ||
                               std::is_same_v<I, model::IdealLens> ||
                               std::is_same_v<I, model::IdealCylinderLens>) {
            kinds[n++] = EventKind::Transmit;
          } else {
            // Fresnel, CoatingRef, IdealAntiReflection (Absorber is handled before).
            if (same_medium) {
              kinds[n++] = EventKind::Transmit;
            } else {
              kinds[n++] = EventKind::Refract;
              kinds[n++] = EventKind::Reflect;
            }
          }
        },
        s.interaction);
  }
  // Diffraction orders (ADR 0057, point 4; #450, point 7): every listed order with every kind;
  // without efficiencies only order 0.
  if (!s.phase_functions.empty() && s.diffraction_efficiency && !inside && !entering) {
    for (std::size_t k = 0; k < n; ++k) {
      for (const model::DiffractionEfficiency& e : *s.diffraction_efficiency) {
        out.push_back({kinds[k], e.order});
      }
    }
    return;
  }
  for (std::size_t k = 0; k < n; ++k) out.push_back({kinds[k], 0});
}

/// Counts a hit on irradiance detector `d` (ADR 0057, point 7): pixel of the local point in the
/// rectangle (points on the rim within kApertureTolerance go to the edge pixel).
void count_irradiance(const Run& run,
                      std::uint32_t d,
                      const SceneHit& hit,
                      const math::Vec3& direction,
                      std::uint16_t wl,
                      double power,
                      double weight,
                      Sink& sink) {
  const auto& det = std::get<IrradianceDetector>(run.detectors[d]);
  const double hx = run.half_width[d][0];
  const double hy = run.half_width[d][1];
  const auto bin = [](double v, double half, std::uint32_t n) {
    const double f = std::floor((v + half) / (2.0 * half) * n);
    return static_cast<std::uint32_t>(std::clamp(f, 0.0, static_cast<double>(n - 1)));
  };
  const std::uint32_t col = bin(hit.local.point.x(), hx, det.nx);
  const std::uint32_t row = bin(hit.local.point.y(), hy, det.ny);
  const std::uint64_t cell = (static_cast<std::uint64_t>(wl) * det.ny + row) * det.nx + col;
  sink.partial.cells.push_back({d, cell, power, weight});
  if (run.options.record_hits) {
    sink.partial.hits.push_back({d, sink.source, sink.index, wl, hit.point, direction, power});
  }
}

/// Counts an escaped ray on every intensity detector (ADR 0057, point 7; Nachtrag PR 4, h).
void count_escaped(
    const Run& run, const RayState& ray, std::uint16_t wl, double power, Partial& partial) {
  for (const std::uint32_t d : run.intensity) {
    const auto& det = std::get<IntensityDetector>(run.detectors[d]);
    const math::Vec3 v = det.pose.rotation().transpose() * ray.dir;
    const double theta = std::atan2(std::hypot(v.x(), v.y()), v.z());
    if (theta > det.max_angle) continue;
    double phi = std::atan2(v.y(), v.x());
    if (phi < 0.0) phi += 2.0 * std::numbers::pi;
    const auto row = static_cast<std::uint32_t>(
        std::min(std::floor(theta / det.max_angle * det.n_theta), det.n_theta - 1.0));
    const auto col = static_cast<std::uint32_t>(
        std::min(std::floor(phi / (2.0 * std::numbers::pi) * det.n_phi), det.n_phi - 1.0));
    const std::uint64_t cell =
        (static_cast<std::uint64_t>(wl) * det.n_theta + row) * det.n_phi + col;
    partial.cells.push_back({d, cell, power, ray.weight});
  }
}

/// Traces one start ray with all its branches into `b`.
void trace_ray(const Run& run, const StartRay& start, Sink& sink) {
  Balance& b = sink.balance;
  std::uint64_t& events = sink.partial.events;
  const std::uint32_t source = sink.source;
  const std::uint64_t index = sink.index;
  const bool report = run.options.path_report;
  const auto& faces = run.scene.faces();
  const auto& media = run.system.media();
  const std::uint16_t wl = start.wavelength;
  const double p0 = start.power;
  const TraceOptions& o = run.options;
  const math::PhiloxKey key{static_cast<std::uint32_t>(o.seed),
                            static_cast<std::uint32_t>(o.seed >> 32)};
  const auto draw = [&](std::uint32_t j) {
    // Counter of ADR 0057, point 5, with the wavelength in the first word (Nachtrag Review PR 1,
    // c); draws 0-2 are reserved for the start ray (PR 4).
    const math::PhiloxCounter counter{(source << 16) | wl, static_cast<std::uint32_t>(index),
                                      static_cast<std::uint32_t>(index >> 32), 3 + j};
    const auto w = math::philox4x32_10(counter, key);
    return math::uniform53(w[0], w[1]);
  };

  std::vector<Branch> stack;
  std::vector<Candidate> cands;
  std::vector<RayState> results;
  Branch first;
  first.ray.pos = start.position;
  first.ray.dir = start.direction.normalized();
  first.ray.wave = first.ray.dir;
  first.ray.last_surface = trace::kNoSurface;
  stack.push_back(first);

  while (!stack.empty()) {
    Branch br = std::move(stack.back());
    stack.pop_back();
    const auto finish = [&](PathEnd end, std::uint32_t face, double power) {
      if (report) sink.partial.paths[source].add(end, face, br.steps, power, 1);
    };
    for (;;) {
      const SceneHit hit = run.scene.closest_hit(br.ray, br.last_face);
      if (hit.status == geom::HitStatus::Missed) {
        const double power = p0 * br.ray.weight;
        b.escaped += power;
        count_escaped(run, br.ray, wl, power, sink.partial);
        finish(PathEnd::Escaped, kNoFace, power);
        break;
      }
      if (hit.status != geom::HitStatus::Hit) {
        b.not_modelled += p0 * br.ray.weight;
        finish(PathEnd::NotModelled, kNoFace, p0 * br.ray.weight);
        break;
      }
      const Face& face = faces[hit.face];
      const CompiledSurface& s = run.scene.surfaces()[face.surface];
      // Media from the side the ray arrives from (ADR 0057, point 2).
      const std::uint32_t m_before = hit.medium_before;
      const std::uint32_t m_beyond = hit.medium_beyond;
      trace::EventMedia base;
      base.before = media[m_before].index[wl];
      base.beyond = media[m_beyond].index[wl];
      base.after = base.before;
      base.wavelength_um = run.system.wavelengths_um()[wl];
      base.crystal_before = crystal_side(media[m_before], wl);
      base.crystal_mode = base.crystal_before ? br.mode : CrystalMode::None;

      const RayState moved = trace::move_to_hit(br.ray, s, hit.local, hit.face, base);
      b.volume += p0 * (br.ray.weight - moved.weight);
      const double w_in = moved.weight;
      const auto absorb = [&] {
        b.absorbed[hit.face] += p0 * w_in;
        finish(PathEnd::Absorbed, hit.face, p0 * w_in);
      };

      if (hit.region == HitRegion::Absorbed) {
        absorb();
        break;
      }
      if (face.kind == FaceKind::Detector) {
        // Irradiance detectors count from both sides (ADR 0057, point 7); a Detector face
        // without one counts only in the balance (Nachtrag PR 4, c).
        for (const std::uint32_t d : run.irradiance_of_face[hit.face]) {
          count_irradiance(run, d, hit, br.ray.dir, wl, p0 * w_in, w_in, sink);
        }
        if (!face.transparent) {
          b.detected[hit.face] += p0 * w_in;
          finish(PathEnd::Detected, hit.face, p0 * w_in);
          break;
        }
        b.passed[hit.face] += p0 * w_in;
        br.ray = moved;
        br.last_face = hit.face;
        continue;
      }
      // The back of a mirror absorbs: the front is the side a ray with d . z_local > 0 comes
      // from, the -z side (ADR 0057, point 2).
      if ((face.kind == FaceKind::Mirror && hit.from_positive) ||
          std::holds_alternative<model::Absorber>(s.interaction)) {
        absorb();
        break;
      }
      const bool same = m_before == m_beyond;
      const bool inside = media[m_before].is_crystal();
      const bool entering = !inside && !same && media[m_beyond].is_crystal();
      candidates(face, s, same, inside, entering, cands);
      // Fresnel, coatings and ideal_anti_reflection between equal media: a pass without event
      // (ADR 0057, point 4); a diffraction order still makes an event.
      const bool interface = std::holds_alternative<model::Fresnel>(s.interaction) ||
                             std::holds_alternative<model::CoatingRef>(s.interaction) ||
                             std::holds_alternative<model::IdealAntiReflection>(s.interaction);
      const bool orders = !s.phase_functions.empty() && s.diffraction_efficiency.has_value();
      if (interface && same && !inside && !orders && face.kind != FaceKind::Mirror &&
          cands.size() == 1 && cands[0].kind == EventKind::Transmit && cands[0].order == 0) {
        br.ray = moved;
        br.last_face = hit.face;
        continue;
      }

      results.clear();
      double total = 0.0;
      bool any_possible = false;
      bool tir = false;
      for (const Candidate& c : cands) {
        trace::EventMedia m = base;
        const bool crossing = c.kind == EventKind::Refract || c.kind == EventKind::Ordinary ||
                              c.kind == EventKind::Extraordinary;
        m.after = crossing ? base.beyond : base.before;
        // A mirror reflects against its substrate, or as an ideal conductor without one.
        if (face.kind == FaceKind::Mirror) {
          m.beyond = face.substrate ? media[*face.substrate].index[wl] : base.before;
        }
        // Into a crystal, and the pass between two segments of the same crystal (both sides).
        if (crossing || (inside && same)) m.crystal_after = crystal_side(media[m_beyond], wl);
        if (s.coating) {
          const bool from_inside = m_before == s.coating->substrate_medium;
          const auto& cs = run.system.coatings()[s.coating->coating];
          m.layers =
              from_inside
                  ? std::span<const coating::Layer<double>>(run.reversed[s.coating->coating][wl])
                  : std::span<const coating::Layer<double>>(cs.layers[wl]);
        }
        RayState r = trace::apply_event(br.ray, s, hit.local, hit.face, c.kind, c.order, m);
        if (!inside && c.kind == EventKind::Refract && r.status == trace::RayStatus::Tir)
          tir = true;
        if (r.status != trace::RayStatus::EventImpossible) any_possible = true;
        if (inside && r.status == trace::RayStatus::Tir) {
          // TIR inside a crystal is not modelled (ADR 0057, Nachtrag Review PR 1, b).
          results.clear();
          total = 0.0;
          any_possible = false;
          break;
        }
        if (r.status != trace::RayStatus::Alive || !(r.weight > 0.0)) r.weight = 0.0;
        total += r.weight;
        results.push_back(r);
      }
      // TIR at an ideal anti-reflection or beam-splitter surface reflects like a Fresnel surface,
      // fully (ADR 0057, Nachtrag PR 3, b): the idealisation holds only below the critical angle.
      if (tir && run.fresnel_twin[hit.face] && results.size() == cands.size()) {
        total = 0.0;
        for (std::size_t i = 0; i < cands.size(); ++i) {
          if (cands[i].kind == EventKind::Reflect) {
            RayState r = trace::apply_event(br.ray, *run.fresnel_twin[hit.face], hit.local,
                                            hit.face, EventKind::Reflect, cands[i].order, base);
            if (r.status != trace::RayStatus::Alive || !(r.weight > 0.0)) r.weight = 0.0;
            results[i] = r;
          }
          total += results[i].weight;
        }
      }
      ++events;
      if (!any_possible) {
        b.not_modelled += p0 * w_in;
        finish(PathEnd::NotModelled, hit.face, p0 * w_in);
        break;
      }
      b.absorbed[hit.face] += p0 * (w_in - total);  // loss of the event
      if (!(total > 0.0)) {
        finish(PathEnd::Absorbed, hit.face, p0 * w_in);
        break;
      }

      const std::uint32_t depth = br.depth + 1;
      const auto next_mode = [&](EventKind k) {
        if (k == EventKind::Ordinary) return CrystalMode::Ordinary;
        if (k == EventKind::Extraordinary) return CrystalMode::Extraordinary;
        if (inside && k == EventKind::Refract) return CrystalMode::None;
        return br.mode;
      };
      const auto cut = [&](double weight) {
        return depth > o.max_depth || weight < o.min_relative_power;
      };
      if (o.splitting == Splitting::MonteCarlo) {
        // Branch i with p_i = w_i / W; the ray keeps P_i and gets the weight W (Nachtrag Review
        // PR 1, a; PBRT 4ed, Eq. (2.16) with one sample, §9.5).
        const double u = draw(br.draws) * total;
        ++br.draws;
        std::size_t pick = results.size();
        double cum = 0.0;
        for (std::size_t i = 0; i < results.size(); ++i) {
          if (!(results[i].weight > 0.0)) continue;
          cum += results[i].weight;
          pick = i;
          if (u < cum) break;
        }
        br.ray = results[pick];
        br.ray.weight = total;
        br.mode = next_mode(cands[pick].kind);
        br.last_face = hit.face;
        br.depth = depth;
        if (report) br.steps.push_back({hit.face, cands[pick].kind, cands[pick].order});
        if (cut(total)) {
          b.cut += p0 * total;
          finish(PathEnd::Cut, kNoFace, p0 * total);
          break;
        }
        continue;
      }
      // Full splitting: push the branches in reverse so that the first candidate runs first.
      for (std::size_t i = results.size(); i-- > 0;) {
        if (!(results[i].weight > 0.0)) continue;
        Branch nb;
        nb.ray = results[i];
        nb.mode = next_mode(cands[i].kind);
        nb.last_face = hit.face;
        nb.depth = depth;
        if (report) {
          nb.steps.reserve(br.steps.size() + 1);
          nb.steps = br.steps;
          nb.steps.push_back({hit.face, cands[i].kind, cands[i].order});
        }
        if (cut(results[i].weight)) {
          b.cut += p0 * results[i].weight;
          if (report)
            sink.partial.paths[source].add(PathEnd::Cut, kNoFace, nb.steps, p0 * results[i].weight,
                                           1);
          continue;
        }
        stack.push_back(std::move(nb));
      }
      break;
    }
  }
}

/// One input source: generated (Source) or an explicit list (RaySource).
struct Input {
  const Source* generated = nullptr;
  const RaySource* list = nullptr;
  std::vector<std::uint16_t> wavelengths;  ///< generated: the wavelengths in run order
  std::vector<double> ray_power;           ///< generated: Phi_lambda / M per wavelength slot
  std::uint64_t m = 0;                     ///< generated: start rays per wavelength
  std::uint64_t count = 0;                 ///< start rays in total
};

void check_options(const TraceOptions& options) {
  if (options.max_depth == 0) throw std::invalid_argument("nonseq::trace: max_depth must be > 0");
  if (!(std::isfinite(options.min_relative_power) && options.min_relative_power >= 0.0)) {
    throw std::invalid_argument("nonseq::trace: min_relative_power must be finite and >= 0");
  }
  if (!(options.report_min_share >= 0.0 && options.report_min_share <= 1.0)) {
    throw std::invalid_argument("nonseq::trace: report_min_share must be in [0, 1]");
  }
}

/// Checks the detectors (ADR 0057, point 7; Nachtrag PR 4, b) and returns the empty fields.
std::vector<DetectorField> prepare_detectors(const compile::CompiledSystem& system,
                                             const Scene& scene,
                                             std::span<const Detector> detectors,
                                             std::size_t wls) {
  std::vector<model::Diagnostic> errors;
  const auto report = [&](diagnostics::DiagnosticCode code, std::string location,
                          std::string message) {
    errors.push_back(
        {code.severity(), std::move(location), std::move(message), std::string(code.str())});
  };
  std::vector<DetectorField> fields;
  for (const Detector& detector : detectors) {
    DetectorField f;
    if (const auto* irr = std::get_if<IrradianceDetector>(&detector)) {
      if (irr->face >= scene.faces().size() ||
          scene.faces()[irr->face].kind != FaceKind::Detector) {
        throw std::invalid_argument("nonseq::trace: irradiance detector on a face of no Detector");
      }
      if (irr->nx == 0 || irr->ny == 0)
        throw std::invalid_argument("nonseq::trace: detector with 0 pixels");
      const CompiledSurface& s = scene.surfaces()[scene.faces()[irr->face].surface];
      const std::string& loc = system.surfaces()[scene.faces()[irr->face].system_surface].location;
      const auto* rect =
          s.aperture ? std::get_if<model::RectangularAperture>(&*s.aperture) : nullptr;
      if (rect == nullptr) {
        report("nonseq.detector_aperture", loc + "/aperture",
               "an irradiance detector needs a rectangular aperture");
      }
      if (!std::holds_alternative<geom::Plane<double>>(s.shape)) {
        report("nonseq.detector_shape", loc + "/shape",
               "an irradiance detector needs a plane surface without terms");
      }
      f.rows = irr->ny;
      f.cols = irr->nx;
      if (rect != nullptr) {
        // A_pixel = (2 h_x / n_x)(2 h_y / n_y).
        f.cell_size = {(2.0 * rect->half_width_x.value / irr->nx) *
                       (2.0 * rect->half_width_y.value / irr->ny)};
      }
    } else {
      const auto& in = std::get<IntensityDetector>(detector);
      if (!(in.max_angle > 0.0 && in.max_angle <= std::numbers::pi))
        throw std::invalid_argument("nonseq::trace: intensity max_angle must be in (0, pi]");
      if (in.n_theta == 0 || in.n_phi == 0)
        throw std::invalid_argument("nonseq::trace: detector with 0 cells");
      if (!in.pose.rotation().allFinite())
        throw std::invalid_argument("nonseq::trace: intensity detector pose not finite");
      f.rows = in.n_theta;
      f.cols = in.n_phi;
      // Delta Omega = Delta phi (cos theta_1 - cos theta_2) (NBS TN 910-1, Eq. (A2-9)), as
      // 2 Delta phi sin((theta_2 - theta_1)/2) sin((theta_1 + theta_2)/2) against
      // cancellation (ADR 0057, point 7).
      const double dt = in.max_angle / in.n_theta;
      const double dp = 2.0 * std::numbers::pi / in.n_phi;
      for (std::uint32_t r = 0; r < in.n_theta; ++r) {
        const double t1 = r * dt;
        const double t2 = (r + 1) * dt;
        f.cell_size.push_back(2.0 * dp * std::sin(0.5 * (t2 - t1)) * std::sin(0.5 * (t1 + t2)));
      }
    }
    const std::size_t cells = wls * f.rows * f.cols;
    f.value.assign(cells, 0.0);
    f.hits.assign(cells, 0);
    f.weights.assign(cells, 0.0);
    fields.push_back(std::move(f));
  }
  if (!errors.empty()) throw SceneError(std::move(errors));
  return fields;
}

/// The report of one source from its merged table (ADR 0057, point 8; Nachtrag PR 4, d).
PathReport make_report(const PathTable& table, double source, double min_share) {
  PathReport rep;
  rep.source = source;
  std::vector<const PathTable::Entry*> entries;
  for (const PathTable::Entry& e : table.entries()) entries.push_back(&e);
  // By end point, then power descending, then the steps (a total order: deterministic).
  const auto step_key = [](const PathStep& st) {
    return std::tuple(st.face, static_cast<int>(st.kind), st.order);
  };
  std::sort(entries.begin(), entries.end(), [&](const auto* a, const auto* b) {
    if (a->end != b->end) return a->end < b->end;
    if (a->face != b->face) return a->face < b->face;
    if (a->power != b->power) return a->power > b->power;
    return std::lexicographical_compare(
        a->steps.begin(), a->steps.end(), b->steps.begin(), b->steps.end(),
        [&](const PathStep& x, const PathStep& y) { return step_key(x) < step_key(y); });
  });
  for (const PathTable::Entry* e : entries) {
    if (rep.endpoints.empty() || rep.endpoints.back().end != e->end ||
        rep.endpoints.back().face != e->face) {
      rep.endpoints.push_back({e->end, e->face, 0.0, 0, {}});
    }
    PathEndpoint& ep = rep.endpoints.back();
    ep.power += e->power;
    ep.rays += e->rays;
    const double share = source > 0.0 ? e->power / source : 0.0;
    if (share >= min_share) ep.paths.push_back({e->steps, e->rays, e->power, share});
  }
  std::stable_sort(rep.endpoints.begin(), rep.endpoints.end(),
                   [](const PathEndpoint& a, const PathEndpoint& b) { return a.power > b.power; });
  return rep;
}

/// Blocks summed per wave: bounds the memory of the block results; the sums run in block order
/// for every wave size, so it does not change the result.
constexpr std::size_t kWave = 64;

NonseqResult run_trace(const compile::CompiledSystem& system,
                       const Scene& scene,
                       std::vector<Input> inputs,
                       std::span<const Detector> detectors,
                       const TraceOptions& options,
                       const trace::RunControl& control) {
  const std::size_t wls = system.wavelengths_um().size();
  const std::size_t faces = scene.faces().size();
  NonseqResult result;
  result.detectors = prepare_detectors(system, scene, detectors, wls);

  std::vector<std::uint64_t> offsets{0};
  for (const Input& in : inputs) offsets.push_back(offsets.back() + in.count);
  const std::uint64_t total = offsets.back();

  Run run{system, scene, options, wls, {}, {}, {}, {}, detectors, {}};
  run.reversed.resize(system.coatings().size());
  for (std::size_t c = 0; c < system.coatings().size(); ++c) {
    for (const auto& layers : system.coatings()[c].layers) {
      run.reversed[c].emplace_back(layers.rbegin(), layers.rend());
    }
  }
  run.fresnel_twin.resize(faces);
  for (std::size_t f = 0; f < faces; ++f) {
    const CompiledSurface& s = scene.surfaces()[scene.faces()[f].surface];
    if (std::holds_alternative<model::IdealAntiReflection>(s.interaction) ||
        std::holds_alternative<model::IdealBeamSplitter>(s.interaction)) {
      CompiledSurface twin = s;
      twin.interaction = model::Fresnel{};
      run.fresnel_twin[f] = std::move(twin);
    }
  }
  run.irradiance_of_face.resize(faces);
  run.half_width.assign(detectors.size(), {0.0, 0.0});
  for (std::uint32_t d = 0; d < detectors.size(); ++d) {
    if (const auto* irr = std::get_if<IrradianceDetector>(&detectors[d])) {
      run.irradiance_of_face[irr->face].push_back(d);
      // prepare_detectors() checked the rectangular aperture.
      const CompiledSurface& s = scene.surfaces()[scene.faces()[irr->face].surface];
      if (s.aperture) {
        const auto& aperture = *s.aperture;
        const auto& rect = std::get<model::RectangularAperture>(aperture);
        run.half_width[d] = {rect.half_width_x.value, rect.half_width_y.value};
      }
    } else {
      run.intensity.push_back(d);
    }
  }

  const std::size_t blocks =
      static_cast<std::size_t>((total + kNonseqBlockSize - 1) / kNonseqBlockSize);
  trace::RunMonitor monitor(control, static_cast<std::size_t>(total), "nonsequential trace");
  std::vector<Partial> partials;
  const auto run_block = [&](std::size_t blk, Partial& p) {
    if (monitor.stop()) return;
    p.balance.assign(inputs.size() * wls, empty_balance(faces));
    if (options.path_report) p.paths.assign(inputs.size(), PathTable{});
    const std::uint64_t begin = blk * kNonseqBlockSize;
    const std::uint64_t end = std::min<std::uint64_t>(begin + kNonseqBlockSize, total);
    auto src = static_cast<std::size_t>(std::upper_bound(offsets.begin(), offsets.end(), begin) -
                                        offsets.begin() - 1);
    for (std::uint64_t g = begin; g < end; ++g) {
      while (g >= offsets[src + 1]) ++src;
      const Input& in = inputs[src];
      const std::uint64_t local = g - offsets[src];
      StartRay start;
      std::uint64_t index = local;
      if (in.generated != nullptr) {
        // Wavelength-major: every wavelength runs the same M start rays (Nachtrag PR 4, g).
        const auto slot = static_cast<std::size_t>(local / in.m);
        index = local % in.m;
        const SourceRay sr =
            source_ray(*in.generated, static_cast<std::uint32_t>(src), index, options.seed);
        start = {sr.position, sr.direction, in.ray_power[slot], in.wavelengths[slot]};
      } else {
        start = in.list->rays[static_cast<std::size_t>(local)];
      }
      Balance& bal = p.balance[src * wls + start.wavelength];
      bal.source += start.power;
      Sink sink{bal, p, static_cast<std::uint32_t>(src), index};
      trace_ray(run, start, sink);
      ++p.rays;
    }
    p.done = true;
    monitor.add(static_cast<std::size_t>(end - begin));
  };

  result.balance.assign(inputs.size() * wls, empty_balance(faces));
  std::vector<PathTable> tables(options.path_report ? inputs.size() : 0);
  for (std::size_t wave = 0; wave < blocks; wave += kWave) {
    const std::size_t n = std::min(kWave, blocks - wave);
    partials.assign(n, Partial{});
    oneapi::tbb::parallel_for(oneapi::tbb::blocked_range<std::size_t>(0, n, 1),
                              [&](const oneapi::tbb::blocked_range<std::size_t>& range) {
                                for (std::size_t k = range.begin(); k != range.end(); ++k) {
                                  run_block(wave + k, partials[k]);
                                }
                              });
    // Sum in block order (rule 7).
    for (const Partial& p : partials) {
      if (!p.done) continue;
      result.rays += p.rays;
      result.events += p.events;
      for (std::size_t k = 0; k < p.balance.size(); ++k) {
        Balance& r = result.balance[k];
        const Balance& q = p.balance[k];
        r.source += q.source;
        r.volume += q.volume;
        r.escaped += q.escaped;
        r.cut += q.cut;
        r.not_modelled += q.not_modelled;
        for (std::size_t f = 0; f < faces; ++f) {
          r.detected[f] += q.detected[f];
          r.absorbed[f] += q.absorbed[f];
          r.passed[f] += q.passed[f];
        }
      }
      for (const CellCount& c : p.cells) {
        DetectorField& f = result.detectors[c.detector];
        f.value[c.cell] += c.power;
        f.hits[c.cell] += 1;
        f.weights[c.cell] += c.weight;
      }
      result.hits.insert(result.hits.end(), p.hits.begin(), p.hits.end());
      for (std::size_t s = 0; s < p.paths.size(); ++s) {
        for (const PathTable::Entry& e : p.paths[s].entries()) {
          tables[s].add(e.end, e.face, e.steps, e.power, e.rays);
        }
      }
    }
  }
  try {
    monitor.finish();
  } catch (const trace::Cancelled&) {
    result.cancelled = true;
  }
  // Power sums to E = sum(Phi) / A_pixel and I = sum(Phi) / Delta Omega.
  for (std::size_t d = 0; d < detectors.size(); ++d) {
    DetectorField& f = result.detectors[d];
    const bool irradiance = std::holds_alternative<IrradianceDetector>(detectors[d]);
    for (std::size_t c = 0; c < f.value.size(); ++c) {
      f.value[c] /= irradiance ? f.cell_size[0] : f.cell_size[(c / f.cols) % f.rows];
    }
  }
  for (std::size_t s = 0; s < tables.size(); ++s) {
    double source = 0.0;
    for (std::size_t w = 0; w < wls; ++w) source += result.balance[s * wls + w].source;
    result.paths.push_back(make_report(tables[s], source, options.report_min_share));
  }
  return result;
}

}  // namespace

NonseqResult trace(const compile::CompiledSystem& system,
                   const Scene& scene,
                   std::span<const Source> sources,
                   std::span<const Detector> detectors,
                   const TraceOptions& options,
                   const trace::RunControl& control) {
  const std::size_t wls = system.wavelengths_um().size();
  if (sources.size() > 0xffff)
    throw std::invalid_argument("nonseq::trace: more than 65535 sources");
  check_options(options);
  const std::vector<double>& g = system.wavelength_weights();
  std::vector<Input> inputs;
  for (const Source& src : sources) {
    check_source(src, wls);
    Input in;
    in.generated = &src;
    if (src.wavelengths.empty()) {
      for (std::size_t w = 0; w < wls; ++w) in.wavelengths.push_back(static_cast<std::uint16_t>(w));
    } else {
      in.wavelengths = src.wavelengths;
    }
    double gsum = 0.0;
    for (const std::uint16_t w : in.wavelengths) gsum += g[w];
    if (!(gsum > 0.0))
      throw std::invalid_argument("nonseq::trace: the source wavelengths have weight 0");
    in.m = source_ray_count(src);
    // Phi_lambda / M with Phi_lambda = Phi g_lambda / sum(g) (Nachtrag PR 4, g).
    for (const std::uint16_t w : in.wavelengths) {
      in.ray_power.push_back(src.power * g[w] / gsum / static_cast<double>(in.m));
    }
    in.count = in.m * in.wavelengths.size();
    inputs.push_back(std::move(in));
  }
  return run_trace(system, scene, std::move(inputs), detectors, options, control);
}

NonseqResult trace(const compile::CompiledSystem& system,
                   const Scene& scene,
                   std::span<const RaySource> sources,
                   std::span<const Detector> detectors,
                   const TraceOptions& options,
                   const trace::RunControl& control) {
  const std::size_t wls = system.wavelengths_um().size();
  if (sources.size() > 0xffff)
    throw std::invalid_argument("nonseq::trace: more than 65535 sources");
  check_options(options);
  std::vector<Input> inputs;
  for (const RaySource& src : sources) {
    for (const StartRay& r : src.rays) {
      if (r.wavelength >= wls)
        throw std::invalid_argument("nonseq::trace: wavelength out of range");
      if (!(std::isfinite(r.power) && r.power >= 0.0)) {
        throw std::invalid_argument("nonseq::trace: power negative or not finite");
      }
      if (!r.direction.allFinite() || !(r.direction.norm() > 0.0) || !r.position.allFinite()) {
        throw std::invalid_argument("nonseq::trace: direction zero, or ray not finite");
      }
    }
    Input in;
    in.list = &src;
    in.count = src.rays.size();
    inputs.push_back(std::move(in));
  }
  return run_trace(system, scene, std::move(inputs), detectors, options, control);
}

NonseqResult trace(const compile::CompiledSystem& system,
                   const Scene& scene,
                   std::span<const RaySource> sources,
                   const TraceOptions& options,
                   const trace::RunControl& control) {
  return trace(system, scene, sources, std::span<const Detector>{}, options, control);
}

}  // namespace rtt::nonseq
