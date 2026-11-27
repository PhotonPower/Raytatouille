#include "rtt/nonseq/trace.hpp"

#include <oneapi/tbb/blocked_range.h>
#include <oneapi/tbb/parallel_for.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include "rtt/coating/transfer_matrix.hpp"
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
};

/// Sums of one block.
struct Partial {
  std::vector<Balance> balance;
  std::uint64_t rays = 0;
  std::uint64_t events = 0;
  bool done = false;
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

/// Traces one start ray with all its branches into `b`.
void trace_ray(const Run& run,
               const StartRay& start,
               std::uint32_t source,
               std::uint64_t index,
               Balance& b,
               std::uint64_t& events) {
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
    Branch br = stack.back();
    stack.pop_back();
    for (;;) {
      const SceneHit hit = run.scene.closest_hit(br.ray, br.last_face);
      if (hit.status == geom::HitStatus::Missed) {
        b.escaped += p0 * br.ray.weight;
        break;
      }
      if (hit.status != geom::HitStatus::Hit) {
        b.not_modelled += p0 * br.ray.weight;
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
      const auto absorb = [&] { b.absorbed[hit.face] += p0 * w_in; };

      if (hit.region == HitRegion::Absorbed) {
        absorb();
        break;
      }
      if (face.kind == FaceKind::Detector) {
        if (!face.transparent) {
          b.detected[hit.face] += p0 * w_in;
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
        break;
      }
      b.absorbed[hit.face] += p0 * (w_in - total);  // loss of the event
      if (!(total > 0.0)) break;

      const std::uint32_t depth = br.depth + 1;
      const auto next_mode = [&](EventKind k) {
        if (k == EventKind::Ordinary) return CrystalMode::Ordinary;
        if (k == EventKind::Extraordinary) return CrystalMode::Extraordinary;
        if (inside && k == EventKind::Refract) return CrystalMode::None;
        return br.mode;
      };
      const auto keep = [&](double weight) {
        if (depth > o.max_depth || weight < o.min_relative_power) {
          b.cut += p0 * weight;
          return false;
        }
        return true;
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
        if (!keep(total)) break;
        continue;
      }
      // Full splitting: push the branches in reverse so that the first candidate runs first.
      for (std::size_t i = results.size(); i-- > 0;) {
        if (!(results[i].weight > 0.0)) continue;
        if (!keep(results[i].weight)) continue;
        Branch nb;
        nb.ray = results[i];
        nb.mode = next_mode(cands[i].kind);
        nb.last_face = hit.face;
        nb.depth = depth;
        stack.push_back(nb);
      }
      break;
    }
  }
}

}  // namespace

NonseqResult trace(const compile::CompiledSystem& system,
                   const Scene& scene,
                   std::span<const RaySource> sources,
                   const TraceOptions& options,
                   const trace::RunControl& control) {
  const std::size_t wls = system.wavelengths_um().size();
  if (sources.size() > 0xffff)
    throw std::invalid_argument("nonseq::trace: more than 65535 sources");
  if (options.max_depth == 0) throw std::invalid_argument("nonseq::trace: max_depth must be > 0");
  if (!(std::isfinite(options.min_relative_power) && options.min_relative_power >= 0.0)) {
    throw std::invalid_argument("nonseq::trace: min_relative_power must be finite and >= 0");
  }
  // Global ray index -> (source, index within the source).
  std::vector<std::uint64_t> offsets{0};
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
    offsets.push_back(offsets.back() + src.rays.size());
  }
  const std::uint64_t total = offsets.back();
  const std::size_t faces = scene.faces().size();

  Run run{system, scene, options, wls, {}, {}};
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

  const std::size_t blocks =
      static_cast<std::size_t>((total + kNonseqBlockSize - 1) / kNonseqBlockSize);
  std::vector<Partial> partials(blocks);
  trace::RunMonitor monitor(control, static_cast<std::size_t>(total), "nonsequential trace");
  const auto run_block = [&](std::size_t blk) {
    if (monitor.stop()) return;
    Partial& p = partials[blk];
    p.balance.assign(sources.size() * wls, empty_balance(faces));
    const std::uint64_t begin = blk * kNonseqBlockSize;
    const std::uint64_t end = std::min<std::uint64_t>(begin + kNonseqBlockSize, total);
    auto src = static_cast<std::size_t>(std::upper_bound(offsets.begin(), offsets.end(), begin) -
                                        offsets.begin() - 1);
    for (std::uint64_t g = begin; g < end; ++g) {
      while (g >= offsets[src + 1]) ++src;
      const std::uint64_t index = g - offsets[src];
      const StartRay& start = sources[src].rays[static_cast<std::size_t>(index)];
      Balance& bal = p.balance[src * wls + start.wavelength];
      bal.source += start.power;
      trace_ray(run, start, static_cast<std::uint32_t>(src), index, bal, p.events);
      ++p.rays;
    }
    p.done = true;
    monitor.add(static_cast<std::size_t>(end - begin));
  };
  oneapi::tbb::parallel_for(oneapi::tbb::blocked_range<std::size_t>(0, blocks, 1),
                            [&](const oneapi::tbb::blocked_range<std::size_t>& range) {
                              for (std::size_t blk = range.begin(); blk != range.end(); ++blk) {
                                run_block(blk);
                              }
                            });
  NonseqResult result;
  try {
    monitor.finish();
  } catch (const trace::Cancelled&) {
    result.cancelled = true;
  }
  // Sum in block order (rule 7).
  result.balance.assign(sources.size() * wls, empty_balance(faces));
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
  }
  return result;
}

}  // namespace rtt::nonseq
