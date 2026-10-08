#include "rtt/trace/sequential.hpp"

#include <oneapi/tbb/blocked_range.h>
#include <oneapi/tbb/parallel_for.h>
#include <oneapi/tbb/partitioner.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#include "rtt/trace/apply_event.hpp"

namespace rtt::trace {
namespace {

/// Diffraction orders and efficiencies are traced from #127 on (ADR 0025). Until then an event
/// with order != 0 ends the ray at its hit point with EventImpossible, as the event kind Diffract
/// did before schema 0.3. So does every event at a surface with diffraction_efficiency: there
/// orders not listed, also order 0, have efficiency 0 (ADR 0025, point 5), so passing with full
/// weight would be wrong. A miss, vignetting and an absorber keep their status, as in
/// sequential_step().
[[nodiscard]] bool order_not_traced(const compile::CompiledEvent& event,
                                    const compile::CompiledSurface& surface) noexcept {
  return event.order != 0 || surface.diffraction_efficiency.has_value();
}

/// The event that order_not_traced() selects: stops at the hit point (see there).
RayState stop_at_order(const RayState& ray,
                       const compile::CompiledSurface& surface,
                       std::uint32_t surface_index,
                       const EventMedia& media) noexcept {
  if (ray.status != RayStatus::Alive) return ray;
  const SurfaceHit hit = intersect_surface(ray, surface);
  RayState out = ray;
  switch (hit.status) {
    case geom::HitStatus::Hit:
      break;
    case geom::HitStatus::Missed:
      out.status = RayStatus::Missed;
      return out;
    case geom::HitStatus::NoConvergence:
      out.status = RayStatus::NoConvergence;
      return out;
  }
  out = move_to_hit(ray, surface, hit, surface_index, media);
  if (!inside_aperture(surface, hit)) {
    out.status = RayStatus::Vignetted;
  } else if (std::holds_alternative<model::Absorber>(surface.interaction)) {
    out.status = RayStatus::Absorbed;
    out.weight = 0.0;
  } else {
    out.status = RayStatus::EventImpossible;
  }
  return out;
}

/// Recorder policy of trace_rays() without recording: every call is an empty inline function,
/// so the ray loop of a plain trace is the same as before #80 (results bitwise unchanged).
struct NoRecord {
  void start(std::size_t /*ray*/, const RayState& /*state*/) const noexcept {}
  void after(std::size_t /*ray*/, std::size_t /*event*/, const RayState& /*state*/) const noexcept {
  }
  void finish(std::size_t /*ray*/,
              std::size_t /*steps*/,
              const RayState& /*state*/) const noexcept {}
};

/// Recorder policy that writes the slots of the selected rays into RayPaths (ray_paths.hpp).
/// Every ray writes only its own row, so the parallel loop stays deterministic (ADR 0004).
/// The arrays are allocated and NaN-filled before the loop; the loop does not allocate.
class Record {
 public:
  /// @param row recorded row of every RayBatch index, -1 if the ray is not recorded
  Record(RayPaths& paths, std::span<const std::int64_t> row) : paths_(paths), row_(row) {}

  /// Slot 0: the start state.
  void start(std::size_t ray, const RayState& state) const noexcept {
    if (row_[ray] >= 0) write(static_cast<std::size_t>(row_[ray]), 0, state);
  }

  /// Slot event + 1: the state after the event.
  void after(std::size_t ray, std::size_t event, const RayState& state) const noexcept {
    if (row_[ray] >= 0) write(static_cast<std::size_t>(row_[ray]), event + 1, state);
  }

  /// After `steps` events: count, lost_at and the status of the unused slots (their values
  /// stay NaN).
  void finish(std::size_t ray, std::size_t steps, const RayState& state) const noexcept {
    if (row_[ray] < 0) return;
    const auto r = static_cast<std::size_t>(row_[ray]);
    paths_.count[r] = static_cast<std::uint32_t>(steps + 1);
    // Stopped at the last executed event if it is no longer Alive; a ray that did not start
    // (steps = 0) or passed every event was not lost on the path.
    paths_.lost_at[r] =
        steps > 0 && state.status != RayStatus::Alive ? static_cast<std::int32_t>(steps - 1) : -1;
    for (std::size_t s = steps + 1; s < paths_.slots; ++s) {
      paths_.status[r * paths_.slots + s] = state.status;
    }
  }

 private:
  void write(std::size_t r, std::size_t slot, const RayState& state) const noexcept {
    const std::size_t i = r * paths_.slots + slot;
    for (int c = 0; c < 3; ++c) {
      paths_.position[3 * i + static_cast<std::size_t>(c)] = state.pos[c];
      paths_.direction[3 * i + static_cast<std::size_t>(c)] = state.dir[c];
    }
    paths_.opl[i] = state.opl;
    paths_.weight[i] = state.weight;
    paths_.status[i] = state.status;
  }

  RayPaths& paths_;
  std::span<const std::int64_t> row_;
};

/// The sequential ray loop shared by both trace() overloads; `recorder` sees the start state,
/// the state after every executed event and the final state of each ray.
/// With an active `control` the batch is split into blocks of control->block_size rays
/// (run_control.hpp, #83); every ray is still traced exactly once by the same code.
template <class Recorder>
TraceStats trace_rays(const compile::CompiledSystem& system,
                      compile::PathId path,
                      RayBatch& rays,
                      const Recorder& recorder,
                      const RunControl* control = nullptr) {
  // Input checks at the API boundary (ADR 0009); the tracing loop itself never throws.
  const compile::CompiledPath& events = system.path(path);
  const std::size_t wavelengths = system.wavelengths_um().size();
  for (std::size_t i = 0; i < rays.size(); ++i) {
    if (rays.wl()[i] >= wavelengths) {
      throw std::invalid_argument("ray " + std::to_string(i) + ": wavelength index " +
                                  std::to_string(rays.wl()[i]) + " is not a system wavelength");
    }
    if (static_cast<std::size_t>(rays.status()[i]) >= kRayStatusCount) {
      throw std::invalid_argument("ray " + std::to_string(i) + ": invalid status");
    }
  }

  // Coating stacks as seen from the substrate side (ADR 0019: light from inside sees the layers
  // in reverse order), prepared once so that the tracing loop does not allocate.
  const auto& coatings = system.coatings();
  std::vector<std::vector<std::vector<coating::Layer<double>>>> reversed(coatings.size());
  for (std::size_t c = 0; c < coatings.size(); ++c) {
    for (const auto& layers : coatings[c].layers) {
      reversed[c].emplace_back(layers.rbegin(), layers.rend());
    }
  }
  const auto event_media = [&](const compile::CompiledEvent& event, std::uint16_t wl) {
    EventMedia m;
    m.before = system.media()[event.medium_before].index[wl];
    m.after = system.media()[event.medium_after].index[wl];
    m.beyond = system.media()[event.medium_beyond].index[wl];
    m.wavelength_um = system.wavelengths_um()[wl];
    if (const auto& c = system.surfaces()[event.surface].coating) {
      // From inside the element = from the substrate (ADR 0019); decided by element, not by
      // medium, since another element of the same glass may lie in front of the surface.
      m.layers = event.from_inside
                     ? std::span<const coating::Layer<double>>(reversed[c->coating][wl])
                     : std::span<const coating::Layer<double>>(coatings[c->coating].layers[wl]);
    }
    return m;
  };

  // Rays are independent; every ray is written by exactly one task, so the result does not
  // depend on the partitioning (ADR 0004).
  const auto trace_range = [&](const oneapi::tbb::blocked_range<std::size_t>& range) {
    for (std::size_t i = range.begin(); i != range.end(); ++i) {
      RayState ray;
      ray.pos = {rays.pos_x()[i], rays.pos_y()[i], rays.pos_z()[i]};
      ray.dir = {rays.dir_x()[i], rays.dir_y()[i], rays.dir_z()[i]};
      ray.opl = rays.opl()[i];
      ray.status = rays.status()[i];
      ray.last_surface = rays.last_surface()[i];
      ray.prt = rays.prt_matrix(i);
      ray.weight = rays.weight()[i];
      const std::uint16_t wl = rays.wl()[i];
      recorder.start(i, ray);
      std::size_t steps = 0;
      for (; steps < events.events.size() && ray.status == RayStatus::Alive; ++steps) {
        const compile::CompiledEvent& event = events.events[steps];
        const compile::CompiledSurface& surface = system.surfaces()[event.surface];
        ray =
            order_not_traced(event, surface)
                ? stop_at_order(ray, surface, event.surface, event_media(event, wl))
                : sequential_step(ray, surface, event.surface, event.kind, event_media(event, wl));
        recorder.after(i, steps, ray);
      }
      recorder.finish(i, steps, ray);
      rays.pos_x()[i] = ray.pos.x();
      rays.pos_y()[i] = ray.pos.y();
      rays.pos_z()[i] = ray.pos.z();
      rays.dir_x()[i] = ray.dir.x();
      rays.dir_y()[i] = ray.dir.y();
      rays.dir_z()[i] = ray.dir.z();
      rays.opl()[i] = ray.opl;
      rays.status()[i] = ray.status;
      rays.last_surface()[i] = ray.last_surface;
      rays.set_prt_matrix(i, ray.prt);
      rays.weight()[i] = ray.weight;
    }
  };
  if (control == nullptr || !control->active()) {
    oneapi::tbb::parallel_for(oneapi::tbb::blocked_range<std::size_t>(0, rays.size()), trace_range,
                              oneapi::tbb::static_partitioner());
  } else {
    // Blocks of at most block_size rays (simple_partitioner splits down to the grain size): no
    // new block after a cancellation request or a failed progress callback; the monitor throws
    // only in finish(), after the parallel part (rule 3). ADR 0004, addendum #83.
    RunMonitor monitor(*control, rays.size(), "trace");
    const std::size_t grain = std::max<std::size_t>(control->block_size, 1);
    oneapi::tbb::parallel_for(
        oneapi::tbb::blocked_range<std::size_t>(0, rays.size(), grain),
        [&](const oneapi::tbb::blocked_range<std::size_t>& range) {
          if (monitor.stop()) return;
          trace_range(range);
          monitor.add(range.size());
        },
        oneapi::tbb::simple_partitioner());
    monitor.finish();
  }

  TraceStats stats;
  for (const RayStatus s : rays.status()) {
    ++stats.rays[static_cast<std::size_t>(s)];
  }
  return stats;
}

}  // namespace

TraceStats SequentialTracer::trace(const compile::CompiledSystem& system,
                                   compile::PathId path,
                                   RayBatch& rays,
                                   const RunControl& control) const {
  return trace_rays(system, path, rays, NoRecord{}, &control);
}

TraceStats SequentialTracer::trace(const compile::CompiledSystem& system,
                                   compile::PathId path,
                                   RayBatch& rays) const {
  return trace_rays(system, path, rays, NoRecord{});
}

TraceStats SequentialTracer::trace(const compile::CompiledSystem& system,
                                   compile::PathId path,
                                   RayBatch& rays,
                                   RayPaths& paths,
                                   std::optional<std::span<const std::size_t>> record_rays,
                                   std::size_t max_recorded_rays) const {
  // Input checks at the API boundary (ADR 0009), before anything is traced.
  const compile::CompiledPath& events = system.path(path);
  const std::size_t slots = events.events.size() + 1;
  const std::size_t n = rays.size();
  std::vector<std::int64_t> row(n, -1);
  std::vector<std::size_t> selected;
  if (record_rays && record_rays->empty()) {
    throw std::invalid_argument(
        "record_rays is empty; pass no selection (std::nullopt, None in Python) to record all "
        "rays");
  }
  if (!record_rays) {
    if (n > max_recorded_rays) {
      throw std::invalid_argument(
          "recording " + std::to_string(n) + " rays of " + std::to_string(slots) +
          " slots needs about " + std::to_string(n * slots * kRecordedBytesPerSlot) +
          " bytes, more than max_recorded_rays = " + std::to_string(max_recorded_rays) +
          " allows; select the rays to record with record_rays (or raise max_recorded_rays)");
    }
    selected.resize(n);
    for (std::size_t i = 0; i < n; ++i) selected[i] = i;
  } else {
    selected.assign(record_rays->begin(), record_rays->end());
  }
  for (std::size_t r = 0; r < selected.size(); ++r) {
    const std::size_t i = selected[r];
    if (i >= n) {
      throw std::invalid_argument("record_rays[" + std::to_string(r) + "] = " + std::to_string(i) +
                                  " is not a ray index (batch of " + std::to_string(n) + " rays)");
    }
    if (row[i] >= 0) {
      throw std::invalid_argument("record_rays lists ray " + std::to_string(i) + " twice");
    }
    row[i] = static_cast<std::int64_t>(r);
  }

  // Recorded into a local object and moved into `paths` only after the trace succeeded, so that
  // `paths` is unchanged if anything throws (strong guarantee; also for bad_alloc).
  constexpr double nan = std::numeric_limits<double>::quiet_NaN();
  const std::size_t recorded = selected.size();
  RayPaths result;
  result.slots = slots;
  result.ray_indices = std::move(selected);
  result.event_surfaces.reserve(events.events.size());
  for (const compile::CompiledEvent& event : events.events) {
    result.event_surfaces.push_back(event.surface);
  }
  result.position.assign(recorded * slots * 3, nan);
  result.direction.assign(recorded * slots * 3, nan);
  result.opl.assign(recorded * slots, nan);
  result.weight.assign(recorded * slots, nan);
  result.status.assign(recorded * slots, RayStatus::Alive);
  result.count.assign(recorded, 0);
  result.lost_at.assign(recorded, -1);
  const TraceStats stats = trace_rays(system, path, rays, Record(result, row));
  paths = std::move(result);
  return stats;
}

}  // namespace rtt::trace
