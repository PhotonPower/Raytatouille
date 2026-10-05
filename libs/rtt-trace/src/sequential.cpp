#include "rtt/trace/sequential.hpp"

#include <oneapi/tbb/blocked_range.h>
#include <oneapi/tbb/parallel_for.h>
#include <oneapi/tbb/partitioner.h>

#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include "rtt/trace/apply_event.hpp"

namespace rtt::trace {

TraceStats SequentialTracer::trace(const compile::CompiledSystem& system,
                                   compile::PathId path,
                                   RayBatch& rays) const {
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
      const bool from_substrate = event.medium_before == c->substrate_medium;
      m.layers = from_substrate
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
      for (const compile::CompiledEvent& event : events.events) {
        if (ray.status != RayStatus::Alive) {
          break;
        }
        ray = sequential_step(ray, system.surfaces()[event.surface], event.surface, event.kind,
                              event_media(event, wl));
      }
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
  oneapi::tbb::parallel_for(oneapi::tbb::blocked_range<std::size_t>(0, rays.size()), trace_range,
                            oneapi::tbb::static_partitioner());

  TraceStats stats;
  for (const RayStatus s : rays.status()) {
    ++stats.rays[static_cast<std::size_t>(s)];
  }
  return stats;
}

}  // namespace rtt::trace
