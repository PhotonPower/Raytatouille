#include "rtt/trace/sequential.hpp"

#include <oneapi/tbb/blocked_range.h>
#include <oneapi/tbb/parallel_for.h>
#include <oneapi/tbb/partitioner.h>

#include <stdexcept>
#include <string>

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
  }

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
      const std::uint16_t wl = rays.wl()[i];
      for (const compile::CompiledEvent& event : events.events) {
        if (ray.status != RayStatus::Alive) {
          break;
        }
        // M1: real part of the complex index (absorption follows in M3).
        const double n_before = system.media()[event.medium_before].index[wl].real();
        const double n_after = system.media()[event.medium_after].index[wl].real();
        ray = apply_event(ray, system.surfaces()[event.surface], event.surface, event.kind,
                          n_before, n_after);
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
