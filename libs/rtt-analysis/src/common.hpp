#pragma once

/// @file common.hpp
/// Internal helpers shared by the analyses of rtt-analysis (not installed, not public API).

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>

#include "rtt/analysis/spot.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/sequential.hpp"
#include "rtt/trace/sources.hpp"

namespace rtt::analysis::detail {

using compile::CompiledSystem;
using compile::PathId;

inline void check_path(const CompiledSystem& system, PathId path) {
  if (path.index >= system.paths().size() || system.path(path).events.empty()) {
    throw std::invalid_argument("analysis: path index " + std::to_string(path.index) +
                                " does not exist or has no events");
  }
}

inline void check_wavelength(const CompiledSystem& system, std::uint16_t wavelength) {
  if (wavelength >= system.wavelengths_um().size()) {
    throw std::invalid_argument("analysis: wavelength index " + std::to_string(wavelength) +
                                " does not exist");
  }
}

/// Image surface: the surface of the last event of the path (decided for #28).
inline std::uint32_t image_surface(const CompiledSystem& system, PathId path) {
  return system.path(path).events.back().surface;
}

/// Traces rays of one field and wavelength through the whole path.
inline trace::RayBatch trace_rays(const CompiledSystem& system,
                                  PathId path,
                                  std::uint16_t field,
                                  std::uint16_t wavelength,
                                  const trace::PupilSampling& sampling,
                                  trace::Aiming aiming) {
  const std::uint16_t fields[] = {field};
  trace::RayBatch rays = trace::make_rays(system, path, fields, wavelength, sampling, aiming);
  [[maybe_unused]] const auto stats = trace::SequentialTracer().trace(system, path, rays);
  return rays;
}

/// Traces one aimed ray through the whole path (batch of size 1).
inline trace::RayBatch trace_aimed(const CompiledSystem& system,
                                   PathId path,
                                   const trace::AimedRay& aimed,
                                   std::uint16_t wavelength) {
  trace::RayBatch rays(1);
  rays.pos_x()[0] = aimed.ray.pos.x();
  rays.pos_y()[0] = aimed.ray.pos.y();
  rays.pos_z()[0] = aimed.ray.pos.z();
  rays.dir_x()[0] = aimed.ray.dir.x();
  rays.dir_y()[0] = aimed.ray.dir.y();
  rays.dir_z()[0] = aimed.ray.dir.z();
  rays.wl()[0] = wavelength;
  rays.status()[0] = aimed.ray.status;
  [[maybe_unused]] const auto stats = trace::SequentialTracer().trace(system, path, rays);
  return rays;
}

/// True if ray i ended Alive on the image surface.
inline bool arrived(const trace::RayBatch& rays, std::size_t i, std::uint32_t image) {
  return rays.status()[i] == trace::RayStatus::Alive && rays.last_surface()[i] == image;
}

/// Local x, y of ray i on the image surface, mm.
inline Point2 local_point(const CompiledSystem& system,
                          const trace::RayBatch& rays,
                          std::size_t i,
                          std::uint32_t image) {
  const math::Vec3 p = system.surfaces()[image].to_local.apply_point(
      math::Vec3(rays.pos_x()[i], rays.pos_y()[i], rays.pos_z()[i]));
  return {p.x(), p.y()};
}

}  // namespace rtt::analysis::detail
