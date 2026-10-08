#pragma once

/// @file common.hpp
/// Internal helpers shared by the analyses of rtt-analysis (not installed, not public API).

#include <array>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <map>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "rtt/analysis/spot.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/diagnostics/codes.hpp"
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

/// Accumulates RayLosses over the batches of one analysis and derives its warnings.
class LossCounter {
 public:
  /// @throws std::invalid_argument if `lost_warning_fraction` is not in [0, 1]
  LossCounter(const CompiledSystem& system, PathId path, double lost_warning_fraction)
      : system_(&system), threshold_(lost_warning_fraction) {
    if (!(lost_warning_fraction >= 0.0 && lost_warning_fraction <= 1.0)) {
      throw std::invalid_argument("analysis: lost_warning_fraction must lie in [0, 1]");
    }
    for (const auto& event : system.path(path).events) {
      if (system.surfaces()[event.surface].element_kind == model::ElementKind::Stop) {
        stop_ = event.surface;
        break;
      }
    }
  }

  /// Counts every ray of `rays`; a ray that ended Alive elsewhere than on `image` counts as
  /// Vignetted (as in the fans).
  void add(const trace::RayBatch& rays, std::uint32_t image) {
    for (std::size_t i = 0; i < rays.size(); ++i) {
      trace::RayStatus s = rays.status()[i];
      if (s == trace::RayStatus::Alive && rays.last_surface()[i] != image) {
        s = trace::RayStatus::Vignetted;
      }
      add_ray(s, rays.last_surface()[i]);
    }
  }

  /// Counts one ray with its final status as the analysis sees it (e.g. Vignetted for an OPD
  /// ray that misses the reference sphere) and its RayBatch::last_surface.
  void add_ray(trace::RayStatus s, std::uint32_t last_surface) {
    ++launched_;
    ++by_status_[static_cast<std::size_t>(s)];
    if (s != trace::RayStatus::Alive && last_surface != trace::kNoSurface) {
      ++lost_at_[last_surface];
    }
    if (s == trace::RayStatus::Vignetted && stop_ && last_surface == *stop_) {
      ++vignetted_at_stop_;
    }
  }

  [[nodiscard]] RayLosses result() const {
    // Built field by field: copying a stored RayLosses with an empty optional made GCC 13 warn
    // "may be used uninitialized" in Release (CI of #104).
    RayLosses l;
    l.launched = launched_;
    l.by_status = by_status_;
    std::uint32_t worst = 0;
    for (const auto& [surface, n] : lost_at_) {  // ascending index: ties keep the lowest
      if (n > l.worst_surface_count) {
        worst = surface;
        l.worst_surface_count = n;
      }
    }
    if (l.worst_surface_count > 0) l.worst_surface = worst;
    return l;
  }

  /// "rays.lost" if the lost fraction exceeds the threshold, at the worst loss surface;
  /// "stop.clips_beam" if rays ended Vignetted at the stop, at the stop surface (ADR 0023).
  [[nodiscard]] std::vector<model::Diagnostic> warnings() const {
    std::vector<model::Diagnostic> out;
    const RayLosses l = result();
    const std::size_t lost = l.launched - l.count(trace::RayStatus::Alive);
    const auto name = [&](std::uint32_t s) { return "'" + system_->surfaces()[s].id.str() + "'"; };
    if (l.launched > 0 &&
        static_cast<double>(lost) > threshold_ * static_cast<double>(l.launched)) {
      std::string message = std::to_string(lost) + " of " + std::to_string(l.launched) +
                            " rays lost, more than the threshold " + percent(threshold_);
      std::string location;
      if (l.worst_surface) {
        message += "; the most (" + std::to_string(l.worst_surface_count) + ") ended at surface " +
                   name(*l.worst_surface);
        location = system_->surfaces()[*l.worst_surface].location;
      }
      add(out, "rays.lost", std::move(location), std::move(message));
    }
    if (stop_ && vignetted_at_stop_ > 0) {  // only counted with a stop
      add(out, "stop.clips_beam", system_->surfaces()[*stop_].location,
          std::to_string(vignetted_at_stop_) + " rays ended Vignetted at the stop " + name(*stop_) +
              ": its aperture clips the beam that the system aperture defines");
    }
    return out;
  }

 private:
  static void add(std::vector<model::Diagnostic>& out,
                  diagnostics::DiagnosticCode code,
                  std::string location,
                  std::string message) {
    out.push_back(
        {code.severity(), std::move(location), std::move(message), std::string(code.str())});
  }

  static std::string percent(double fraction) {
    std::ostringstream text;
    text << std::setprecision(3) << fraction * 100.0 << " %";
    return text.str();
  }

  const CompiledSystem* system_;
  double threshold_;
  std::optional<std::uint32_t> stop_;
  std::size_t vignetted_at_stop_ = 0;
  std::size_t launched_ = 0;
  std::array<std::size_t, trace::kRayStatusCount> by_status_{};
  std::map<std::uint32_t, std::size_t> lost_at_;
};

/// AnalysisError for ray i that did not arrive, with its surface, status, field and wavelength.
[[noreturn]] inline void throw_lost(const CompiledSystem& system,
                                    const trace::RayBatch& rays,
                                    std::size_t i,
                                    std::optional<std::uint16_t> field,
                                    const std::string& message) {
  AnalysisError::LostRay ray;
  const std::uint32_t last = rays.last_surface()[i];
  if (last != trace::kNoSurface) {
    ray.surface = system.surfaces()[last].id;
    ray.location = system.surfaces()[last].location;
  }
  ray.ray_status = rays.status()[i];
  ray.field = field;
  ray.wavelength = rays.wl()[i];
  throw AnalysisError(message, std::move(ray));
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
