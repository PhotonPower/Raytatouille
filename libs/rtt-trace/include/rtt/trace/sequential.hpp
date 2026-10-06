#pragma once

/// @file sequential.hpp
/// Sequential and multi-sequential tracing of a RayBatch along one compiled path
/// (docs/architecture.md, Engine 2).

#include <array>
#include <cstddef>
#include <span>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/ray_paths.hpp"

namespace rtt::trace {

/// Number of values of RayStatus.
inline constexpr std::size_t kRayStatusCount = 7;
static_assert(kRayStatusCount == static_cast<std::size_t>(RayStatus::EventImpossible) + 1,
              "kRayStatusCount must match RayStatus");

/// Result summary of a trace: number of rays per status after the trace.
struct TraceStats {
  std::array<std::size_t, kRayStatusCount> rays{};  ///< indexed by RayStatus

  /// Number of rays with status `s`.
  [[nodiscard]] std::size_t count(RayStatus s) const noexcept {
    return rays[static_cast<std::size_t>(s)];
  }
};

/// Common interface of the trace engines (docs/architecture.md).
class Tracer {
 public:
  virtual ~Tracer() = default;

  /// Traces `rays` along path `path` of `system` in place.
  [[nodiscard]] virtual TraceStats trace(const compile::CompiledSystem& system,
                                         compile::PathId path,
                                         RayBatch& rays) const = 0;

 protected:
  Tracer() = default;
  Tracer(const Tracer&) = default;
  Tracer(Tracer&&) noexcept = default;
  Tracer& operator=(const Tracer&) = default;
  Tracer& operator=(Tracer&&) noexcept = default;
};

/// Sequential tracer: every ray visits the events of the path in order via sequential_step()
/// (intersect_surface, aperture check -> Vignetted, apply_event).
///
/// Positions (mm) and unit directions of the batch are in global coordinates (right-handed,
/// optical axis +z); |dir| = 1 is required so that OPL is in mm. After the trace, pos is the last
/// valid point of each ray, dir its direction there, opl the accumulated optical path in mm and
/// last_surface the last surface reached. Rays whose status is not Alive at the start are
/// skipped. Each ray is traced independently;
/// the batch is split with oneTBB parallel_for and static_partitioner (ADR 0004), so the result
/// does not depend on the number of threads. The medium indices are the real parts of the
/// compiled media at the ray's wavelength.
class SequentialTracer final : public Tracer {
 public:
  /// @throws std::out_of_range if `path` does not belong to `system`
  /// @throws std::invalid_argument if a ray's wavelength index is not a system wavelength or
  ///         its status is not a valid RayStatus
  [[nodiscard]] TraceStats trace(const compile::CompiledSystem& system,
                                 compile::PathId path,
                                 RayBatch& rays) const override;

  /// trace() that also records the path of selected rays into `paths` (rtt/trace/ray_paths.hpp,
  /// #80). The rays are traced exactly as without recording: the batch ends bitwise the same.
  /// `paths` is replaced on success and unchanged if an exception is thrown; `rays` is unchanged
  /// if an input check fails (all checks run before tracing).
  /// @param record_rays       RayBatch indices of the rays to record, in this order; empty:
  ///                          all rays (the Python API rejects an empty selection instead)
  /// @param max_recorded_rays limit on the number of rays recorded without a selection
  /// @throws std::invalid_argument as trace(), or if `record_rays` has an index >= rays.size()
  ///         or a duplicate, or if it is empty and rays.size() > max_recorded_rays (the message
  ///         gives the memory need and points to record_rays)
  [[nodiscard]] TraceStats trace(const compile::CompiledSystem& system,
                                 compile::PathId path,
                                 RayBatch& rays,
                                 RayPaths& paths,
                                 std::span<const std::size_t> record_rays = {},
                                 std::size_t max_recorded_rays = kDefaultMaxRecordedRays) const;
};

}  // namespace rtt::trace
