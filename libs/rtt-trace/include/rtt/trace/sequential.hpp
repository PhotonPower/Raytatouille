#pragma once

/// @file sequential.hpp
/// Sequential and multi-sequential tracing of a RayBatch along one compiled path
/// (docs/architecture.md, Engine 2).

#include <array>
#include <cstddef>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/trace/ray_batch.hpp"

namespace rtt::trace {

/// Number of values of RayStatus.
inline constexpr std::size_t kRayStatusCount = 7;

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

/// Sequential tracer: every ray visits the events of the path in order via apply_event().
///
/// Rays whose status is not Alive at the start are skipped. Each ray is traced independently;
/// the batch is split with oneTBB parallel_for and static_partitioner (ADR 0004), so the result
/// does not depend on the number of threads. The medium indices are the real parts of the
/// compiled media at the ray's wavelength.
class SequentialTracer final : public Tracer {
 public:
  /// @throws std::out_of_range if `path` does not belong to `system`
  /// @throws std::invalid_argument if a ray's wavelength index is not a system wavelength
  [[nodiscard]] TraceStats trace(const compile::CompiledSystem& system,
                                 compile::PathId path,
                                 RayBatch& rays) const override;
};

}  // namespace rtt::trace
