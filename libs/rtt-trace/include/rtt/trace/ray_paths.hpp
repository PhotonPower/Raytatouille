#pragma once

/// @file ray_paths.hpp
/// Recorded paths of traced rays (#80, docs/architecture.md, Engine 2): the state of every
/// recorded ray before the first event and after each event of the path, for drawing rays in a
/// layout. Filled by SequentialTracer::trace(system, path, rays, paths, ...); a trace without a
/// RayPaths argument records nothing and gives bitwise the same rays.
///
/// Layout (structure of arrays, C order, so that Python gets views without a copy): S = slots()
/// = number of path events + 1; slot 0 is the start state, slot k the state after event k - 1.
/// For recorded ray r (r < ray_count()) and slot s, position and direction are at
/// [(r * S + s) * 3 + c], opl, weight and status at [r * S + s].
///
/// Lost rays: a ray that stops at event j (status no longer Alive) has count = j + 2 valid slots
/// and lost_at = j; slot j + 1 holds the state returned for that event (the hit point for
/// Vignetted, Tir, Absorbed, EventImpossible, Evanescent, the unchanged state for Missed and
/// NoConvergence) with the loss status. The slots after it have NaN in position, direction, opl and
/// weight and repeat the loss status. A ray that is not Alive at the start has count = 1 and
/// lost_at = -1; a ray that passes all events has count = S and lost_at = -1. The last valid slot,
/// count - 1, is bitwise the state that the trace leaves in the RayBatch.
///
/// Units and frames as RayBatch: positions in mm and unit directions in global coordinates,
/// OPL in mm, weight dimensionless (ADR 0021).

#include <cstddef>
#include <cstdint>
#include <vector>

#include "rtt/math/types.hpp"
#include "rtt/trace/ray_batch.hpp"

namespace rtt::trace {

/// Default limit of SequentialTracer::trace() on the number of rays recorded without an
/// explicit selection: 10000 rays of S slots need about 10000 * S * 65 bytes.
inline constexpr std::size_t kDefaultMaxRecordedRays = 10000;

/// Bytes recorded per ray and slot: position, direction (6 doubles), opl, weight (2 doubles),
/// status (1 byte).
inline constexpr std::size_t kRecordedBytesPerSlot = 8 * sizeof(double) + 1;

/// Recorded path states of selected rays of a RayBatch (see the file comment).
struct RayPaths {
  std::size_t slots = 0;                      ///< S = number of path events + 1
  std::vector<std::size_t> ray_indices;       ///< RayBatch index of recorded ray r
  std::vector<std::uint32_t> event_surfaces;  ///< surface index of event k, S - 1 entries
  std::vector<double> position;               ///< mm, global, [(r * S + s) * 3 + c]
  std::vector<double> direction;              ///< unit, global, [(r * S + s) * 3 + c]
  std::vector<double> opl;                    ///< mm, [r * S + s]
  std::vector<double> weight;                 ///< power, unpolarized source = 1, [r * S + s]
  std::vector<RayStatus> status;              ///< [r * S + s]
  std::vector<std::uint32_t> count;           ///< number of valid slots of ray r
  std::vector<std::int32_t> lost_at;          ///< event where ray r stopped, -1 if none

  /// Number of recorded rays.
  [[nodiscard]] std::size_t ray_count() const noexcept { return ray_indices.size(); }

  /// Position of recorded ray r in slot s, mm, global.
  [[nodiscard]] math::Vec3 position_at(std::size_t r, std::size_t s) const {
    const std::size_t i = (r * slots + s) * 3;
    return {position[i], position[i + 1], position[i + 2]};
  }

  /// Unit direction of recorded ray r in slot s, global.
  [[nodiscard]] math::Vec3 direction_at(std::size_t r, std::size_t s) const {
    const std::size_t i = (r * slots + s) * 3;
    return {direction[i], direction[i + 1], direction[i + 2]};
  }
};

}  // namespace rtt::trace
