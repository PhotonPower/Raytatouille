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
/// Polarization (ADR 0034, point 3): with `polarization` set before the trace, prt and q hold the
/// accumulated P and Q of every slot, 9 values each at [(r * S + s) * 9 + 3 * row + col]
/// (row-major), NaN in the unused slots like the other values. Slot 0 is the state of the input
/// batch (I for rays from a source); slot s the product after event s - 1; P and Q of the last
/// valid slot are bitwise those of the RayBatch. The P of one event follows from two slots. Without
/// the request prt and q stay empty and everything else is bitwise the same.
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

/// Additional bytes per ray and slot with RayPaths::polarization: P (9 complex) and Q (9 real),
/// 216 bytes (ADR 0034, point 3).
inline constexpr std::size_t kRecordedPolarizationBytesPerSlot =
    9 * sizeof(math::Complex) + 9 * sizeof(double);

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
  /// Request to record P and Q per slot, set by the caller before the trace; kept by it.
  bool polarization = false;
  std::vector<math::Complex> prt;  ///< P per slot (RayBatch::prt), only with `polarization`
  std::vector<double> q;           ///< Q per slot (RayBatch::q), only with `polarization`

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

  /// Accumulated P of recorded ray r in slot s (power-normalised, global; ADR 0021).
  /// @pre polarization was set for the trace (prt is not empty)
  [[nodiscard]] math::CMat3 prt_at(std::size_t r, std::size_t s) const {
    const std::size_t i = (r * slots + s) * 9;
    math::CMat3 p;
    for (int row = 0; row < 3; ++row) {
      for (int col = 0; col < 3; ++col) {
        p(row, col) = prt[i + static_cast<std::size_t>(3 * row + col)];
      }
    }
    return p;
  }

  /// Accumulated Q of recorded ray r in slot s (global; ADR 0034, point 1).
  /// @pre polarization was set for the trace (q is not empty)
  [[nodiscard]] math::Mat3 q_at(std::size_t r, std::size_t s) const {
    const std::size_t i = (r * slots + s) * 9;
    math::Mat3 m;
    for (int row = 0; row < 3; ++row) {
      for (int col = 0; col < 3; ++col)
        m(row, col) = q[i + static_cast<std::size_t>(3 * row + col)];
    }
    return m;
  }
};

}  // namespace rtt::trace
