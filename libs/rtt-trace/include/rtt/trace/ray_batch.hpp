#pragma once

/// @file ray_batch.hpp
/// Batch of real rays as Structure-of-Arrays (ADR 0005; docs/architecture.md, "Strahl-Batch").
///
/// Conventions (docs/architecture.md, "Konventionen"): positions in mm and directions as unit
/// vectors, both in global coordinates (right-handed, optical axis +z); optical path length in
/// mm; weight as power with the source normalised to 1.

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <vector>

#include "rtt/math/types.hpp"

namespace rtt::trace {

/// State of a ray. Problems while tracing are status flags, never exceptions (ADR 0009).
enum class RayStatus : std::uint8_t {
  Alive,            ///< still being traced
  Missed,           ///< did not intersect the next surface
  Vignetted,        ///< hit a surface outside its aperture
  Tir,              ///< total internal reflection where refraction was requested
  NoConvergence,    ///< iterative intersection did not converge
  Absorbed,         ///< stopped by an absorbing surface
  EventImpossible,  ///< the requested event cannot happen (e.g. evanescent order)
  Evanescent,  // STUB for the red run of #127 PR A
};

/// Value of `last_surface` before a ray has hit any surface.
inline constexpr std::uint32_t kNoSurface = std::numeric_limits<std::uint32_t>::max();

/// Rays stored component-wise: every field is one contiguous array of length size().
///
/// A newly created or grown ray starts at the global origin travelling along +z with
/// wavelength index 0, OPL 0, polarization ray-tracing matrix P = identity, weight 1, status
/// Alive, field index 0, pupil coordinates (0, 0) and last_surface = kNoSurface.
/// Callers set the fields they need; the batch never checks physical consistency
/// (e.g. |dir| = 1), which is the job of whoever fills it.
class RayBatch {
 public:
  /// Empty batch.
  RayBatch() = default;

  /// Batch of `size` rays in the default state described above.
  explicit RayBatch(std::size_t size);

  /// Number of rays.
  [[nodiscard]] std::size_t size() const noexcept { return status_.size(); }

  /// Changes the number of rays; new rays get the default state, existing rays are kept.
  void resize(std::size_t size);

  // Position in global coordinates, mm.
  [[nodiscard]] std::span<double> pos_x() noexcept { return pos_x_; }
  [[nodiscard]] std::span<double> pos_y() noexcept { return pos_y_; }
  [[nodiscard]] std::span<double> pos_z() noexcept { return pos_z_; }
  [[nodiscard]] std::span<const double> pos_x() const noexcept { return pos_x_; }
  [[nodiscard]] std::span<const double> pos_y() const noexcept { return pos_y_; }
  [[nodiscard]] std::span<const double> pos_z() const noexcept { return pos_z_; }

  // Unit direction in global coordinates.
  [[nodiscard]] std::span<double> dir_x() noexcept { return dir_x_; }
  [[nodiscard]] std::span<double> dir_y() noexcept { return dir_y_; }
  [[nodiscard]] std::span<double> dir_z() noexcept { return dir_z_; }
  [[nodiscard]] std::span<const double> dir_x() const noexcept { return dir_x_; }
  [[nodiscard]] std::span<const double> dir_y() const noexcept { return dir_y_; }
  [[nodiscard]] std::span<const double> dir_z() const noexcept { return dir_z_; }

  /// Index into rtt::compile::CompiledSystem::wavelengths_um().
  [[nodiscard]] std::span<std::uint16_t> wl() noexcept { return wl_; }
  [[nodiscard]] std::span<const std::uint16_t> wl() const noexcept { return wl_; }

  /// Accumulated optical path length, mm.
  [[nodiscard]] std::span<double> opl() noexcept { return opl_; }
  [[nodiscard]] std::span<const double> opl() const noexcept { return opl_; }

  /// Element (row, col) of the accumulated 3x3 polarization ray-tracing matrix P
  /// (Yun, McClain, Chipman, "Three-dimensional polarization ray-tracing calculus I",
  /// Applied Optics 50, 2011), global coordinates, POWER-NORMALISED (ADR 0021): for an incident
  /// field E_in with |E_in| = 1 transverse to the initial direction, |P E_in|^2 is the transmitted
  /// power fraction WITHOUT the polarization-independent factors s (volume absorption, absorber),
  /// which only weight() carries, and the phases are those of the field PRT matrix; the field
  /// amplitude differs
  /// by the factors sqrt(c) of the interfaces passed (Byrnes, Eqs. (21), (22)). P k_0 = k.
  /// @pre 0 <= row < 3 and 0 <= col < 3 (not checked)
  [[nodiscard]] std::span<math::Complex> prt(int row, int col) noexcept;
  [[nodiscard]] std::span<const math::Complex> prt(int row, int col) const noexcept;

  /// P of ray `i` as a matrix (copy).
  [[nodiscard]] math::CMat3 prt_matrix(std::size_t i) const;

  /// Sets P of ray `i`.
  void set_prt_matrix(std::size_t i, const math::CMat3& p);

  /// Power for an unpolarized source, dimensionless (source normalised to 1): weight =
  /// s ||P_T||_F^2 / 2 with P_T = P - k k_0^T and s the polarization-independent factors (volume
  /// absorption, absorber). Power for a polarized state E: weight |P E|^2 / (||P_T||_F^2 / 2)
  /// (ADR 0021).
  [[nodiscard]] std::span<double> weight() noexcept { return weight_; }
  [[nodiscard]] std::span<const double> weight() const noexcept { return weight_; }

  /// Ray status.
  [[nodiscard]] std::span<RayStatus> status() noexcept { return status_; }
  [[nodiscard]] std::span<const RayStatus> status() const noexcept { return status_; }

  /// Index of the field point the ray belongs to (for analyses).
  [[nodiscard]] std::span<std::uint16_t> field() noexcept { return field_; }
  [[nodiscard]] std::span<const std::uint16_t> field() const noexcept { return field_; }

  // Normalised pupil coordinates of the ray's origin (for analyses), dimensionless.
  [[nodiscard]] std::span<double> pupil_x() noexcept { return pupil_x_; }
  [[nodiscard]] std::span<double> pupil_y() noexcept { return pupil_y_; }
  [[nodiscard]] std::span<const double> pupil_x() const noexcept { return pupil_x_; }
  [[nodiscard]] std::span<const double> pupil_y() const noexcept { return pupil_y_; }

  /// Index into rtt::compile::CompiledSystem::surfaces() of the last surface hit, or kNoSurface.
  [[nodiscard]] std::span<std::uint32_t> last_surface() noexcept { return last_surface_; }
  [[nodiscard]] std::span<const std::uint32_t> last_surface() const noexcept {
    return last_surface_;
  }

 private:
  std::vector<double> pos_x_, pos_y_, pos_z_;
  std::vector<double> dir_x_, dir_y_, dir_z_;
  std::vector<std::uint16_t> wl_;
  std::vector<double> opl_;
  std::array<std::vector<math::Complex>, 9> prt_;  // row-major: prt_[3 * row + col]
  std::vector<double> weight_;
  std::vector<RayStatus> status_;
  std::vector<std::uint16_t> field_;
  std::vector<double> pupil_x_, pupil_y_;
  std::vector<std::uint32_t> last_surface_;
};

}  // namespace rtt::trace
