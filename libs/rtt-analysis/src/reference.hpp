#pragma once

/// @file reference.hpp
/// Reference sphere of the OPD and the OPL of a ray up to it (internal, not installed): shared by
/// opd.cpp and pupil.cpp, so that the pupil function and its refocus() give W bit for bit as
/// opd_points() (ADR 0032, points 2, 3 and 5).

#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

#include "common.hpp"
#include "rtt/analysis/opd.hpp"
#include "rtt/math/types.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/run_control.hpp"
#include "rtt/trace/sources.hpp"

namespace rtt::analysis::detail {

/// Reference sphere and the data needed to measure OPL up to it (decided for #29; exit pupil at
/// infinity #102; through focus ADR 0032, point 5).
struct Reference {
  ReferenceSphere sphere;  ///< radius +infinity if the exit pupil is at infinity
  /// Crossing of a ray line with the sphere on the side of the exit pupil: +1 if the exit
  /// pupil lies upstream of C along the chief ray (s = +R + ..., the usual case), -1 if it lies
  /// downstream (virtual exit pupil behind the image). Unused for R = infinity.
  double branch = 1.0;
  double lambda_mm = 0.0;  ///< reference wavelength, mm
  /// C: the chief ray of the reference wavelength on the image surface, before focus_shift.
  math::Vec3 chief_pos = math::Vec3::Zero();
  math::Vec3 chief_dir = math::Vec3::UnitZ();  ///< d_c, its unit direction
  std::optional<math::Vec3> exit_pupil;        ///< E; none at infinity
};

/// The reference about C' = chief_pos + focus_shift chief_dir (only for focus_shift != 0, so that
/// 0 keeps C bit for bit) with radius |C' - E| and the branch of E.
/// @throws AnalysisError if C' lies on E
[[nodiscard]] Reference reference_at(const math::Vec3& chief_pos,
                                     const math::Vec3& chief_dir,
                                     const std::optional<math::Vec3>& exit_pupil,
                                     double focus_shift,
                                     double lambda_mm);

/// The reference of `field`: chief ray of the reference wavelength, paraxial exit pupil.
/// @throws as opd_points() (no stop, chief ray lost, exit pupil on the image surface)
[[nodiscard]] Reference make_reference(const CompiledSystem& system,
                                       PathId path,
                                       std::uint16_t field,
                                       trace::Aiming aiming,
                                       double focus_shift);

/// OPL of a ray at p with direction d and OPL `opl` on the image surface up to the reference
/// sphere, up to a constant n' R common to all rays; none if the line misses the sphere.
[[nodiscard]] std::optional<double> opl_to_sphere(
    const math::Vec3& p, const math::Vec3& d, double opl, const Reference& ref, double n_image);

/// n' = |Re n| of the image space of `path` at `wavelength` (isotropic image medium).
[[nodiscard]] double image_index(const CompiledSystem& system,
                                 PathId path,
                                 std::uint16_t wavelength);

/// The chief ray of `field` at `wavelength`, traced to the image surface.
/// @throws AnalysisError if it does not arrive
[[nodiscard]] trace::RayBatch chief_ray(const CompiledSystem& system,
                                        PathId path,
                                        std::uint16_t field,
                                        std::uint16_t wavelength,
                                        trace::Aiming aiming);

/// OPL of the chief ray up to the sphere, the reference of W.
/// @throws AnalysisError if it misses the sphere
[[nodiscard]] double chief_opl(
    const math::Vec3& p, const math::Vec3& d, double opl, const Reference& ref, double n_image);

/// W and status of one ray against `ref` (px, py are not set): W = (OPL_chief - OPL_ray) /
/// lambda_ref if it arrived on the image surface and crosses the sphere, otherwise W = 0 and its
/// status, Vignetted if it was Alive.
[[nodiscard]] OpdPoint point_against(const math::Vec3& p,
                                     const math::Vec3& d,
                                     double opl,
                                     bool arrived_on_image,
                                     trace::RayStatus status,
                                     double opl_chief,
                                     const Reference& ref,
                                     double n_image);

/// OPD of every ray of `sampling` at `wavelength` against `ref`, in waves at the reference
/// wavelength; every ray counted in `losses` with the status of its point.
[[nodiscard]] std::vector<OpdPoint> points_against(const CompiledSystem& system,
                                                   PathId path,
                                                   std::uint16_t field,
                                                   std::uint16_t wavelength,
                                                   const trace::PupilSampling& sampling,
                                                   trace::Aiming aiming,
                                                   const Reference& ref,
                                                   LossCounter& losses,
                                                   const trace::RunControl& control);

}  // namespace rtt::analysis::detail
