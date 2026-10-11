#pragma once

/// @file pupil.hpp
/// Pupil function and its cells (ADR 0032; #202).
///
/// Conventions:
/// - Normalised pupil coordinates p = (px, py) as in rtt-trace (sources.hpp): the unit circle is
///   the rim of the paraxial entrance pupil, +y meridional.
/// - The cells divide the square [-1, 1]^2 into N x N squares of edge a = 2/N, in the order of
///   trace::GridPupil: rows from py = -1 upwards, px fastest; cell (i, k), i along x and k along
///   y, has index k N + i and the centre (-1 + (i + 1/2) a, -1 + (k + 1/2) a).
/// - W in waves at the reference wavelength as in opd.hpp (W = OPL_ref - OPL, reference sphere
///   about the chief ray of the reference wavelength, the chief ray of the evaluated wavelength
///   as reference ray). The phase of a pupil point relative to the chief ray at the evaluated
///   wavelength lambda is -2 pi W lambda_ref / lambda (ADR 0032, "Bezeichnungen": fields
///   ~ exp(i (k . r - omega t)), so the phase grows with the OPL).
/// - Directions in the image space: s = n' (d . e_x, d . e_y), the direction sines in the plane
///   perpendicular to the chief ray of the reference wavelength (ADR 0032, point 4).

#include <array>
#include <complex>
#include <cstdint>
#include <optional>
#include <vector>

#include "rtt/analysis/opd.hpp"
#include "rtt/analysis/spot.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/math/types.hpp"
#include "rtt/model/validate.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/run_control.hpp"
#include "rtt/trace/sources.hpp"

namespace rtt::analysis {

/// Largest grid N of the pupil function. Memory bound: each traced cell is a ray of a
/// trace::RayBatch (13 doubles, 9 complex PRT elements, ids and status: about 260 bytes) plus up
/// to about 390 bytes in the pupil function (cells, end states, PRT copy), so N = 2048 needs
/// about 0.79 N^2 = 3.3 million traced cells (pi/4 of N^2) or 2.7 GB at most. A typo such as
/// 20480 would ask for 100 times that; it is an input error instead. The acceptance of M6 needs
/// N <= 256.
inline constexpr int kMaxPupilGrid = 2048;

/// Cells of the pupil function (ADR 0032, point 1), N x N entries each in the order above.
struct PupilCells {
  int grid = 0;  ///< N
  /// Area fraction f in [0, 1] of the cell inside the unit circle, in closed form: exactly 1 for
  /// a cell inside the circle, exactly 0 for a cell outside.
  std::vector<double> fraction;
  /// Sample point: the centre of a full cell (f = 1) and of an empty one (f = 0, not traced);
  /// the centroid of the part inside the circle for a boundary cell (0 < f < 1), which lies
  /// strictly inside the circle.
  std::vector<double> px;
  std::vector<double> py;
};

/// Cells of an N x N grid, N = `grid`.
/// @throws std::invalid_argument for grid < 2 or grid > kMaxPupilGrid
[[nodiscard]] PupilCells pupil_cells(int grid);

/// Options of the pupil function (ADR 0032, point 3).
struct PupilFunctionOptions {
  int grid = 64;  ///< N in [2, kMaxPupilGrid]
  /// Gaussian apodisation factor G >= 0 (finite): A_apod = exp(-G rho^2) with rho = |p| of the
  /// sample point, G = 1 gives 1/e at the rim of the entrance pupil (Ansys Zemax OpticStudio User
  /// Guide, Release 2025 R1, "Apodization Type"). Analysis option only, not in the file format.
  double apodization = 0.0;
  /// Through focus (ADR 0032, point 5): the reference sphere about C' = C + focus_shift d_c, mm,
  /// as OpdOptions::focus_shift. Finite.
  double focus_shift = 0.0;
  trace::Aiming aiming = trace::Aiming::Real;  ///< aiming of all rays, chief rays included
  bool with_prt = false;  ///< copy the accumulated PRT matrix P of each cell (ADR 0034)
  /// Warning "rays.lost" if more than this fraction of the traced cells is lost; in [0, 1].
  double lost_warning_fraction = 0.5;
  /// Warning "pupil.mapping_nonlinear" if delta_max of the mapping exceeds this; in [0, 1].
  double mapping_warning = 1e-2;
};

/// Affine mapping of the pupil coordinates to the direction sines, s ~ s0 + M p (ADR 0032,
/// point 4), at the evaluated wavelength. M comes from central differences at the four parabasal
/// rays p = (+-1e-3, 0), (0, +-1e-3) of the same wavelength and aiming (error O(1e-6)).
struct PupilMapping {
  /// M row-major: rows s_x, s_y, columns px, py; dimensionless (direction sines per pupil unit).
  std::array<double, 4> m{};
  /// s of the chief ray of the evaluated wavelength; exactly 0 at the reference wavelength (by
  /// construction: the same ray as d_c).
  std::array<double, 2> s0{};
  /// Origin of the PSF plane: C' (the centre of the reference sphere after focus_shift), global,
  /// mm.
  math::Vec3 centre = math::Vec3::Zero();
  /// e_x: the local x axis of the image surface made perpendicular to d_c and normalised, global.
  math::Vec3 e_x = math::Vec3::UnitX();
  math::Vec3 e_y = math::Vec3::UnitY();  ///< d_c x e_x, global
  /// Unit direction of the chief ray of the reference wavelength on the image surface, global.
  math::Vec3 d_c = math::Vec3::UnitZ();
  /// max |s - s0 - M p| / ||M||_2 over the arrived cells with f > 0, pupil units (||M||_2 the
  /// largest singular value).
  double delta_max = 0.0;
  /// The same, root mean square weighted with the area fractions f.
  double delta_rms = 0.0;
};

/// End states of the traced cells on the image surface, kept for refocus() (an addition to
/// ADR 0032 for ADR 0033, through focus without tracing again; #202). Per cell in cell order;
/// cells with f = 0 are not traced and hold NaN, the status Vignetted and kNoSurface.
struct PupilEndStates {
  std::vector<double> pos_x, pos_y, pos_z;    ///< global, mm
  std::vector<double> dir_x, dir_y, dir_z;    ///< global unit direction
  std::vector<double> opl;                    ///< mm
  std::vector<double> weight;                 ///< power of the ray (ADR 0021)
  std::vector<trace::RayStatus> ray_status;   ///< status of the traced ray
  std::vector<std::uint32_t> last_surface;    ///< RayBatch::last_surface
  math::Vec3 chief_pos = math::Vec3::Zero();  ///< chief ray of the evaluated wavelength
  math::Vec3 chief_dir = math::Vec3::UnitZ();
  double chief_opl = 0.0;
  math::Vec3 reference_pos = math::Vec3::Zero();   ///< C: chief ray of lambda_ref, focus_shift 0
  math::Vec3 reference_dir = math::Vec3::UnitZ();  ///< d_c
  std::optional<math::Vec3>
      exit_pupil;                   ///< E, the centre of the paraxial exit pupil; none at infinity
  double image_index = 1.0;         ///< n' = |Re n| of the image space at the evaluated wavelength
  std::uint32_t image_surface = 0;  ///< index into CompiledSystem::surfaces()
};

/// Pupil function of one field at one wavelength on N x N cells (ADR 0032, point 3).
///
/// **Who counts or evaluates cells checks f > 0 first, then status**: cells outside the circle
/// (f = 0) are not traced and carry status Vignetted with A = W = 0; they do not count in
/// losses. A traced cell is Alive if its ray arrived on the image surface and crossed the
/// reference sphere, otherwise it keeps the status of its ray (Vignetted if it arrived but
/// misses the sphere) with A = W = 0, as in opd_points().
struct PupilFunction {
  compile::PathId path{};
  std::uint16_t field = 0;               ///< index into CompiledSystem::fields().points
  std::uint16_t wavelength = 0;          ///< index into CompiledSystem::wavelengths_um()
  double wavelength_um = 0.0;            ///< evaluated wavelength in vacuum, um
  double reference_wavelength_um = 0.0;  ///< reference wavelength in vacuum, um (unit of W)
  int grid = 0;                          ///< N
  double focus_shift = 0.0;              ///< mm (PupilFunctionOptions::focus_shift)
  double apodization = 0.0;              ///< G
  double lost_warning_fraction = 0.5;    ///< threshold of "rays.lost" (options, for refocus)
  double mapping_warning = 1e-2;         ///< threshold of "pupil.mapping_nonlinear" (options)
  ReferenceSphere sphere;                ///< after focus_shift
  // Per cell, N * N entries in the order of PupilCells:
  std::vector<double> fraction;  ///< f
  std::vector<double> px, py;    ///< sample point (centre or centroid)
  /// sqrt(weight) exp(-G rho^2); 0 if the cell did not arrive or f = 0.
  std::vector<double> amplitude;
  /// W = OPL_ref - OPL in waves at lambda_ref, > 0 where the wave leads (opd.hpp); 0 if the
  /// cell did not arrive or f = 0.
  std::vector<double> w;
  std::vector<trace::RayStatus> status;
  /// With PupilFunctionOptions::with_prt: 9 elements per cell, global, row-major (3 row + col)
  /// as RayBatch::prt; zero for cells with f = 0. Empty otherwise.
  std::vector<std::complex<double>> prt;
  /// None for an afocal image space (|det M| <= 1e-12 n'^2); analyses that need the mapping
  /// reject it (ADR 0033).
  std::optional<PupilMapping> mapping;
  /// sum f over the arrived cells / sum f over all cells.
  double arrived_fraction = 0.0;
  RayLosses losses;  ///< the traced cells (f > 0) by final status (ADR 0023)
  /// Warnings with stable codes: "rays.lost", "stop.clips_beam" as opd_points(),
  /// "pupil.mapping_nonlinear" (location /paths/i).
  std::vector<model::Diagnostic> warnings;
  PupilEndStates end;  ///< for refocus()
};

/// Pupil function of `field` at `wavelength` (ADR 0032, point 3). The traced cells are the
/// cells with f > 0 as a trace::PointListPupil, so W and status are bit for bit those of
/// opd_points() with that list and the same options; bit for bit at every thread count.
/// @throws std::invalid_argument for an invalid path or wavelength, options out of range (grid,
///         apodization G < 0 or not finite, focus_shift not finite, lost_warning_fraction or
///         mapping_warning outside [0, 1]) and as trace::make_rays (e.g. an invalid field)
/// @throws AnalysisError with location /paths/i for a path with an ordinary or extraordinary
///         event (crystal modes, #35), before any ray is traced; with the location of the image
///         surface if the chief ray grazes it (|x_B x d_c| <= 1e-6, no frame of the PSF plane);
///         and as opd_points() (chief ray lost or missing the reference sphere)
/// @throws rtt::paraxial::ParaxialError, rtt::compile::NoStopError as opd_points()
[[nodiscard]] PupilFunction pupil_function(const compile::CompiledSystem& system,
                                           compile::PathId path,
                                           std::uint16_t field,
                                           std::uint16_t wavelength,
                                           const PupilFunctionOptions& options = {});

/// pupil_function() with cancellation and progress (#83), stages "aim" and "trace"; with an
/// empty control the same result bit for bit.
[[nodiscard]] PupilFunction pupil_function(const compile::CompiledSystem& system,
                                           compile::PathId path,
                                           std::uint16_t field,
                                           std::uint16_t wavelength,
                                           const PupilFunctionOptions& options,
                                           const trace::RunControl& control);

/// The pupil function `pupil` for another focus_shift, from its end states without tracing
/// again (ADR 0032, point 5; ADR 0033, through focus): W, status, amplitude, sphere,
/// arrived_fraction, losses, warnings, mapping->centre and the deltas of the mapping are
/// recomputed. Bit for bit the same as pupil_function() with that focus_shift. `system` must be
/// the system that produced `pupil` (for the surface names and pointers of the warnings).
/// @throws std::invalid_argument if focus_shift is not finite or `system` cannot be the system of
///         `pupil` (path, field or wavelength out of range, another image surface, or a last
///         surface of a ray outside its surfaces)
/// @throws AnalysisError as pupil_function() for the shifted sphere
[[nodiscard]] PupilFunction refocus(const compile::CompiledSystem& system,
                                    const PupilFunction& pupil,
                                    double focus_shift);

}  // namespace rtt::analysis
