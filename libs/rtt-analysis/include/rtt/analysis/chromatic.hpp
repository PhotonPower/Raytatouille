#pragma once

/// @file chromatic.hpp
/// Longitudinal (axial) and lateral (transverse) colour as data objects (docs/architecture.md,
/// Analyse; decided for #31).
///
/// Conventions:
/// - Wavelength pair (first, second): differences are value(first) - value(second). Default:
///   first = the first system wavelength, second = the last one (F - C for an F, d, C list).
/// - Longitudinal colour, paraxial: focus per wavelength from rtt-paraxial (rear focal point for
///   an object at infinity, paraxial image of the axial object point otherwise).
/// - Longitudinal colour, real: axis crossing of the real ray at normalised pupil coordinates
///   (0, zone) of the on-axis field per wavelength.
/// - Both differences are measured along the image-space propagation direction, mm.
/// - Lateral colour: real chief ray (pupil centre) of each wavelength on the image surface
///   (local x, y as in spot.hpp), and its offset from the chief ray of the reference wavelength.

#include <cstdint>
#include <optional>
#include <vector>

#include "rtt/analysis/spot.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/paraxial/seidel.hpp"
#include "rtt/trace/sources.hpp"

namespace rtt::analysis {

/// Options of the longitudinal-colour analysis.
struct ChromaticOptions {
  /// Wavelength pair; std::nullopt: first and last system wavelength.
  std::optional<paraxial::ChromaticPair> pair;
  double zone = 1.0;  ///< normalised pupil height of the real ray, in (0, 1]
  trace::Aiming aiming = trace::Aiming::Real;  ///< aiming of the real rays
};

/// Focus positions of one wavelength.
struct FocusPosition {
  std::uint16_t wavelength = 0;  ///< index into CompiledSystem::wavelengths_um()
  double paraxial_z = 0.0;       ///< global z of the paraxial focus, mm
  double real_z = 0.0;           ///< global z of the axis crossing of the real zone ray, mm
};

/// Longitudinal colour of a path.
struct LongitudinalColour {
  paraxial::ChromaticPair pair;     ///< pair used
  std::vector<FocusPosition> foci;  ///< every system wavelength, in system order
  double paraxial = 0.0;  ///< paraxial focus(first) - focus(second) along the propagation, mm
  double real = 0.0;      ///< real focus(first) - focus(second) along the propagation, mm
};

/// Longitudinal colour of `path`.
/// @throws std::invalid_argument for an invalid path, wavelength pair or zone (and as
///         rtt::trace::aim_ray)
/// @throws rtt::compile::NoStopError (a std::invalid_argument) if the path has no stop; checked
///         by rtt::compile::require_stop before any ray is traced (ADR 0022)
/// @throws rtt::paraxial::ParaxialError if the path is not rotationally symmetric
/// @throws AnalysisError if a paraxial focus does not exist (afocal), a real ray does not reach
///         the image surface or does not cross the axis
[[nodiscard]] LongitudinalColour longitudinal_colour(const compile::CompiledSystem& system,
                                                     compile::PathId path,
                                                     const ChromaticOptions& options = {});

/// Lateral colour of one field.
struct LateralColour {
  std::uint16_t field = 0;     ///< index into CompiledSystem::fields().points
  std::vector<Point2> chief;   ///< chief ray per system wavelength on the image surface, mm
  std::vector<Point2> offset;  ///< chief[i] - chief[reference wavelength], mm
};

/// Lateral colour of `field`.
/// @throws std::invalid_argument for an invalid path or field (and as rtt::trace::aim_ray)
/// @throws rtt::compile::NoStopError (a std::invalid_argument) if the path has no stop; checked
///         by rtt::compile::require_stop before any ray is traced (ADR 0022)
/// @throws rtt::paraxial::ParaxialError if the path is not rotationally symmetric
/// @throws AnalysisError if a chief ray does not reach the image surface
[[nodiscard]] LateralColour lateral_colour(const compile::CompiledSystem& system,
                                           compile::PathId path,
                                           std::uint16_t field,
                                           trace::Aiming aiming = trace::Aiming::Real);

/// longitudinal_colour() and lateral_colour(), overloads with cancellation and progress (#83): the
/// same analysis with a trace::RunControl (rtt/trace/run_control.hpp) as last parameter.
/// Cancellation and progress are run-time control, never part of the options or of the result: with
/// an empty control the result is the same as with the overload without it, bit for bit, and it
/// never depends on the control. A cancelled call throws trace::Cancelled, an exception of the
/// progress callback is rethrown; both only after the parallel parts. Progress stages: "aim" and
/// "trace" for every ray bundle (e.g. once per wavelength of a polychromatic spot, so a stage name
/// can repeat), "field" for the points of a field sweep and "wavelength" for the colour analyses.
/// @throws as the overloads without a control, trace::Cancelled, or the callback's exception
/// @{
[[nodiscard]] LongitudinalColour longitudinal_colour(const compile::CompiledSystem& system,
                                                     compile::PathId path,
                                                     const ChromaticOptions& options,
                                                     const trace::RunControl& control);
[[nodiscard]] LateralColour lateral_colour(const compile::CompiledSystem& system,
                                           compile::PathId path,
                                           std::uint16_t field,
                                           trace::Aiming aiming,
                                           const trace::RunControl& control);
/// @}

}  // namespace rtt::analysis
