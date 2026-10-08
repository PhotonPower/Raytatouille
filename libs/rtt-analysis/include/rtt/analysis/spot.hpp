#pragma once

/// @file spot.hpp
/// Spot diagram and ray fans as data objects, never plots (docs/architecture.md, Analyse).
///
/// Conventions (decided for #28):
/// - Image surface: the surface of the last event of the path. Spot and fan coordinates x, y
///   are in mm in the local coordinate system of that surface, with its axes as given by its
///   pose (for a tilted or curved image surface: the local x, y of the hit point). A tilted or
///   decentred image surface currently makes the path non rotationally symmetric, so the ray
///   aiming of #8 (rtt::paraxial::first_order) throws rtt::paraxial::ParaxialError.
/// - Reference ray ("chief ray"): the ray at normalised pupil coordinates (0, 0), aimed with the
///   chosen Aiming at the reference wavelength of the system. It is the reference for every
///   spot (monochromatic and polychromatic) and every fan, so lateral colour stays visible.
/// - Weights: polychromatic spots use the wavelength weights of the model normalised to a sum of
///   1, times the power weight of each ray. Field weights (model::Field::weight) play no role
///   here, because every spot belongs to one field; they matter for multi-field merit functions.
/// - A ray counts as arrived if its status is Alive after the trace and its last surface is the
///   image surface; every other ray (vignetted, missed, TIR, aiming failed, ...) counts as
///   vignetted.

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/run_control.hpp"
#include "rtt/trace/sequential.hpp"
#include "rtt/trace/sources.hpp"

namespace rtt::analysis {

/// Thrown when an analysis has no defined result, e.g. the chief ray does not reach the image
/// surface or no ray arrives. For a single lost ray (chief or zone ray) the accessors give the
/// place and the cause (ADR 0022); otherwise they are empty.
class AnalysisError : public std::runtime_error {
 public:
  /// The lost ray of an AnalysisError.
  struct LostRay {
    /// Last surface the ray reached (RayBatch::last_surface): for Vignetted, Absorbed, Tir and
    /// EventImpossible the surface where it stopped; for Missed and NoConvergence the surface
    /// before the one it did not reach. None if the ray reached no surface.
    std::optional<model::SurfaceId> surface;
    std::optional<std::string> location;  ///< JSON pointer of that surface in the system file
    trace::RayStatus ray_status = trace::RayStatus::Alive;
    std::optional<std::uint16_t> field;  ///< index into CompiledSystem::fields().points
    std::uint16_t wavelength = 0;        ///< index into CompiledSystem::wavelengths_um()
  };

  using std::runtime_error::runtime_error;

  AnalysisError(const std::string& message, LostRay ray)
      : std::runtime_error(message), ray_(std::move(ray)) {}

  [[nodiscard]] std::optional<model::SurfaceId> surface() const {
    return ray_ ? ray_->surface : std::nullopt;
  }
  [[nodiscard]] std::optional<std::string> location() const {
    return ray_ ? ray_->location : std::nullopt;
  }
  /// Status of the lost ray; Alive if it ended alive on another surface than the image.
  [[nodiscard]] std::optional<trace::RayStatus> ray_status() const {
    return ray_ ? std::optional(ray_->ray_status) : std::nullopt;
  }
  [[nodiscard]] std::optional<std::uint16_t> field() const {
    return ray_ ? ray_->field : std::nullopt;
  }
  [[nodiscard]] std::optional<std::uint16_t> wavelength() const {
    return ray_ ? std::optional(ray_->wavelength) : std::nullopt;
  }

 private:
  std::optional<LostRay> ray_;
};

/// Where the rays of an analysis went (ADR 0023): the rays of its pupil sampling (not the
/// chief ray) by final status, and the surface where most of the lost rays ended.
struct RayLosses {
  std::size_t launched = 0;  ///< rays traced
  /// Rays per final RayStatus, indexed by its value. Alive counts the rays that arrived at the
  /// image surface; a ray that ended Alive elsewhere counts as Vignetted.
  std::array<std::size_t, trace::kRayStatusCount> by_status{};
  /// Index into CompiledSystem::surfaces() of the surface that most lost rays reached last
  /// (RayBatch::last_surface: the surface where a Vignetted, Absorbed, Tir or EventImpossible
  /// ray stopped, the surface before the one a Missed or NoConvergence ray did not reach);
  /// ties go to the lowest index. None if no lost ray reached a surface.
  std::optional<std::uint32_t> worst_surface;
  std::size_t worst_surface_count = 0;  ///< lost rays whose last surface is worst_surface

  /// Number of rays that ended with status `s`.
  [[nodiscard]] std::size_t count(trace::RayStatus s) const noexcept {
    return by_status[static_cast<std::size_t>(s)];
  }
};

/// Point in the local x, y plane of the image surface, mm.
struct Point2 {
  double x = 0.0;
  double y = 0.0;
};

/// One arrived ray of a spot diagram.
struct SpotPoint {
  double x = 0.0;                ///< local x on the image surface, mm
  double y = 0.0;                ///< local y on the image surface, mm
  std::uint16_t wavelength = 0;  ///< index into CompiledSystem::wavelengths_um()
  double weight = 0.0;           ///< statistical weight (see file comment), dimensionless
};

/// Statistics of a weighted point set.
struct SpotStatistics {
  Point2 centroid;            ///< c = sum w r / sum w, mm
  double rms_centroid = 0.0;  ///< sqrt(sum w |r - c|^2 / sum w), mm
  double rms_chief = 0.0;     ///< sqrt(sum w |r - chief|^2 / sum w), mm
  double geo_centroid = 0.0;  ///< max |r - c| over the points with weight > 0, mm
  double geo_chief = 0.0;     ///< max |r - chief| over the points with weight > 0, mm
};

/// Weighted statistics of `points` about their centroid and about `chief`.
/// @throws AnalysisError if the points are empty or their weights do not sum to a positive value
[[nodiscard]] SpotStatistics spot_statistics(std::span<const SpotPoint> points, Point2 chief);

/// Options of a spot diagram.
struct SpotOptions {
  trace::PupilSampling sampling = trace::HexapolarPupil{6};  ///< pupil sampling per wavelength
  trace::Aiming aiming = trace::Aiming::Real;  ///< aiming of all rays, chief ray included
  /// Warning "rays.lost" if more than this fraction of the launched rays is lost (ADR 0023);
  /// in [0, 1]. Vignetting at the field edge is intended, hence the default of one half.
  double lost_warning_fraction = 0.5;
};

/// Spot diagram of one field at one wavelength or polychromatic.
struct SpotDiagram {
  std::uint16_t field = 0;                  ///< index into CompiledSystem::fields().points
  std::optional<std::uint16_t> wavelength;  ///< none: polychromatic over all system wavelengths
  std::uint32_t image_surface = 0;          ///< index into CompiledSystem::surfaces()
  std::vector<SpotPoint> points;            ///< arrived rays, in trace order
  Point2 chief;                             ///< reference ray on the image surface, mm
  SpotStatistics stats;             ///< statistics of `points` about their centroid and `chief`
  std::size_t rays_launched = 0;    ///< rays started, over all wavelengths of the spot
  std::size_t rays_arrived = 0;     ///< rays that reached the image surface (= points.size())
  double vignetted_fraction = 0.0;  ///< (launched - arrived) / launched, unweighted
  RayLosses losses;  ///< rays of the sampling by final status, worst loss surface (ADR 0023)
  /// Warnings with stable codes (ADR 0022, 0023): "rays.lost" above lost_warning_fraction,
  /// "stop.clips_beam" if rays end Vignetted at the stop surface.
  std::vector<model::Diagnostic> warnings;
};

/// Spot diagram of `field`. Rays come from rtt::trace::make_rays and are traced with the
/// SequentialTracer; all sums run serially in trace order, so the result is bitwise independent
/// of the number of threads.
/// @param wavelength index of one system wavelength, or std::nullopt for all wavelengths with
///                   the model's wavelength weights
/// @throws std::invalid_argument for an invalid path, field or wavelength, a wavelength weight
///         that is negative or not finite, weights that do not sum to a positive value, or
///         options.lost_warning_fraction outside [0, 1] (and as rtt::trace::make_rays)
/// @throws rtt::compile::NoStopError (a std::invalid_argument) if the path has no stop; checked
///         by rtt::compile::require_stop before any ray is traced (ADR 0022)
/// @throws rtt::paraxial::ParaxialError if the path is not rotationally symmetric (aiming)
/// @throws AnalysisError if the chief ray does not reach the image surface, or no ray with a
///         positive weight arrives (e.g. all rays of the weighted wavelengths are vignetted)
[[nodiscard]] SpotDiagram spot(const compile::CompiledSystem& system,
                               compile::PathId path,
                               std::uint16_t field,
                               std::optional<std::uint16_t> wavelength,
                               const SpotOptions& options = {});

/// Options of ray fans.
struct FanOptions {
  int points = 21;                             ///< points per fan, evenly spaced on [-1, 1]
  trace::Aiming aiming = trace::Aiming::Real;  ///< aiming of all rays, chief ray included
  /// Warning "rays.lost" if more than this fraction of the launched rays is lost (ADR 0023);
  /// in [0, 1]. Vignetting at the field edge is intended, hence the default of one half.
  double lost_warning_fraction = 0.5;
};

/// One point of a ray fan: transverse aberration relative to the chief ray.
struct FanPoint {
  double p = 0.0;   ///< normalised pupil coordinate (py for tangential, px for sagittal)
  double ex = 0.0;  ///< x - chief.x on the image surface, mm; 0 if the ray did not arrive
  double ey = 0.0;  ///< y - chief.y on the image surface, mm; 0 if the ray did not arrive
  /// Alive if the ray arrived; otherwise its status (Vignetted if it stopped on another
  /// surface while Alive). Consumers must check it before using ex, ey.
  trace::RayStatus status = trace::RayStatus::Alive;
};

/// Tangential (px = 0, py in [-1, 1]) and sagittal (py = 0, px in [-1, 1]) ray fans.
struct RayFan {
  std::uint16_t field = 0;           ///< index into CompiledSystem::fields().points
  std::uint16_t wavelength = 0;      ///< index into CompiledSystem::wavelengths_um()
  std::uint32_t image_surface = 0;   ///< index into CompiledSystem::surfaces()
  Point2 chief;                      ///< reference ray on the image surface, mm
  std::vector<FanPoint> tangential;  ///< epsilon_y(py) is the tangential aberration
  std::vector<FanPoint> sagittal;    ///< epsilon_x(px) is the sagittal aberration
  RayLosses losses;  ///< rays of the sampling by final status, worst loss surface (ADR 0023)
  /// Warnings with stable codes (ADR 0022, 0023): "rays.lost" above lost_warning_fraction,
  /// "stop.clips_beam" if rays end Vignetted at the stop surface.
  std::vector<model::Diagnostic> warnings;
};

/// Ray fans of `field` at `wavelength`, relative to the chief ray (see file comment).
/// @throws std::invalid_argument for an invalid path, field, wavelength, points < 1 or
///         options.lost_warning_fraction outside [0, 1] (and as rtt::trace::make_rays)
/// @throws rtt::compile::NoStopError (a std::invalid_argument) if the path has no stop; checked
///         by rtt::compile::require_stop before any ray is traced (ADR 0022)
/// @throws rtt::paraxial::ParaxialError if the path is not rotationally symmetric (aiming)
/// @throws AnalysisError if the chief ray does not reach the image surface
[[nodiscard]] RayFan ray_fan(const compile::CompiledSystem& system,
                             compile::PathId path,
                             std::uint16_t field,
                             std::uint16_t wavelength,
                             const FanOptions& options = {});

/// Overloads with cancellation and progress (#83): the same analysis with a trace::RunControl
/// (rtt/trace/run_control.hpp) as last parameter. Cancellation and progress are run-time
/// control, never part of the options or of the result: with an empty control the result is the
/// same as with the overload without it, bit for bit, and it never depends on the control.
/// A cancelled call throws trace::Cancelled, an exception of the progress callback is rethrown;
/// both only after the parallel parts. Progress stages: "aim" and "trace" for every ray bundle
/// (e.g. once per wavelength of a polychromatic spot, so a stage name can repeat), "field" for
/// the points of a field sweep and "wavelength" for the colour analyses.
/// @throws as the overloads without a control, trace::Cancelled, or the callback's exception
/// @{
[[nodiscard]] SpotDiagram spot(const compile::CompiledSystem& system,
                               compile::PathId path,
                               std::uint16_t field,
                               std::optional<std::uint16_t> wavelength,
                               const SpotOptions& options,
                               const trace::RunControl& control);
[[nodiscard]] RayFan ray_fan(const compile::CompiledSystem& system,
                             compile::PathId path,
                             std::uint16_t field,
                             std::uint16_t wavelength,
                             const FanOptions& options,
                             const trace::RunControl& control);
/// @}

}  // namespace rtt::analysis
