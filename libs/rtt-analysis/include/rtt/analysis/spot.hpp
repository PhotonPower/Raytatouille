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

#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <stdexcept>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/sources.hpp"

namespace rtt::analysis {

/// Thrown when an analysis has no defined result, e.g. the chief ray does not reach the image
/// surface or no ray arrives.
class AnalysisError : public std::runtime_error {
 public:
  using std::runtime_error::runtime_error;
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
};

/// Spot diagram of `field`. Rays come from rtt::trace::make_rays and are traced with the
/// SequentialTracer; all sums run serially in trace order, so the result is bitwise independent
/// of the number of threads.
/// @param wavelength index of one system wavelength, or std::nullopt for all wavelengths with
///                   the model's wavelength weights
/// @throws std::invalid_argument for an invalid path, field or wavelength, a wavelength weight
///         that is negative or not finite, or weights that do not sum to a positive value (and
///         as rtt::trace::make_rays)
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
};

/// Ray fans of `field` at `wavelength`, relative to the chief ray (see file comment).
/// @throws std::invalid_argument for an invalid path, field, wavelength or points < 1 (and as
///         rtt::trace::make_rays)
/// @throws rtt::paraxial::ParaxialError if the path is not rotationally symmetric (aiming)
/// @throws AnalysisError if the chief ray does not reach the image surface
[[nodiscard]] RayFan ray_fan(const compile::CompiledSystem& system,
                             compile::PathId path,
                             std::uint16_t field,
                             std::uint16_t wavelength,
                             const FanOptions& options = {});

}  // namespace rtt::analysis
