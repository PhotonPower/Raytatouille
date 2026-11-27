#pragma once

/// @file source.hpp
/// Sources of the non-sequential tracer (ADR 0057, point 6, Nachtrag Review PR 1, c and d,
/// Nachtrag PR 4, e-g): point, collimated and Lambert emitters, sampled at random with a seed or
/// on a regular grid. A start ray depends only on the source, its index and the seed, so trace()
/// creates it in its block instead of storing all start rays.
///
/// Conventions: positions in mm, directions as unit vectors; local source coordinates with the
/// axis +z, mapped to global coordinates by Source::pose; angles in rad; power in W.

#include <cstddef>
#include <cstdint>
#include <numbers>
#include <variant>
#include <vector>

#include "rtt/math/isometry.hpp"
#include "rtt/math/types.hpp"

namespace rtt::nonseq {

/// Emitting disc of radius `radius` about the local origin in the plane z = 0, mm.
struct Disc {
  double radius = 0.0;  ///< mm, finite and > 0
};

/// Emitting rectangle |x| <= half_x, |y| <= half_y in the plane z = 0, mm.
struct Rectangle {
  double half_x = 0.0;  ///< mm, finite and > 0
  double half_y = 0.0;  ///< mm, finite and > 0
};

/// Shape of an extended emitter.
using EmitterShape = std::variant<Disc, Rectangle>;

/// Point source at the local origin, directions uniform in the cone theta <= max_angle about +z
/// (PBRT 4ed, §A.5.4); intensity I = power / (2 pi (1 - cos max_angle)).
struct PointEmitter {
  double max_angle = std::numbers::pi;  ///< cone half-angle, rad, in (0, pi]
};

/// Collimated source: direction +z, positions uniform over the shape (disc: PBRT 4ed, §A.5.1).
struct CollimatedEmitter {
  EmitterShape shape = Disc{1.0};
};

/// Lambert source: positions as CollimatedEmitter, directions cosine-weighted (Malley, PBRT 4ed,
/// §A.5.3) and limited to theta <= max_angle; radiance L = power / (pi A sin^2 max_angle)
/// (NBS TN 910-1, Eq. (A2-14)), A the area of the shape.
struct LambertEmitter {
  EmitterShape shape = Disc{1.0};
  double max_angle = std::numbers::pi / 2.0;  ///< rad, in (0, pi/2]
};

/// How the unit square (xi_1, xi_2) (and (xi_3, xi_4) for an extended Lambert source) is sampled.
enum class Sampling : std::uint8_t {
  Random,  ///< Philox-4x32-10 (ADR 0057, point 5; Nachtrag PR 4, f)
  Grid,    ///< cell centres of a regular grid (ADR 0057, Nachtrag Review PR 1, d)
};

/// A source (ADR 0057, point 6). Each start ray is unpolarized and carries Phi_lambda / M, with
/// Phi_lambda = power * g_lambda / sum(g) over `wavelengths` (g the system wavelength weights)
/// and M = source_ray_count(); every wavelength runs with the same M start rays.
struct Source {
  std::variant<PointEmitter, CollimatedEmitter, LambertEmitter> emitter = PointEmitter{};
  math::Isometry3 pose;    ///< local source coordinates -> global coordinates
  double power = 0.0;      ///< W, finite and > 0
  std::uint64_t rays = 0;  ///< N, start rays per wavelength; > 0 and <= 2^40
  Sampling sampling = Sampling::Random;
  /// Indices into CompiledSystem::wavelengths_um(), distinct; empty = all.
  std::vector<std::uint16_t> wavelengths;
};

/// Position and direction of a start ray, global, mm and unit vector.
struct SourceRay {
  math::Vec3 position = math::Vec3::Zero();
  math::Vec3 direction = math::Vec3::UnitZ();
};

/// Number M of start rays per wavelength (ADR 0057, Nachtrag Review PR 1, d): N for Random; n^2
/// with n = ceil(sqrt(N)) on a grid of two dimensions (point and collimated sources); n^4 with
/// n = ceil(N^(1/4)) for a Lambert source on a grid (n^2 positions times n^2 directions).
[[nodiscard]] std::uint64_t source_ray_count(const Source& source) noexcept;

/// Start ray `index` (0 <= index < source_ray_count()) of source number `source_number`
/// (ADR 0057, Nachtrag PR 4, f). Random: xi_1, xi_2 = uniform53 of words 0, 1 and 2, 3 of
/// Philox-4x32-10 with key = seed (low, high 32 bits) and counter (source_number 2^16, index low,
/// index high, 0); xi_3, xi_4 the same with draw 1. Grid: the cell centres (i + 1/2)/n, with
/// index = i_1 + n i_2 on a grid of two dimensions and index = p n^2 + q (position cell p,
/// direction cell q, each as before) for a Lambert source.
/// - point: 1 - cos theta = 2 xi_1 sin^2(max_angle / 2), phi = 2 pi xi_2 (PBRT 4ed, §A.5.4);
/// - disc: r = R sqrt(xi_1), phi = 2 pi xi_2 (§A.5.1); rectangle: x = (2 xi_1 - 1) half_x,
///   y = (2 xi_2 - 1) half_y;
/// - Lambert direction: the point (rho cos phi, rho sin phi) of the disc of radius sin max_angle,
///   rho = sin(max_angle) sqrt(xi_3), phi = 2 pi xi_4, lifted to z = sqrt(1 - rho^2) (Malley,
///   §A.5.3).
/// The result is mapped to global coordinates with Source::pose.
[[nodiscard]] SourceRay source_ray(const Source& source,
                                   std::uint32_t source_number,
                                   std::uint64_t index,
                                   std::uint64_t seed) noexcept;

/// Checks a source against the rules of the struct comments.
/// @param wavelengths number of system wavelengths
/// @throws std::invalid_argument naming the broken rule
void check_source(const Source& source, std::size_t wavelengths);

}  // namespace rtt::nonseq
