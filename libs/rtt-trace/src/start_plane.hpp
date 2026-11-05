#pragma once

// Internal to rtt-trace: the lowest z a bundle can reach on one surface, for the start plane of
// the aiming (sources.cpp). A template over the concrete shape so that the tests reach it with a
// shape whose vertex is not at the origin (ADR 0041, point 6). Not installed.

#include <algorithm>
#include <cmath>
#include <numbers>
#include <optional>

#include "rtt/math/isometry.hpp"
#include "rtt/math/types.hpp"

namespace rtt::trace::detail {

/// Global z of the vertex of a surface: the local point (0, 0, sag(0, 0)), not the origin
/// (ADR 0041, point 6). For sag(0, 0) == 0 it is the z of the origin bit for bit, also its sign.
/// @param to_global local surface coordinates -> global coordinates (mm)
/// @param sag0      sag(0, 0) of the shape in mm (the `value` of its vertex expansion)
[[nodiscard]] inline double vertex_global_z(const math::Isometry3& to_global, double sag0) {
  if (sag0 == 0.0) return to_global.translation().z();
  return to_global.apply_point(math::Vec3(0.0, 0.0, sag0)).z();
}

/// Lowest global z of the part of a surface that the bundle of an infinite-object field can
/// reach. The bundle is bounded by a cylinder of radius r_bundle around the chief ray through
/// the EP centre (slope tan_field against the axis), so at global z it is at most
/// r_bundle + |z - z_ep| tan_field off the axis. The surface is sampled in local coordinates up
/// to that radius, capped by `cap` (its aperture and the domain of its shape); since a surface
/// that curves back upstream meets the bundle where it is wider, the radius is re-evaluated at
/// the lowest z found until it no longer drops (decided for #8, start plane). The sampling starts
/// at the vertex (vertex_global_z()).
/// Valid for the rotationally symmetric paths that rtt::paraxial::first_order accepts (the
/// sampled disk is centred on the vertex). The bound drops monotonically. With a finite cap
/// (outer radius of the aperture, or shape domain) it is bounded below by the lowest z of the
/// whole cap disk (up to the sampling, hence the final std::min), but it may converge only
/// linearly (factor tan_field |dsag/dr| close to 1); if it has not settled after kMaxRounds and
/// the cap is below kDivergentRadius, the lowest z of the whole cap disk is returned, a safe
/// lower bound since the bundle cannot meet the surface outside its cap (decided for #72).
/// Without such a cap the bound can drop without limit (an unbounded surface without aperture,
/// or with an aperture of 1 km or more, that curves back against a steep field, e.g. a concave
/// paraboloid); then there is no safe start plane and the result is empty.
/// @param shape     the shape (sag in mm over local x, y in mm; NaN outside its domain)
/// @param to_global local surface coordinates -> global coordinates
/// @param cap       largest local radius the bundle can meet, mm (infinity if unbounded)
/// @param z_ep      global z of the entrance pupil, mm
/// @param r_bundle  bundle radius at the entrance pupil, mm
/// @param tan_field tangent of the field angle against the global z axis
/// @return the lowest global z in mm, or std::nullopt if there is no safe start plane
template <class S>
[[nodiscard]] std::optional<double> lowest_z(const S& shape,
                                             const math::Isometry3& to_global,
                                             double cap,
                                             double z_ep,
                                             double r_bundle,
                                             double tan_field) {
  // 2 pi rad; not named kTwoPi, which would hide that of sources.cpp (MSVC C4459).
  constexpr double kFullTurn = 2.0 * std::numbers::pi;
  constexpr int kRadialSamples = 64;
  constexpr int kAzimuthSamples = 8;
  constexpr int kMaxRounds = 50;
  constexpr double kDivergentRadius = 1e6;  // mm
  const double z_vertex = vertex_global_z(to_global, shape.sag(0.0, 0.0));
  // Lowest global z of the surface on the local disk of radius rho around the vertex.
  const auto lowest_on_disk = [&](double rho) {
    double low = z_vertex;
    for (int i = 1; i <= kRadialSamples; ++i) {
      const double r = rho * i / kRadialSamples;
      for (int j = 0; j < kAzimuthSamples; ++j) {
        const double phi = kFullTurn * j / kAzimuthSamples;
        const double x = r * std::sin(phi);
        const double y = r * std::cos(phi);
        // NaN only at the rim of the domain (rounding); skipped.
        const double h = shape.sag(x, y);
        if (!std::isfinite(h)) continue;
        low = std::min(low, to_global.apply_point(math::Vec3(x, y, h)).z());
      }
    }
    return low;
  };
  double z = z_vertex;
  for (int round = 0; round < kMaxRounds; ++round) {
    const double rho = std::min(cap, r_bundle + std::abs(z - z_ep) * tan_field);
    // A bundle radius beyond any optical size means the bound diverges (overflow guard).
    if (!(rho < kDivergentRadius)) break;
    const double low = lowest_on_disk(rho);
    if (low >= z - 1e-9) return std::min(z, low);
    z = low;
  }
  // Not settled: with a cap of optical size, the whole cap disk bounds the bundle's reach (#72).
  // A cap beyond the divergence guard (e.g. a placeholder aperture) would place the start plane
  // where positions and OPL lose their precision; that stays an error (review of #72).
  if (cap < kDivergentRadius) return std::min(z, lowest_on_disk(cap));
  return std::nullopt;
}

}  // namespace rtt::trace::detail
