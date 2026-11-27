#pragma once

/// @file aabb.hpp
/// Axis-aligned bounding boxes for the BVH of the non-sequential scene (ADR 0057, point 3).
/// Coordinates in mm, global unless stated otherwise.

#include <array>
#include <limits>

#include "rtt/math/types.hpp"

namespace rtt::nonseq {

/// Axis-aligned box [lo, hi] in mm. The default box is empty (lo = +inf, hi = -inf), so that
/// growing it by a point gives that point.
struct Aabb {
  math::Vec3 lo = math::Vec3::Constant(std::numeric_limits<double>::infinity());
  math::Vec3 hi = math::Vec3::Constant(-std::numeric_limits<double>::infinity());

  /// True if lo <= hi on every axis (a single point is not empty).
  [[nodiscard]] bool empty() const noexcept {
    return !(lo.x() <= hi.x() && lo.y() <= hi.y() && lo.z() <= hi.z());
  }
  /// True if the box is not empty and every coordinate is finite.
  [[nodiscard]] bool finite() const noexcept {
    return !empty() && lo.allFinite() && hi.allFinite();
  }
  /// Grows the box to contain `p`.
  void grow(const math::Vec3& p) noexcept {
    lo = lo.cwiseMin(p);
    hi = hi.cwiseMax(p);
  }
  /// Grows the box to contain `b`.
  void grow(const Aabb& b) noexcept {
    lo = lo.cwiseMin(b.lo);
    hi = hi.cwiseMax(b.hi);
  }
  /// Centre (lo + hi) / 2, mm.
  [[nodiscard]] math::Vec3 centroid() const noexcept { return 0.5 * (lo + hi); }
  /// Surface area 2 (dx dy + dy dz + dz dx) in mm^2; 0 for an empty box.
  [[nodiscard]] double surface_area() const noexcept {
    if (empty()) return 0.0;
    const math::Vec3 d = hi - lo;
    return 2.0 * (d.x() * d.y() + d.y() * d.z() + d.z() * d.x());
  }
};

/// gamma_n = n eps_m / (1 - n eps_m) with eps_m = epsilon / 2 the unit roundoff (PBRT 4ed, §6.8.1,
/// definition after Eq. (6.23); docs/quellen.md): bound of the relative error of n rounded
/// floating-point operations.
[[nodiscard]] constexpr double gamma_bound(int n) noexcept {
  constexpr double eps_m = std::numeric_limits<double>::epsilon() * 0.5;
  return (n * eps_m) / (1.0 - n * eps_m);
}

/// Precomputed values of a ray for the slab test (PBRT 4ed, §7.3.5): origin, reciprocal
/// direction (infinite for a zero component) and whether each component is negative.
struct SlabRay {
  math::Vec3 origin = math::Vec3::Zero();
  math::Vec3 inv_dir = math::Vec3::Ones();
  std::array<bool, 3> dir_is_neg{};

  /// Prepares the slab test for the ray origin + t direction (global, mm).
  [[nodiscard]] static SlabRay from(const math::Vec3& origin, const math::Vec3& direction) noexcept;
};

/// Slab test of a ray against a box (PBRT 4ed, §6.1.2 and §6.8.2): true if the ray meets the box
/// for some t in [0, t_max) or starts inside it; `t_near` receives the entry parameter (0 if the
/// origin is inside). The far parameter of every slab is enlarged by 2 gamma_3 (§6.8.2), so that
/// rounding never loses a box the ray really meets. NaN from 0 * inf (origin on a slab plane of
/// an axis-parallel ray) fails every comparison and leaves the interval unchanged on every axis.
[[nodiscard]] bool slab_test(const Aabb& box,
                             const SlabRay& ray,
                             double t_max,
                             double& t_near) noexcept;

}  // namespace rtt::nonseq
