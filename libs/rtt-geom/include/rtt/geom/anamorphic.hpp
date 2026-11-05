#pragma once

/// @file anamorphic.hpp
/// Anamorphic shapes (ADR 0042, point 1): biconic, toroid and cylinder in local coordinates
/// (lengths in mm). Their sections x-z and y-z have different curvatures; the orientation of the
/// sections is that of the local x and y axes, set by the rotation of the surface about z.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <utility>
#include <vector>

#include "rtt/geom/asphere.hpp"
#include "rtt/geom/shape.hpp"

namespace rtt::geom {

namespace detail {

/// +-infinity with the sign of v, 0 for v = 0: a slope component at a vertical rim.
template <rtt::math::Real T>
T vertical_slope(T v) noexcept {
  const T inf = std::numeric_limits<T>::infinity();
  return v == T(0) ? T(0) : (v > T(0) ? inf : -inf);
}

}  // namespace detail

/// Biconic (Ansys Zemax OpticStudio User Guide, Release 2025 R1, "Biconic"; equation without
/// number):
///   z = (c_x x^2 + c_y y^2) / (1 + S),  S = sqrt(1 - a x^2 - b y^2),
///   a = (1 + k_x) c_x^2,  b = (1 + k_y) c_y^2,
/// with c_x = 1/R_x, c_y = 1/R_y in 1/mm and the conic constants k_x, k_y. Sign conventions of
/// c and k in each section as for Conic. Defined where 1 - a x^2 - b y^2 >= 0.
template <rtt::math::Real T>
class Biconic final {
 public:
  /// @param curvature_x c_x in 1/mm, section x-z
  /// @param curvature_y c_y in 1/mm, section y-z
  /// @param conic_x     k_x, dimensionless
  /// @param conic_y     k_y, dimensionless
  Biconic(T curvature_x, T curvature_y, T conic_x, T conic_y)
      : cx_(curvature_x), cy_(curvature_y), kx_(conic_x), ky_(conic_y) {}

  [[nodiscard]] T curvature_x() const noexcept { return cx_; }
  [[nodiscard]] T curvature_y() const noexcept { return cy_; }
  [[nodiscard]] T conic_x() const noexcept { return kx_; }
  [[nodiscard]] T conic_y() const noexcept { return ky_; }

  /// Sag z in mm at (x, y) in mm; NaN outside the domain.
  [[nodiscard]] T sag(T x, T y) const noexcept {
    using std::sqrt;
    const T s2 = radicand(x, y);
    if (s2 < T(0)) return std::numeric_limits<T>::quiet_NaN();
    return (cx_ * (x * x) + cy_ * (y * y)) / (T(1) + sqrt(s2));
  }

  /// Slopes (dz/dx, dz/dy), derived from the sag (ADR 0042, point 1):
  ///   dz/dx = 2 c_x x / (1 + S) + a x N / (S (1 + S)^2),  N = c_x x^2 + c_y y^2,
  /// and dz/dy likewise with c_y and b. For equal sections this is c x / S, the conic's.
  /// NaN outside the domain, +-infinity on its rim (S = 0).
  [[nodiscard]] std::pair<T, T> grad(T x, T y) const noexcept {
    using std::sqrt;
    const T s2 = radicand(x, y);
    if (s2 < T(0)) {
      const T nan = std::numeric_limits<T>::quiet_NaN();
      return {nan, nan};
    }
    const T s = sqrt(s2);
    const T n = cx_ * (x * x) + cy_ * (y * y);
    if (s == T(0)) {
      return {detail::vertical_slope(a() * x * n), detail::vertical_slope(b() * y * n)};
    }
    const T p = T(1) + s;
    const T q = s * p * p;
    return {T(2) * cx_ * x / p + a() * x * n / q, T(2) * cy_ * y / p + b() * y * n / q};
  }

  /// Start value of the intersection: (c, k) of the section with the larger (1 + k) c^2, whose
  /// disk of definition is the largest disk inside the ellipse; y for equal values.
  [[nodiscard]] std::pair<T, T> base_conic() const {
    return a() > b() ? std::pair{cx_, kx_} : std::pair{cy_, ky_};
  }

  /// Outer extent of the domain: the major semi-axis of the ellipse a x^2 + b y^2 <= 1, as
  /// 1 / (|c| sqrt(1 + k)) of the section with the smaller (1 + k) c^2 (the formula of Conic, so
  /// equal sections give the conic's value bit for bit); unbounded unless a > 0 and b > 0.
  [[nodiscard]] std::optional<T> max_radius() const {
    using std::abs;
    using std::sqrt;
    if (!(a() > T(0) && b() > T(0))) return std::nullopt;
    return a() < b() ? T(1) / (abs(cx_) * sqrt(T(1) + kx_)) : T(1) / (abs(cy_) * sqrt(T(1) + ky_));
  }

  /// True where the sag is defined: 1 - a x^2 - b y^2 >= 0.
  [[nodiscard]] bool in_domain(T x, T y) const noexcept { return radicand(x, y) >= T(0); }

  /// Vertex at the origin, Hessian (c_x, 0, c_y), order 4 (ADR 0042, point 1). Rotational to
  /// fourth order exactly when both sections have the same c and k; then c is assigned to both
  /// diagonal elements once and the cap is (k, A4 = 0), the conic's.
  [[nodiscard]] VertexExpansion<T> vertex_expansion() const noexcept {
    VertexExpansion<T> e;
    e.order = 4;
    if (cx_ == cy_ && kx_ == ky_) {
      e.hessian = {cy_, T(0), cy_};
      e.rotational4 = true;
      e.conic = ky_;
    } else {
      e.hessian = {cx_, T(0), cy_};
    }
    return e;
  }

 private:
  [[nodiscard]] T a() const noexcept { return (T(1) + kx_) * cx_ * cx_; }
  [[nodiscard]] T b() const noexcept { return (T(1) + ky_) * cy_ * cy_; }
  [[nodiscard]] T radicand(T x, T y) const noexcept { return T(1) - a() * (x * x) - b() * (y * y); }

  T cx_;
  T cy_;
  T kx_;
  T ky_;
};

/// Cylinder (ADR 0042, point 1): the y-z profile of the toroid, flat in x,
///   z = c y^2 / (1 + sqrt(1 - (1 + k) c^2 y^2)) + A4 y^4 + A6 y^6 + ...,
/// the even asphere of the same parameters evaluated at (0, y). Defined on the strip where the
/// profile is.
template <rtt::math::Real T>
class Cylinder final {
 public:
  /// @param curvature      c = 1/R of the profile in 1/mm
  /// @param conic_constant k of the profile, dimensionless
  /// @param coefficients   A4, A6, ... (coefficients[0] = A4), A_n in mm^(1 - n)
  Cylinder(T curvature, T conic_constant, std::vector<T> coefficients)
      : profile_(curvature, conic_constant, std::move(coefficients)) {}

  /// The profile as an even asphere; its sag at (0, y) is z_y(y).
  [[nodiscard]] const EvenAsphere<T>& profile() const noexcept { return profile_; }

  /// Sag z = z_y(y) in mm; NaN outside the strip of the profile.
  [[nodiscard]] T sag(T /*x*/, T y) const noexcept { return profile_.sag(T(0), y); }

  /// Slopes (0, z_y'(y)).
  [[nodiscard]] std::pair<T, T> grad(T /*x*/, T y) const noexcept {
    return {T(0), profile_.grad(T(0), y).second};
  }

  /// (c, k) of the profile.
  [[nodiscard]] std::pair<T, T> base_conic() const { return profile_.base_conic(); }

  /// Unbounded: the strip extends in x without limit; the aperture bounds the surface.
  [[nodiscard]] std::optional<T> max_radius() const { return std::nullopt; }

  /// True where the sag is defined.
  [[nodiscard]] bool in_domain(T x, T y) const noexcept {
    using std::isnan;
    return !isnan(sag(x, y));
  }

  /// Vertex at the origin, Hessian (0, 0, c), order 4, never rotational (ADR 0042, point 1).
  [[nodiscard]] VertexExpansion<T> vertex_expansion() const noexcept {
    VertexExpansion<T> e;
    e.hessian = {T(0), T(0), profile_.base_conic().first};
    e.order = 4;
    return e;
  }

 private:
  EvenAsphere<T> profile_;
};

/// Toroid (Ansys Zemax OpticStudio User Guide, Release 2025 R1, "Toroidal"): the profile
/// z_y(y) of Cylinder rotated about the axis parallel to y through (x, z) = (0, R), R the
/// signed radius of rotation, c_R = 1/R. The guide gives the surface only in words; ADR 0042,
/// point 1 derives
///   z = z_y + c' x^2 / (1 + sqrt(1 - c'^2 x^2)),  c' = c_R / (1 - c_R z_y(y)),
/// from x^2 + (z - R)^2 = (R - z_y)^2 on the branch through the vertex. c_R = 0 is the
/// cylinder. Defined where the profile is, where the profile has not crossed the axis
/// (1 - c_R z_y > 0) and where |x| <= |R - z_y|.
template <rtt::math::Real T>
class Toroid final {
 public:
  /// @param curvature             c = 1/R_y of the profile in 1/mm
  /// @param conic_constant        k of the profile, dimensionless
  /// @param coefficients          A4, A6, ... of the profile, A_n in mm^(1 - n)
  /// @param curvature_of_rotation c_R = 1/R in 1/mm, R the radius of rotation (sign as for c:
  ///                              c_R > 0 puts the axis on the +z side)
  Toroid(T curvature, T conic_constant, std::vector<T> coefficients, T curvature_of_rotation)
      : profile_(curvature, conic_constant, std::move(coefficients)),
        cr_(curvature_of_rotation),
        max_radius_(outer_extent()) {}

  /// The profile as an even asphere; its sag at (0, y) is z_y(y).
  [[nodiscard]] const EvenAsphere<T>& profile() const noexcept { return profile_; }

  /// c_R = 1/R in 1/mm.
  [[nodiscard]] T curvature_of_rotation() const noexcept { return cr_; }

  /// Sag z in mm at (x, y) in mm; NaN outside the domain.
  [[nodiscard]] T sag(T x, T y) const noexcept {
    using std::sqrt;
    const T zy = profile_.sag(T(0), y);
    const T g = T(1) - cr_ * zy;
    if (!(g > T(0))) return std::numeric_limits<T>::quiet_NaN();
    const T cp = cr_ / g;
    const T s2 = T(1) - cp * cp * (x * x);
    if (s2 < T(0)) return std::numeric_limits<T>::quiet_NaN();
    return zy + cp * (x * x) / (T(1) + sqrt(s2));
  }

  /// Slopes (ADR 0042, point 1): dz/dx = c' x / s, dz/dy = z_y'(y) / s, s = sqrt(1 - c'^2 x^2).
  /// NaN outside the domain, +-infinity where the surface is vertical (s = 0).
  [[nodiscard]] std::pair<T, T> grad(T x, T y) const noexcept {
    using std::sqrt;
    const T nan = std::numeric_limits<T>::quiet_NaN();
    const T zy = profile_.sag(T(0), y);
    const T g = T(1) - cr_ * zy;
    if (!(g > T(0))) return {nan, nan};
    const T cp = cr_ / g;
    const T s2 = T(1) - cp * cp * (x * x);
    if (s2 < T(0)) return {nan, nan};
    const T s = sqrt(s2);
    const T zyp = profile_.grad(T(0), y).second;
    if (s == T(0)) {
      return {detail::vertical_slope(cp * x), detail::vertical_slope(zyp)};
    }
    return {cp * x / s, zyp / s};
  }

  /// Start value of the intersection: the section with the larger (1 + k) c^2, the circle of
  /// rotation (c_R, 0) or the profile (c, k); the profile for equal values.
  [[nodiscard]] std::pair<T, T> base_conic() const {
    const auto [c, k] = profile_.base_conic();
    return cr_ * cr_ > (T(1) + k) * c * c ? std::pair{cr_, T(0)} : std::pair{c, k};
  }

  /// Outer extent of the domain (ADR 0042, point 1): sup sqrt(y^2 + d(y)^2) over the profile up
  /// to the first point where it meets the axis of rotation, d = (1 - c_R z_y) / |c_R|;
  /// unbounded for c_R = 0 or an unbounded profile that never meets the axis. Computed once.
  [[nodiscard]] std::optional<T> max_radius() const { return max_radius_; }

  /// True where the sag is defined.
  [[nodiscard]] bool in_domain(T x, T y) const noexcept {
    using std::isnan;
    return !isnan(sag(x, y));
  }

  /// Vertex at the origin, Hessian (c_R, 0, c), order 4 (ADR 0042, point 1). Rotational only for
  /// the sphere: c_R = c, k = 0 and no non-zero coefficient; then the cap is (0, 0).
  [[nodiscard]] VertexExpansion<T> vertex_expansion() const noexcept {
    const auto [c, k] = profile_.base_conic();
    VertexExpansion<T> e;
    e.order = 4;
    const auto& a = profile_.coefficients();
    if (cr_ == c && k == T(0) && std::all_of(a.begin(), a.end(), [](T v) { return v == T(0); })) {
      e.hessian = {c, T(0), c};
      e.rotational4 = true;
    } else {
      e.hessian = {cr_, T(0), c};
    }
    return e;
  }

 private:
  /// Samples of the profile for the outer extent (ADR 0042, point 1).
  static constexpr int kSamples = 4096;
  /// Relative margin by which the sampled supremum is enlarged.
  static constexpr double kExtentMargin = 1e-6;

  /// 1 - c_R z_y(y): positive while the profile has not reached the axis; NaN outside the
  /// profile.
  [[nodiscard]] T gap(T y) const noexcept { return T(1) - cr_ * profile_.sag(T(0), y); }

  /// Largest y in [lo, hi] with gap > 0, given gap(lo) > 0 and not gap(hi) > 0 (bisection).
  [[nodiscard]] T last_inside(T lo, T hi) const noexcept {
    for (int i = 0; i < 200 && hi - lo > std::numeric_limits<T>::epsilon() * hi; ++i) {
      const T mid = lo + (hi - lo) / T(2);
      if (gap(mid) > T(0)) {
        lo = mid;
      } else {
        hi = mid;
      }
    }
    return lo;
  }

  /// First y >= 0 where the profile ends or meets the axis, searched on [0, limit] with
  /// kSamples steps and refined by bisection; limit if there is none.
  [[nodiscard]] T end_of_profile(T limit) const noexcept {
    T previous = T(0);
    for (int i = 1; i <= kSamples; ++i) {
      const T y = limit * static_cast<T>(i) / static_cast<T>(kSamples);
      if (!(gap(y) > T(0))) return last_inside(previous, y);
      previous = y;
    }
    return limit;
  }

  [[nodiscard]] std::optional<T> outer_extent() const {
    using std::abs;
    using std::sqrt;
    if (cr_ == T(0)) return std::nullopt;
    T y_end = T(0);
    if (const std::optional<T> y_lim = profile_.max_radius()) {
      y_end = end_of_profile(*y_lim);
    } else {
      // Unbounded profile: double y until the profile meets the axis; never meeting it leaves
      // the toroid unbounded.
      T y = T(1) / abs(cr_);
      while (gap(y) > T(0)) {
        y *= T(2);
        if (!std::isfinite(y)) return std::nullopt;
      }
      y_end = end_of_profile(y);
    }
    T sup = T(0);
    for (int i = 0; i <= kSamples; ++i) {
      const T y = y_end * static_cast<T>(i) / static_cast<T>(kSamples);
      const T g = gap(y);
      const T d = g > T(0) ? g / abs(cr_) : T(0);
      sup = std::max(sup, y * y + d * d);
    }
    return sqrt(sup) * (T(1) + T(kExtentMargin));
  }

  EvenAsphere<T> profile_;
  T cr_;
  std::optional<T> max_radius_;
};

}  // namespace rtt::geom
