#include "rtt/compile/layout.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

namespace rtt::compile {
namespace {

using math::Vec3;

/// Tolerance of the check that the section plane is parallel to the local z axis.
constexpr double kParallel = 1e-12;

template <class>
inline constexpr bool kUnhandledAperture = false;

const CompiledSurface& surface_at(const CompiledSystem& system, std::uint32_t surface) {
  if (surface >= system.surfaces().size()) {
    throw std::out_of_range("surface index " + std::to_string(surface) + " out of range (" +
                            std::to_string(system.surfaces().size()) + " surfaces)");
  }
  return system.surfaces()[surface];
}

double sag(const CompiledShape& shape, double x, double y) {
  return std::visit([&](const auto& s) { return s.sag(x, y); }, shape);
}

std::pair<double, double> grad(const CompiledShape& shape, double x, double y) {
  return std::visit([&](const auto& s) { return s.grad(x, y); }, shape);
}

std::optional<double> max_radius(const CompiledShape& shape) {
  return std::visit([](const auto& s) { return s.max_radius(); }, shape);
}

/// A parameter interval [lo, hi] on the section line.
struct Interval {
  double lo = 0.0;
  double hi = 0.0;
};

/// The section line in local coordinates: p(s) = q0 + s t in the local x-y plane, with q0 the
/// point closest to the local z axis (|q0| = |d|) and t the unit direction z_local x n.
struct Line {
  double q0x = 0.0;
  double q0y = 0.0;
  double tx = 0.0;
  double ty = 0.0;
  double d2 = 0.0;  ///< squared distance of the line from the local z axis

  [[nodiscard]] double x(double s) const noexcept { return q0x + s * tx; }
  [[nodiscard]] double y(double s) const noexcept { return q0y + s * ty; }
};

/// Parameters where |p(s)| <= r: s^2 <= r^2 - d^2 (q0 is perpendicular to t). None if the line
/// misses the circle or only touches it.
std::optional<Interval> disc(const Line& line, double r) {
  const double h2 = r * r - line.d2;
  if (!(h2 > 0.0)) return std::nullopt;
  const double h = std::sqrt(h2);
  return Interval{-h, h};
}

/// Intersection of two sets of disjoint intervals, sorted by lo.
std::vector<Interval> intersect(const std::vector<Interval>& a, const std::vector<Interval>& b) {
  std::vector<Interval> out;
  for (const Interval& p : a) {
    for (const Interval& q : b) {
      const Interval i{std::max(p.lo, q.lo), std::min(p.hi, q.hi)};
      if (i.hi > i.lo) out.push_back(i);
    }
  }
  std::sort(out.begin(), out.end(),
            [](const Interval& p, const Interval& q) { return p.lo < q.lo; });
  return out;
}

/// Parameter intervals of the line inside the aperture (local x, y, centred on the axis).
std::vector<Interval> aperture_pieces(const model::Aperture& aperture, const Line& line) {
  return std::visit(
      [&](const auto& a) -> std::vector<Interval> {
        using A = std::decay_t<decltype(a)>;
        if constexpr (std::is_same_v<A, model::CircularAperture>) {
          const std::optional<Interval> outer = disc(line, a.radius);
          if (!outer) return {};
          const std::optional<Interval> inner =
              a.inner_radius > 0.0 ? disc(line, a.inner_radius) : std::nullopt;
          if (!inner) return {*outer};
          return {{outer->lo, inner->lo}, {inner->hi, outer->hi}};
        } else if constexpr (std::is_same_v<A, model::RectangularAperture>) {
          // Slab clipping: |q0 + s t| <= half width in x and in y.
          Interval range{-std::numeric_limits<double>::infinity(),
                         std::numeric_limits<double>::infinity()};
          const auto slab = [&](double q, double t, double half) {
            if (t == 0.0) {
              if (std::abs(q) > half) range = {0.0, 0.0};
              return;
            }
            const double s1 = (-half - q) / t;
            const double s2 = (half - q) / t;
            range.lo = std::max(range.lo, std::min(s1, s2));
            range.hi = std::min(range.hi, std::max(s1, s2));
          };
          slab(line.q0x, line.tx, a.half_width_x);
          slab(line.q0y, line.ty, a.half_width_y);
          if (!(range.hi > range.lo)) return {};
          return {range};
        } else if constexpr (std::is_same_v<A, model::EllipticalAperture>) {
          // (x/a)^2 + (y/b)^2 <= 1 is A s^2 + B s + C <= 0 with A > 0.
          const double ia = 1.0 / (a.semi_axis_x * a.semi_axis_x);
          const double ib = 1.0 / (a.semi_axis_y * a.semi_axis_y);
          const double qa = line.tx * line.tx * ia + line.ty * line.ty * ib;
          const double qb = 2.0 * (line.q0x * line.tx * ia + line.q0y * line.ty * ib);
          const double qc = line.q0x * line.q0x * ia + line.q0y * line.q0y * ib - 1.0;
          const double disc2 = qb * qb - 4.0 * qa * qc;
          if (!(disc2 > 0.0)) return {};
          const double root = std::sqrt(disc2);
          return {{(-qb - root) / (2.0 * qa), (-qb + root) / (2.0 * qa)}};
        } else {
          static_assert(kUnhandledAperture<A>, "aperture type not handled");
        }
      },
      aperture);
}

/// The section line of `plane` in the local coordinates of `s`.
/// @throws std::invalid_argument as surface_profile()
Line section_line(const CompiledSurface& s, const SectionPlane& plane) {
  if (!plane.point.allFinite() || !plane.normal.allFinite() || !(plane.normal.norm() > 0.0)) {
    throw std::invalid_argument("section plane: point and normal must be finite, normal != 0");
  }
  const Vec3 n = s.to_local.apply_vector(plane.normal.normalized());
  if (!(std::abs(n.z()) <= kParallel)) {
    throw std::invalid_argument("section plane is not parallel to the local z axis of surface " +
                                s.id.str() +
                                "; general surface sections come with the 3D meshes (M10)");
  }
  const double norm = std::hypot(n.x(), n.y());
  const double mx = n.x() / norm;
  const double my = n.y() / norm;
  const Vec3 p = s.to_local.apply_point(plane.point);
  const double d = mx * p.x() + my * p.y();
  // t = z x n: rotates (mx, my) by +90 degrees in the local x-y plane.
  return {d * mx, d * my, -my, mx, d * d};
}

/// Parameter intervals of the profile of `s` on `line`: aperture and shape domain.
std::vector<Interval> profile_pieces(const CompiledSurface& s, const Line& line) {
  const std::optional<double> r_max = max_radius(s.shape);
  if (!s.aperture && !r_max) {
    throw std::invalid_argument("surface " + s.id.str() +
                                " is unbounded (no aperture and an unbounded shape): no profile");
  }
  std::vector<Interval> pieces;
  if (s.aperture) {
    pieces = aperture_pieces(*s.aperture, line);
  } else {
    pieces = {{-std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()}};
  }
  if (r_max) {
    const std::optional<Interval> domain = disc(line, *r_max);
    pieces = domain ? intersect(pieces, {*domain}) : std::vector<Interval>{};
  }
  return pieces;
}

/// Point of the profile at parameter s, local coordinates.
Vec3 profile_point(const CompiledSurface& s, const Line& line, double param) {
  const double x = line.x(param);
  const double y = line.y(param);
  return {x, y, sag(s.shape, x, y)};
}

/// Global points of a piece: `samples` points equally spaced in the parameter, both ends exact.
Polyline sample_piece(const CompiledSurface& s,
                      const Line& line,
                      const Interval& piece,
                      std::size_t samples) {
  Polyline points;
  points.reserve(samples);
  const auto last = static_cast<double>(samples - 1);
  for (std::size_t k = 0; k < samples; ++k) {
    const double param = k + 1 == samples
                             ? piece.hi
                             : piece.lo + (piece.hi - piece.lo) * static_cast<double>(k) / last;
    points.push_back(s.to_global.apply_point(profile_point(s, line, param)));
  }
  return points;
}

/// Edge of an outline from rim point `from` to rim point `to` (global; one on surface a, one on
/// the next surface): a step corner if the rims have different distances from the axis of `a`,
/// in a's section coordinates (parameter on a's line, local z of a).
void append_edge(Polyline& outline,
                 const CompiledSurface& a,
                 const Line& line,
                 const Vec3& from,
                 const Vec3& to) {
  const auto coordinates = [&](const Vec3& global) {
    const Vec3 local = a.to_local.apply_point(global);
    const double param = (local.x() - line.q0x) * line.tx + (local.y() - line.q0y) * line.ty;
    return std::pair<double, double>{param, local.z()};
  };
  const auto [s_from, z_from] = coordinates(from);
  const auto [s_to, z_to] = coordinates(to);
  if (std::abs(std::abs(s_from) - std::abs(s_to)) > 1e-12) {
    // Parallel to the axis at the larger rim, then along the section line to the smaller one.
    const bool from_larger = std::abs(s_from) > std::abs(s_to);
    const double s_corner = from_larger ? s_from : s_to;
    const double z_corner = from_larger ? z_to : z_from;
    const Vec3 corner(line.x(s_corner), line.y(s_corner), z_corner);
    outline.push_back(a.to_global.apply_point(corner));
  }
  outline.push_back(to);
}

void check_points(std::span<const double> x,
                  std::span<const double> y,
                  std::size_t out_size,
                  std::size_t per_point) {
  if (x.size() != y.size() || out_size != per_point * x.size()) {
    throw std::invalid_argument("x, y and the output must have matching lengths");
  }
}

}  // namespace

void surface_sag(const CompiledSystem& system,
                 std::uint32_t surface,
                 std::span<const double> x,
                 std::span<const double> y,
                 std::span<double> z) {
  const CompiledSurface& s = surface_at(system, surface);
  check_points(x, y, z.size(), 1);
  for (std::size_t i = 0; i < x.size(); ++i) z[i] = sag(s.shape, x[i], y[i]);
}

void surface_normal(const CompiledSystem& system,
                    std::uint32_t surface,
                    std::span<const double> x,
                    std::span<const double> y,
                    Frame frame,
                    std::span<double> normals) {
  const CompiledSurface& s = surface_at(system, surface);
  check_points(x, y, normals.size(), 3);
  for (std::size_t i = 0; i < x.size(); ++i) {
    const auto [gx, gy] = grad(s.shape, x[i], y[i]);
    Vec3 n = Vec3(-gx, -gy, 1.0).normalized();
    if (frame == Frame::Global) n = s.to_global.apply_vector(n);
    for (std::size_t c = 0; c < 3; ++c) normals[3 * i + c] = n[static_cast<Eigen::Index>(c)];
  }
}

std::vector<Polyline> surface_profile(const CompiledSystem& system,
                                      std::uint32_t surface,
                                      const SectionPlane& plane,
                                      std::size_t samples) {
  const CompiledSurface& s = surface_at(system, surface);
  if (samples < 2) throw std::invalid_argument("a profile needs at least 2 samples");
  const Line line = section_line(s, plane);
  std::vector<Polyline> profile;
  for (const Interval& piece : profile_pieces(s, line)) {
    profile.push_back(sample_piece(s, line, piece, samples));
  }
  return profile;
}

std::vector<Polyline> element_outlines(const CompiledSystem& system,
                                       std::uint32_t element,
                                       const SectionPlane& plane,
                                       std::size_t samples) {
  if (element >= system.elements().size()) {
    throw std::out_of_range("element index " + std::to_string(element) + " out of range (" +
                            std::to_string(system.elements().size()) + " elements)");
  }
  const CompiledElement& e = system.elements()[element];
  // Only bodies whose surfaces bound the segments pairwise (ADR 0017); the others are drawn by
  // their surface profiles.
  if (!e.segmented) return {};
  std::vector<Polyline> outlines;
  for (std::uint32_t j = 0; j + 1 < e.surface_count; ++j) {
    const std::uint32_t ia = e.first_surface + j;
    const CompiledSurface& a = system.surfaces()[ia];
    const std::vector<Polyline> front = surface_profile(system, ia, plane, samples);
    const std::vector<Polyline> back = surface_profile(system, ia + 1, plane, samples);
    if (front.size() != back.size()) {
      throw std::invalid_argument(
          "element " + e.name + ", segment " + std::to_string(j) + ": the section plane cuts " +
          std::to_string(front.size()) + " piece(s) of " + a.id.str() + " but " +
          std::to_string(back.size()) + " of " + system.surfaces()[ia + 1].id.str());
    }
    const Line line = section_line(a, plane);
    for (std::size_t k = 0; k < front.size(); ++k) {
      Polyline outline(front[k]);
      append_edge(outline, a, line, front[k].back(), back[k].back());
      outline.insert(outline.end(), back[k].rbegin() + 1, back[k].rend());
      append_edge(outline, a, line, back[k].front(), front[k].front());
      outlines.push_back(std::move(outline));
    }
  }
  return outlines;
}

}  // namespace rtt::compile
