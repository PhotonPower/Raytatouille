#include "rtt/nonseq/scene.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numbers>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include "rtt/diagnostics/codes.hpp"
#include "rtt/geom/intersect.hpp"
#include "rtt/math/isometry.hpp"
#include "rtt/model/element.hpp"
#include "rtt/model/surface.hpp"

namespace rtt::nonseq {

namespace {

using compile::CompiledSurface;
using math::Vec3;
using model::ElementKind;

constexpr double kInf = std::numeric_limits<double>::infinity();
/// Tolerance of the coaxial and equal-rim checks of a mantle (ADR 0057, point 2), relative.
constexpr double kRimRelTol = 1e-12;
/// Points per axis of the bounding grid and points on a rim (ADR 0057, point 3).
constexpr int kGridPoints = 64;
constexpr int kRimPoints = 256;
/// Relative enlargement of every box (ADR 0057, point 3).
constexpr double kBoxRelMargin = 1e-9;
/// Crossings of one face tried per query: a conic has at most two; the bound only stops a
/// Newton surface that keeps finding further roots.
constexpr int kMaxCrossings = 8;
/// Rim points of a circle or ellipse that stand for its corners in the prism checks.
constexpr int kCornerPoints = 64;
/// Samples of F = z - sag along the ray when the next crossing of a Newton surface is
/// bracketed (ADR 0057, Nachtrag PR 2, i), and the bound of the bisection steps (2^-200 of
/// the interval is far below the rounding limit).
constexpr int kBracketSamples = 64;
constexpr int kMaxBisections = 200;

double sag(const CompiledSurface& s, double x, double y) noexcept {
  return std::visit([x, y](const auto& shape) { return shape.sag(x, y); }, s.shape);
}

std::pair<double, double> grad(const CompiledSurface& s, double x, double y) noexcept {
  return std::visit([x, y](const auto& shape) { return shape.grad(x, y); }, s.shape);
}

/// Shapes whose sag depends only on r (ADR 0041, ADR 0043): the rim height of a circle is then
/// constant (ADR 0057, point 2, step at unequal radii).
bool rotational(const CompiledSurface& s) noexcept {
  return std::holds_alternative<geom::Plane<double>>(s.shape) ||
         std::holds_alternative<geom::Conic<double>>(s.shape) ||
         std::holds_alternative<geom::EvenAsphere<double>>(s.shape) ||
         std::holds_alternative<geom::OddAsphere<double>>(s.shape) ||
         std::holds_alternative<geom::Axicon<double>>(s.shape) ||
         std::holds_alternative<geom::QCon<double>>(s.shape) ||
         std::holds_alternative<geom::QBfs<double>>(s.shape);
}

/// Shapes with an analytic intersection (rtt-geom: plane, conic, axicon): with t_min they find
/// every crossing; every other shape intersects by Newton from its base conic.
bool analytic(const CompiledSurface& s) noexcept {
  return std::holds_alternative<geom::Plane<double>>(s.shape) ||
         std::holds_alternative<geom::Conic<double>>(s.shape) ||
         std::holds_alternative<geom::Axicon<double>>(s.shape);
}

/// JSON pointer of the element of a surface: the surface pointer without "/surfaces/<k>".
std::string element_location(const CompiledSurface& s) {
  const std::size_t pos = s.location.rfind("/surfaces/");
  return pos == std::string::npos ? s.location : s.location.substr(0, pos);
}

Vec3 axis_z(const CompiledSurface& s) {
  return s.to_global.rotation().col(2);
}
Vec3 axis_x(const CompiledSurface& s) {
  return s.to_global.rotation().col(0);
}
Vec3 vertex(const CompiledSurface& s) {
  return s.to_global.translation();
}

double sign_of(double v) noexcept {
  return v > 0.0 ? 1.0 : (v < 0.0 ? -1.0 : 0.0);
}

/// Closed polygon of a rim in local x, y, counter-clockwise (positive signed area).
std::vector<std::array<double, 2>> ccw(std::vector<std::array<double, 2>> v) {
  double area2 = 0.0;
  for (std::size_t i = 0; i < v.size(); ++i) {
    const auto& p = v[i];
    const auto& q = v[(i + 1) % v.size()];
    area2 += p[0] * q[1] - q[0] * p[1];
  }
  if (area2 < 0.0) std::reverse(v.begin(), v.end());
  return v;
}

/// Polygon of an aperture whose rim is straight (rectangle and polygons); empty otherwise.
std::vector<std::array<double, 2>> rim_polygon(const model::Aperture& a) {
  if (const auto* r = std::get_if<model::RectangularAperture>(&a)) {
    const double x = r->half_width_x.value;
    const double y = r->half_width_y.value;
    return {{{x, -y}}, {{x, y}}, {{-x, y}}, {{-x, -y}}};
  }
  return ccw(model::aperture_vertices(a));
}

/// Points on the outer rim of an aperture (local x, y): `n` points around a circle or ellipse;
/// the vertices of a straight rim plus points spread over its edges by length.
std::vector<std::array<double, 2>> rim_points(const model::Aperture& a, int n) {
  std::vector<std::array<double, 2>> pts;
  const auto ellipse = [&](double ax, double by) {
    for (int k = 0; k < n; ++k) {
      const double phi = 2.0 * std::numbers::pi * k / n;
      pts.push_back({ax * std::cos(phi), by * std::sin(phi)});
    }
  };
  if (const auto* c = std::get_if<model::CircularAperture>(&a)) {
    ellipse(c->radius.value, c->radius.value);
  } else if (const auto* e = std::get_if<model::EllipticalAperture>(&a)) {
    ellipse(e->semi_axis_x.value, e->semi_axis_y.value);
  } else {
    const auto poly = rim_polygon(a);
    double perimeter = 0.0;
    for (std::size_t i = 0; i < poly.size(); ++i) {
      const auto& p = poly[i];
      const auto& q = poly[(i + 1) % poly.size()];
      perimeter += std::hypot(q[0] - p[0], q[1] - p[1]);
    }
    for (std::size_t i = 0; i < poly.size(); ++i) {
      const auto& p = poly[i];
      const auto& q = poly[(i + 1) % poly.size()];
      const double len = std::hypot(q[0] - p[0], q[1] - p[1]);
      const int m = std::max(1, static_cast<int>(std::ceil(n * len / perimeter)));
      for (int k = 0; k < m; ++k) {
        const double u = static_cast<double>(k) / m;
        pts.push_back({p[0] + u * (q[0] - p[0]), p[1] + u * (q[1] - p[1])});
      }
    }
  }
  return pts;
}

/// Inner radius of a circular aperture, 0 for none or another kind.
double inner_radius(const model::Aperture& a) {
  const auto* c = std::get_if<model::CircularAperture>(&a);
  return c != nullptr ? c->inner_radius.value : 0.0;
}

/// Corners of an aperture for the prism checks (ADR 0057, point 2): the vertices of a straight
/// rim, kCornerPoints rim points of a circle or ellipse.
std::vector<std::array<double, 2>> corners(const model::Aperture& a) {
  if (std::holds_alternative<model::CircularAperture>(a) ||
      std::holds_alternative<model::EllipticalAperture>(a)) {
    return rim_points(a, kCornerPoints);
  }
  return rim_polygon(a);
}

/// The aperture scaled by `s` about the local origin, without inner radius (screen of a stop).
model::Aperture scaled(const model::Aperture& a, double s) {
  return std::visit(
      [s](const auto& ap) -> model::Aperture {
        using A = std::remove_cvref_t<decltype(ap)>;
        A out = ap;
        if constexpr (std::is_same_v<A, model::CircularAperture>) {
          out.radius.value *= s;
          out.inner_radius.value = 0.0;
        } else if constexpr (std::is_same_v<A, model::RectangularAperture>) {
          out.half_width_x.value *= s;
          out.half_width_y.value *= s;
        } else if constexpr (std::is_same_v<A, model::EllipticalAperture>) {
          out.semi_axis_x.value *= s;
          out.semi_axis_y.value *= s;
        } else if constexpr (std::is_same_v<A, model::RegularPolygonAperture>) {
          out.circumradius.value *= s;
        } else {
          for (auto& v : out.vertices) {
            v[0] *= s;
            v[1] *= s;
          }
        }
        return out;
      },
      a);
}

/// Box in local x, y of an aperture: [x0, x1] x [y0, y1].
std::array<double, 4> aperture_xy(const model::Aperture& a) {
  if (const auto* c = std::get_if<model::CircularAperture>(&a)) {
    const double r = c->radius.value;
    return {-r, r, -r, r};
  }
  if (const auto* e = std::get_if<model::EllipticalAperture>(&a)) {
    return {-e->semi_axis_x.value, e->semi_axis_x.value, -e->semi_axis_y.value,
            e->semi_axis_y.value};
  }
  std::array<double, 4> b{kInf, -kInf, kInf, -kInf};
  for (const auto& v : rim_polygon(a)) {
    b[0] = std::min(b[0], v[0]);
    b[1] = std::max(b[1], v[0]);
    b[2] = std::min(b[2], v[1]);
    b[3] = std::max(b[3], v[1]);
  }
  return b;
}

/// Adds the relative margin of ADR 0057, point 3.
Aabb with_margin(Aabb b) {
  if (!b.finite()) return b;
  const double scale = 1.0 + std::max(b.lo.cwiseAbs().maxCoeff(), b.hi.cwiseAbs().maxCoeff());
  const Vec3 m = Vec3::Constant(kBoxRelMargin * scale);
  b.lo -= m;
  b.hi += m;
  return b;
}

/// Global box of a local box: the box of its eight transformed corners.
Aabb to_global_box(const Aabb& local, const math::Isometry3& to_global) {
  if (!local.finite()) return local;
  Aabb g;
  for (int k = 0; k < 8; ++k) {
    const Vec3 c((k & 1) != 0 ? local.hi.x() : local.lo.x(),
                 (k & 2) != 0 ? local.hi.y() : local.lo.y(),
                 (k & 4) != 0 ? local.hi.z() : local.lo.z());
    g.grow(to_global.apply_point(c));
  }
  return with_margin(g);
}

/// Local box of `surface` over its aperture (ADR 0057, point 3; ADR 0057, Nachtrag PR 2, d): a
/// grid of 64 x 64 points over the aperture box (points inside the aperture) and 256 points on
/// each rim; the sag range grows by h max|grad sag| with h the larger grid step. Points outside
/// the domain of the shape are skipped. Not finite without aperture or with an infinite slope.
Aabb local_surface_box(const CompiledSurface& surface) {
  Aabb box;
  if (!surface.aperture) return box;
  const model::Aperture& ap = *surface.aperture;
  const auto xy = aperture_xy(ap);
  const double hx = (xy[1] - xy[0]) / (kGridPoints - 1);
  const double hy = (xy[3] - xy[2]) / (kGridPoints - 1);
  double z_lo = kInf;
  double z_hi = -kInf;
  double g_max = 0.0;
  const auto take = [&](double x, double y) {
    const double z = sag(surface, x, y);
    if (!std::isfinite(z)) return;
    const auto [gx, gy] = grad(surface, x, y);
    z_lo = std::min(z_lo, z);
    z_hi = std::max(z_hi, z);
    g_max = std::max(g_max, std::hypot(gx, gy));
  };
  trace::SurfaceHit probe;
  for (int i = 0; i < kGridPoints; ++i) {
    for (int k = 0; k < kGridPoints; ++k) {
      probe.point = Vec3(xy[0] + i * hx, xy[2] + k * hy, 0.0);
      if (trace::inside_aperture(surface, probe)) take(probe.point.x(), probe.point.y());
    }
  }
  for (const auto& p : rim_points(ap, kRimPoints)) take(p[0], p[1]);
  if (const double ri = inner_radius(ap); ri > 0.0) {
    for (int k = 0; k < kRimPoints; ++k) {
      const double phi = 2.0 * std::numbers::pi * k / kRimPoints;
      take(ri * std::cos(phi), ri * std::sin(phi));
    }
  }
  if (!(z_lo <= z_hi)) return box;
  const double grow = std::max(hx, hy) * g_max;
  box.lo = Vec3(xy[0], xy[2], z_lo - grow);
  box.hi = Vec3(xy[1], xy[3], z_hi + grow);
  return box;
}

/// z in the frame of `frame` of the point of surface `other` above the frame point (x, y).
double height_in_frame(const CompiledSurface& frame,
                       const CompiledSurface& other,
                       double x,
                       double y) {
  const Vec3 q = other.to_local.apply_point(frame.to_global.apply_point(Vec3(x, y, 0.0)));
  const Vec3 on = other.to_global.apply_point(Vec3(q.x(), q.y(), sag(other, q.x(), q.y())));
  return frame.to_local.apply_point(on).z();
}

}  // namespace

SceneError::SceneError(std::vector<model::Diagnostic> diagnostics)
    : std::invalid_argument([&] {
        std::string m = "nonseq: the scene cannot be built:";
        for (const auto& d : diagnostics)
          m += "\n  [" + d.code + "] " + d.location + ": " + d.message;
        return m;
      }()),
      diagnostics_(std::move(diagnostics)) {}

/// Builds the faces of a Scene (ADR 0057, point 2).
class SceneBuilder {
 public:
  SceneBuilder(const compile::CompiledSystem& system, const SceneOptions& options)
      : system_(system), options_(options) {}

  Scene run() {
    check_options();
    scene_.environment_ = system_.environment_medium();
    for (std::uint32_t e = 0; e < system_.elements().size(); ++e) {
      const compile::CompiledElement& el = system_.elements()[e];
      if (!el.active || el.surface_count == 0) continue;
      switch (el.kind) {
        case ElementKind::Lens:
          add_segmented(e);
          break;
        case ElementKind::Plate:
          if (el.segmented || el.surface_count == 2) {
            add_segmented(e);
          } else {
            add_prism(e);
          }
          break;
        case ElementKind::Mirror:
        case ElementKind::ThinElement:
        case ElementKind::Stop:
        case ElementKind::Detector:
          add_thin(e);
          break;
      }
    }
    if (!errors_.empty()) throw SceneError(std::move(errors_));
    finish();
    return std::move(scene_);
  }

 private:
  void check_options() const {
    if (!(std::isfinite(options_.screen_scale) && options_.screen_scale > 1.0)) {
      throw std::invalid_argument("nonseq: screen_scale must be finite and > 1");
    }
    const auto& els = system_.elements();
    for (const std::uint32_t e : options_.absorbing_rims) {
      if (e >= els.size() ||
          (els[e].kind != ElementKind::Lens && els[e].kind != ElementKind::Plate)) {
        throw std::invalid_argument("nonseq: absorbing_rims names no lens or plate");
      }
    }
    for (const std::uint32_t e : options_.transparent_detectors) {
      if (e >= els.size() || els[e].kind != ElementKind::Detector) {
        throw std::invalid_argument("nonseq: transparent_detectors names no detector");
      }
    }
  }

  void report(diagnostics::DiagnosticCode code, std::string location, std::string message) {
    errors_.push_back(
        {code.severity(), std::move(location), std::move(message), std::string(code.str())});
  }

  const CompiledSurface& sys(std::uint32_t s) const { return system_.surfaces()[s]; }

  /// Copies a system surface into the scene and returns its scene index.
  std::uint32_t copy_surface(std::uint32_t s) {
    scene_.surfaces_.push_back(sys(s));
    return static_cast<std::uint32_t>(scene_.surfaces_.size() - 1);
  }

  /// Adds a face with its box in the frame of its scene surface (margin included).
  void add_face(const Face& f, const Aabb& local_box) {
    scene_.faces_.push_back(f);
    scene_.local_bounds_.push_back(local_box);
    scene_.bounds_.push_back(to_global_box(local_box, scene_.surfaces_[f.surface].to_global));
  }

  /// Local box of a surface face with margin.
  static Aabb surface_box(const CompiledSurface& s) { return with_margin(local_surface_box(s)); }

  /// Lens, segmented Plate, Plate of one material with two surfaces (segment rule, ADR 0017).
  void add_segmented(std::uint32_t e) {
    const compile::CompiledElement& el = system_.elements()[e];
    const std::uint32_t f = el.first_surface;
    const std::uint32_t n = el.surface_count;
    // Side of segment j on surface j (low) and on surface j + 1 (up), ADR 0057, point 2.
    std::vector<double> low(n - 1, 0.0);
    std::vector<double> up(n - 1, 0.0);
    bool ok = true;
    for (std::uint32_t j = 0; j + 1 < n; ++j) {
      const Vec3 d = vertex(sys(f + j + 1)) - vertex(sys(f + j));
      low[j] = sign_of(axis_z(sys(f + j)).dot(d));
      up[j] = sign_of(-axis_z(sys(f + j + 1)).dot(d));
      if (low[j] == 0.0 || up[j] == 0.0) {
        report("nonseq.segment_thickness", element_location(sys(f)),
               "segment " + std::to_string(j) + " of element '" + el.name +
                   "' has centre thickness 0 along the surface axes");
        ok = false;
      }
    }
    for (std::uint32_t i = 1; ok && i + 1 < n; ++i) {
      if (up[i - 1] == low[i]) {
        report("nonseq.segment_thickness", sys(f + i).location,
               "the segments on both sides of surface '" + sys(f + i).id.str() +
                   "' lie on the same side: they overlap");
        ok = false;
      }
    }
    std::vector<std::uint32_t> scene_index(n, 0);
    for (std::uint32_t i = 0; i < n; ++i) {
      if (!sys(f + i).aperture) {
        report("nonseq.aperture_missing", sys(f + i).location,
               "surface '" + sys(f + i).id.str() + "' of a body has no aperture");
        ok = false;
      }
    }
    if (!ok) return;
    const std::uint32_t env = system_.environment_medium();
    for (std::uint32_t i = 0; i < n; ++i) {
      Face face;
      face.kind = FaceKind::Body;
      face.surface = copy_surface(f + i);
      face.system_surface = f + i;
      face.element = e;
      face.medium_positive = env;
      face.medium_negative = env;
      const auto assign = [&](double side, std::uint32_t medium) {
        (side > 0.0 ? face.medium_positive : face.medium_negative) = medium;
      };
      if (i >= 1) assign(up[i - 1], el.media[i - 1]);
      if (i + 1 < n) assign(low[i], el.media[i]);
      scene_index[i] = face.surface;
      add_face(face, surface_box(sys(f + i)));
    }
    for (std::uint32_t j = 0; j + 1 < n; ++j) {
      add_mantle(e, j, f + j, scene_index[j], scene_index[j + 1], low[j], up[j]);
    }
  }

  /// Synthetic scene surface of a generated face: the frame of surface j, plane shape, no
  /// aperture, the rim interaction (Fresnel, or Absorber per element; ADR 0057, point 2).
  std::uint32_t generated_surface(std::uint32_t e, std::uint32_t sj) {
    CompiledSurface s = sys(sj);
    s.shape = geom::Plane<double>{};
    s.aperture.reset();
    s.obscurations.clear();
    s.phases.clear();
    s.phase_functions.clear();
    s.diffraction_efficiency.reset();
    s.coating.reset();
    s.ideal_axis.reset();
    s.ideal_lens.reset();
    s.ideal_conductor = false;
    s.location = element_location(sys(sj));
    const bool absorbs = std::find(options_.absorbing_rims.begin(), options_.absorbing_rims.end(),
                                   e) != options_.absorbing_rims.end();
    s.interaction =
        absorbs ? model::Interaction{model::Absorber{}} : model::Interaction{model::Fresnel{}};
    scene_.surfaces_.push_back(std::move(s));
    return static_cast<std::uint32_t>(scene_.surfaces_.size() - 1);
  }

  /// Mantle of segment j between the scene surfaces lo and hi (ADR 0057, point 2).
  /// `sj` is the system index of surface j.
  void add_mantle(std::uint32_t e,
                  std::uint32_t j,
                  std::uint32_t sj,
                  std::uint32_t lo,
                  std::uint32_t hi,
                  double lo_side,
                  double hi_side) {
    const compile::CompiledElement& el = system_.elements()[e];
    // Copies: generated_surface() grows scene_.surfaces_.
    const CompiledSurface a = scene_.surfaces_[lo];
    const CompiledSurface b = scene_.surfaces_[hi];
    const std::string where = element_location(a);
    const std::string seg = "segment " + std::to_string(j) + " of element '" + el.name + "'";
    // Both apertures exist: add_segmented() reported nonseq.aperture_missing before.
    if (!a.aperture || !b.aperture) return;
    const model::Aperture& ap_a = *a.aperture;
    const model::Aperture& ap_b = *b.aperture;
    if (ap_a.index() != ap_b.index()) {
      report("nonseq.rim_unsupported", where, seg + ": the two apertures differ in kind");
      return;
    }
    const auto xy = aperture_xy(ap_a);
    const double size = std::max({std::abs(xy[0]), xy[1], std::abs(xy[2]), xy[3]});
    const Vec3 d = vertex(b) - vertex(a);
    const Vec3 lateral = d - d.dot(axis_z(a)) * axis_z(a);
    // The vertex coordinates enter through their rounding (ADR 0057, Nachtrag PR 2, h).
    const double scale = std::max({d.norm(), size, vertex(a).norm(), vertex(b).norm()});
    if ((axis_z(a) - axis_z(b)).norm() > kRimRelTol || lateral.norm() > kRimRelTol * scale) {
      report("nonseq.rim_unsupported", where,
             seg + ": the surface axes are not parallel and coaxial");
      return;
    }
    const bool circle = std::holds_alternative<model::CircularAperture>(ap_a);
    const auto equal = [&](double u, double v) {
      return std::abs(u - v) <= kRimRelTol * std::max(std::abs(u), std::abs(v));
    };
    if (!circle) {
      if ((axis_x(a) - axis_x(b)).norm() > kRimRelTol) {
        report("nonseq.rim_unsupported", where,
               seg + ": the apertures are turned against each other");
        return;
      }
      const auto pa = rim_polygon(ap_a);
      const auto pb = rim_polygon(ap_b);
      bool same = pa.size() == pb.size();
      for (std::size_t i = 0; same && i < pa.size(); ++i) {
        same = equal(pa[i][0], pb[i][0]) && equal(pa[i][1], pb[i][1]);
      }
      if (const auto* ea = std::get_if<model::EllipticalAperture>(&ap_a)) {
        const auto& eb = std::get<model::EllipticalAperture>(ap_b);
        same = equal(ea->semi_axis_x.value, eb.semi_axis_x.value) &&
               equal(ea->semi_axis_y.value, eb.semi_axis_y.value);
      }
      if (!same) {
        report("nonseq.rim_unsupported", where,
               seg + ": unequal non-circular apertures have no mantle");
        return;
      }
    }
    Face base;
    base.kind = FaceKind::Mantle;
    base.surface = generated_surface(e, sj);
    base.system_surface = sj;
    base.element = e;
    base.segment = j;
    base.lower = lo;
    base.upper = hi;
    base.lower_side = lo_side;
    base.upper_side = hi_side;
    base.medium_positive = system_.environment_medium();
    base.medium_negative = el.media[j];

    // z range in the frame of surface j of both surfaces over a rim curve (for the boxes); none
    // if the rim leaves the domain of a surface (the body would be open there).
    const auto rim_z = [&](const std::vector<std::array<double, 2>>& pts,
                           double h) -> std::optional<std::array<double, 2>> {
      double z0 = kInf;
      double z1 = -kInf;
      double g = 0.0;
      for (const auto& p : pts) {
        const double za = sag(a, p[0], p[1]);
        const double zb = height_in_frame(a, b, p[0], p[1]);
        if (!std::isfinite(za) || !std::isfinite(zb)) {
          report("nonseq.rim_unsupported", where,
                 seg + ": the rim lies outside the domain of a surface");
          return std::nullopt;
        }
        const auto [gx, gy] = grad(a, p[0], p[1]);
        z0 = std::min({z0, za, zb});
        z1 = std::max({z1, za, zb});
        g = std::max(g, std::hypot(gx, gy));
      }
      // The slope of surface j + 1 in frame j: estimated from its own gradient (frames parallel).
      for (const auto& p : pts) {
        const Vec3 q = b.to_local.apply_point(a.to_global.apply_point(Vec3(p[0], p[1], 0.0)));
        const auto [gx, gy] = grad(b, q.x(), q.y());
        g = std::max(g, std::hypot(gx, gy));
      }
      return std::array<double, 2>{z0 - h * g, z1 + h * g};
    };
    const auto wall_box = [&](double x0, double x1, double y0, double y1, double z0, double z1) {
      Aabb local;
      local.lo = Vec3(x0, y0, z0);
      local.hi = Vec3(x1, y1, z1);
      return with_margin(local);
    };

    if (!circle) {
      if (const auto* ea = std::get_if<model::EllipticalAperture>(&ap_a)) {
        Face f = base;
        f.generated.wall = WallShape::Ellipse;
        f.generated.a = ea->semi_axis_x.value;
        f.generated.b = ea->semi_axis_y.value;
        const double h =
            2.0 * std::numbers::pi * std::max(f.generated.a, f.generated.b) / kRimPoints;
        const auto z = rim_z(rim_points(ap_a, kRimPoints), h);
        if (!z) return;
        add_face(f, wall_box(-f.generated.a, f.generated.a, -f.generated.b, f.generated.b, (*z)[0],
                             (*z)[1]));
        return;
      }
      const auto poly = rim_polygon(ap_a);
      for (std::size_t i = 0; i < poly.size(); ++i) {
        Face f = base;
        f.generated.wall = WallShape::Edge;
        f.generated.p0 = poly[i];
        f.generated.p1 = poly[(i + 1) % poly.size()];
        std::vector<std::array<double, 2>> pts;
        constexpr int kEdgePoints = 64;
        for (int k = 0; k <= kEdgePoints; ++k) {
          const double u = static_cast<double>(k) / kEdgePoints;
          pts.push_back({f.generated.p0[0] + u * (f.generated.p1[0] - f.generated.p0[0]),
                         f.generated.p0[1] + u * (f.generated.p1[1] - f.generated.p0[1])});
        }
        const double len = std::hypot(f.generated.p1[0] - f.generated.p0[0],
                                      f.generated.p1[1] - f.generated.p0[1]);
        const auto z = rim_z(pts, len / kEdgePoints);
        if (!z) return;
        add_face(f, wall_box(std::min(f.generated.p0[0], f.generated.p1[0]),
                             std::max(f.generated.p0[0], f.generated.p1[0]),
                             std::min(f.generated.p0[1], f.generated.p1[1]),
                             std::max(f.generated.p0[1], f.generated.p1[1]), (*z)[0], (*z)[1]));
      }
      return;
    }

    // Circles: outer rim, and the bore of an inner radius (ADR 0057, point 2).
    const auto& ca = std::get<model::CircularAperture>(ap_a);
    const auto& cb = std::get<model::CircularAperture>(ap_b);
    const double ri_a = ca.inner_radius.value;
    const double ri_b = cb.inner_radius.value;
    if ((ri_a > 0.0) != (ri_b > 0.0)) {
      report("nonseq.rim_unsupported", where, seg + ": an inner radius on one surface only");
      return;
    }
    // Ring at radius r_a on surface j and r_b on surface j + 1; outward +1 for the outer rim
    // (the larger radius owns the wall), -1 for the bore (the smaller radius owns it).
    const auto ring = [&](double r_a, double r_b, double outward) {
      if (equal(r_a, r_b)) {
        Face f = base;
        f.generated.wall = WallShape::Circle;
        f.generated.a = std::max(r_a, r_b);
        f.generated.orientation = outward;
        const double r = f.generated.a;
        const auto z = rim_z(rim_points(model::CircularAperture{r, 0.0}, kRimPoints),
                             2.0 * std::numbers::pi * r / kRimPoints);
        if (!z) return;
        add_face(f, wall_box(-r, r, -r, r, (*z)[0], (*z)[1]));
        return;
      }
      if (!rotational(a) || !rotational(b)) {
        report("nonseq.rim_unsupported", where,
               seg + ": a step between unequal circles needs rotationally symmetric surfaces");
        return;
      }
      const bool a_owns = outward > 0.0 ? r_a > r_b : r_a < r_b;
      const double r_own = a_owns ? r_a : r_b;
      const double r_oth = a_owns ? r_b : r_a;
      const auto z_of = [&](bool on_a, double r) {
        return on_a ? sag(a, r, 0.0) : height_in_frame(a, b, r, 0.0);
      };
      const double h_own = z_of(a_owns, r_own);
      const double h_oth = z_of(!a_owns, r_oth);
      const double z_mid = z_of(a_owns, 0.5 * (r_own + r_oth));
      const double glass = sign_of(z_mid - h_oth);
      if (!std::isfinite(h_own) || !std::isfinite(h_oth) || !std::isfinite(z_mid) || glass == 0.0) {
        report("nonseq.rim_unsupported", where,
               seg + ": the rim lies outside the domain of a surface, or the shoulder is flat");
        return;
      }
      if (h_own != h_oth) {
        Face w = base;
        w.generated.wall = WallShape::Circle;
        w.generated.a = r_own;
        w.generated.orientation = outward;
        w.generated.constant_z = true;
        w.generated.z_lo = std::min(h_own, h_oth);
        w.generated.z_hi = std::max(h_own, h_oth);
        add_face(w, wall_box(-r_own, r_own, -r_own, r_own, w.generated.z_lo, w.generated.z_hi));
      }
      Face s = base;
      s.kind = FaceKind::Shoulder;
      s.generated.orientation = -glass;  // away from the glass between owner and ring
      s.generated.z_lo = h_oth;
      s.generated.z_hi = h_oth;
      s.generated.r_in = std::min(r_own, r_oth);
      s.generated.r_out = std::max(r_own, r_oth);
      const double r = s.generated.r_out;
      add_face(s, wall_box(-r, r, -r, r, h_oth, h_oth));
    };
    ring(ca.radius.value, cb.radius.value, 1.0);
    if (ri_a > 0.0) ring(ri_a, ri_b, -1.0);
  }

  /// Plate of one material with more than two surfaces (ADR 0057, point 2).
  void add_prism(std::uint32_t e) {
    const compile::CompiledElement& el = system_.elements()[e];
    const std::uint32_t f = el.first_surface;
    const std::uint32_t n = el.surface_count;
    bool ok = true;
    std::vector<std::vector<Vec3>> pts(n);
    Vec3 centroid = Vec3::Zero();
    std::size_t count = 0;
    for (std::uint32_t i = 0; i < n; ++i) {
      const CompiledSurface& s = sys(f + i);
      if (!s.aperture) {
        report("nonseq.aperture_missing", s.location,
               "surface '" + s.id.str() + "' of a body has no aperture");
        ok = false;
        continue;
      }
      for (const auto& c : corners(*s.aperture)) {
        pts[i].push_back(s.to_global.apply_point(Vec3(c[0], c[1], sag(s, c[0], c[1]))));
        centroid += pts[i].back();
        ++count;
      }
    }
    if (!ok) return;
    centroid /= static_cast<double>(count);
    std::vector<double> side(n, 0.0);
    for (std::uint32_t i = 0; i < n; ++i) {
      side[i] = sign_of(axis_z(sys(f + i)).dot(centroid - vertex(sys(f + i))));
      if (side[i] == 0.0) {
        report("nonseq.plate_not_convex", sys(f + i).location,
               "the centroid of the corners lies on surface '" + sys(f + i).id.str() + "'");
        return;
      }
    }
    for (std::uint32_t k = 0; k < n; ++k) {
      const CompiledSurface& plane = sys(f + k);
      if (!std::holds_alternative<geom::Plane<double>>(plane.shape)) continue;
      for (std::uint32_t i = 0; i < n; ++i) {
        if (i == k) continue;
        for (const Vec3& p : pts[i]) {
          if (side[k] * axis_z(plane).dot(p - vertex(plane)) < -trace::kApertureTolerance) {
            report("nonseq.plate_not_convex", plane.location,
                   "a corner of surface '" + sys(f + i).id.str() + "' lies outside the plane '" +
                       plane.id.str() + "' of element '" + el.name + "'");
            return;
          }
        }
      }
    }
    const std::uint32_t env = system_.environment_medium();
    for (std::uint32_t i = 0; i < n; ++i) {
      Face face;
      face.kind = FaceKind::Body;
      face.surface = copy_surface(f + i);
      face.system_surface = f + i;
      face.element = e;
      face.medium_positive = side[i] > 0.0 ? el.media[0] : env;
      face.medium_negative = side[i] > 0.0 ? env : el.media[0];
      add_face(face, surface_box(sys(f + i)));
    }
  }

  /// Mirror, ThinElement, Stop, Detector: thin faces (ADR 0057, point 2).
  void add_thin(std::uint32_t e) {
    const compile::CompiledElement& el = system_.elements()[e];
    const std::uint32_t env = system_.environment_medium();
    for (std::uint32_t i = 0; i < el.surface_count; ++i) {
      const std::uint32_t si = el.first_surface + i;
      Face face;
      face.system_surface = si;
      face.element = e;
      face.medium_positive = env;
      face.medium_negative = env;
      switch (el.kind) {
        case ElementKind::Mirror:
          face.kind = FaceKind::Mirror;
          if (!el.media.empty()) face.substrate = el.media[0];
          break;
        case ElementKind::Stop:
          face.kind = FaceKind::Stop;
          break;
        case ElementKind::Detector:
          face.kind = FaceKind::Detector;
          face.transparent = std::find(options_.transparent_detectors.begin(),
                                       options_.transparent_detectors.end(),
                                       e) != options_.transparent_detectors.end();
          break;
        default:
          face.kind = FaceKind::Thin;
          break;
      }
      if (face.kind == FaceKind::Stop) {
        const auto& aperture = sys(si).aperture;
        if (!aperture) continue;  // open everywhere (ADR 0057, Nachtrag PR 2, c)
        face.surface = copy_surface(si);
        CompiledSurface screen = sys(si);
        screen.aperture = scaled(*aperture, options_.screen_scale);
        screen.obscurations.clear();
        const Aabb box = surface_box(screen);
        scene_.surfaces_.push_back(std::move(screen));
        face.screen = static_cast<std::uint32_t>(scene_.surfaces_.size() - 1);
        add_face(face, box);
      } else {
        face.surface = copy_surface(si);
        add_face(face, surface_box(sys(si)));
      }
    }
  }

  void finish() {
    std::vector<Aabb> boxes;
    for (std::uint32_t f = 0; f < scene_.faces_.size(); ++f) {
      if (scene_.bounds_[f].finite()) {
        scene_.bounded_.push_back(f);
        boxes.push_back(scene_.bounds_[f]);
      } else {
        scene_.unbounded_.push_back(f);
      }
    }
    scene_.bvh_ = Bvh(boxes);
  }

  const compile::CompiledSystem& system_;
  const SceneOptions& options_;
  Scene scene_;
  std::vector<model::Diagnostic> errors_;
};

Scene Scene::build(const compile::CompiledSystem& system, const SceneOptions& options) {
  return SceneBuilder(system, options).run();
}

namespace {

/// Fills the global and media part of a hit on `face` from its local part.
SceneHit finish_hit(const Face& face,
                    const CompiledSurface& s,
                    std::uint32_t index,
                    const trace::SurfaceHit& local,
                    HitRegion region) noexcept {
  SceneHit h;
  h.status = geom::HitStatus::Hit;
  h.face = index;
  h.region = region;
  h.local = local;
  h.point = s.to_global.apply_point(local.point);
  h.normal = s.to_global.apply_vector(local.normal);
  h.from_positive = local.direction.dot(local.normal) < 0.0;
  h.medium_before = h.from_positive ? face.medium_positive : face.medium_negative;
  h.medium_beyond = h.from_positive ? face.medium_negative : face.medium_positive;
  return h;
}

}  // namespace

trace::SurfaceHit Scene::next_crossing(const trace::RayState& ray,
                                       std::uint32_t index,
                                       double t1,
                                       double t_max) const noexcept {
  const Face& face = faces_[index];
  const CompiledSurface& s = surfaces_[face.surface];
  const Aabb& box = local_bounds_[index];
  if (!box.finite()) return trace::intersect_surface(ray, s, t1);  // no interval to bracket
  const Vec3 o = s.to_local.apply_point(ray.pos);
  const Vec3 d = s.to_local.apply_vector(ray.dir);
  // Exit of the ray from the local box: no crossing of the face lies beyond it.
  double t_exit = kInf;
  for (int a = 0; a < 3; ++a) {
    if (d[a] > 0.0) t_exit = std::min(t_exit, (box.hi[a] - o[a]) / d[a]);
    if (d[a] < 0.0) t_exit = std::min(t_exit, (box.lo[a] - o[a]) / d[a]);
  }
  const double t_end = std::min(t_exit, t_max);
  trace::SurfaceHit none;
  none.direction = d;
  if (!(t_end > t1)) return none;
  // F(t) = z(t) - sag(x(t), y(t)) along the ray; a sign change between two samples brackets a
  // crossing, found by bisection (Press et al., Numerical Recipes, 3rd ed., Sec. 9.1), until
  // |F| <= kNewtonTolerance or the bracket reaches the rounding limit of geom::intersect.
  const auto f_at = [&](double t) {
    const Vec3 p = o + t * d;
    return p.z() - sag(s, p.x(), p.y());
  };
  const double h = (t_end - t1) / kBracketSamples;
  // The first sample sits just behind the crossing at t1 (F(t1) = 0), so that a near second
  // crossing in [t1, t1 + h] is bracketed as well; a NaN sample (outside the domain of the shape)
  // is skipped and the last finite one kept.
  double t_a = std::min(t1 + geom::kDefaultTMin, t1 + 0.5 * h);
  double f_a = f_at(t_a);
  for (int k = 1; k <= kBracketSamples; ++k) {
    const double t_b = k == kBracketSamples ? t_end : t1 + k * h;
    const double f_b = f_at(t_b);
    if (std::isfinite(f_a) && std::isfinite(f_b) && (f_a < 0.0) != (f_b < 0.0)) {
      double lo = t_a;
      double hi = t_b;
      double f_lo = f_a;
      double t = t_b;
      for (int i = 0; i < kMaxBisections; ++i) {
        t = 0.5 * (lo + hi);
        const double f = f_at(t);
        if (!std::isfinite(f)) {
          none.status = geom::HitStatus::NoConvergence;
          return none;
        }
        if (std::abs(f) <= geom::kNewtonTolerance ||
            hi - lo <= geom::kNewtonStepFactor * std::numeric_limits<double>::epsilon() *
                           (1.0 + std::abs(t) + o.norm())) {
          break;
        }
        if ((f < 0.0) == (f_lo < 0.0)) {
          lo = t;
          f_lo = f;
        } else {
          hi = t;
        }
      }
      trace::SurfaceHit hit;
      hit.status = geom::HitStatus::Hit;
      hit.t = t;
      hit.point = o + t * d;
      const auto [gx, gy] = grad(s, hit.point.x(), hit.point.y());
      hit.normal = Vec3(-gx, -gy, 1.0).normalized();  // as rtt-geom: +z at the vertex
      hit.direction = d;
      return hit;
    }
    if (std::isfinite(f_b) || !std::isfinite(f_a)) {
      t_a = t_b;
      f_a = f_b;
    }
  }
  return none;  // Missed: no further crossing inside the box
}

SceneHit Scene::intersect_face(const trace::RayState& ray,
                               std::uint32_t index,
                               double t_min,
                               double t_max) const noexcept {
  const Face& face = faces_[index];
  const CompiledSurface& s = surfaces_[face.surface];
  constexpr double tol = trace::kApertureTolerance;
  if (face.kind != FaceKind::Mantle && face.kind != FaceKind::Shoulder) {
    // NoConvergence on this face; `bound` is a lower bound of the t of a crossing it may hide.
    const auto no_convergence = [index](double bound) {
      SceneHit h;
      h.status = geom::HitStatus::NoConvergence;
      h.face = index;
      h.local.t = bound;
      return h;
    };
    // Where the face is there (ADR 0057, points 2 and 3; Nachtrag PR 2, e): a stop is open in its
    // opening outside its obscurations and absorbs elsewhere on its screen (inner radius and
    // obscurations included), and is not there beyond the screen; every other face is there
    // inside its aperture and absorbs in its obscurations.
    const auto region = [&](const trace::SurfaceHit& hit) -> std::optional<HitRegion> {
      const bool obscured = trace::in_obscuration(s, hit.point.x(), hit.point.y());
      if (face.kind == FaceKind::Stop) {
        if (trace::inside_aperture(s, hit) && !obscured) return std::nullopt;
        if (trace::inside_aperture(surfaces_[face.screen], hit)) return HitRegion::Absorbed;
        return std::nullopt;
      }
      if (!trace::inside_aperture(s, hit)) return std::nullopt;
      return obscured ? HitRegion::Absorbed : HitRegion::Clear;
    };
    trace::SurfaceHit hit = trace::intersect_surface(ray, s, t_min);
    for (int k = 0;; ++k) {
      if (hit.status == geom::HitStatus::Missed) return {};
      if (hit.status == geom::HitStatus::NoConvergence) return no_convergence(k == 0 ? 0.0 : hit.t);
      if (!(hit.t < t_max)) return {};
      if (const auto r = region(hit)) return finish_hit(face, s, index, hit, *r);
      if (k == kMaxCrossings) return no_convergence(hit.t);
      // The face is not there: its next crossing. Analytic shapes find it with t_min; a Newton
      // shape would restart at its base conic and fall back onto this crossing, so the next one
      // is bracketed instead (ADR 0057, Nachtrag PR 2, i).
      const double t_prev = hit.t;
      hit = analytic(s) ? trace::intersect_surface(ray, s, t_prev)
                        : next_crossing(ray, index, t_prev, t_max);
      if (hit.status == geom::HitStatus::NoConvergence) hit.t = t_prev;
    }
  }

  const Vec3 o = s.to_local.apply_point(ray.pos);
  const Vec3 d = s.to_local.apply_vector(ray.dir);
  const GeneratedGeometry& g = face.generated;
  // Pointwise bounds of a mantle: on the segment side of both surfaces (ADR 0057, point 2).
  const auto within = [&](const Vec3& p) {
    if (g.constant_z) return p.z() >= g.z_lo - tol && p.z() <= g.z_hi + tol;
    const CompiledSurface& lo = surfaces_[face.lower];
    const CompiledSurface& hi = surfaces_[face.upper];
    const Vec3 q = hi.to_local.apply_point(s.to_global.apply_point(p));
    return face.lower_side * (p.z() - sag(lo, p.x(), p.y())) >= -tol &&
           face.upper_side * (q.z() - sag(hi, q.x(), q.y())) >= -tol;
  };
  const auto accept = [&](double t, const Vec3& normal) -> SceneHit {
    trace::SurfaceHit hit;
    hit.status = geom::HitStatus::Hit;
    hit.t = t;
    hit.point = o + t * d;
    hit.normal = normal;
    hit.direction = d;
    return finish_hit(face, s, index, hit, HitRegion::Clear);
  };
  if (face.kind == FaceKind::Shoulder) {
    if (d.z() == 0.0) return {};
    const double t = (g.z_lo - o.z()) / d.z();
    if (!(t > t_min && t < t_max)) return {};
    const Vec3 p = o + t * d;
    const double r2 = p.x() * p.x() + p.y() * p.y();
    const double r_in = std::max(0.0, g.r_in - tol);
    const double r_out = g.r_out + tol;
    if (r2 < r_in * r_in || r2 > r_out * r_out) return {};
    return accept(t, Vec3(0.0, 0.0, g.orientation));
  }
  if (g.wall == WallShape::Edge) {
    const double ex = g.p1[0] - g.p0[0];
    const double ey = g.p1[1] - g.p0[1];
    const double len = std::hypot(ex, ey);
    const double nx = ey / len;  // outward normal of a counter-clockwise edge
    const double ny = -ex / len;
    const double denom = nx * d.x() + ny * d.y();
    if (denom == 0.0) return {};
    const double t = (nx * (g.p0[0] - o.x()) + ny * (g.p0[1] - o.y())) / denom;
    if (!(t > t_min && t < t_max)) return {};
    const Vec3 p = o + t * d;
    const double u = ((p.x() - g.p0[0]) * ex + (p.y() - g.p0[1]) * ey) / len;
    if (u < -tol || u > len + tol || !within(p)) return {};
    return accept(t, Vec3(nx, ny, 0.0));
  }
  // Circle or ellipse: |(x / a, y / b)|^2 = 1 along the ray, roots without cancellation (Press et
  // al., Numerical Recipes, 3rd ed., Sec. 5.6, Eqs. 5.6.4-5.6.5, as in rtt::geom::intersect_conic).
  const double a = g.a;
  const double b = g.wall == WallShape::Circle ? g.a : g.b;
  const double ox = o.x() / a;
  const double oy = o.y() / b;
  const double dx = d.x() / a;
  const double dy = d.y() / b;
  const double qa = dx * dx + dy * dy;
  if (qa == 0.0) return {};
  const double qb = ox * dx + oy * dy;
  const double qc = ox * ox + oy * oy - 1.0;
  const double disc = qb * qb - qa * qc;
  if (!(disc >= 0.0)) return {};
  const double q = -(qb + (qb < 0.0 ? -std::sqrt(disc) : std::sqrt(disc)));
  std::array<double, 2> roots{q / qa, q != 0.0 ? qc / q : q / qa};
  if (roots[1] < roots[0]) std::swap(roots[0], roots[1]);
  for (const double t : roots) {
    if (!(t > t_min && t < t_max)) continue;
    const Vec3 p = o + t * d;
    if (!within(p)) continue;
    const Vec3 n = Vec3(p.x() / (a * a), p.y() / (b * b), 0.0).normalized() * g.orientation;
    return accept(t, n);
  }
  return {};
}

SceneHit Scene::closest_hit(const trace::RayState& ray, std::uint32_t last_face) const noexcept {
  SceneHit best;
  double best_t = kInf;
  double noconv_t = kInf;
  std::uint32_t noconv_face = kNoFace;
  const SlabRay slab = SlabRay::from(ray.pos, ray.dir);
  const auto test = [&](std::uint32_t face) {
    const double t_min = face == last_face ? geom::kDefaultTMin : 0.0;
    SceneHit h = intersect_face(ray, face, t_min, best_t);
    if (h.status == geom::HitStatus::Hit) {
      best_t = h.local.t;
      best = std::move(h);
    } else if (h.status == geom::HitStatus::NoConvergence) {
      // A hidden crossing lies inside the face's own box and beyond the bound of the search.
      double t_box = 0.0;
      if (bounds_[face].finite() && !slab_test(bounds_[face], slab, kInf, t_box)) t_box = kInf;
      const double bound = std::max(t_box, h.local.t);
      if (bound < noconv_t) {
        noconv_t = bound;
        noconv_face = face;
      }
    }
    return best_t;
  };
  for (const std::uint32_t f : unbounded_) test(f);
  bvh_.traverse(slab, kInf, [&](std::uint32_t prim, double) { return test(bounded_[prim]); });
  // A face whose Newton intersection failed may hold a nearer hit, unless that hit could only lie
  // beyond the best one (ADR 0057, Nachtrag PR 1, e, and Nachtrag PR 2, g).
  if (noconv_face != kNoFace && noconv_t <= best_t) {
    SceneHit h;
    h.status = geom::HitStatus::NoConvergence;
    h.face = noconv_face;
    return h;
  }
  return best;
}

}  // namespace rtt::nonseq
