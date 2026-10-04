#include "rtt/trace/sources.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <optional>
#include <random>
#include <stdexcept>
#include <string>
#include <type_traits>

#include "rtt/paraxial/paraxial.hpp"

namespace rtt::trace {

namespace {

using compile::CompiledSystem;
using compile::PathId;
using math::Vec3;

constexpr double kTwoPi = 2.0 * std::numbers::pi;

/// Maximum number of step halvings per Newton step of real aiming.
constexpr int kMaxAimHalvings = 10;

/// Draw in [0, 1) from the upper 53 bits of a 64-bit engine output; identical on all platforms
/// (std::uniform_real_distribution is implementation-defined).
double unit_draw(std::mt19937_64& engine) {
  return static_cast<double>(engine() >> 11) * 0x1.0p-53;
}

/// n points evenly spaced on [-1, 1]; n = 1 gives {0}.
std::vector<double> even_points(int n) {
  if (n < 1) throw std::invalid_argument("pupil sampling: n must be >= 1");
  if (n == 1) return {0.0};
  std::vector<double> v(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    v[static_cast<std::size_t>(i)] = -1.0 + 2.0 * static_cast<double>(i) / (n - 1);
  }
  return v;
}

/// Everything that is the same for all rays of one path and wavelength.
struct Context {
  const CompiledSystem* system = nullptr;
  PathId path;
  std::uint16_t wavelength = 0;
  paraxial::FirstOrder first_order;
  std::size_t stop_event = 0;  ///< index of the first Stop event in the path
  double z_ep = 0.0;           ///< global z of the paraxial entrance pupil, mm
  double r_ep = 0.0;           ///< radius of the paraxial entrance pupil, mm
  /// Signed paraxial stop height of a ray through the EP point at unit height: the stop target
  /// of pupil point (px, py) is (px, py) * r_ep * stop_scale (negative for an inverted image).
  double stop_scale = 1.0;
  double z_min = 0.0;  ///< start limit: rays of an infinite object start at least 1 mm before it
};

/// Ray start for one field.
struct FieldStart {
  bool infinite = true;
  Vec3 direction = Vec3::UnitZ();  ///< object at infinity: field direction
  double plane_offset = 0.0;       ///< object at infinity: L of the start plane, mm
  Vec3 object = Vec3::Zero();      ///< finite object: object point, mm
};

Context make_context(const CompiledSystem& system, PathId path, std::uint16_t wavelength) {
  if (path.index >= system.paths().size()) {
    throw std::invalid_argument("sources: path index " + std::to_string(path.index) +
                                " does not exist");
  }
  if (wavelength >= system.wavelengths_um().size()) {
    throw std::invalid_argument("sources: wavelength index " + std::to_string(wavelength) +
                                " does not exist");
  }
  Context c;
  c.system = &system;
  c.path = path;
  c.wavelength = wavelength;
  const auto& events = system.path(path).events;
  const auto stop = std::find_if(events.begin(), events.end(), [&](const auto& e) {
    return system.surfaces()[e.surface].element_kind == model::ElementKind::Stop;
  });
  if (stop == events.end()) {
    throw std::invalid_argument("sources: the path has no stop to aim at");
  }
  c.stop_event = static_cast<std::size_t>(stop - events.begin());
  c.first_order = paraxial::first_order(system, path, wavelength);
  const auto& ep = c.first_order.entrance_pupil;
  if (!ep || !ep->z || !ep->diameter || !(*ep->diameter > 0.0)) {
    throw std::invalid_argument(
        "sources: the entrance pupil is not defined (pupil at infinity, or the system aperture "
        "gives no diameter for this object)");
  }
  c.z_ep = *ep->z;
  c.r_ep = *ep->diameter / 2.0;
  // The EP plane and the stop plane are conjugate, so the stop height of a paraxial ray from
  // the EP plane does not depend on its slope: y_stop = stop_scale * y_ep.
  c.stop_scale = paraxial::trace_ray(system, path, wavelength, c.z_ep, 1.0, 0.0)[c.stop_event].y;

  // Lowest z the bundle may reach before the system: every surface of the path is bounded by a
  // sphere of radius rho around its vertex, rho = outer extent of its aperture (or 3 EP radii
  // without aperture); this holds for planes, spheres and ellipsoids within their aperture.
  double z_min = c.z_ep;
  for (const auto& e : events) {
    const compile::CompiledSurface& s = system.surfaces()[e.surface];
    double rho = 3.0 * c.r_ep;
    if (s.aperture) {
      rho = std::visit(
          [](const auto& a) {
            using A = std::decay_t<decltype(a)>;
            if constexpr (std::is_same_v<A, model::CircularAperture>) {
              return a.radius;
            } else if constexpr (std::is_same_v<A, model::RectangularAperture>) {
              return std::hypot(a.half_width_x, a.half_width_y);
            } else {
              return std::max(a.semi_axis_x, a.semi_axis_y);
            }
          },
          *s.aperture);
    }
    z_min = std::min(z_min, s.to_global.translation().z() - rho);
  }
  c.z_min = z_min - 1.0;
  return c;
}

/// Paraxial image height (at the paraxial image plane) of the chief ray with unit field value:
/// unit slope through the EP centre (object at infinity) or unit object height (finite object).
double unit_image_height(const Context& c) {
  const auto& fo = c.first_order;
  if (!fo.image_z) {
    throw std::invalid_argument("sources: paraxial image height needs a finite paraxial image");
  }
  const CompiledSystem& system = *c.system;
  std::vector<paraxial::RayAtEvent> ray;
  if (system.object().at_infinity) {
    ray = paraxial::trace_ray(system, c.path, c.wavelength, c.z_ep, 0.0, 1.0);
  } else {
    const double z_obj = -system.object().distance.value;
    ray = paraxial::trace_ray(system, c.path, c.wavelength, z_obj, 1.0, -1.0 / (c.z_ep - z_obj));
  }
  const auto& last = ray.back();
  const double y = last.y + (*fo.image_z - last.z) * last.u;
  if (y == 0.0) throw std::invalid_argument("sources: chief ray does not reach the image");
  return y;
}

FieldStart make_field(const Context& c, std::uint16_t field) {
  const CompiledSystem& system = *c.system;
  const auto& points = system.fields().points;
  if (field >= points.size()) {
    throw std::invalid_argument("sources: field index " + std::to_string(field) +
                                " does not exist");
  }
  const model::Field& f = points[field];
  const bool infinite = system.object().at_infinity;
  const double z_obj = infinite ? 0.0 : -system.object().distance.value;

  // Field value as tangent of the angle (infinite object) or object height (finite object).
  double tx = 0.0;
  double ty = 0.0;
  double hx = 0.0;
  double hy = 0.0;
  switch (system.fields().type) {
    case model::FieldType::AngleDeg: {
      if (!(std::abs(f.x) < 90.0 && std::abs(f.y) < 90.0)) {
        throw std::invalid_argument("sources: field angles must lie in (-90, 90) degree");
      }
      tx = std::tan(f.x * std::numbers::pi / 180.0);
      ty = std::tan(f.y * std::numbers::pi / 180.0);
      // Finite object: the object point on the chief ray through the EP centre.
      hx = (z_obj - c.z_ep) * tx;
      hy = (z_obj - c.z_ep) * ty;
      break;
    }
    case model::FieldType::ObjectHeight:
      if (infinite) {
        throw std::invalid_argument("sources: object height needs a finite object distance");
      }
      hx = f.x;
      hy = f.y;
      break;
    case model::FieldType::ParaxialImageHeight: {
      // Linear in the field value: scale the unit chief ray (decided for #8).
      const double unit = unit_image_height(c);
      tx = f.x / unit;
      ty = f.y / unit;
      hx = tx;
      hy = ty;
      break;
    }
  }

  FieldStart start;
  start.infinite = infinite;
  if (!infinite) {
    if (c.z_ep == z_obj) throw std::invalid_argument("sources: entrance pupil in the object");
    start.object = Vec3(hx, hy, z_obj);
    return start;
  }
  start.direction = Vec3(tx, ty, 1.0).normalized();
  // Start plane perpendicular to d through E - L d (E = EP centre). A ray through (a, b, z_ep)
  // starts at z_ep - (a d_x + b d_y + L) d_z; with |(a, b)| <= 3 r_ep this is <= z_min for
  // L >= (z_ep - z_min) / d_z + 3 r_ep |d_xy| (decided for #8, same L for all rays of a field).
  const Vec3& d = start.direction;
  start.plane_offset = (c.z_ep - c.z_min) / d.z() + 3.0 * c.r_ep * std::hypot(d.x(), d.y());
  return start;
}

/// Ray through the point (a, b, z_ep) on the EP plane.
RayState ray_through(const Context& c, const FieldStart& f, double a, double b) {
  const Vec3 q(a, b, c.z_ep);
  RayState ray;
  if (f.infinite) {
    const Vec3 e(0.0, 0.0, c.z_ep);
    const double s = f.direction.dot(q - e) + f.plane_offset;
    ray.pos = q - s * f.direction;
    ray.dir = f.direction;
  } else {
    ray.pos = f.object;
    ray.dir = (q - f.object).normalized();
  }
  return ray;
}

/// Local (x, y) of the hit on the stop surface, apertures ignored; nullopt if the ray does not
/// get there.
std::optional<std::pair<double, double>> stop_hit(const Context& c, RayState ray) {
  const CompiledSystem& system = *c.system;
  const auto& events = system.path(c.path).events;
  for (std::size_t i = 0; i < c.stop_event; ++i) {
    const auto& e = events[i];
    const compile::CompiledSurface& s = system.surfaces()[e.surface];
    const SurfaceHit hit = intersect_surface(ray, s);
    ray = apply_event(ray, s, hit, e.surface, e.kind,
                      system.media()[e.medium_before].index[c.wavelength].real(),
                      system.media()[e.medium_after].index[c.wavelength].real());
    if (ray.status != RayStatus::Alive) return std::nullopt;
  }
  const SurfaceHit hit = intersect_surface(ray, system.surfaces()[events[c.stop_event].surface]);
  if (hit.status != geom::HitStatus::Hit) return std::nullopt;
  return std::pair{hit.point.x(), hit.point.y()};
}

AimedRay aim(const Context& c, const FieldStart& f, double px, double py, Aiming aiming) {
  // Paraxial solution: straight through the EP point (px, py) * r_ep.
  double a = px * c.r_ep;
  double b = py * c.r_ep;
  AimedRay out;
  out.ray = ray_through(c, f, a, b);
  if (aiming == Aiming::Paraxial) return out;

  // Real aiming: Newton on the residual at the stop, central-difference Jacobian (ADR 0007;
  // Newton's method in 2D, e.g. Press et al., Numerical Recipes, 3rd ed., Sec. 9.6).
  const double r_s = c.r_ep * c.stop_scale;
  const double tx = px * r_s;
  const double ty = py * r_s;
  const double h = kAimStepRelative * std::abs(r_s);
  const auto residual = [&](double u, double v) -> std::optional<std::pair<double, double>> {
    const auto hit = stop_hit(c, ray_through(c, f, u, v));
    if (!hit) return std::nullopt;
    return std::pair{hit->first - tx, hit->second - ty};
  };
  for (int iteration = 0;; ++iteration) {
    out.iterations = iteration;
    const auto r = residual(a, b);
    if (!r) break;
    out.residual = std::hypot(r->first, r->second);
    if (out.residual < kAimTolerance) {
      out.ray = ray_through(c, f, a, b);
      return out;
    }
    if (iteration == kMaxAimIterations) break;
    const auto ap = residual(a + h, b);
    const auto am = residual(a - h, b);
    const auto bp = residual(a, b + h);
    const auto bm = residual(a, b - h);
    if (!ap || !am || !bp || !bm) break;
    const double j11 = (ap->first - am->first) / (2.0 * h);
    const double j21 = (ap->second - am->second) / (2.0 * h);
    const double j12 = (bp->first - bm->first) / (2.0 * h);
    const double j22 = (bp->second - bm->second) / (2.0 * h);
    const double det = j11 * j22 - j12 * j21;
    if (!(det != 0.0) || !std::isfinite(det)) break;
    const double da = -(j22 * r->first - j12 * r->second) / det;
    const double db = -(-j21 * r->first + j11 * r->second) / det;
    // Backtracking: halve the Newton step while it does not reduce the residual or the ray
    // does not reach the stop (damped Newton, cf. Numerical Recipes, Sec. 9.7).
    double lambda = 1.0;
    bool accepted = false;
    for (int k = 0; k < kMaxAimHalvings; ++k, lambda *= 0.5) {
      const auto trial = residual(a + lambda * da, b + lambda * db);
      if (trial && std::hypot(trial->first, trial->second) < out.residual) {
        accepted = true;
        break;
      }
    }
    if (!accepted) break;
    a += lambda * da;
    b += lambda * db;
    out.ray = ray_through(c, f, a, b);
    if (!out.ray.pos.allFinite() || !out.ray.dir.allFinite()) break;
  }
  out.ray.status = RayStatus::NoConvergence;
  if (!out.ray.pos.allFinite() || !out.ray.dir.allFinite()) {
    out.ray = ray_through(c, f, px * c.r_ep, py * c.r_ep);
    out.ray.status = RayStatus::NoConvergence;
  }
  return out;
}

}  // namespace

std::vector<PupilPoint> pupil_points(const PupilSampling& sampling) {
  return std::visit(
      [](const auto& s) -> std::vector<PupilPoint> {
        using S = std::decay_t<decltype(s)>;
        std::vector<PupilPoint> p;
        if constexpr (std::is_same_v<S, SinglePupilPoint>) {
          p.push_back({s.px, s.py});
        } else if constexpr (std::is_same_v<S, HexapolarPupil>) {
          if (s.rings < 0) throw std::invalid_argument("pupil sampling: rings must be >= 0");
          p.push_back({0.0, 0.0});
          for (int k = 1; k <= s.rings; ++k) {
            const double r = static_cast<double>(k) / s.rings;
            for (int j = 0; j < 6 * k; ++j) {
              const double phi = kTwoPi * j / (6.0 * k);
              p.push_back({r * std::sin(phi), r * std::cos(phi)});
            }
          }
        } else if constexpr (std::is_same_v<S, GridPupil>) {
          const auto v = even_points(s.n);
          for (const double y : v) {
            for (const double x : v) {
              if (x * x + y * y <= 1.0 + 1e-12) p.push_back({x, y});
            }
          }
        } else if constexpr (std::is_same_v<S, FanXPupil>) {
          for (const double x : even_points(s.n)) p.push_back({x, 0.0});
        } else if constexpr (std::is_same_v<S, FanYPupil>) {
          for (const double y : even_points(s.n)) p.push_back({0.0, y});
        } else {
          std::mt19937_64 engine(s.seed);
          p.reserve(s.count);
          for (std::size_t i = 0; i < s.count; ++i) {
            const double u1 = unit_draw(engine);
            const double u2 = unit_draw(engine);
            const double r = std::sqrt(u1);
            const double phi = kTwoPi * u2;
            p.push_back({r * std::sin(phi), r * std::cos(phi)});
          }
        }
        return p;
      },
      sampling);
}

AimedRay aim_ray(const compile::CompiledSystem& system,
                 compile::PathId path,
                 std::uint16_t field,
                 std::uint16_t wavelength,
                 double px,
                 double py,
                 Aiming aiming) {
  const Context c = make_context(system, path, wavelength);
  return aim(c, make_field(c, field), px, py, aiming);
}

RayBatch make_rays(const compile::CompiledSystem& system,
                   compile::PathId path,
                   std::span<const std::uint16_t> fields,
                   std::uint16_t wavelength,
                   const PupilSampling& sampling,
                   Aiming aiming) {
  const Context c = make_context(system, path, wavelength);
  const std::vector<PupilPoint> points = pupil_points(sampling);
  std::vector<FieldStart> starts;
  starts.reserve(fields.size());
  for (const std::uint16_t f : fields) starts.push_back(make_field(c, f));

  RayBatch rays(fields.size() * points.size());
  std::size_t i = 0;
  for (std::size_t f = 0; f < fields.size(); ++f) {
    for (const PupilPoint& p : points) {
      const AimedRay aimed = aim(c, starts[f], p.px, p.py, aiming);
      rays.pos_x()[i] = aimed.ray.pos.x();
      rays.pos_y()[i] = aimed.ray.pos.y();
      rays.pos_z()[i] = aimed.ray.pos.z();
      rays.dir_x()[i] = aimed.ray.dir.x();
      rays.dir_y()[i] = aimed.ray.dir.y();
      rays.dir_z()[i] = aimed.ray.dir.z();
      rays.wl()[i] = wavelength;
      rays.field()[i] = fields[f];
      rays.pupil_x()[i] = p.px;
      rays.pupil_y()[i] = p.py;
      rays.status()[i] = aimed.ray.status;
      ++i;
    }
  }
  return rays;
}

}  // namespace rtt::trace
