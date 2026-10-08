#include "rtt/trace/sources.hpp"

#include <oneapi/tbb/blocked_range.h>
#include <oneapi/tbb/parallel_for.h>
#include <oneapi/tbb/partitioner.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <optional>
#include <random>
#include <stdexcept>
#include <string>
#include <type_traits>

#include "rtt/compile/errors.hpp"
#include "rtt/paraxial/paraxial.hpp"

namespace rtt::trace {

namespace {

using compile::CompiledSystem;
using compile::PathId;
using math::Vec3;

constexpr double kTwoPi = 2.0 * std::numbers::pi;

/// Maximum number of step halvings per Newton step of real aiming.
constexpr int kMaxAimHalvings = 10;

/// Half-width of the bundle of one field in units of the EP radius, around the chief ray
/// through the EP centre: reserve for the pupil aberration found by real aiming.
constexpr double kBundleReserve = 3.0;

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
  /// Object-space telecentric (#96): finite object and the paraxial entrance pupil at infinity
  /// at the ray's wavelength. Pupil coordinates are then object-space slopes, not EP points.
  bool telecentric = false;
  double z_ep = 0.0;  ///< global z of the paraxial entrance pupil, mm (unused if telecentric)
  /// Radius of the paraxial entrance pupil, mm; if telecentric, the paraxial marginal slope u_m
  /// in object space (pupil coordinates are slopes).
  double r_ep = 0.0;
  /// Signed paraxial stop height of a ray through the EP point at unit height (telecentric: of
  /// the ray from the axial object point with unit slope): the stop target of pupil point
  /// (px, py) is (px, py) * r_ep * stop_scale (negative for an inverted image).
  double stop_scale = 1.0;
  /// First-order data for the field conversion, so a field point is the same physical direction
  /// or object point for every wavelength (#31). Set at the reference wavelength only for the
  /// field types that need it (paraxial image height, angle with a finite object); otherwise a
  /// copy of the ray-wavelength data that make_field does not read.
  paraxial::FirstOrder first_order_ref;
  double z_ep_ref = 0.0;         ///< global z of the paraxial entrance pupil for the conversion, mm
  bool telecentric_ref = false;  ///< EP at infinity for the conversion (finite object, #96)
  /// Orientation of the pupil labels at the reference wavelength (#96): -1 if the entrance
  /// pupil of this wavelength lies on the other side of the object than that of the reference
  /// wavelength (finite object), else +1. Pupil point (px, py) is then the EP point
  /// label_sign * (px, py) * r_ep, so it means the same side of the stop for every wavelength.
  double label_sign = 1.0;
};

/// Ray start for one field.
struct FieldStart {
  bool infinite = true;
  Vec3 direction = Vec3::UnitZ();  ///< object at infinity: field direction
  double plane_offset = 0.0;       ///< object at infinity: L of the start plane, mm
  Vec3 object = Vec3::Zero();      ///< finite object: object point, mm
};

/// Object-space telecentric context (#96): finite object, paraxial entrance pupil at infinity
/// (the front matrix [[a, b], [c, d]] from the first vertex to the stop has a = 0; y-nu
/// equations of Greivenkamp, OPTI-201/202 lecture notes, Sec. 9, p. 9-2, matrix form as in
/// rtt-paraxial first_order). A paraxial
/// ray then reaches the stop at y_s = a y + b n1 u = b n1 u, whatever its height: the stop point
/// depends on the object-space slope only, so pupil coordinates are slopes, d ~ (a, b, 1) from
/// the object point, and the chief ray is parallel to the axis (decided for #96).
/// - stop_scale = s, the stop height of the paraxial ray from the axial object point with unit
///   slope (= b n1; s != 0 because a d - b c = 1).
/// - r_ep = u_m, the paraxial marginal slope: NA / n1 for object_na (the paraxial reading of
///   first_order), r_stop / |s| for stop_size. An entrance pupil diameter or an image-space
///   F-number (EPD = EFL / F#) does not define a bundle for a pupil at infinity: error.
void set_telecentric(Context& c) {
  const CompiledSystem& system = *c.system;
  c.telecentric = true;
  const double z_obj = -system.object().distance.value;
  c.stop_scale = paraxial::trace_ray(system, c.path, c.wavelength, z_obj, 0.0, 1.0)[c.stop_event].y;
  if (c.stop_scale == 0.0) {
    throw std::invalid_argument("sources: telecentric object space, but no ray reaches the stop");
  }
  const double value = system.aperture().value.value;
  switch (system.aperture().type) {
    case model::SystemApertureType::ObjectSpaceNA:
      c.r_ep = value / c.first_order.object_index;
      break;
    case model::SystemApertureType::StopSize: {
      const auto& events = system.path(c.path).events;
      const auto& aperture = system.surfaces()[events[c.stop_event].surface].aperture;
      // first_order() has checked that the stop aperture is circular.
      const auto* circle = aperture ? std::get_if<model::CircularAperture>(&*aperture) : nullptr;
      if (circle == nullptr) {
        throw std::invalid_argument("sources: the stop needs a circular aperture");
      }
      c.r_ep = circle->radius / std::abs(c.stop_scale);
      break;
    }
    case model::SystemApertureType::EntrancePupilDiameter:
    case model::SystemApertureType::ImageSpaceFNumber:
      throw std::invalid_argument(
          "sources: the entrance pupil is at infinity (object-space telecentric); give the "
          "system aperture as object_na or stop_size");
  }
  if (!(c.r_ep > 0.0)) {
    throw std::invalid_argument("sources: the system aperture gives no bundle (u_m <= 0)");
  }
}

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
  compile::require_stop(system, path);  // the one check for all users of the stop (ADR 0022)
  const auto& events = system.path(path).events;
  const auto stop = std::find_if(events.begin(), events.end(), [&](const auto& e) {
    return system.surfaces()[e.surface].element_kind == model::ElementKind::Stop;
  });
  c.stop_event = static_cast<std::size_t>(stop - events.begin());
  c.first_order = paraxial::first_order(system, path, wavelength);
  const auto& ep = c.first_order.entrance_pupil;
  const bool infinite = system.object().at_infinity;
  if (ep && !ep->z && !infinite) {
    set_telecentric(c);
  } else {
    if (!ep || !ep->z || !ep->diameter || !(*ep->diameter > 0.0)) {
      throw std::invalid_argument(
          "sources: the entrance pupil is not defined (pupil at infinity with the object at "
          "infinity, or the system aperture gives no diameter for this object)");
    }
    c.z_ep = *ep->z;
    c.r_ep = *ep->diameter / 2.0;
    // The EP plane and the stop plane are conjugate, so the stop height of a paraxial ray from
    // the EP plane does not depend on its slope: y_stop = stop_scale * y_ep.
    c.stop_scale = paraxial::trace_ray(system, path, wavelength, c.z_ep, 1.0, 0.0)[c.stop_event].y;
  }
  // Field conversion at the reference wavelength (decided for #31, fix of #8). Only paraxial
  // image heights and angles with a finite object need paraxial data for the conversion; other
  // field types skip the extra first-order computation.
  const bool needs_ref =
      system.fields().type == model::FieldType::ParaxialImageHeight ||
      (system.fields().type == model::FieldType::AngleDeg && !system.object().at_infinity);
  if (wavelength == system.reference_wavelength() || !needs_ref) {
    c.first_order_ref = c.first_order;
    c.z_ep_ref = c.z_ep;
    c.telecentric_ref = c.telecentric;
  } else {
    c.first_order_ref = paraxial::first_order(system, path, system.reference_wavelength());
    const auto& ep_ref = c.first_order_ref.entrance_pupil;
    if (ep_ref && !ep_ref->z && !infinite) {
      c.telecentric_ref = true;
    } else if (!ep_ref || !ep_ref->z) {
      throw std::invalid_argument(
          "sources: the entrance pupil at the reference wavelength is at infinity");
    } else {
      c.z_ep_ref = *ep_ref->z;
    }
  }
  // Pupil labels oriented at the reference wavelength (decided for #96, refining #93): the
  // image of the pupil turns over when the EP passes through infinity, so a wavelength whose
  // EP lies on the other side of the object (z_ep - z_obj of the other sign) would otherwise
  // see (px, py) mirrored on the stop. sigma = side(ref) side(lambda), side = sign(z_ep - z_obj),
  // and an EP at infinity (object-space telecentric) counts as the + side (py > 0 = slope up).
  // An object at infinity has no such side: sigma = +1. Without a change of side over the
  // wavelengths sigma = +1 and the rays are bit for bit the same as before.
  if (!infinite && wavelength != system.reference_wavelength()) {
    const double z_obj = -system.object().distance.value;
    const auto side = [z_obj](const paraxial::FirstOrder& fo) {
      const auto& e = fo.entrance_pupil;
      if (!e || !e->z) return 1.0;
      return *e->z < z_obj ? -1.0 : 1.0;
    };
    const paraxial::FirstOrder fo_ref =
        needs_ref ? c.first_order_ref
                  : paraxial::first_order(system, path, system.reference_wavelength());
    c.label_sign = side(fo_ref) * side(c.first_order);
  }
  return c;
}

/// Outer radius of an aperture in local x, y, mm.
double aperture_extent(const model::Aperture& aperture) {
  return std::visit(
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
      aperture);
}

/// Lowest global z of the part of `surface` that the bundle of an infinite-object field can
/// reach. The bundle is bounded by a cylinder of radius r_bundle around the chief ray through
/// the EP centre (slope tan_field against the axis), so at global z it is at most
/// r_bundle + |z - z_ep| tan_field off the axis. The surface is sampled in local coordinates up
/// to that radius, capped by its aperture and by the domain of its shape; since a surface that
/// curves back upstream meets the bundle where it is wider, the radius is re-evaluated at the
/// lowest z found until it no longer drops (decided for #8, start plane).
/// Valid for the rotationally symmetric paths that rtt::paraxial::first_order accepts (the
/// sampled disk is centred on the vertex). The bound drops monotonically. With a finite cap
/// (outer radius of the aperture, or shape domain) it is bounded below by the lowest z of the
/// whole cap disk (up to the sampling, hence the final std::min), but it may converge only
/// linearly (factor tan_field |dsag/dr| close to 1); if it has not settled after kMaxRounds and
/// the cap is below kDivergentRadius, the lowest z of the whole cap disk is returned, a safe
/// lower bound since the bundle cannot meet the surface outside its cap (decided for #72).
/// Without such a cap the bound can drop without limit (an unbounded surface without aperture,
/// or with an aperture of 1 km or more, that curves back against a steep field, e.g. a concave
/// paraboloid); then there is no safe start plane: std::invalid_argument, raised before any ray
/// is traced.
double lowest_z(const compile::CompiledSurface& surface,
                double z_ep,
                double r_bundle,
                double tan_field) {
  constexpr int kRadialSamples = 64;
  constexpr int kAzimuthSamples = 8;
  constexpr int kMaxRounds = 50;
  constexpr double kDivergentRadius = 1e6;  // mm
  double cap = std::numeric_limits<double>::infinity();
  if (surface.aperture) cap = aperture_extent(*surface.aperture);
  const std::optional<double> domain =
      std::visit([](const auto& shape) { return shape.max_radius(); }, surface.shape);
  if (domain) cap = std::min(cap, *domain);
  const double z_vertex = surface.to_global.translation().z();
  // Lowest global z of the surface on the local disk of radius rho around the vertex.
  const auto lowest_on_disk = [&](double rho) {
    double low = z_vertex;
    for (int i = 1; i <= kRadialSamples; ++i) {
      const double r = rho * i / kRadialSamples;
      for (int j = 0; j < kAzimuthSamples; ++j) {
        const double phi = kTwoPi * j / kAzimuthSamples;
        const double x = r * std::sin(phi);
        const double y = r * std::cos(phi);
        // NaN only at the rim of the domain (rounding); skipped.
        const double h =
            std::visit([x, y](const auto& shape) { return shape.sag(x, y); }, surface.shape);
        if (!std::isfinite(h)) continue;
        low = std::min(low, surface.to_global.apply_point(Vec3(x, y, h)).z());
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
  throw std::invalid_argument("sources: no start plane before surface '" + surface.id.str() +
                              "' for this field (unbounded surface without aperture, or with an "
                              "aperture of 1 km or more, that curves back against the field); "
                              "give it an aperture of optical size");
}

/// Paraxial image height (at the paraxial image plane) of the chief ray with unit field value:
/// unit slope through the EP centre (object at infinity) or unit object height (finite object),
/// at the reference wavelength (#31).
double unit_image_height(const Context& c) {
  const auto& fo = c.first_order_ref;
  if (!fo.image_z) {
    throw std::invalid_argument("sources: paraxial image height needs a finite paraxial image");
  }
  const CompiledSystem& system = *c.system;
  const std::uint16_t ref = system.reference_wavelength();
  std::vector<paraxial::RayAtEvent> ray;
  if (system.object().at_infinity) {
    ray = paraxial::trace_ray(system, c.path, ref, c.z_ep_ref, 0.0, 1.0);
  } else {
    const double z_obj = -system.object().distance.value;
    if (c.telecentric_ref) {
      // EP at infinity: the chief ray of unit object height is parallel to the axis (#96).
      ray = paraxial::trace_ray(system, c.path, ref, z_obj, 1.0, 0.0);
    } else {
      if (c.z_ep_ref == z_obj) {
        throw std::invalid_argument("sources: entrance pupil in the object (reference wavelength)");
      }
      ray = paraxial::trace_ray(system, c.path, ref, z_obj, 1.0, -1.0 / (c.z_ep_ref - z_obj));
    }
  }
  const auto& last = ray.back();
  const double y = last.y + (*fo.image_z - last.z) * last.u;
  if (y == 0.0) throw std::invalid_argument("sources: chief ray does not reach the image");
  return y;
}

/// Field point `field` of the system.
const model::Field& field_point(const Context& c, std::uint16_t field) {
  const auto& points = c.system->fields().points;
  if (field >= points.size()) {
    throw std::invalid_argument("sources: field index " + std::to_string(field) +
                                " does not exist");
  }
  return points[field];
}

/// Ray start for the field value `f`, interpreted with the system's field type.
FieldStart make_field(const Context& c, const model::Field& f) {
  const CompiledSystem& system = *c.system;
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
      // Finite object: the object point on the chief ray through the EP centre at the
      // reference wavelength (#31). With the EP at infinity every chief ray is parallel to the
      // axis, so a field angle does not define an object point (#96).
      if (!infinite && c.telecentric_ref) {
        throw std::invalid_argument(
            "sources: a field angle with a finite object needs a finite entrance pupil (object "
            "space telecentric: give object heights or paraxial image heights)");
      }
      if (!infinite && c.z_ep_ref == z_obj) {
        throw std::invalid_argument("sources: entrance pupil in the object (reference wavelength)");
      }
      hx = (z_obj - c.z_ep_ref) * tx;
      hy = (z_obj - c.z_ep_ref) * ty;
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
    if (!c.telecentric && c.z_ep == z_obj) {
      throw std::invalid_argument("sources: entrance pupil in the object");
    }
    start.object = Vec3(hx, hy, z_obj);
    return start;
  }
  start.direction = Vec3(tx, ty, 1.0).normalized();
  const Vec3& d = start.direction;
  const double d_xy = std::hypot(d.x(), d.y());
  // The bundle reaches at most 3 EP radii off the chief ray (reserve for real aiming); every
  // ray must start at least 1 mm before the EP and before every surface of the path it can
  // reach (lowest_z).
  const double r_bundle = kBundleReserve * c.r_ep;
  double z_min = c.z_ep;
  for (const auto& e : system.path(c.path).events) {
    z_min = std::min(z_min, lowest_z(system.surfaces()[e.surface], c.z_ep, r_bundle, d_xy / d.z()));
  }
  z_min -= 1.0;
  // Start plane perpendicular to d through E - L d (E = EP centre). A ray through (a, b, z_ep)
  // starts at z_ep - (a d_x + b d_y + L) d_z; with |(a, b)| <= r_bundle this is <= z_min for
  // L >= (z_ep - z_min) / d_z + r_bundle |d_xy| (decided for #8, same L for all rays of a field).
  start.plane_offset = (c.z_ep - z_min) / d.z() + r_bundle * d_xy;
  return start;
}

/// Ray through the point (a, b, z_ep) on the EP plane; object-space telecentric: the ray from
/// the object point with slopes (a, b).
RayState ray_through(const Context& c, const FieldStart& f, double a, double b) {
  const Vec3 q(a, b, c.z_ep);
  RayState ray;
  if (f.infinite) {
    const Vec3 e(0.0, 0.0, c.z_ep);
    const double s = f.direction.dot(q - e) + f.plane_offset;
    ray.pos = q - s * f.direction;
    ray.dir = f.direction;
  } else if (c.telecentric) {
    // EP at infinity (#96): (a, b) are the object-space slopes of the ray from P, into +z.
    ray.pos = f.object;
    ray.dir = Vec3(a, b, 1.0).normalized();
  } else {
    // The ray is the line through the object point P and Q and travels into the system, in +z
    // (light starts towards +z in object space, rtt/paraxial/paraxial.hpp). With a virtual EP
    // on the far side of the object (z_ep < z_obj, e.g. a stop behind the rear focal point of
    // the group before it) that is away from Q (#93). The line is unchanged, so Q still maps to
    // the conjugate stop point stop_scale * (a, b). P and all Q lie on the planes z_obj and
    // z_ep, so the sign is the same for every ray of every field (paths are rotationally
    // symmetric about z, make_context).
    ray.pos = f.object;
    ray.dir = c.z_ep < f.object.z() ? (f.object - q).normalized() : (q - f.object).normalized();
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
  // Pupil labels oriented at the reference wavelength (label_sign, #96); 1 * p == p exactly.
  const double sx = c.label_sign * px;
  const double sy = c.label_sign * py;
  // Paraxial solution: straight through the EP point (sx, sy) * r_ep, or with the slopes
  // (sx, sy) * u_m for an object-space telecentric system.
  double a = sx * c.r_ep;
  double b = sy * c.r_ep;
  AimedRay out;
  out.ray = ray_through(c, f, a, b);
  if (aiming == Aiming::Paraxial) return out;

  // Real aiming: Newton on the residual at the stop, central-difference Jacobian (ADR 0007;
  // Newton's method in 2D, e.g. Press et al., Numerical Recipes, 3rd ed., Sec. 9.6).
  out.residual = std::numeric_limits<double>::infinity();  // until the stop is reached
  const double r_s = c.r_ep * c.stop_scale;
  const double tx = sx * r_s;
  const double ty = sy * r_s;
  // Step in pupil units (ADR 0007, addendum #96): sized in stop units, but with a floor in
  // pupil units for a very distant entrance pupil, where |r_s / r_ep| = |stop_scale| -> 0.
  // Object-space telecentric (#96): the pupil coordinates are slopes, so the step is taken
  // relative to u_m directly (new path, no earlier results to keep).
  const double h = c.telecentric
                       ? kAimStepRelative * c.r_ep
                       : std::max(kAimStepRelative * std::abs(r_s), kAimStepPupilFloor * c.r_ep);
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
    out.ray = ray_through(c, f, sx * c.r_ep, sy * c.r_ep);
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
  if (!std::isfinite(px) || !std::isfinite(py)) {
    throw std::invalid_argument("sources: pupil coordinates must be finite");
  }
  const Context c = make_context(system, path, wavelength);
  return aim(c, make_field(c, field_point(c, field)), px, py, aiming);
}

AimedRay aim_ray(const compile::CompiledSystem& system,
                 compile::PathId path,
                 const model::Field& field,
                 std::uint16_t wavelength,
                 double px,
                 double py,
                 Aiming aiming) {
  if (!std::isfinite(px) || !std::isfinite(py)) {
    throw std::invalid_argument("sources: pupil coordinates must be finite");
  }
  if (!std::isfinite(field.x) || !std::isfinite(field.y)) {
    throw std::invalid_argument("sources: field values must be finite");
  }
  const Context c = make_context(system, path, wavelength);
  return aim(c, make_field(c, field), px, py, aiming);
}

namespace {

/// make_rays() with an optional control (run_control.hpp, #83). The rays are aimed in parallel
/// (#119): every ray from its own paraxial start value, written only to its own slot, with no
/// reduction over rays, so the result does not depend on the number of threads or on the
/// partitioning (ADR 0004, addendum #119). The preparation (context, pupil points, field starts)
/// stays serial and is the only part that throws. With an active control the rays are aimed in
/// blocks of control->block_size: progress (stage "aim") after each block, no new block after a
/// cancellation request or a failed callback, and the monitor throws only after the parallel
/// part (rule 3), exactly as in SequentialTracer::trace().
RayBatch make_rays_impl(const compile::CompiledSystem& system,
                        compile::PathId path,
                        std::span<const std::uint16_t> fields,
                        std::uint16_t wavelength,
                        const PupilSampling& sampling,
                        Aiming aiming,
                        const RunControl* control) {
  const Context c = make_context(system, path, wavelength);
  const std::vector<PupilPoint> points = pupil_points(sampling);
  for (const PupilPoint& p : points) {
    if (!std::isfinite(p.px) || !std::isfinite(p.py)) {
      throw std::invalid_argument("sources: pupil coordinates must be finite");
    }
  }
  std::vector<FieldStart> starts;
  starts.reserve(fields.size());
  for (const std::uint16_t f : fields) starts.push_back(make_field(c, field_point(c, f)));

  // Ray i belongs to field i / P and pupil point i % P (fields as the outer loop).
  RayBatch rays(fields.size() * points.size());
  const std::size_t per_field = points.size();
  const auto aim_range = [&](const oneapi::tbb::blocked_range<std::size_t>& range) {
    for (std::size_t i = range.begin(); i != range.end(); ++i) {
      const std::size_t f = i / per_field;
      const PupilPoint& p = points[i % per_field];
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
    }
  };
  if (control == nullptr || !control->active()) {
    oneapi::tbb::parallel_for(oneapi::tbb::blocked_range<std::size_t>(0, rays.size()), aim_range,
                              oneapi::tbb::static_partitioner());
    return rays;
  }
  RunMonitor monitor(*control, rays.size(), "aim");
  const std::size_t grain = std::max<std::size_t>(control->block_size, 1);
  oneapi::tbb::parallel_for(
      oneapi::tbb::blocked_range<std::size_t>(0, rays.size(), grain),
      [&](const oneapi::tbb::blocked_range<std::size_t>& range) {
        if (monitor.stop()) return;
        aim_range(range);
        monitor.add(range.size());
      },
      oneapi::tbb::simple_partitioner());
  monitor
      .finish();  // after the parallel part: Cancelled, the callback's exception, or done == total
  return rays;
}

}  // namespace

RayBatch make_rays(const compile::CompiledSystem& system,
                   compile::PathId path,
                   std::span<const std::uint16_t> fields,
                   std::uint16_t wavelength,
                   const PupilSampling& sampling,
                   Aiming aiming) {
  return make_rays_impl(system, path, fields, wavelength, sampling, aiming, nullptr);
}

RayBatch make_rays(const compile::CompiledSystem& system,
                   compile::PathId path,
                   std::span<const std::uint16_t> fields,
                   std::uint16_t wavelength,
                   const PupilSampling& sampling,
                   Aiming aiming,
                   const RunControl& control) {
  return make_rays_impl(system, path, fields, wavelength, sampling, aiming, &control);
}

}  // namespace rtt::trace
