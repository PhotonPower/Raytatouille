#include "parabasal.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "rtt/math/types.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/ray_paths.hpp"
#include "rtt/trace/sequential.hpp"
#include "rtt/trace/sources.hpp"

namespace rtt::optim::detail {

namespace {

using Vec2 = std::array<double, 2>;
/// 2x2 matrix, row-major: {a, b, c, d} = [[a, b], [c, d]].
using Mat2 = std::array<double, 4>;

/// An orthonormal pair (e1, e2) perpendicular to the unit vector d, from the global x axis (the
/// global y axis if d is close to x). Only sums of singular values and dot products are taken in
/// it, so the choice does not enter the results.
struct Basis {
  math::Vec3 e1;
  math::Vec3 e2;
};

Basis basis_perpendicular(const math::Vec3& d) {
  const math::Vec3 a = std::abs(d.x()) < 0.9 ? math::Vec3::UnitX() : math::Vec3::UnitY();
  const math::Vec3 e1 = (a - a.dot(d) * d).normalized();
  // e2 = d x e1, written out (rtt-optim does not include Eigen's geometry module).
  const math::Vec3 e2(d.y() * e1.z() - d.z() * e1.y(), d.z() * e1.x() - d.x() * e1.z(),
                      d.x() * e1.y() - d.y() * e1.x());
  return {e1, e2};
}

/// Transverse slope of the direction `d` against the chief ray direction `c` in the basis: the
/// tangent of the angle in the planes (c, e1) and (c, e2), u = dy/dz of the source with z along c
/// (Greivenkamp, OPTI-201/202, Sec. 9, p. 9-2: u is the slope).
Vec2 slope(const math::Vec3& d, const math::Vec3& c, const Basis& b) {
  const double dc = d.dot(c);
  return {d.dot(b.e1) / dc, d.dot(b.e2) / dc};
}

Vec2 minus(const Vec2& a, const Vec2& b) {
  return {a[0] - b[0], a[1] - b[1]};
}

/// sigma_1 + sigma_2 of a 2x2 matrix A (addendum #256): sigma_1^2 + sigma_2^2 = ||A||_F^2 (the
/// Frobenius norm is the root of the sum of the squared singular values) and sigma_1 sigma_2 =
/// |det A| (the determinant is the product of the singular values up to its sign), so
/// (sigma_1 + sigma_2)^2 = ||A||_F^2 + 2 |det A|. Smooth also where sigma_1 = sigma_2, unlike the
/// single singular values, and independent of the bases on both sides (rotation about z,
/// reflections of an axis).
double singular_sum(const Mat2& a) {
  const double frobenius = a[0] * a[0] + a[1] * a[1] + a[2] * a[2] + a[3] * a[3];
  const double det = a[0] * a[3] - a[1] * a[2];
  return std::sqrt(frobenius + 2.0 * std::abs(det));
}

/// A = Y X^-1 for the columns X = (x_x, x_y) and Y = (y_x, y_y); none (NaN) if X is singular.
Mat2 map_of(const Vec2& xx, const Vec2& xy, const Vec2& yx, const Vec2& yy) {
  const double det = xx[0] * xy[1] - xy[0] * xx[1];
  const double nan = std::numeric_limits<double>::quiet_NaN();
  if (det == 0.0) return {nan, nan, nan, nan};
  // X^-1 = [[xy1, -xy0], [-xx1, xx0]] / det; A = Y X^-1 with Y = [[yx0, yy0], [yx1, yy1]].
  return {(yx[0] * xy[1] - yy[0] * xx[1]) / det, (yy[0] * xx[0] - yx[0] * xy[0]) / det,
          (yx[1] * xy[1] - yy[1] * xx[1]) / det, (yy[1] * xx[0] - yx[1] * xy[0]) / det};
}

/// The mean principal value (sigma_1 + sigma_2)/2 of A = Y X^-1 (map_of, singular_sum); none if
/// X is singular, i.e. the input differences of the two pairs do not span the plane (coincident
/// or collinear parabasal rays).
std::optional<double> mean_principal(const Vec2& xx,
                                     const Vec2& xy,
                                     const Vec2& yx,
                                     const Vec2& yy) {
  const Mat2 a = map_of(xx, xy, yx, yy);
  if (std::isnan(a[0])) return std::nullopt;
  return singular_sum(a) / 2.0;
}

constexpr const char* kDegenerate =
    "degenerate parabasal rays (their offsets do not span the plane)";

/// The five rays traced with recording: 0 chief, 1 (+delta, 0), 2 (-delta, 0), 3 (0, +delta),
/// 4 (0, -delta).
struct Bundle {
  trace::RayBatch rays;
  trace::RayPaths paths;
  std::string undefined;  ///< why the bundle is not usable; empty if all rays arrive
};

Bundle trace_bundle(const compile::CompiledSystem& system,
                    compile::PathId path,
                    std::uint16_t wavelength,
                    const trace::AimReference& reference,
                    double delta) {
  const std::array<std::array<double, 2>, 5> pupil = {
      {{0.0, 0.0}, {delta, 0.0}, {-delta, 0.0}, {0.0, delta}, {0.0, -delta}}};
  // The axis point of the object (ADR 0036, point 2): field value 0 in the field type of the
  // system, which the nominal aiming data convert as in the nominal system.
  const model::Field axis{0.0, 0.0, 1.0};
  Bundle b;
  b.rays = trace::RayBatch(pupil.size());
  for (std::size_t i = 0; i < pupil.size(); ++i) {
    const trace::AimedRay aimed = trace::aim_ray(system, path, axis, wavelength, pupil[i][0],
                                                 pupil[i][1], trace::Aiming::Real, reference);
    b.rays.pos_x()[i] = aimed.ray.pos.x();
    b.rays.pos_y()[i] = aimed.ray.pos.y();
    b.rays.pos_z()[i] = aimed.ray.pos.z();
    b.rays.dir_x()[i] = aimed.ray.dir.x();
    b.rays.dir_y()[i] = aimed.ray.dir.y();
    b.rays.dir_z()[i] = aimed.ray.dir.z();
    b.rays.wl()[i] = wavelength;
    b.rays.pupil_x()[i] = pupil[i][0];
    b.rays.pupil_y()[i] = pupil[i][1];
    b.rays.status()[i] = aimed.ray.status;
  }
  static_cast<void>(trace::SequentialTracer().trace(system, path, b.rays, b.paths));
  for (std::size_t r = 0; r < pupil.size(); ++r) {
    const std::size_t last = b.paths.slots - 1;
    if (b.paths.count[r] != b.paths.slots ||
        b.paths.status[r * b.paths.slots + last] != trace::RayStatus::Alive) {
      b.undefined = "a parabasal ray does not reach the image surface (lost, or not aimed)";
      break;
    }
  }
  return b;
}

/// Pair differences of a per-ray quantity: (q(+delta e_x) - q(-delta e_x), q(+delta e_y) -
/// q(-delta e_y)).
template <class F>
std::array<Vec2, 2> pairs(F q) {
  return {minus(q(1), q(2)), minus(q(3), q(4))};
}

ParabasalValue undefined(std::string why) {
  ParabasalValue v;
  v.undefined = std::move(why);
  return v;
}

double index_of(const compile::CompiledSystem& cs, std::uint32_t medium, std::uint16_t wl) {
  return std::abs(cs.media()[medium].index[wl].real());
}

}  // namespace

ParabasalValue parabasal_value(const model::FirstOrderOperand& op,
                               const compile::CompiledSystem& trial,
                               const compile::CompiledSystem* trial_infinite,
                               const WavelengthData& nominal,
                               const analysis::NominalReference& reference,
                               bool object_at_infinity,
                               compile::PathId path,
                               std::uint16_t wavelength,
                               double delta) {
  using Q = model::FirstOrderQuantity;
  const bool parallel = op.quantity == Q::Efl || op.quantity == Q::Bfl;
  if (parallel && !nominal.efl) return undefined("the path is afocal (no EFL)");
  if (op.quantity == Q::Magnification && !nominal.magnification) {
    return undefined("no paraxial magnification");
  }
  if (op.quantity == Q::ImageFNumber && !nominal.f_number) {
    // As rtt::paraxial::prescription in the nominal system (addendum #256).
    return undefined("no paraxial working F-number (no marginal ray, or collimated output)");
  }
  // efl and bfl: the axis-parallel bundle (ADR 0036, point 2), of the variant for a finite
  // object; image_fnumber and magnification: the bundle of the axis point of the object.
  const bool variant = parallel && !object_at_infinity;
  if (variant && trial_infinite == nullptr) {
    throw std::logic_error("optimize: the variant with the object at infinity is missing");
  }
  const compile::CompiledSystem& system = variant ? *trial_infinite : trial;
  const trace::AimReference* aim = nullptr;
  if (variant) {
    if (!nominal.infinite_aim.has_value()) {
      throw std::logic_error("optimize: no aiming data of the variant with the object at infinity");
    }
    aim = &nominal.infinite_aim.value();
  } else {
    aim = &reference.aim(wavelength);
  }
  const Bundle b = trace_bundle(system, path, wavelength, *aim, delta);
  if (!b.undefined.empty()) return undefined(b.undefined);

  const std::vector<compile::CompiledEvent>& events = system.path(path).events;
  const std::size_t last = b.paths.slots - 1;
  const double n_image = index_of(system, events.back().medium_after, wavelength);
  const math::Vec3 d_image = b.paths.direction_at(0, last);
  const Basis image = basis_perpendicular(d_image);
  const auto image_slope = [&](std::size_t r) {
    return slope(b.paths.direction_at(r, last), d_image, image);
  };

  switch (op.quantity) {
    case Q::Efl: {
      // phi = -n' u'_k / y_1 (Greivenkamp, OPTI-201/202, Sec. 9, p. 9-12): the map from the
      // transverse start offsets y_1 of the parallel rays to n' u' in image space, its mean
      // principal value (sigma_1 + sigma_2)/2 = |phi|; the sign from the nominal EFL.
      const math::Vec3 d_object = b.paths.direction_at(0, 0);
      const Basis object = basis_perpendicular(d_object);
      const auto offset = [&](std::size_t r) {
        const math::Vec3 o = b.paths.position_at(r, 0) - b.paths.position_at(0, 0);
        return Vec2{o.dot(object.e1), o.dot(object.e2)};
      };
      const auto nu = [&](std::size_t r) {
        const Vec2 u = image_slope(r);
        return Vec2{n_image * u[0], n_image * u[1]};
      };
      const std::array<Vec2, 2> x = pairs(offset);
      const std::array<Vec2, 2> y = pairs(nu);
      const std::optional<double> phi = mean_principal(x[0], x[1], y[0], y[1]);
      if (!phi.has_value()) return undefined(kDegenerate);
      if (!(phi.value() > 0.0)) return undefined("the trial is afocal (no EFL)");
      const double efl = 1.0 / phi.value();
      return {nominal.efl.has_value() && nominal.efl.value() < 0.0 ? -efl : efl, {}};
    }
    case Q::Bfl: {
      // BFD = -y_k / u'_k (Greivenkamp, p. 9-12), from the chief ray's point q_c on the last
      // vertex along the propagation: each parabasal ray is drawn from its point on that surface
      // to the plane through q_c perpendicular to the chief ray (offset h), with its slope w
      // against the chief ray; the focus at distance s has h + s w = 0, least squares over the
      // two pair differences (addendum #256).
      if (!nominal.last_vertex) return undefined("no event of the path changes the ray");
      const std::size_t slot = *nominal.last_vertex + 1;
      const math::Vec3 qc = b.paths.position_at(0, slot);
      const math::Vec3 dc = b.paths.direction_at(0, slot);
      const Basis t = basis_perpendicular(dc);
      const auto height = [&](std::size_t r) {
        const math::Vec3 q = b.paths.position_at(r, slot);
        const math::Vec3 d = b.paths.direction_at(r, slot);
        const math::Vec3 h = q - ((q - qc).dot(dc) / d.dot(dc)) * d - qc;
        return Vec2{h.dot(t.e1), h.dot(t.e2)};
      };
      const auto w = [&](std::size_t r) { return slope(b.paths.direction_at(r, slot), dc, t); };
      const std::array<Vec2, 2> h = pairs(height);
      const std::array<Vec2, 2> a = pairs(w);
      const double num =
          h[0][0] * a[0][0] + h[0][1] * a[0][1] + h[1][0] * a[1][0] + h[1][1] * a[1][1];
      const double den =
          a[0][0] * a[0][0] + a[0][1] * a[0][1] + a[1][0] * a[1][0] + a[1][1] * a[1][1];
      if (!(den > 0.0)) return undefined("the trial is afocal (no BFL)");
      return {-num / den, {}};
    }
    case Q::ImageFNumber: {
      // F/#_W = 1 / (2 n' |u'_m|) (Greivenkamp, OPTI-502, Sec. 9, pp. 9-35 to 9-37), u'_m the
      // paraxial marginal slope, the marginal ray at p = 1. The pupil coordinate of a ray is its
      // hit on the stop over the nominal R_s (the aim target is sigma p R_s, ADR 0007); the
      // measured hits are taken, not the requested p, because real aiming stops at a residual
      // of up to kAimTolerance = 1e-9 mm, which is 2e-7 of a parabasal offset delta R_s = 5e-3
      // mm. The mean principal value of the map p -> u' is |u'_m| (linear in p). The hits are
      // taken over the paraxial R_s (u = NA / n for object_na, as first_order), also where
      // object_na is read real and the rays aim at the real cone (ADR 0040, point 5).
      const std::optional<std::size_t> stop = compile::stop_event(system, path);
      if (!stop) return undefined("the path has no stop");
      const compile::CompiledSurface& s =
          system.surfaces()[system.path(path).events[*stop].surface];
      const double r_s = trace::paraxial_stop_radius(*aim);
      const auto pupil_of = [&](std::size_t r) {
        const math::Vec3 local = s.to_local.apply_point(b.paths.position_at(r, *stop + 1));
        return Vec2{local.x() / r_s, local.y() / r_s};
      };
      const std::array<Vec2, 2> x = pairs(pupil_of);
      const std::array<Vec2, 2> y = pairs(image_slope);
      const std::optional<double> u_image = mean_principal(x[0], x[1], y[0], y[1]);
      if (!u_image.has_value()) return undefined(kDegenerate);
      if (!(u_image.value() > 0.0)) return undefined("collimated output (no marginal slope)");
      return {1.0 / (2.0 * n_image * u_image.value()), {}};
    }
    case Q::Magnification: {
      // m = n u / (n' u') (docs/quellen.md, "Pupillen, Winkel- und Lateralvergrößerung"), u and
      // u' the marginal slopes in object and image space. Both from the same rays: the map from
      // the object-space slopes to the image-space slopes has the mean principal value u'/u.
      // The sign from the nominal magnification.
      const math::Vec3 d_object = b.paths.direction_at(0, 0);
      const Basis object = basis_perpendicular(d_object);
      const std::array<Vec2, 2> x =
          pairs([&](std::size_t r) { return slope(b.paths.direction_at(r, 0), d_object, object); });
      const std::array<Vec2, 2> y = pairs(image_slope);
      const std::optional<double> ratio = mean_principal(x[0], x[1], y[0], y[1]);  // u' / u
      if (!ratio.has_value()) return undefined(kDegenerate);
      if (!(ratio.value() > 0.0)) return undefined("collimated output (no marginal slope)");
      const double n_object = index_of(system, events.front().medium_before, wavelength);
      const double m = n_object / (n_image * ratio.value());
      return {nominal.magnification.has_value() && nominal.magnification.value() < 0.0 ? -m : m,
              {}};
    }
  }
  throw std::invalid_argument("optimize: invalid first-order quantity");
}

}  // namespace rtt::optim::detail
