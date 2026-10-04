#include "rtt/paraxial/seidel.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <string>
#include <variant>

namespace rtt::paraxial {
namespace {

using compile::CompiledSystem;
using compile::PathId;

/// Curvature, conic constant and A4 of a surface in global orientation.
struct ShapeData {
  double c = 0.0;
  double conic = 0.0;
  double a4 = 0.0;
};

/// The paraxial trace (paraxial.cpp) uses c in global orientation: 1/R times the sign of the
/// local z axis along global z. The fourth-order sag c^3 r^4 / 8 (K + 1) + A4 r^4 changes sign
/// with the axis as well, so A4 flips with c while K, a ratio, stays.
ShapeData shape_data(const compile::CompiledSurface& s) {
  ShapeData d;
  const auto [c_local, k] =
      std::visit([](const auto& shape) { return shape.base_conic(); }, s.shape);
  double a4_local = 0.0;
  if (const auto* asphere = std::get_if<geom::EvenAsphere<double>>(&s.shape)) {
    if (!asphere->coefficients().empty()) a4_local = asphere->coefficients().front();
  }
  const double sign = s.to_global.apply_vector(math::Vec3::UnitZ()).z() > 0.0 ? 1.0 : -1.0;
  d.c = sign * c_local;
  d.conic = k;
  d.a4 = sign * a4_local;
  return d;
}

std::string field_error(const std::string& what) {
  return "paraxial: seidel: " + what;
}

/// Chief ray of the maximum field (see seidel.hpp), in object space.
RayStart chief_ray(const CompiledSystem& system,
                   PathId path,
                   std::uint16_t wavelength,
                   const FirstOrder& fo,
                   double z_ep,
                   double z_obj) {
  const bool infinite = system.object().at_infinity;
  const auto& points = system.fields().points;
  // Largest radial field value: tan theta for angles (chief ray d ~ (tan theta_x, tan theta_y,
  // 1), docs/architecture.md "Feldwinkel und Pupille"), mm for heights.
  double value = 0.0;
  for (const auto& f : points) {
    double radial = std::hypot(f.x, f.y);
    if (system.fields().type == model::FieldType::AngleDeg) {
      if (!(std::abs(f.x) < 90.0 && std::abs(f.y) < 90.0)) {
        throw ParaxialError(field_error("field angles must lie in (-90, 90) degree"));
      }
      radial = std::hypot(std::tan(f.x * std::numbers::pi / 180.0),
                          std::tan(f.y * std::numbers::pi / 180.0));
    }
    value = std::max(value, radial);
  }

  // Ray through the EP centre with unit field value: unit slope (angle; also for a finite
  // object, whose object point then lies on this ray as in rtt-trace, #8) or unit object height.
  const auto unit_ray = [&](bool angle) -> RayStart {
    if (angle) return {z_ep, 0.0, 1.0};
    return {z_obj, 1.0, -1.0 / (z_ep - z_obj)};
  };
  switch (system.fields().type) {
    case model::FieldType::AngleDeg: {
      const RayStart r = unit_ray(true);
      return {r.z, 0.0, value};
    }
    case model::FieldType::ObjectHeight: {
      if (infinite) {
        throw ParaxialError(field_error("an object height needs a finite object distance"));
      }
      const RayStart r = unit_ray(false);
      return {r.z, value, value * r.u};
    }
    case model::FieldType::ParaxialImageHeight: {
      // Linear in the field value: scale the unit chief ray to the requested paraxial image
      // height (as rtt-trace does for #8).
      if (!fo.image_z) {
        throw ParaxialError(field_error("a paraxial image height needs a finite paraxial image"));
      }
      const RayStart r = unit_ray(infinite);
      const auto ray = trace_ray(system, path, wavelength, r.z, r.y, r.u);
      const RayAtEvent& last = ray.back();
      const double unit = last.y + (*fo.image_z - last.z) * last.u;
      if (unit == 0.0) throw ParaxialError(field_error("the chief ray does not reach the image"));
      const double scale = value / unit;
      return {r.z, scale * r.y, scale * r.u};
    }
  }
  return {};
}

}  // namespace

/// Seidel sums of one surface.
///
/// Source: J. Sasian, OPTI 517 Lens Design, lecture notes L4 "Seidel aberration coefficients",
/// https://wp.optics.arizona.edu/jsasian/wp-content/uploads/sites/33/2016/02/
/// L4-OPTI-517-Aberration-Coefficients.pdf, read in full text (slides without equation numbers):
/// - p. 24: A = n i = n u + n y c, A_bar = n u_bar + n y_bar c, Lagrange invariant
///   H (Zhe in the source) = n u_bar y - n u y_bar = A_bar y - A y_bar, c = 1/r,
///   P = c Delta(1/n).
/// - p. 23: S_I = -sum A^2 y Delta(u/n), S_II = -sum A A_bar y Delta(u/n),
///   S_III = -sum A_bar^2 y Delta(u/n), S_IV = -H^2 sum P, and S_V in the second form
///   S_V = -sum A_bar [A_bar^2 Delta(1/n^2) y - (H + A_bar y) y_bar P], which needs no division
///   by A (the first form (A_bar/A)[H^2 P + A_bar^2 y Delta(u/n)] fails for A = 0, an axial
///   object at the centre of curvature). Both are equal: with u/n = A/n^2 - y c/n,
///   Delta(u/n) = A Delta(1/n^2) - y P, and H^2 - A_bar^2 y^2 = -A y_bar (H + A_bar y).
///   C_L = sum A y Delta(dn/n), C_T = sum A_bar y Delta(dn/n); W040 = S_I/8, W131 = S_II/2,
///   W222 = S_III/2, W220 = (S_IV + S_III)/4, W311 = S_V/2, d W020 = C_L/2, d W111 = C_T.
/// - p. 7: u = -y/s, the slope of the ray (s < 0 for an object on the left), the same as our
///   u = dy/dz. Also arXiv:2203.02302 (EUV lithography objective), Eq. (3) and (4): the same
///   formulas for mirrors with n_(i+1) = n'_i = -n_i, our signed-index convention.
/// Chromatic terms also in OPTI 518 L6 "Chromatic aberrations", p. 4 (single surface,
/// d W020 = 1/2 A Delta(dn/n) y, Delta(dn/n) = dn'/n' - dn/n) and p. 14 (sum over surfaces).
/// Aspheric cap: OPTI 518 L14 "Aspheric surfaces and stop shifting", p. 16/17: a = 8 A4 y^4
/// Delta(n) for A4, a = -eps^2 c^3 y^4 Delta(n) for a conic (eps^2 = -K); dS_I = a,
/// dS_II = (y_bar/y) a, dS_III = (y_bar/y)^2 a, dS_IV = 0, dS_V = (y_bar/y)^3 a,
/// dC_L = dC_T = 0. Also arXiv:2202.11378 (TianQin telescope), Eq. (7)-(9):
/// dS_I = 8 (K c^3/8 + c4) y^4 Delta(n).
/// Conventions checked in #30 against a real ray trace; see seidel.hpp and docs/quellen.md.
SeidelTerms surface_seidel(const SeidelSurfaceInput& in) {
  const double n = in.n;
  const double n2 = in.n_after;
  const double y = in.y;
  const double yb = in.y_bar;
  // Slope after the surface, y-nu trace (Greivenkamp, OPTI-201/202, p. 9-2).
  const double u2 = (n * in.u - y * in.c * (n2 - n)) / n2;

  const double a = n * (in.u + y * in.c);
  const double ab = n * (in.u_bar + yb * in.c);
  const double h = n * (in.u_bar * y - in.u * yb);
  const double d_u_n = u2 / n2 - in.u / n;
  const double p = in.c * (1.0 / n2 - 1.0 / n);
  const double d_inv_n2 = 1.0 / (n2 * n2) - 1.0 / (n * n);

  SeidelTerms t;
  t.s1 = -a * a * y * d_u_n;
  t.s2 = -a * ab * y * d_u_n;
  t.s3 = -ab * ab * y * d_u_n;
  t.s4 = -h * h * p;
  t.s5 = -ab * (ab * ab * d_inv_n2 * y - (h + ab * y) * yb * p);

  // Aspheric cap, (y_bar/y)^k a written without division by y.
  const double c3 = in.c * in.c * in.c;
  const double cap = (in.conic * c3 + 8.0 * in.a4) * (n2 - n);
  const double y2 = y * y;
  t.s1 += cap * y2 * y2;
  t.s2 += cap * y2 * y * yb;
  t.s3 += cap * y2 * yb * yb;
  t.s5 += cap * y * yb * yb * yb;

  const double d_dn_n = in.dn_after / n2 - in.dn / n;
  t.c_l = a * y * d_dn_n;
  t.c_t = ab * y * d_dn_n;
  return t;
}

Seidel seidel(const compile::CompiledSystem& system,
              compile::PathId path,
              std::uint16_t wavelength,
              std::optional<ChromaticPair> chromatic) {
  // Validates path, wavelength and rotational symmetry (ParaxialError).
  const FirstOrder fo = first_order(system, path, wavelength);
  const std::size_t n_wl = system.wavelengths_um().size();
  if (chromatic && (chromatic->first >= n_wl || chromatic->second >= n_wl)) {
    throw ParaxialError(field_error("wavelength index of the chromatic pair does not exist"));
  }
  if (!fo.entrance_pupil) throw ParaxialError(field_error("the path needs a stop"));
  const Pupil& ep = *fo.entrance_pupil;
  if (!ep.z) throw ParaxialError(field_error("the entrance pupil lies at infinity"));
  if (!ep.diameter) {
    throw ParaxialError(field_error("the system aperture does not define the entrance pupil"));
  }
  const double z_ep = *ep.z;
  const double r_ep = 0.5 * *ep.diameter;
  const bool infinite = system.object().at_infinity;
  const double z_obj = infinite ? 0.0 : -system.object().distance.value;
  if (!infinite && z_ep == z_obj) {
    throw ParaxialError(field_error("the entrance pupil lies in the object plane"));
  }

  Seidel out;
  out.chromatic = chromatic;
  // Marginal ray: axial object point to the rim of the entrance pupil at +y.
  out.marginal = infinite ? RayStart{z_ep, r_ep, 0.0} : RayStart{z_obj, 0.0, r_ep / (z_ep - z_obj)};
  out.chief = chief_ray(system, path, wavelength, fo, z_ep, z_obj);

  const auto marginal =
      trace_ray(system, path, wavelength, out.marginal.z, out.marginal.y, out.marginal.u);
  const auto chief = trace_ray(system, path, wavelength, out.chief.z, out.chief.y, out.chief.u);

  const auto& events = system.path(path).events;
  const auto dn_of = [&](std::uint32_t medium, double signed_n) {
    if (!chromatic) return 0.0;
    const auto& index = system.media()[medium].index;
    const double dn = index[chromatic->first].real() - index[chromatic->second].real();
    return signed_n > 0.0 ? dn : -dn;
  };

  double n = fo.object_index;  // signed index before the event; light starts towards +z
  double u = out.marginal.u;
  double u_bar = out.chief.u;
  out.surfaces.reserve(events.size());
  for (std::size_t i = 0; i < events.size(); ++i) {
    const ShapeData shape = shape_data(system.surfaces()[events[i].surface]);
    SeidelSurfaceInput in;
    in.c = shape.c;
    in.conic = shape.conic;
    in.a4 = shape.a4;
    in.n = n;
    in.n_after = marginal[i].n;
    in.y = marginal[i].y;
    in.u = u;
    in.y_bar = chief[i].y;
    in.u_bar = u_bar;
    in.dn = dn_of(events[i].medium_before, n);
    in.dn_after = dn_of(events[i].medium_after, marginal[i].n);

    SeidelSurface s;
    s.surface = events[i].surface;
    s.y = in.y;
    s.y_bar = in.y_bar;
    s.a = n * (u + in.y * in.c);
    s.a_bar = n * (u_bar + in.y_bar * in.c);
    s.terms = surface_seidel(in);
    s.lagrange = marginal[i].n * (chief[i].u * marginal[i].y - marginal[i].u * chief[i].y);
    if (i == 0) out.lagrange = n * (u_bar * in.y - u * in.y_bar);  // object space

    out.sum.s1 += s.terms.s1;
    out.sum.s2 += s.terms.s2;
    out.sum.s3 += s.terms.s3;
    out.sum.s4 += s.terms.s4;
    out.sum.s5 += s.terms.s5;
    out.sum.c_l += s.terms.c_l;
    out.sum.c_t += s.terms.c_t;
    out.surfaces.push_back(s);

    n = marginal[i].n;
    u = marginal[i].u;
    u_bar = chief[i].u;
  }
  return out;
}

}  // namespace rtt::paraxial
