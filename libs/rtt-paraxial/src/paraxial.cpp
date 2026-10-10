#include "rtt/paraxial/paraxial.hpp"

#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <variant>

namespace rtt::paraxial {
namespace {

using compile::CompiledSystem;
using compile::PathId;
using model::EventKind;

/// Largest off-axis vertex position (mm) and largest transverse component of a surface's local
/// z axis that still count as rotationally symmetric (decided for #7).
constexpr double kSymmetryTolerance = 1e-12;

/// Factor k of the afocal threshold |Phi| <= k N u S (power_scale(), #35 B9).
constexpr double kAfocalFactor = 16.0;

/// One event of the path reduced to what the y-nu trace needs.
struct Step {
  double z = 0.0;        ///< global z of the vertex, mm
  double c = 0.0;        ///< vertex curvature in global orientation, 1/mm
  double n_after = 1.0;  ///< |n| of the medium after the event
  EventKind kind = EventKind::Transmit;
  bool stop = false;  ///< surface of a Stop element
  std::uint32_t surface = 0;
};

/// Paraxial ray in reduced form: height y and optical direction cosine nu = n * u (signed n).
struct State {
  double z = 0.0;
  double y = 0.0;
  double nu = 0.0;
  double n = 1.0;
};

/// Reduced 2x2 matrix (y, nu)_out = [[a, b], [c, d]] (y, nu)_in with det = 1, plus the signed
/// index at the output plane.
struct Matrix {
  double a = 1.0, b = 0.0, c = 0.0, d = 1.0;
  double n_out = 1.0;
};

std::string surface_name(const CompiledSystem& cs, std::uint32_t surface) {
  return "'" + cs.surfaces()[surface].id.str() + "'";
}

std::vector<Step> prepare(const CompiledSystem& cs, PathId path, std::uint16_t wavelength) {
  if (path.index >= cs.paths().size()) {
    throw ParaxialError("paraxial: path index " + std::to_string(path.index) + " does not exist");
  }
  if (wavelength >= cs.wavelengths_um().size()) {
    throw ParaxialError("paraxial: wavelength index " + std::to_string(wavelength) +
                        " does not exist");
  }
  const auto& events = cs.path(path).events;
  if (events.empty()) throw ParaxialError("paraxial: path has no events");

  std::vector<Step> steps;
  steps.reserve(events.size());
  for (const auto& event : events) {
    const compile::CompiledSurface& s = cs.surfaces()[event.surface];
    const auto fail = [&](const std::string& what) {
      throw ParaxialError("paraxial: surface " + surface_name(cs, event.surface) + ": " + what, s);
    };
    if (event.kind != EventKind::Refract && event.kind != EventKind::Reflect &&
        event.kind != EventKind::Transmit) {
      fail("ordinary and extraordinary events have no paraxial model");
    }
    // ADR 0025, point 4: order 0 at a surface with phase layers is the surface without them
    // (point 2); other orders have no paraxial model in M4.
    if (event.order != 0) fail("diffraction orders other than 0 have no paraxial model");
    // Interim state of #178 until rtt-paraxial knows the ideal lens (ADR 0031, point 8): an
    // error, never an EFL without the lens. The cylinder lens stays an error afterwards.
    if (s.ideal_lens) fail("ideal lenses have no paraxial model yet (#178)");
    const math::Vec3 vertex = s.to_global.translation();
    const math::Vec3 axis = s.to_global.apply_vector(math::Vec3::UnitZ());
    if (std::abs(vertex.x()) > kSymmetryTolerance || std::abs(vertex.y()) > kSymmetryTolerance) {
      fail("decentred from the z axis; the path is not rotationally symmetric");
    }
    if (std::abs(axis.x()) > kSymmetryTolerance || std::abs(axis.y()) > kSymmetryTolerance) {
      fail("tilted against the z axis; the path is not rotationally symmetric");
    }
    // The paraxial curvature is the vertex curvature of the base conic; its sign follows the
    // orientation of the local z axis in global coordinates.
    const double c_local =
        std::visit([](const auto& shape) { return shape.base_conic().first; }, s.shape);
    Step step;
    step.z = vertex.z();
    step.c = axis.z() > 0.0 ? c_local : -c_local;
    step.n_after = cs.media()[event.medium_after].index[wavelength].real();
    step.kind = event.kind;
    step.stop = s.element_kind == model::ElementKind::Stop;
    step.surface = event.surface;
    steps.push_back(step);
  }
  return steps;
}

/// |n| in object space: the medium before the first event.
double object_index(const CompiledSystem& cs, PathId path, std::uint16_t wavelength) {
  return cs.media()[cs.path(path).events.front().medium_before].index[wavelength].real();
}

/// Transfer to the next vertex, then refraction or reflection there (y-nu trace).
///
/// Transfer: y' = y + (z' - z) u with u = nu / n (global slope, global z).
/// Refraction: n'u' = nu - y phi, phi = c (n' - n). A reflection is a refraction with
/// n' = -n. Source: J. E. Greivenkamp, OPTI-201/202 Geometrical and Instrumental Optics,
/// lecture notes (2018), Sec. 9 "Paraxial Raytracing", read in full text: p. 9-2 (YNU raytrace:
/// omega = n u, phi = (n' - n) C, refraction or reflection n'u' = nu - y phi, transfer
/// y' = y + u't'), pp. 9-26/9-27 (Cassegrain example: n = 1, -1, 1, t = -80 mm after the
/// primary, C1 = -0.005 for R1 = -200 mm, u = dy/dz in the global frame). The notes have no
/// equation numbers; see docs/quellen.md.
/// With the signed index n * u is |n| times the direction cosine along the propagation
/// direction, so both equations hold unchanged after mirrors with global z and c.
void apply(State& st, const Step& step) {
  st.y += (step.z - st.z) * st.nu / st.n;
  st.z = step.z;
  const double direction = st.n > 0.0 ? 1.0 : -1.0;
  const double n_after = step.kind == EventKind::Reflect ? -st.n : direction * step.n_after;
  const double phi = step.c * (n_after - st.n);
  st.nu -= st.y * phi;
  st.n = n_after;
}

/// Reduced matrix from the plane z_in (index n_in, signed) through steps [first, last) to the
/// plane z_out. Traced with the two basis rays (y, nu) = (1, 0) and (0, 1).
///
/// Both y-nu equations are linear in (y, nu): transfer [[1, dz/n], [0, 1]], refraction
/// [[1, 0], [-phi, 1]]. Each has determinant 1, so the product has a d - b c = 1.
Matrix propagate(const std::vector<Step>& steps,
                 std::size_t first,
                 std::size_t last,
                 double z_in,
                 double n_in,
                 double z_out) {
  State r1{z_in, 1.0, 0.0, n_in};
  State r2{z_in, 0.0, 1.0, n_in};
  for (std::size_t i = first; i < last; ++i) {
    apply(r1, steps[i]);
    apply(r2, steps[i]);
  }
  r1.y += (z_out - r1.z) * r1.nu / r1.n;
  r2.y += (z_out - r2.z) * r2.nu / r2.n;
  return {r1.y, r2.y, r1.nu, r2.nu, r1.n};
}

/// Magnitudes (y, nu) of a y-nu trace for the infinity tests of first_order (#35, rest of B9):
/// the trace of propagate() for the start ray (y0, nu0) from the plane z_in (index n_in, signed)
/// through steps [first, last) to the plane z_out, with every term by its magnitude as in
/// power_scale(): |y| += (|z'| + |z|) |nu| / |n|, |nu| += |y| |phi|. The rounding model of
/// power_scale() applies unchanged: each step rounds y and nu at most 8 times relative to these
/// magnitudes, the final transfer y 4 more times, so an element of the matrix that is zero in
/// exact arithmetic comes out with at most 8 (N + 1) u times its magnitude for N steps.
struct Magnitudes {
  double y = 0.0;
  double nu = 0.0;
};

Magnitudes magnitude_trace(const std::vector<Step>& steps,
                           std::size_t first,
                           std::size_t last,
                           double z_in,
                           double n_in,
                           double z_out,
                           double y0,
                           double nu0) {
  double y = std::abs(y0);
  double nu = std::abs(nu0);
  double z = z_in;
  double n = n_in;
  for (std::size_t i = first; i < last; ++i) {
    const Step& step = steps[i];
    y += (std::abs(step.z) + std::abs(z)) * nu / std::abs(n);
    z = step.z;
    const double direction = n > 0.0 ? 1.0 : -1.0;
    const double n_after = step.kind == EventKind::Reflect ? -n : direction * step.n_after;
    nu += y * std::abs(step.c * (n_after - n));
    n = n_after;
  }
  y += (std::abs(z_out) + std::abs(z)) * nu / std::abs(n);
  return {y, nu};
}

/// Scale S of the power for the afocal test (#35, item B9): the y-nu trace of the ray
/// (y, nu) = (1, 0) of propagate() with every term by its magnitude, |y| += (|z'| + |z|) |nu| / |n|
/// and |nu| += |y| |phi|, so S is the sum of the magnitudes of the terms that make up Phi.
///
/// Derivation (no literature source; docs/quellen.md): in the rounding model of IEEE double,
/// fl(a op b) = (a op b)(1 + delta) with |delta| <= u = 2^-53, apply() rounds four times for y
/// (difference z' - z, division by n, product, sum) and four times for nu (difference n' - n,
/// two products, difference); the inputs z and c are themselves known only to relative u
/// (z comes from composed poses, so its error is about u |z|, not u |z' - z|, hence |z'| + |z|).
/// The index n counts as exact: its value from the dispersion formula defines the system. Then
/// n' - n is exact for 1/2 <= n'/n <= 2 (Sterbenz: both operands lie within a factor 2, so the
/// difference is a multiple of the last place of the smaller one and no larger than it, hence
/// representable), and -2 n at a mirror is exact. Counting n with an error u |n| instead would
/// put u |n| into |n' - n|, which S does not hold for nearly index-matched cemented surfaces.
/// To first order every term of the computed Phi then carries a relative error of at most 8 u
/// per step, so |fl(Phi) - Phi| <= 8 N u S for N steps. The threshold k = 16 is twice that
/// bound. Measured with the same operations for 1999 nearly afocal thick lenses
/// (R = +-r, n = 1.5, d = 6 r, r from 0.7 um to 1.4 mm, at z = 0, 10 and 1000 mm): at most
/// |Phi| = 0.30 N u S, so the threshold lies more than 50 times above the rounding.
/// Computed as the nu of magnitude_trace() for (1, 0) through all steps: the same operations in
/// the same order as before magnitude_trace() existed, so the afocal threshold is bitwise
/// unchanged; the final transfer (z_out) only touches y.
double power_scale(const std::vector<Step>& steps, double z_in, double n_in) {
  return magnitude_trace(steps, 0, steps.size(), z_in, n_in, z_in, 1.0, 0.0).nu;
}

/// True if `value`, an element of a matrix traced through `steps` events, is zero up to the
/// rounding of the trace: |value| <= 16 (N + 1) u S with the magnitude S of magnitude_trace(),
/// twice the bound of its rounding model, as the afocal threshold (#35, rest of B9).
bool zero_up_to_rounding(double value, std::size_t steps, double scale) {
  return std::abs(value) <= kAfocalFactor * static_cast<double>(steps + 1) *
                                std::numeric_limits<double>::epsilon() / 2.0 * scale;
}

bool changes_ray(const Step& step, double n_before) {
  return step.kind == EventKind::Reflect || step.n_after != n_before;
}

}  // namespace

std::vector<RayAtEvent> trace_ray(const compile::CompiledSystem& system,
                                  compile::PathId path,
                                  std::uint16_t wavelength,
                                  double z0,
                                  double y0,
                                  double u0) {
  const std::vector<Step> steps = prepare(system, path, wavelength);
  const double n0 = object_index(system, path, wavelength);
  State st{z0, y0, n0 * u0, n0};
  std::vector<RayAtEvent> out;
  out.reserve(steps.size());
  for (const Step& step : steps) {
    apply(st, step);
    out.push_back({st.z, st.y, st.nu / st.n, st.n});
  }
  return out;
}

FirstOrder first_order(const compile::CompiledSystem& system,
                       compile::PathId path,
                       std::uint16_t wavelength) {
  const std::vector<Step> steps = prepare(system, path, wavelength);
  const double n1 = object_index(system, path, wavelength);
  const double z_first = steps.front().z;
  const double z_last = steps.back().z;

  FirstOrder fo;
  fo.object_index = n1;

  // Vertices: first and last events that change the ray.
  std::optional<std::size_t> v1;
  std::optional<std::size_t> vk;
  double n_abs = n1;
  for (std::size_t i = 0; i < steps.size(); ++i) {
    if (changes_ray(steps[i], n_abs)) {
      if (!v1) v1 = i;
      vk = i;
    }
    n_abs = steps[i].n_after;
  }

  // Whole path between the vertices (or between first and last event if nothing changes the
  // ray). Transfers outside [V1, Vk] cancel because object and image space are homogeneous.
  const double z_v1 = v1 ? steps[*v1].z : z_first;
  const double z_vk = vk ? steps[*vk].z : z_last;
  const Matrix m = propagate(steps, 0, steps.size(), z_v1, n1, z_vk);
  const double n_img = m.n_out;  // signed
  fo.image_index = std::abs(n_img);
  fo.image_direction = n_img > 0.0 ? 1 : -1;

  // Cardinal points. Source: Greivenkamp, OPTI-201/202 lecture notes, Sec. 9, p. 9-12 (rear:
  // phi = -n'u'_k / y_1, f_E = 1/phi, f'_R = n'/phi, BFD = -y_k/u'_k, d' = BFD - f'_R) and
  // p. 9-14 (front: f_F = -n/phi, FFD = -y_1/u_1, d = FFD - f_F); BFD, FFD, d, d' are directed
  // distances along z from the vertex. Mapping to our outputs (conventions agreed in #7):
  // ffl = -FFD and front_focal_length = -f_F (Hecht: positive for converging);
  // rear_focal_length = |n'|/phi and bfl = BFD * sign(n') (measured along the propagation
  // direction, differs from f'_R and BFD only after an odd number of reflections);
  // H = V1 + d, H' = Vk + d' exactly as in the source.
  // In terms of the reduced matrix M = [[a, b], [c, d]] from the plane of V1 (object space,
  // index n1 > 0) to the plane of Vk (image space, signed index n'):
  // - Ray from infinity (y, nu) = (1, 0) leaves with (y, nu) = (a, c), so the power (definition
  //   agreed in #7) is Phi = -n'u'/y_1 = -c.
  // - Rear focal point F': that ray meets the axis at z_Vk - y/u' = z_Vk - a n'/c
  //   = z_Vk + a n'/Phi. Along the image-space direction (sign of n') this is
  //   bfl = a n'/Phi * sign(n') = a |n'|/Phi.
  // - Front focal point F: an axial ray from z_F with slope u reaches V1 with
  //   (y, nu) = ((z_V1 - z_F) u, n1 u) and must leave parallel: c (z_V1 - z_F) u + d n1 u = 0,
  //   so z_V1 - z_F = d n1 / Phi = ffl (light travels +z in object space).
  // - Principal points: the focal lengths f = n1/Phi and f' = |n'|/Phi are the distances H -> F
  //   and H' -> F' along the propagation direction: z_H = z_F + n1/Phi, z_H' = z_F' - n'/Phi.
  fo.power = -m.c;
  // Afocal if Phi is zero up to the rounding of the trace and of the coordinates (#35, B9):
  // |Phi| <= 16 N u S, see power_scale().
  const double afocal_threshold = kAfocalFactor * static_cast<double>(steps.size()) *
                                  std::numeric_limits<double>::epsilon() / 2.0 *
                                  power_scale(steps, z_v1, n1);
  if (std::abs(fo.power) > afocal_threshold) {
    const double phi = fo.power;
    fo.efl = 1.0 / phi;
    fo.front_focal_length = n1 / phi;
    fo.rear_focal_length = fo.image_index / phi;
    const double ffl = m.d * n1 / phi;
    const double front_focal_z = z_v1 - ffl;
    const double rear_focal_z = z_vk + m.a * n_img / phi;
    fo.ffl = ffl;
    fo.bfl = m.a * fo.image_index / phi;
    fo.front_focal_z = front_focal_z;
    fo.rear_focal_z = rear_focal_z;
    fo.front_principal_z = front_focal_z + n1 / phi;
    fo.rear_principal_z = rear_focal_z - n_img / phi;
  }

  // Paraxial image of the axial object point: the axial ray (y, nu) = (0, n1) from the object
  // leaves the last vertex with (y, nu) = (b n1, d n1) and meets the axis at z - y n'/nu.
  // Lateral magnification m = n u / (n' u') (definition agreed in #7) = n1 / (d n1) per unit
  // object slope.
  if (system.object().at_infinity) {
    fo.image_z = fo.rear_focal_z;
  } else {
    const double z_obj = -system.object().distance.value;
    const Matrix o = propagate(steps, 0, steps.size(), z_obj, n1, z_last);
    // Axial ray (y, nu) = (0, n1 * 1) at the object.
    const double y_out = o.b * n1;
    const double nu_out = o.d * n1;
    // Image at infinity if d is zero up to rounding (#35, rest of B9): the axial ray (0, 1).
    const double d_scale = magnitude_trace(steps, 0, steps.size(), z_obj, n1, z_last, 0.0, 1.0).nu;
    if (!zero_up_to_rounding(o.d, steps.size(), d_scale)) {
      fo.image_z = z_last - y_out * n_img / nu_out;
      fo.lateral_magnification = n1 / nu_out;
    }
  }

  // Pupils: images of the stop through the steps before and after it.
  std::size_t stop = steps.size();  // index of the stop event; steps.size() = no stop
  for (std::size_t i = 0; i < steps.size(); ++i) {
    if (steps[i].stop) {
      stop = i;
      break;
    }
  }
  if (stop == steps.size()) return fo;

  const compile::CompiledSurface& stop_surface = system.surfaces()[steps[stop].surface];
  const model::CircularAperture* circle = nullptr;
  if (stop_surface.aperture.has_value()) {
    circle = std::get_if<model::CircularAperture>(&stop_surface.aperture.value());
  }
  if (circle == nullptr) {
    throw ParaxialError("paraxial: the stop " + surface_name(system, steps[stop].surface) +
                            " needs a circular aperture for the pupils",
                        stop_surface);
  }
  const double z_stop = steps[stop].z;
  // Object space -> stop plane, and stop plane -> image space (the stop itself only transmits).
  const Matrix front = propagate(steps, 0, stop, z_first, n1, z_stop);
  const Matrix back = propagate(steps, stop + 1, steps.size(), z_stop, front.n_out, z_last);
  // Pupils at infinity if front.a or back.d is zero up to rounding (#35, rest of B9): a with the
  // magnitudes of the ray (1, 0) from z_first to the stop, d with those of (0, 1) from the stop.
  const bool ep_finite = !zero_up_to_rounding(
      front.a, stop, magnitude_trace(steps, 0, stop, z_first, n1, z_stop, 1.0, 0.0).y);
  const std::size_t after_stop = steps.size() - stop - 1;
  const bool xp_finite = !zero_up_to_rounding(
      back.d, after_stop,
      magnitude_trace(steps, stop + 1, steps.size(), z_stop, front.n_out, z_last, 0.0, 1.0).nu);

  // Entrance pupil, derived from the y-nu equations with front = [[a, b], [c, d]] from z_first
  // (object space) to the stop plane:
  // - Position: an object-space ray (y, nu) = (y0, n1 u) at z_first passes the stop centre if
  //   a y0 + b n1 u = 0, i.e. y0 = -b n1 u / a; it meets the axis at
  //   z_first - y0/u = z_first + b n1 / a.
  // - Size: transferring the start plane to the pupil plane changes b but not a, and there
  //   b = 0 (conjugate planes), so y_stop = a y_pupil: the stop radius r appears as r / |a|.
  Pupil ep;
  Pupil xp;
  if (ep_finite) ep.z = z_first + front.b * n1 / front.a;
  // Exit pupil with back = [[a, b], [c, d]] from the stop plane to z_last (image space, n'):
  // - Position: the ray from the stop centre (0, nu_s) leaves with (b nu_s, d nu_s) and meets
  //   the axis at z_last - b n' / d.
  // - Size: in the conjugate pupil plane the matrix is [[1/d, 0], [c, d]] (det = 1), so the
  //   stop radius r appears as r / |d|.
  if (xp_finite) xp.z = z_last - back.b * n_img / back.d;

  const double value = system.aperture().value.value;
  switch (system.aperture().type) {
    case model::SystemApertureType::EntrancePupilDiameter:
      ep.diameter = value;
      break;
    case model::SystemApertureType::StopSize:
      if (ep_finite) ep.diameter = 2.0 * circle->radius / std::abs(front.a);
      break;
    case model::SystemApertureType::ImageSpaceFNumber:
      // F-number at infinite conjugates, EPD = EFL / F#, also for finite objects (#7).
      if (fo.efl) ep.diameter = std::abs(*fo.efl) / value;
      break;
    case model::SystemApertureType::ObjectSpaceNA:
      // Paraxial marginal slope u = NA / n from the axial object point to the pupil edge.
      if (!system.object().at_infinity && ep.z) {
        const double z_obj = -system.object().distance.value;
        ep.diameter = 2.0 * value / n1 * std::abs(*ep.z - z_obj);
      }
      break;
  }
  // Stop radius actually used = EP radius * |a|; the exit pupil shows it magnified by 1 / |d|.
  if (ep.diameter && xp_finite) {
    xp.diameter = *ep.diameter * std::abs(front.a) / std::abs(back.d);
  }
  fo.entrance_pupil = ep;
  fo.exit_pupil = xp;

  // Chief ray (through the stop centre), from the y-nu equations: with y0 = -b n1 u / a (see
  // above) the stop is reached with nu_s = c y0 + d n1 u = n1 u (a d - b c) / a = n1 u / a, and
  // the system leaves with nu' = d_back nu_s. Slopes measured along the propagation direction
  // (u' sign(n') / u) give the ratio n1 d_back / (|n'| a).
  if (ep_finite) {
    fo.angular_magnification = n1 * back.d / (fo.image_index * front.a);
  }
  return fo;
}

}  // namespace rtt::paraxial
