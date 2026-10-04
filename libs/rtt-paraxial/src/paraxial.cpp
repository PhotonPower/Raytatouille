#include "rtt/paraxial/paraxial.hpp"

#include <cmath>
#include <cstddef>
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

/// |Phi| at or below this value (1/mm) counts as afocal (EFL > 1e14 mm).
constexpr double kAfocalPower = 1e-14;

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
      throw ParaxialError("paraxial: surface " + surface_name(cs, event.surface) + ": " + what);
    };
    if (event.kind != EventKind::Refract && event.kind != EventKind::Reflect &&
        event.kind != EventKind::Transmit) {
      fail("diffraction and birefringent events have no paraxial model before M4");
    }
    if (!s.phases.empty()) fail("phase layers are not supported before M4");
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
/// lecture notes (2018), Sec. 9 "Paraxial Raytracing": p. 9-2 (YNU raytrace: refraction or
/// reflection n'u' = nu - y phi with phi = (n' - n) C, transfer y' = y + u' t'), pp. 9-26/9-27
/// (mirror system: n' = -n after each reflection, negative distances after a mirror).
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

  // Cardinal points, derived from the y-nu equations with the reduced matrix M = [[a, b],
  // [c, d]] from the plane of V1 (object space, index n1 > 0) to the plane of Vk (image space,
  // signed index n'):
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
  if (std::abs(fo.power) > kAfocalPower) {
    const double phi = fo.power;
    fo.efl = 1.0 / phi;
    fo.front_focal_length = n1 / phi;
    fo.rear_focal_length = fo.image_index / phi;
    fo.ffl = m.d * n1 / phi;
    fo.bfl = m.a * fo.image_index / phi;
    fo.front_focal_z = z_v1 - *fo.ffl;
    fo.rear_focal_z = z_vk + m.a * n_img / phi;
    fo.front_principal_z = *fo.front_focal_z + n1 / phi;
    fo.rear_principal_z = *fo.rear_focal_z - n_img / phi;
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
    if (nu_out != 0.0) {
      fo.image_z = z_last - y_out * n_img / nu_out;
      fo.lateral_magnification = n1 / nu_out;
    }
  }

  // Pupils: images of the stop through the steps before and after it.
  std::optional<std::size_t> stop;
  for (std::size_t i = 0; i < steps.size(); ++i) {
    if (steps[i].stop) {
      stop = i;
      break;
    }
  }
  if (!stop) return fo;

  const compile::CompiledSurface& stop_surface = system.surfaces()[steps[*stop].surface];
  const auto* circle = stop_surface.aperture
                           ? std::get_if<model::CircularAperture>(&*stop_surface.aperture)
                           : nullptr;
  if (circle == nullptr) {
    throw ParaxialError("paraxial: the stop " + surface_name(system, steps[*stop].surface) +
                        " needs a circular aperture for the pupils");
  }
  const double z_stop = steps[*stop].z;
  // Object space -> stop plane, and stop plane -> image space (the stop itself only transmits).
  const Matrix front = propagate(steps, 0, *stop, z_first, n1, z_stop);
  const Matrix back = propagate(steps, *stop + 1, steps.size(), z_stop, front.n_out, z_last);

  // Entrance pupil, derived from the y-nu equations with front = [[a, b], [c, d]] from z_first
  // (object space) to the stop plane:
  // - Position: an object-space ray (y, nu) = (y0, n1 u) at z_first passes the stop centre if
  //   a y0 + b n1 u = 0, i.e. y0 = -b n1 u / a; it meets the axis at
  //   z_first - y0/u = z_first + b n1 / a.
  // - Size: transferring the start plane to the pupil plane changes b but not a, and there
  //   b = 0 (conjugate planes), so y_stop = a y_pupil: the stop radius r appears as r / |a|.
  Pupil ep;
  Pupil xp;
  if (front.a != 0.0) ep.z = z_first + front.b * n1 / front.a;
  // Exit pupil with back = [[a, b], [c, d]] from the stop plane to z_last (image space, n'):
  // - Position: the ray from the stop centre (0, nu_s) leaves with (b nu_s, d nu_s) and meets
  //   the axis at z_last - b n' / d.
  // - Size: in the conjugate pupil plane the matrix is [[1/d, 0], [c, d]] (det = 1), so the
  //   stop radius r appears as r / |d|.
  if (back.d != 0.0) xp.z = z_last - back.b * n_img / back.d;

  const double value = system.aperture().value.value;
  switch (system.aperture().type) {
    case model::SystemApertureType::EntrancePupilDiameter:
      ep.diameter = value;
      break;
    case model::SystemApertureType::StopSize:
      if (front.a != 0.0) ep.diameter = 2.0 * circle->radius / std::abs(front.a);
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
  if (ep.diameter && back.d != 0.0) {
    xp.diameter = *ep.diameter * std::abs(front.a) / std::abs(back.d);
  }
  fo.entrance_pupil = ep;
  fo.exit_pupil = xp;

  // Chief ray (through the stop centre), from the y-nu equations: with y0 = -b n1 u / a (see
  // above) the stop is reached with nu_s = c y0 + d n1 u = n1 u (a d - b c) / a = n1 u / a, and
  // the system leaves with nu' = d_back nu_s. Slopes measured along the propagation direction
  // (u' sign(n') / u) give the ratio n1 d_back / (|n'| a).
  if (front.a != 0.0) {
    fo.angular_magnification = n1 * back.d / (fo.image_index * front.a);
  }
  return fo;
}

}  // namespace rtt::paraxial
