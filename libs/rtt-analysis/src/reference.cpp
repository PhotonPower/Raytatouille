#include "reference.hpp"

#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "rtt/compile/errors.hpp"
#include "rtt/paraxial/paraxial.hpp"

namespace rtt::analysis::detail {

namespace {

math::Vec3 global_pos(const trace::RayBatch& rays, std::size_t i) {
  return {rays.pos_x()[i], rays.pos_y()[i], rays.pos_z()[i]};
}

math::Vec3 global_dir(const trace::RayBatch& rays, std::size_t i) {
  return {rays.dir_x()[i], rays.dir_y()[i], rays.dir_z()[i]};
}

}  // namespace

Reference reference_at(const math::Vec3& chief_pos,
                       const math::Vec3& chief_dir,
                       const std::optional<math::Vec3>& exit_pupil,
                       double focus_shift,
                       double lambda_mm) {
  Reference r;
  r.chief_pos = chief_pos;
  r.chief_dir = chief_dir;
  r.exit_pupil = exit_pupil;
  r.sphere.centre = chief_pos;
  // Through focus (ADR 0032, point 5): C' = C + Delta d_c on the chief ray, so W(0, 0) = 0 stays.
  // Only for Delta != 0, so that Delta = 0 keeps C bit for bit (C + 0 d could turn -0.0 into +0.0).
  if (focus_shift != 0.0) r.sphere.centre += focus_shift * chief_dir;
  if (!exit_pupil) {
    // Image-space telecentric (#102): the limit R -> infinity of the reference sphere.
    r.sphere.radius = std::numeric_limits<double>::infinity();
  } else {
    r.sphere.radius = (r.sphere.centre - *exit_pupil).norm();
    if (!(r.sphere.radius > 0.0)) {
      throw AnalysisError(focus_shift != 0.0
                              ? "analysis: the shifted reference point (focus_shift) lies on the "
                                "exit pupil"
                              : "analysis: the exit pupil lies on the image surface");
    }
    // The exit pupil lies at C - s d of the chief ray, s = +R upstream, s = -R downstream.
    r.branch = chief_dir.dot(r.sphere.centre - *exit_pupil) >= 0.0 ? 1.0 : -1.0;
  }
  r.lambda_mm = lambda_mm;
  return r;
}

Reference make_reference(const CompiledSystem& system,
                         PathId path,
                         std::uint16_t field,
                         trace::Aiming aiming,
                         double focus_shift,
                         const NominalReference* nominal) {
  compile::require_stop(system, path);  // before any ray is traced (ADR 0022)
  const std::uint16_t ref = system.reference_wavelength();
  const std::uint32_t image = image_surface(system, path);
  // A virtual reference ray (ADR 0032, addendum #229): a central obscuration does not stop it.
  const trace::RayBatch chief = trace_reference_rays(
      system, path, field, ref, trace::SinglePupilPoint{0.0, 0.0}, aiming, nominal);
  if (!arrived(chief, 0, image)) {
    throw_lost(system, chief, 0, field,
               "analysis: the chief ray of field " + std::to_string(field) +
                   " does not reach the image surface");
  }
  // The exit pupil of the nominal system (ADR 0036, point 2), computed there as here.
  std::optional<math::Vec3> exit_pupil;
  if (nominal != nullptr) {
    exit_pupil = nominal->exit_pupil();
  } else {
    const paraxial::FirstOrder fo = paraxial::first_order(system, path, ref);
    // With a stop on the path (require_stop above) first_order() always gives the exit pupil.
    if (!fo.exit_pupil) throw std::logic_error("analysis: no exit pupil despite a stop");
    const paraxial::Pupil& xp = *fo.exit_pupil;
    if (xp.z) exit_pupil = math::Vec3(0.0, 0.0, *xp.z);
  }
  return reference_at(global_pos(chief, 0), global_dir(chief, 0), exit_pupil, focus_shift,
                      system.wavelengths_um()[ref] * 1e-3);
}

/// OPL of a ray up to the reference sphere, up to a constant n' R common to all rays (it cancels
/// in W = OPL_ref - OPL_ray): OPL at the image surface minus |n'| (s - branch R), where s is the
/// signed path from the sphere to the image surface along the ray, on the crossing on the side
/// of the exit pupil. None if the line misses the sphere.
///
/// Derivation (definition of Wyant & Creath, Sec. I; #102): with p the hit on the image
/// surface, d the unit direction and C the sphere centre, |p - s d - C|^2 = R^2 gives
/// s^2 - 2 b s + |p - C|^2 - R^2 = 0, b = d . (p - C), so s = b +- root with
/// root = sqrt(q + R^2), q = b^2 - |p - C|^2. Then
///   s - branch R = b + branch (root - R) = b + branch q / (root + R),
/// which is free of the cancellation of root - R and also exact for R -> infinity:
/// s - branch R -> b. The OPD against the sphere of infinite radius about C is therefore
/// W_inf = OPL_chief - OPL_ray + n' d . (p - C) (the path to the foot of the perpendicular from C
/// on the ray), and W(R) = W_inf + O(1 / R), continuous from both sides of infinity.
std::optional<double> opl_to_sphere(
    const math::Vec3& p, const math::Vec3& d, double opl, const Reference& ref, double n_image) {
  const math::Vec3 pc = p - ref.sphere.centre;
  const double b = d.dot(pc);
  if (std::isinf(ref.sphere.radius)) return opl - n_image * b;
  const double radius = ref.sphere.radius;
  const double q = b * b - pc.squaredNorm();
  const double disc = q + radius * radius;
  if (!(disc >= 0.0)) return std::nullopt;
  const double root = std::sqrt(disc);
  return opl - n_image * (b + ref.branch * q / (root + radius));
}

double image_index(const CompiledSystem& system, PathId path, std::uint16_t wavelength) {
  // Isotropic image medium only. A crystal (ADR 0026) cannot be the image medium here:
  // make_reference and make_rays call paraxial::first_order, which rejects the ordinary and
  // extraordinary events that lead into a crystal (guard test in test_opd.cpp). An OPD inside a
  // crystal needs the index of the mode and is a follow-up after M4 (#35).
  const auto& last = system.path(path).events.back();
  return std::abs(system.media()[last.medium_after].index[wavelength].real());
}

trace::RayBatch chief_ray(const CompiledSystem& system,
                          PathId path,
                          std::uint16_t field,
                          std::uint16_t wavelength,
                          trace::Aiming aiming,
                          const NominalReference* nominal) {
  // Reference ray: chief ray of the same wavelength (W(0, 0) = 0, decided for #29).
  trace::RayBatch chief = trace_reference_rays(system, path, field, wavelength,
                                               trace::SinglePupilPoint{0.0, 0.0}, aiming, nominal);
  if (!arrived(chief, 0, image_surface(system, path))) {
    throw_lost(system, chief, 0, field,
               "analysis: the chief ray of field " + std::to_string(field) +
                   " does not reach the image surface at wavelength " + std::to_string(wavelength));
  }
  return chief;
}

double chief_opl(
    const math::Vec3& p, const math::Vec3& d, double opl, const Reference& ref, double n_image) {
  const std::optional<double> o = opl_to_sphere(p, d, opl, ref, n_image);
  if (!o) throw AnalysisError("analysis: the chief ray misses the reference sphere");
  return *o;
}

OpdPoint point_against(const math::Vec3& p,
                       const math::Vec3& d,
                       double opl,
                       bool arrived_on_image,
                       trace::RayStatus status,
                       double opl_chief,
                       const Reference& ref,
                       double n_image) {
  OpdPoint q;
  std::optional<double> to_sphere;
  if (arrived_on_image) to_sphere = opl_to_sphere(p, d, opl, ref, n_image);
  if (to_sphere) {
    // Wyant & Creath, Sec. I: W > 0 if the wavefront leads, W = OPL_ref - OPL_ray.
    q.w = (opl_chief - *to_sphere) / ref.lambda_mm;
  } else {
    q.status = status == trace::RayStatus::Alive ? trace::RayStatus::Vignetted : status;
  }
  return q;
}

std::vector<OpdPoint> points_against(const CompiledSystem& system,
                                     PathId path,
                                     std::uint16_t field,
                                     std::uint16_t wavelength,
                                     const trace::PupilSampling& sampling,
                                     trace::Aiming aiming,
                                     const Reference& ref,
                                     LossCounter& losses,
                                     const trace::RunControl& control,
                                     const NominalReference* nominal) {
  const std::uint32_t image = image_surface(system, path);
  const double n_image = image_index(system, path, wavelength);
  const trace::RayBatch chief = chief_ray(system, path, field, wavelength, aiming, nominal);
  const double opl_chief =
      chief_opl(global_pos(chief, 0), global_dir(chief, 0), chief.opl()[0], ref, n_image);

  const trace::RayBatch rays =
      trace_rays(system, path, field, wavelength, sampling, aiming, control, nominal);
  std::vector<OpdPoint> points;
  points.reserve(rays.size());
  for (std::size_t i = 0; i < rays.size(); ++i) {
    OpdPoint q = point_against(global_pos(rays, i), global_dir(rays, i), rays.opl()[i],
                               arrived(rays, i, image), rays.status()[i], opl_chief, ref, n_image);
    q.px = rays.pupil_x()[i];
    q.py = rays.pupil_y()[i];
    // Counted with the status of the point, so that losses agree with arrived and vignetted
    // (a ray on the image surface that misses the reference sphere is Vignetted there).
    losses.add_ray(q.status, rays.last_surface()[i], global_pos(rays, i));
    points.push_back(q);
  }
  return points;
}

}  // namespace rtt::analysis::detail
