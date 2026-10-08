#include "rtt/paraxial/prescription.hpp"

#include <cmath>
#include <cstddef>
#include <optional>
#include <vector>

#include "paraxial_rays.hpp"

namespace rtt::paraxial {
namespace {

/// Entrance pupil radius from the system aperture alone, for an object at infinity without a
/// stop: the same formulas as first_order() (EPD directly, EPD = EFL / F#); the other aperture
/// types need the stop or a finite object.
std::optional<double> radius_without_stop(const compile::CompiledSystem& system,
                                          const FirstOrder& fo) {
  const double value = system.aperture().value.value;
  switch (system.aperture().type) {
    case model::SystemApertureType::EntrancePupilDiameter:
      return 0.5 * value;
    case model::SystemApertureType::ImageSpaceFNumber:
      if (fo.efl) return 0.5 * (std::abs(*fo.efl) / value);
      return std::nullopt;
    case model::SystemApertureType::StopSize:
    case model::SystemApertureType::ObjectSpaceNA:
      return std::nullopt;
  }
  return std::nullopt;
}

}  // namespace

/// Sources and conventions in prescription.hpp: angle of incidence i = u + y c and Lagrange
/// invariant H = n (u_bar y - u y_bar) from Sasian, OPTI 517 L4, p. 24 (the same expressions
/// as in seidel(), so that both agree bit for bit); NA ~ n' |u'| and F/#_W = 1 / (2 n' |u'|)
/// from Greivenkamp, OPTI-502 Sec. 9, p. 9-34 to 9-37.
Prescription prescription(const compile::CompiledSystem& system,
                          compile::PathId path,
                          std::uint16_t wavelength) {
  Prescription out;
  // Validates path, wavelength and rotational symmetry (ParaxialError).
  out.first_order = first_order(system, path, wavelength);
  const FirstOrder& fo = out.first_order;
  const bool infinite = system.object().at_infinity;
  const double z_obj = infinite ? 0.0 : -system.object().distance.value;

  // Vertex z and signed index of every event from the axial ray; they do not depend on the
  // ray, so they exist without a marginal or chief ray.
  const std::vector<RayAtEvent> axis = trace_ray(system, path, wavelength, 0.0, 0.0, 0.0);

  std::optional<double> z_ep;
  std::optional<double> r_ep;
  if (fo.entrance_pupil) {
    z_ep = fo.entrance_pupil->z;
    if (fo.entrance_pupil->diameter) r_ep = 0.5 * *fo.entrance_pupil->diameter;
  } else if (infinite) {
    r_ep = radius_without_stop(system, fo);
  }
  const bool ep_usable = z_ep && (infinite || *z_ep != z_obj);

  if (r_ep && infinite) {
    // Parallel to the axis: the start plane only matters for the record (seidel() uses the
    // entrance pupil).
    const double z_start = z_ep ? *z_ep : (axis.empty() ? 0.0 : axis.front().z);
    out.marginal_start = detail::marginal_start(true, z_obj, z_start, *r_ep);
  } else if (r_ep && ep_usable) {
    out.marginal_start = detail::marginal_start(false, z_obj, *z_ep, *r_ep);
  }
  // The field conversion at the reference wavelength (#35) needs a usable pupil there as well;
  // without one the chief ray is none, as for an unusable pupil at `wavelength`.
  if (fo.entrance_pupil && ep_usable &&
      detail::reference_pupil_usable(system, path, wavelength, z_obj)) {
    out.chief_start =
        detail::chief_start(system, path, wavelength, fo, *z_ep, z_obj, "prescription");
  }

  std::vector<RayAtEvent> marginal;
  std::vector<RayAtEvent> chief;
  if (out.marginal_start) {
    const RayStart& m = *out.marginal_start;
    marginal = trace_ray(system, path, wavelength, m.z, m.y, m.u);
  }
  if (out.chief_start) {
    const RayStart& c = *out.chief_start;
    chief = trace_ray(system, path, wavelength, c.z, c.y, c.u);
  }

  const auto& events = system.path(path).events;
  double n = fo.object_index;  // signed index before the event; light starts towards +z
  double u = out.marginal_start ? out.marginal_start->u : 0.0;
  double u_bar = out.chief_start ? out.chief_start->u : 0.0;
  out.surfaces.reserve(events.size());
  for (std::size_t k = 0; k < events.size(); ++k) {
    const double c = detail::shape_data(system.surfaces()[events[k].surface]).c;
    PrescriptionSurface s;
    s.surface = events[k].surface;
    s.z = axis[k].z;
    s.n = axis[k].n;
    if (!marginal.empty())
      s.marginal = PrescriptionRay{marginal[k].y, marginal[k].u, u + marginal[k].y * c};
    if (!chief.empty()) s.chief = PrescriptionRay{chief[k].y, chief[k].u, u_bar + chief[k].y * c};
    if (!marginal.empty() && !chief.empty()) {
      s.lagrange = marginal[k].n * (chief[k].u * marginal[k].y - marginal[k].u * chief[k].y);
      if (k == 0) out.lagrange_invariant = n * (u_bar * marginal[k].y - u * chief[k].y);
    }
    if (k > 0) out.total_track += std::abs(axis[k].z - axis[k - 1].z);
    out.surfaces.push_back(s);

    n = axis[k].n;
    if (!marginal.empty()) u = marginal[k].u;
    if (!chief.empty()) u_bar = chief[k].u;
  }

  if (!infinite && !axis.empty()) out.object_distance = axis.front().z - z_obj;
  if (!marginal.empty()) {
    // The working F/# is defined for any cone of light (Greivenkamp p. 9-36), so also for an
    // afocal system with a finite object. An afocal system with the object at infinity leaves
    // the marginal ray parallel: u'_K = 0 in theory, decided by the afocal criterion of
    // first_order() (efl none) rather than by the rounding of u'_K.
    const bool collimated = !fo.efl && infinite;
    const double nu = collimated ? 0.0 : std::abs(marginal.back().n * marginal.back().u);
    out.paraxial_image_na = nu;
    if (nu != 0.0) out.paraxial_working_f_number = 1.0 / (2.0 * nu);
  }
  return out;
}

}  // namespace rtt::paraxial
