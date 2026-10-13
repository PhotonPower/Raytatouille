#include "rtt/analysis/jones.hpp"

#include <Eigen/Geometry>
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numbers>
#include <numeric>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#include "common.hpp"
#include "reference.hpp"
#include "rtt/diagnostics/codes.hpp"
#include "rtt/math/types.hpp"
#include "rtt/polar/interface.hpp"
#include "rtt/polar/prt_analysis.hpp"
#include "rtt/trace/apply_event.hpp"
#include "rtt/trace/ray_paths.hpp"
#include "rtt/trace/run_control.hpp"
#include "rtt/trace/sequential.hpp"

namespace rtt::analysis {
namespace {

using compile::CompiledSystem;
using compile::PathId;
using math::Complex;
using math::Vec3;

/// A basis is not defined within this angle of a pole, rad (ADR 0034, point 5).
constexpr double kPoleAngle = 1e-6;
/// |x_B x k_c| at or below this: the chief ray grazes the image surface (ADR 0034, point 5, as
/// the frame of the PSF plane in ADR 0032, point 4).
constexpr double kGrazingLimit = 1e-6;
/// |det J(0, 0)| at or below this times ||J(0, 0)||_F^2: the global phase from the largest
/// element instead of the determinant (ADR 0034, point 5).
constexpr double kSingularDeterminant = 1e-12;
/// |c_s - c_p| above this times c_s: the power factors of a transmission differ, so J_field is
/// no scalar multiple of J and becomes NaN (ADR 0034, point 5; addendum #208 B). The threshold
/// comes from the accuracy (review of #242): below it J is divided by sqrt(c_s), which leaves the
/// p part off by |c_s - c_p| / (2 c_s) <= 5e-10 relative, below every other tolerance of J. For
/// real indices c_s = c_p = Re(q_t) / Re(q_i) analytically (Byrnes, Eqs. (21), (22), with
/// Re(n^2 q*) / |n|^2 = Re(q) for real n), up to rounding. For a weakly absorbing medium
/// n + i kappa they differ in second order: |c_s - c_p| / c_s ~ 2 kappa^2 xi^2 / (n^2 Re(q)^2)
/// (derived in the review of #242, checked against transmission_power_factors to 2 %), e.g. at
/// n = 1.5 and 80 deg 6.7e-13 for kappa = 1e-6 and 6.7e-11 for kappa = 1e-5, the order of
/// catalogue glasses at their UV edge: they keep a finite J.
constexpr double kEqualPowerFactors = 1e-9;

std::string path_pointer(PathId path) {
  return "/paths/" + std::to_string(path.index);
}

/// ADR 0034, H6: paths with crystal modes are rejected before any ray is traced (as the pupil
/// function, pupil.cpp).
void reject_crystal_modes(const CompiledSystem& system, PathId path) {
  for (const auto& event : system.path(path).events) {
    if (event.kind == model::EventKind::Ordinary || event.kind == model::EventKind::Extraordinary) {
      throw AnalysisError("analysis: paths with crystal modes are not supported (#35)",
                          path_pointer(path));
    }
  }
}

Vec3 position(const trace::RayBatch& r, std::size_t i) {
  return {r.pos_x()[i], r.pos_y()[i], r.pos_z()[i]};
}

Vec3 direction(const trace::RayBatch& r, std::size_t i) {
  return {r.dir_x()[i], r.dir_y()[i], r.dir_z()[i]};
}

/// Start directions k0 of a freshly made batch, before it is traced.
std::vector<Vec3> start_directions(const trace::RayBatch& rays) {
  std::vector<Vec3> k0(rays.size());
  for (std::size_t i = 0; i < rays.size(); ++i) k0[i] = direction(rays, i);
  return k0;
}

/// b^T P a for real unit vectors a, b: one element of B_out^T P B_in.
Complex project(const math::CMat3& p, const Vec3& b, const Vec3& a) {
  const Eigen::Vector3cd pa = p * a.cast<Complex>();
  return b.cast<Complex>().dot(pa);  // dot conjugates its first argument, b is real
}

/// Components of a transverse complex vector in the real orthonormal basis (e_1, e_2).
PolarizationEllipse ellipse_in(const Eigen::Vector3cd& v, const std::array<Vec3, 2>& basis) {
  return polarization_ellipse(basis[0].cast<Complex>().dot(v), basis[1].cast<Complex>().dot(v));
}

/// J = B_out^T P B_in, without the phase factor, row-major.
std::array<Complex, 4> jones_matrix(const math::CMat3& p,
                                    const std::array<Vec3, 2>& in,
                                    const std::array<Vec3, 2>& out) {
  return {project(p, out[0], in[0]), project(p, out[0], in[1]), project(p, out[1], in[0]),
          project(p, out[1], in[1])};
}

/// Global phase psi of J(0, 0) (ADR 0034, point 5).
double global_phase_of(const std::array<Complex, 4>& j) {
  const Complex det = j[0] * j[3] - j[1] * j[2];
  double frobenius2 = 0.0;
  for (const Complex& z : j) frobenius2 += std::norm(z);
  if (std::abs(det) > kSingularDeterminant * frobenius2) return std::arg(det) / 2.0;
  const auto largest = std::max_element(j.begin(), j.end(), [](const Complex& a, const Complex& b) {
    return std::abs(a) < std::abs(b);
  });
  return std::arg(*largest);
}

/// Throws if the chief ray grazes the image surface: |x_B x k_c| <= 1e-6 (ADR 0034, point 5).
void check_grazing(const CompiledSystem& system,
                   std::uint32_t image,
                   const Vec3& x_b,
                   const Vec3& k_c) {
  if (x_b.cross(k_c).norm() <= kGrazingLimit) {
    throw AnalysisError(
        "analysis: the chief ray grazes the image surface (its x axis is parallel to the chief "
        "ray); the output basis of the Jones pupil has no frame",
        system.surfaces()[image].location);
  }
}

/// Both chief frames (ADR 0034, point 5): the input frame from k_in, the output frame of the PSF
/// plane from x_B and k_out (ADR 0032, point 4).
void set_frames(JonesPupil& jp, const Vec3& k_in, const Vec3& k_out, const Vec3& x_b) {
  jp.k_in = k_in;
  const std::array<Vec3, 2> frame_in = input_chief_frame(k_in);
  jp.x_in = frame_in[0];
  jp.y_in = frame_in[1];
  jp.k_out = k_out;
  jp.x_out = (x_b - x_b.dot(k_out) * k_out).normalized();
  jp.y_out = k_out.cross(jp.x_out);
}

/// The output basis of the chosen kind at direction k; none at a pole.
std::optional<std::array<Vec3, 2>> output_basis(const JonesPupil& jp, const Vec3& k) {
  return jp.basis == JonesBasis::DoublePole ? double_pole_basis(jp.x_out, jp.y_out, jp.k_out, k)
                                            : dipole_basis(jp.y_out, k);
}

/// Sizes the per-point columns for n points: NaN in j and in every map.
void resize_points(JonesPupil& jp, std::size_t n) {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  jp.px.resize(n);
  jp.py.resize(n);
  jp.status.resize(n);
  jp.w.resize(n);
  jp.j.assign(4 * n, Complex(nan, nan));
  for (auto* column : {&jp.transmission, &jp.diattenuation, &jp.diattenuation_orientation,
                       &jp.diattenuation_ellipticity, &jp.retardance, &jp.retardance_orientation,
                       &jp.retardance_ellipticity}) {
    column->assign(n, nan);
  }
}

/// J and the maps of point i that arrived (ADR 0034, points 5 and 6): ray i of the traced
/// `rays` with start direction k0, jp.w[i] already set; J divided by `field_factor` (option
/// field amplitude, 1 without it; NaN: j stays NaN, the maps are filled). False at a pole of a
/// basis, where the point keeps NaN.
bool fill_point(JonesPupil& jp,
                const trace::RayBatch& rays,
                std::size_t i,
                const Vec3& k0,
                double field_factor = 1.0) {
  const Vec3 k = direction(rays, i);
  const std::optional<std::array<Vec3, 2>> in = double_pole_basis(jp.x_in, jp.y_in, jp.k_in, k0);
  const std::optional<std::array<Vec3, 2>> out = output_basis(jp, k);
  if (!in || !out) return false;
  // Dereferenced once, right after the check (bugprone-unchecked-optional-access).
  const std::array<Vec3, 2>& b_in = *in;
  const std::array<Vec3, 2>& b_out = *out;
  const math::CMat3 p = rays.prt_matrix(i);
  const std::array<Complex, 4> j = jones_matrix(p, b_in, b_out);
  // Phase factor exp(-2 pi i W lambda_ref / lambda) (ADR 0034, point 5).
  const double phase_scale =
      -2.0 * std::numbers::pi * jp.reference_wavelength_um / jp.wavelength_um;
  if (!std::isnan(field_factor)) {
    const Complex factor = std::polar(1.0 / field_factor, phase_scale * jp.w[i]);
    for (std::size_t e = 0; e < 4; ++e) jp.j[4 * i + e] = factor * j[e];
  }
  jp.transmission[i] = rays.weight()[i];
  const polar::Diattenuation d = polar::diattenuation(p, k0, k);
  jp.diattenuation[i] = d.value;
  const PolarizationEllipse de = ellipse_in(d.axis, b_in);
  jp.diattenuation_orientation[i] = de.orientation;
  jp.diattenuation_ellipticity[i] = de.ellipticity;
  const polar::Retardance r = polar::physical_retardance(p, rays.q_matrix(i), k0);
  jp.retardance[i] = r.value;
  const PolarizationEllipse re = ellipse_in(r.fast_axis, b_in);
  jp.retardance_orientation[i] = re.orientation;
  jp.retardance_ellipticity[i] = re.ellipticity;
  return true;
}

/// Removes the global phase psi of J(0, 0) from every j (option remove_global_phase).
void remove_global_phase(JonesPupil& jp, const std::array<Complex, 4>& j00) {
  jp.global_phase = global_phase_of(j00);
  const Complex turn = std::polar(1.0, -jp.global_phase);
  for (Complex& z : jp.j) z *= turn;
}

/// Traces `rays` (Q requested). With the field-amplitude option it records the paths of all
/// rays for the power factors: rtt-trace has no recording trace with RunControl, so the trace is
/// one block then, cancellation is checked before it and progress reported after it; the batch
/// ends bitwise as without recording (sequential.hpp).
std::optional<trace::RayPaths> trace_points(const CompiledSystem& system,
                                            PathId path,
                                            trace::RayBatch& rays,
                                            bool record,
                                            const trace::RunControl& control) {
  if (!record) {
    [[maybe_unused]] const auto stats =
        trace::SequentialTracer().trace(system, path, rays, control);
    return std::nullopt;
  }
  trace::RunMonitor monitor(control, rays.size(), "trace");
  if (monitor.stop()) monitor.finish();  // throws Cancelled
  std::vector<std::size_t> all(rays.size());
  std::iota(all.begin(), all.end(), std::size_t{0});
  trace::RayPaths paths;
  [[maybe_unused]] const auto stats =
      trace::SequentialTracer().trace(system, path, rays, paths, std::span<const std::size_t>(all));
  monitor.add(rays.size());
  monitor.finish();
  return paths;
}

/// prod_q sqrt(c_q) of recorded ray r (ADR 0034, point 5, F4): c_q = c_s of every Refract at a
/// Fresnel or coated surface, the only events whose P carries sqrt(c) (apply_event.cpp; a
/// diffraction order inherits P_0 of order 0 with the same tangential invariant; reflections and
/// ideal elements carry none). The same computation as the tracer: the hit from the recorded
/// state before the event (slot q, bitwise the tracer's), the global normal of that hit,
/// polar::tangential_invariant and polar::transmission_power_factors (Byrnes, Eqs. (21), (22)).
/// NaN if c_s and c_p differ by more than kEqualPowerFactors (absorbing media).
double field_factor(const CompiledSystem& system,
                    PathId path,
                    std::uint16_t wavelength,
                    const trace::RayPaths& paths,
                    std::size_t r) {
  const auto& events = system.path(path).events;
  double product = 1.0;
  for (std::size_t q = 0; q < events.size(); ++q) {
    const compile::CompiledEvent& e = events[q];
    if (e.kind != model::EventKind::Refract) continue;
    const compile::CompiledSurface& surface = system.surfaces()[e.surface];
    const bool fresnel = std::holds_alternative<model::Fresnel>(surface.interaction);
    const bool coated =
        std::holds_alternative<model::CoatingRef>(surface.interaction) && surface.coating;
    if (!fresnel && !coated) continue;
    trace::RayState in;
    in.pos = paths.position_at(r, q);
    in.dir = paths.direction_at(r, q);
    const trace::SurfaceHit hit = trace::intersect_surface(in, surface);
    const Vec3 normal = surface.to_global.apply_vector(hit.normal);
    const Complex n_i = system.media()[e.medium_before].index[wavelength];
    const Complex n_t = system.media()[e.medium_after].index[wavelength];
    const double xi = polar::tangential_invariant(in.dir, normal, n_i);
    const polar::PowerFactors<double> c = polar::transmission_power_factors(n_i, n_t, xi);
    if (!(std::abs(c.s - c.p) <= kEqualPowerFactors * std::abs(c.s))) {
      return std::numeric_limits<double>::quiet_NaN();
    }
    product *= c.s;
  }
  return std::sqrt(product);
}

/// "pupil.field_amplitude_absorbing" with the number of points, if any (ADR 0034, point 5).
void warn_absorbing(JonesPupil& jp, std::size_t absorbing) {
  if (absorbing == 0) return;
  constexpr diagnostics::DiagnosticCode kCode = "pupil.field_amplitude_absorbing";
  jp.warnings.push_back({kCode.severity(), path_pointer(jp.path),
                         std::to_string(absorbing) +
                             " point(s) of the Jones pupil cross an absorbing interface "
                             "(c_s != c_p): the field amplitude is no scalar multiple of J, j is "
                             "NaN there",
                         std::string(kCode.str())});
}

/// "pupil.chief_phase_unavailable": the global phase was requested, but the chief ray has no J.
void warn_chief_phase(JonesPupil& jp) {
  constexpr diagnostics::DiagnosticCode kCode = "pupil.chief_phase_unavailable";
  jp.warnings.push_back({kCode.severity(), path_pointer(jp.path),
                         "the chief ray has no J (field amplitude NaN): the global phase is not "
                         "removed",
                         std::string(kCode.str())});
}

/// "pupil.basis_pole" with the number of points at a pole, if any (ADR 0034, point 5).
void warn_poles(JonesPupil& jp, std::size_t poles) {
  if (poles == 0) return;
  constexpr diagnostics::DiagnosticCode kCode = "pupil.basis_pole";
  jp.warnings.push_back({kCode.severity(), path_pointer(jp.path),
                         std::to_string(poles) +
                             " point(s) of the Jones pupil within 1e-6 rad of a pole of the "
                             "basis: J and the maps are NaN there",
                         std::string(kCode.str())});
}

}  // namespace

std::array<Vec3, 2> input_chief_frame(const Vec3& k0_c) noexcept {
  if (k0_c == Vec3(0.0, 0.0, -1.0)) {
    // R(z -> -z) is not unique; the rotation by pi about x (ADR 0034, point 5).
    return {Vec3::UnitX(), Vec3(0.0, -1.0, 0.0)};
  }
  const math::Mat3 r = trace::rotation_between(Vec3::UnitZ(), k0_c);
  return {r * Vec3::UnitX(), r * Vec3::UnitY()};
}

std::optional<std::array<Vec3, 2>> double_pole_basis(const Vec3& x_c,
                                                     const Vec3& y_c,
                                                     const Vec3& k_c,
                                                     const Vec3& k) {
  // Angle between k and the pole -k_c: |k + k_c| = 2 sin(angle / 2).
  if ((k + k_c).norm() <= 2.0 * std::sin(kPoleAngle / 2.0)) return std::nullopt;
  const math::Mat3 r = trace::rotation_between(k_c, k);
  return std::array<Vec3, 2>{r * x_c, r * y_c};
}

std::optional<std::array<Vec3, 2>> dipole_basis(const Vec3& y_c, const Vec3& k) {
  const Vec3 c = y_c.cross(k);
  const double s = c.norm();  // sine of the angle between k and the pole axis
  if (!(s > std::sin(kPoleAngle))) return std::nullopt;
  const Vec3 x = c / s;
  return std::array<Vec3, 2>{x, k.cross(x)};
}

PolarizationEllipse polarization_ellipse(Complex a, Complex b) noexcept {
  // Stokes parameters of the project (docs/architecture.md, "Händigkeit und Stokes").
  const double s0 = std::norm(a) + std::norm(b);
  if (!(s0 > 0.0)) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    return {nan, nan};
  }
  const double s1 = std::norm(a) - std::norm(b);
  const double s2 = 2.0 * (std::conj(a) * b).real();
  const double s3 = 2.0 * (a * std::conj(b)).imag();
  double psi = 0.5 * std::atan2(s2, s1);  // in (-pi/2, pi/2]
  if (psi < 0.0) psi += std::numbers::pi;
  if (psi >= std::numbers::pi) psi -= std::numbers::pi;
  const double ratio = std::clamp(s3 / s0, -1.0, 1.0);
  return {psi, 0.5 * std::asin(ratio)};
}

JonesPupil jones_pupil(const CompiledSystem& system,
                       PathId path,
                       std::uint16_t field,
                       std::uint16_t wavelength,
                       const JonesPupilOptions& options) {
  return jones_pupil(system, path, field, wavelength, options, trace::RunControl{});
}

JonesPupil jones_pupil(const CompiledSystem& system,
                       PathId path,
                       std::uint16_t field,
                       std::uint16_t wavelength,
                       const JonesPupilOptions& options,
                       const trace::RunControl& control) {
  detail::check_path(system, path);
  detail::check_wavelength(system, wavelength);
  detail::LossCounter losses(system, path, options.lost_warning_fraction);
  reject_crystal_modes(system, path);  // before any ray is traced

  // Reference sphere and chief ray of the reference wavelength, as opd_points().
  const detail::Reference ref =
      detail::make_reference(system, path, field, options.aiming, 0.0, nullptr);
  const std::uint32_t image = detail::image_surface(system, path);
  const Vec3 x_b = system.surfaces()[image].to_global.apply_vector(Vec3::UnitX());
  check_grazing(system, image, x_b, ref.chief_dir);

  JonesPupil jp;
  jp.path = path;
  jp.field = field;
  jp.wavelength = wavelength;
  const std::uint16_t reference = system.reference_wavelength();
  jp.wavelength_um = system.wavelengths_um()[wavelength];
  jp.reference_wavelength_um = system.wavelengths_um()[reference];
  jp.basis = options.basis;
  jp.remove_global_phase = options.remove_global_phase;
  jp.field_amplitude = options.field_amplitude;
  jp.sphere = ref.sphere;

  // Chief frames, both from the chief ray of the reference wavelength (addendum #208).
  const std::uint16_t fields[] = {field};
  const trace::RayBatch chief_ref_start = trace::make_rays(
      system, path, fields, reference, trace::SinglePupilPoint{0.0, 0.0}, options.aiming);
  set_frames(jp, direction(chief_ref_start, 0), ref.chief_dir, x_b);

  // The chief ray of the evaluated wavelength: the reference of W (as opd_points()) and J(0, 0).
  const double n_image = detail::image_index(system, path, wavelength);
  const trace::RayBatch chief_start = trace::make_rays(
      system, path, fields, wavelength, trace::SinglePupilPoint{0.0, 0.0}, options.aiming);
  const trace::RayBatch chief =
      detail::chief_ray(system, path, field, wavelength, options.aiming, nullptr);
  const double opl_chief =
      detail::chief_opl(position(chief, 0), direction(chief, 0), chief.opl()[0], ref, n_image);

  // The rays of opd_points(), traced with Q requested (bitwise the same otherwise, #207).
  trace::RayBatch rays =
      trace::make_rays(system, path, fields, wavelength, options.sampling, options.aiming, control);
  const std::vector<Vec3> k0 = start_directions(rays);
  rays.enable_q();
  const std::optional<trace::RayPaths> paths =
      trace_points(system, path, rays, options.field_amplitude, control);

  const std::size_t n = rays.size();
  resize_points(jp, n);
  std::size_t poles = 0;
  std::size_t absorbing = 0;
  for (std::size_t i = 0; i < n; ++i) {
    const OpdPoint q = detail::point_against(position(rays, i), direction(rays, i), rays.opl()[i],
                                             detail::arrived(rays, i, image), rays.status()[i],
                                             opl_chief, ref, n_image);
    losses.add_ray(q.status, rays.last_surface()[i]);
    jp.px[i] = rays.pupil_x()[i];
    jp.py[i] = rays.pupil_y()[i];
    jp.status[i] = q.status;
    jp.w[i] = q.w;
    if (q.status != trace::RayStatus::Alive) continue;
    const double ff = paths ? field_factor(system, path, wavelength, *paths, i) : 1.0;
    if (std::isnan(ff)) ++absorbing;
    if (!fill_point(jp, rays, i, k0[i], ff)) ++poles;  // unreachable in this form, see jones.hpp
  }

  bool chief_phase_unavailable = false;
  if (options.remove_global_phase) {
    // J(0, 0) of the chief ray of the evaluated wavelength (W = 0 there by construction). The
    // field amplitude scales it by a positive factor, which leaves psi; it only matters whether
    // the chief ray has a J at all.
    bool chief_has_j = true;
    if (options.field_amplitude) {
      trace::RayBatch again = chief_start;
      // The same virtual reference ray as `chief` (ADR 0032, addendum #229): a central
      // obscuration does not take its J away.
      again.set_ignore_obscurations(true);
      const std::optional<trace::RayPaths> chief_path =
          trace_points(system, path, again, true, trace::RunControl{});
      // field_factor() reads the slots of an arrived ray only: after a loss they are NaN.
      chief_has_j = chief_path && detail::arrived(again, 0, image) &&
                    !std::isnan(field_factor(system, path, wavelength, *chief_path, 0));
    }
    const Vec3 k0_c = direction(chief_start, 0);
    const Vec3 k_c = direction(chief, 0);
    const auto in = double_pole_basis(jp.x_in, jp.y_in, jp.k_in, k0_c);
    const auto out = output_basis(jp, k_c);
    if (!chief_has_j) {
      chief_phase_unavailable = true;
    } else if (in && out) {
      remove_global_phase(jp, jones_matrix(chief.prt_matrix(0), *in, *out));
    }
  }

  jp.losses = losses.result();
  jp.warnings = losses.warnings();
  warn_poles(jp, poles);
  warn_absorbing(jp, absorbing);
  if (chief_phase_unavailable) warn_chief_phase(jp);
  return jp;
}

JonesPupil jones_pupil(const CompiledSystem& system,
                       PathId path,
                       const trace::RayBatch& start,
                       std::size_t chief,
                       const JonesPupilOptions& options) {
  return jones_pupil(system, path, start, chief, options, trace::RunControl{});
}

JonesPupil jones_pupil(const CompiledSystem& system,
                       PathId path,
                       const trace::RayBatch& start,
                       std::size_t chief,
                       const JonesPupilOptions& options,
                       const trace::RunControl& control) {
  detail::check_path(system, path);
  if (start.size() == 0) throw std::invalid_argument("analysis: no start rays for the Jones pupil");
  if (chief >= start.size()) {
    throw std::invalid_argument("analysis: chief ray index " + std::to_string(chief) +
                                " outside the " + std::to_string(start.size()) + " start rays");
  }
  const std::uint16_t wavelength = start.wl()[0];
  detail::check_wavelength(system, wavelength);
  if (std::any_of(start.wl().begin(), start.wl().end(),
                  [wavelength](std::uint16_t wl) { return wl != wavelength; })) {
    throw std::invalid_argument("analysis: the start rays of a Jones pupil need one wavelength");
  }
  detail::LossCounter losses(system, path, options.lost_warning_fraction);
  reject_crystal_modes(system, path);  // before any ray is traced

  trace::RayBatch rays = start;
  const std::vector<Vec3> k0 = start_directions(rays);
  rays.enable_q();
  const std::optional<trace::RayPaths> paths =
      trace_points(system, path, rays, options.field_amplitude, control);

  // The chief ray gives both frames and the reference of the phase (addendum #208, PR C).
  const std::uint32_t image = detail::image_surface(system, path);
  if (!detail::arrived(rays, chief, image)) {
    detail::throw_lost(system, rays, chief, rays.field()[chief],
                       "analysis: the chief ray of the Jones pupil does not arrive at the image "
                       "surface");
  }
  const Vec3 k_c = direction(rays, chief);
  const Vec3 x_b = system.surfaces()[image].to_global.apply_vector(Vec3::UnitX());
  check_grazing(system, image, x_b, k_c);

  JonesPupil jp;
  jp.path = path;
  jp.field = rays.field()[chief];
  jp.wavelength = wavelength;
  jp.wavelength_um = system.wavelengths_um()[wavelength];
  jp.reference_wavelength_um = system.wavelengths_um()[system.reference_wavelength()];
  jp.basis = options.basis;
  jp.remove_global_phase = options.remove_global_phase;
  jp.field_amplitude = options.field_amplitude;
  jp.sphere.centre = position(rays, chief);
  jp.sphere.radius = std::numeric_limits<double>::quiet_NaN();  // no sphere in this form
  set_frames(jp, k0[chief], k_c, x_b);

  // phi = 2 pi (OPL - OPL_c) / lambda up to the image surface, so w = (OPL_c - OPL) / lambda_ref
  // in waves at lambda_ref with the sign of W (addendum #208, PR C).
  const double lambda_ref_mm = 1e-3 * jp.reference_wavelength_um;
  const double opl_chief = rays.opl()[chief];
  const std::size_t n = rays.size();
  resize_points(jp, n);
  std::size_t poles = 0;
  std::size_t absorbing = 0;
  for (std::size_t i = 0; i < n; ++i) {
    trace::RayStatus status = rays.status()[i];
    if (status == trace::RayStatus::Alive && rays.last_surface()[i] != image) {
      status = trace::RayStatus::Vignetted;  // ended Alive elsewhere (as the fans)
    }
    losses.add_ray(status, rays.last_surface()[i]);
    jp.px[i] = rays.pupil_x()[i];
    jp.py[i] = rays.pupil_y()[i];
    jp.status[i] = status;
    if (status != trace::RayStatus::Alive) {
      jp.w[i] = 0.0;  // as opd_points(): masked by the status
      continue;
    }
    jp.w[i] = (opl_chief - rays.opl()[i]) / lambda_ref_mm;
    const double ff = paths ? field_factor(system, path, wavelength, *paths, i) : 1.0;
    if (std::isnan(ff)) ++absorbing;
    if (!fill_point(jp, rays, i, k0[i], ff)) ++poles;
  }

  // J(0, 0) is the chief ray's: w = 0 there exactly, never at a pole (k0 = k_in, k = k_c); with
  // the field amplitude it is NaN if the chief ray crosses an absorbing interface.
  const bool chief_phase_unavailable =
      options.remove_global_phase && std::isnan(jp.j[4 * chief].real());
  if (options.remove_global_phase && !chief_phase_unavailable) {
    remove_global_phase(
        jp, {jp.j[4 * chief], jp.j[4 * chief + 1], jp.j[4 * chief + 2], jp.j[4 * chief + 3]});
  }

  jp.losses = losses.result();
  jp.warnings = losses.warnings();
  warn_poles(jp, poles);
  warn_absorbing(jp, absorbing);
  if (chief_phase_unavailable) warn_chief_phase(jp);
  return jp;
}

}  // namespace rtt::analysis
