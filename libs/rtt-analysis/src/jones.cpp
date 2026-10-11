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
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "common.hpp"
#include "reference.hpp"
#include "rtt/math/types.hpp"
#include "rtt/polar/prt_analysis.hpp"
#include "rtt/trace/apply_event.hpp"
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
  const detail::Reference ref = detail::make_reference(system, path, field, options.aiming, 0.0);
  const std::uint32_t image = detail::image_surface(system, path);
  const Vec3 x_b = system.surfaces()[image].to_global.apply_vector(Vec3::UnitX());
  if (x_b.cross(ref.chief_dir).norm() <= kGrazingLimit) {
    throw AnalysisError(
        "analysis: the chief ray grazes the image surface (its x axis is parallel to the chief "
        "ray); the output basis of the Jones pupil has no frame",
        system.surfaces()[image].location);
  }

  JonesPupil jp;
  jp.path = path;
  jp.field = field;
  jp.wavelength = wavelength;
  const std::uint16_t reference = system.reference_wavelength();
  jp.wavelength_um = system.wavelengths_um()[wavelength];
  jp.reference_wavelength_um = system.wavelengths_um()[reference];
  jp.basis = options.basis;
  jp.remove_global_phase = options.remove_global_phase;
  jp.sphere = ref.sphere;

  // Chief frames, both from the chief ray of the reference wavelength (addendum #208).
  const std::uint16_t fields[] = {field};
  const trace::RayBatch chief_ref_start = trace::make_rays(
      system, path, fields, reference, trace::SinglePupilPoint{0.0, 0.0}, options.aiming);
  jp.k_in = direction(chief_ref_start, 0);
  const std::array<Vec3, 2> frame_in = input_chief_frame(jp.k_in);
  jp.x_in = frame_in[0];
  jp.y_in = frame_in[1];
  jp.k_out = ref.chief_dir;
  jp.x_out = (x_b - x_b.dot(ref.chief_dir) * ref.chief_dir).normalized();
  jp.y_out = ref.chief_dir.cross(jp.x_out);
  const auto output_basis = [&](const Vec3& k) {
    return options.basis == JonesBasis::DoublePole
               ? double_pole_basis(jp.x_out, jp.y_out, jp.k_out, k)
               : dipole_basis(jp.y_out, k);
  };

  // The chief ray of the evaluated wavelength: the reference of W (as opd_points()) and J(0, 0).
  const double n_image = detail::image_index(system, path, wavelength);
  const trace::RayBatch chief_start = trace::make_rays(
      system, path, fields, wavelength, trace::SinglePupilPoint{0.0, 0.0}, options.aiming);
  const trace::RayBatch chief = detail::chief_ray(system, path, field, wavelength, options.aiming);
  const double opl_chief =
      detail::chief_opl(position(chief, 0), direction(chief, 0), chief.opl()[0], ref, n_image);

  // The rays of opd_points(), traced with Q requested (bitwise the same otherwise, #207).
  trace::RayBatch rays =
      trace::make_rays(system, path, fields, wavelength, options.sampling, options.aiming, control);
  const std::vector<Vec3> k0 = start_directions(rays);
  rays.enable_q();
  [[maybe_unused]] const auto stats = trace::SequentialTracer().trace(system, path, rays, control);

  const double nan = std::numeric_limits<double>::quiet_NaN();
  const std::size_t n = rays.size();
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
  // Phase factor exp(-2 pi i W lambda_ref / lambda) (ADR 0034, point 5).
  const double phase_scale =
      -2.0 * std::numbers::pi * jp.reference_wavelength_um / jp.wavelength_um;
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
    const Vec3 k = direction(rays, i);
    const std::optional<std::array<Vec3, 2>> in =
        double_pole_basis(jp.x_in, jp.y_in, jp.k_in, k0[i]);
    const std::optional<std::array<Vec3, 2>> out = output_basis(k);
    if (!in || !out) continue;  // at a pole: NaN (unreachable in this form, see jones.hpp)
    // Dereferenced once, right after the check (bugprone-unchecked-optional-access).
    const std::array<Vec3, 2>& b_in = *in;
    const std::array<Vec3, 2>& b_out = *out;
    const math::CMat3 p = rays.prt_matrix(i);
    const std::array<Complex, 4> j = jones_matrix(p, b_in, b_out);
    const Complex factor = std::polar(1.0, phase_scale * q.w);
    for (std::size_t e = 0; e < 4; ++e) jp.j[4 * i + e] = factor * j[e];
    jp.transmission[i] = rays.weight()[i];
    const polar::Diattenuation d = polar::diattenuation(p, k0[i], k);
    jp.diattenuation[i] = d.value;
    const PolarizationEllipse de = ellipse_in(d.axis, b_in);
    jp.diattenuation_orientation[i] = de.orientation;
    jp.diattenuation_ellipticity[i] = de.ellipticity;
    const polar::Retardance r = polar::physical_retardance(p, rays.q_matrix(i), k0[i]);
    jp.retardance[i] = r.value;
    const PolarizationEllipse re = ellipse_in(r.fast_axis, b_in);
    jp.retardance_orientation[i] = re.orientation;
    jp.retardance_ellipticity[i] = re.ellipticity;
  }

  if (options.remove_global_phase) {
    // J(0, 0) of the chief ray of the evaluated wavelength (W = 0 there by construction).
    const Vec3 k0_c = direction(chief_start, 0);
    const Vec3 k_c = direction(chief, 0);
    const auto in = double_pole_basis(jp.x_in, jp.y_in, jp.k_in, k0_c);
    const auto out = output_basis(k_c);
    if (in && out) {
      jp.global_phase = global_phase_of(jones_matrix(chief.prt_matrix(0), *in, *out));
      const Complex turn = std::polar(1.0, -jp.global_phase);
      for (Complex& z : jp.j) z *= turn;
    }
  }

  jp.losses = losses.result();
  jp.warnings = losses.warnings();
  return jp;
}

}  // namespace rtt::analysis
