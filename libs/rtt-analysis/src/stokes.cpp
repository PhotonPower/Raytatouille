#include "rtt/analysis/stokes.hpp"

#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <stdexcept>

namespace rtt::analysis {

StokesPupil stokes_pupil(const JonesPupil& jones,
                         std::complex<double> e1,
                         std::complex<double> e2) {
  if (jones.field_amplitude) {
    throw std::invalid_argument(
        "analysis: the Stokes map needs a power-normalised Jones pupil, not field amplitudes "
        "(its power scale weight / (||J||_F^2 / 2) has no meaning there; #211)");
  }
  const double norm = std::sqrt(std::norm(e1) + std::norm(e2));
  if (!std::isfinite(norm) || !(norm > 0.0)) {
    throw std::invalid_argument(
        "analysis: the input state of the Stokes map must be finite and not zero");
  }
  StokesPupil out;
  out.path = jones.path;
  out.field = jones.field;
  out.wavelength = jones.wavelength;
  out.wavelength_um = jones.wavelength_um;
  out.e1 = e1 / norm;
  out.e2 = e2 / norm;
  out.x_in = jones.x_in;
  out.y_in = jones.y_in;
  out.x_out = jones.x_out;
  out.y_out = jones.y_out;
  out.px = jones.px;
  out.py = jones.py;
  out.status = jones.status;
  out.warnings = jones.warnings;

  const std::size_t n = jones.status.size();
  const double nan = std::numeric_limits<double>::quiet_NaN();
  for (auto* column : {&out.s0, &out.s1, &out.s2, &out.s3, &out.orientation, &out.ellipticity}) {
    column->assign(n, nan);
  }
  for (std::size_t i = 0; i < n; ++i) {
    const std::array<std::complex<double>, 4> j = {jones.j[4 * i], jones.j[4 * i + 1],
                                                   jones.j[4 * i + 2], jones.j[4 * i + 3]};
    double frobenius2 = 0.0;
    for (const std::complex<double>& z : j) frobenius2 += std::norm(z);
    if (!std::isfinite(frobenius2)) continue;  // no J: not arrived, at a pole (ADR 0034)
    // E_out = J E_in in the output basis of the point.
    const std::complex<double> a = j[0] * out.e1 + j[1] * out.e2;
    const std::complex<double> b = j[2] * out.e1 + j[3] * out.e2;
    // Power scale (ADR 0021, point 2): polarized power = weight |P E|^2 / (||P_T||_F^2 / 2) with
    // ||P_T||_F = ||J||_F; the factor is the polarization-independent s of the ray. A path that
    // blocks every state (J = 0) transmits nothing.
    const double scale = frobenius2 > 0.0 ? jones.transmission[i] / (0.5 * frobenius2) : 0.0;
    // Stokes parameters of the project (polar::stokes; Lam, Table 2.1).
    out.s0[i] = scale * (std::norm(a) + std::norm(b));
    out.s1[i] = scale * (std::norm(a) - std::norm(b));
    out.s2[i] = scale * 2.0 * std::real(std::conj(a) * b);
    out.s3[i] = scale * 2.0 * std::imag(a * std::conj(b));
    const PolarizationEllipse ellipse = polarization_ellipse(a, b);
    out.orientation[i] = ellipse.orientation;
    out.ellipticity[i] = ellipse.ellipticity;
  }
  return out;
}

}  // namespace rtt::analysis
