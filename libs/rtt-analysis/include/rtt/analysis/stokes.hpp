#pragma once

/// @file stokes.hpp
/// Stokes parameters and the polarization map of the pupil for a given input state (ADR 0034,
/// point 6; #211), from a Jones pupil (jones.hpp), so for both of its forms.
///
/// Conventions:
/// - Input state E_in = e1 x_in + e2 y_in in the input chief frame of the Jones pupil,
///   normalised to |E_in| = 1. Output field E_out = J E_in, components (a, b) in the output
///   basis B_out (double pole or dipole) of each point.
/// - Stokes parameters of the project (docs/architecture.md, "Händigkeit und Stokes";
///   polar::stokes, Lam, Table 2.1): s0 = |a|^2 + |b|^2, s1 = |a|^2 - |b|^2, s2 = 2 Re(a* b),
///   s3 = 2 Im(a b*); s3 > 0 is right circular ((i e1 + e2) / sqrt(2)).
/// - Scale (coordinator's decision for #211): power for the polarized input, with every factor
///   of the ray (ADR 0021, point 2): polarized power = weight |P E|^2 / (||P_T||_F^2 / 2), and
///   ||P_T||_F = ||J||_F, so S = (weight / (||J||_F^2 / 2)) Stokes(J E_in). weight / (||J||_F^2 /
///   2) is the polarization-independent factor s of the ray (volume absorption, efficiencies);
///   for a lossless path S0 = |J E_in|^2. Only for a power-normalised Jones pupil: with field
///   amplitudes the factor has no meaning.
/// - J is deterministic, so the output is fully polarized (degree of polarization 1); the
///   unpolarized case is the transmission map of the Jones pupil.

#include <complex>
#include <cstdint>
#include <vector>

#include "rtt/analysis/jones.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/math/types.hpp"
#include "rtt/model/validate.hpp"
#include "rtt/trace/ray_batch.hpp"

namespace rtt::analysis {

/// Stokes parameters and polarization ellipse of the output per pupil point.
struct StokesPupil {
  compile::PathId path{};
  std::uint16_t field = 0;                ///< as JonesPupil::field
  std::uint16_t wavelength = 0;           ///< as JonesPupil::wavelength
  double wavelength_um = 0.0;             ///< evaluated wavelength in vacuum, um
  std::complex<double> e1{1.0, 0.0};      ///< input component along x_in, normalised
  std::complex<double> e2{0.0, 0.0};      ///< input component along y_in, normalised
  math::Vec3 x_in = math::Vec3::UnitX();  ///< input chief frame (as JonesPupil)
  math::Vec3 y_in = math::Vec3::UnitY();
  math::Vec3 x_out = math::Vec3::UnitX();  ///< output chief frame (as JonesPupil)
  math::Vec3 y_out = math::Vec3::UnitY();
  // Per point, in the order of the Jones pupil:
  std::vector<double> px, py;
  std::vector<trace::RayStatus> status;  ///< as JonesPupil::status
  /// Stokes parameters, power (see the file comment); NaN where the point has no J.
  std::vector<double> s0, s1, s2, s3;
  /// Orientation psi in [0, pi) of the output ellipse from the first towards the second output
  /// basis vector, rad; NaN without J or without output power.
  std::vector<double> orientation;
  /// Ellipticity chi in [-pi/4, pi/4] of the output ellipse, > 0 right-handed, rad.
  std::vector<double> ellipticity;
  /// The warnings of the Jones pupil (ADR 0023, 0034); the Stokes map adds none.
  std::vector<model::Diagnostic> warnings;
};

/// Stokes parameters of E_out = J E_in at every point of `jones` for the input state
/// E_in = (e1 x_in + e2 y_in) / |(e1, e2)| (ADR 0034, point 6; #211).
/// @throws std::invalid_argument if (e1, e2) is zero or not finite, or if `jones` holds field
///         amplitudes (JonesPupil::field_amplitude; the power scale needs the power-normalised J)
[[nodiscard]] StokesPupil stokes_pupil(const JonesPupil& jones,
                                       std::complex<double> e1,
                                       std::complex<double> e2);

}  // namespace rtt::analysis
