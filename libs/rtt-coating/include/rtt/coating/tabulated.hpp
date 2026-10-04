#pragma once

/// @file tabulated.hpp
/// Coating given as a table of amplitudes r_s, r_p, t_s, t_p over angle of incidence and vacuum
/// wavelength (e.g. computed elsewhere or measured), interpolated bilinearly.

#include <cstddef>
#include <vector>

#include "rtt/coating/transfer_matrix.hpp"

namespace rtt::coating {

/// Tabulated coating for fixed incident and exit media. Amplitudes as in transfer_matrix.hpp
/// (basis [s, p, k] of rtt-polar, Convention A).
///
/// The angle grid covers the whole range [0, pi/2] (checked at construction), so every ray finds
/// a value; the wavelength range is checked against the system wavelengths when the system is
/// compiled (#61), like for materials. Evaluation is noexcept: between grid points the real and
/// imaginary parts are interpolated bilinearly; outside the wavelength range the nearest
/// wavelength column is used (clamped).
class TabulatedCoating {
 public:
  /// @param angles_rad     angles of incidence in the incident medium, rad, strictly
  ///                       increasing, first 0, last pi/2 (both within 1e-12 rad)
  /// @param wavelengths_um vacuum wavelengths, um, strictly increasing, at least one, > 0
  /// @param values         amplitudes for every (angle i, wavelength j) at index
  ///                       i * wavelengths_um.size() + j
  /// @throws std::invalid_argument if a grid is too short, not finite, not strictly increasing,
  ///         the angles do not cover [0, pi/2], a wavelength is not > 0, the number of values is
  ///         not angles x wavelengths, or a value is not finite
  TabulatedCoating(std::vector<double> angles_rad,
                   std::vector<double> wavelengths_um,
                   std::vector<Amplitudes<double>> values);

  /// Interpolated amplitudes.
  /// @param angle_rad     angle of incidence, rad (clamped to [0, pi/2])
  /// @param wavelength_um vacuum wavelength, um (clamped to the table range, see class comment)
  [[nodiscard]] Amplitudes<double> amplitudes(double angle_rad,
                                              double wavelength_um) const noexcept;

  /// Smallest tabulated wavelength, um.
  [[nodiscard]] double min_wavelength_um() const noexcept { return wavelengths_um_.front(); }
  /// Largest tabulated wavelength, um.
  [[nodiscard]] double max_wavelength_um() const noexcept { return wavelengths_um_.back(); }

 private:
  [[nodiscard]] const Amplitudes<double>& at(std::size_t angle, std::size_t wavelength) const {
    return values_[angle * wavelengths_um_.size() + wavelength];
  }

  std::vector<double> angles_rad_;
  std::vector<double> wavelengths_um_;
  std::vector<Amplitudes<double>> values_;
};

}  // namespace rtt::coating
