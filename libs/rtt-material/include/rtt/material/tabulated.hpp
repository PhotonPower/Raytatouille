#pragma once

/// @file tabulated.hpp
/// Material given by tabulated n and kappa over the vacuum wavelength.

#include <optional>
#include <vector>

#include "rtt/material/material.hpp"

namespace rtt::material {

/// One measured point: vacuum wavelength in um, real index n and extinction coefficient kappa.
struct IndexSample {
  double wavelength_um = 0.0;
  double n = 1.0;
  double kappa = 0.0;
  bool operator==(const IndexSample&) const = default;
};

/// Complex index n + i*kappa from tabulated samples, interpolated linearly in the wavelength,
/// n and kappa separately. The valid range is [first sample, last sample]; outside it index()
/// returns the value of the nearest end (rtt-compile rejects such wavelengths anyway).
/// index() ignores temperature and pressure.
class TabulatedMaterial final : public Material {
 public:
  /// @param samples at least 2 samples with finite values, strictly increasing wavelength > 0,
  ///                n > 0 and kappa >= 0
  /// @throws std::invalid_argument otherwise
  explicit TabulatedMaterial(std::vector<IndexSample> samples);

  /// Linear interpolation n + i*kappa at the wavelength, see the class comment.
  [[nodiscard]] math::Complex index(double wavelength_um,
                                    double temperature_c,
                                    double pressure_atm) const override;

  [[nodiscard]] std::optional<WavelengthRange> wavelength_range_um() const override;

  [[nodiscard]] const std::vector<IndexSample>& samples() const noexcept { return samples_; }

 private:
  std::vector<IndexSample> samples_;
};

}  // namespace rtt::material
