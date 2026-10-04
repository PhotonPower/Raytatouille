#include "rtt/material/tabulated.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace rtt::material {

TabulatedMaterial::TabulatedMaterial(std::vector<IndexSample> samples)
    : samples_(std::move(samples)) {
  if (samples_.size() < 2) {
    throw std::invalid_argument("TabulatedMaterial: at least 2 samples are required");
  }
  for (std::size_t i = 0; i < samples_.size(); ++i) {
    const IndexSample& s = samples_[i];
    if (!std::isfinite(s.wavelength_um) || !std::isfinite(s.n) || !std::isfinite(s.kappa)) {
      throw std::invalid_argument("TabulatedMaterial: samples must be finite");
    }
    if (!(s.wavelength_um > 0.0) || !(s.n > 0.0) || s.kappa < 0.0) {
      throw std::invalid_argument(
          "TabulatedMaterial: samples need wavelength > 0 um, n > 0 and kappa >= 0");
    }
    if (i > 0 && !(samples_[i - 1].wavelength_um < s.wavelength_um)) {
      throw std::invalid_argument("TabulatedMaterial: wavelengths must increase strictly");
    }
    // kappa = -0 counts as 0; store +0 so that the sign of the imaginary part stays >= 0.
    samples_[i].kappa = s.kappa + 0.0;
  }
}

math::Complex TabulatedMaterial::index(double wavelength_um,
                                       double /*temperature_c*/,
                                       double /*pressure_atm*/) const {
  const IndexSample& first = samples_.front();
  const IndexSample& last = samples_.back();
  if (!(wavelength_um > first.wavelength_um)) return {first.n, first.kappa};
  if (!(wavelength_um < last.wavelength_um)) return {last.n, last.kappa};
  // First sample with a wavelength above the argument; its predecessor lies at or below it.
  const auto upper =
      std::upper_bound(samples_.begin(), samples_.end(), wavelength_um,
                       [](double wl, const IndexSample& s) { return wl < s.wavelength_um; });
  const IndexSample& b = *upper;
  const IndexSample& a = *(upper - 1);
  // Linear interpolation in the wavelength (convention of #23), n and kappa separately.
  const double t = (wavelength_um - a.wavelength_um) / (b.wavelength_um - a.wavelength_um);
  return {a.n + t * (b.n - a.n), a.kappa + t * (b.kappa - a.kappa)};
}

std::optional<WavelengthRange> TabulatedMaterial::wavelength_range_um() const {
  return WavelengthRange{samples_.front().wavelength_um, samples_.back().wavelength_um};
}

}  // namespace rtt::material
