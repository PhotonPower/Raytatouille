#include "rtt/material/dispersion.hpp"

#include <cmath>
#include <stdexcept>

namespace rtt::material {

DispersionMaterial::DispersionMaterial(DispersionFormula formula, WavelengthRange range)
    : formula_(formula), range_(range) {
  if (!std::isfinite(range.min_um) || !std::isfinite(range.max_um) || !(range.min_um > 0.0) ||
      !(range.min_um < range.max_um)) {
    throw std::invalid_argument(
        "DispersionMaterial: the wavelength range must be finite with 0 < min_um < max_um");
  }
}

math::Complex DispersionMaterial::index(double wavelength_um,
                                        double /*temperature_c*/,
                                        double /*pressure_atm*/) const {
  // Temperature and pressure enter with dn/dT and air (#25); the formula value is absolute
  // until then (convention of #23).
  return {refractive_index(formula_, wavelength_um), 0.0};
}

}  // namespace rtt::material
