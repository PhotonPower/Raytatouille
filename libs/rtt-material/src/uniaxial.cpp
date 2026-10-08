#include "rtt/material/uniaxial.hpp"

#include <algorithm>
#include <stdexcept>
#include <string>
#include <utility>

namespace rtt::material {

UniaxialMaterial::UniaxialMaterial(std::shared_ptr<const Material> ordinary,
                                   std::shared_ptr<const Material> extraordinary)
    : ordinary_(std::move(ordinary)), extraordinary_(std::move(extraordinary)) {
  if (ordinary_ == nullptr) throw std::invalid_argument("UniaxialMaterial: ordinary is null");
  if (extraordinary_ == nullptr) {
    throw std::invalid_argument("UniaxialMaterial: extraordinary is null");
  }
}

math::Complex UniaxialMaterial::n_ordinary(double wavelength_um,
                                           double temperature_c,
                                           double pressure_atm) const {
  return ordinary_->index(wavelength_um, temperature_c, pressure_atm);
}

math::Complex UniaxialMaterial::n_extraordinary(double wavelength_um,
                                                double temperature_c,
                                                double pressure_atm) const {
  return extraordinary_->index(wavelength_um, temperature_c, pressure_atm);
}

std::optional<WavelengthRange> UniaxialMaterial::wavelength_range_um() const {
  const std::optional<WavelengthRange> o = ordinary_->wavelength_range_um();
  const std::optional<WavelengthRange> e = extraordinary_->wavelength_range_um();
  if (!o) return e;
  if (!e) return o;
  return WavelengthRange{std::max(o->min_um, e->min_um), std::min(o->max_um, e->max_um)};
}

UniaxialMaterial MaterialLibrary::resolve_uniaxial(std::string_view ordinary,
                                                   std::string_view extraordinary) const {
  // Each part through resolve() (same cache and errors); the message names the failing part.
  const auto part = [this](std::string_view reference, const char* name) {
    try {
      return resolve(reference);
    } catch (const UnknownMaterial& e) {
      throw UnknownMaterial(std::string(name) + ": " + e.what());
    }
  };
  // Ordinary first, as documented; two statements make the order explicit.
  auto o = part(ordinary, "ordinary");
  auto e = part(extraordinary, "extraordinary");
  return {std::move(o), std::move(e)};
}

}  // namespace rtt::material
