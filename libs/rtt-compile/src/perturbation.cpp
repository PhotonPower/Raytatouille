#include "rtt/compile/perturbation.hpp"

#include <array>
#include <cstddef>
#include <optional>

#include "rtt/material/air.hpp"

namespace rtt::compile {

bool RigidMotion::is_identity() const {
  return (translation_mm.array() == 0.0).all() && (rotation_rad.array() == 0.0).all();
}

double line_vacuum_wavelength_um(double air_um) {
  // lambda_air = lambda_vac / n_air(lambda_vac, 20 degC, 1 atm), the wavelength scale of the
  // catalogues (docs/quellen.md, "Umrechnung von Katalogdaten"; CatalogMaterial in agf.hpp, #25).
  // The map v -> air_um * n_air(v) contracts strongly (|d n_air / d lambda| * lambda ~ 1e-5).
  constexpr int kMaxSteps = 10;
  double vacuum = air_um;
  for (int i = 0; i < kMaxSteps; ++i) {
    const double next = air_um * material::ciddor_air_index(vacuum, 20.0, 1.0);
    if (next == vacuum) break;
    vacuum = next;
  }
  return vacuum;
}

std::optional<double> FraunhoferIndices::abbe_number() const noexcept {
  const double d = dispersion();
  if (d == 0.0) return std::nullopt;
  return (nd - 1.0) / d;  // TIE-29, Eq. (2.1-1)
}

FraunhoferIndices fraunhofer_indices(const material::Material& material,
                                     double temperature_c,
                                     double pressure_atm) {
  constexpr std::array<double, 3> kLines{kLineDAirUm, kLineFAirUm, kLineCAirUm};
  const std::optional<material::WavelengthRange> range = material.wavelength_range_um();
  FraunhoferIndices out;
  std::array<double, 3> n{};
  for (std::size_t k = 0; k < kLines.size(); ++k) {
    const double wl = line_vacuum_wavelength_um(kLines[k]);
    if (range && !range->contains(wl)) {
      out.status = FraunhoferIndices::Status::OutOfRange;
      out.line = k;
      out.wavelength_um = wl;
      return out;
    }
    const math::Complex value = material.index(wl, temperature_c, pressure_atm);
    if (value.imag() != 0.0) {
      out.status = FraunhoferIndices::Status::Complex;
      out.line = k;
      out.wavelength_um = wl;
      return out;
    }
    n[k] = value.real();
  }
  out.nd = n[0];
  out.nf = n[1];
  out.nc = n[2];
  return out;
}

}  // namespace rtt::compile
