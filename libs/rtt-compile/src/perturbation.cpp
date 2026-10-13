#include "rtt/compile/perturbation.hpp"

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

// Stubs for the red run of the tests (#260).
std::optional<double> FraunhoferIndices::abbe_number() const noexcept {
  return std::nullopt;
}

FraunhoferIndices fraunhofer_indices(const material::Material&, double, double) {
  return {};
}

}  // namespace rtt::compile
