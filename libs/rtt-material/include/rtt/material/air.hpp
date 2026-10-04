#pragma once

/// @file air.hpp
/// Refractive index of air after Ciddor and the AIR material.
///
/// Source (docs/quellen.md): NIST Engineering Metrology Toolbox, "Refractive Index of Air –
/// Documentation", Appendix A-III "Ciddor Calculation of the Index of Refraction",
/// Eq. (A21)-(A41). Decided for #25: dry air (water vapour mole fraction x_v = 0) with a CO2
/// content of 450 umol/mol. The vacuum wavelength is in um, the temperature in degC and the
/// pressure in atm (NIST computes in Pa). The result is the phase refractive index of air
/// (absolute).

#include "rtt/material/material.hpp"
#include "rtt/math/real.hpp"

namespace rtt::material {

/// One standard atmosphere in Pa (NIST Eq. (A27), p_R1).
inline constexpr double kStandardAtmospherePa = 101325.0;

/// CO2 content of the AIR material in umol/mol (#25).
inline constexpr double kAirCo2Ppm = 450.0;

/// Phase refractive index of dry air with kAirCo2Ppm CO2 after Ciddor (NIST Eq. (A21)-(A41)
/// with x_v = 0, so all water vapour terms vanish).
/// Ciddor's range of validity (NIST calculator): 0.3 um to 1.7 um, -40 degC to 100 degC,
/// 10 kPa to 140 kPa; outside the formula is extrapolated without warning (#25, P6).
/// @param wavelength_um vacuum wavelength in um
/// @param temperature_c air temperature in degC
/// @param pressure_atm  air pressure in atm (converted to Pa as NIST requires)
/// @return phase refractive index of air, absolute (dimensionless)
template <math::Real T>
[[nodiscard]] T ciddor_air_index(T wavelength_um, T temperature_c, T pressure_atm) {
  // NIST works in Pa: 1 atm = p_R1 = 101325 Pa (A27).
  const T pressure_pa = pressure_atm * T(kStandardAtmospherePa);
  // (A22) dispersion constants of standard air
  const T k0 = T(238.0185);
  const T k1 = T(5792105.0);
  const T k2 = T(57.362);
  const T k3 = T(167917.0);
  // (A23) compressibility constants (dry air: b and c terms vanish with x_v = 0)
  const T a0 = T(1.58123e-6);
  const T a1 = T(-2.9331e-8);
  const T a2 = T(1.1043e-10);
  // (A26) second virial term (e vanishes with x_v = 0)
  const T d = T(1.83e-11);
  // (A27) reference conditions of standard air, (A28) its compressibility, (A30) gas constant
  const T p_r1 = T(kStandardAtmospherePa);
  const T t_r1 = T(288.15);
  const T z_a = T(0.9995922115);
  const T r_gas = T(8.314472);

  // (A31) S = 1/lambda^2 in um^-2
  const T s = T(1) / (wavelength_um * wavelength_um);
  // (A32) refractivity of standard air
  const T r_as = T(1e-8) * (k1 / (k0 - s) + k3 / (k2 - s));
  // (A34) molar mass of dry air and (A35) refractivity for the CO2 content
  const T x_co2 = T(kAirCo2Ppm);
  const T m_a = T(0.0289635) + T(1.2011e-8) * (x_co2 - T(400));
  const T r_axs = r_as * (T(1) + T(5.34e-7) * (x_co2 - T(450)));
  // (A36) Kelvin temperature, (A37) compressibility with x_v = 0
  const T t = temperature_c;
  const T t_k = t + T(273.15);
  const T p_over_t = pressure_pa / t_k;
  const T z_m = T(1) - p_over_t * (a0 + a1 * t + a2 * t * t) + p_over_t * p_over_t * d;
  // (A38) density of standard dry air, (A40) density of the dry air component (x_v = 0)
  const T rho_axs = p_r1 * m_a / (z_a * r_gas * t_r1);
  const T rho_a = pressure_pa * m_a / (z_m * r_gas * t_k);
  // (A41) without the water vapour term
  return T(1) + (rho_a / rho_axs) * r_axs;
}

/// The AIR material: dry air after Ciddor at the temperature and pressure of the medium
/// (ciddor_air_index). No wavelength range is reported (#25, P6): outside 0.3 um to 1.7 um the
/// formula is extrapolated.
class AirMaterial final : public Material {
 public:
  /// Ciddor index n + 0i; see Material::index.
  [[nodiscard]] math::Complex index(double wavelength_um,
                                    double temperature_c,
                                    double pressure_atm) const override {
    return {ciddor_air_index(wavelength_um, temperature_c, pressure_atm), 0.0};
  }
};

}  // namespace rtt::material
