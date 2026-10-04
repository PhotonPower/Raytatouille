#pragma once

/// @file thermal.hpp
/// Temperature dependence of the absolute refractive index of glass (SCHOTT model).
///
/// Source (docs/quellen.md): SCHOTT, "TIE-19: Temperature Coefficient of the Refractive Index",
/// Eq. (3) and (4); denominator lambda^2 - S_tk lambda_tk^2 with S_tk = sign(lambda_tk) as in the
/// Ansys Zemax OpticStudio User Guide, "Index of Refraction Computation" (identical to TIE-19 for
/// the positive lambda_tk of SCHOTT glasses).

#include "rtt/math/real.hpp"

namespace rtt::material {

/// Coefficients of the TD record of an AGF catalogue (in this order), see agf.hpp.
struct SchottThermalCoefficients {
  double d0 = 0.0;                      ///< 1/K
  double d1 = 0.0;                      ///< 1/K^2
  double d2 = 0.0;                      ///< 1/K^3
  double e0 = 0.0;                      ///< um^2/K
  double e1 = 0.0;                      ///< um^2/K^2
  double lambda_tk = 0.0;               ///< um
  double reference_temperature_c = 20;  ///< T_ref in degC
  bool operator==(const SchottThermalCoefficients&) const = default;
};

/// Change of the absolute index from T_ref to T (TIE-19 Eq. (3)):
///   Delta n_abs = (n^2 - 1) / (2 n) * (D0 dT + D1 dT^2 + D2 dT^3
///                                      + (E0 dT + E1 dT^2) / (lambda^2 - S_tk lambda_tk^2)),
/// dT = T - T_ref. TIE-19 allows the relative index at T_ref for n.
/// @param n_ref         index at the reference temperature (relative to air, see TIE-19)
/// @param wavelength_um vacuum wavelength in um (TIE-19: "wavelength ... in a vacuum")
/// @param delta_t       T - T_ref in K
template <math::Real T>
[[nodiscard]] T schott_delta_n_abs(T n_ref,
                                   T wavelength_um,
                                   T delta_t,
                                   const SchottThermalCoefficients& c) {
  const T s_tk = c.lambda_tk < 0.0 ? T(-1) : T(1);
  const T denominator = wavelength_um * wavelength_um - s_tk * T(c.lambda_tk) * T(c.lambda_tk);
  const T polynomial = delta_t * (T(c.d0) + delta_t * (T(c.d1) + delta_t * T(c.d2)));
  const T resonance = delta_t * (T(c.e0) + delta_t * T(c.e1)) / denominator;
  return (n_ref * n_ref - T(1)) / (T(2) * n_ref) * (polynomial + resonance);
}

}  // namespace rtt::material
