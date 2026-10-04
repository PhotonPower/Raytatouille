#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>

#include "rtt/material/air.hpp"
#include "rtt/material/material.hpp"
#include "rtt/material/thermal.hpp"

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using rtt::material::ciddor_air_index;
using rtt::material::kStandardAtmospherePa;
using rtt::material::MaterialLibrary;
using rtt::material::schott_delta_n_abs;
using rtt::material::SchottThermalCoefficients;
using rtt::math::Complex;

namespace {

/// Reference values of the NIST Engineering Metrology Toolbox Ciddor calculator
/// (emtoolbox.nist.gov/Wavelength/Ciddor.asp), queried 2026-10-04 for 0 % relative humidity and
/// 450 umol/mol CO2: vacuum wavelength in nm, temperature in degC, pressure in kPa, index.
/// The calculator prints 9 decimals, so the tolerance is 1e-9 (docs/quellen.md, #25).
struct NistCase {
  double wavelength_nm, temperature_c, pressure_kpa, index;
};
constexpr NistCase kNist[] = {
    {633, 20, 101.325, 1.0002718},    {400, 20, 101.325, 1.000277925},
    {1550, 20, 101.325, 1.000268586}, {300, 20, 101.325, 1.000286581},
    {1700, 20, 101.325, 1.000268479}, {633, 0, 101.325, 1.000291767},
    {633, 40, 101.325, 1.000254398},  {633, -20, 101.325, 1.000314912},
    {633, 20, 50, 1.000134099},       {633, 20, 120, 1.000321916},
    {500, 35, 80, 1.000205918},       {1000, 10, 110, 1.000302923},
};

}  // namespace

TEST_CASE("Ciddor air index matches the NIST calculator", "[air]") {
  for (const NistCase& c : kNist) {
    INFO(c.wavelength_nm << " nm, " << c.temperature_c << " degC, " << c.pressure_kpa << " kPa");
    REQUIRE_THAT(
        ciddor_air_index(c.wavelength_nm / 1000.0, c.temperature_c, c.pressure_kpa / 101.325),
        WithinAbs(c.index, 1e-9));
  }
}

TEST_CASE("AIR is Ciddor air at the temperature and pressure of the medium", "[air]") {
  const MaterialLibrary lib;
  const auto air = lib.resolve("AIR");
  // 1 atm = 101.325 kPa.
  REQUIRE_THAT(air->index(0.633, 20.0, 1.0).real(), WithinAbs(1.0002718, 1e-9));
  REQUIRE(air->index(0.633, 20.0, 1.0).imag() == 0.0);
  REQUIRE_THAT(air->index(0.633, 40.0, 1.0).real(), WithinAbs(1.000254398, 1e-9));
  REQUIRE_THAT(air->index(0.633, 20.0, 120.0 / 101.325).real(), WithinAbs(1.000321916, 1e-9));
  // No pressure, no air: the refractivity vanishes.
  REQUIRE(air->index(0.633, 20.0, 0.0) == Complex(1.0, 0.0));
  // No wavelength range (#25, P6): Ciddor is extrapolated outside 0.3 um to 1.7 um.
  REQUIRE_FALSE(air->wavelength_range_um().has_value());
  // VACUUM stays exactly 1.
  REQUIRE(lib.resolve("VACUUM")->index(0.633, 40.0, 2.0) == Complex(1.0, 0.0));
  REQUIRE(kStandardAtmospherePa == 101325.0);
}

TEST_CASE("Schott Delta n_abs (TIE-19 Eq. (3)) by hand", "[air][thermal]") {
  // Invented coefficients; n = 1.5 so (n^2 - 1)/(2n) = 1.25/3, lambda = 1 um, dT = 10 K.
  const SchottThermalCoefficients c{1e-6, 1e-8, 1e-11, 1e-6, 1e-9, 0.2, 20.0};
  const double polynomial = 1e-6 * 10 + 1e-8 * 100 + 1e-11 * 1000;   // D0 dT + D1 dT^2 + D2 dT^3
  const double resonance = (1e-6 * 10 + 1e-9 * 100) / (1.0 - 0.04);  // / (l^2 - l_tk^2)
  REQUIRE_THAT(schott_delta_n_abs(1.5, 1.0, 10.0, c),
               WithinRel(1.25 / 3.0 * (polynomial + resonance), 1e-12));
  // dT = 0: no change.
  REQUIRE(schott_delta_n_abs(1.5, 1.0, 0.0, c) == 0.0);
  // Negative lambda_tk: S_tk = -1, the denominator becomes l^2 + l_tk^2 (Ansys convention).
  SchottThermalCoefficients neg = c;
  neg.lambda_tk = -0.2;
  const double resonance_neg = (1e-6 * 10 + 1e-9 * 100) / (1.0 + 0.04);
  REQUIRE_THAT(schott_delta_n_abs(1.5, 1.0, 10.0, neg),
               WithinRel(1.25 / 3.0 * (polynomial + resonance_neg), 1e-12));
}
