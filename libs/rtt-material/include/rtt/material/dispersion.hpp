#pragma once

/// @file dispersion.hpp
/// Dispersion formulas n(lambda) and a Material that evaluates one of them.
///
/// All formulas take a wavelength lambda in micrometre (um); their coefficients carry the
/// matching powers of um. They return the real index n. DispersionMaterial passes the vacuum
/// wavelength and returns the formula value as absolute index; CatalogMaterial (agf.hpp) passes
/// the wavelength in air because catalogue data are relative to air (#25). Temperature and
/// pressure are not part of these formulas (thermal.hpp, air.hpp). The formula kernels are free
/// function templates over rtt::math::Real (ADR 0006, 0014, 0015); DispersionMaterial calls them
/// with double.
///
/// Sources (docs/quellen.md): Ansys Zemax OpticStudio User Guide, Release 2025 R1, "The Glass
/// Dispersion Formulas" (one page per formula, each with a single unnumbered equation, lambda in
/// um); Cauchy: J. M. Palmer, "Models for fitting refractive index n vs. lambda", University of
/// Arizona, p. 1, equation CAUCHY.

#include <array>
#include <cmath>
#include <cstddef>
#include <optional>
#include <variant>

#include "rtt/material/material.hpp"
#include "rtt/math/real.hpp"

namespace rtt::material {

/// Schott formula: n^2 = a0 + a1 l^2 + a2 l^-2 + a3 l^-4 + a4 l^-6 + a5 l^-8
/// (OpticStudio User Guide, "The Schott Formula").
struct SchottCoefficients {
  /// a0 dimensionless, a1 in um^-2, a2 in um^2, a3 in um^4, a4 in um^6, a5 in um^8
  std::array<double, 6> a{};
  bool operator==(const SchottCoefficients&) const = default;
};

/// Sellmeier sum with N terms: n^2 - 1 = sum_i K_i l^2 / (l^2 - L_i)
/// (OpticStudio User Guide, "The Sellmeier 1 Formula" with N = 3, "The Sellmeier 3 Formula"
/// with N = 4, "The Sellmeier 5 Formula" with N = 5). L_i in um^2.
template <std::size_t N>
struct SellmeierCoefficients {
  std::array<double, N> k{};  ///< K_1 .. K_N, dimensionless
  std::array<double, N> l{};  ///< L_1 .. L_N in um^2
  bool operator==(const SellmeierCoefficients&) const = default;
};
using Sellmeier1Coefficients = SellmeierCoefficients<3>;
using Sellmeier3Coefficients = SellmeierCoefficients<4>;
using Sellmeier5Coefficients = SellmeierCoefficients<5>;

/// Sellmeier 2: n^2 - 1 = A + B1 l^2 / (l^2 - lambda1^2) + B2 / (l^2 - lambda2^2)
/// (OpticStudio User Guide, "The Sellmeier 2 Formula"). Note the squared lambda1, lambda2 and
/// the missing l^2 in the numerator of the second term.
struct Sellmeier2Coefficients {
  double a = 0.0;        ///< dimensionless
  double b1 = 0.0;       ///< dimensionless
  double lambda1 = 0.0;  ///< um
  double b2 = 0.0;       ///< um^2
  double lambda2 = 0.0;  ///< um
  bool operator==(const Sellmeier2Coefficients&) const = default;
};

/// Sellmeier 4: n^2 = A + B l^2 / (l^2 - C) + D l^2 / (l^2 - E)
/// (OpticStudio User Guide, "The Sellmeier 4 Formula"). Note n^2, not n^2 - 1.
struct Sellmeier4Coefficients {
  double a = 0.0;  ///< dimensionless
  double b = 0.0;  ///< dimensionless
  double c = 0.0;  ///< um^2
  double d = 0.0;  ///< dimensionless
  double e = 0.0;  ///< um^2
  bool operator==(const Sellmeier4Coefficients&) const = default;
};

/// Herzberger: n = A + B L + C L^2 + D l^2 + E l^4 + F l^6 with L = 1 / (l^2 - 0.028)
/// (OpticStudio User Guide, "The Herzberger Formula"; 0.028 in um^2).
struct HerzbergerCoefficients {
  double a = 0.0;  ///< dimensionless
  double b = 0.0;  ///< um^2
  double c = 0.0;  ///< um^4
  double d = 0.0;  ///< um^-2
  double e = 0.0;  ///< um^-4
  double f = 0.0;  ///< um^-6
  bool operator==(const HerzbergerCoefficients&) const = default;
};

/// Conrady: n = n0 + A / l + B / l^3.5 (OpticStudio User Guide, "The Conrady Formula").
struct ConradyCoefficients {
  double n0 = 0.0;  ///< dimensionless
  double a = 0.0;   ///< um
  double b = 0.0;   ///< um^3.5
  bool operator==(const ConradyCoefficients&) const = default;
};

/// Cauchy: n = A + B / l^2 + C / l^4 (Palmer, "Models for fitting refractive index n vs.
/// lambda", p. 1, CAUCHY).
struct CauchyCoefficients {
  double a = 0.0;  ///< dimensionless
  double b = 0.0;  ///< um^2
  double c = 0.0;  ///< um^4
  bool operator==(const CauchyCoefficients&) const = default;
};

/// n from n^2; NaN if n^2 < 0. With measured catalogue data this happens only outside the
/// material's wavelength range; coefficients with a pole (l^2 = L_i) inside the range are not
/// detected and give inf or NaN there.
template <math::Real T>
[[nodiscard]] T index_from_square(T n_squared) {
  using std::sqrt;
  return sqrt(n_squared);
}

/// Schott formula, see SchottCoefficients. @param wavelength_um vacuum wavelength in um
template <math::Real T>
[[nodiscard]] T refractive_index(const SchottCoefficients& c, T wavelength_um) {
  const T l2 = wavelength_um * wavelength_um;
  const T inv = T(1) / l2;
  // Horner form of a2 l^-2 + a3 l^-4 + a4 l^-6 + a5 l^-8.
  const T negative = inv * (c.a[2] + inv * (c.a[3] + inv * (c.a[4] + inv * c.a[5])));
  return index_from_square(c.a[0] + c.a[1] * l2 + negative);
}

/// Sellmeier sum, see SellmeierCoefficients. @param wavelength_um vacuum wavelength in um
template <math::Real T, std::size_t N>
[[nodiscard]] T refractive_index(const SellmeierCoefficients<N>& c, T wavelength_um) {
  const T l2 = wavelength_um * wavelength_um;
  T n2 = T(1);
  for (std::size_t i = 0; i < N; ++i) n2 += c.k[i] * l2 / (l2 - c.l[i]);
  return index_from_square(n2);
}

/// Sellmeier 2, see Sellmeier2Coefficients. @param wavelength_um vacuum wavelength in um
template <math::Real T>
[[nodiscard]] T refractive_index(const Sellmeier2Coefficients& c, T wavelength_um) {
  const T l2 = wavelength_um * wavelength_um;
  return index_from_square(T(1) + c.a + c.b1 * l2 / (l2 - c.lambda1 * c.lambda1) +
                           c.b2 / (l2 - c.lambda2 * c.lambda2));
}

/// Sellmeier 4, see Sellmeier4Coefficients. @param wavelength_um vacuum wavelength in um
template <math::Real T>
[[nodiscard]] T refractive_index(const Sellmeier4Coefficients& c, T wavelength_um) {
  const T l2 = wavelength_um * wavelength_um;
  return index_from_square(c.a + c.b * l2 / (l2 - c.c) + c.d * l2 / (l2 - c.e));
}

/// Herzberger, see HerzbergerCoefficients. @param wavelength_um vacuum wavelength in um
template <math::Real T>
[[nodiscard]] T refractive_index(const HerzbergerCoefficients& c, T wavelength_um) {
  const T l2 = wavelength_um * wavelength_um;
  const T big_l = T(1) / (l2 - T(0.028));
  return c.a + c.b * big_l + c.c * big_l * big_l + l2 * (c.d + l2 * (c.e + l2 * c.f));
}

/// Conrady, see ConradyCoefficients. @param wavelength_um vacuum wavelength in um
template <math::Real T>
[[nodiscard]] T refractive_index(const ConradyCoefficients& c, T wavelength_um) {
  using std::pow;
  return c.n0 + c.a / wavelength_um + c.b / pow(wavelength_um, T(3.5));
}

/// Cauchy, see CauchyCoefficients. @param wavelength_um vacuum wavelength in um
template <math::Real T>
[[nodiscard]] T refractive_index(const CauchyCoefficients& c, T wavelength_um) {
  const T inv2 = T(1) / (wavelength_um * wavelength_um);
  return c.a + inv2 * (c.b + inv2 * c.c);
}

/// One of the supported dispersion formulas with its coefficients.
using DispersionFormula = std::variant<SchottCoefficients,
                                       Sellmeier1Coefficients,
                                       Sellmeier2Coefficients,
                                       Sellmeier3Coefficients,
                                       Sellmeier4Coefficients,
                                       Sellmeier5Coefficients,
                                       HerzbergerCoefficients,
                                       ConradyCoefficients,
                                       CauchyCoefficients>;

/// Real index n(lambda) of a formula. @param wavelength_um vacuum wavelength in um
template <math::Real T>
[[nodiscard]] T refractive_index(const DispersionFormula& formula, T wavelength_um) {
  return std::visit([&](const auto& c) { return refractive_index(c, wavelength_um); }, formula);
}

/// Non-absorbing material given by a dispersion formula on a wavelength range.
/// index() is absolute and independent of temperature and pressure (user-defined data, #25 P5).
class DispersionMaterial final : public Material {
 public:
  /// @param formula formula and coefficients (lambda in um)
  /// @param range   valid vacuum wavelengths in um
  /// @throws std::invalid_argument if the range is not finite with 0 < min_um < max_um
  DispersionMaterial(DispersionFormula formula, WavelengthRange range);

  /// n(lambda) + 0i; see Material::index. NaN if n^2 < 0 (outside the range only).
  [[nodiscard]] math::Complex index(double wavelength_um,
                                    double temperature_c,
                                    double pressure_atm) const override;

  [[nodiscard]] std::optional<WavelengthRange> wavelength_range_um() const override {
    return range_;
  }

  /// The formula and coefficients given at construction (lambda in um).
  [[nodiscard]] const DispersionFormula& formula() const noexcept { return formula_; }

 private:
  DispersionFormula formula_;
  WavelengthRange range_;
};

}  // namespace rtt::material
