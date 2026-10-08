#pragma once

/// @file uniaxial.hpp
/// Uniaxial crystal media (ADR 0014, point 4; ADR 0026, point 1): the two principal indices n_O
/// (ordinary) and n_E (extraordinary) from two ordinary materials, e.g. the Zemax convention of
/// two catalogue glasses X and X-E. The optic axis is not part of the material; it belongs to
/// the element (ADR 0026, point 2).
///
/// Conventions as in material.hpp: vacuum wavelengths in um, temperatures in degree Celsius,
/// pressures in atm, absolute complex indices n + i*kappa with kappa >= 0 for absorption.

#include <memory>
#include <optional>

#include "rtt/material/material.hpp"
#include "rtt/math/types.hpp"

namespace rtt::material {

/// Uniaxial crystal: principal indices from two shared, immutable materials. A value type:
/// copies share the two materials; safe to use from several threads like Material.
///
/// Absorption is passed through unchanged: a crystal with kappa != 0 at a system wavelength is
/// an error of rtt-compile (code crystal.absorbing, ADR 0026), not of this class.
class UniaxialMaterial {
 public:
  /// @param ordinary      material of the ordinary principal index n_O
  /// @param extraordinary material of the extraordinary principal index n_E
  /// @throws std::invalid_argument if a pointer is null
  UniaxialMaterial(std::shared_ptr<const Material> ordinary,
                   std::shared_ptr<const Material> extraordinary);

  // Copies only: a move would leave null parts behind and break the invariant of the
  // constructor; with the copy operations declared, a "move" copies the two shared pointers.
  UniaxialMaterial(const UniaxialMaterial&) = default;
  UniaxialMaterial& operator=(const UniaxialMaterial&) = default;
  ~UniaxialMaterial() = default;

  /// Ordinary principal index n_O + i*kappa_O, absolute: ordinary().index(...).
  /// @param wavelength_um vacuum wavelength in micrometre
  /// @param temperature_c medium temperature in degree Celsius
  /// @param pressure_atm  medium pressure in atm
  [[nodiscard]] math::Complex n_ordinary(double wavelength_um,
                                         double temperature_c,
                                         double pressure_atm) const;

  /// Extraordinary principal index n_E + i*kappa_E, absolute: extraordinary().index(...).
  /// @param wavelength_um vacuum wavelength in micrometre
  /// @param temperature_c medium temperature in degree Celsius
  /// @param pressure_atm  medium pressure in atm
  [[nodiscard]] math::Complex n_extraordinary(double wavelength_um,
                                              double temperature_c,
                                              double pressure_atm) const;

  /// Wavelengths on which both parts are valid, in um: the intersection of their ranges;
  /// std::nullopt if both are unbounded. For disjoint ranges the result is empty
  /// (min_um > max_um), so that WavelengthRange::contains() is false for every wavelength.
  [[nodiscard]] std::optional<WavelengthRange> wavelength_range_um() const;

  /// The material of n_O.
  [[nodiscard]] const Material& ordinary() const noexcept { return *ordinary_; }
  /// The material of n_E.
  [[nodiscard]] const Material& extraordinary() const noexcept { return *extraordinary_; }

 private:
  std::shared_ptr<const Material> ordinary_;
  std::shared_ptr<const Material> extraordinary_;
};

}  // namespace rtt::material
