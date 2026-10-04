#pragma once

/// @file thickness.hpp
/// Layer thickness as physical thickness in um or as quarter-wave optical thickness (QWOT) at a
/// design wavelength.

#include <variant>

#include "rtt/math/real.hpp"

namespace rtt::coating {

/// Physical thickness.
struct PhysicalThickness {
  double um = 0.0;  ///< um, >= 0
};

/// Quarter-wave optical thickness: `count` quarter waves at the vacuum wavelength
/// `design_wavelength_um` and normal incidence, i.e. optical thickness n d = count lambda0 / 4.
struct QuarterWaves {
  double count = 1.0;                 ///< number of quarter waves, >= 0 (1 = lambda/4)
  double design_wavelength_um = 0.0;  ///< design vacuum wavelength lambda0, um, > 0
};

/// Thickness of a layer in a coating design.
using Thickness = std::variant<PhysicalThickness, QuarterWaves>;

/// Physical thickness of `count` quarter waves: d = count lambda0 / (4 n) (optical thickness
/// n d = count lambda0 / 4 at normal incidence; then delta = 2 pi n d / lambda0 = count pi / 2,
/// Byrnes Eq. (8)).
/// @param count         number of quarter waves, dimensionless
/// @param wavelength_um design vacuum wavelength lambda0, um
/// @param index         real part of the layer index at lambda0
/// @return thickness, um
template <math::Real T>
[[nodiscard]] T quarter_wave_thickness_um(T count, T wavelength_um, T index) noexcept {
  return count * wavelength_um / (T(4) * index);
}

/// Physical thickness in um of a design thickness.
/// @param thickness       physical thickness or QWOT
/// @param index_at_design real part of the layer index at the design wavelength (QWOT only)
/// @throws std::invalid_argument for a negative or non-finite thickness or count, a design
///         wavelength that is not finite and > 0, or (QWOT) an index that is not finite and > 0
[[nodiscard]] double physical_thickness_um(const Thickness& thickness, double index_at_design);

}  // namespace rtt::coating
