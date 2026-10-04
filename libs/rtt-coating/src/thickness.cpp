#include "rtt/coating/thickness.hpp"

#include <cmath>
#include <stdexcept>
#include <type_traits>
#include <variant>

namespace rtt::coating {

double physical_thickness_um(const Thickness& thickness, double index_at_design) {
  return std::visit(
      [&](const auto& t) -> double {
        if constexpr (std::is_same_v<std::decay_t<decltype(t)>, PhysicalThickness>) {
          if (!std::isfinite(t.um) || t.um < 0.0) {
            throw std::invalid_argument("coating: a physical thickness must be finite and >= 0");
          }
          return t.um;
        } else {
          if (!std::isfinite(t.count) || t.count < 0.0) {
            throw std::invalid_argument("coating: a QWOT count must be finite and >= 0");
          }
          if (!std::isfinite(t.design_wavelength_um) || !(t.design_wavelength_um > 0.0)) {
            throw std::invalid_argument("coating: a design wavelength must be finite and > 0 um");
          }
          if (!std::isfinite(index_at_design) || !(index_at_design > 0.0)) {
            throw std::invalid_argument("coating: a QWOT needs a finite index > 0");
          }
          return quarter_wave_thickness_um(t.count, t.design_wavelength_um, index_at_design);
        }
      },
      thickness);
}

}  // namespace rtt::coating
