#pragma once

/// @file event_media.hpp
/// Private: the physical inputs of a compiled event, shared by the sequential tracer and ray
/// aiming (stop_hit), so that further event data (crystals, #132) are filled in one place.

#include <cstdint>
#include <optional>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/trace/apply_event.hpp"

namespace rtt::trace::detail {

/// Crystal data of `medium` at wavelength index `wl`: n_O (real part of `index`), n_E and the
/// global optic axis (ADR 0026, point 3); none for an isotropic medium.
[[nodiscard]] inline std::optional<CrystalSide> crystal_side(const compile::CompiledMedium& medium,
                                                             std::uint16_t wl) noexcept {
  if (!medium.is_crystal()) return std::nullopt;
  return CrystalSide{medium.index[wl].real(), medium.index_extraordinary[wl],
                     medium.optic_axis.value_or(math::Vec3::UnitZ())};
}

/// Complex indices of the media before, after and beyond `event`, the vacuum wavelength at
/// wavelength index `wl` (ADR 0021) and the crystal data of the media before and after with the
/// mode of the event (ADR 0026); no coating layers (the caller adds them where needed).
[[nodiscard]] inline EventMedia event_media(const compile::CompiledSystem& system,
                                            const compile::CompiledEvent& event,
                                            std::uint16_t wl) noexcept {
  EventMedia m;
  m.before = system.media()[event.medium_before].index[wl];
  m.after = system.media()[event.medium_after].index[wl];
  m.beyond = system.media()[event.medium_beyond].index[wl];
  m.wavelength_um = system.wavelengths_um()[wl];
  m.crystal_before = crystal_side(system.media()[event.medium_before], wl);
  m.crystal_after = crystal_side(system.media()[event.medium_after], wl);
  m.crystal_mode = event.crystal_mode;
  return m;
}

}  // namespace rtt::trace::detail
