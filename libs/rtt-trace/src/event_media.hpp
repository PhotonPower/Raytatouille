#pragma once

/// @file event_media.hpp
/// Private: the physical inputs of a compiled event, shared by the sequential tracer and ray
/// aiming (stop_hit), so that further event data (crystals, #132) are filled in one place.

#include <cstdint>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/trace/apply_event.hpp"

namespace rtt::trace::detail {

/// Complex indices of the media before, after and beyond `event` and the vacuum wavelength at
/// wavelength index `wl` (ADR 0021); no coating layers (the caller adds them where needed).
[[nodiscard]] inline EventMedia event_media(const compile::CompiledSystem& system,
                                            const compile::CompiledEvent& event,
                                            std::uint16_t wl) noexcept {
  EventMedia m;
  m.before = system.media()[event.medium_before].index[wl];
  m.after = system.media()[event.medium_after].index[wl];
  m.beyond = system.media()[event.medium_beyond].index[wl];
  m.wavelength_um = system.wavelengths_um()[wl];
  return m;
}

}  // namespace rtt::trace::detail
