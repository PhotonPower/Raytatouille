#pragma once

/// @file path.hpp
/// Explicit ray paths (sequences) through the surface tree.

#include <cstdint>
#include <string>
#include <vector>

#include "rtt/model/ids.hpp"

namespace rtt::model {

enum class EventKind : std::uint8_t {
  Refract,        ///< refraction into the next medium
  Reflect,        ///< reflection
  Transmit,       ///< pass without refraction (thin element, detector)
  Ordinary,       ///< ordinary ray in a uniaxial crystal
  Extraordinary,  ///< extraordinary ray in a uniaxial crystal
  Diffract,       ///< diffraction into `order`
};

struct Event {
  SurfaceId surface;
  EventKind kind = EventKind::Refract;
  int order = 0;  ///< diffraction order, only meaningful for Diffract
  bool operator==(const Event&) const = default;
};

/// Named ordered list of events. `automatic` = all surfaces in tree order, refracting
/// (reflecting at mirrors); the explicit list is then empty.
struct Path {
  std::string name;
  bool automatic = false;
  std::vector<Event> events;
  bool operator==(const Path&) const = default;
};

}  // namespace rtt::model
