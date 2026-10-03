#pragma once

/// @file param.hpp
/// Numeric design parameter (docs/architecture.md, "Konfigurationen und Variablen").

#include <optional>
#include <string>

namespace rtt::model {

/// A numeric value that the optimizer may vary or that may be derived from another
/// parameter (pickup). Units are those of the field that holds the parameter.
struct Param {
  double value = 0.0;
  /// True if the optimizer may change this value.
  bool variable = false;
  /// Pickup expression (evaluated from M5 on); empty means an independent value.
  std::optional<std::string> pickup;

  Param() = default;
  /// Implicit on purpose so that `Param r = 51.68;` reads naturally.
  Param(double v) : value(v) {}  // NOLINT(google-explicit-constructor)

  /// True if this parameter carries nothing but a plain value.
  [[nodiscard]] bool is_plain() const noexcept { return !variable && !pickup.has_value(); }

  bool operator==(const Param&) const = default;
};

}  // namespace rtt::model
