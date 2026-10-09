#pragma once

/// @file param.hpp
/// Numeric design parameter (docs/architecture.md, "Parameter"; ADR 0029).

#include <optional>
#include <string>
#include <utility>

namespace rtt::model {

/// A numeric value that the optimizer may vary, or a reference to a row of the parameter table
/// (System::parameters, ADR 0029). Units are those of the field that holds the parameter.
///
/// - Unbound (`param` empty): `value` is the value; `variable` marks an unknown of the
///   optimizer; `min`/`max` bound it (same unit as `value`; bounds on a value that is not
///   variable are allowed and have no effect).
/// - Bound (`param` set): the value comes from the row of that name. `value` is meaningless
///   until the table is evaluated (the reader sets 0; validate and the writers do not read it),
///   and a bound Param has neither `variable` nor bounds (validate: param.bound_conflict).
struct Param {
  double value = 0.0;
  /// True if the optimizer may change this value (unbound Params only).
  bool variable = false;
  /// Lower and upper bound for the optimizer (unbound Params only), unit of `value`.
  std::optional<double> min;
  std::optional<double> max;
  /// Name of a row of System::parameters; replaces the pickup of schema 0.3 (ADR 0029).
  std::optional<std::string> param;

  Param() = default;
  /// Implicit on purpose so that `Param r = 51.68;` reads naturally.
  Param(double v) : value(v) {}  // NOLINT(google-explicit-constructor)

  /// A Param bound to the row `name` of the parameter table (value 0, not variable).
  [[nodiscard]] static Param bound(std::string name) {
    Param p;
    p.param = std::move(name);
    return p;
  }

  /// True if the value comes from a row of the parameter table.
  [[nodiscard]] bool is_bound() const noexcept { return param.has_value(); }

  /// True if this parameter carries nothing but a plain value.
  [[nodiscard]] bool is_plain() const noexcept {
    return !variable && !min.has_value() && !max.has_value() && !param.has_value();
  }

  bool operator==(const Param&) const = default;
};

}  // namespace rtt::model
