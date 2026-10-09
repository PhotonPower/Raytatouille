#pragma once

/// @file variables.hpp
/// The unknowns of an optimization run (ADR 0030, point 5; ADR 0029): variable rows of the
/// parameter table and variable model Params, with their pointers into the edit form for the
/// result patch (ADR 0024) and their bounds.

#include <cstddef>
#include <optional>
#include <span>
#include <string>
#include <vector>

#include "rtt/model/system.hpp"
#include "rtt/optim/bounds.hpp"

namespace rtt::optim {

/// One unknown of the optimizer.
struct Variable {
  /// JSON pointer of the value in the edit form (ADR 0024) that the result patch replaces:
  /// "<pointer of the Param>/value" for a model Param, "/parameters/i/value" or
  /// "/parameters/i/values/k" for a table row.
  std::string pointer;
  std::string row;                           ///< name of the table row; empty for a model Param
  std::optional<std::size_t> configuration;  ///< column k of a `values` row; empty otherwise
  /// Position of the Param in the order of model::for_each_param; empty for a table row.
  std::optional<std::size_t> param_index;
  double start = 0.0;  ///< value in the input system, in the unit of the field or row
  Bounds bounds;       ///< min and max of the Param or the row
};

/// The variables of `system` in the order of ADR 0030, point 5: first the variable rows of the
/// parameter table in table order (a `value` row is one unknown for all configurations, a
/// `values` row one unknown per configuration in column order), then the variable model Params
/// in the order of model::for_each_param. Expression rows and bound Params are never variable
/// (validate: parameters.variable_expression, param.bound_conflict) and are skipped.
[[nodiscard]] std::vector<Variable> collect_variables(const model::System& system);

/// A copy of `system` with `variables` set to `values`: the `value` or `values[k]` of a table
/// row, the `value` of a model Param. Nothing else changes; variable, bounds and bindings stay.
/// @param variables the result of collect_variables(system)
/// @param values    one value per variable, in the unit of the field or row
/// @throws std::invalid_argument if the sizes differ
[[nodiscard]] model::System with_values(const model::System& system,
                                        std::span<const Variable> variables,
                                        std::span<const double> values);

}  // namespace rtt::optim
