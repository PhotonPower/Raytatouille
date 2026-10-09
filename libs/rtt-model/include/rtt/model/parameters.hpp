#pragma once

/// @file parameters.hpp
/// Evaluation of the parameter table and its configurations (ADR 0029): row values per
/// configuration, the resolved system of one configuration, and a visitor over every Param.

#include <cstddef>
#include <functional>
#include <optional>
#include <string_view>
#include <vector>

#include "rtt/model/param.hpp"
#include "rtt/model/system.hpp"
#include "rtt/model/validate.hpp"

namespace rtt::model {

/// Number C of configurations (ADR 0029, point 1): the size of System::configurations, or 1 (the
/// nominal configuration) if there are none. Always at least 1.
[[nodiscard]] std::size_t configuration_count(const System& system) noexcept;

/// Index of the configuration called `name`, if any (the nominal configuration has no name).
[[nodiscard]] std::optional<std::size_t> find_configuration(const System& system,
                                                            std::string_view name);

/// Index of the first row of System::parameters called `name`, if any.
[[nodiscard]] std::optional<std::size_t> find_parameter(const System& system,
                                                        std::string_view name);

/// Values of all rows of the parameter table in all configurations, in the unit of each row.
struct ParameterValues {
  std::size_t rows = 0;            ///< number of rows (System::parameters)
  std::size_t configurations = 1;  ///< C, see configuration_count()
  /// Row-major: values[row * configurations + configuration]. NaN in every configuration for a
  /// row that is not defined: a syntax or name error in its expression, a `values` row with the
  /// wrong length, or an expression that uses a row with an error or with a non-finite value. A
  /// non-finite value itself (from the file or an expression) stays as it is.
  std::vector<double> values;

  /// The value of `row` in `configuration`.
  /// @pre row < rows and configuration < configurations
  [[nodiscard]] double at(std::size_t row, std::size_t configuration) const {
    return values[row * configurations + configuration];
  }
};

/// Evaluates the parameter table for every configuration (ADR 0029, point 2): rows in table
/// order, an expression only over earlier rows, each operation in double in the order of the
/// grammar without FMA contraction, so that the result is the same on every platform.
///
/// Appends to `diagnostics` (if not null) the evaluating codes, at /parameters/i/expression:
/// parameters.expression_syntax, parameters.unknown_name, parameters.forward_reference (with the
/// position in the expression) and parameters.not_finite (once per configuration in which the
/// result is not finite, with the configuration in the message). A row that uses an undefined
/// row is undefined too, without a diagnostic of its own. The structural codes of the table
/// (names, values_count, bounds, value.not_finite) come from validate(), which also calls this.
/// Never throws (except std::bad_alloc).
[[nodiscard]] ParameterValues evaluate_parameters(const System& system,
                                                  std::vector<Diagnostic>* diagnostics = nullptr);

/// The system of one configuration (ADR 0029, point 3): a copy in which every Param bound to a
/// row is replaced by an unbound Param with the row's value in `configuration` (variable false,
/// no bounds, `param` empty). The table and the configurations stay in the copy; unbound Params
/// are unchanged. Without bound Params the copy equals the input.
/// @throws std::invalid_argument if `configuration` >= configuration_count(system), or if the
///         table or a binding has an error (any diagnostic of validate() with a code
///         parameters.*, param.* or configurations.*, or value.not_finite in the table);
///         compile() runs validate() first, so that it reports these as CompileError
[[nodiscard]] System resolve_parameters(const System& system, std::size_t configuration);

/// Calls `visit` for every Param of the model with its JSON pointer into the edit form
/// (rtt::io::to_edit_json, ADR 0024), in the order in which the edit form writes them: the
/// object distance, the system aperture, then the tree in pre-order (pose position 0..2 and
/// rotation_deg 0..2 of each node; per surface its pose, the shape base radius, conic and
/// coefficients, the Zernike terms, then the phase layers). Rows of the parameter table are not
/// Params and are not visited. The order is the order of the optimizer's model variables
/// (ADR 0030, point 5).
void for_each_param(System& system,
                    const std::function<void(std::string_view pointer, Param& param)>& visit);

/// for_each_param() for reading.
void for_each_param(const System& system,
                    const std::function<void(std::string_view pointer, const Param& param)>& visit);

}  // namespace rtt::model
