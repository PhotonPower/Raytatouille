#include "rtt/optim/variables.hpp"

#include <cstddef>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <variant>
#include <vector>

#include "rtt/diagnostics/codes.hpp"
#include "rtt/model/parameters.hpp"

namespace rtt::optim {

std::vector<Variable> collect_variables(const model::System& system) {
  std::vector<Variable> out;
  // Table rows first, in table order (ADR 0030, point 5).
  for (std::size_t i = 0; i < system.parameters.size(); ++i) {
    const model::ParameterRow& row = system.parameters[i];
    if (!row.variable) continue;
    const Bounds bounds{row.min, row.max};
    const std::string base = "/parameters/" + std::to_string(i);
    if (const auto* value = std::get_if<double>(&row.form)) {
      out.push_back({base + "/value", row.name, std::nullopt, std::nullopt, *value, bounds});
    } else if (const auto* values = std::get_if<std::vector<double>>(&row.form)) {
      for (std::size_t k = 0; k < values->size(); ++k) {
        out.push_back({base + "/values/" + std::to_string(k), row.name, k, std::nullopt,
                       (*values)[k], bounds});
      }
    }
    // An expression row is never variable (validate: parameters.variable_expression).
  }
  // Then the model Params, in the order of for_each_param (the order of the edit form); a Param
  // in an inactive node is no variable (ADR 0049, point 4).
  const std::vector<std::string> inactive = model::inactive_nodes(system);
  std::size_t index = 0;
  model::for_each_param(system, [&](std::string_view pointer, const model::Param& param) {
    if (param.variable && !param.is_bound() && !model::in_inactive_node(inactive, pointer)) {
      out.push_back({std::string(pointer) + "/value", std::string(), std::nullopt, index,
                     param.value, Bounds{param.min, param.max}});
    }
    ++index;
  });
  return out;
}

std::vector<model::Diagnostic> inactive_variables(const model::System& system) {
  constexpr diagnostics::DiagnosticCode kCode = "optim.variable_inactive";
  std::vector<model::Diagnostic> out;
  const std::vector<std::string> inactive = model::inactive_nodes(system);
  if (inactive.empty()) return out;
  model::for_each_param(system, [&](std::string_view pointer, const model::Param& param) {
    if (param.variable && !param.is_bound() && model::in_inactive_node(inactive, pointer)) {
      out.push_back({kCode.severity(), std::string(pointer),
                     "the Param is variable but in an inactive node: it is not a variable and "
                     "keeps its value (ADR 0049)",
                     std::string(kCode.str())});
    }
  });
  return out;
}

model::System with_values(const model::System& system,
                          std::span<const Variable> variables,
                          std::span<const double> values) {
  if (variables.size() != values.size()) {
    throw std::invalid_argument("with_values: one value per variable is needed");
  }
  model::System out = system;
  std::vector<std::optional<std::size_t>> by_param;  // Param index -> variable index
  for (std::size_t v = 0; v < variables.size(); ++v) {
    const Variable& var = variables[v];
    if (var.param_index) {
      if (by_param.size() <= *var.param_index) by_param.resize(*var.param_index + 1);
      by_param[*var.param_index] = v;
      continue;
    }
    const std::optional<std::size_t> row = model::find_parameter(out, var.row);
    if (!row) throw std::invalid_argument("with_values: unknown parameter row '" + var.row + "'");
    model::ParameterForm& form = out.parameters[*row].form;
    if (var.configuration) {
      std::get<std::vector<double>>(form).at(*var.configuration) = values[v];
    } else {
      form = values[v];
    }
  }
  if (!by_param.empty()) {
    std::size_t index = 0;
    model::for_each_param(out, [&](std::string_view /*pointer*/, model::Param& param) {
      if (index < by_param.size() && by_param[index]) param.value = values[*by_param[index]];
      ++index;
    });
  }
  return out;
}

}  // namespace rtt::optim
