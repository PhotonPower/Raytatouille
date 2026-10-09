// Evaluation of the parameter table and its configurations (ADR 0029).

#include "rtt/model/parameters.hpp"

#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include "expression.hpp"
#include "rtt/diagnostics/codes.hpp"

namespace rtt::model {
namespace {

constexpr double kUndefined = std::numeric_limits<double>::quiet_NaN();

std::string row_pointer(std::size_t row) {
  return "/parameters/" + std::to_string(row);
}

std::string configuration_text(const System& system, std::size_t k) {
  if (system.configurations.empty()) return "the nominal configuration";
  return "configuration '" + system.configurations[k].name + "'";
}

void add(std::vector<Diagnostic>* out,
         diagnostics::DiagnosticCode code,
         std::string location,
         std::string message) {
  if (out == nullptr) return;
  out->push_back(
      {code.severity(), std::move(location), std::move(message), std::string(code.str())});
}

/// Visits every Param of the model in edit-form order (see for_each_param()); `S` is System or
/// const System, `P` the matching Param type.
template <typename S, typename Visit>
class ParamWalker {
 public:
  explicit ParamWalker(const Visit& visit) : visit_(visit) {}

  void run(S& s) {
    visit_("/object/distance", s.object.distance);
    visit_("/aperture/value", s.aperture.value);
    assembly(s.root, "/root");
  }

 private:
  template <typename PoseT>
  void pose(PoseT& p, const std::string& loc) {
    for (std::size_t k = 0; k < 3; ++k) {
      visit_(loc + "/pose/position/" + std::to_string(k), p.position[k]);
    }
    for (std::size_t k = 0; k < 3; ++k) {
      visit_(loc + "/pose/rotation_deg/" + std::to_string(k), p.rotation_deg[k]);
    }
  }

  template <typename List>
  void list(List& params, const std::string& loc) {
    for (std::size_t k = 0; k < params.size(); ++k)
      visit_(loc + "/" + std::to_string(k), params[k]);
  }

  template <typename A>
  void assembly(A& a, const std::string& loc) {
    pose(a.pose, loc);
    for (std::size_t i = 0; i < a.children.size(); ++i) {
      const std::string child = loc + "/children/" + std::to_string(i);
      auto& value = a.children[i].value;
      if (auto* sub = std::get_if<Assembly>(&value)) {
        assembly(*sub, child);
      } else {
        element(std::get<Element>(value), child);
      }
    }
  }

  template <typename E>
  void element(E& e, const std::string& loc) {
    pose(e.pose, loc);
    for (std::size_t i = 0; i < e.surfaces.size(); ++i) {
      surface(e.surfaces[i], loc + "/surfaces/" + std::to_string(i));
    }
  }

  template <typename Sf>
  void surface(Sf& s, const std::string& loc) {
    pose(s.pose, loc);
    const std::string base = loc + "/shape/base";
    if (auto* c = std::get_if<Conic>(&s.shape.base)) {
      visit_(base + "/radius", c->radius);
      visit_(base + "/conic", c->conic);
    } else if (auto* a = std::get_if<EvenAsphere>(&s.shape.base)) {
      visit_(base + "/radius", a->radius);
      visit_(base + "/conic", a->conic);
      list(a->coefficients, base + "/coefficients");
    }
    for (std::size_t t = 0; t < s.shape.terms.size(); ++t) {
      const std::string term = loc + "/shape/terms/" + std::to_string(t);
      auto& z = std::get<ZernikeSag>(s.shape.terms[t]);
      visit_(term + "/normalization_radius", z.normalization_radius);
      list(z.coefficients, term + "/coefficients");
    }
    for (std::size_t k = 0; k < s.phases.size(); ++k) {
      const std::string phase = loc + "/phases/" + std::to_string(k);
      if (auto* g = std::get_if<LinearGrating>(&s.phases[k])) {
        visit_(phase + "/lines_per_mm", g->lines_per_mm);
      } else if (auto* r = std::get_if<RadialPhase>(&s.phases[k])) {
        visit_(phase + "/normalization_radius", r->normalization_radius);
        list(r->coefficients, phase + "/coefficients");
      }
    }
  }

  const Visit& visit_;
};

/// The registered code of a parse error (expression.hpp names one of these three).
diagnostics::DiagnosticCode code_of(const detail::ExpressionError& error) {
  if (error.code == "parameters.unknown_name") return "parameters.unknown_name";
  if (error.code == "parameters.forward_reference") return "parameters.forward_reference";
  return "parameters.expression_syntax";
}

/// A diagnostic of validate() that makes the table or a binding unusable for resolve.
bool table_error(const Diagnostic& d) {
  if (d.severity != Severity::Error) return false;
  return d.code.starts_with("parameters.") || d.code.starts_with("param.") ||
         d.code.starts_with("configurations.") ||
         (d.code == "value.not_finite" && d.location.starts_with("/parameters/"));
}

}  // namespace

std::size_t configuration_count(const System& system) noexcept {
  return system.configurations.empty() ? 1 : system.configurations.size();
}

std::optional<std::size_t> find_configuration(const System& system, std::string_view name) {
  for (std::size_t k = 0; k < system.configurations.size(); ++k) {
    if (system.configurations[k].name == name) return k;
  }
  return std::nullopt;
}

std::optional<std::size_t> find_parameter(const System& system, std::string_view name) {
  for (std::size_t i = 0; i < system.parameters.size(); ++i) {
    if (system.parameters[i].name == name) return i;
  }
  return std::nullopt;
}

ParameterValues evaluate_parameters(const System& system, std::vector<Diagnostic>* diagnostics) {
  const std::size_t columns = configuration_count(system);
  ParameterValues out;
  out.rows = system.parameters.size();
  out.configurations = columns;
  out.values.assign(out.rows * columns, kUndefined);
  // A row is usable as an input if it is defined and finite in every configuration; a row that
  // uses an unusable one is undefined without a diagnostic of its own (no follow-up errors).
  std::vector<bool> usable(out.rows, false);
  std::vector<double> column(out.rows, 0.0);

  for (std::size_t i = 0; i < out.rows; ++i) {
    const ParameterRow& row = system.parameters[i];
    double* values = out.values.data() + i * columns;
    if (const auto* value = std::get_if<double>(&row.form)) {
      for (std::size_t k = 0; k < columns; ++k) values[k] = *value;
      usable[i] = std::isfinite(*value);
      continue;
    }
    if (const auto* list = std::get_if<std::vector<double>>(&row.form)) {
      // A wrong length is parameters.values_count from validate(); the row stays undefined.
      if (list->size() != columns) continue;
      usable[i] = true;
      for (std::size_t k = 0; k < columns; ++k) {
        values[k] = (*list)[k];
        usable[i] = usable[i] && std::isfinite(values[k]);
      }
      continue;
    }
    const std::string& text = std::get<ParameterExpression>(row.form).text;
    const std::string location = row_pointer(i) + "/expression";
    const auto parsed = detail::parse_expression(text, [&](std::string_view name) {
      for (std::size_t j = 0; j < out.rows; ++j) {
        if (system.parameters[j].name == name) {
          return j < i ? detail::NameLookup{detail::NameStatus::Found, j}
                       : detail::NameLookup{detail::NameStatus::Later, 0};
        }
      }
      return detail::NameLookup{detail::NameStatus::Unknown, 0};
    });
    if (const auto* error = std::get_if<detail::ExpressionError>(&parsed)) {
      add(diagnostics, code_of(*error), location, "row '" + row.name + "': " + error->message);
      continue;
    }
    const detail::Expression& expression = std::get<detail::Expression>(parsed);
    bool inputs_usable = true;
    for (const detail::Instruction& in : expression.code) {
      if (in.op == detail::OpCode::Row && !usable[in.row]) inputs_usable = false;
    }
    if (!inputs_usable) continue;  // undefined, reported at its input
    usable[i] = true;
    for (std::size_t k = 0; k < columns; ++k) {
      for (std::size_t j = 0; j < i; ++j) column[j] = out.values[j * columns + k];
      values[k] = detail::evaluate(expression, column);
      if (!std::isfinite(values[k])) {
        usable[i] = false;
        add(diagnostics, "parameters.not_finite", location,
            "row '" + row.name + "' is not finite (" + std::to_string(values[k]) + ") in " +
                configuration_text(system, k));
      }
    }
  }
  return out;
}

System resolve_parameters(const System& system, std::size_t configuration) {
  const std::size_t columns = configuration_count(system);
  if (configuration >= columns) {
    throw std::invalid_argument(
        "resolve_parameters: configuration " + std::to_string(configuration) + " does not exist (" +
        std::to_string(columns) + " configuration" + (columns == 1 ? "" : "s") + ")");
  }
  for (const Diagnostic& d : validate(system)) {
    if (table_error(d)) {
      throw std::invalid_argument("resolve_parameters: " + to_string(d));
    }
  }
  const ParameterValues values = evaluate_parameters(system);
  System out = system;
  for_each_param(out, [&](std::string_view /*pointer*/, Param& p) {
    if (!p.is_bound()) return;
    // validate() found the row (param.unknown_parameter) and a finite value in every
    // configuration (parameters.not_finite, value.not_finite of the row).
    const std::size_t row = *find_parameter(system, *p.param);
    p = Param(values.at(row, configuration));
  });
  return out;
}

void for_each_param(System& system,
                    const std::function<void(std::string_view pointer, Param& param)>& visit) {
  ParamWalker<System, std::function<void(std::string_view, Param&)>>(visit).run(system);
}

void for_each_param(
    const System& system,
    const std::function<void(std::string_view pointer, const Param& param)>& visit) {
  ParamWalker<const System, std::function<void(std::string_view, const Param&)>>(visit).run(system);
}

}  // namespace rtt::model
