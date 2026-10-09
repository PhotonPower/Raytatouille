#include "rtt/optim/optimize.hpp"

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <string>
#include <system_error>
#include <utility>
#include <variant>
#include <vector>

#include "rtt/diagnostics/codes.hpp"
#include "rtt/model/optimization.hpp"
#include "rtt/optim/bounds.hpp"
#include "rtt/optim/merit.hpp"
#include "rtt/optim/variables.hpp"

namespace rtt::optim {

namespace {

model::Diagnostic diagnostic(diagnostics::DiagnosticCode code,
                             std::string location,
                             std::string message) {
  return {code.severity(), std::move(location), std::move(message), std::string(code.str())};
}

/// Shortest round-trip form of a finite double (std::to_chars, ADR 0030, point 13).
std::string number(double v) {
  std::array<char, 32> buffer{};
  const std::to_chars_result r = std::to_chars(buffer.data(), buffer.data() + buffer.size(), v);
  if (r.ec != std::errc()) throw std::logic_error("optimize: number does not fit");
  return {buffer.data(), r.ptr};
}

/// A JSON string; the pointers come from for_each_param and table indices, the escape covers
/// the two characters a JSON string cannot hold literally there.
std::string json_string(const std::string& s) {
  std::string out = "\"";
  for (const char c : s) {
    if (c == '"' || c == '\\') out += '\\';
    out += c;
  }
  return out + "\"";
}

}  // namespace

OptimResult optimize(const model::System& system,
                     const material::MaterialLibrary& materials,
                     const coating::CoatingLibrary* coatings,
                     const OptimizeOptions& options,
                     const trace::RunControl& control) {
  const MeritFunction merit(system, materials, coatings);
  const std::vector<Variable>& variables = merit.variables();
  if (variables.empty()) {
    throw OptimError({diagnostic("optim.no_variables", "",
                                 "no variable Param and no variable row of the parameter table")});
  }
  const std::vector<double> p0 = merit.start();
  std::vector<Bounds> bounds;
  bounds.reserve(variables.size());
  for (const Variable& v : variables) bounds.push_back(v.bounds);

  // Start (ADR 0030, point 10): the exception of an analysis comes through unchanged; an
  // undefined value is an input error naming the operand.
  const MeritEvaluation start = merit.evaluate(p0);
  for (std::size_t i = 0; i < start.residuals.size(); ++i) {
    if (!std::isfinite(start.residuals[i])) {
      const std::string& why = start.undefined[i];
      throw std::invalid_argument("optimize: operand /optimization/operands/" + std::to_string(i) +
                                  " has no finite residual at the start" +
                                  (why.empty() ? std::string() : ": " + why));
    }
  }

  LmResult lm;
  if (merit.size() == 0) {
    // No operand: F = 0 for every p, so g = 0 and the gradient test holds at the start
    // (ADR 0030, point 8, step 1); the solver needs m >= 1.
    lm.status = LmStatus::ConvergedGradient;
    lm.p = p0;
    lm.changed.assign(p0.size(), 0);
    lm.at_bound.resize(p0.size());
    for (std::size_t j = 0; j < p0.size(); ++j) {
      lm.at_bound[j] = at_bound(p0[j], bounds[j]) ? 1 : 0;
    }
    lm.F = 0.0;
  } else {
    const LmOptions lm_options{options.max_iterations, options.ftol, options.xtol,
                               options.gtol,           options.tau,  options.function_precision};
    lm = levenberg_marquardt(
        [&merit](std::span<const double> p, std::span<double> f) {
          // Exceptions and non-finite residuals are invalid evaluations of the solver.
          const MeritEvaluation e = merit.evaluate(p);
          std::copy(e.residuals.begin(), e.residuals.end(), f.begin());
          return true;
        },
        merit.size(), p0, bounds, lm_options, control);
  }

  OptimResult result;
  result.status = lm.status;
  result.iterations = lm.iterations;
  result.failed_evaluations = lm.failed_evaluations;
  result.system = with_values(system, variables, lm.p);

  // The final state: operand table and warnings (the last accepted state was valid).
  const MeritEvaluation end = merit.evaluate(lm.p);
  result.evaluations = lm.evaluations + 2;  // the start check and the final state

  double weights = 0.0;
  double squares = 0.0;
  for (const model::Operand& op : system.optimization.operands) {
    weights += std::visit([](const auto& o) { return o.common.weight; }, op);
  }
  for (const double r : end.residuals) squares += r * r;
  for (const LmIteration& it : lm.history) {
    result.history.push_back({it.k, weights > 0.0 ? 2.0 * it.F / weights : 0.0, it.mu, it.rho,
                              it.step_norm, it.accepted, it.evaluations});
  }
  for (std::size_t i = 0; i < end.values.size(); ++i) {
    const model::Operand& op = system.optimization.operands[i];
    const double target = std::visit([](const auto& o) { return o.common.target; }, op);
    const double weight = std::visit([](const auto& o) { return o.common.weight; }, op);
    const double r = end.residuals[i];
    result.operands.push_back({"/optimization/operands/" + std::to_string(i), end.values[i], target,
                               weight, squares > 0.0 ? 100.0 * r * r / squares : 0.0});
  }

  std::string patch = "[";
  for (std::size_t j = 0; j < variables.size(); ++j) {
    const Variable& v = variables[j];
    const bool changed = lm.changed[j] != 0;
    result.variables.push_back(
        {v.pointer, v.row, v.configuration, v.start, lm.p[j], changed, lm.at_bound[j] != 0});
    if (changed) {
      if (patch.size() > 1) patch += ",";
      patch += R"({"op":"replace","path":)" + json_string(v.pointer) + R"(,"value":)" +
               number(lm.p[j]) + "}";
    }
    if (lm.at_bound[j] != 0) {
      result.diagnostics.push_back(diagnostic("optim.parameter_at_bound", v.pointer,
                                              "the result lies at a bound of the variable"));
    }
  }
  result.patch = patch + "]";

  if (lm.failed_evaluations > 0) {
    result.diagnostics.push_back(
        diagnostic("optim.evaluation_failed", "/optimization",
                   std::to_string(lm.failed_evaluations) +
                       " evaluation(s) failed (analysis error or undefined value); the trials "
                       "counted as rejected"));
  }
  if (lm.status == LmStatus::Failed && lm.failed_variable) {
    result.diagnostics.push_back(
        diagnostic("optim.jacobian_failed", variables[*lm.failed_variable].pointer,
                   "no valid evaluation for the Jacobian column of this variable"));
  }
  for (const model::Diagnostic& d : end.warnings) result.diagnostics.push_back(d);
  return result;
}

}  // namespace rtt::optim
