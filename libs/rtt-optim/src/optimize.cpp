#include "rtt/optim/optimize.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <cstdint>
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
#include "rtt/optim/generators.hpp"
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

/// optimize() with or without the overlay of a trial (ADR 0030, addendum #256).
OptimResult run(const model::System& system,
                const material::MaterialLibrary& materials,
                const coating::CoatingLibrary* coatings,
                const OptimizeOptions& options,
                const trace::RunControl& control,
                const TrialOverlay* overlay) {
  // ADR 0030, point 11 (addendum #167): a cancellation requested before the start ends the run
  // before any compile or evaluation, with the input.
  if (control.cancel && control.cancel->cancelled()) {
    OptimResult result;
    result.status = LmStatus::Cancelled;
    result.system = system;
    result.patch = "[]";
    for (const Variable& v : collect_variables(system)) {
      result.variables.push_back({v.pointer, v.row, v.configuration, v.start, v.start, false,
                                  at_bound(v.start, v.bounds)});
    }
    return result;
  }
  const MeritFunction merit(system, materials, coatings, overlay);
  const std::vector<Variable>& variables = merit.variables();
  // ADR 0030, point 10 (addendum #167): a run without variables or without a merit function is
  // an input error, never a silent "converged".
  std::vector<model::Diagnostic> errors;
  if (variables.empty()) {
    errors.push_back(diagnostic("optim.no_variables", "",
                                "no variable Param and no variable row of the parameter table"));
  }
  if (system.optimization.empty()) {
    errors.push_back(
        diagnostic("optim.no_operands", "",
                   "the system has no operand and no generator (section optimization)"));
  }
  if (!errors.empty()) throw OptimError(std::move(errors));
  const std::vector<double> p0 = merit.start();
  std::vector<Bounds> bounds;
  bounds.reserve(variables.size());
  for (const Variable& v : variables) bounds.push_back(v.bounds);

  // Start (ADR 0030, point 10): the exception of an analysis comes through unchanged; an
  // undefined value is an input error naming the operand.
  const MeritEvaluation start = merit.evaluate(p0);
  const std::size_t n_operands = system.optimization.operands.size();
  for (std::size_t i = 0; i < n_operands; ++i) {
    if (!std::isfinite(start.residuals[i])) {
      const std::string& why = start.undefined[i];
      throw std::invalid_argument("optimize: operand /optimization/operands/" + std::to_string(i) +
                                  " has no finite residual at the start" +
                                  (why.empty() ? std::string() : ": " + why));
    }
  }
  // Generator residuals follow in blocks of generator_sizes() (ADR 0030, addendum #168).
  for (std::size_t g = 0, offset = n_operands; g < merit.generator_sizes().size(); ++g) {
    const std::size_t block_end = offset + merit.generator_sizes()[g];
    for (std::size_t i = offset; i < block_end; ++i) {
      if (std::isfinite(start.residuals[i])) continue;
      const std::string& why = start.generators[g].undefined;
      throw std::invalid_argument("optimize: generator /optimization/generators/" +
                                  std::to_string(g) + " has no finite residual at the start" +
                                  (why.empty() ? std::string() : ": " + why));
    }
    offset = block_end;
  }

  // ADR 0030, addendum #289: the set of generator rays that arrive at the start holds for the
  // run. An evaluation with another set (a ray of the start lost, or a ray lost at the start
  // arrived) is invalid, so the merit is always a function of p over the same rays (point 11).
  std::vector<std::vector<std::uint8_t>> start_arrived;
  start_arrived.reserve(start.generators.size());
  for (const GeneratorStats& g : start.generators) start_arrived.push_back(g.arrived);
  // Evaluations refused for this reason; an integer sum, independent of the order of the
  // concurrent Jacobian evaluations.
  std::atomic<std::size_t> changed_sets{0};

  const LmOptions lm_options{options.max_iterations, options.ftol, options.xtol,
                             options.gtol,           options.tau,  options.function_precision};
  const LmResult lm = levenberg_marquardt(
      [&merit, &start_arrived, &changed_sets](std::span<const double> p, std::span<double> f) {
        // Exceptions and non-finite residuals are invalid evaluations of the solver.
        const MeritEvaluation e = merit.evaluate(p);
        for (std::size_t g = 0; g < e.generators.size(); ++g) {
          if (e.generators[g].arrived != start_arrived[g]) {
            changed_sets.fetch_add(1, std::memory_order_relaxed);
            return false;
          }
        }
        std::copy(e.residuals.begin(), e.residuals.end(), f.begin());
        return true;
      },
      merit.size(), p0, bounds, lm_options, control);

  OptimResult result;
  result.status = lm.status;
  result.iterations = lm.iterations;
  result.failed_evaluations = lm.failed_evaluations;
  result.system = with_values(system, variables, lm.p);

  // The final state: operand table and warnings (the last accepted state was valid). Not after
  // a cancellation: no evaluation after the request (ADR 0030, point 11, addendum #167).
  MeritEvaluation end;
  result.evaluations = lm.evaluations + 1;  // the start check
  if (lm.status != LmStatus::Cancelled) {
    end = merit.evaluate(lm.p);
    ++result.evaluations;
  }

  double weights = 0.0;
  double squares = 0.0;
  for (const model::Operand& op : system.optimization.operands) {
    weights += std::visit([](const auto& o) { return o.common.weight; }, op);
  }
  for (const model::Generator& g : system.optimization.generators) {
    weights += std::visit([](const auto& o) { return o.weight; }, g);
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
  // Generators (ADR 0030, addendum #168): only with an evaluated final state.
  for (std::size_t g = 0, offset = end.values.size(); g < end.generators.size(); ++g) {
    const std::size_t n = merit.generator_sizes()[g];
    double own = 0.0;
    for (std::size_t i = offset; i < offset + n; ++i) own += end.residuals[i] * end.residuals[i];
    offset += n;
    const GeneratorStats& st = end.generators[g];
    result.generators.push_back(
        {"/optimization/generators/" + std::to_string(g), std::sqrt(st.mean_square),
         std::visit([](const auto& o) { return o.weight; }, system.optimization.generators[g]),
         squares > 0.0 ? 100.0 * own / squares : 0.0, st.rays_launched, st.rays_lost});
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
    const std::size_t changed = changed_sets.load();
    result.diagnostics.push_back(diagnostic(
        "optim.evaluation_failed", "/optimization",
        std::to_string(lm.failed_evaluations) +
            " evaluation(s) failed (analysis error, undefined value or changed generator rays); "
            "the trials counted as rejected" +
            (changed > 0 ? "; " + std::to_string(changed) +
                               " of them because the set of arrived generator rays differed "
                               "from the start (ADR 0030, addendum #289)"
                         : std::string())));
  }
  if (lm.status == LmStatus::Failed && lm.failed_variable) {
    result.diagnostics.push_back(
        diagnostic("optim.jacobian_failed", variables[*lm.failed_variable].pointer,
                   "no valid evaluation for the Jacobian column of this variable"));
  }
  // F7: lost rays of a generator in the final state (ADR 0030, point 4, addenda #168 and #289):
  // the same rays as at the start, since every accepted state has the start set.
  for (const GeneratorValue& g : result.generators) {
    if (g.rays_lost == 0) continue;
    result.diagnostics.push_back(diagnostic(
        "optim.rays_lost", g.pointer,
        std::to_string(g.rays_lost) + " of " + std::to_string(g.rays_launched) +
            " rays of the generator were lost at the start and stay out of the run (residuals "
            "0)"));
  }
  for (const model::Diagnostic& d : end.warnings) result.diagnostics.push_back(d);
  return result;
}

}  // namespace

OptimResult optimize(const model::System& system,
                     const material::MaterialLibrary& materials,
                     const coating::CoatingLibrary* coatings,
                     const OptimizeOptions& options,
                     const trace::RunControl& control) {
  return run(system, materials, coatings, options, control, nullptr);
}

OptimResult optimize(const model::System& system,
                     const material::MaterialLibrary& materials,
                     const coating::CoatingLibrary* coatings,
                     const OptimizeOptions& options,
                     const trace::RunControl& control,
                     const TrialOverlay& overlay) {
  return run(system, materials, coatings, options, control, &overlay);
}

}  // namespace rtt::optim
