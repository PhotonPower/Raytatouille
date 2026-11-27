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

/// RFC 6902 patch from the input to the external values p: one "replace" per changed variable
/// (ADR 0030, point 12), for OptimResult::patch and IterationState::patch alike (#378).
std::string patch_of(const std::vector<Variable>& variables,
                     std::span<const double> p,
                     std::span<const std::uint8_t> changed) {
  std::string patch = "[";
  for (std::size_t j = 0; j < variables.size(); ++j) {
    if (changed[j] == 0) continue;
    if (patch.size() > 1) patch += ",";
    patch += R"({"op":"replace","path":)" + json_string(variables[j].pointer) + R"(,"value":)" +
             number(p[j]) + "}";
  }
  return patch + "]";
}

/// optimize() with or without the overlay of a trial (ADR 0030, addendum #256), with or without
/// the live state (#378).
OptimResult run(const model::System& system,
                const material::MaterialLibrary& materials,
                const coating::CoatingLibrary* coatings,
                const OptimizeOptions& options,
                const trace::RunControl& control,
                const TrialOverlay* overlay,
                const IterationCallback* on_iteration) {
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
    result.diagnostics = inactive_variables(system);
    return result;
  }
  const MeritFunction merit(system, materials, coatings, overlay, options.aiming);
  const std::vector<Variable>& variables = merit.variables();
  // ADR 0030, point 10 (addendum #167): a run without variables or without a merit function is
  // an input error, never a silent "converged".
  std::vector<model::Diagnostic> errors;
  if (variables.empty()) {
    errors.push_back(diagnostic("optim.no_variables", "",
                                "no variable Param and no variable row of the parameter table"));
  }
  // ADR 0049, point 4: operands left out because their surface is inactive do not count.
  const std::vector<bool>& inactive = merit.inactive_operands();
  // Exact constraints (ADR 0056): the active operands with a constraint, in file order. An
  // inactive one is left out like any operand (its residual stays 0).
  const std::size_t n_operands = system.optimization.operands.size();
  std::vector<std::size_t> constraints;
  for (std::size_t i = 0; i < n_operands; ++i) {
    const bool constrained =
        std::visit([](const auto& o) { return o.common.constraint.has_value(); },
                   system.optimization.operands[i]);
    if (constrained && !inactive[i]) constraints.push_back(i);
  }
  const auto is_constraint = [&constraints](std::size_t i) {
    return std::find(constraints.begin(), constraints.end(), i) != constraints.end();
  };
  // Declared as a constraint in the file, active or not: such an operand is never part of F, so
  // it counts neither in the weights of phi nor in the table with its weight (ADR 0056, point 7;
  // review of #196). Only the active ones go to the solver.
  const auto declared_constraint = [&system, n_operands](std::size_t i) {
    return i < n_operands &&
           std::visit([](const auto& o) { return o.common.constraint.has_value(); },
                      system.optimization.operands[i]);
  };
  std::size_t objective_operands = 0;
  for (std::size_t i = 0; i < n_operands; ++i) {
    if (!inactive[i] && !is_constraint(i)) ++objective_operands;
  }
  if (system.optimization.generators.empty() && objective_operands == 0) {
    errors.push_back(diagnostic(
        "optim.no_operands", "",
        system.optimization.empty()
            ? "the system has no operand and no generator (section optimization)"
        : constraints.empty()
            ? "every operand is at a surface of an inactive node and no generator is left"
            : "only exact constraints and no operand or generator to minimize (ADR 0056)"));
  }
  if (!variables.empty() && constraints.size() > variables.size()) {
    errors.push_back(diagnostic("optim.constraints_too_many", "/optimization",
                                std::to_string(constraints.size()) + " exact constraints for " +
                                    std::to_string(variables.size()) +
                                    " variables; at most one per variable (ADR 0056, point 2)"));
  }
  if (!errors.empty()) {
    // Why there may be no variable: the variable Params of inactive nodes (ADR 0049, point 4).
    if (variables.empty()) {
      const std::vector<model::Diagnostic> left_out = inactive_variables(system);
      errors.insert(errors.end(), left_out.begin(), left_out.end());
    }
    throw OptimError(std::move(errors));
  }
  const std::vector<double> p0 = merit.start();
  std::vector<Bounds> bounds;
  bounds.reserve(variables.size());
  for (const Variable& v : variables) bounds.push_back(v.bounds);

  // Start (ADR 0030, point 10): the exception of an analysis comes through unchanged; an
  // undefined value is an input error naming the operand.
  const MeritEvaluation start = merit.evaluate(p0);
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

  // ADR 0030, addenda #289 and #299: the set of rays that arrive at the start holds for the run,
  // for each generator and each spot_rms and opd_rms operand. An evaluation with another set (a
  // ray of the start lost, or a ray lost at the start arrived) is invalid, so the merit is always
  // a function of p over the same rays (point 11).
  std::vector<std::vector<std::uint8_t>> start_arrived;
  start_arrived.reserve(start.generators.size());
  for (const GeneratorStats& g : start.generators) start_arrived.push_back(g.arrived);
  const std::vector<std::vector<std::uint8_t>>& start_operand_arrived = start.operand_arrived;
  // Evaluations refused for this reason; an integer sum, independent of the order of the
  // concurrent Jacobian evaluations.
  std::atomic<std::size_t> changed_sets{0};

  // Sum of the weights for phi = 2 F / weights (ADR 0030, point 12), in the history and the live
  // state alike.
  double weights = 0.0;
  for (std::size_t i = 0; i < n_operands; ++i) {
    if (declared_constraint(i)) continue;  // not part of F (ADR 0056, point 7)
    weights +=
        std::visit([](const auto& o) { return o.common.weight; }, system.optimization.operands[i]);
  }
  for (const model::Generator& g : system.optimization.generators) {
    weights += std::visit([](const auto& o) { return o.weight; }, g);
  }
  const auto phi_of = [weights](double F) { return weights > 0.0 ? 2.0 * F / weights : 0.0; };
  LmAcceptedCallback accepted;
  if (on_iteration != nullptr && *on_iteration) {
    accepted = [&variables, &phi_of, on_iteration](const LmIteration& it, std::span<const double> p,
                                                   std::span<const std::uint8_t> changed) {
      (*on_iteration)(IterationState{it.k, phi_of(it.F), std::vector<double>(p.begin(), p.end()),
                                     patch_of(variables, p, changed)});
    };
  }

  LmOptions lm_options;
  lm_options.max_iterations = options.max_iterations;
  lm_options.ftol = options.ftol;
  lm_options.xtol = options.xtol;
  lm_options.gtol = options.gtol;
  lm_options.tau = options.tau;
  lm_options.function_precision = options.function_precision;
  // The solver takes the constraints as the last residuals (ADR 0056, point 8): `order` lists the
  // merit residuals in the solver's order, the others first in their order, then the
  // constraints; each with the tolerance eps_c max(|t|, 1).
  std::vector<std::size_t> order;
  order.reserve(merit.size());
  for (std::size_t i = 0; i < merit.size(); ++i) {
    if (!is_constraint(i)) order.push_back(i);
  }
  std::vector<double> tolerances;
  tolerances.reserve(constraints.size());
  for (const std::size_t i : constraints) {
    order.push_back(i);
    const double t =
        std::visit([](const auto& o) { return o.common.target; }, system.optimization.operands[i]);
    tolerances.push_back(options.constraint_tolerance * std::max(std::abs(t), 1.0));
  }
  lm_options.constraint_tolerances = tolerances;
  const LmResult lm = levenberg_marquardt(
      [&merit, &start_arrived, &start_operand_arrived, &changed_sets, &order](
          std::span<const double> p, std::span<double> f) {
        // Exceptions and non-finite residuals are invalid evaluations of the solver.
        const MeritEvaluation e = merit.evaluate(p);
        if (e.operand_arrived != start_operand_arrived) {
          changed_sets.fetch_add(1, std::memory_order_relaxed);
          return false;
        }
        for (std::size_t g = 0; g < e.generators.size(); ++g) {
          if (e.generators[g].arrived != start_arrived[g]) {
            changed_sets.fetch_add(1, std::memory_order_relaxed);
            return false;
          }
        }
        for (std::size_t k = 0; k < order.size(); ++k) f[k] = e.residuals[order[k]];
        return true;
      },
      merit.size(), p0, bounds, lm_options, control, accepted);

  // ADR 0056, point 5: a start that cannot be projected onto the constraints is an input error at
  // the constraint with the largest relative violation.
  if (lm.status == LmStatus::Failed && lm.infeasible_constraint) {
    const std::size_t i = constraints[*lm.infeasible_constraint];
    throw OptimError({diagnostic(
        "optim.constraint_infeasible", "/optimization/operands/" + std::to_string(i),
        "the start cannot be projected onto the exact constraints; this one has the largest "
        "relative violation (value " +
            number(start.values[i]) + ")")});
  }

  OptimResult result;
  result.status = lm.status;
  result.iterations = lm.iterations;
  result.diagnostics = inactive_variables(system);  // ADR 0049, point 4: first
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

  double squares = 0.0;  // of F: without the constraints
  for (std::size_t i = 0; i < end.residuals.size(); ++i) {
    if (!is_constraint(i)) squares += end.residuals[i] * end.residuals[i];
  }
  for (const LmIteration& it : lm.history) {
    result.history.push_back(
        {it.k, phi_of(it.F), it.mu, it.rho, it.step_norm, it.accepted, it.evaluations});
  }
  for (std::size_t i = 0; i < end.values.size(); ++i) {
    const model::Operand& op = system.optimization.operands[i];
    const double target = std::visit([](const auto& o) { return o.common.target; }, op);
    const bool constraint = declared_constraint(i);
    const double weight =
        constraint ? 0.0 : std::visit([](const auto& o) { return o.common.weight; }, op);
    const double r = end.residuals[i];
    result.operands.push_back({"/optimization/operands/" + std::to_string(i), end.values[i], target,
                               weight, !constraint && squares > 0.0 ? 100.0 * r * r / squares : 0.0,
                               constraint});
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

  for (std::size_t j = 0; j < variables.size(); ++j) {
    const Variable& v = variables[j];
    result.variables.push_back({v.pointer, v.row, v.configuration, v.start, lm.p[j],
                                lm.changed[j] != 0, lm.at_bound[j] != 0});
    if (lm.at_bound[j] != 0) {
      result.diagnostics.push_back(diagnostic("optim.parameter_at_bound", v.pointer,
                                              "the result lies at a bound of the variable"));
    }
  }
  result.patch = patch_of(variables, lm.p, lm.changed);

  if (lm.failed_evaluations > 0) {
    const std::size_t changed = changed_sets.load();
    result.diagnostics.push_back(diagnostic(
        "optim.evaluation_failed", "/optimization",
        std::to_string(lm.failed_evaluations) +
            " evaluation(s) failed (analysis error, undefined value or changed ray set); the "
            "trials counted as rejected" +
            (changed > 0 ? "; " + std::to_string(changed) +
                               " of them because the set of arrived rays of a generator or a "
                               "spot_rms/opd_rms operand differed from the start (ADR 0030, "
                               "addenda #289 and #299)"
                         : std::string())));
  }
  if (lm.status == LmStatus::Failed && lm.dependent_constraints) {
    result.diagnostics.push_back(diagnostic(
        "optim.constraints_dependent", "/optimization",
        "the exact constraints are linearly dependent at the last accepted state (ADR 0056, "
        "point 2)"));
  }
  // ADR 0056, point 7: a constraint that misses its tolerance in the final state. Every accepted
  // state and the projected start are feasible, so this happens only after status Failed before
  // the start was projected (a failed Jacobian column at the start).
  for (std::size_t k = 0; k < constraints.size() && !end.values.empty(); ++k) {
    const std::size_t i = constraints[k];
    const double target =
        std::visit([](const auto& o) { return o.common.target; }, system.optimization.operands[i]);
    if (!(std::abs(end.values[i] - target) <= tolerances[k])) {
      result.diagnostics.push_back(
          diagnostic("optim.constraint_violated", "/optimization/operands/" + std::to_string(i),
                     "the exact constraint is not met in the final state: value " +
                         number(end.values[i]) + ", target " + number(target)));
    }
  }
  if (lm.status == LmStatus::Failed && lm.failed_variable) {
    result.diagnostics.push_back(
        diagnostic("optim.jacobian_failed", variables[*lm.failed_variable].pointer,
                   "no valid evaluation for the Jacobian column of this variable"));
  }
  // F7: lost rays of a spot_rms or opd_rms operand and of a generator in the final state
  // (ADR 0030, point 4, addenda #168, #289 and #299): the same rays as at the start, since every
  // accepted state has the start set. None after a cancellation (no final evaluation).
  for (std::size_t i = 0; i < end.operand_arrived.size(); ++i) {
    const std::vector<std::uint8_t>& mask = end.operand_arrived[i];
    const auto lost =
        static_cast<std::size_t>(std::count(mask.begin(), mask.end(), std::uint8_t{0}));
    if (lost == 0) continue;
    result.diagnostics.push_back(
        diagnostic("optim.rays_lost", "/optimization/operands/" + std::to_string(i),
                   std::to_string(lost) + " of " + std::to_string(mask.size()) +
                       " rays of the operand were lost at the start and stay out of the run"));
  }
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
                     const trace::RunControl& control,
                     const IterationCallback& on_iteration) {
  return run(system, materials, coatings, options, control, nullptr, &on_iteration);
}

OptimResult optimize(const model::System& system,
                     const material::MaterialLibrary& materials,
                     const coating::CoatingLibrary* coatings,
                     const OptimizeOptions& options,
                     const trace::RunControl& control,
                     const TrialOverlay& overlay) {
  return run(system, materials, coatings, options, control, &overlay, nullptr);
}

}  // namespace rtt::optim
