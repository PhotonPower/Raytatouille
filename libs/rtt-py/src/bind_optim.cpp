// Optimization from Python (#169, ADR 0030; rtt/optim/*.hpp): the variables of a system, its
// merit function and optimize() with the result classes. raytatouille.optim wraps the
// functions (default libraries, warnings).

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "bindings.hpp"
#include "rtt/coating/catalog.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/system.hpp"
#include "rtt/optim/generators.hpp"
#include "rtt/optim/levenberg_marquardt.hpp"
#include "rtt/optim/merit.hpp"
#include "rtt/optim/optimize.hpp"
#include "rtt/optim/variables.hpp"
#include "rtt/trace/run_control.hpp"

namespace nb = nanobind;
using namespace nb::literals;

namespace rtt::py {
namespace {

/// The history of a run, bound with columns().
struct OptimIterations {
  std::vector<optim::OptimIteration> points;
};

/// Read-only NumPy bool copy of the accepted flags (read_only_array needs vector::data(), which
/// std::vector<bool> lacks).
ReadOnlyArray<bool> accepted_flags(const OptimIterations& h) {
  const std::size_t n = h.points.size();
  auto data = std::make_unique<bool[]>(n);
  for (std::size_t i = 0; i < n; ++i) data[i] = h.points[i].accepted;
  const std::size_t shape[1] = {n};
  const bool* values = data.get();
  // The capsule owns the copy from here on and frees it with the array.
  const nb::capsule owner(data.get(), [](void* p) noexcept {
    const std::unique_ptr<bool[]> owned(static_cast<bool*>(p));
  });
  static_cast<void>(data.release());
  return ReadOnlyArray<bool>(values, 1, shape, owner);
}

/// Read-only NumPy copy of a vector of doubles.
ReadOnlyArray<double> doubles(const std::vector<double>& v) {
  return read_only_array<double>(v, [](double x) { return x; });
}

/// A merit function for Python: owns a copy of the system (optim::MeritFunction keeps a
/// reference to it) and, without a library from Python, the default libraries. The libraries
/// from Python are kept alive by the binding (keep_alive).
class Merit {
 public:
  Merit(const model::System& system,
        const material::MaterialLibrary* materials,
        const coating::CoatingLibrary* coatings)
      : system_(std::make_unique<model::System>(system)),
        own_materials_(materials != nullptr ? nullptr
                                            : std::make_unique<material::MaterialLibrary>()),
        merit_(std::make_unique<optim::MeritFunction>(
            *system_, materials != nullptr ? *materials : *own_materials_, coatings)) {}

  [[nodiscard]] const optim::MeritFunction& get() const noexcept { return *merit_; }

 private:
  std::unique_ptr<model::System> system_;
  std::unique_ptr<material::MaterialLibrary> own_materials_;
  std::unique_ptr<optim::MeritFunction> merit_;
};

}  // namespace

void bind_optim(nb::module_& m) {
  using namespace optim;

  nb::enum_<LmStatus>(m, "OptimStatus", "How an optimization run ended (ADR 0030, point 8).")
      .value("CONVERGED_GRADIENT", LmStatus::ConvergedGradient,
             "scaled gradient below gtol (off by default)")
      .value("CONVERGED_STEP", LmStatus::ConvergedStep, "step below xtol; the step is not applied")
      .value("CONVERGED_MERIT", LmStatus::ConvergedMerit,
             "small relative actual and predicted decrease of the merit (ftol)")
      .value("MAX_ITERATIONS", LmStatus::MaxIterations, "max_iterations solves done")
      .value("CANCELLED", LmStatus::Cancelled,
             "cancelled through the CancelToken; the last accepted state")
      .value("FAILED", LmStatus::Failed,
             "a Jacobian column without valid evaluations; the last accepted state");

  nb::class_<OptimizeOptions>(
      m, "OptimizeOptions",
      "Settings of the solver (ADR 0030, point 2: run parameters, not part of the file).")
      .def(
          "__init__",
          [](OptimizeOptions* o, int max_iterations, double ftol, double xtol, double gtol,
             double tau, double function_precision) {
            new (o) OptimizeOptions{max_iterations, ftol, xtol, gtol, tau, function_precision};
          },
          nb::kw_only(), "max_iterations"_a = OptimizeOptions{}.max_iterations,
          "ftol"_a = OptimizeOptions{}.ftol, "xtol"_a = OptimizeOptions{}.xtol,
          "gtol"_a = OptimizeOptions{}.gtol, "tau"_a = OptimizeOptions{}.tau,
          "function_precision"_a = OptimizeOptions{}.function_precision)
      .def_rw("max_iterations", &OptimizeOptions::max_iterations,
              "k_max, the number of solves (ADR 0030, point 8).")
      .def_rw("ftol", &OptimizeOptions::ftol, "Merit test, default sqrt(eps_M).")
      .def_rw("xtol", &OptimizeOptions::xtol, "Step test, default sqrt(eps_M).")
      .def_rw("gtol", &OptimizeOptions::gtol, "Scaled gradient test; 0 switches it off.")
      .def_rw("tau", &OptimizeOptions::tau, "mu_0 = tau max diag(J^T J) (point 7).")
      .def_rw("function_precision", &OptimizeOptions::function_precision,
              "Relative precision eps_f of a merit evaluation for the difference step "
              "(point 9); default MERIT_PRECISION.");
  m.attr("MERIT_PRECISION") = kMeritPrecision;

  nb::class_<Variable>(m, "Variable", "One unknown of the optimizer (ADR 0030, point 5).")
      .def_ro("pointer", &Variable::pointer,
              "JSON pointer of the value in the edit form that the result patch replaces.")
      .def_ro("row", &Variable::row, "Name of the table row; empty for a model Param.")
      .def_ro("configuration", &Variable::configuration,
              "Column of a `values` row; None otherwise.")
      .def_ro("param_index", &Variable::param_index,
              "Position of the Param in the order of the model's Params; None for a row.")
      .def_ro("start", &Variable::start, "Value in the input system, unit of the field or row.")
      .def_prop_ro(
          "min", [](const Variable& v) { return v.bounds.min; }, "Lower bound; None if not set.")
      .def_prop_ro(
          "max", [](const Variable& v) { return v.bounds.max; }, "Upper bound; None if not set.");

  nb::class_<GeneratorStats>(m, "GeneratorStats",
                             "Rays and mean square of one generator in one evaluation.")
      .def_ro("rays_launched", &GeneratorStats::rays_launched, "Rays of the generator.")
      .def_ro("rays_lost", &GeneratorStats::rays_lost,
              "Rays with residuals 0 because they did not arrive.")
      .def_ro("mean_square", &GeneratorStats::mean_square,
              "Weighted mean square of the deviations (mm^2 for rms_spot, waves^2 for "
              "rms_wavefront), independent of the weight; NaN if a residual is NaN.")
      .def_ro("undefined", &GeneratorStats::undefined,
              "Why the residuals are NaN (empty if they are finite).");

  nb::class_<MeritEvaluation>(m, "MeritEvaluation",
                              "One evaluation of the merit function (ADR 0030, points 1, 10).")
      .def_prop_ro(
          "values", [](const MeritEvaluation& e) { return doubles(e.values); },
          nb::rv_policy::reference,
          "Value of each operand in file order, unit of the operand; NaN if not defined "
          "(read-only copy).")
      .def_prop_ro(
          "residuals", [](const MeritEvaluation& e) { return doubles(e.residuals); },
          nb::rv_policy::reference,
          "sqrt(weight) (value - target) per operand, then the residuals of the generators in "
          "file order (read-only copy).")
      .def_ro("undefined", &MeritEvaluation::undefined,
              "Why an operand has no value (empty string if it has one), operand order.")
      .def_ro("warnings", &MeritEvaluation::warnings,
              "Warnings of the analyses (rays.lost, stop.clips_beam), each once.")
      .def_ro("generators", &MeritEvaluation::generators, "One GeneratorStats per generator.")
      .def("valid", &MeritEvaluation::valid, "True if every residual is finite.");

  nb::class_<Merit>(m, "MeritFunction",
                    "The merit function of a system (ADR 0030): its variables and its operands and "
                    "generators. Keeps a copy of the system.")
      .def(nb::init<const model::System&, const material::MaterialLibrary*,
                    const coating::CoatingLibrary*>(),
           "system"_a, "materials"_a.none() = nb::none(), "coatings"_a.none() = nb::none(),
           nb::keep_alive<1, 3>(), nb::keep_alive<1, 4>(),
           "Prepares the evaluation and checks the operands (ADR 0030, point 10). Without "
           "`materials` only VACUUM, AIR and CONST: resolve.")
      .def_prop_ro(
          "variables", [](const Merit& f) { return f.get().variables(); },
          "The variables, in the order of ADR 0030, point 5.")
      .def(
          "start", [](const Merit& f) { return doubles(f.get().start()); },
          nb::rv_policy::reference, "Start values of the variables (read-only copy).")
      .def_prop_ro(
          "size", [](const Merit& f) { return f.get().size(); },
          "Number of residuals: one per operand, then generator_sizes.")
      .def_prop_ro(
          "generator_sizes", [](const Merit& f) { return f.get().generator_sizes(); },
          "Number of residuals of each generator, file order.")
      .def(
          "evaluate",
          [](const Merit& f, const std::vector<double>& values, std::optional<int> threads) {
            if (values.size() != f.get().variables().size()) {
              throw std::invalid_argument(
                  "MeritFunction.evaluate: " + std::to_string(values.size()) + " values for " +
                  std::to_string(f.get().variables().size()) + " variables");
            }
            return released(threads, [&] { return f.get().evaluate(values); });
          },
          "values"_a, "threads"_a.none() = nb::none(),
          "Evaluates the operands and generators with the variables at `values` (one per "
          "variable, unit of the field or row).");

  m.def(
      "collect_variables", [](const model::System& system) { return collect_variables(system); },
      "system"_a, "The variables of `system` in the order of ADR 0030, point 5.");

  auto history = columns<OptimIterations>(
      m, "OptimIterations",
      "One entry per solve of the run (ADR 0030, point 12); read-only NumPy copies.");
  column<std::int32_t>(
      history, "k", [](const OptimIteration& i) { return i.k; },
      "Number of the solve, 1, 2, ... (copy).");
  column<double>(
      history, "phi", [](const OptimIteration& i) { return i.phi; },
      "Normalised merit 2 F / sum of the weights of the current state after the solve "
      "(copy).");
  column<double>(history, "mu", [](const OptimIteration& i) { return i.mu; }, "Damping (copy).");
  column<double>(
      history, "rho", [](const OptimIteration& i) { return i.rho; },
      "Gain ratio; NaN without a valid trial (copy).");
  column<double>(
      history, "step_norm", [](const OptimIteration& i) { return i.step_norm; },
      "||D h|| in scaled internal variables; NaN if no step (copy).");
  history.def_prop_ro("accepted", &accepted_flags, nb::rv_policy::reference,
                      "True if the trial was accepted (copy).");
  column<std::uint64_t>(
      history, "evaluations", [](const OptimIteration& i) { return i.evaluations; },
      "Evaluations of the solver so far, cumulative (copy).");

  nb::class_<OperandValue>(m, "OperandValue", "One operand in the final state.")
      .def_ro("pointer", &OperandValue::pointer, "/optimization/operands/i")
      .def_ro("value", &OperandValue::value, "Value, unit of the operand.")
      .def_ro("target", &OperandValue::target, "Target, unit of the operand.")
      .def_ro("weight", &OperandValue::weight, "Weight >= 0.")
      .def_ro("contribution", &OperandValue::contribution, "Share of the merit, percent.");

  nb::class_<GeneratorValue>(m, "GeneratorValue", "One generator in the final state.")
      .def_ro("pointer", &GeneratorValue::pointer, "/optimization/generators/g")
      .def_ro("rms", &GeneratorValue::rms,
              "Weighted RMS: mm (rms_spot) or waves at the reference wavelength "
              "(rms_wavefront); also for weight 0.")
      .def_ro("weight", &GeneratorValue::weight, "Weight >= 0.")
      .def_ro("contribution", &GeneratorValue::contribution, "Share of the merit, percent.")
      .def_ro("rays_launched", &GeneratorValue::rays_launched, "Rays of the generator.")
      .def_ro("rays_lost", &GeneratorValue::rays_lost, "Rays that did not arrive.");

  nb::class_<VariableValue>(m, "VariableValue", "One variable in the final state.")
      .def_ro("pointer", &VariableValue::pointer, "JSON pointer of the value in the edit form.")
      .def_ro("row", &VariableValue::row, "Name of the table row; empty for a model Param.")
      .def_ro("configuration", &VariableValue::configuration,
              "Column of a `values` row; None otherwise.")
      .def_ro("start", &VariableValue::start, "Input value.")
      .def_ro("end", &VariableValue::end, "Result value (bitwise the input if unchanged).")
      .def_ro("changed", &VariableValue::changed, "The run changed this variable.")
      .def_ro("at_bound", &VariableValue::at_bound, "The result lies at a bound.");

  nb::class_<OptimResult>(
      m, "OptimResult",
      "Result of optimize() (ADR 0030, point 12). `system` is a new System; the input is not "
      "changed. `patch` turns the input into `system` (RFC 6902, for an Editor with undo).")
      .def_ro("status", &OptimResult::status, "How the run ended.")
      .def_ro("system", &OptimResult::system,
              "The input with the result values of the changed variables.")
      .def_ro("patch", &OptimResult::patch,
              "RFC 6902 patch from the input to `system` (JSON text): one replace per changed "
              "variable; \"[]\" if nothing changed.")
      .def_prop_ro(
          "history", [](const OptimResult& r) { return OptimIterations{r.history}; },
          "One entry per solve.")
      .def_ro("operands", &OptimResult::operands, "Final state of the operands, file order.")
      .def_ro("generators", &OptimResult::generators, "Final state of the generators, file order.")
      .def_ro("variables", &OptimResult::variables,
              "Final state of the variables, order of ADR 0030, point 5.")
      .def_ro("diagnostics", &OptimResult::diagnostics,
              "optim.parameter_at_bound, optim.evaluation_failed, optim.jacobian_failed, "
              "optim.rays_lost, then the warnings of the analyses in the final state.")
      .def_ro("iterations", &OptimResult::iterations, "Number of solves.")
      .def_ro("evaluations", &OptimResult::evaluations,
              "Evaluations of the merit function in this call.")
      .def_ro("failed_evaluations", &OptimResult::failed_evaluations,
              "Trials that were invalid (ADR 0030, point 10).");

  m.def(
      "optimize",
      [](const model::System& system, const material::MaterialLibrary* materials,
         const coating::CoatingLibrary* coatings, const OptimizeOptions& options,
         std::optional<int> threads, const std::optional<trace::CancelToken>& cancel,
         const std::optional<nb::callable>& progress) {
        const trace::RunControl control = run_control(cancel, progress);
        const material::MaterialLibrary default_materials;
        const material::MaterialLibrary& library =
            materials != nullptr ? *materials : default_materials;
        return released(threads,
                        [&] { return optimize(system, library, coatings, options, control); });
      },
      "system"_a, "materials"_a.none(), "coatings"_a.none(), "options"_a, "threads"_a.none(),
      "cancel"_a.none(), "progress"_a.none(),
      "Optimizes the variables of `system` against its merit function (ADR 0030).");
}

}  // namespace rtt::py
