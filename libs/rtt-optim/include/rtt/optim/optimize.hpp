#pragma once

/// @file optimize.hpp
/// Optimization of a system against its merit function (ADR 0030, points 10-12; #167):
/// Levenberg-Marquardt (levenberg_marquardt.hpp) over the variables of the system
/// (variables.hpp) with the residuals of its operands (merit.hpp).

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

#include "rtt/coating/catalog.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/system.hpp"
#include "rtt/model/validate.hpp"
#include "rtt/optim/levenberg_marquardt.hpp"
#include "rtt/trace/run_control.hpp"

namespace rtt::optim {

/// Relative precision eps_f of a merit evaluation (ADR 0030, point 9 and addendum #167), the
/// default of OptimizeOptions::function_precision; a run option, so a caller with a smoother
/// problem may set it smaller. Measured in #167 (test_merit_noise.cpp, tag [.noise], GCC 16
/// Release) over 26 paths of the reference systems, the M5 acceptance systems included: the
/// largest relative noise is 8.17e-11 (m2/telecentric_4f, opd_rms on axis: rounding noise of a
/// perfectly corrected afocal relay, whose ||f(theta_0)||_inf is small without an EFL), the
/// second 3.11e-12 (m1/telecentric_singlet), all others <= 1.5e-13. Times 10 and rounded up to a
/// power of ten: 1e-9, so the difference step is h = 1e-3 max(|theta|, 1). Pairs with a step or
/// a kink (telecentric aiming threshold, an RMS at its minimum) are listed, not counted.
inline constexpr double kMeritPrecision = 1e-9;

/// Settings of the solver (ADR 0030, point 2: run parameters, not part of the file); as
/// LmOptions, with the measured precision of a merit evaluation as default.
struct OptimizeOptions {
  int max_iterations = 100;                     ///< k_max, number of solves (point 8)
  double ftol = 1.4901161193847656e-8;          ///< merit test, sqrt(eps_M) (point 8)
  double xtol = 1.4901161193847656e-8;          ///< step test, sqrt(eps_M) (point 8)
  double gtol = 0.0;                            ///< scaled gradient test, off (point 8)
  double tau = 1e-3;                            ///< mu_0 = tau max diag(J^T J) (point 7)
  double function_precision = kMeritPrecision;  ///< eps_f for the difference step (point 9)
};

/// One solve of the run (ADR 0030, point 12).
struct OptimIteration {
  int k = 0;  ///< number of this solve (1, 2, ...)
  /// Normalised merit phi = 2 F / sum of the operand and generator weights, of the current
  /// (last accepted) state after this solve; the mean of w (v - t)^2, a generator counting with
  /// w times its mean square. 0 if all weights are 0.
  double phi = 0.0;
  double mu = 0.0;              ///< damping of this solve
  double rho = 0.0;             ///< gain ratio; NaN without a valid trial
  double step_norm = 0.0;       ///< ||D h|| in scaled internal variables; NaN if no step
  bool accepted = false;        ///< the trial was accepted
  std::size_t evaluations = 0;  ///< evaluations of the solver so far (cumulative)
};

/// One operand in the final state.
struct OperandValue {
  std::string pointer;  ///< /optimization/operands/i
  double value = 0.0;   ///< value of the operand, unit of the operand
  double target = 0.0;  ///< target of the operand, unit of the operand
  double weight = 0.0;  ///< weight w >= 0, dimensionless
  /// Share w (v - t)^2 / sum_j w_j (v_j - t_j)^2 of the merit, percent; 0 if the sum is 0.
  double contribution = 0.0;
};

/// One generator in the final state (ADR 0030, addendum #168).
struct GeneratorValue {
  std::string pointer;  ///< /optimization/generators/g
  /// Weighted RMS sqrt(GeneratorStats::mean_square): mm (rms_spot) or waves at the reference
  /// wavelength (rms_wavefront); independent of the weight, so also reported for weight 0.
  double rms = 0.0;
  double weight = 0.0;  ///< weight w >= 0, dimensionless
  /// Share of its residuals in the merit, sum r^2 / sum_j f_j^2, percent; 0 if the sum is 0.
  double contribution = 0.0;
  std::size_t rays_launched = 0;  ///< rays of the generator
  std::size_t rays_lost = 0;      ///< rays that did not arrive (residuals 0)
};

/// One variable in the final state.
struct VariableValue {
  std::string pointer;  ///< JSON pointer of the value in the edit form (Variable::pointer)
  std::string row;      ///< name of the table row; empty for a model Param
  std::optional<std::size_t> configuration;  ///< column of a `values` row
  double start = 0.0;                        ///< input value
  double end = 0.0;                          ///< result value (bitwise the input if unchanged)
  bool changed = false;                      ///< the run changed this variable
  bool at_bound = false;                     ///< the result lies at a bound (bounds.hpp)
};

/// Result of optimize() (ADR 0030, point 12).
struct OptimResult {
  LmStatus status = LmStatus::MaxIterations;
  /// The input with the result values of the changed variables; everything else, bound
  /// Params and unchanged variables included, is bitwise the input.
  model::System system;
  /// RFC 6902 patch from the input to `system`: one "replace" per changed variable, numbers in
  /// the shortest round-trip form; "[]" if nothing changed.
  std::string patch;
  std::vector<OptimIteration> history;     ///< one entry per solve
  std::vector<OperandValue> operands;      ///< final state, file order
  std::vector<GeneratorValue> generators;  ///< final state, file order
  std::vector<VariableValue> variables;    ///< order of ADR 0030, point 5
  /// optim.parameter_at_bound (per variable), optim.evaluation_failed, optim.jacobian_failed
  /// (status Failed, error), optim.rays_lost (per generator with lost rays in the final state),
  /// then the warnings of the analyses in the final state (rays.lost, stop.clips_beam).
  std::vector<model::Diagnostic> diagnostics;
  int iterations = 0;  ///< number of solves
  /// Evaluations of the merit function in this call: the solver's, the start check and the
  /// final state (not after a cancellation; 0 if cancelled before the start).
  std::size_t evaluations = 0;
  std::size_t failed_evaluations = 0;  ///< trials that were invalid (ADR 0030, point 10)
};

/// Optimizes the variables of `system` against its merit function System::optimization
/// (ADR 0030). Every evaluation sets the values in a copy, compiles every used configuration
/// and evaluates the operands (MeritFunction). The input is not changed.
///
/// Known limit (ADR 0030, addendum #170): an equality condition stated as an operand with a
/// large weight together with a residual that cannot vanish (e.g. a fixed EFL and a minimal
/// RMS spot) makes a narrow, curved valley; the run may then end with ConvergedStep far from the
/// optimum. Hold such conditions exactly through the parameter table (an expression row binding
/// the dependent value, ADR 0029), as in tests/reference/m5/singlet_solve.rtt.json; exact
/// equality constraints in the solver are #196.
///
/// @param system    a valid system with at least one variable and one operand or generator
/// @param materials material library for compile
/// @param coatings  coating library for compile, or nullptr
/// @param options   solver settings
/// @param control   cancellation and progress (stages "jacobian" and "optimize"); a
///                  cancellation ends the run with status Cancelled and the last accepted state.
///                  Requested before the start, nothing is evaluated and the result is the
///                  input; after a cancellation the final state is not evaluated again, so
///                  `operands`, `generators`, optim.rays_lost and the warnings of the analyses
///                  stay empty (ADR 0030, addendum #167)
/// @throws compile::CompileError if the system is invalid or does not compile at the start
/// @throws OptimError with optim.no_variables and/or optim.no_operands (both if both apply), or
///         with merit.operand_unsupported
/// @throws the exception of an analysis at the start (ParaxialError, AnalysisError,
///         NoStopError), or std::invalid_argument naming the operand or generator if its value
///         or a residual is not defined at the start (ADR 0030, point 10)
[[nodiscard]] OptimResult optimize(const model::System& system,
                                   const material::MaterialLibrary& materials,
                                   const coating::CoatingLibrary* coatings = nullptr,
                                   const OptimizeOptions& options = {},
                                   const trace::RunControl& control = {});

}  // namespace rtt::optim
