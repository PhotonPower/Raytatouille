#pragma once

/// @file trial.hpp
/// One tolerance trial (ADR 0036, points 1 to 4; #257): the nominal system with given
/// perturbation values for its items, the compensators run, the criteria evaluated with the
/// nominal reference. Sensitivity (#258) and Monte Carlo (#259) are built only on this interface;
/// it knows no seed and draws nothing.
///
/// Units: a perturbation value is in the unit of its item's kind (mm, degree, dimensionless; ADR
/// 0035, point 2), a criterion value in the unit of its operand (mm, waves, dimensionless; ADR
/// 0030, point 3), a compensator value in the unit of its target.

#include <cstdint>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "rtt/coating/catalog.hpp"
#include "rtt/compile/perturbation.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/system.hpp"
#include "rtt/model/tolerance.hpp"
#include "rtt/model/validate.hpp"
#include "rtt/optim/levenberg_marquardt.hpp"
#include "rtt/optim/optimize.hpp"
#include "rtt/trace/run_control.hpp"

namespace rtt::tolerance {

/// Options of the trials of a context.
struct TrialOptions {
  /// The Levenberg-Marquardt run of the compensators (ADR 0036, point 3; ADR 0030, point 8).
  optim::OptimizeOptions compensators;
};

/// How a trial ended.
enum class TrialStatus : std::uint8_t {
  Valid,      ///< every criterion has a value
  Invalid,    ///< the trial counts as invalid in Monte Carlo (ADR 0036, point 7); see TrialResult
  Cancelled,  ///< cancelled through the RunControl; the result is not usable
};

/// The value of one criterion in one trial.
struct CriterionValue {
  /// The value Phi of the criterion's operand, raw (not a change), in the unit of the operand;
  /// none if it is not defined in this trial (an analysis error, NaN, a lost chief ray).
  std::optional<double> value;
  std::string reason;  ///< why `value` is none; empty if it has a value
  /// satisfies(criterion, value, nominal value); none if that cannot be decided.
  std::optional<bool> satisfied;
};

/// One compensator in one trial (ADR 0036, point 3).
struct CompensatorValue {
  double start = 0.0;       ///< the perturbed value before the run, unit of the target
  double end = 0.0;         ///< the value after the run (start if not run or infeasible)
  bool at_bound = false;    ///< the end lies at a bound (travel or Param.min/max)
  bool infeasible = false;  ///< start not inside (min, max): kept at start, with a warning
};

/// The result of one trial.
///
/// criteria[k].value is what counts; status Invalid only means that Monte Carlo counts the
/// trial as invalid. It is set
/// - for a radius item whose value makes the radius 0, infinite or of the other sign (ADR 0036,
///   addendum #257) and for an invalid compensator run (LM failed, an analysis error or an
///   undefined value of the compensator merit): then no criterion is evaluated and every value
///   is none;
/// - for a single criterion without value: then the other criteria keep their values.
struct TrialResult {
  TrialStatus status = TrialStatus::Valid;
  std::string reason;                          ///< for Invalid and Cancelled: the first reason
  std::vector<CriterionValue> criteria;        ///< order of TrialContext::criteria()
  std::vector<CompensatorValue> compensators;  ///< order of TrialContext::compensators()
  std::optional<optim::LmStatus> lm_status;    ///< none without a compensator run
  std::optional<double> mu;                    ///< damping at the end of the run (last history)
  int lm_iterations = 0;
  /// tolerances.compensator_infeasible (one per compensator), then the warnings of the
  /// compensator run and of the criteria (rays.lost, stop.clips_beam, ...), as each reports
  /// them: the same code can appear several times (rays.lost of the run and of the criteria).
  std::vector<model::Diagnostic> warnings;
};

namespace detail {
struct ContextData;    // trial.cpp
struct ContextAccess;  // trial.cpp
}  // namespace detail

/// Everything that stays the same over the trials of a tolerance run: the nominal system, the
/// items, criteria and compensators, the nominal basis of the overlays and the nominal result
/// Phi_0. Immutable with shared content: cheap to copy, safe to use from many threads at once.
class TrialContext {
 public:
  /// The items, in the order of the values of evaluate().
  [[nodiscard]] std::span<const model::ToleranceItem> items() const noexcept;
  /// The criteria of the nominal system's tolerance set.
  [[nodiscard]] std::span<const model::Criterion> criteria() const noexcept;
  /// The compensators of the nominal system's tolerance set.
  [[nodiscard]] std::span<const model::Compensator> compensators() const noexcept;
  /// Phi_0: the trial with all values 0, computed once by prepare() with the same procedure
  /// (compensators included; ADR 0036, points 2 and 6). Status Cancelled if prepare() was
  /// cancelled.
  [[nodiscard]] const TrialResult& nominal() const noexcept;

 private:
  friend TrialContext prepare(const model::System& nominal,
                              const material::MaterialLibrary& materials,
                              const coating::CoatingLibrary* coatings,
                              std::span<const model::ToleranceItem> items,
                              const TrialOptions& options,
                              const trace::RunControl* control);
  friend struct detail::ContextAccess;
  explicit TrialContext(std::shared_ptr<const detail::ContextData> data) : data_(std::move(data)) {}
  std::shared_ptr<const detail::ContextData> data_;
};

/// Prepares the trials of `items` on `nominal` (ADR 0036, point 4): the criteria, compensators,
/// compensator merit (compensator_merit, else optimization) and pivots come from
/// `nominal.tolerances`; the items are those of expand() (#260) or any subset. Computes the
/// nominal basis (rtt::optim::nominal_basis) and Phi_0 (nominal()).
/// @param nominal   a valid system; copied
/// @param materials material library; must outlive the context
/// @param coatings  coating library, or nullptr; must outlive the context
/// @param items     the items, explicit or generated; copied
/// @param control   cancellation of the computation of Phi_0 (nominal() is then Cancelled), or
///                  nullptr
/// @throws compile::CompileError if `nominal` is not valid (model::validate), does not compile,
///         or the glass of a material item gives no index or Abbe number at the lines d, F, C
///         (tolerances.material_unsupported, material.unknown, tolerances.abbe_on_const)
/// @throws std::invalid_argument for an item whose target does not exist or does not fit its
///         kind, or a material item whose segment does not fit the element (as compile: the
///         shorthand takes none, a list needs one below its size, a list of equal entries only
///         0; ADR 0035, point 2)
/// @throws optim::OptimError with merit.operand_unsupported for a criterion or an operand of
///         the compensator merit that the nominal system does not support
[[nodiscard]] TrialContext prepare(const model::System& nominal,
                                   const material::MaterialLibrary& materials,
                                   const coating::CoatingLibrary* coatings,
                                   std::span<const model::ToleranceItem> items,
                                   const TrialOptions& options = {},
                                   const trace::RunControl* control = nullptr);

/// One trial (ADR 0036, point 4): the values applied (offsets on Params and rows, the rigid and
/// material overlay), the compensators run (ADR 0036, point 3), the criteria evaluated with the
/// nominal reference (point 2). Reentrant and deterministic: bitwise the same for every number
/// of threads and every repetition.
/// @param values  one perturbation per item, in the unit of the item's kind
/// @param control cancellation: checked before and inside the compensator run and before the
///                criteria; then status Cancelled. Progress is not reported (the caller reports
///                per trial).
/// @throws std::invalid_argument if `values` has not one finite value per item; everything that
///         fails in the trial itself is a status
[[nodiscard]] TrialResult evaluate(const TrialContext& context,
                                   std::span<const double> values,
                                   const trace::RunControl* control = nullptr);

/// Whether `value` satisfies `criterion` (ADR 0035, point 7): min <= value, value <= max and
/// |value - nominal_value| <= max_change, each where given, bounds included. None if `value` is
/// none, or if max_change is given and `nominal_value` is none.
[[nodiscard]] std::optional<bool> satisfies(const model::Criterion& criterion,
                                            std::optional<double> value,
                                            std::optional<double> nominal_value);

/// The perturbed trial before the compensators (for acceptance tests and hand calculations).
struct Applied {
  model::System system;                ///< the nominal model with the offsets (ADR 0036, 1a)
  compile::Perturbation perturbation;  ///< the rigid and material overlay (1b, 1c)
};

/// The values applied to the nominal system, without compensators and criteria.
/// @throws as evaluate() for `values`
[[nodiscard]] Applied apply(const TrialContext& context, std::span<const double> values);

}  // namespace rtt::tolerance
