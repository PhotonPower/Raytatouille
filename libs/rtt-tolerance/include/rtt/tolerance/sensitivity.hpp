#pragma once

/// @file sensitivity.hpp
/// Sensitivity, RSS and inverse sensitivity of a tolerance set (ADR 0036, point 6; #258).
///
/// Every evaluation is one trial of the TrialContext (trial.hpp, #257): the item at +T_i or
/// -T_i, all other items at 0, after the compensators, with the nominal reference. Phi_0 is the
/// nominal evaluation of the context (TrialContext::nominal(), the same procedure). Values are
/// in the unit of the criterion's operand (mm, waves, dimensionless), limits in the unit of the
/// item's kind (mm, degree, dimensionless).
///
/// Sources: J. H. Burge, "Tolerancing Optical Systems", OPTI 521, slides 5 to 11 (RSS, sensitivity
/// by perturbation, inverse sensitivity, tolerance as +-2 sigma); docs/quellen.md.

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "rtt/coating/catalog.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/optimization.hpp"
#include "rtt/model/system.hpp"
#include "rtt/model/tolerance.hpp"
#include "rtt/model/validate.hpp"
#include "rtt/tolerance/expand.hpp"
#include "rtt/tolerance/trial.hpp"
#include "rtt/trace/run_control.hpp"

namespace rtt::tolerance {

/// How the contributions of a criterion add (ADR 0036, point 6): spot_rms and opd_rms are
/// quadratic (non-negative RMS quantities), all other operands linear (signed quantities).
enum class CriterionKind : std::uint8_t { Linear, Quadratic };

/// The kind of a criterion with this operand.
[[nodiscard]] CriterionKind criterion_kind(const model::Operand& operand) noexcept;

/// c = sigma / T of a distribution (ADR 0036, points 5 and 6): truncated normal
/// sqrt(1 - 2 a phi(a) / Z) / 2 with a = 2, Z = erf(a / sqrt(2)) (0.43981), uniform 1/sqrt(3),
/// endpoint 1. Dimensionless.
/// @throws std::invalid_argument for a value outside the enumeration
[[nodiscard]] double sigma_ratio(model::ToleranceDistribution distribution);

/// Options of sensitivity().
struct SensitivityOptions {
  /// Confirmation run at +-T'_i after the inverse sensitivity (ADR 0036, point 6).
  bool confirm = true;
  /// The trials (compensator run); only for the overload with a System, which prepares the
  /// context itself. The overload with a TrialContext uses the options of its prepare().
  TrialOptions trial;
};

/// Whether the run finished (ADR 0036, point 4; addendum #258, F2).
enum class SensitivityStatus : std::uint8_t {
  Completed,
  /// Partial: only items with both sides; no RSS, no inverse, no confirmation. A cancellation
  /// during the confirmation keeps the complete first part and drops only the confirmation.
  Cancelled,
};

/// One item for one criterion. All values in the unit of the criterion.
struct ItemContribution {
  std::optional<double> plus;   ///< Phi(+T_i); none: not defined (reason in `excluded`)
  std::optional<double> minus;  ///< Phi(-T_i); none: not defined (reason in `excluded`)
  /// Not empty if the item is excluded for this criterion (addendum #258, F1): a side is not
  /// defined (its reason), Phi_0 is not defined, or the run was cancelled before both sides
  /// were done. Then all values below are none.
  std::string excluded;
  /// Phi(+T_i) - Phi_0 or Phi(-T_i) - Phi_0, whichever has the larger magnitude (the plus side
  /// on a tie); for display only.
  std::optional<double> worst;
  std::optional<double> delta;  ///< linear: (Phi(+T_i) - Phi(-T_i)) / 2, signed
  std::optional<double> shift;  ///< linear: (Phi(+T_i) + Phi(-T_i)) / 2 - Phi_0, second order
  std::optional<double> q;      ///< quadratic: (Phi(+T_i)^2 + Phi(-T_i)^2) / 2 - Phi_0^2, signed
};

/// Inverse sensitivity of one criterion (ADR 0036, point 6; Burge, slide 6).
struct InverseSensitivity {
  /// Allowance A in the unit of the criterion: linear the smallest distance of Phi_0 to min, max
  /// and max_change; quadratic sqrt(U^2 - Phi_0^2) with U = min(max, Phi_0 + max_change). 0 if
  /// Phi_0 lies on or outside a bound; none if there is no bound to distribute or Phi_0 is not
  /// defined. Only a positive A is distributed (`share`, `limits`).
  std::optional<double> allowance;
  /// Why the criterion is not distributed: "no upper bound" (quadratic), "no bound" (linear),
  /// "nominal value on or outside a bound" (A = 0), "nominal value not defined",
  /// "no items".
  std::string note;
  /// Share of one item: linear A / sqrt(n) in the unit of the criterion, quadratic A^2 / n in its
  /// square; n the number of items of the context (excluded items count, addendum #258).
  std::optional<double> share;
  /// T'_i per item in the unit of the item: linear T_i share / |delta_i|, quadratic
  /// T_i sqrt(share / q_i); may exceed T_i. None: the item limits nothing or is excluded.
  /// Empty if the criterion is not distributed (see `note`).
  std::vector<std::optional<double>> limits;
  /// Items that limit nothing, ascending (ADR 0036, addendum #258): linear
  /// |delta_i| <= eps_M max(|Phi_0|, |Phi(+T_i)|, |Phi(-T_i)|, A), quadratic
  /// q_i <= eps_M max(Phi_0^2, Phi(+T_i)^2, Phi(-T_i)^2, A^2), eps_M = optim::kMeritPrecision.
  std::vector<std::size_t> unlimited;
};

/// Sensitivity of one criterion.
struct CriterionSensitivity {
  std::string id;  ///< the criterion's id
  CriterionKind kind = CriterionKind::Linear;
  /// Phi_0; none: not defined, then no contribution, no RSS and no inverse (addendum #258, F1).
  std::optional<double> nominal;
  std::string nominal_reason;           ///< reason if `nominal` is none
  std::vector<ItemContribution> items;  ///< per item, in the order of the context
  /// Items excluded for this criterion (F1), ascending; RSS and prediction leave them out and
  /// are marked incomplete.
  std::vector<std::size_t> excluded;
  bool complete = true;  ///< `excluded` is empty
  /// RSS (Burge, slides 5 and 9): linear sqrt(sum delta_i^2); quadratic
  /// Phi_RSS = sqrt(max(0, Phi_0^2 + sum q_i)) (warning tolerances.rss_negative if the sum is
  /// negative). None if cancelled or Phi_0 is not defined.
  std::optional<double> rss;
  /// Linear only: sigma_RSS = sqrt(sum (c_i delta_i)^2), the predicted standard deviation, unit
  /// of the criterion. None for a quadratic criterion and as `rss`.
  std::optional<double> predicted_sigma;
  /// Quadratic only: E[Phi^2]_RSS = max(0, Phi_0^2 + sum c_i^2 q_i), the predicted mean of
  /// Phi^2, in the unit of the criterion squared (mm^2, waves^2); a negative sum is reported as 0
  /// with the warning tolerances.rss_negative (it can be negative although Phi_0^2 + sum q_i is
  /// not). None for a linear criterion and as `rss`.
  std::optional<double> predicted_mean_square;
  InverseSensitivity inverse;  ///< empty if cancelled or Phi_0 is not defined
};

/// Confirmation run at +-T'_i (ADR 0036, point 6).
struct Confirmation {
  /// The limits used, per item: the smallest T'_i over the criteria that limit it (addendum
  /// #258, F3), possibly larger than T_i; T_i for an item that no criterion limits, also if it
  /// is excluded in every criterion or its criteria are not distributed.
  std::vector<double> limits;
  std::vector<TrialResult> plus;   ///< per item, the trial at +limits[i]
  std::vector<TrialResult> minus;  ///< per item, the trial at -limits[i]
  /// Contributions and RSS at the limits (`inverse` empty), per criterion.
  std::vector<CriterionSensitivity> criteria;
};

/// The trials of one item.
struct ItemTrials {
  model::ToleranceItem item;
  double sigma_ratio = 1.0;          ///< c_i of the item's distribution
  std::optional<TrialResult> plus;   ///< at +T_i; none if cancelled before it was done
  std::optional<TrialResult> minus;  ///< at -T_i; none if cancelled before it was done
};

/// Result of sensitivity() (ADR 0036, point 10).
struct ToleranceSensitivity {
  TrialResult nominal;                         ///< Phi_0 (TrialContext::nominal())
  std::vector<ItemTrials> items;               ///< per item, order of the context
  std::vector<CriterionSensitivity> criteria;  ///< per criterion, order of the context
  /// Per item the smallest T'_i over the criteria (F3); none if no criterion limits the item,
  /// and none for every item if cancelled.
  std::vector<std::optional<double>> limits;
  std::optional<Confirmation> confirmation;  ///< with SensitivityOptions::confirm, not cancelled
  /// Warnings of the run itself: tolerances.rss_negative at /tolerances/criteria/<k> (of the
  /// first part, then of the confirmation; the message names the part). The warnings of each trial
  /// are in its TrialResult.
  std::vector<model::Diagnostic> warnings;
  /// Only from the overload with a System: the result of expand(), whose items are those of
  /// `items` (ADR 0036, addendum #260).
  std::optional<ToleranceItems> expansion;
  SensitivityStatus status = SensitivityStatus::Completed;
};

/// Sensitivity, RSS, inverse sensitivity and the confirmation run of the items of `context`
/// (ADR 0036, point 6). The 2n evaluations (and 2n of the confirmation) run in parallel in fixed
/// slots and are reduced in item order: bitwise the same for every number of threads.
/// RunControl: stages "sensitivity" and "confirmation", one work item per trial; a cancellation
/// stops between trials (and inside a compensator run) and gives a partial result (status
/// Cancelled).
/// @throws whatever the progress callback throws (rtt::trace::RunControl)
[[nodiscard]] ToleranceSensitivity sensitivity(const TrialContext& context,
                                               const SensitivityOptions& options = {},
                                               const trace::RunControl* control = nullptr);

/// The same for all items of `system` (ADR 0036, addendum #260): expand(system, materials),
/// then prepare() with these items and options.trial, then the run above; the result carries
/// the expansion. The cancellation of `control` also covers Phi_0 in prepare().
/// @param materials material library (also for the Abbe limits of expand)
/// @param coatings  coating library, or nullptr
/// @throws compile::CompileError, std::invalid_argument, optim::OptimError as expand() and
///         prepare(); whatever the progress callback throws (rtt::trace::RunControl)
[[nodiscard]] ToleranceSensitivity sensitivity(const model::System& system,
                                               const material::MaterialLibrary& materials,
                                               const coating::CoatingLibrary* coatings,
                                               const SensitivityOptions& options = {},
                                               const trace::RunControl* control = nullptr);

}  // namespace rtt::tolerance
