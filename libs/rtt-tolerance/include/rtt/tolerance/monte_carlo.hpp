#pragma once

/// @file monte_carlo.hpp
/// Monte Carlo of a tolerance set (ADR 0036, point 7; #259).
///
/// Trial n (0-based) perturbs every item of the TrialContext (trial.hpp, #257) with
/// perturbations(seed, n, items) (random.hpp, ADR 0036, point 5) and is one evaluate(): after the
/// compensators, with the nominal reference. Criterion values are in the unit of the criterion's
/// operand (mm, waves, dimensionless), perturbations in the unit of the item's kind (mm, degree,
/// dimensionless).
///
/// Sources: NIST/SEMATECH e-Handbook, Sec. 7.2.4.1 (Wilson interval); docs/quellen.md.

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "rtt/coating/catalog.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/system.hpp"
#include "rtt/model/tolerance.hpp"
#include "rtt/model/validate.hpp"
#include "rtt/tolerance/expand.hpp"
#include "rtt/tolerance/trial.hpp"
#include "rtt/trace/run_control.hpp"

namespace rtt::tolerance {

/// Options of monte_carlo().
struct MonteCarloOptions {
  /// The seed (ADR 0035, point 8; the full range of uint64). None: tolerances.monte_carlo.seed of
  /// the nominal system, else 0.
  std::optional<std::uint64_t> seed;
  /// The number of trials N in [1, 2^32]. None: tolerances.monte_carlo.trials of the nominal
  /// system; there is no default. The run holds one slot per trial, so its memory grows with N.
  std::optional<std::uint64_t> trials;
  /// The percentiles p in (0, 1] of each criterion, dimensionless (nearest rank, ADR 0036,
  /// point 7).
  std::vector<double> percentiles = {0.5, 0.9, 0.95, 0.99};
  /// z = z_{1 - alpha/2} > 0 of the Wilson interval, dimensionless; the default is 95 %.
  double z = 1.959963984540054;
  /// Keep every trial in MonteCarloResult::trials. Off: only the trials that are not Valid (they
  /// are always listed, ADR 0036, point 7).
  bool keep_trials = true;
  /// The trials (compensator run); only for the overload with a System, which prepares the
  /// context itself. The overload with a TrialContext uses the options of its prepare().
  TrialOptions trial;
};

/// Whether the run finished (ADR 0036, point 4).
enum class MonteCarloStatus : std::uint8_t {
  Completed,
  /// Only the trials that finished before the cancellation; which ones depends on the timing, so
  /// a cancelled run is not reproducible (ADR 0036, point 7).
  Cancelled,
};

/// One trial.
struct TrialRecord {
  std::uint32_t trial = 0;     ///< its number n, 0-based
  std::vector<double> values;  ///< the perturbation per item, in the unit of the item
  /// The evaluation: status (Valid or Invalid), reason, criterion values with `satisfied`,
  /// compensators, LM status, mu and the warnings of this trial.
  TrialResult result;
  /// The trial is good: Valid and every criterion satisfied (ADR 0035, point 7).
  bool good = false;
};

/// A yield with its Wilson interval (ADR 0036, point 7; NIST e-Handbook, Sec. 7.2.4.1).
struct Yield {
  std::size_t good = 0;   ///< trials counted as good
  std::size_t total = 0;  ///< N, the finished trials
  /// good / total and the interval, dimensionless in [0, 1]; none for total = 0. lower is exactly
  /// 0 for good = 0, upper exactly 1 for good = total.
  std::optional<double> fraction;
  std::optional<double> lower;
  std::optional<double> upper;
};

/// One percentile: Q(p) = x_(r) of the ascending values, r = ceil(p n) (ADR 0036, point 7).
struct Percentile {
  double p = 0.0;        ///< the requested p in (0, 1], dimensionless
  std::size_t rank = 0;  ///< r in [1, n], 1-based
  double value = 0.0;    ///< unit of the criterion
};

/// The statistics of one criterion (ADR 0036, point 7; addendum #259).
struct CriterionStatistics {
  std::string id;  ///< the criterion's id
  /// Finished trials in which the criterion has a value, and those in which it has none. Both
  /// count over all N trials, also the Invalid ones (addendum #259, D5).
  std::size_t with_value = 0;
  std::size_t without_value = 0;
  /// n: the values in the statistics, those of the Valid trials (ADR 0036, point 7).
  std::size_t count = 0;
  /// Mean, standard deviation (n - 1, two passes in trial order, from n = 2), minimum, maximum;
  /// unit of the criterion. None if n is too small.
  std::optional<double> mean;
  std::optional<double> standard_deviation;
  std::optional<double> min;
  std::optional<double> max;
  std::vector<Percentile> percentiles;  ///< order of MonteCarloOptions::percentiles; empty if n = 0
  /// Good: the criterion is satisfied (CriterionValue::satisfied); a trial without a value counts
  /// as bad. Denominator N.
  Yield yield;
  std::string note;  ///< "no valid values" for n = 0
};

/// Result of monte_carlo() (ADR 0036, point 10).
struct MonteCarloResult {
  std::uint64_t seed = 0;                     ///< the seed used (option, file or 0)
  std::uint64_t requested = 0;                ///< the number of trials asked for
  TrialResult nominal;                        ///< Phi_0 (TrialContext::nominal())
  std::vector<model::ToleranceItem> items;    ///< order of the values of each trial
  std::vector<CriterionStatistics> criteria;  ///< per criterion, order of the context
  /// Good: Valid and every criterion satisfied; an Invalid trial counts as bad. Denominator N.
  Yield yield;
  /// The finished trials in ascending order: all with keep_trials, otherwise those that are not
  /// Valid.
  std::vector<TrialRecord> trials;
  /// The warnings of the trials, one per code and location in the order of their first trial:
  /// that trial's diagnostic, with the number of trials and the first trial in the message
  /// (addendum #259, D7). Each trial's own warnings are in its TrialResult.
  std::vector<model::Diagnostic> warnings;
  /// Only from the overload with a System: the result of expand() (ADR 0036, addendum #260).
  std::optional<ToleranceItems> expansion;
  MonteCarloStatus status = MonteCarloStatus::Completed;
};

/// Monte Carlo over the items of `context` (ADR 0036, point 7). The trials run in parallel in
/// fixed slots and are reduced in trial order: bitwise the same for every number of threads and
/// every repetition (point 8). RunControl: stage "monte_carlo", one work item per trial; a
/// cancellation stops between trials (and inside a compensator run) and keeps the finished trials
/// (status Cancelled). If nominal() of the context is Cancelled, no trial runs.
/// @throws std::invalid_argument if the number of trials is neither given nor in the file, or
///         lies outside [1, 2^32]; for a percentile outside (0, 1] or a z that is not finite and
///         positive; for an item whose limit is not finite and > 0 or whose distribution is
///         unknown (perturbation()); before any trial
/// @throws whatever the progress callback throws (rtt::trace::RunControl)
[[nodiscard]] MonteCarloResult monte_carlo(const TrialContext& context,
                                           const MonteCarloOptions& options = {},
                                           const trace::RunControl* control = nullptr);

/// The same for all items of `system` (ADR 0036, addendum #260): expand(system, materials), then
/// prepare() with these items and options.trial, then the run above; the result carries the
/// expansion. The cancellation of `control` also covers Phi_0 in prepare().
/// @param materials material library (also for the Abbe limits of expand)
/// @param coatings  coating library, or nullptr
/// @throws compile::CompileError, std::invalid_argument, optim::OptimError as expand() and
///         prepare(); std::invalid_argument as above; whatever the progress callback throws
[[nodiscard]] MonteCarloResult monte_carlo(const model::System& system,
                                           const material::MaterialLibrary& materials,
                                           const coating::CoatingLibrary* coatings,
                                           const MonteCarloOptions& options = {},
                                           const trace::RunControl* control = nullptr);

}  // namespace rtt::tolerance
