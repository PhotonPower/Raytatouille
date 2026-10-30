#include "rtt/tolerance/monte_carlo.hpp"

#include <oneapi/tbb/blocked_range.h>
#include <oneapi/tbb/parallel_for.h>
#include <oneapi/tbb/partitioner.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "monte_carlo_detail.hpp"
#include "rtt/tolerance/random.hpp"

namespace rtt::tolerance {

namespace detail {

std::size_t nearest_rank(double p, std::size_t n) noexcept {
  const double x = p * static_cast<double>(n);
  const double k = std::round(x);
  const double r =
      std::abs(x - k) <= 4.0 * std::numeric_limits<double>::epsilon() * x ? k : std::ceil(x);
  if (r < 1.0) return 1;
  if (r >= static_cast<double>(n)) return n;
  return static_cast<std::size_t>(r);
}

Yield wilson(std::size_t good, std::size_t total, double z) noexcept {
  Yield y;
  y.good = good;
  y.total = total;
  if (total == 0) return y;
  // NIST e-Handbook, Sec. 7.2.4.1, with z = z_{1 - alpha/2} > 0 in both bounds.
  const auto n = static_cast<double>(total);
  const double p = static_cast<double>(good) / n;
  const double z2 = z * z;
  const double denominator = 1.0 + z2 / n;
  const double centre = p + z2 / (2.0 * n);
  const double half = z * std::sqrt(p * (1.0 - p) / n + z2 / (4.0 * n * n));
  y.fraction = p;
  y.lower = good == 0 ? 0.0 : std::clamp((centre - half) / denominator, 0.0, 1.0);
  y.upper = good == total ? 1.0 : std::clamp((centre + half) / denominator, 0.0, 1.0);
  return y;
}

Moments moments(std::span<const double> x) noexcept {
  Moments m;
  if (x.empty()) return m;
  const auto n = static_cast<double>(x.size());
  double sum = 0.0;
  for (const double v : x) sum += v;
  const double mean = sum / n;
  m.mean = mean;
  if (x.size() < 2) return m;
  double squares = 0.0;
  for (const double v : x) squares += (v - mean) * (v - mean);
  m.standard_deviation = std::sqrt(squares / (n - 1.0));
  return m;
}

}  // namespace detail

namespace {

constexpr std::uint64_t kMaxTrials = std::uint64_t{1} << 32;  // ADR 0035, point 8

/// The checked options: seed and number of trials resolved (option before file, D6).
struct Resolved {
  std::uint64_t seed = 0;
  std::uint64_t trials = 0;
};

/// @throws std::invalid_argument as monte_carlo()
Resolved resolve(const std::optional<model::MonteCarloDefaults>& file,
                 const MonteCarloOptions& options) {
  Resolved r;
  if (options.seed) {
    r.seed = *options.seed;
  } else if (file && file->seed) {
    r.seed = *file->seed;
  }
  std::optional<std::uint64_t> trials = options.trials;
  if (!trials && file) trials = file->trials;
  if (!trials) {
    throw std::invalid_argument(
        "monte_carlo: the number of trials is neither an option nor in tolerances.monte_carlo");
  }
  if (*trials < 1 || *trials > kMaxTrials) {
    throw std::invalid_argument("monte_carlo: trials must lie in [1, 2^32], got " +
                                std::to_string(*trials));
  }
  r.trials = *trials;
  for (const double p : options.percentiles) {
    if (!(p > 0.0 && p <= 1.0)) {
      throw std::invalid_argument("monte_carlo: a percentile must lie in (0, 1], got " +
                                  std::to_string(p));
    }
  }
  if (!(std::isfinite(options.z) && options.z > 0.0)) {
    throw std::invalid_argument("monte_carlo: z must be finite and positive, got " +
                                std::to_string(options.z));
  }
  return r;
}

/// The trials in parallel, each in its own slot (ADR 0004; ADR 0036, point 4): bitwise the same
/// for every number of threads. A cancellation leaves the remaining slots empty.
/// @throws whatever the progress callback throws (RunMonitor::finish)
std::vector<std::optional<TrialRecord>> run_trials(const TrialContext& context,
                                                   std::uint64_t seed,
                                                   std::size_t count,
                                                   const trace::RunControl& control,
                                                   bool& cancelled) {
  std::vector<std::optional<TrialRecord>> slots(count);
  trace::RunMonitor monitor(control, count, "monte_carlo");
  std::atomic<bool> trial_cancelled{false};
  const auto body = [&](const oneapi::tbb::blocked_range<std::size_t>& range) {
    for (std::size_t idx = range.begin(); idx != range.end(); ++idx) {
      if (monitor.stop() || trial_cancelled.load()) continue;  // before every trial
      TrialRecord record;
      record.trial = static_cast<std::uint32_t>(idx);  // idx < 2^32 (resolve)
      record.values = perturbations(seed, record.trial, context.items());
      record.result = evaluate(context, record.values, &control);
      if (record.result.status == TrialStatus::Cancelled) {
        trial_cancelled = true;
        continue;
      }
      record.good =
          record.result.status == TrialStatus::Valid &&
          std::all_of(record.result.criteria.begin(), record.result.criteria.end(),
                      [](const CriterionValue& c) { return c.satisfied.value_or(false); });
      slots[idx] = std::move(record);
      monitor.add(1);
    }
  };
  const oneapi::tbb::blocked_range<std::size_t> all(0, count, 1);
  if (control.active()) {
    // Grain 1: the cancellation is checked before every trial (ADR 0004, addendum #83).
    oneapi::tbb::parallel_for(all, body, oneapi::tbb::simple_partitioner());
  } else {
    oneapi::tbb::parallel_for(all, body, oneapi::tbb::static_partitioner());
  }
  cancelled = trial_cancelled.load() ||
              std::any_of(slots.begin(), slots.end(),
                          [](const std::optional<TrialRecord>& t) { return !t.has_value(); });
  try {
    monitor.finish();  // rethrows an exception of the progress callback
  } catch (const trace::Cancelled&) {
    cancelled = true;  // also a request after the last trial
  }
  return slots;
}

/// The statistics of criterion k over the finished trials (ADR 0036, point 7; addendum #259, D5).
CriterionStatistics statistics(const model::Criterion& criterion,
                               std::size_t k,
                               std::span<const TrialRecord> trials,
                               const MonteCarloOptions& options) {
  CriterionStatistics c;
  c.id = criterion.id;
  std::vector<double> values;
  std::size_t good = 0;
  for (const TrialRecord& t : trials) {
    const CriterionValue& v = t.result.criteria[k];
    if (v.value) {
      ++c.with_value;
      if (t.result.status == TrialStatus::Valid) values.push_back(*v.value);
    } else {
      ++c.without_value;
    }
    good += v.satisfied.value_or(false) ? 1 : 0;
  }
  c.count = values.size();
  c.yield = detail::wilson(good, trials.size(), options.z);
  if (values.empty()) {
    c.note = "no valid values";
    return c;
  }
  const detail::Moments m = detail::moments(values);
  c.mean = m.mean;
  c.standard_deviation = m.standard_deviation;
  std::sort(values.begin(), values.end());
  c.min = values.front();
  c.max = values.back();
  for (const double p : options.percentiles) {
    const std::size_t rank = detail::nearest_rank(p, values.size());
    c.percentiles.push_back({p, rank, values[rank - 1]});
  }
  return c;
}

/// The warnings of the trials, once per code and location in the order of their first
/// appearance (trial order, then within the trial), with the count (addendum #259, D7).
std::vector<model::Diagnostic> collect_warnings(std::span<const TrialRecord> trials) {
  struct Entry {
    model::Diagnostic first;
    std::uint32_t trial = 0;
    std::size_t count = 0;
  };
  std::vector<Entry> entries;
  const auto find = [&](const model::Diagnostic& w) {
    return std::find_if(entries.begin(), entries.end(), [&](const Entry& e) {
      return e.first.code == w.code && e.first.location == w.location;
    });
  };
  for (const TrialRecord& t : trials) {
    std::vector<const model::Diagnostic*> in_trial;
    for (const model::Diagnostic& w : t.result.warnings) {
      const bool repeated = std::any_of(in_trial.begin(), in_trial.end(), [&](const auto* o) {
        return o->code == w.code && o->location == w.location;
      });
      if (repeated) continue;
      in_trial.push_back(&w);
      const auto it = find(w);
      if (it == entries.end()) {
        entries.push_back({w, t.trial, 1});
      } else {
        ++it->count;
      }
    }
  }
  std::vector<model::Diagnostic> out;
  out.reserve(entries.size());
  for (Entry& e : entries) {
    e.first.message += " (in " + std::to_string(e.count) + " of " + std::to_string(trials.size()) +
                       " trials, first trial " + std::to_string(e.trial) + ")";
    out.push_back(std::move(e.first));
  }
  return out;
}

}  // namespace

MonteCarloResult monte_carlo(const TrialContext& context,
                             const MonteCarloOptions& options,
                             const trace::RunControl* control) {
  const Resolved resolved = resolve(context.monte_carlo_defaults(), options);
  const trace::RunControl none;
  const trace::RunControl& run = control != nullptr ? *control : none;
  MonteCarloResult out;
  out.seed = resolved.seed;
  out.requested = resolved.trials;
  out.nominal = context.nominal();
  out.items.assign(context.items().begin(), context.items().end());

  std::vector<TrialRecord> finished;
  if (out.nominal.status == TrialStatus::Cancelled) {
    out.status = MonteCarloStatus::Cancelled;
  } else {
    bool cancelled = false;
    std::vector<std::optional<TrialRecord>> slots = run_trials(
        context, resolved.seed, static_cast<std::size_t>(resolved.trials), run, cancelled);
    if (cancelled) out.status = MonteCarloStatus::Cancelled;
    for (std::optional<TrialRecord>& t : slots) {
      if (t) finished.push_back(std::move(*t));
    }
  }

  const std::span<const model::Criterion> criteria = context.criteria();
  for (std::size_t k = 0; k < criteria.size(); ++k) {
    out.criteria.push_back(statistics(criteria[k], k, finished, options));
  }
  const auto good = static_cast<std::size_t>(
      std::count_if(finished.begin(), finished.end(), [](const TrialRecord& t) { return t.good; }));
  out.yield = detail::wilson(good, finished.size(), options.z);
  out.warnings = collect_warnings(finished);
  for (TrialRecord& t : finished) {
    if (options.keep_trials || t.result.status != TrialStatus::Valid) {
      out.trials.push_back(std::move(t));
    }
  }
  return out;
}

MonteCarloResult monte_carlo(const model::System& system,
                             const material::MaterialLibrary& materials,
                             const coating::CoatingLibrary* coatings,
                             const MonteCarloOptions& options,
                             const trace::RunControl* control) {
  resolve(system.tolerances.monte_carlo, options);  // throws before the expensive part
  ToleranceItems expansion = expand(system, materials);
  const TrialContext context =
      prepare(system, materials, coatings, expansion.items, options.trial, control);
  MonteCarloResult out = monte_carlo(context, options, control);
  out.expansion = std::move(expansion);
  return out;
}

}  // namespace rtt::tolerance
