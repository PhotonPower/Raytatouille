#include "rtt/tolerance/sensitivity.hpp"

#include <oneapi/tbb/blocked_range.h>
#include <oneapi/tbb/parallel_for.h>
#include <oneapi/tbb/partitioner.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

#include "rtt/diagnostics/codes.hpp"
#include "rtt/optim/optimize.hpp"
#include "sensitivity_detail.hpp"

namespace rtt::tolerance {

namespace {

/// The trials of one stage: slot 2 i at +limits[i], slot 2 i + 1 at -limits[i]; none if the
/// slot was not done (cancelled).
struct Slots {
  std::vector<std::optional<TrialResult>> trials;
  bool cancelled = false;
};

/// The 2n trials of one stage in parallel, each in its own slot (ADR 0004): bitwise the same for
/// every number of threads. A cancellation (RunControl or a trial that ends Cancelled) leaves the
/// remaining slots empty.
/// @throws whatever the progress callback throws (RunMonitor::finish)
Slots run_trials(const TrialContext& context,
                 std::span<const double> limits,
                 const trace::RunControl& control,
                 std::string_view stage) {
  const std::size_t n = limits.size();
  const std::size_t count = 2 * n;
  Slots out;
  out.trials.resize(count);
  trace::RunMonitor monitor(control, count, stage);
  std::atomic<bool> trial_cancelled{false};
  const auto body = [&](const oneapi::tbb::blocked_range<std::size_t>& range) {
    std::vector<double> values(n, 0.0);
    for (std::size_t idx = range.begin(); idx != range.end(); ++idx) {
      if (monitor.stop() || trial_cancelled.load()) continue;  // before every evaluation
      const std::size_t i = idx / 2;
      values[i] = idx % 2 == 0 ? limits[i] : -limits[i];
      TrialResult r = evaluate(context, values, &control);
      values[i] = 0.0;
      if (r.status == TrialStatus::Cancelled) {
        trial_cancelled = true;
        continue;
      }
      out.trials[idx] = std::move(r);
      monitor.add(1);
    }
  };
  const oneapi::tbb::blocked_range<std::size_t> all(0, count, 1);
  if (control.active()) {
    // Grain 1: the cancellation is checked before every evaluation (ADR 0004, addendum #83).
    oneapi::tbb::parallel_for(all, body, oneapi::tbb::simple_partitioner());
  } else {
    oneapi::tbb::parallel_for(all, body, oneapi::tbb::static_partitioner());
  }
  out.cancelled = trial_cancelled.load() ||
                  std::any_of(out.trials.begin(), out.trials.end(),
                              [](const std::optional<TrialResult>& t) { return !t.has_value(); });
  try {
    monitor.finish();  // rethrows an exception of the progress callback
  } catch (const trace::Cancelled&) {
    out.cancelled = true;  // also a request after the last evaluation
  }
  return out;
}

/// The contribution of one item to criterion k (ADR 0036, point 6; addendum #258, F1).
ItemContribution contribution(CriterionKind kind,
                              const std::optional<double>& phi0,
                              bool nominal_cancelled,
                              const std::optional<TrialResult>& plus,
                              const std::optional<TrialResult>& minus,
                              std::size_t k) {
  ItemContribution c;
  if (plus) c.plus = plus->criteria[k].value;
  if (minus) c.minus = minus->criteria[k].value;
  if (!phi0) {
    c.excluded = nominal_cancelled ? "cancelled" : "nominal value not defined";
    return c;
  }
  if (!plus || !minus) {
    c.excluded = "cancelled before both sides were done";
    return c;
  }
  if (!c.plus) {
    c.excluded = "+T: " + plus->criteria[k].reason;
    return c;
  }
  if (!c.minus) {
    c.excluded = "-T: " + minus->criteria[k].reason;
    return c;
  }
  const double p = *c.plus;
  const double m = *c.minus;
  const double p0 = *phi0;
  const double up = p - p0;
  const double down = m - p0;
  c.worst = std::abs(down) > std::abs(up) ? down : up;
  if (kind == CriterionKind::Linear) {
    c.delta = (p - m) / 2.0;       // slope times tolerance (Burge, slide 6)
    c.shift = (p + m) / 2.0 - p0;  // second order
  } else {
    c.q = (p * p + m * m) / 2.0 - p0 * p0;  // Burge, slide 7, averaged over both sides
  }
  return c;
}

/// tolerances.rss_negative of criterion k; `stage` names the part of the run ("sensitivity" or
/// "confirmation"), so that the warnings of both parts can be told apart.
model::Diagnostic rss_warning(const model::Criterion& criterion,
                              std::size_t k,
                              std::string_view stage) {
  const diagnostics::DiagnosticCode code("tolerances.rss_negative");
  return {code.severity(), "/tolerances/criteria/" + std::to_string(k),
          "criterion '" + criterion.id + "' (" + std::string(stage) +
              "): Phi_0^2 + sum q_i is negative; the RSS is reported as 0",
          std::string(code.str())};
}

/// A distributed allowance: A and the share of one item, units of the criterion (squared for a
/// quadratic share).
struct Distribution {
  double allowance = 0.0;
  double share = 0.0;
};

/// The allowance A and the share of one item (ADR 0036, point 6). Sets allowance, note and share
/// of `inverse`; returns them if A > 0 is distributed.
std::optional<Distribution> distribute(const model::Criterion& criterion,
                                       CriterionKind kind,
                                       double phi0,
                                       std::size_t n,
                                       InverseSensitivity& inverse) {
  if (n == 0) {
    inverse.note = "no items";
    return std::nullopt;
  }
  const auto items = static_cast<double>(n);
  if (kind == CriterionKind::Linear) {
    std::optional<double> a;
    const auto take = [&](double x) { a = a ? std::min(*a, x) : x; };
    if (criterion.min) take(phi0 - *criterion.min);
    if (criterion.max) take(*criterion.max - phi0);
    if (criterion.max_change) take(*criterion.max_change);
    if (!a) {
      inverse.note = "no bound";
      return std::nullopt;
    }
    const double allowance = *a;
    if (!(allowance > 0.0)) {
      inverse.allowance = 0.0;
      inverse.note = "nominal value on or outside a bound";
      return std::nullopt;
    }
    const Distribution d{allowance, allowance / std::sqrt(items)};
    inverse.allowance = d.allowance;
    inverse.share = d.share;
    return d;
  }
  // Quadratic: the upper bound U = min(max, Phi_0 + max_change), A^2 = U^2 - Phi_0^2.
  std::optional<double> upper = criterion.max;
  if (criterion.max_change) {
    const double x = phi0 + *criterion.max_change;
    upper = upper ? std::min(*upper, x) : x;
  }
  if (!upper) {
    inverse.note = "no upper bound";
    return std::nullopt;
  }
  const double u = *upper;
  const double a2 = u * u - phi0 * phi0;
  if (!(u > phi0) || !(a2 > 0.0)) {
    inverse.allowance = 0.0;
    inverse.note = "nominal value on or outside a bound";
    return std::nullopt;
  }
  const Distribution d{std::sqrt(a2), a2 / items};
  inverse.allowance = d.allowance;
  inverse.share = d.share;
  return d;
}

/// The sensitivity of criterion k over all items: contributions, RSS, prediction and, with
/// `with_inverse`, the inverse sensitivity (ADR 0036, point 6; addendum #258).
CriterionSensitivity criterion_sensitivity(const TrialContext& context,
                                           std::size_t k,
                                           std::span<const double> limits,
                                           std::span<const double> ratios,
                                           const Slots& slots,
                                           bool with_inverse,
                                           std::string_view stage,
                                           std::vector<model::Diagnostic>& warnings) {
  const model::Criterion& criterion = context.criteria()[k];
  const CriterionValue& nominal = context.nominal().criteria[k];
  const bool nominal_cancelled = context.nominal().status == TrialStatus::Cancelled;
  CriterionSensitivity out;
  out.id = criterion.id;
  out.kind = criterion_kind(criterion.operand);
  out.nominal = nominal.value;
  if (!nominal.value) out.nominal_reason = nominal.reason;
  const std::size_t n = limits.size();
  for (std::size_t i = 0; i < n; ++i) {
    out.items.push_back(contribution(out.kind, nominal.value, nominal_cancelled,
                                     slots.trials[2 * i], slots.trials[2 * i + 1], k));
    if (!out.items.back().excluded.empty()) out.excluded.push_back(i);
  }
  out.complete = out.excluded.empty();
  if (!nominal.value) {
    if (with_inverse && !nominal_cancelled) out.inverse.note = "nominal value not defined";
    return out;
  }
  if (slots.cancelled) return out;
  const double phi0 = *nominal.value;

  // RSS and prediction over the items that are not excluded, in item order (point 8).
  if (out.kind == CriterionKind::Linear) {
    double sum = 0.0;
    double predicted = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
      const ItemContribution& c = out.items[i];
      if (!c.delta) continue;
      sum += *c.delta * *c.delta;
      const double s = ratios[i] * *c.delta;
      predicted += s * s;
    }
    out.rss = std::sqrt(sum);
    out.predicted_sigma = std::sqrt(predicted);
  } else {
    std::vector<double> q;
    double predicted = phi0 * phi0;
    for (std::size_t i = 0; i < n; ++i) {
      const ItemContribution& c = out.items[i];
      if (!c.q) continue;
      q.push_back(*c.q);
      predicted += ratios[i] * ratios[i] * *c.q;
    }
    const detail::QuadraticRss rss = detail::quadratic_rss(phi0, q);
    out.rss = rss.value;
    if (rss.negative) warnings.push_back(rss_warning(criterion, k, stage));
    out.predicted_mean_square = predicted;
  }
  if (!with_inverse) return out;

  // Inverse sensitivity: A distributed equally on all n items (addendum #258), T'_i per item.
  out.inverse.limits.assign(n, std::nullopt);
  const std::optional<Distribution> distribution =
      distribute(criterion, out.kind, phi0, n, out.inverse);
  if (!distribution) return out;
  const double share = distribution->share;
  const double a = distribution->allowance;
  constexpr double eps = optim::kMeritPrecision;
  for (std::size_t i = 0; i < n; ++i) {
    const ItemContribution& c = out.items[i];
    if (!c.excluded.empty() || !c.plus || !c.minus) continue;
    const double p = *c.plus;
    const double m = *c.minus;
    if (c.delta) {
      const double scale = std::max({std::abs(phi0), std::abs(p), std::abs(m), a});
      const double d = std::abs(*c.delta);
      if (d <= eps * scale) {
        out.inverse.unlimited.push_back(i);
      } else {
        out.inverse.limits[i] = limits[i] * share / d;
      }
    } else if (c.q) {
      const double q = *c.q;
      const double scale = std::max({phi0 * phi0, p * p, m * m, a * a});
      if (q <= eps * scale) {
        out.inverse.unlimited.push_back(i);
      } else {
        out.inverse.limits[i] = limits[i] * std::sqrt(share / q);
      }
    }
  }
  return out;
}

}  // namespace

CriterionKind criterion_kind(const model::Operand& operand) noexcept {
  return std::holds_alternative<model::SpotRmsOperand>(operand) ||
                 std::holds_alternative<model::OpdRmsOperand>(operand)
             ? CriterionKind::Quadratic
             : CriterionKind::Linear;
}

double sigma_ratio(model::ToleranceDistribution distribution) {
  switch (distribution) {
    case model::ToleranceDistribution::Normal: {
      // Variance of the normal truncated at +-a sigma: (1 - 2 a phi(a) / Z) sigma^2, derived in
      // ADR 0036, point 5 (docs/quellen.md, "Varianz der bei +-a sigma abgeschnittenen
      // Normalverteilung"), with Z = erf(a / sqrt(2)) the probability of |x| <= a sigma (DLMF,
      // Eq. (7.2.1)); a = 2 and sigma = T / 2 (point 5).
      const double a = 2.0;
      const double phi = std::exp(-a * a / 2.0) / std::sqrt(2.0 * std::numbers::pi);
      const double z = std::erf(a / std::numbers::sqrt2);
      return std::sqrt(1.0 - 2.0 * a * phi / z) / 2.0;
    }
    case model::ToleranceDistribution::Uniform:
      return 1.0 / std::sqrt(3.0);  // NIST e-Handbook, Sec. 1.3.6.6.2
    case model::ToleranceDistribution::Endpoint:
      return 1.0;
  }
  throw std::invalid_argument("sigma_ratio: unknown distribution");
}

ToleranceSensitivity sensitivity(const TrialContext& context,
                                 const SensitivityOptions& options,
                                 const trace::RunControl* control) {
  const trace::RunControl none;
  const trace::RunControl& run = control != nullptr ? *control : none;
  const std::span<const model::ToleranceItem> items = context.items();
  const std::size_t n = items.size();
  const std::size_t criteria = context.criteria().size();

  ToleranceSensitivity out;
  out.nominal = context.nominal();
  std::vector<double> limits(n);
  std::vector<double> ratios(n);
  for (std::size_t i = 0; i < n; ++i) {
    limits[i] = items[i].limit;
    ratios[i] = sigma_ratio(items[i].distribution);
    out.items.push_back({items[i], ratios[i], std::nullopt, std::nullopt});
  }

  // Phi_0 cancelled in prepare (#257): nothing to evaluate against.
  Slots slots;
  if (out.nominal.status == TrialStatus::Cancelled) {
    slots.trials.resize(2 * n);
    slots.cancelled = true;
  } else {
    slots = run_trials(context, limits, run, "sensitivity");
  }
  for (std::size_t i = 0; i < n; ++i) {
    out.items[i].plus = slots.trials[2 * i];
    out.items[i].minus = slots.trials[2 * i + 1];
  }
  for (std::size_t k = 0; k < criteria; ++k) {
    out.criteria.push_back(criterion_sensitivity(context, k, limits, ratios, slots, true,
                                                 "sensitivity", out.warnings));
  }
  out.limits.assign(n, std::nullopt);
  if (slots.cancelled) {
    out.status = SensitivityStatus::Cancelled;
    return out;
  }

  // Per item the smallest T'_i over the criteria that limit it (addendum #258, F3).
  for (const CriterionSensitivity& c : out.criteria) {
    for (std::size_t i = 0; i < c.inverse.limits.size(); ++i) {
      const std::optional<double>& t = c.inverse.limits[i];
      if (!t) continue;
      const double limit = *t;
      std::optional<double>& smallest = out.limits[i];
      smallest = smallest ? std::min(*smallest, limit) : limit;
    }
  }
  if (!options.confirm) return out;

  // Confirmation at +-T'_i, T_i for an item that no criterion limits (F3). A trial that is not
  // valid there excludes the item for that criterion as in the first part (F1; ADR 0036,
  // addendum #257: a radius of the other sign).
  Confirmation confirmation;
  for (std::size_t i = 0; i < n; ++i)
    confirmation.limits.push_back(out.limits[i].value_or(limits[i]));
  const Slots second = run_trials(context, confirmation.limits, run, "confirmation");
  if (second.cancelled) {
    out.status = SensitivityStatus::Cancelled;
    return out;
  }
  for (std::size_t i = 0; i < n; ++i) {
    const std::optional<TrialResult>& plus = second.trials[2 * i];
    const std::optional<TrialResult>& minus = second.trials[2 * i + 1];
    if (!plus || !minus) throw std::logic_error("sensitivity: a confirmation slot is empty");
    confirmation.plus.push_back(*plus);
    confirmation.minus.push_back(*minus);
  }
  for (std::size_t k = 0; k < criteria; ++k) {
    confirmation.criteria.push_back(criterion_sensitivity(
        context, k, confirmation.limits, ratios, second, false, "confirmation", out.warnings));
  }
  out.confirmation = std::move(confirmation);
  return out;
}

ToleranceSensitivity sensitivity(const model::System& system,
                                 const material::MaterialLibrary& materials,
                                 const coating::CoatingLibrary* coatings,
                                 const SensitivityOptions& options,
                                 const trace::RunControl* control) {
  ToleranceItems expansion = expand(system, materials);
  const TrialContext context =
      prepare(system, materials, coatings, expansion.items, options.trial, control);
  ToleranceSensitivity out = sensitivity(context, options, control);
  out.expansion = std::move(expansion);
  return out;
}

}  // namespace rtt::tolerance
