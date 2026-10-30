// Tolerance cases of rtt_py_reference (#261): the same calls as CASES in
// test_bitwise_tolerance.py, on the M7 reference achromat with the SCHOTT catalogue: expand of
// the full set, and sensitivity of its three explicit items (class none).

#include "tolerance_cases.hpp"

#include <oneapi/tbb/task_arena.h>

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "npy_writer.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/system.hpp"
#include "rtt/tolerance/expand.hpp"
#include "rtt/tolerance/monte_carlo.hpp"
#include "rtt/tolerance/sensitivity.hpp"
#include "rtt/tolerance/trial.hpp"

namespace rtt::py::reference {
namespace {

namespace fs = std::filesystem;

/// Writes the arrays of one case as <out>/<case>.<name>.npy.
class Writer {
 public:
  Writer(fs::path out, std::string name) : out_(std::move(out)), name_(std::move(name)) {}

  void f8(const std::string& array, const std::vector<double>& v) const {
    write_column<double>(file(array), "<f8", v);
  }
  void u8(const std::string& array, const std::vector<std::uint64_t>& v) const {
    write_column<std::uint64_t>(file(array), "<u8", v);
  }
  /// Lines joined with '\n' as UTF-8 bytes.
  void text(const std::string& array, const std::vector<std::string>& lines) const {
    std::string joined;
    for (std::size_t i = 0; i < lines.size(); ++i) {
      if (i > 0) joined += '\n';
      joined += lines[i];
    }
    write_column<std::uint8_t>(file(array), "|u1",
                               std::vector<std::uint8_t>(joined.begin(), joined.end()));
  }

 private:
  [[nodiscard]] fs::path file(const std::string& array) const {
    return out_ / (name_ + "." + array + ".npy");
  }

  fs::path out_;
  std::string name_;
};

template <typename E>
std::uint64_t code(E e) {
  return static_cast<std::uint64_t>(e);
}

/// 0 without, else the index + 1.
std::uint64_t optional_index(const std::optional<std::size_t>& i) {
  return i ? *i + 1 : 0;
}

// Flattening; the names match flatten_items() in test_bitwise_tolerance.py.

void write_items(const Writer& w, const tolerance::ToleranceItems& r) {
  std::vector<std::string> ids;
  std::vector<double> limit;
  std::vector<std::uint64_t> kind;
  std::vector<std::uint64_t> distribution;
  for (const model::ToleranceItem& i : r.items) {
    ids.push_back(i.id);
    limit.push_back(i.limit);
    kind.push_back(code(i.kind));
    distribution.push_back(code(i.distribution));
  }
  std::vector<std::uint64_t> origin;
  for (const tolerance::ItemOrigin& o : r.origins) {
    origin.push_back(code(o.source));
    origin.push_back(o.index);
  }
  std::vector<std::string> class_ids;
  std::vector<std::uint64_t> class_ints;
  for (const model::ClassTolerance& c : r.classes.items) {
    class_ids.push_back(c.id);
    class_ints.insert(class_ints.end(),
                      {code(c.kind), code(c.level), code(c.tolerance_class), code(c.skip),
                       optional_index(c.value_override), optional_index(c.replaced_by)});
  }
  std::vector<std::uint64_t> overrides;
  for (const model::OverrideStatus& v : r.classes.value_overrides) {
    overrides.insert(overrides.end(), {v.index, code(v.effect), code(v.skip)});
  }
  w.text("ids", ids);
  w.f8("limit", limit);
  w.u8("kind", kind);
  w.u8("distribution", distribution);
  w.u8("origin", origin);
  w.text("class_ids", class_ids);
  w.u8("class_ints", class_ints);
  w.u8("overrides", overrides);
}

/// NaN for none (the same bits as math.nan in Python).
double value(const std::optional<double>& x) {
  return x.value_or(std::numeric_limits<double>::quiet_NaN());
}

/// Criterion values, compensator ends, mu and the iterations of one trial (NaN for none).
void add_trial(std::vector<double>& v, const std::optional<tolerance::TrialResult>& t) {
  if (!t) return;
  for (const tolerance::CriterionValue& c : t->criteria) v.push_back(value(c.value));
  for (const tolerance::CompensatorValue& c : t->compensators) v.push_back(c.end);
  v.push_back(value(t->mu));
  v.push_back(static_cast<double>(t->lm_iterations));
}

/// Per criterion: Phi_0, rss, the predictions, allowance and share, then per item the six
/// contributions and the inverse limit (NaN for none).
std::vector<double> criterion_values(const std::vector<tolerance::CriterionSensitivity>& cs) {
  std::vector<double> v;
  for (const tolerance::CriterionSensitivity& c : cs) {
    for (const auto* x : {&c.nominal, &c.rss, &c.predicted_sigma, &c.predicted_mean_square,
                          &c.inverse.allowance, &c.inverse.share}) {
      v.push_back(value(*x));
    }
    for (const tolerance::ItemContribution& i : c.items) {
      for (const auto* x : {&i.plus, &i.minus, &i.worst, &i.delta, &i.shift, &i.q}) {
        v.push_back(value(*x));
      }
    }
    for (const std::optional<double>& x : c.inverse.limits) v.push_back(value(x));
  }
  return v;
}

/// Per criterion: kind, complete, the count of excluded and unlimited items, then their indices.
std::vector<std::uint64_t> criterion_ints(const std::vector<tolerance::CriterionSensitivity>& cs) {
  std::vector<std::uint64_t> v;
  for (const tolerance::CriterionSensitivity& c : cs) {
    v.insert(v.end(),
             {code(c.kind), c.complete ? 1U : 0U, c.excluded.size(), c.inverse.unlimited.size()});
    v.insert(v.end(), c.excluded.begin(), c.excluded.end());
    v.insert(v.end(), c.inverse.unlimited.begin(), c.inverse.unlimited.end());
  }
  return v;
}

void write_sensitivity(const Writer& w, const tolerance::ToleranceSensitivity& r) {
  std::vector<std::string> ids;
  std::vector<double> sigma_ratio;
  std::vector<double> trials;
  add_trial(trials, r.nominal);
  for (const tolerance::ItemTrials& i : r.items) {
    ids.push_back(i.item.id);
    sigma_ratio.push_back(i.sigma_ratio);
    add_trial(trials, i.plus);
    add_trial(trials, i.minus);
  }
  std::vector<double> limits;
  limits.reserve(r.limits.size());
  for (const std::optional<double>& x : r.limits) limits.push_back(value(x));
  std::vector<std::uint64_t> status{code(r.status), code(r.nominal.status),
                                    r.confirmation ? 1U : 0U};
  for (const tolerance::ItemTrials& i : r.items) {
    for (const auto* t : {&i.plus, &i.minus}) status.push_back(*t ? code((*t)->status) + 1 : 0);
  }
  w.text("ids", ids);
  w.f8("sigma_ratio", sigma_ratio);
  w.f8("trials", trials);
  w.f8("criteria", criterion_values(r.criteria));
  w.u8("criterion_ints", criterion_ints(r.criteria));
  w.f8("limits", limits);
  w.u8("status", status);
  std::vector<double> confirmation;
  if (r.confirmation) {
    confirmation = r.confirmation->limits;
    for (std::size_t i = 0; i < r.confirmation->plus.size(); ++i) {
      add_trial(confirmation, r.confirmation->plus[i]);
      add_trial(confirmation, r.confirmation->minus[i]);
    }
    const std::vector<double> c = criterion_values(r.confirmation->criteria);
    confirmation.insert(confirmation.end(), c.begin(), c.end());
  }
  w.f8("confirmation", confirmation);
}

/// good, total, fraction, lower, upper of a yield (NaN for none).
void add_yield(std::vector<double>& v, const tolerance::Yield& y) {
  v.insert(v.end(), {static_cast<double>(y.good), static_cast<double>(y.total), value(y.fraction),
                     value(y.lower), value(y.upper)});
}

void write_monte_carlo(const Writer& w, const tolerance::MonteCarloResult& r) {
  std::vector<std::uint64_t> ints{r.seed, r.requested, code(r.status), r.trials.size()};
  std::vector<double> statistics;
  add_yield(statistics, r.yield);
  for (const tolerance::CriterionStatistics& c : r.criteria) {
    ints.insert(ints.end(), {c.with_value, c.without_value, c.count, c.percentiles.size()});
    for (const auto* x : {&c.mean, &c.standard_deviation, &c.min, &c.max}) {
      statistics.push_back(value(*x));
    }
    for (const tolerance::Percentile& p : c.percentiles) {
      ints.push_back(p.rank);
      statistics.insert(statistics.end(), {p.p, p.value});
    }
    add_yield(statistics, c.yield);
  }
  // The columns of TrialTable, row by row. Flags: lm_status + 1 with 0 for no run, satisfied
  // as 2/1/0 for yes/no/undecidable (as flatten_monte_carlo in test_bitwise_tolerance.py).
  std::vector<double> rows;
  std::vector<std::uint64_t> flags;
  for (const tolerance::TrialRecord& t : r.trials) {
    flags.insert(flags.end(), {t.trial, code(t.result.status), t.good ? 1U : 0U,
                               t.result.lm_status ? code(*t.result.lm_status) + 1 : 0});
    rows.insert(rows.end(), t.values.begin(), t.values.end());
    for (const tolerance::CriterionValue& c : t.result.criteria) {
      rows.push_back(value(c.value));
      flags.push_back(c.satisfied ? (*c.satisfied ? 2U : 1U) : 0U);
    }
    for (const tolerance::CompensatorValue& c : t.result.compensators) {
      rows.insert(rows.end(), {c.start, c.end});
      flags.insert(flags.end(), {c.at_bound ? 1U : 0U, c.infeasible ? 1U : 0U});
    }
    rows.push_back(value(t.result.mu));
    rows.push_back(static_cast<double>(t.result.lm_iterations));
  }
  w.u8("ints", ints);
  w.f8("statistics", statistics);
  w.f8("rows", rows);
  w.u8("flags", flags);
}

/// The achromat with the class switched off: only its three explicit items (as
/// test_tolerance_sensitivity.py, explicit_only()).
model::System explicit_only(model::System s) {
  s.tolerances.tolerance_class = model::ToleranceClass::None;
  auto& overrides = s.tolerances.overrides;
  overrides.erase(overrides.begin() + 4);  // image none
  overrides.erase(overrides.begin());      // objective precise
  return s;
}

}  // namespace

void run_tolerance_cases(const fs::path& reference_dir,
                         const fs::path& catalog_dir,
                         const fs::path& out,
                         int threads) {
  material::MaterialLibrary schott;
  schott.add_catalog(catalog_dir / "schott.agf");
  const model::System achromat = io::load_system(reference_dir / "m7/achromat_tolerances.rtt.json");
  oneapi::tbb::task_arena arena(threads);
  arena.execute(
      [&] { write_items(Writer(out, "achromat_expand"), tolerance::expand(achromat, schott)); });
  arena.execute([&] {
    write_sensitivity(Writer(out, "achromat_sensitivity"),
                      tolerance::sensitivity(explicit_only(achromat), schott, nullptr));
  });
  arena.execute([&] {
    tolerance::MonteCarloOptions options;
    options.trials = 12;  // the seed of the file
    write_monte_carlo(Writer(out, "achromat_monte_carlo"),
                      tolerance::monte_carlo(explicit_only(achromat), schott, nullptr, options));
  });
}

}  // namespace rtt::py::reference
