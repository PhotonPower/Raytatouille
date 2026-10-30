// M7 acceptance cases of rtt_py_reference (#262): the same calls as CASES in test_bitwise_m7.py,
// on the reference systems of test_m7_acceptance.cpp (tests/reference/m7). Python binds only the
// overloads with a System, so the item subsets of the C++ acceptance are made in the system here
// (keep_items) and with apply_patch there.

#include "m7_cases.hpp"

#include <oneapi/tbb/task_arena.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <initializer_list>
#include <limits>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "npy_writer.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/system.hpp"
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

/// NaN for none (the same bits as math.nan in Python).
double value(const std::optional<double>& x) {
  return x.value_or(std::numeric_limits<double>::quiet_NaN());
}

// Flattening; the names match flatten_sensitivity() and flatten_monte_carlo() in
// test_bitwise_m7.py.

/// Per criterion: Phi_0, rss, predicted sigma and mean square, then per item delta, q, Phi(+T),
/// Phi(-T) (NaN for none); the status of the run and of each trial.
void write_sensitivity(const Writer& w, const tolerance::ToleranceSensitivity& r) {
  std::vector<double> v;
  for (const tolerance::CriterionSensitivity& c : r.criteria) {
    for (const auto* x : {&c.nominal, &c.rss, &c.predicted_sigma, &c.predicted_mean_square}) {
      v.push_back(value(*x));
    }
    for (const tolerance::ItemContribution& i : c.items) {
      for (const auto* x : {&i.delta, &i.q, &i.plus, &i.minus}) v.push_back(value(*x));
    }
  }
  std::vector<std::uint64_t> status{code(r.status), code(r.nominal.status)};
  for (const tolerance::ItemTrials& i : r.items) {
    for (const auto* t : {&i.plus, &i.minus}) status.push_back(*t ? code((*t)->status) + 1 : 0);
  }
  w.f8("criteria", v);
  w.u8("status", status);
}

/// The yield as fraction, lower, upper (NaN for none) and good, total.
void add_yield(std::vector<double>& f, std::vector<std::uint64_t>& u, const tolerance::Yield& y) {
  f.insert(f.end(), {value(y.fraction), value(y.lower), value(y.upper)});
  u.insert(u.end(), {y.good, y.total});
}

/// Per trial (ascending): number, status, good, the perturbations, the criterion values (NaN for
/// none) and the compensator ends; per criterion the statistics and the yield; the total yield.
void write_monte_carlo(const Writer& w, const tolerance::MonteCarloResult& r) {
  std::vector<std::uint64_t> trial_ints;
  std::vector<double> values;
  std::vector<double> criteria;
  std::vector<double> compensators;
  for (const tolerance::TrialRecord& t : r.trials) {
    trial_ints.insert(trial_ints.end(), {t.trial, code(t.result.status), t.good ? 1U : 0U});
    values.insert(values.end(), t.values.begin(), t.values.end());
    for (const tolerance::CriterionValue& c : t.result.criteria) criteria.push_back(value(c.value));
    for (const tolerance::CompensatorValue& c : t.result.compensators) {
      compensators.push_back(c.end);
    }
  }
  std::vector<double> statistics;
  std::vector<std::uint64_t> counts;
  for (const tolerance::CriterionStatistics& c : r.criteria) {
    for (const auto* x : {&c.mean, &c.standard_deviation, &c.min, &c.max}) {
      statistics.push_back(value(*x));
    }
    for (const tolerance::Percentile& p : c.percentiles) {
      statistics.push_back(p.value);
      counts.push_back(p.rank);
    }
    counts.insert(counts.end(), {c.with_value, c.without_value, c.count});
    add_yield(statistics, counts, c.yield);
  }
  add_yield(statistics, counts, r.yield);
  counts.insert(counts.end(), {code(r.status), r.seed, r.requested});
  w.u8("trial_ints", trial_ints);
  w.f8("values", values);
  w.f8("criteria", criteria);
  w.f8("compensators", compensators);
  w.f8("statistics", statistics);
  w.u8("counts", counts);
}

/// The system with only the explicit items `ids`, in file order (as keep_items() in
/// test_bitwise_m7.py).
model::System keep_items(model::System s, std::initializer_list<const char*> ids) {
  auto& items = s.tolerances.items;
  std::erase_if(items, [&](const model::ToleranceItem& i) {
    return std::none_of(ids.begin(), ids.end(), [&](const char* id) { return i.id == id; });
  });
  return s;
}

}  // namespace

void run_m7_cases(const fs::path& reference_dir,
                  const fs::path& catalog_dir,
                  const fs::path& out,
                  int threads) {
  material::MaterialLibrary schott;
  schott.add_catalog(catalog_dir / "schott.agf");
  const material::MaterialLibrary no_glass;
  const auto load = [&](const char* name) { return io::load_system(reference_dir / "m7" / name); };
  tolerance::SensitivityOptions sensitivity;
  sensitivity.confirm = false;
  oneapi::tbb::task_arena arena(threads);
  arena.execute([&] {
    // Case 6: decenters of lens and detector at m = -0.5.
    write_sensitivity(Writer(out, "m7_sensitivity_lens"),
                      tolerance::sensitivity(keep_items(load("decentred_lens_finite.rtt.json"),
                                                        {"lens-dy", "image-dy"}),
                                             no_glass, nullptr, sensitivity));
    // Case 4: index and Abbe number of N-BK7, efl at F, d and C.
    write_sensitivity(
        Writer(out, "m7_sensitivity_bk7"),
        tolerance::sensitivity(load("planoconvex_bk7.rtt.json"), schott, nullptr, sensitivity));
    // Case 8: the uniform decenter, seed and trials of the file.
    write_monte_carlo(
        Writer(out, "m7_mc_lens"),
        tolerance::monte_carlo(load("decentred_lens.rtt.json"), no_glass, nullptr, {}));
    // Case 7: LM in every trial, 20 trials (the bitwise path, not the statistics).
    tolerance::MonteCarloOptions twenty;
    twenty.trials = 20;
    write_monte_carlo(
        Writer(out, "m7_mc_compensator"),
        tolerance::monte_carlo(load("paraboloid_compensator.rtt.json"), no_glass, nullptr, twenty));
    // Case 13: the items of class, overrides and explicit items, 10 trials.
    tolerance::MonteCarloOptions ten;
    ten.trials = 10;
    write_monte_carlo(
        Writer(out, "m7_mc_classes"),
        tolerance::monte_carlo(load("achromat_tolerances.rtt.json"), schott, nullptr, ten));
  });
}

}  // namespace rtt::py::reference
