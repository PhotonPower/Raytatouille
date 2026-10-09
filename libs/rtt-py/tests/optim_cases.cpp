// Optimization cases of rtt_py_reference (#169): the same calls as CASES in
// test_bitwise_optim.py, on the M5 reference systems with the merit functions of the M5
// acceptance (#170) in their section "optimization".

#include "optim_cases.hpp"

#include <oneapi/tbb/task_arena.h>

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <iterator>
#include <string>
#include <utility>
#include <vector>

#include "npy_writer.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/system.hpp"
#include "rtt/optim/merit.hpp"
#include "rtt/optim/optimize.hpp"

namespace rtt::py::reference {
namespace {

namespace fs = std::filesystem;
using namespace rtt::optim;

/// Writes the arrays of one case as <out>/<case>.<name>.npy.
class Writer {
 public:
  Writer(fs::path out, std::string name) : out_(std::move(out)), name_(std::move(name)) {}

  template <typename T, typename Items, typename Get>
  void column(const std::string& array, const char* descr, const Items& items, Get get) const {
    std::vector<T> v;
    v.reserve(std::size(items));
    for (const auto& item : items) v.push_back(static_cast<T>(get(item)));
    write_column<T>(out_ / (name_ + "." + array + ".npy"), descr, v);
  }
  void f8(const std::string& array, const std::vector<double>& v) const {
    column<double>(array, "<f8", v, [](double x) { return x; });
  }
  void u8(const std::string& array, const std::vector<std::uint64_t>& v) const {
    column<std::uint64_t>(array, "<u8", v, [](std::uint64_t x) { return x; });
  }

 private:
  fs::path out_;
  std::string name_;
};

// Flattening; the names match flatten_result() and flatten_evaluation() in
// test_bitwise_optim.py.

void write_result(const Writer& w, const OptimResult& r) {
  const auto& h = r.history;
  w.column<std::int32_t>("k", "<i4", h, [](const OptimIteration& i) { return i.k; });
  w.column<double>("phi", "<f8", h, [](const OptimIteration& i) { return i.phi; });
  w.column<double>("mu", "<f8", h, [](const OptimIteration& i) { return i.mu; });
  w.column<double>("rho", "<f8", h, [](const OptimIteration& i) { return i.rho; });
  w.column<double>("step_norm", "<f8", h, [](const OptimIteration& i) { return i.step_norm; });
  w.column<std::uint8_t>("accepted", "|b1", h,
                         [](const OptimIteration& i) { return i.accepted ? 1 : 0; });
  w.column<std::uint64_t>("evaluations", "<u8", h,
                          [](const OptimIteration& i) { return i.evaluations; });
  const auto& v = r.variables;
  w.column<double>("var_start", "<f8", v, [](const VariableValue& x) { return x.start; });
  w.column<double>("var_end", "<f8", v, [](const VariableValue& x) { return x.end; });
  w.column<std::uint8_t>("var_changed", "|b1", v,
                         [](const VariableValue& x) { return x.changed ? 1 : 0; });
  w.column<std::uint8_t>("var_at_bound", "|b1", v,
                         [](const VariableValue& x) { return x.at_bound ? 1 : 0; });
  w.column<double>("op_value", "<f8", r.operands, [](const OperandValue& x) { return x.value; });
  w.column<double>("op_contribution", "<f8", r.operands,
                   [](const OperandValue& x) { return x.contribution; });
  w.column<double>("gen_rms", "<f8", r.generators, [](const GeneratorValue& x) { return x.rms; });
  w.column<double>("gen_contribution", "<f8", r.generators,
                   [](const GeneratorValue& x) { return x.contribution; });
  std::vector<std::uint64_t> rays;
  for (const GeneratorValue& g : r.generators) {
    rays.push_back(g.rays_launched);
    rays.push_back(g.rays_lost);
  }
  w.u8("gen_rays", rays);
  w.u8("ints", {static_cast<std::uint64_t>(r.status), static_cast<std::uint64_t>(r.iterations),
                r.evaluations, r.failed_evaluations, r.diagnostics.size()});
  w.column<std::uint8_t>("patch", "|u1", r.patch,
                         [](char c) { return static_cast<std::uint8_t>(c); });
}

void write_evaluation(const Writer& w, const MeritEvaluation& e) {
  w.f8("values", e.values);
  w.f8("residuals", e.residuals);
  w.column<double>("gen_mean_square", "<f8", e.generators,
                   [](const GeneratorStats& g) { return g.mean_square; });
  std::vector<std::uint64_t> rays;
  for (const GeneratorStats& g : e.generators) {
    rays.push_back(g.rays_launched);
    rays.push_back(g.rays_lost);
  }
  w.u8("gen_rays", rays);
  w.u8("ints", {e.valid() ? 1U : 0U, e.warnings.size()});
}

}  // namespace

void run_optim_cases(const fs::path& reference_dir, const fs::path& out, int threads) {
  const auto w = [&](const std::string& name) { return Writer(out, name); };
  const material::MaterialLibrary plain;
  const model::System singlet = io::load_system(reference_dir / "m5/singlet_optim.rtt.json");
  const model::System gap = io::load_system(reference_dir / "m5/two_lens_gap.rtt.json");
  const MeritFunction merit(singlet, plain, nullptr);
  const std::vector<double> start = merit.start();

  oneapi::tbb::task_arena arena(threads);
  arena.execute([&] {
    write_result(w("singlet_optimize"), optimize(singlet, plain));
    write_result(w("gap_optimize"), optimize(gap, plain));
    write_evaluation(w("singlet_merit_start"), merit.evaluate(start));
  });
}

}  // namespace rtt::py::reference
