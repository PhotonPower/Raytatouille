// Intermediate state per accepted iteration for a live display (#378; ADR 0030, addendum #378).
// The reference is the run itself: the callback must report exactly the accepted entries of the
// history, its last patch must give the result system, and the run must not change with it.

#include <tbb/global_control.h>
#include <tbb/info.h>

#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#include "rtt/io/edit.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/model/optimization.hpp"
#include "rtt/optim/levenberg_marquardt.hpp"
#include "rtt/optim/optimize.hpp"
#include "rtt/trace/run_control.hpp"

using rtt::material::MaterialLibrary;
using rtt::model::Param;
using rtt::model::System;
using rtt::optim::IterationState;
using rtt::optim::LmStatus;
using rtt::optim::OptimIteration;
using rtt::optim::optimize;
using rtt::optim::OptimResult;

namespace {

System load(const std::string& relative) {
  return rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/" + relative);
}

/// m5/singlet_optim: three variables, an EFL operand and an RMS spot generator.
System singlet_optim() {
  return load("m5/singlet_optim.rtt.json");
}

/// m1/singlet_const with two configurations, a variable table row bound to R1 and a variable
/// detector position; EFL 100 and 80 in the two configurations (as test M6 in
/// test_optimize.cpp).
System zoom_singlet() {
  System s = load("m1/singlet_const.rtt.json");
  s.configurations = {{"long"}, {"short"}};
  rtt::model::ParameterRow r;
  r.name = "R";
  r.form = std::vector<double>{60.0, 60.0};
  r.variable = true;
  s.parameters = {r};
  auto& lens = std::get<rtt::model::Element>(s.root.children[1].value);
  std::get<rtt::model::Conic>(lens.surfaces[0].shape.base).radius = Param::bound("R");
  std::get<rtt::model::Element>(s.root.children[2].value).pose.position[2].variable = true;
  const auto efl = [](double target, const char* configuration) {
    rtt::model::FirstOrderOperand op;
    op.common = {target, 1.0, std::string(configuration)};
    op.quantity = rtt::model::FirstOrderQuantity::Efl;
    op.path = "main";
    return op;
  };
  s.optimization.operands = {efl(100.0, "long"), efl(80.0, "short")};
  return s;
}

std::vector<OptimIteration> accepted(const OptimResult& r) {
  std::vector<OptimIteration> out;
  for (const OptimIteration& it : r.history) {
    if (it.accepted) out.push_back(it);
  }
  return out;
}

}  // namespace

TEST_CASE("#378: one call per accepted solve, phi bitwise the history", "[optim][optimize][378]") {
  const MaterialLibrary lib;
  for (const char* which : {"m5/singlet_optim", "zoom_singlet"}) {
    INFO(which);
    const System s = std::string(which) == "zoom_singlet" ? zoom_singlet() : singlet_optim();
    std::vector<IterationState> calls;
    const OptimResult r = optimize(s, lib, nullptr, {}, {},
                                   [&calls](const IterationState& st) { calls.push_back(st); });
    const std::vector<OptimIteration> acc = accepted(r);
    REQUIRE(acc.size() >= 2);  // content guard: the run takes several steps
    REQUIRE(calls.size() == acc.size());
    for (std::size_t i = 0; i < calls.size(); ++i) {
      CHECK(calls[i].k == acc[i].k);
      CHECK(calls[i].phi == acc[i].phi);  // bitwise
      CHECK(calls[i].values.size() == r.variables.size());
    }
    // phi falls from call to call: every call is an accepted step (rho > 0).
    for (std::size_t i = 1; i < calls.size(); ++i) CHECK(calls[i].phi < calls[i - 1].phi);
  }
}

TEST_CASE("#378: the patch of a call gives that state; the last one gives the result system",
          "[optim][optimize][378]") {
  const MaterialLibrary lib;
  for (const char* which : {"m5/singlet_optim", "zoom_singlet"}) {
    INFO(which);
    const System s = std::string(which) == "zoom_singlet" ? zoom_singlet() : singlet_optim();
    std::vector<IterationState> calls;
    const OptimResult r = optimize(s, lib, nullptr, {}, {},
                                   [&calls](const IterationState& st) { calls.push_back(st); });
    REQUIRE(!calls.empty());
    // The last accepted state is the result state: same values, same patch.
    CHECK(calls.back().patch == r.patch);
    for (std::size_t j = 0; j < r.variables.size(); ++j) {
      CHECK(calls.back().values[j] == r.variables[j].end);  // bitwise
    }
    const rtt::io::PatchResult last = rtt::io::apply_patch(s, calls.back().patch);
    CHECK(rtt::io::to_json(last.system) == rtt::io::to_json(r.system));
    // Every intermediate patch applies to the input (a valid system) and differs from the next.
    for (std::size_t i = 0; i < calls.size(); ++i) {
      INFO("call " << i);
      const rtt::io::PatchResult applied = rtt::io::apply_patch(s, calls[i].patch);
      CHECK(applied.system.optimization == s.optimization);
      if (i + 1 < calls.size()) CHECK(calls[i].patch != calls[i + 1].patch);
    }
  }
}

TEST_CASE("#378: results are bitwise the same with and without the callback, any thread count",
          "[optim][optimize][378]") {
  const System s = zoom_singlet();
  const MaterialLibrary lib;
  const OptimResult plain = [&] {
    const tbb::global_control limit(tbb::global_control::max_allowed_parallelism, 1);
    return optimize(s, lib);
  }();
  REQUIRE(plain.patch != "[]");  // content guard
  for (const std::size_t threads : {std::size_t{1}, std::size_t{4},
                                    static_cast<std::size_t>(tbb::info::default_concurrency())}) {
    INFO("threads " << threads);
    const tbb::global_control limit(tbb::global_control::max_allowed_parallelism, threads);
    std::vector<IterationState> calls;
    const OptimResult r = optimize(s, lib, nullptr, {}, {},
                                   [&calls](const IterationState& st) { calls.push_back(st); });
    CHECK(r.status == plain.status);
    CHECK(r.patch == plain.patch);
    CHECK(r.evaluations == plain.evaluations);
    REQUIRE(r.history.size() == plain.history.size());
    for (std::size_t k = 0; k < r.history.size(); ++k) {
      CHECK(r.history[k].phi == plain.history[k].phi);
      CHECK(r.history[k].accepted == plain.history[k].accepted);
    }
    CHECK(calls.size() == accepted(plain).size());
  }
}

TEST_CASE("#378: an exception in the callback ends the run like a cancellation",
          "[optim][optimize][378]") {
  const System s = zoom_singlet();
  const MaterialLibrary lib;
  std::vector<IterationState> calls;
  OptimResult r;
  REQUIRE_NOTHROW(r = optimize(s, lib, nullptr, {}, {}, [&calls](const IterationState& st) {
                    calls.push_back(st);
                    if (calls.size() == 2) throw std::runtime_error("stop the display");
                  }));
  CHECK(r.status == LmStatus::Cancelled);
  REQUIRE(calls.size() == 2);
  // The last accepted state is the one of the throwing call (ADR 0030, addendum #378).
  CHECK(r.patch == calls[1].patch);
  for (std::size_t j = 0; j < r.variables.size(); ++j) {
    CHECK(r.variables[j].end == calls[1].values[j]);
  }
  // The history ends with that solve; no final evaluation after a cancellation (addendum #167).
  REQUIRE(!r.history.empty());
  CHECK(r.history.back().k == calls[1].k);
  CHECK(r.history.back().accepted);
  CHECK(accepted(r).size() == 2);
  CHECK(r.operands.empty());
  // The same state as a run limited to that solve (ADR 0030, point 8, step 7).
  rtt::optim::OptimizeOptions k_max;
  k_max.max_iterations = calls[1].k;
  const OptimResult limited = optimize(s, lib, nullptr, k_max);
  CHECK(limited.patch == r.patch);
}

TEST_CASE("#378: no call without an accepted step", "[optim][optimize][378]") {
  const System s = zoom_singlet();
  const MaterialLibrary lib;
  rtt::optim::OptimizeOptions none;
  none.max_iterations = 0;
  int calls = 0;
  const OptimResult r =
      optimize(s, lib, nullptr, none, {}, [&calls](const IterationState&) { ++calls; });
  CHECK(r.status == LmStatus::MaxIterations);
  CHECK(calls == 0);
}

TEST_CASE("#378: levenberg_marquardt reports every accepted state", "[optim][lm][378]") {
  // f(p) = (p0 - 3, 2 (p1 + 1)): minimum at (3, -1); the solver takes accepted steps towards it.
  const rtt::optim::Residuals f = [](std::span<const double> p, std::span<double> out) {
    out[0] = p[0] - 3.0;
    out[1] = 2.0 * (p[1] + 1.0);
    return true;
  };
  std::vector<std::vector<double>> states;
  std::vector<int> ks;
  const rtt::optim::LmResult r = rtt::optim::levenberg_marquardt(
      f, 2, {0.0, 0.0}, {rtt::optim::Bounds{}, rtt::optim::Bounds{}}, {}, {},
      [&](const rtt::optim::LmIteration& it, std::span<const double> p,
          std::span<const std::uint8_t> changed) {
        // CHECK, not REQUIRE: an exception from here would end the run as Cancelled (#378).
        CHECK(it.accepted);
        CHECK(changed.size() == 2);
        ks.push_back(it.k);
        states.emplace_back(p.begin(), p.end());
      });
  std::size_t n_accepted = 0;
  for (const auto& it : r.history) n_accepted += it.accepted ? 1 : 0;
  REQUIRE(n_accepted >= 1);  // content guard
  CHECK(states.size() == n_accepted);
  CHECK(states.back() == r.p);  // bitwise
}
