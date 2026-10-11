// A variable focal length of an ideal lens that reaches f = 0 during optimization (#178, review
// of #217 by f2): compile rejects f = 0 (interaction.focal_length_invalid), so that evaluation is
// invalid (ADR 0030, point 10) and the run ends without an exception.

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <optional>
#include <string>
#include <variant>

#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
#include "rtt/optim/optimize.hpp"

using rtt::model::Element;
using rtt::model::IdealLens;
using rtt::model::Param;
using rtt::model::System;
using rtt::optim::LmStatus;
using rtt::optim::OptimResult;

namespace {

bool has(const OptimResult& r, const std::string& code, const std::string& location) {
  for (const rtt::model::Diagnostic& d : r.diagnostics) {
    if (d.code == code && d.location == location) return true;
  }
  return false;
}

}  // namespace

TEST_CASE("ideal lens: an evaluation at f = 0 is invalid, the run ends without exception",
          "[optim][ideal]") {
  // Start f = h mm, unbounded (theta = p): the Jacobian steps by h = cbrt(eps_f) max(|theta|, 1)
  // (levenberg_marquardt.cpp); for theta = h < 1 that is h itself, so theta - h = 0 exactly,
  // provided the test has the same h as the solver to the last bit. Both therefore call cbrt at
  // run time on the same value (volatile: no folding by the compiler, which rounds correctly while
  // a C library need not; GCC 13 -O3 folded cbrt(kMeritPrecision) in CI). No exactness of cbrt is
  // assumed: glibc's cbrt is not correctly rounded, even for an exact cube such as 2^-30.
  // The operand is a ray of a path that bypasses the lenses (STO -> IMG): every evaluation
  // compiles the whole system, the lens included, so only f = 0 fails.
  System s = rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/r2/ideal_lens.rtt.json");
  rtt::optim::OptimizeOptions options;
  volatile double eps_f = rtt::optim::kMeritPrecision;  // not a constant for the compiler
  options.function_precision = eps_f;
  const double h = std::cbrt(static_cast<double>(eps_f));
  // Precondition, checked loudly: theta - h_solver = 0 with the solver's formula.
  REQUIRE(h < 1.0);
  REQUIRE(h - h * std::max(std::abs(h), 1.0) == 0.0);
  Param f(h);
  f.variable = true;
  std::get<IdealLens>(std::get<Element>(s.root.children[1].value).surfaces[0].interaction)
      .focal_length = f;
  using rtt::model::EventKind;
  using rtt::model::SurfaceId;
  s.paths.push_back(
      {"free",
       false,
       {{SurfaceId{"STO"}, EventKind::Transmit, 0}, {SurfaceId{"IMG"}, EventKind::Transmit, 0}}});
  rtt::model::RayOperand op;
  op.path = "free";
  op.surface = SurfaceId{"IMG"};
  op.py = 1.0;
  op.common.target = 0.0;
  s.optimization.operands = {op};

  OptimResult r;
  REQUIRE_NOTHROW(r = rtt::optim::optimize(s, rtt::material::MaterialLibrary{}, nullptr, options));
  CHECK(r.status == LmStatus::Failed);
  CHECK(r.failed_evaluations >= 1);
  CHECK(has(r, "optim.evaluation_failed", "/optimization"));
  CHECK(has(r, "optim.jacobian_failed",
            "/root/children/1/surfaces/0/interaction/focal_length/value"));
  // The result is the last accepted state, the input.
  CHECK(r.system == s);
  CHECK(r.patch == "[]");
}
