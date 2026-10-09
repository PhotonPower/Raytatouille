// A variable focal length of an ideal lens that reaches f = 0 during optimization (#178, review
// of #217 by f2): compile rejects f = 0 (interaction.focal_length_invalid), so that evaluation is
// invalid (ADR 0030, point 10) and the run ends without an exception.

#include <catch2/catch_test_macros.hpp>
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
  // Start f = 1e-3 mm, unbounded (theta = p): the Jacobian steps by h = cbrt(eps_f) max(|theta|,
  // 1) = 1e-3 (kMeritPrecision = 1e-9), so theta - h = 0 exactly. The merit is the value of a
  // table row, so every other evaluation is valid and only f = 0 fails.
  System s = rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/r2/ideal_lens.rtt.json");
  Param f(1e-3);
  f.variable = true;
  std::get<IdealLens>(std::get<Element>(s.root.children[1].value).surfaces[0].interaction)
      .focal_length = f;
  rtt::model::ParameterRow row;
  row.name = "X";
  row.form = 1.0;
  s.parameters = {row};
  rtt::model::ParamValueOperand op;
  op.parameter = "X";
  op.common.target = 0.0;
  s.optimization.operands = {op};

  OptimResult r;
  REQUIRE_NOTHROW(r = rtt::optim::optimize(s, rtt::material::MaterialLibrary{}));
  CHECK(r.status == LmStatus::Failed);
  CHECK(r.failed_evaluations >= 1);
  CHECK(has(r, "optim.evaluation_failed", "/optimization"));
  CHECK(has(r, "optim.jacobian_failed",
            "/root/children/1/surfaces/0/interaction/focal_length/value"));
  // The result is the last accepted state, the input.
  CHECK(r.system == s);
  CHECK(r.patch == "[]");
}
