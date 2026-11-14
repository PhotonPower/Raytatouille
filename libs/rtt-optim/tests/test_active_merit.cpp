// Merit function with nodes taken out of the computation (#380, ADR 0049, points 2 and 4): ray
// operands at inactive surfaces are left out in their place, and computations on a path without
// events throw at the API boundary instead of reading past the events.

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <string_view>
#include <variant>
#include <vector>

#include "rtt/analysis/opd.hpp"
#include "rtt/analysis/spot.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/compile/errors.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/model/optimization.hpp"
#include "rtt/model/parameters.hpp"
#include "rtt/optim/generators.hpp"
#include "rtt/optim/merit.hpp"
#include "rtt/optim/optimize.hpp"
#include "rtt/paraxial/paraxial.hpp"
#include "rtt/trace/sources.hpp"

using Catch::Matchers::ContainsSubstring;
using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::material::MaterialLibrary;
using rtt::model::CircularAperture;
using rtt::model::Element;
using rtt::model::ElementKind;
using rtt::model::Node;
using rtt::model::Param;
using rtt::model::Pose;
using rtt::model::Surface;
using rtt::model::SurfaceId;
using rtt::model::System;
using rtt::optim::MeritEvaluation;
using rtt::optim::MeritFunction;
using rtt::optim::OptimError;

namespace {

/// m1/singlet_const: stop STO (child 0), lens L1 in CONST:1.5168 (child 1), detector IMG
/// (child 2), automatic path "main"; every node placed absolutely.
System singlet() {
  System s = rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m1/singlet_const.rtt.json");
  // The front radius of the lens is the variable of every case.
  bool found = false;
  rtt::model::for_each_param(s, [&](std::string_view pointer, Param& p) {
    if (!found && pointer == "/root/children/1/surfaces/0/shape/base/radius") {
      p.variable = true;
      found = true;
    }
  });
  REQUIRE(found);
  return s;
}

/// The singlet with an inactive plane-parallel window (CONST:1.5, 2 mm) at z = 1 mm, between stop
/// and lens; the lens moves to child 2.
System with_inactive_window() {
  System s = singlet();
  Surface a;
  a.id = SurfaceId("W.S1");
  a.aperture = CircularAperture{12.7, 0.0};
  Surface b = a;
  b.id = SurfaceId("W.S2");
  b.pose = Pose::along_z(2.0);
  Element window{"window", ElementKind::Plate, Pose::along_z(1.0), "CONST:1.5", {a, b}};
  window.active = false;
  s.root.children.insert(s.root.children.begin() + 1, Node{window});
  return s;
}

rtt::model::RayOperand ray_y(const std::string& surface, double py) {
  rtt::model::RayOperand op;
  op.common = {0.0, 1.0, std::nullopt};
  op.coordinate = rtt::model::RayCoordinate::Y;
  op.path = "main";
  op.surface = SurfaceId(surface);
  op.py = py;
  return op;
}

rtt::model::FirstOrderOperand efl(double target) {
  rtt::model::FirstOrderOperand op;
  op.common = {target, 1.0, std::nullopt};
  op.quantity = rtt::model::FirstOrderQuantity::Efl;
  op.path = "main";
  return op;
}

bool has_code(const OptimError& e, const std::string& code) {
  const auto& d = e.diagnostics();
  return std::any_of(d.begin(), d.end(), [&](const auto& x) { return x.code == code; });
}

}  // namespace

TEST_CASE("active: a ray operand at an inactive surface is left out in its place", "[active]") {
  const MaterialLibrary lib;
  System a = with_inactive_window();
  a.optimization.operands = {ray_y("W.S1", 1.0), efl(100.0), ray_y("IMG", 1.0)};
  // B: the same system without the window and without the left-out operand.
  System b = singlet();
  b.optimization.operands = {efl(100.0), ray_y("IMG", 1.0)};

  const MeritFunction ma(a, lib, nullptr);
  const MeritFunction mb(b, lib, nullptr);
  CHECK(ma.inactive_operands() == std::vector<bool>{true, false, false});
  CHECK(mb.inactive_operands() == std::vector<bool>{false, false});
  REQUIRE(ma.size() == 3);
  const MeritEvaluation ea = ma.evaluate(ma.start());
  const MeritEvaluation eb = mb.evaluate(mb.start());
  // Left out: value NaN, residual 0, a reason.
  CHECK(std::isnan(ea.values[0]));
  CHECK(ea.residuals[0] == 0.0);
  CHECK(ea.undefined[0].find("inactive") != std::string::npos);
  CHECK(ea.valid());
  // The others as without the window, bit for bit (the paths agree, ADR 0049, point 2).
  REQUIRE(std::isfinite(eb.values[0]));
  REQUIRE(std::isfinite(eb.values[1]));
  CHECK(eb.values[1] != 0.0);  // content guard: the marginal ray at the image is off axis
  CHECK(ea.values[1] == eb.values[0]);
  CHECK(ea.values[2] == eb.values[1]);
  CHECK(ea.residuals[1] == eb.residuals[0]);
  CHECK(ea.residuals[2] == eb.residuals[1]);
}

TEST_CASE("active: only left-out operands and no generator is optim.no_operands", "[active]") {
  const MaterialLibrary lib;
  System s = with_inactive_window();
  s.optimization.operands = {ray_y("W.S1", 1.0), ray_y("W.S2", 0.5)};
  try {
    static_cast<void>(rtt::optim::optimize(s, lib));
    FAIL("no OptimError");
  } catch (const OptimError& e) {
    CHECK(has_code(e, "optim.no_operands"));
  }
  // With a generator left, the run is no input error of this kind.
  rtt::model::SpotGenerator g;
  g.path = "main";
  s.optimization.generators = {g};
  rtt::optim::OptimizeOptions options;
  options.max_iterations = 1;
  CHECK_NOTHROW(rtt::optim::optimize(s, lib, nullptr, options));
}

TEST_CASE("active: computations on a path without events throw", "[active]") {
  // ADR 0049, point 2: all elements inactive (paths.all_inactive). Every computation on the
  // path throws at the API boundary, never reads past its events.
  const MaterialLibrary lib;
  System s = singlet();
  for (Node& child : s.root.children) std::get<Element>(child.value).active = false;
  const CompiledSystem cs = rtt::compile::compile(s, lib);
  REQUIRE(cs.path(PathId{0}).events.empty());
  CHECK_THROWS_AS(rtt::paraxial::first_order(cs, PathId{0}, 0), rtt::paraxial::ParaxialError);
  CHECK_THROWS_AS(rtt::analysis::spot(cs, PathId{0}, 0, 0), std::invalid_argument);
  CHECK_THROWS_AS(rtt::analysis::opd_map(cs, PathId{0}, 0, 0), std::invalid_argument);
  CHECK_THROWS_AS(rtt::trace::aim_reference(cs, PathId{0}, 0), std::invalid_argument);

  // A generator evaluated directly (public generator_residuals): without the check it would
  // read the last event of the empty path.
  rtt::model::SpotGenerator g;
  g.path = "main";
  const rtt::model::Generator generator = g;
  std::vector<double> out(rtt::optim::generator_size(generator, cs));
  rtt::optim::GeneratorStats stats;
  CHECK_THROWS_WITH(rtt::optim::generator_residuals(generator, cs, out, stats),
                    ContainsSubstring("all its elements are inactive"));

  // The merit function: an operand on the path and a generator on it. The message is that of
  // the check, not merit.operand_unsupported from first_order (an OptimError).
  System op = s;
  op.optimization.operands = {efl(100.0)};
  CHECK_THROWS_WITH(MeritFunction(op, lib, nullptr),
                    ContainsSubstring("all its elements are inactive"));
  System gen = s;
  gen.optimization.generators = {g};
  CHECK_THROWS_WITH(MeritFunction(gen, lib, nullptr),
                    ContainsSubstring("all its elements are inactive"));
  CHECK_THROWS_AS(rtt::optim::optimize(gen, lib), std::invalid_argument);
}
