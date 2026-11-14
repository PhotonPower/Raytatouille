// Tolerancing with nodes taken out of the computation (#380, ADR 0049, point 5): no class items
// in inactive subtrees, explicit items, single values and compensators at inactive nodes are
// skipped with tolerances.target_inactive.

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <optional>
#include <string>
#include <variant>
#include <vector>

#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
#include "rtt/model/tolerance_classes.hpp"
#include "rtt/tolerance/expand.hpp"
#include "rtt/tolerance/monte_carlo.hpp"
#include "rtt/tolerance/trial.hpp"

using rtt::model::ClassSkip;
using rtt::model::Criterion;
using rtt::model::Diagnostic;
using rtt::model::Element;
using rtt::model::ElementKind;
using rtt::model::FieldType;
using rtt::model::Node;
using rtt::model::Param;
using rtt::model::Pose;
using rtt::model::Surface;
using rtt::model::SurfaceId;
using rtt::model::System;
using rtt::model::ToleranceAnchor;
using rtt::model::ToleranceItem;
using rtt::model::ToleranceKind;
using rtt::model::ToleranceTarget;

namespace {

Surface surface(const std::string& id, double z_mm = 0.0, std::optional<double> radius = {}) {
  Surface s;
  s.id = SurfaceId(id);
  s.pose = Pose::along_z(z_mm);
  s.aperture = rtt::model::CircularAperture{10.0, 0.0};
  if (radius) s.shape.base = rtt::model::Conic{Param(*radius), Param(0.0)};
  return s;
}

ToleranceTarget at(ToleranceAnchor anchor,
                   const std::string& name,
                   std::optional<std::string> param = {}) {
  ToleranceTarget t;
  t.anchor = anchor;
  t.name = name;
  t.param = std::move(param);
  return t;
}

/// A singlet in CONST:1.5168 (R = 50, -50, 5 mm) with a stop and a detector, object at infinity,
/// EPD 10, class normal, one criterion; every node placed absolutely.
System singlet() {
  System s;
  s.name = "singlet";
  s.wavelengths = {{0.5876, 1.0, true}};
  s.aperture = {rtt::model::SystemApertureType::EntrancePupilDiameter, Param(10.0)};
  s.fields = {FieldType::AngleDeg, {{0.0, 0.0, 1.0}}};
  s.root.name = "root";
  s.root.children = {
      {Element{"stop", ElementKind::Stop, Pose::along_z(-2.0), std::nullopt, {surface("STO")}}},
      {Element{"L",
               ElementKind::Lens,
               Pose::along_z(0.0),
               "CONST:1.5168",
               {surface("L.S1", 0.0, 50.0), surface("L.S2", 5.0, -50.0)}}},
      {Element{"D", ElementKind::Detector, Pose::along_z(50.0), std::nullopt, {surface("IMG")}}}};
  s.paths = {{"main", true, {}}};
  s.tolerances.tolerance_class = rtt::model::ToleranceClass::Normal;
  rtt::model::SpotRmsOperand spot;
  spot.path = "main";
  Criterion c;
  c.id = "spot";
  c.operand = spot;
  c.max = 1.0;
  s.tolerances.criteria = {c};
  return s;
}

/// The singlet with an inactive window W (plate, CONST:1.5, planes W.S1 and W.S2 at z = 20 and
/// 22 mm) between lens and detector.
System with_inactive_window() {
  System s = singlet();
  Element w{"W",
            ElementKind::Plate,
            Pose::along_z(20.0),
            "CONST:1.5",
            {surface("W.S1"), surface("W.S2", 2.0)}};
  w.active = false;
  s.root.children.insert(s.root.children.begin() + 2, Node{w});
  return s;
}

std::size_t count(const std::vector<Diagnostic>& d,
                  const std::string& code,
                  const std::string& location) {
  return static_cast<std::size_t>(std::count_if(d.begin(), d.end(), [&](const Diagnostic& x) {
    return x.code == code && x.location == location;
  }));
}

}  // namespace

TEST_CASE("inactive: no class items in an inactive subtree", "[tolerance][active]") {
  const rtt::model::ClassTolerances classes = rtt::model::class_tolerances(with_inactive_window());
  std::size_t at_window = 0;
  for (const rtt::model::ClassTolerance& c : classes.items) {
    const bool window = c.target.name == "W" || c.target.name.starts_with("W.");
    INFO(c.id);
    if (window) {
      ++at_window;
      CHECK(c.skip == ClassSkip::Inactive);
      CHECK_FALSE(c.generated());
    } else {
      CHECK(c.skip != ClassSkip::Inactive);
    }
  }
  CHECK(at_window > 0);
  // The generated items are those of the singlet without the window: ids, limits, order.
  const rtt::material::MaterialLibrary lib;
  const rtt::tolerance::ToleranceItems a = rtt::tolerance::expand(with_inactive_window(), lib);
  const rtt::tolerance::ToleranceItems b = rtt::tolerance::expand(singlet(), lib);
  REQUIRE_FALSE(b.items.empty());
  REQUIRE(a.items.size() == b.items.size());
  for (std::size_t i = 0; i < a.items.size(); ++i) {
    INFO(b.items[i].id);
    CHECK(a.items[i] == b.items[i]);
  }
  CHECK(a.warnings.empty());
}

TEST_CASE("inactive: an inactive assembly takes the class items of its subtree",
          "[tolerance][active]") {
  System s = singlet();
  rtt::model::Assembly group;
  group.name = "group";
  group.active = false;
  group.children.push_back(s.root.children[1]);
  s.root.children[1] = Node{group};
  const rtt::model::ClassTolerances classes = rtt::model::class_tolerances(s);
  std::size_t inactive = 0;
  for (const rtt::model::ClassTolerance& c : classes.items) {
    INFO(c.id);
    const bool below =
        c.target.name == "group" || c.target.name == "L" || c.target.name.starts_with("L.");
    CHECK((c.skip == ClassSkip::Inactive) == below);
    if (below) ++inactive;
  }
  CHECK(inactive > 0);
}

TEST_CASE("inactive: explicit items, single values and compensators are skipped with a warning",
          "[tolerance][active]") {
  System s = with_inactive_window();
  ToleranceItem on_window;
  on_window.id = "w-decenter";
  on_window.target = at(ToleranceAnchor::Node, "W");
  on_window.kind = ToleranceKind::DecenterX;
  on_window.limit = 0.1;
  ToleranceItem on_lens = on_window;
  on_lens.id = "l-decenter";
  on_lens.target = at(ToleranceAnchor::Node, "L");
  ToleranceItem on_window_surface;
  on_window_surface.id = "w-surface";
  on_window_surface.target = at(ToleranceAnchor::Surface, "W.S2");
  on_window_surface.kind = ToleranceKind::Thickness;
  on_window_surface.limit = 0.05;
  s.tolerances.items = {on_window, on_lens, on_window_surface};
  rtt::model::ToleranceValueOverride single;
  single.target = at(ToleranceAnchor::Node, "W");
  single.kind = ToleranceKind::DecenterY;
  single.limit = 0.2;
  rtt::model::ToleranceValueOverride single_lens = single;
  single_lens.target = at(ToleranceAnchor::Node, "L");
  s.tolerances.overrides = {single_lens, single};
  s.tolerances.compensators = {{at(ToleranceAnchor::Node, "W", "/pose/position/2"), -1.0, 1.0},
                               {at(ToleranceAnchor::Node, "D", "/pose/position/2"), -1.0, 1.0}};
  rtt::model::SpotGenerator focus;
  focus.path = "main";
  s.tolerances.compensator_merit = rtt::model::Optimization{{}, {focus}};
  REQUIRE_FALSE(rtt::model::has_errors(rtt::model::validate(s)));

  const rtt::material::MaterialLibrary lib;
  const rtt::tolerance::ToleranceItems x = rtt::tolerance::expand(s, lib);
  const auto explicit_ids = [&] {
    std::vector<std::string> ids;
    for (std::size_t i = 0; i < x.items.size(); ++i) {
      if (x.origins[i].source == rtt::tolerance::ItemOrigin::Source::Explicit) {
        ids.push_back(x.items[i].id);
      }
    }
    return ids;
  }();
  CHECK(explicit_ids == std::vector<std::string>{"l-decenter"});
  CHECK(count(x.warnings, "tolerances.target_inactive", "/tolerances/items/0") == 1);
  CHECK(count(x.warnings, "tolerances.target_inactive", "/tolerances/items/2") == 1);
  CHECK(count(x.warnings, "tolerances.target_inactive", "/tolerances/overrides/1") == 1);
  CHECK(count(x.warnings, "tolerances.target_inactive", "/tolerances/compensators/0") == 1);
  CHECK(x.warnings.size() == 4);
  for (const Diagnostic& d : x.warnings) CHECK(d.severity == rtt::model::Severity::Warning);
  // The single value at the lens applies; the one at the window is skipped with the class.
  REQUIRE(x.classes.value_overrides.size() == 2);
  CHECK(x.classes.value_overrides[0].effect == rtt::model::OverrideEffect::Applied);
  CHECK(x.classes.value_overrides[1].effect == rtt::model::OverrideEffect::Skipped);
  CHECK(x.classes.value_overrides[1].skip == ClassSkip::Inactive);
}

TEST_CASE("inactive: a compensator at an inactive node is not run", "[tolerance][active]") {
  // Only a compensator at the inactive window: it keeps its value and the trial runs no
  // optimization; a compensator at the detector next to it still runs.
  System s = with_inactive_window();
  rtt::model::SpotGenerator focus;
  focus.path = "main";
  s.tolerances.compensator_merit = rtt::model::Optimization{{}, {focus}};
  s.tolerances.compensators = {{at(ToleranceAnchor::Node, "W", "/pose/position/2"), -1.0, 1.0}};
  const rtt::material::MaterialLibrary lib;
  ToleranceItem shift;
  shift.id = "lens-z";
  shift.target = at(ToleranceAnchor::Node, "L");
  shift.kind = ToleranceKind::Thickness;
  shift.limit = 0.5;
  const std::vector<ToleranceItem> items{shift};
  const rtt::tolerance::TrialContext alone = rtt::tolerance::prepare(s, lib, nullptr, items);
  const std::vector<double> values{0.5};
  const rtt::tolerance::TrialResult r = rtt::tolerance::evaluate(alone, values);
  REQUIRE(r.compensators.size() == 1);
  CHECK(r.compensators[0].start == 20.0);
  CHECK(r.compensators[0].end == 20.0);
  CHECK(r.lm_iterations == 0);

  s.tolerances.compensators.push_back(
      {at(ToleranceAnchor::Node, "D", "/pose/position/2"), -2.0, 2.0});
  const rtt::tolerance::TrialContext both = rtt::tolerance::prepare(s, lib, nullptr, items);
  const rtt::tolerance::TrialResult q = rtt::tolerance::evaluate(both, values);
  REQUIRE(q.compensators.size() == 2);
  CHECK(q.compensators[0].end == 20.0);
  CHECK(q.compensators[1].end != q.compensators[1].start);  // the detector follows the lens
}

TEST_CASE("inactive: a target below an inactive assembly or reached through param",
          "[tolerance][active]") {
  // Review of PR 4: one rule with rtt-optim, the pointer the target acts on. The lens L sits in
  // an assembly G. (a) G inactive: an item and a compensator at L are inactive through the
  // ancestor. (b) G active, L inactive: a compensator at G whose param points into L is
  // inactive, one whose param points at G itself is not.
  System base = singlet();
  rtt::model::Assembly g;
  g.name = "G";
  g.children.push_back(base.root.children[1]);
  base.root.children[1] = Node{g};
  rtt::model::SpotGenerator focus;
  focus.path = "main";
  base.tolerances.compensator_merit = rtt::model::Optimization{{}, {focus}};
  const rtt::material::MaterialLibrary lib;

  System a = base;
  std::get<rtt::model::Assembly>(a.root.children[1].value).active = false;
  ToleranceItem on_lens;
  on_lens.id = "l-decenter";
  on_lens.target = at(ToleranceAnchor::Node, "L");
  on_lens.kind = ToleranceKind::DecenterX;
  on_lens.limit = 0.1;
  a.tolerances.items = {on_lens};
  a.tolerances.compensators = {{at(ToleranceAnchor::Node, "L", "/pose/position/2"), -1.0, 1.0}};
  const rtt::tolerance::ToleranceItems xa = rtt::tolerance::expand(a, lib);
  CHECK(count(xa.warnings, "tolerances.target_inactive", "/tolerances/items/0") == 1);
  CHECK(count(xa.warnings, "tolerances.target_inactive", "/tolerances/compensators/0") == 1);

  System b = base;
  std::get<Element>(std::get<rtt::model::Assembly>(b.root.children[1].value).children[0].value)
      .active = false;
  b.tolerances.compensators = {
      {at(ToleranceAnchor::Node, "G", "/children/0/pose/position/2"), -1.0, 1.0},
      {at(ToleranceAnchor::Node, "D", "/pose/position/2"), -2.0, 2.0}};
  REQUIRE_FALSE(rtt::model::has_errors(rtt::model::validate(b)));
  const rtt::tolerance::ToleranceItems xb = rtt::tolerance::expand(b, lib);
  CHECK(count(xb.warnings, "tolerances.target_inactive", "/tolerances/compensators/0") == 1);
  CHECK(count(xb.warnings, "tolerances.target_inactive", "/tolerances/compensators/1") == 0);
  // The trial does not try to run the inactive one: Phi_0 is valid.
  ToleranceItem detector_z;
  detector_z.id = "d-z";
  detector_z.target = at(ToleranceAnchor::Node, "D");
  detector_z.kind = ToleranceKind::Thickness;
  detector_z.limit = 0.5;
  const std::vector<ToleranceItem> items{detector_z};
  const rtt::tolerance::TrialContext ctx = rtt::tolerance::prepare(b, lib, nullptr, items);
  CHECK(ctx.nominal().status == rtt::tolerance::TrialStatus::Valid);
  REQUIRE(ctx.nominal().compensators.size() == 2);
  CHECK(ctx.nominal().compensators[0].start == ctx.nominal().compensators[0].end);
  CHECK_FALSE(ctx.nominal().compensators[0].infeasible);
  for (const Diagnostic& w : ctx.nominal().warnings) CHECK(w.code != "optim.variable_inactive");
}

TEST_CASE("inactive: a single value at an inactive node without class still warns",
          "[tolerance][active]") {
  System s = with_inactive_window();
  s.tolerances.overrides = {
      rtt::model::ToleranceClassOverride{"W", rtt::model::ToleranceClass::None},
      rtt::model::ToleranceValueOverride{at(ToleranceAnchor::Node, "W"), ToleranceKind::DecenterY,
                                         0.2, rtt::model::ToleranceDistribution::Normal}};
  const rtt::material::MaterialLibrary lib;
  const rtt::tolerance::ToleranceItems x = rtt::tolerance::expand(s, lib);
  REQUIRE(x.classes.value_overrides.size() == 1);
  CHECK(x.classes.value_overrides[0].effect == rtt::model::OverrideEffect::ClassNone);
  CHECK(count(x.warnings, "tolerances.target_inactive", "/tolerances/overrides/1") == 1);
}

TEST_CASE("inactive: a criterion at an inactive surface is left out, the trial stays valid",
          "[tolerance][active]") {
  // Decision of the maintainer in #380: as a merit operand, value none, not decided, not counted
  // for the yield; validate warns with merit.surface_inactive.
  System s = with_inactive_window();
  rtt::model::RayOperand at_window;
  at_window.path = "main";
  at_window.surface = SurfaceId("W.S1");
  at_window.py = 1.0;
  Criterion c;
  c.id = "at-window";
  c.operand = at_window;
  c.max = 0.0;  // would never be satisfied
  s.tolerances.criteria.push_back(c);
  const std::vector<Diagnostic> v = rtt::model::validate(s);
  CHECK_FALSE(rtt::model::has_errors(v));
  CHECK(count(v, "merit.surface_inactive", "/tolerances/criteria/1/operand/surface") == 1);

  const rtt::material::MaterialLibrary lib;
  rtt::tolerance::MonteCarloOptions options;
  options.seed = 7;
  options.trials = 8;
  const rtt::tolerance::MonteCarloResult r = rtt::tolerance::monte_carlo(s, lib, nullptr, options);
  CHECK(r.nominal.status == rtt::tolerance::TrialStatus::Valid);
  REQUIRE(r.nominal.criteria.size() == 2);
  CHECK(r.nominal.criteria[1].inactive);
  CHECK_FALSE(r.nominal.criteria[1].value.has_value());
  CHECK_FALSE(r.nominal.criteria[1].satisfied.has_value());
  CHECK_FALSE(r.nominal.criteria[0].inactive);
  REQUIRE(r.trials.size() == 8);
  for (const rtt::tolerance::TrialRecord& t : r.trials) {
    CHECK(t.result.status == rtt::tolerance::TrialStatus::Valid);
    CHECK(t.result.criteria[1].inactive);
    // Good exactly if the active criterion is satisfied.
    CHECK(t.good == t.result.criteria[0].satisfied.value_or(false));
  }
  CHECK(r.yield.good > 0);
  REQUIRE(r.criteria.size() == 2);
  CHECK(r.criteria[1].count == 0);
  CHECK(r.criteria[1].note.find("inactive") != std::string::npos);
}

TEST_CASE("inactive: Monte Carlo draws and criteria are those of the system without the window",
          "[tolerance][active]") {
  // ADR 0049, addendum PR 4 (corrects point 5): the draws depend only on seed, trial and id, and
  // with absolute placement (point 2) the inactive window changes nothing.
  const rtt::material::MaterialLibrary lib;
  rtt::tolerance::MonteCarloOptions options;
  options.seed = 11;
  options.trials = 6;
  const rtt::tolerance::MonteCarloResult a =
      rtt::tolerance::monte_carlo(with_inactive_window(), lib, nullptr, options);
  const rtt::tolerance::MonteCarloResult b =
      rtt::tolerance::monte_carlo(singlet(), lib, nullptr, options);
  REQUIRE(a.items == b.items);
  REQUIRE(a.trials.size() == b.trials.size());
  for (std::size_t n = 0; n < a.trials.size(); ++n) {
    INFO("trial " << n);
    CHECK(a.trials[n].values == b.trials[n].values);
    REQUIRE(a.trials[n].result.criteria.size() == 1);
    CHECK(a.trials[n].result.criteria[0].value == b.trials[n].result.criteria[0].value);
  }
  // Content guard: the trials differ from each other and have values.
  REQUIRE(a.trials[0].result.criteria[0].value.has_value());
  CHECK(a.trials[0].values != a.trials[1].values);
}
