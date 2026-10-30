// A base radius that a trial turns to zero, to infinity or to the other sign, whatever changed it
// (#293; ADR 0036, addendum #293): a row of the parameter table (also through an expression row,
// also in one configuration only) or a compensator. The radius items themselves are tested in
// test_trial.cpp (addendum #257).
//
// The system is a singlet of N-BK7 (R1 = 50 mm, R2 = -50 mm) as in test_trial.cpp. The values of
// the invalid cases put the radius clearly on the other side of 0, those of the content guards
// clearly on the same side, so that a check of the wrong sign or of the unresolved table fails one
// of them.

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
#include "rtt/model/tolerance.hpp"
#include "rtt/tolerance/monte_carlo.hpp"
#include "rtt/tolerance/sensitivity.hpp"
#include "rtt/tolerance/trial.hpp"

using rtt::model::Criterion;
using rtt::model::Element;
using rtt::model::ElementKind;
using rtt::model::FieldType;
using rtt::model::Param;
using rtt::model::ParameterExpression;
using rtt::model::ParameterRow;
using rtt::model::Pose;
using rtt::model::Surface;
using rtt::model::SurfaceId;
using rtt::model::System;
using rtt::model::ToleranceAnchor;
using rtt::model::ToleranceDistribution;
using rtt::model::ToleranceItem;
using rtt::model::ToleranceKind;
using rtt::model::ToleranceTarget;
using rtt::tolerance::TrialContext;
using rtt::tolerance::TrialResult;
using rtt::tolerance::TrialStatus;

namespace {

const rtt::material::MaterialLibrary& library() {
  struct Schott {  // MaterialLibrary is not movable
    rtt::material::MaterialLibrary lib;
    Schott() { lib.add_catalog(std::string(RTT_CATALOG_DIR) + "/schott.agf"); }
  };
  static const Schott s;
  return s.lib;
}

Surface surface(const std::string& id, double z_mm, std::optional<Param> radius = {}) {
  Surface s;
  s.id = SurfaceId(id);
  s.pose = Pose::along_z(z_mm);
  if (radius) s.shape.base = rtt::model::Conic{*radius, Param(0.0)};
  return s;
}

/// The singlet of test_trial.cpp with the radii `r1` and `r2` (plain or bound to a row): stop at
/// z = -2, lens "L" (L.S1 at 0, L.S2 at 5 mm, N-BK7), detector at 50; object at infinity, EPD 10;
/// one criterion, the RMS spot of field 0.
System singlet(Param r1, Param r2) {
  System s;
  s.name = "singlet";
  s.wavelengths = {{0.5876, 1.0, true}};
  s.aperture = {rtt::model::SystemApertureType::EntrancePupilDiameter, Param(10.0)};
  s.fields = {FieldType::AngleDeg, {{0.0, 0.0, 1.0}}};
  s.root.name = "root";
  Surface stop = surface("STO", 0.0);
  stop.aperture = rtt::model::CircularAperture{8.0, 0.0};
  s.root.children = {
      {Element{"stop", ElementKind::Stop, Pose::along_z(-2.0), std::nullopt, {stop}}},
      {Element{"L",
               ElementKind::Lens,
               Pose::along_z(0.0),
               "SCHOTT:N-BK7",
               {surface("L.S1", 0.0, std::move(r1)), surface("L.S2", 5.0, std::move(r2))}}},
      {Element{
          "D", ElementKind::Detector, Pose::along_z(50.0), std::nullopt, {surface("IMG", 0.0)}}}};
  s.paths = {{"main", true, {}}};
  rtt::model::SpotRmsOperand spot;
  spot.path = "main";
  Criterion c;
  c.id = "spot";
  c.operand = spot;
  c.max = 1.0;
  s.tolerances.criteria = {c};
  return s;
}

ParameterRow row(const std::string& name, rtt::model::ParameterForm form) {
  ParameterRow r;
  r.name = name;
  r.form = std::move(form);
  return r;
}

ToleranceItem parameter_item(const std::string& id,
                             const std::string& row_name,
                             double limit,
                             ToleranceDistribution distribution = ToleranceDistribution::Normal) {
  ToleranceItem i;
  i.id = id;
  i.target.anchor = ToleranceAnchor::Parameter;
  i.target.name = row_name;
  i.kind = ToleranceKind::Parameter;
  i.limit = limit;
  i.distribution = distribution;
  return i;
}

TrialContext context(const System& s, const std::vector<ToleranceItem>& items) {
  return rtt::tolerance::prepare(s, library(), nullptr, items);
}

bool contains(const std::string& text, const std::string& part) {
  return text.find(part) != std::string::npos;
}

/// An Invalid trial without criterion values, with a reason naming `surface` and `row`.
void check_invalid(const TrialResult& r, const std::string& surface, const std::string& row) {
  INFO("reason: " << r.reason);
  CHECK(r.status == TrialStatus::Invalid);
  CHECK(contains(r.reason, "surface '" + surface + "'"));
  CHECK(contains(r.reason, "row '" + row + "'"));
  CHECK(contains(r.reason, "other sign"));
  REQUIRE(r.criteria.size() == 1);
  CHECK(!r.criteria[0].value.has_value());
  CHECK(r.criteria[0].reason == r.reason);
}

void check_valid(const TrialResult& r) {
  INFO("reason: " << r.reason);
  CHECK(r.status == TrialStatus::Valid);
  REQUIRE(r.criteria.size() == 1);
  CHECK(r.criteria[0].value.has_value());
}

}  // namespace

TEST_CASE("radius sign: a row of the parameter table", "[tolerance][trial][radius]") {
  // L.S1 bound to the row R1 = 50. Offset -70: R' = -20 (other sign); offset -20: R' = 30.
  System s = singlet(Param::bound("R1"), Param(-50.0));
  s.parameters = {row("R1", 50.0)};
  const TrialContext ctx = context(s, {parameter_item("r1", "R1", 100.0)});
  REQUIRE(ctx.nominal().status == TrialStatus::Valid);
  check_invalid(rtt::tolerance::evaluate(ctx, std::vector<double>{-70.0}), "L.S1", "R1");
  check_valid(rtt::tolerance::evaluate(ctx, std::vector<double>{-20.0}));
}

TEST_CASE("radius sign: a row through an expression", "[tolerance][trial][radius]") {
  // L.S2 bound to the expression row R2 = A - 70 with A = 20, so R2 = -50; L.S1 stays 50. Offset
  // +80 on A: R2 = 30 (other sign); offset +10: R2 = -40. Only the expression changes sign, A
  // stays positive in both.
  System s = singlet(Param(50.0), Param::bound("R2"));
  s.parameters = {row("A", 20.0), row("R2", ParameterExpression{"A - 70"})};
  const TrialContext ctx = context(s, {parameter_item("a", "A", 100.0)});
  REQUIRE(ctx.nominal().status == TrialStatus::Valid);
  check_invalid(rtt::tolerance::evaluate(ctx, std::vector<double>{80.0}), "L.S2", "R2");
  check_valid(rtt::tolerance::evaluate(ctx, std::vector<double>{10.0}));
}

TEST_CASE("radius sign: a row that changes sign in one configuration only",
          "[tolerance][trial][radius]") {
  // R1 = 50 (near) and 20 (far). Offset -30: 20 and -10, so only "far" changes sign; offset -10:
  // 40 and 10.
  System s = singlet(Param::bound("R1"), Param(-50.0));
  s.configurations = {{"near"}, {"far"}};
  s.parameters = {row("R1", std::vector<double>{50.0, 20.0})};
  const TrialContext ctx = context(s, {parameter_item("r1", "R1", 100.0)});
  REQUIRE(ctx.nominal().status == TrialStatus::Valid);
  const TrialResult r = rtt::tolerance::evaluate(ctx, std::vector<double>{-30.0});
  check_invalid(r, "L.S1", "R1");
  CHECK(contains(r.reason, "configuration 'far'"));
  check_valid(rtt::tolerance::evaluate(ctx, std::vector<double>{-10.0}));
}

TEST_CASE("radius sign: a compensator that ends on the other side of 0",
          "[tolerance][trial][radius]") {
  // The compensator is the row R1 (L.S1, nominal 50), travel +-200; its merit is the value of the
  // row with a target. Target -40: LM ends near -40, on the other side of 0, so already the
  // nominal trial Phi_0 (all items at 0) is Invalid after the compensator run, and so is every
  // trial; the compensator values stay in the result. Target 30: valid. The item only gives the
  // context a perturbation; it does not touch R1.
  const auto with_target = [](double target) {
    System s = singlet(Param::bound("R1"), Param(-50.0));
    s.parameters = {row("R1", 50.0)};
    ToleranceTarget t;
    t.anchor = ToleranceAnchor::Parameter;
    t.name = "R1";
    s.tolerances.compensators = {{t, -200.0, 200.0}};
    rtt::model::ParamValueOperand value;
    value.parameter = "R1";
    value.common.target = target;
    s.tolerances.compensator_merit = rtt::model::Optimization{{value}, {}};
    ToleranceItem shift;
    shift.id = "D.dz";
    shift.target.anchor = ToleranceAnchor::Node;
    shift.target.name = "D";
    shift.kind = ToleranceKind::ShiftZ;
    shift.limit = 0.1;
    return context(s, {shift});
  };
  const TrialContext flipped = with_target(-40.0);
  const TrialResult r = rtt::tolerance::evaluate(flipped, std::vector<double>{0.05});
  check_invalid(r, "L.S1", "R1");
  CHECK(contains(r.reason, "after the compensator run"));
  REQUIRE(r.compensators.size() == 1);
  CHECK(r.compensators[0].start == 50.0);
  CHECK(r.compensators[0].end < -10.0);  // the run went across 0
  CHECK(r.lm_status.has_value());
  CHECK(flipped.nominal().status == TrialStatus::Invalid);
  CHECK(contains(flipped.nominal().reason, "after the compensator run"));

  const TrialContext kept = with_target(30.0);
  const TrialResult ok = rtt::tolerance::evaluate(kept, std::vector<double>{0.05});
  check_valid(ok);
  REQUIRE(ok.compensators.size() == 1);
  CHECK(ok.compensators[0].end > 10.0);  // content guard: the compensator moved, same side
}

TEST_CASE("radius sign: sensitivity and Monte Carlo exclude the flipped trials",
          "[tolerance][trial][radius]") {
  // R1 = 50 with a parameter item of limit 60: +T gives 110, -T gives -10 (other sign). The
  // sensitivity excludes the item for the criterion with that reason (addendum #258); in the
  // Monte Carlo (uniform, 100 trials, P(offset < -50) = 10/120) the flipped trials are Invalid
  // with the reason and stay out of the statistics (addendum #259, D3 and D5).
  System s = singlet(Param::bound("R1"), Param(-50.0));
  s.parameters = {row("R1", 50.0)};
  const TrialContext ctx =
      context(s, {parameter_item("r1", "R1", 60.0, ToleranceDistribution::Uniform)});

  rtt::tolerance::SensitivityOptions so;
  so.confirm = false;
  const rtt::tolerance::ToleranceSensitivity sens = rtt::tolerance::sensitivity(ctx, so);
  REQUIRE(sens.criteria.size() == 1);
  const rtt::tolerance::CriterionSensitivity& c = sens.criteria[0];
  CHECK(!c.complete);
  CHECK(c.excluded == std::vector<std::size_t>{0});
  REQUIRE(c.items.size() == 1);
  INFO("excluded: " << c.items[0].excluded);
  CHECK(contains(c.items[0].excluded, "row 'R1'"));
  REQUIRE(sens.items.size() == 1);
  REQUIRE(sens.items[0].plus.has_value());
  REQUIRE(sens.items[0].minus.has_value());
  CHECK(sens.items[0].plus->status == TrialStatus::Valid);
  CHECK(sens.items[0].minus->status == TrialStatus::Invalid);

  rtt::tolerance::MonteCarloOptions mo;
  mo.trials = 100;
  mo.seed = 293;
  const rtt::tolerance::MonteCarloResult mc = rtt::tolerance::monte_carlo(ctx, mo);
  REQUIRE(mc.trials.size() == 100);
  // A radius between 0 and the ray height 5 mm keeps its sign but loses rays (another reason);
  // those trials are only counted, not classified.
  std::size_t flipped = 0;
  std::size_t valid = 0;
  for (const rtt::tolerance::TrialRecord& t : mc.trials) {
    const double r1 = 50.0 + t.values[0];
    if (r1 <= 0.0) {
      ++flipped;
      CHECK(t.result.status == TrialStatus::Invalid);
      CHECK(contains(t.result.reason, "row 'R1'"));
    } else {
      CHECK(!contains(t.result.reason, "other sign"));
      if (r1 >= 10.0) CHECK(t.result.status == TrialStatus::Valid);
    }
    if (t.result.status == TrialStatus::Valid) ++valid;
  }
  CHECK(flipped > 0);  // content guard: some trials changed the sign
  CHECK(valid > 0);
  REQUIRE(mc.criteria.size() == 1);
  CHECK(mc.criteria[0].count == valid);  // Invalid trials stay out of the statistics
}
