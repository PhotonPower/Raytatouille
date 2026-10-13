// One tolerance trial (ADR 0036, points 1 to 4; #257): prepare() and evaluate(). The values are
// applied as offsets on Params and rows (1a) and as the rigid and material overlay (1b, 1c), the
// compensators run with Levenberg-Marquardt (point 3), the criteria are evaluated with the
// nominal reference (point 2).

#include <oneapi/tbb/info.h>
#include <oneapi/tbb/parallel_for.h>
#include <oneapi/tbb/task_arena.h>

#include <bit>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numbers>
#include <optional>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
#include "rtt/optim/bounds.hpp"
#include "rtt/optim/merit.hpp"
#include "rtt/tolerance/trial.hpp"

using rtt::model::Compensator;
using rtt::model::Criterion;
using rtt::model::Element;
using rtt::model::ElementKind;
using rtt::model::FieldType;
using rtt::model::Param;
using rtt::model::Pose;
using rtt::model::Surface;
using rtt::model::SurfaceId;
using rtt::model::System;
using rtt::model::ToleranceAnchor;
using rtt::model::ToleranceItem;
using rtt::model::ToleranceKind;
using rtt::model::ToleranceTarget;
using rtt::tolerance::TrialContext;
using rtt::tolerance::TrialResult;
using rtt::tolerance::TrialStatus;

namespace {

bool same(double a, double b) {
  return std::bit_cast<std::uint64_t>(a) == std::bit_cast<std::uint64_t>(b);
}

bool same(const std::optional<double>& a, const std::optional<double>& b) {
  return a.has_value() == b.has_value() && (!a || same(*a, *b));
}

void require_same(const TrialResult& a, const TrialResult& b) {
  CHECK(a.status == b.status);
  CHECK(a.reason == b.reason);
  REQUIRE(a.criteria.size() == b.criteria.size());
  for (std::size_t k = 0; k < a.criteria.size(); ++k) {
    CHECK(same(a.criteria[k].value, b.criteria[k].value));
    CHECK(a.criteria[k].satisfied == b.criteria[k].satisfied);
  }
  REQUIRE(a.compensators.size() == b.compensators.size());
  for (std::size_t c = 0; c < a.compensators.size(); ++c) {
    CHECK(same(a.compensators[c].start, b.compensators[c].start));
    CHECK(same(a.compensators[c].end, b.compensators[c].end));
  }
  CHECK(a.lm_status == b.lm_status);
  CHECK(same(a.mu, b.mu));
  CHECK(a.lm_iterations == b.lm_iterations);
  CHECK(a.warnings.size() == b.warnings.size());
}

/// The library of the tests: SCHOTT for N-BK7 (not movable, so built in place once).
const rtt::material::MaterialLibrary& library() {
  static rtt::material::MaterialLibrary materials;
  static const bool loaded = [] {
    materials.add_catalog(std::string(RTT_CATALOG_DIR) + "/schott.agf");
    return true;
  }();
  static_cast<void>(loaded);
  return materials;
}

Surface plane(const std::string& id, double z_mm = 0.0, std::optional<double> radius = {}) {
  Surface s;
  s.id = SurfaceId(id);
  s.pose = Pose::along_z(z_mm);
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

ToleranceItem item(const std::string& id, ToleranceTarget target, ToleranceKind kind, double t) {
  ToleranceItem i;
  i.id = id;
  i.target = std::move(target);
  i.kind = kind;
  i.limit = t;
  return i;
}

/// Ideal lens f = 100 in vacuum (ADR 0031): stop r = 10 at z = -5, lens "lens" at z = 0 (aperture
/// 30), detector "D" at z = 100 (the rear focal plane); object at infinity, EPD 10, fields 0 and
/// 1 degree. Compensator: the detector position z, travel +-2 mm, against an rms_spot generator
/// of field 0 about the chief ray. Criteria: spot_rms of field 0 about the chief ray ("spot"),
/// ray_y of the chief ray of field 1 on the detector ("chief", max_change 0.05 mm).
System ideal_lens() {
  System s;
  s.name = "ideal lens with a detector compensator";
  s.environment.medium = "VACUUM";
  s.wavelengths = {{0.5876, 1.0, true}};
  s.aperture = {rtt::model::SystemApertureType::EntrancePupilDiameter, Param(10.0)};
  s.fields = {FieldType::AngleDeg, {{0.0, 0.0, 1.0}, {0.0, 1.0, 1.0}}};
  s.root.name = "root";
  Surface stop = plane("STO");
  stop.aperture = rtt::model::CircularAperture{10.0, 0.0};
  Surface lens = plane("IL");
  lens.aperture = rtt::model::CircularAperture{30.0, 0.0};
  lens.interaction = rtt::model::IdealLens{Param(100.0), std::nullopt};
  s.root.children = {
      {Element{"stop", ElementKind::Stop, Pose::along_z(-5.0), std::nullopt, {stop}}},
      {Element{"lens", ElementKind::ThinElement, Pose::along_z(0.0), std::nullopt, {lens}}},
      {Element{"D", ElementKind::Detector, Pose::along_z(100.0), std::nullopt, {plane("IMG")}}}};
  s.paths = {{"main", true, {}}};

  rtt::model::SpotGenerator focus;
  focus.path = "main";
  focus.fields = std::vector<std::uint16_t>{0};
  focus.reference = rtt::model::SpotReference::Chief;
  s.tolerances.compensator_merit = rtt::model::Optimization{{}, {focus}};
  s.tolerances.compensators = {{at(ToleranceAnchor::Node, "D", "/pose/position/2"), -2.0, 2.0}};

  rtt::model::SpotRmsOperand spot;
  spot.path = "main";
  spot.reference = rtt::model::SpotReference::Chief;
  Criterion c_spot;
  c_spot.id = "spot";
  c_spot.operand = spot;
  c_spot.max = 1e-3;
  rtt::model::RayOperand chief;
  chief.path = "main";
  chief.surface = SurfaceId("IMG");
  chief.field = 1;
  Criterion c_chief;
  c_chief.id = "chief";
  c_chief.operand = chief;
  c_chief.max_change = 0.05;
  s.tolerances.criteria = {c_spot, c_chief};
  return s;
}

/// A singlet of N-BK7 (R = 50, -50, 5 mm) with a stop and a detector, object at infinity, EPD 10;
/// for the offsets of apply() (no compensator, one criterion).
System singlet() {
  System s;
  s.name = "singlet";
  s.wavelengths = {{0.5876, 1.0, true}};
  s.aperture = {rtt::model::SystemApertureType::EntrancePupilDiameter, Param(10.0)};
  s.fields = {FieldType::AngleDeg, {{0.0, 0.0, 1.0}}};
  s.root.name = "root";
  Surface stop = plane("STO");
  stop.aperture = rtt::model::CircularAperture{8.0, 0.0};
  s.root.children = {
      {Element{"stop", ElementKind::Stop, Pose::along_z(-2.0), std::nullopt, {stop}}},
      {Element{"L",
               ElementKind::Lens,
               Pose::along_z(0.0),
               "SCHOTT:N-BK7",
               {plane("L.S1", 0.0, 50.0), plane("L.S2", 5.0, -50.0)}}},
      {Element{"D", ElementKind::Detector, Pose::along_z(50.0), std::nullopt, {plane("IMG")}}}};
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

TrialContext context(const System& s, const std::vector<ToleranceItem>& items) {
  return rtt::tolerance::prepare(s, library(), nullptr, items);
}

}  // namespace

TEST_CASE("trial: the detector compensator takes back a detector shift", "[tolerance][trial]") {
  // Acceptance #257: the detector moved by shift_z = Dz (rigid overlay) is brought back by its
  // own position Param: compensator end - start = -Dz, and spot_rms is nominal again (0 for the
  // ideal lens: every ray of field 0 meets the focus, ADR 0031).
  // Bound (ADR 0036, point 3): the residuals of the generator are linear in p (transverse errors
  // of rays converging on the focus), so the run stops on the step test ||D h|| <= xtol (||D
  // theta|| + xtol) (MNT eq. (3.15b)) with ||D e|| <= ||D h|| (1 + mu / lambda_min). The two
  // bounds a, b = p0 -+ 2 give p = a + (b - a)/2 (sin theta + 1) (MINUIT, eq. (1.2); bounds.hpp):
  // theta_0 = 0 and dp/dtheta = (b - a)/2 cos theta <= 2; the scaling D is the column norm at the
  // start, so lambda_min = cos^2 theta_end. With D >= xtol (generous): |e_p| <= 2 xtol
  // (|theta_end| + 1) (1 + mu / cos^2 theta_end).
  const double dz = 0.2;
  const System s = ideal_lens();
  const std::vector<ToleranceItem> items = {
      item("D.shift", at(ToleranceAnchor::Node, "D"), ToleranceKind::ShiftZ, 0.5)};
  const TrialContext ctx = context(s, items);
  const std::vector<double> values = {dz};
  const TrialResult r = rtt::tolerance::evaluate(ctx, values);
  REQUIRE(r.status == TrialStatus::Valid);
  REQUIRE(r.compensators.size() == 1);
  REQUIRE(r.lm_status.has_value());
  REQUIRE(r.mu.has_value());
  const double xtol = rtt::tolerance::TrialOptions{}.compensators.xtol;
  const rtt::optim::Bounds bounds{98.0, 102.0};
  const double theta = rtt::optim::to_internal(r.compensators[0].end, bounds);
  const double c2 = std::cos(theta) * std::cos(theta);
  const double bound = 2.0 * xtol * (std::abs(theta) + 1.0) * (1.0 + *r.mu / c2);
  INFO("start " << r.compensators[0].start << ", end " << r.compensators[0].end << ", mu " << *r.mu
                << ", bound " << bound);
  CHECK(r.compensators[0].start == 100.0);
  CHECK(std::abs((r.compensators[0].end - r.compensators[0].start) + dz) <= bound);
  CHECK_FALSE(r.compensators[0].infeasible);
  // spot_rms: with the detector |e_p| from the focus the spot radius is |e_p| times the marginal
  // slope 5/100 at most, so rms <= bound.
  REQUIRE(r.criteria[0].value.has_value());
  CHECK(*r.criteria[0].value <= bound);
  CHECK(r.criteria[0].satisfied == std::optional<bool>(true));
  // Content guard: without the compensator the trial is defocused by dz (rms > 0.1 dz 5/100).
  System fixed = s;
  fixed.tolerances.compensators.clear();
  fixed.tolerances.compensator_merit.reset();
  const TrialResult open = rtt::tolerance::evaluate(context(fixed, items), values);
  REQUIRE(open.criteria[0].value.has_value());
  CHECK(*open.criteria[0].value > 0.1 * dz * 0.05);
  CHECK(open.compensators.empty());
  CHECK(!open.lm_status.has_value());
}

TEST_CASE("trial: a lens shift is followed by the detector", "[tolerance][trial]") {
  // The lens moved by shift_z = Dz moves the focus by Dz: the detector follows by +Dz (same
  // bound as above).
  const double dz = -0.3;
  const System s = ideal_lens();
  const std::vector<ToleranceItem> items = {
      item("lens.shift", at(ToleranceAnchor::Node, "lens"), ToleranceKind::ShiftZ, 0.5)};
  const TrialResult r = rtt::tolerance::evaluate(context(s, items), std::vector<double>{dz});
  REQUIRE(r.status == TrialStatus::Valid);
  REQUIRE(r.mu.has_value());
  const rtt::optim::Bounds bounds{98.0, 102.0};
  const double theta = rtt::optim::to_internal(r.compensators.at(0).end, bounds);
  const double c2 = std::cos(theta) * std::cos(theta);
  const double xtol = rtt::tolerance::TrialOptions{}.compensators.xtol;
  const double bound = 2.0 * xtol * (std::abs(theta) + 1.0) * (1.0 + *r.mu / c2);
  CHECK(std::abs((r.compensators[0].end - r.compensators[0].start) - dz) <= bound);
}

TEST_CASE("trial: Phi_0 is the trial at zero, compensators included", "[tolerance][trial]") {
  // ADR 0036, points 2 and 6: nominal() is computed with the same procedure, and evaluate() at
  // all zeros gives it bitwise again.
  const System s = ideal_lens();
  const std::vector<ToleranceItem> items = {
      item("lens.dx", at(ToleranceAnchor::Node, "lens"), ToleranceKind::DecenterX, 0.1),
      item("D.shift", at(ToleranceAnchor::Node, "D"), ToleranceKind::ShiftZ, 0.5)};
  const TrialContext ctx = context(s, items);
  const TrialResult& nominal = ctx.nominal();
  REQUIRE(nominal.status == TrialStatus::Valid);
  REQUIRE(nominal.criteria.size() == 2);
  REQUIRE(nominal.criteria[1].value.has_value());
  // The chief ray of 1 degree meets the detector at 100 tan(1 deg) (collinear mapping).
  CHECK(std::abs(*nominal.criteria[1].value - 100.0 * std::tan(std::numbers::pi / 180.0)) < 1e-9);
  CHECK(nominal.criteria[1].satisfied == std::optional<bool>(true));  // |Phi - Phi_0| = 0
  require_same(rtt::tolerance::evaluate(ctx, std::vector<double>{0.0, 0.0}), nominal);
  // The nominal system is in focus: the compensator ends where it started.
  CHECK(std::abs(nominal.compensators.at(0).end - 100.0) < 1e-6);
}

TEST_CASE("trial: the values are applied as offsets and as the overlay", "[tolerance][trial]") {
  // ADR 0036, point 1: (a) offsets on Params and rows, (b) one rigid motion per target with the
  // pivot of the target, (c) material changes per element and segment.
  System s = singlet();
  rtt::model::TolerancePivotOverride pivot;
  pivot.target = at(ToleranceAnchor::Node, "L");
  pivot.pivot = {0.0, 0.0, 2.5};
  s.tolerances.overrides = {pivot};
  const std::vector<ToleranceItem> items = {
      item("r1", at(ToleranceAnchor::Surface, "L.S1"), ToleranceKind::CurvatureRelative, 0.01),
      item("r2", at(ToleranceAnchor::Surface, "L.S2"), ToleranceKind::Radius, 1.0),
      item("k1", at(ToleranceAnchor::Surface, "L.S1"), ToleranceKind::Conic, 0.1),
      item("t", at(ToleranceAnchor::Surface, "L.S2"), ToleranceKind::Thickness, 0.1),
      item("dx", at(ToleranceAnchor::Node, "L"), ToleranceKind::DecenterX, 0.1),
      item("tx", at(ToleranceAnchor::Node, "L"), ToleranceKind::TiltXDeg, 0.1),
      item("sx", at(ToleranceAnchor::Surface, "L.S2"), ToleranceKind::DecenterY, 0.1),
      item("n", at(ToleranceAnchor::Node, "L"), ToleranceKind::Index, 1e-3),
      item("v", at(ToleranceAnchor::Node, "L"), ToleranceKind::Abbe, 0.5),
      item("x", at(ToleranceAnchor::Node, "D", "/pose/position/0"), ToleranceKind::Param, 0.1),
  };
  const TrialContext ctx = context(s, items);
  const std::vector<double> values = {0.004, -0.5, 0.05, 0.02, 0.03, -0.07, 0.01, 2e-4, -0.3, 0.06};
  const rtt::tolerance::Applied a = rtt::tolerance::apply(ctx, values);
  const auto& lens = std::get<Element>(a.system.root.children[1].value);
  const auto& s1 = std::get<rtt::model::Conic>(lens.surfaces[0].shape.base);
  const auto& s2 = std::get<rtt::model::Conic>(lens.surfaces[1].shape.base);
  CHECK(same(s1.radius.value, 50.0 / (1.0 + 0.004)));  // R' = R / (1 + delta)
  CHECK(same(s2.radius.value, -50.0 - 0.5));           // R' = R + dR
  CHECK(same(s1.conic.value, 0.05));
  CHECK(same(lens.surfaces[1].pose.position[2].value, 5.0 + 0.02));
  const auto& detector = std::get<Element>(a.system.root.children[2].value);
  CHECK(same(detector.pose.position[0].value, 0.06));
  // The rigid motions: one per target, the pivot of the target (override).
  REQUIRE(a.perturbation.nodes.size() == 1);
  const rtt::compile::RigidMotion& m = a.perturbation.nodes.at("L");
  CHECK(same(m.translation_mm.x(), 0.03));
  CHECK(m.translation_mm.y() == 0.0);
  CHECK(same(m.rotation_rad.x(), -0.07 * std::numbers::pi / 180.0));
  REQUIRE(m.pivot_mm.has_value());
  CHECK(*m.pivot_mm == rtt::math::Vec3(0.0, 0.0, 2.5));
  REQUIRE(a.perturbation.surfaces.size() == 1);
  CHECK(same(a.perturbation.surfaces.at(SurfaceId("L.S2")).translation_mm.y(), 0.01));
  CHECK(!a.perturbation.surfaces.at(SurfaceId("L.S2")).pivot_mm.has_value());
  REQUIRE(a.perturbation.materials.size() == 1);
  const rtt::compile::MaterialChange& mc =
      a.perturbation.materials.at(rtt::compile::MaterialTarget{"L", std::nullopt});
  CHECK(same(mc.delta_nd, 2e-4));
  CHECK(same(mc.delta_vd, -0.3));
  // The nominal model is unchanged.
  const auto& nominal = std::get<Element>(s.root.children[1].value);
  CHECK(std::get<rtt::model::Conic>(nominal.surfaces[0].shape.base).radius.value == 50.0);
}

TEST_CASE("trial: a parameter item shifts its row in every configuration", "[tolerance][trial]") {
  // ADR 0036, point 1a (K6).
  System s = singlet();
  s.configurations = {{"near"}, {"far"}};
  rtt::model::ParameterRow gap;
  gap.name = "GAP";
  gap.form = std::vector<double>{50.0, 52.0};
  s.parameters = {gap};
  std::get<Element>(s.root.children[2].value).pose.position[2] = Param::bound("GAP");
  const std::vector<ToleranceItem> items = {
      item("gap", at(ToleranceAnchor::Parameter, "GAP"), ToleranceKind::Parameter, 0.5)};
  const rtt::tolerance::Applied a =
      rtt::tolerance::apply(context(s, items), std::vector<double>{0.25});
  REQUIRE(a.system.parameters.size() == 1);
  const auto& values = std::get<std::vector<double>>(a.system.parameters[0].form);
  REQUIRE(values.size() == 2);
  CHECK(same(values[0], 50.25));
  CHECK(same(values[1], 52.25));
}

TEST_CASE("trial: a criterion without value invalidates only itself", "[tolerance][trial]") {
  // The detector gets a clear aperture of 1 mm: the chief ray of field 1 (1.745 mm) misses it.
  // Its criterion has no value with the reason, status Invalid (Monte Carlo); spot_rms of field 0
  // keeps its value.
  System s = ideal_lens();
  auto& d = std::get<Element>(s.root.children[2].value);
  d.surfaces[0].aperture = rtt::model::CircularAperture{1.0, 0.0};
  const TrialContext ctx = context(s, {});
  const TrialResult& r = ctx.nominal();
  CHECK(r.status == TrialStatus::Invalid);
  CHECK(!r.reason.empty());
  REQUIRE(r.criteria.size() == 2);
  CHECK(r.criteria[0].value.has_value());
  CHECK(!r.criteria[1].value.has_value());
  CHECK(!r.criteria[1].reason.empty());
  CHECK(!r.criteria[1].satisfied.has_value());
}

TEST_CASE("trial: a compensator run that fails invalidates the trial", "[tolerance][trial]") {
  // The lens decentered by 40 mm: every ray misses its aperture of 30 mm, the generator of the
  // compensator merit has no value at the start. No criterion is evaluated.
  const System s = ideal_lens();
  const std::vector<ToleranceItem> items = {
      item("lens.dy", at(ToleranceAnchor::Node, "lens"), ToleranceKind::DecenterY, 50.0)};
  const TrialResult r = rtt::tolerance::evaluate(context(s, items), std::vector<double>{40.0});
  CHECK(r.status == TrialStatus::Invalid);
  CHECK(!r.reason.empty());
  for (const auto& c : r.criteria) CHECK(!c.value.has_value());
}

TEST_CASE("trial: a compensator without inner start stays and warns", "[tolerance][trial]") {
  // ADR 0036, point 3: the detector position has max = 100, its perturbed value: no inner start,
  // the compensator stays at 100 with tolerances.compensator_infeasible; the trial stays valid.
  System s = ideal_lens();
  std::get<Element>(s.root.children[2].value).pose.position[2].max = 100.0;
  const std::vector<ToleranceItem> items = {
      item("D.shift", at(ToleranceAnchor::Node, "D"), ToleranceKind::ShiftZ, 0.5)};
  const TrialResult r = rtt::tolerance::evaluate(context(s, items), std::vector<double>{0.2});
  CHECK(r.status == TrialStatus::Valid);
  REQUIRE(r.compensators.size() == 1);
  CHECK(r.compensators[0].infeasible);
  CHECK(r.compensators[0].end == 100.0);
  CHECK(!r.lm_status.has_value());
  REQUIRE(!r.warnings.empty());
  CHECK(r.warnings[0].code == "tolerances.compensator_infeasible");
  CHECK(r.warnings[0].location == "/tolerances/compensators/0");
}

TEST_CASE("trial: bitwise the same for 1, 4 and all threads and in parallel",
          "[tolerance][trial]") {
  const System s = ideal_lens();
  const std::vector<ToleranceItem> items = {
      item("lens.dx", at(ToleranceAnchor::Node, "lens"), ToleranceKind::DecenterX, 0.1),
      item("lens.tx", at(ToleranceAnchor::Node, "lens"), ToleranceKind::TiltXDeg, 0.2),
      item("D.shift", at(ToleranceAnchor::Node, "D"), ToleranceKind::ShiftZ, 0.5)};
  const TrialContext ctx = context(s, items);
  const std::vector<std::vector<double>> trials = {
      {0.05, 0.1, -0.2}, {-0.08, -0.15, 0.3}, {0.0, 0.2, 0.1}, {0.1, 0.0, -0.4}};
  std::vector<TrialResult> serial;
  {
    oneapi::tbb::task_arena one(1);
    one.execute([&] {
      for (const auto& v : trials) serial.push_back(rtt::tolerance::evaluate(ctx, v));
    });
  }
  for (const int threads : {4, oneapi::tbb::info::default_concurrency()}) {
    INFO(threads << " threads");
    std::vector<TrialResult> parallel(trials.size());
    oneapi::tbb::task_arena arena(threads);
    arena.execute([&] {
      oneapi::tbb::parallel_for(std::size_t{0}, trials.size(), [&](std::size_t i) {
        parallel[i] = rtt::tolerance::evaluate(ctx, trials[i]);
      });
    });
    for (std::size_t i = 0; i < trials.size(); ++i) require_same(parallel[i], serial[i]);
  }
}

TEST_CASE("trial: input errors are exceptions, cancellation a status", "[tolerance][trial]") {
  const System s = ideal_lens();
  const std::vector<ToleranceItem> items = {
      item("lens.dx", at(ToleranceAnchor::Node, "lens"), ToleranceKind::DecenterX, 0.1)};
  const TrialContext ctx = context(s, items);
  CHECK_THROWS_AS(rtt::tolerance::evaluate(ctx, std::vector<double>{}), std::invalid_argument);
  CHECK_THROWS_AS(rtt::tolerance::evaluate(ctx, std::vector<double>{0.1, 0.2}),
                  std::invalid_argument);
  CHECK_THROWS_AS(
      rtt::tolerance::evaluate(ctx, std::vector<double>{std::numeric_limits<double>::quiet_NaN()}),
      std::invalid_argument);
  // An item whose target does not exist.
  CHECK_THROWS_AS(
      context(s, {item("x", at(ToleranceAnchor::Node, "nowhere"), ToleranceKind::DecenterX, 0.1)}),
      std::invalid_argument);
  // A material item at a metal: no index at the lines (compile::fraunhofer_indices, #260).
  System metal = s;
  auto& lens = std::get<Element>(metal.root.children[1].value);
  lens.kind = ElementKind::Plate;
  lens.material = "CONST:0.2,3.0";
  lens.surfaces = {plane("IL"), plane("IL2", 1.0)};
  CHECK_THROWS_AS(
      context(metal, {item("n", at(ToleranceAnchor::Node, "lens"), ToleranceKind::Index, 1e-3)}),
      rtt::compile::CompileError);
  // A criterion the nominal system does not support: magnification with the object at infinity.
  System unsupported = s;
  rtt::model::FirstOrderOperand m;
  m.quantity = rtt::model::FirstOrderQuantity::Magnification;
  m.path = "main";
  Criterion c;
  c.id = "m";
  c.operand = m;
  c.max = 1.0;
  unsupported.tolerances.criteria.push_back(c);
  CHECK_THROWS_AS(context(unsupported, items), rtt::optim::OptimError);

  // Cancellation requested before: the trial is Cancelled; prepare's Phi_0 too.
  rtt::trace::RunControl control;
  control.cancel.emplace();
  control.cancel->request_cancel();
  const TrialResult r = rtt::tolerance::evaluate(ctx, std::vector<double>{0.05}, &control);
  CHECK(r.status == TrialStatus::Cancelled);
  for (const auto& v : r.criteria) CHECK(!v.value.has_value());
  const TrialContext cancelled =
      rtt::tolerance::prepare(s, library(), nullptr, items, {}, &control);
  CHECK(cancelled.nominal().status == TrialStatus::Cancelled);
}

TEST_CASE("trial: satisfies takes every given bound, bounds included", "[tolerance][trial]") {
  using rtt::tolerance::satisfies;
  Criterion c;
  c.min = 1.0;
  c.max = 2.0;
  CHECK(satisfies(c, 1.0, std::nullopt) == std::optional<bool>(true));
  CHECK(satisfies(c, 2.0, std::nullopt) == std::optional<bool>(true));
  CHECK(satisfies(c, 2.5, std::nullopt) == std::optional<bool>(false));
  CHECK(satisfies(c, 0.5, 1.0) == std::optional<bool>(false));
  CHECK(!satisfies(c, std::nullopt, 1.0).has_value());
  c.max_change = 0.125;  // binary fractions: the differences below are exact
  CHECK(satisfies(c, 1.625, 1.5) == std::optional<bool>(true));
  CHECK(satisfies(c, 1.625, 1.375) == std::optional<bool>(false));
  CHECK(!satisfies(c, 1.625, std::nullopt).has_value());
}
