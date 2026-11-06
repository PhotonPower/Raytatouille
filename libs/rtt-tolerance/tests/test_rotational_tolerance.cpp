// Tolerances on the odd asphere and the axicon (#314; ADR 0035 and ADR 0036, addenda #314):
// curvature_relative scales the curvature, C' = C (1 + delta); a curvature that reaches 0 or
// changes its sign leaves the trial valid on every path (item, param, row), because c = 0 is a
// plane base and c changes continuously through 0. The axicon angle makes a trial invalid when
// it changes its sign or reaches 90 deg in magnitude.
//
// The system is the singlet of test_radius_sign.cpp (N-BK7, R2 = -50 mm) with L.S1 replaced.

#include <catch2/catch_test_macros.hpp>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
#include "rtt/model/tolerance.hpp"
#include "rtt/tolerance/trial.hpp"

using rtt::model::Criterion;
using rtt::model::Element;
using rtt::model::ElementKind;
using rtt::model::FieldType;
using rtt::model::Param;
using rtt::model::ParameterRow;
using rtt::model::Pose;
using rtt::model::Surface;
using rtt::model::SurfaceId;
using rtt::model::System;
using rtt::model::ToleranceAnchor;
using rtt::model::ToleranceDistribution;
using rtt::model::ToleranceItem;
using rtt::model::ToleranceKind;
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

Surface surface(const std::string& id, double z_mm, rtt::model::BaseShape base) {
  Surface s;
  s.id = SurfaceId(id);
  s.pose = Pose::along_z(z_mm);
  s.shape.base = std::move(base);
  return s;
}

/// The singlet of test_radius_sign.cpp with `front` as L.S1: stop at z = -2, lens "L" (L.S1 at
/// 0, L.S2 at 5 mm with R = -50, N-BK7), detector at 50; object at infinity, EPD 10; with
/// `criterion` the RMS spot of field 0 as the one criterion.
System singlet(rtt::model::BaseShape front, bool criterion) {
  System s;
  s.name = "singlet";
  s.wavelengths = {{0.5876, 1.0, true}};
  s.aperture = {rtt::model::SystemApertureType::EntrancePupilDiameter, Param(10.0)};
  s.fields = {FieldType::AngleDeg, {{0.0, 0.0, 1.0}}};
  s.root.name = "root";
  Surface stop = surface("STO", 0.0, rtt::model::Plane{});
  stop.aperture = rtt::model::CircularAperture{8.0, 0.0};
  s.root.children = {
      {Element{"stop", ElementKind::Stop, Pose::along_z(-2.0), std::nullopt, {stop}}},
      {Element{"L",
               ElementKind::Lens,
               Pose::along_z(0.0),
               "SCHOTT:N-BK7",
               {surface("L.S1", 0.0, std::move(front)),
                surface("L.S2", 5.0, rtt::model::Conic{Param(-50.0), Param(0.0)})}}},
      {Element{"D",
               ElementKind::Detector,
               Pose::along_z(50.0),
               std::nullopt,
               {surface("IMG", 0.0, rtt::model::Plane{})}}}};
  s.paths = {{"main", true, {}}};
  if (criterion) {
    rtt::model::SpotRmsOperand spot;
    spot.path = "main";
    Criterion c;
    c.id = "spot";
    c.operand = spot;
    c.max = 1.0;
    s.tolerances.criteria = {c};
  }
  return s;
}

/// An odd asphere with base curvature `c` (a Param) and an A4 term: order 4, a paraxial model.
rtt::model::OddAsphere odd_asphere(Param c) {
  return {std::move(c), Param(0.0), {Param(0.0), Param(0.0), Param(0.0), Param(1e-7)}};
}

ToleranceItem surface_item(ToleranceKind kind, double limit) {
  ToleranceItem i;
  i.id = "i";
  i.target.anchor = ToleranceAnchor::Surface;
  i.target.name = "L.S1";
  i.kind = kind;
  i.limit = limit;
  i.distribution = ToleranceDistribution::Uniform;
  return i;
}

ToleranceItem param_item(const std::string& pointer, double limit) {
  ToleranceItem i = surface_item(ToleranceKind::Param, limit);
  i.target.param = pointer;
  return i;
}

ToleranceItem row_item(const std::string& row, double limit) {
  ToleranceItem i;
  i.id = "row";
  i.target.anchor = ToleranceAnchor::Parameter;
  i.target.name = row;
  i.kind = ToleranceKind::Parameter;
  i.limit = limit;
  i.distribution = ToleranceDistribution::Uniform;
  return i;
}

TrialContext context(const System& s, const ToleranceItem& item) {
  return rtt::tolerance::prepare(s, library(), nullptr, std::vector<ToleranceItem>{item});
}

/// The base curvature of L.S1 in an applied trial.
double curvature(const System& s) {
  const auto& lens = std::get<Element>(s.root.children[1].value);
  return std::get<rtt::model::OddAsphere>(lens.surfaces[0].shape.base).curvature.value;
}

bool contains(const std::string& text, const std::string& part) {
  return text.find(part) != std::string::npos;
}

/// A valid trial with a value of the spot criterion.
void check_valid(const TrialResult& r) {
  INFO("reason: " << r.reason);
  CHECK(r.status == TrialStatus::Valid);
  REQUIRE(r.criteria.size() == 1);
  CHECK(r.criteria[0].value.has_value());
}

}  // namespace

TEST_CASE("curvature_relative scales the curvature of an odd asphere",
          "[tolerance][trial][rotational]") {
  // ADR 0035, addendum #314: C' = C (1 + delta) on /shape/base/curvature.
  const System s = singlet(odd_asphere(Param(0.02)), true);
  const TrialContext ctx = context(s, surface_item(ToleranceKind::CurvatureRelative, 0.5));
  REQUIRE(ctx.nominal().status == TrialStatus::Valid);
  REQUIRE(ctx.nominal().criteria[0].value.has_value());
  for (const double delta : {0.25, -0.1, -1.0, -1.5}) {
    INFO("delta " << delta);
    const double expected = 0.02 * (1.0 + delta);
    CHECK(curvature(rtt::tolerance::apply(ctx, std::vector<double>{delta}).system) == expected);
  }
  // A curvature that reaches 0 (delta = -1) or changes its sign (delta = -1.5) is a valid trial:
  // the base becomes a plane or bends the other way (ADR 0036, addendum #314).
  for (const double delta : {0.25, -1.0, -1.5}) {
    INFO("delta " << delta);
    check_valid(rtt::tolerance::evaluate(ctx, std::vector<double>{delta}));
  }
  // Content guard: the criterion follows the curvature.
  const TrialResult flat = rtt::tolerance::evaluate(ctx, std::vector<double>{-1.0});
  REQUIRE(flat.criteria[0].value.has_value());
  CHECK(*flat.criteria[0].value != *ctx.nominal().criteria[0].value);
}

TEST_CASE("curvature_relative leaves a plane base of an odd asphere at 0",
          "[tolerance][trial][rotational]") {
  // ADR 0035, addendum #314: a nominal c = 0 is not changed by curvature_relative.
  const System s = singlet(odd_asphere(Param(0.0)), true);
  const TrialContext ctx = context(s, surface_item(ToleranceKind::CurvatureRelative, 0.5));
  for (const double delta : {0.3, -2.0}) {
    CHECK(curvature(rtt::tolerance::apply(ctx, std::vector<double>{delta}).system) == 0.0);
  }
}

TEST_CASE("a param item that turns the curvature through 0 leaves the trial valid",
          "[tolerance][trial][rotational]") {
  const System s = singlet(odd_asphere(Param(0.02)), true);
  const TrialContext ctx = context(s, param_item("/shape/base/curvature", 0.05));
  REQUIRE(ctx.nominal().status == TrialStatus::Valid);
  for (const double offset : {-0.02, -0.03, 0.01}) {
    INFO("offset " << offset);
    CHECK(curvature(rtt::tolerance::apply(ctx, std::vector<double>{offset}).system) ==
          0.02 + offset);
    check_valid(rtt::tolerance::evaluate(ctx, std::vector<double>{offset}));
  }
}

TEST_CASE("a row that turns a bound curvature through 0 leaves the trial valid",
          "[tolerance][trial][rotational]") {
  // ADR 0036, addendum #314: curvature is not one of the radius fields of the sign check #293.
  System s = singlet(odd_asphere(Param::bound("C")), true);
  ParameterRow row;
  row.name = "C";
  row.form = 0.02;
  s.parameters = {row};
  const TrialContext ctx = context(s, row_item("C", 0.05));
  REQUIRE(ctx.nominal().status == TrialStatus::Valid);
  for (const double offset : {-0.02, -0.03}) {
    INFO("offset " << offset);
    check_valid(rtt::tolerance::evaluate(ctx, std::vector<double>{offset}));
  }
}

TEST_CASE("an axicon angle that changes its sign or reaches 90 deg makes the trial invalid",
          "[tolerance][trial][rotational]") {
  // ADR 0035 and ADR 0036, addenda #314. The axicon has no paraxial model (ADR 0043, F6), so
  // the trials run without a criterion; only the status counts.
  const System s = singlet(rtt::model::Axicon{Param(-5.0)}, false);
  const TrialContext ctx = context(s, param_item("/shape/base/angle_deg", 100.0));
  REQUIRE(ctx.nominal().status == TrialStatus::Valid);
  for (const double offset : {10.0, 5.5, -85.0, -95.0, -180.0}) {
    INFO("offset " << offset);
    const TrialResult r = rtt::tolerance::evaluate(ctx, std::vector<double>{offset});
    INFO("reason: " << r.reason);
    CHECK(r.status == TrialStatus::Invalid);
    CHECK(contains(r.reason, "surface 'L.S1'"));
    CHECK(contains(r.reason, "angle_deg"));
  }
  // Content guard: the same sign and below 90 deg in magnitude stays valid, also at the plane
  // alpha = 0, which is no change of sign.
  for (const double offset : {4.0, 5.0, -84.0}) {
    INFO("offset " << offset);
    const TrialResult r = rtt::tolerance::evaluate(ctx, std::vector<double>{offset});
    INFO("reason: " << r.reason);
    CHECK(r.status == TrialStatus::Valid);
  }
}

TEST_CASE("an axicon with a flat nominal: any sign below 90 deg is valid",
          "[tolerance][trial][rotational]") {
  // ADR 0043, addendum PR B: a nominal alpha = 0 (the plane) has no sign to change; only
  // |alpha| >= 90 deg makes the trial invalid.
  const System s = singlet(rtt::model::Axicon{Param(0.0)}, false);
  const TrialContext ctx = context(s, param_item("/shape/base/angle_deg", 100.0));
  REQUIRE(ctx.nominal().status == TrialStatus::Valid);
  for (const double offset : {3.0, -3.0}) {
    const TrialResult r = rtt::tolerance::evaluate(ctx, std::vector<double>{offset});
    INFO("offset " << offset << ", reason: " << r.reason);
    CHECK(r.status == TrialStatus::Valid);
  }
  for (const double offset : {95.0, -95.0}) {
    const TrialResult r = rtt::tolerance::evaluate(ctx, std::vector<double>{offset});
    INFO("offset " << offset << ", reason: " << r.reason);
    CHECK(r.status == TrialStatus::Invalid);
    CHECK(contains(r.reason, "angle_deg"));
  }
}

TEST_CASE("an axicon angle bound to a row is checked per row value",
          "[tolerance][trial][rotational]") {
  System s = singlet(rtt::model::Axicon{Param::bound("A")}, false);
  ParameterRow row;
  row.name = "A";
  row.form = 7.0;
  s.parameters = {row};
  const TrialContext ctx = context(s, row_item("A", 100.0));
  REQUIRE(ctx.nominal().status == TrialStatus::Valid);
  const TrialResult flipped = rtt::tolerance::evaluate(ctx, std::vector<double>{-10.0});
  INFO("reason: " << flipped.reason);
  CHECK(flipped.status == TrialStatus::Invalid);
  CHECK(contains(flipped.reason, "angle_deg"));
  CHECK(contains(flipped.reason, "row 'A'"));
  CHECK(rtt::tolerance::evaluate(ctx, std::vector<double>{3.0}).status == TrialStatus::Valid);
  CHECK(rtt::tolerance::evaluate(ctx, std::vector<double>{83.0}).status == TrialStatus::Invalid);
}
