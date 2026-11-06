// The radius sign check of ADR 0036, addendum #293, on the radii of the anamorphic shapes (ADR
// 0042): radius_x and radius_y of the biconic and radius_of_rotation of the toroid. A trial that
// turns one of them to zero, to infinity or to the other sign is Invalid, as for the base radius
// of conic and even asphere (test_radius_sign.cpp).
//
// The system is the singlet of test_radius_sign.cpp (N-BK7, R2 = -50 mm) with L.S1 anamorphic in
// form but rotationally symmetric in the nominal (all radii 50 mm), so the nominal criterion has a
// value. The invalid cases put the radius clearly on the other side of 0, the content guards
// clearly on the same side.

#include <catch2/catch_test_macros.hpp>
#include <string>
#include <utility>
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

/// The singlet of test_radius_sign.cpp with `front` as L.S1 and the row R = 50: stop at z = -2,
/// lens "L" (L.S1 at 0, L.S2 at 5 mm with R = -50, N-BK7), detector at 50; object at infinity,
/// EPD 10; one criterion, the RMS spot of field 0.
System singlet(rtt::model::BaseShape front) {
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
  ParameterRow r;
  r.name = "R";
  r.form = 50.0;
  s.parameters = {r};
  rtt::model::SpotRmsOperand spot;
  spot.path = "main";
  Criterion c;
  c.id = "spot";
  c.operand = spot;
  c.max = 1.0;
  s.tolerances.criteria = {c};
  return s;
}

/// A parameter item on the row R, uniform within +-100 mm.
TrialContext context(const System& s) {
  ToleranceItem i;
  i.id = "r";
  i.target.anchor = ToleranceAnchor::Parameter;
  i.target.name = "R";
  i.kind = ToleranceKind::Parameter;
  i.limit = 100.0;
  i.distribution = ToleranceDistribution::Uniform;
  return rtt::tolerance::prepare(s, library(), nullptr, std::vector<ToleranceItem>{i});
}

bool contains(const std::string& text, const std::string& part) {
  return text.find(part) != std::string::npos;
}

/// An Invalid trial on L.S1 and the row R, whose reason names `field`.
void check_invalid(const TrialResult& r, const std::string& field) {
  INFO("reason: " << r.reason);
  CHECK(r.status == TrialStatus::Invalid);
  CHECK(contains(r.reason, "surface 'L.S1'"));
  CHECK(contains(r.reason, "the " + field + " "));
  CHECK(contains(r.reason, "row 'R'"));
  CHECK(contains(r.reason, "other sign"));
}

/// A trial that the radius check lets through. Its surface may be anamorphic, which the criterion
/// (aiming) does not support yet, so only the status is checked.
void check_not_flipped(const TrialResult& r) {
  INFO("reason: " << r.reason);
  CHECK(r.status == TrialStatus::Valid);
  CHECK(!contains(r.reason, "other sign"));
}

}  // namespace

TEST_CASE("radius sign: radius_x and radius_y of a biconic on one row",
          "[tolerance][trial][radius][anamorphic]") {
  // Both radii bound to R = 50, so L.S1 stays a sphere. Offset -70: R' = -20 for both, radius_x
  // is reported first; offset -20: R' = 30, a sphere, the criterion has a value.
  const System s =
      singlet(rtt::model::Biconic{Param::bound("R"), Param::bound("R"), Param(0.0), Param(0.0)});
  const TrialContext ctx = context(s);
  REQUIRE(ctx.nominal().status == TrialStatus::Valid);
  REQUIRE(ctx.nominal().criteria.size() == 1);
  REQUIRE(ctx.nominal().criteria[0].value.has_value());
  check_invalid(rtt::tolerance::evaluate(ctx, std::vector<double>{-70.0}), "radius_x");
  const TrialResult ok = rtt::tolerance::evaluate(ctx, std::vector<double>{-20.0});
  check_not_flipped(ok);
  REQUIRE(ok.criteria.size() == 1);
  CHECK(ok.criteria[0].value.has_value());
}

TEST_CASE("radius sign: radius_y of a biconic", "[tolerance][trial][radius][anamorphic]") {
  // radius_x = 50 fixed, radius_y bound to R = 50. Offset -70: R_y' = -20; offset -20: 30.
  const System s =
      singlet(rtt::model::Biconic{Param(50.0), Param::bound("R"), Param(0.0), Param(0.0)});
  const TrialContext ctx = context(s);
  REQUIRE(ctx.nominal().status == TrialStatus::Valid);
  check_invalid(rtt::tolerance::evaluate(ctx, std::vector<double>{-70.0}), "radius_y");
  check_not_flipped(rtt::tolerance::evaluate(ctx, std::vector<double>{-20.0}));
}

TEST_CASE("radius sign: radius_of_rotation of a toroid", "[tolerance][trial][radius][anamorphic]") {
  // Profile radius 50 fixed, radius_of_rotation bound to R = 50: a sphere in the nominal. Offset
  // -70: R' = -20; offset -20: 30.
  const System s = singlet(rtt::model::Toroid{Param(50.0), Param(0.0), {}, Param::bound("R")});
  const TrialContext ctx = context(s);
  REQUIRE(ctx.nominal().status == TrialStatus::Valid);
  check_invalid(rtt::tolerance::evaluate(ctx, std::vector<double>{-70.0}), "radius_of_rotation");
  check_not_flipped(rtt::tolerance::evaluate(ctx, std::vector<double>{-20.0}));
}
