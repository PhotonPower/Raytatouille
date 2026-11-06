// curvature_relative on Q-con and Q-bfs (#314; ADR 0035 and ADR 0036, addenda #314): C' = C
// (1 + delta) on /shape/base/curvature, and a curvature that reaches 0 or changes its sign leaves
// the trial valid, also for a limit >= 1 (ADR 0043, addendum PR C).
//
// The system is the singlet of test_radius_sign.cpp (N-BK7, R2 = -50 mm) with L.S1 replaced.

#include <catch2/catch_test_macros.hpp>
#include <string>
#include <type_traits>
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

TrialContext context(const System& s, const ToleranceItem& item) {
  return rtt::tolerance::prepare(s, library(), nullptr, std::vector<ToleranceItem>{item});
}

/// The base curvature of L.S1 (Q-con or Q-bfs) in an applied trial.
double curvature(const System& s) {
  const auto& lens = std::get<Element>(s.root.children[1].value);
  return std::visit(
      [](const auto& b) -> double {
        using B = std::decay_t<decltype(b)>;
        if constexpr (std::is_same_v<B, rtt::model::QCon> || std::is_same_v<B, rtt::model::QBfs>) {
          return b.curvature.value;
        } else {
          FAIL("not a Forbes base");
          return 0.0;
        }
      },
      lens.surfaces[0].shape.base);
}

/// A valid trial with a value of the spot criterion.
void check_valid(const TrialResult& r) {
  INFO("reason: " << r.reason);
  CHECK(r.status == TrialStatus::Valid);
  REQUIRE(r.criteria.size() == 1);
  CHECK(r.criteria[0].value.has_value());
}

}  // namespace

TEST_CASE("curvature_relative scales the curvature of Q-con and Q-bfs, through 0 as well",
          "[tolerance][trial][forbes]") {
  for (const rtt::model::BaseShape& base :
       {rtt::model::BaseShape{
            rtt::model::QCon{Param(0.02), Param(0.0), Param(5.0), {Param(1e-3), Param(-2e-4)}}},
        rtt::model::BaseShape{
            rtt::model::QBfs{Param(0.02), Param(0.0), Param(5.0), {Param(1e-3), Param(-2e-4)}}}}) {
    const System s = singlet(base, true);
    // A limit >= 1 is allowed on a curvature (ADR 0043, addendum PR C).
    const TrialContext ctx = context(s, surface_item(ToleranceKind::CurvatureRelative, 1.5));
    REQUIRE(ctx.nominal().status == TrialStatus::Valid);
    REQUIRE(ctx.nominal().criteria[0].value.has_value());
    for (const double delta : {0.25, -1.0, -1.5}) {
      INFO("delta " << delta);
      CHECK(curvature(rtt::tolerance::apply(ctx, std::vector<double>{delta}).system) ==
            0.02 * (1.0 + delta));
      check_valid(rtt::tolerance::evaluate(ctx, std::vector<double>{delta}));
    }
  }
}
