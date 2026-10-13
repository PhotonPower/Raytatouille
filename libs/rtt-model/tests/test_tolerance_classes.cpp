// The items a tolerance class generates (ADR 0035, point 5; #253 part B): targets, kinds, ids,
// order, skips with their reasons, single-value overrides and the precedence of explicit items
// by perturbed quantity; and the FNV-1a hash of the ids (ADR 0036, point 5).

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "rtt/model/model.hpp"
#include "rtt/model/tolerance_classes.hpp"
#include "test_support.hpp"
#include "tolerance_fixture.hpp"

using namespace rtt::model;
using namespace rtt::model::test;

namespace {

using K = ToleranceKind;
using A = ToleranceAnchor;

/// (id, skip, effective class) of every listed item.
struct Row {
  std::string id;
  ClassSkip skip;
  ToleranceClass cls;
  bool operator==(const Row&) const = default;
};

std::vector<Row> rows(const ClassTolerances& r) {
  std::vector<Row> out;
  out.reserve(r.items.size());
  for (const ClassTolerance& c : r.items) out.push_back({c.id, c.skip, c.tolerance_class});
  return out;
}

/// The item with this id, by value: a reference returned from a call with a string literal is
/// a dangling temporary to GCC 13 (-Wdangling-reference).
ClassTolerance find(const ClassTolerances& r, const std::string& id) {
  for (const ClassTolerance& c : r.items) {
    if (c.id == id) return c;
  }
  FAIL("no generated item " << id);
  return {};
}

}  // namespace

TEST_CASE("FNV-1a 64 of tolerance ids (RFC 9923, Sec. 2 and 5)", "[model][tolerances]") {
  // The empty input gives offset_basis; "a" and "foobar" from the algorithm of Sec. 2 with the
  // constants of Sec. 5, octet by octet.
  CHECK(tolerance_id_hash("") == 0xCBF29CE484222325ULL);
  CHECK(tolerance_id_hash("a") == 0xAF63DC4C8601EC8CULL);
  CHECK(tolerance_id_hash("foobar") == 0x85944171F73967E8ULL);
  // UTF-8 bytes, not characters: "Ä" is the two octets C3 84.
  CHECK(tolerance_id_hash("\xC3\x84") == 0x0AC1F207B717DF22ULL);
  CHECK(tolerance_id_hash("class:node:L1:index:0") == 0xCA314017AC320210ULL);
}

TEST_CASE("tolerance kind names are the file names", "[model][tolerances]") {
  CHECK(tolerance_kind_name(K::Radius) == "radius");
  CHECK(tolerance_kind_name(K::CurvatureRelative) == "curvature_relative");
  CHECK(tolerance_kind_name(K::Conic) == "conic");
  CHECK(tolerance_kind_name(K::Thickness) == "thickness");
  CHECK(tolerance_kind_name(K::Param) == "param");
  CHECK(tolerance_kind_name(K::Parameter) == "parameter");
  CHECK(tolerance_kind_name(K::DecenterX) == "decenter_x");
  CHECK(tolerance_kind_name(K::DecenterY) == "decenter_y");
  CHECK(tolerance_kind_name(K::ShiftZ) == "shift_z");
  CHECK(tolerance_kind_name(K::TiltXDeg) == "tilt_x_deg");
  CHECK(tolerance_kind_name(K::TiltYDeg) == "tilt_y_deg");
  CHECK(tolerance_kind_name(K::TiltZDeg) == "tilt_z_deg");
  CHECK(tolerance_kind_name(K::Index) == "index");
  CHECK(tolerance_kind_name(K::Abbe) == "abbe");
}

TEST_CASE("class: targets, kinds, ids, order and skips", "[model][tolerances]") {
  System s = toleranced();
  s.tolerances.overrides = {ToleranceClassOverride{"objective", ToleranceClass::Precise},
                            ToleranceClassOverride{"image", ToleranceClass::None}};
  const ClassTolerances r = class_tolerances(s);
  using C = ClassSkip;
  const ToleranceClass n = ToleranceClass::Normal;
  const ToleranceClass p = ToleranceClass::Precise;
  const std::vector<Row> expected = {
      // stop (element): no material in the beam; its plane surface has no radius
      {"class:node:stop:thickness", C::None, n},
      {"class:node:stop:decenter_x", C::None, n},
      {"class:node:stop:decenter_y", C::None, n},
      {"class:node:stop:tilt_x_deg", C::None, n},
      {"class:node:stop:tilt_y_deg", C::None, n},
      {"class:node:stop:index", C::NoMaterialInBeam, n},
      {"class:node:stop:abbe", C::NoMaterialInBeam, n},
      {"class:surface:STO:curvature_relative", C::Plane, n},
      {"class:surface:STO:conic", C::Plane, n},
      // the assembly, with its class override
      {"class:node:objective:decenter_x", C::None, p},
      {"class:node:objective:decenter_y", C::None, p},
      {"class:node:objective:shift_z", C::None, p},
      {"class:node:objective:tilt_x_deg", C::None, p},
      {"class:node:objective:tilt_y_deg", C::None, p},
      // L1: two segments with their own glasses
      {"class:node:L1:thickness", C::None, p},
      {"class:node:L1:decenter_x", C::None, p},
      {"class:node:L1:decenter_y", C::None, p},
      {"class:node:L1:tilt_x_deg", C::None, p},
      {"class:node:L1:tilt_y_deg", C::None, p},
      {"class:node:L1:index:0", C::None, p},
      {"class:node:L1:index:1", C::None, p},
      {"class:node:L1:abbe:0", C::None, p},
      {"class:node:L1:abbe:1", C::None, p},
      {"class:surface:L1.S1:curvature_relative", C::None, p},
      {"class:surface:L1.S1:conic", C::Sphere, p},
      {"class:surface:L1.S2:curvature_relative", C::None, p},
      {"class:surface:L1.S2:conic", C::Sphere, p},
      {"class:surface:L1.S2:thickness", C::None, p},
      {"class:surface:L1.S3:curvature_relative", C::None, p},
      {"class:surface:L1.S3:conic", C::None, p},
      {"class:surface:L1.S3:thickness", C::None, p},
      // L2: one real CONST material (no Abbe number), an even asphere, a bound thickness
      {"class:node:L2:thickness", C::None, p},
      {"class:node:L2:decenter_x", C::None, p},
      {"class:node:L2:decenter_y", C::None, p},
      {"class:node:L2:tilt_x_deg", C::None, p},
      {"class:node:L2:tilt_y_deg", C::None, p},
      {"class:node:L2:index", C::None, p},
      {"class:node:L2:abbe", C::ConstAbbe, p},
      {"class:surface:L2.S1:curvature_relative", C::None, p},
      {"class:surface:L2.S1:conic", C::None, p},
      {"class:surface:L2.S2:curvature_relative", C::Plane, p},
      {"class:surface:L2.S2:conic", C::Plane, p},
      {"class:surface:L2.S2:thickness", C::Bound, p},
      // the mirror: no material in the beam
      {"class:node:fold:thickness", C::None, n},
      {"class:node:fold:decenter_x", C::None, n},
      {"class:node:fold:decenter_y", C::None, n},
      {"class:node:fold:tilt_x_deg", C::None, n},
      {"class:node:fold:tilt_y_deg", C::None, n},
      {"class:node:fold:index", C::NoMaterialInBeam, n},
      {"class:node:fold:abbe", C::NoMaterialInBeam, n},
      {"class:surface:FOLD:curvature_relative", C::Plane, n},
      {"class:surface:FOLD:conic", C::Plane, n},
      // the detector has the class none: nothing listed
  };
  CHECK(rows(r) == expected);

  const ClassTolerance index1 = find(r, "class:node:L1:index:1");
  CHECK(index1.target == [] {
    ToleranceTarget t = at(A::Node, "L1");
    t.segment = 1;
    return t;
  }());
  CHECK(index1.kind == K::Index);
  CHECK(index1.level == ToleranceLevel::Element);
  CHECK(find(r, "class:node:L2:index").target == at(A::Node, "L2"));  // one material: no segment
  CHECK(find(r, "class:surface:L1.S2:thickness").target == at(A::Surface, "L1.S2"));
  CHECK(find(r, "class:surface:L1.S2:thickness").level == ToleranceLevel::Surface);
  CHECK(find(r, "class:node:objective:shift_z").level == ToleranceLevel::Assembly);
  CHECK(r.value_overrides.empty());
}

TEST_CASE("class: a uniform segment list is one material (ADR 0035, addendum #254)",
          "[model][tolerances]") {
  // A list whose entries are all the same string is one body (ADR 0017): the class generates
  // index and abbe only for segment 0, with the id suffix :0.
  System s = toleranced();
  Element prism{"prism",
                ElementKind::Plate,
                Pose::along_z(30.0),
                std::nullopt,
                {surface("P.S1", Plane{}, 0.0), surface("P.S2", Plane{}, 5.0),
                 surface("P.S3", Plane{}, 10.0), surface("P.S4", Plane{}, 15.0)}};
  prism.segment_materials = {"SCHOTT:N-BK7", "SCHOTT:N-BK7", "SCHOTT:N-BK7"};
  s.root.children.insert(s.root.children.begin() + 2, Node{prism});
  const ClassTolerances r = class_tolerances(s);
  std::vector<std::string> material;
  for (const ClassTolerance& c : r.items) {
    if (c.id.starts_with("class:node:prism:index") || c.id.starts_with("class:node:prism:abbe")) {
      material.push_back(c.id);
    }
  }
  CHECK(material ==
        std::vector<std::string>{"class:node:prism:index:0", "class:node:prism:abbe:0"});
  CHECK(find(r, "class:node:prism:index:0").target.segment == 0U);
  // L1 with two glasses stays segmented.
  CHECK(find(r, "class:node:L1:index:1").skip == ClassSkip::None);
}

TEST_CASE("class: crystals, elements without material and a bound element position",
          "[model][tolerances]") {
  System s = toleranced();
  Element& stop = test::element(s, 0);
  stop.pose.position[2] = Param::bound("T2");
  Element& fold = test::element(s, 2);
  fold.kind = ElementKind::Plate;
  fold.material = std::nullopt;
  fold.crystal = CrystalMaterial{"BIREFRINGENT:CALCITE", "BIREFRINGENT:CALCITE-E"};
  fold.optic_axis = std::array<double, 3>{0.0, 0.0, 1.0};
  Element plate = fold;
  plate.name = "plate";
  plate.crystal = std::nullopt;
  plate.optic_axis = std::nullopt;
  plate.surfaces[0].id = SurfaceId("P.S1");
  Element metal = plate;
  metal.name = "metal";
  metal.material = "CONST:0.2,3.5";
  metal.surfaces[0].id = SurfaceId("M.S1");
  Element vacuum = plate;
  vacuum.name = "vacuum";
  vacuum.material = "VACUUM";
  vacuum.surfaces[0].id = SurfaceId("V.S1");
  s.root.children.insert(s.root.children.begin() + 3, {Node{plate}, Node{metal}, Node{vacuum}});
  const ClassTolerances r = class_tolerances(s);
  CHECK(find(r, "class:node:stop:thickness").skip == ClassSkip::Bound);
  CHECK(find(r, "class:node:fold:index").skip == ClassSkip::Crystal);
  CHECK(find(r, "class:node:fold:abbe").skip == ClassSkip::Crystal);
  CHECK(find(r, "class:node:plate:index").skip == ClassSkip::NoMaterial);
  CHECK(find(r, "class:node:plate:abbe").skip == ClassSkip::NoMaterial);
  CHECK(find(r, "class:node:metal:index").skip == ClassSkip::Metal);
  CHECK(find(r, "class:node:metal:abbe").skip == ClassSkip::Metal);
  CHECK(find(r, "class:node:vacuum:index").skip == ClassSkip::None);
  CHECK(find(r, "class:node:vacuum:abbe").skip == ClassSkip::ConstAbbe);
  // The detector has the system class here.
  CHECK(find(r, "class:node:image:thickness").skip == ClassSkip::Bound);
}

TEST_CASE("class: none, missing class and the nearest ancestor", "[model][tolerances]") {
  System s = toleranced();
  SECTION("no class: nothing listed") {
    s.tolerances.tolerance_class = ToleranceClass::None;
    CHECK(class_tolerances(s).items.empty());
  }
  SECTION("a class only below the root") {
    s.tolerances.tolerance_class = ToleranceClass::None;
    s.tolerances.overrides = {ToleranceClassOverride{"objective", ToleranceClass::Coarse},
                              ToleranceClassOverride{"L2", ToleranceClass::None}};
    const ClassTolerances r = class_tolerances(s);
    REQUIRE_FALSE(r.items.empty());
    for (const ClassTolerance& c : r.items) {
      INFO(c.id);
      CHECK(c.tolerance_class == ToleranceClass::Coarse);
      CHECK(
          (c.id.find(":objective:") != std::string::npos || c.id.find(":L1") != std::string::npos));
    }
  }
  SECTION("the target's own override before its parent's") {
    s.tolerances.overrides = {ToleranceClassOverride{"L1", ToleranceClass::Coarse},
                              ToleranceClassOverride{"objective", ToleranceClass::Precise}};
    const ClassTolerances r = class_tolerances(s);
    CHECK(find(r, "class:node:objective:tilt_x_deg").tolerance_class == ToleranceClass::Precise);
    CHECK(find(r, "class:node:L1:tilt_x_deg").tolerance_class == ToleranceClass::Coarse);
    CHECK(find(r, "class:surface:L1.S2:thickness").tolerance_class == ToleranceClass::Coarse);
    CHECK(find(r, "class:node:L2:tilt_x_deg").tolerance_class == ToleranceClass::Precise);
    CHECK(find(r, "class:node:fold:tilt_x_deg").tolerance_class == ToleranceClass::Normal);
  }
}

TEST_CASE("class: single values and explicit items by perturbed quantity", "[model][tolerances]") {
  System s = toleranced();
  ToleranceValueOverride index1 = value(A::Node, "L1", K::Index);
  index1.target.segment = 1;
  s.tolerances.overrides = {
      ToleranceClassOverride{"image", ToleranceClass::None},              // 0
      value(A::Surface, "L1.S2", K::CurvatureRelative),                   // 1
      index1,                                                             // 2
      value(A::Node, "image", K::Thickness),                              // 3
      value(A::Surface, "L2.S2", K::Thickness),                           // 4
      value(A::Surface, "L1.S1", K::Conic),                               // 5
      value(A::Node, "L1", K::ShiftZ),                                    // 6
      TolerancePivotOverride{at(A::Node, "objective"), {0.0, 0.0, 9.5}},  // 7
      value(A::Surface, "L1.S3", K::Conic),                               // 8
      value(A::Node, "nowhere", K::DecenterX),                            // 9
  };
  ToleranceTarget s3_conic = at(A::Surface, "L1.S3");
  s3_conic.param = "/shape/base/conic";
  ToleranceTarget l2_gap = at(A::Node, "L2");
  l2_gap.param = "/pose/position/2";
  ToleranceTarget l1_abbe0 = at(A::Node, "L1");
  l1_abbe0.segment = 0;
  // L1.S2 through its element: the same Param as the thickness of L1.S2.
  ToleranceTarget via_node = at(A::Node, "L1");
  via_node.param = "/surfaces/1/pose/position/2";
  s.tolerances.items = {
      item("S3-conic", s3_conic, K::Param),                   // 0
      item("S1-radius", at(A::Surface, "L1.S1"), K::Radius),  // 1
      item("gap", l2_gap, K::Param),                          // 2
      item("L1-dx", at(A::Node, "L1"), K::DecenterX),         // 3
      item("glass", l1_abbe0, K::Abbe),                       // 4
      item("fold-roll", at(A::Node, "fold"), K::TiltZDeg),    // 5
      item("via-node", via_node, K::Param),                   // 6
  };
  const ClassTolerances r = class_tolerances(s);

  CHECK(find(r, "class:surface:L1.S2:curvature_relative").value_override == 1U);
  CHECK(find(r, "class:node:L1:index:1").value_override == 2U);
  CHECK_FALSE(find(r, "class:node:L1:index:0").value_override);

  CHECK(find(r, "class:surface:L1.S3:conic").replaced_by == 0U);
  CHECK(find(r, "class:surface:L1.S1:curvature_relative").replaced_by == 1U);
  CHECK(find(r, "class:node:L2:thickness").replaced_by == 2U);
  CHECK(find(r, "class:node:L1:decenter_x").replaced_by == 3U);
  CHECK(find(r, "class:node:L1:abbe:0").replaced_by == 4U);
  CHECK_FALSE(find(r, "class:node:L1:abbe:1").replaced_by);
  CHECK(find(r, "class:surface:L1.S2:thickness").replaced_by == 6U);
  CHECK_FALSE(find(r, "class:node:fold:tilt_x_deg").replaced_by);
  CHECK_FALSE(find(r, "class:surface:L1.S3:conic").generated());
  CHECK(find(r, "class:surface:L1.S3:curvature_relative").generated());

  using E = OverrideEffect;
  const std::vector<OverrideStatus> expected = {
      {1, E::Applied, ClassSkip::None},   {2, E::Applied, ClassSkip::None},
      {3, E::ClassNone, ClassSkip::None}, {4, E::Skipped, ClassSkip::Bound},
      {5, E::Skipped, ClassSkip::Sphere}, {6, E::NotGenerated, ClassSkip::None},
      {8, E::Replaced, ClassSkip::None},  {9, E::NotGenerated, ClassSkip::None},
  };
  CHECK(r.value_overrides == expected);
}

TEST_CASE("class: the result does not depend on the limits or on other items",
          "[model][tolerances]") {
  // Ids and order depend only on the tree, the classes and the overrides (ADR 0035, point 5):
  // an explicit item elsewhere and other limits change nothing but the replaced entry.
  const System a = toleranced();
  System b = a;
  b.tolerances.items = {item("roll", at(A::Node, "fold"), K::TiltZDeg)};
  b.tolerances.overrides = {value(A::Node, "L1", K::Thickness)};
  std::get<ToleranceValueOverride>(b.tolerances.overrides[0]).limit = 0.5;
  const ClassTolerances ra = class_tolerances(a);
  const ClassTolerances rb = class_tolerances(b);
  REQUIRE_FALSE(ra.items.empty());  // content guard: the comparison below is not vacuous
  REQUIRE(ra.items.size() == rb.items.size());
  for (std::size_t i = 0; i < ra.items.size(); ++i) {
    INFO(ra.items[i].id);
    CHECK(ra.items[i].id == rb.items[i].id);
    CHECK(ra.items[i].skip == rb.items[i].skip);
  }
}
