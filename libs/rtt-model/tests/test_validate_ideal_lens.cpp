// validate() for the ideal lens and the ideal cylinder lens (#178, ADR 0031, points 1 and 7):
// only at the single plane surface of a thin_element without shape terms or phases; f and
// object_distance non-zero and finite; axis_deg finite.

#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <limits>
#include <optional>
#include <string>
#include <vector>

#include "rtt/model/validate.hpp"
#include "test_support.hpp"

using rtt::model::Diagnostic;
using rtt::model::Element;
using rtt::model::ElementKind;
using rtt::model::has_errors;
using rtt::model::IdealCylinderLens;
using rtt::model::IdealLens;
using rtt::model::Param;
using rtt::model::Severity;
using rtt::model::System;
using rtt::model::validate;
using rtt::model::test::element;
using rtt::model::test::make_singlet;

namespace {

constexpr double kInf = std::numeric_limits<double>::infinity();
const std::string kLens = "/root/children/2/surfaces/0/interaction";

/// The singlet with a thin element at z = 50 (child 2) whose only surface carries `interaction`.
System with_thin(const rtt::model::Interaction& interaction) {
  System s = make_singlet();
  Element thin;
  thin.name = "IL";
  thin.kind = ElementKind::ThinElement;
  thin.pose.position[2] = Param(50.0);
  rtt::model::Surface surface;
  surface.id = rtt::model::SurfaceId{"IL.S1"};
  surface.aperture = rtt::model::CircularAperture{10.0, 0.0};
  surface.interaction = interaction;
  thin.surfaces = {surface};
  s.root.children.insert(s.root.children.begin() + 2, rtt::model::Node{thin});
  return s;
}

/// True if `d` has an error with this code at this location.
bool has(const std::vector<Diagnostic>& d, const std::string& code, const std::string& location) {
  for (const Diagnostic& x : d) {
    if (x.severity == Severity::Error && x.code == code && x.location == location) return true;
  }
  return false;
}

}  // namespace

TEST_CASE("ideal lenses at the plane surface of a thin element are valid", "[validate][ideal]") {
  CHECK_FALSE(has_errors(validate(with_thin(IdealLens{Param(50.0), std::nullopt}))));
  CHECK_FALSE(has_errors(validate(with_thin(IdealLens{Param(-50.0), Param(-30.0)}))));
  CHECK_FALSE(has_errors(validate(with_thin(IdealCylinderLens{Param(50.0), 30.0, Param(75.0)}))));
}

TEST_CASE("interaction.ideal_lens_not_allowed", "[validate][ideal]") {
  SECTION("at a lens surface") {
    System s = make_singlet();
    element(s, 1).surfaces[0].interaction = IdealLens{Param(50.0), std::nullopt};
    CHECK(has(validate(s), "interaction.ideal_lens_not_allowed",
              "/root/children/1/surfaces/0/interaction"));
  }
  SECTION("at a curved surface of a thin element") {
    System s = with_thin(IdealLens{Param(50.0), std::nullopt});
    element(s, 2).surfaces[0].shape.base = rtt::model::Conic{Param(100.0), Param(0.0)};
    CHECK(has(validate(s), "interaction.ideal_lens_not_allowed", kLens));
  }
  SECTION("with shape terms") {
    System s = with_thin(IdealCylinderLens{Param(50.0), 0.0, std::nullopt});
    rtt::model::ZernikeSag z;
    z.normalization_radius = Param(10.0);
    z.coefficients = {Param(0.0), Param(0.0), Param(0.0), Param(0.001)};
    element(s, 2).surfaces[0].shape.terms = {z};
    CHECK(has(validate(s), "interaction.ideal_lens_not_allowed", kLens));
  }
  SECTION("with a phase layer") {
    System s = with_thin(IdealLens{Param(50.0), std::nullopt});
    rtt::model::LinearGrating g;
    g.lines_per_mm = Param(500.0);
    element(s, 2).surfaces[0].phases = {g};
    CHECK(has(validate(s), "interaction.ideal_lens_not_allowed", kLens));
  }
}

TEST_CASE("interaction.focal_length_invalid", "[validate][ideal]") {
  for (const double f : {0.0, kInf, -kInf}) {
    INFO("f = " << f);
    CHECK(has(validate(with_thin(IdealLens{Param(f), std::nullopt})),
              "interaction.focal_length_invalid", kLens + "/focal_length"));
    CHECK(has(validate(with_thin(IdealCylinderLens{Param(f), 0.0, std::nullopt})),
              "interaction.focal_length_invalid", kLens + "/focal_length"));
  }
  // A bound focal length is checked by compile after the resolution (ADR 0029).
  System s = with_thin(IdealLens{Param::bound("F"), std::nullopt});
  rtt::model::ParameterRow row;
  row.name = "F";
  row.form = 0.0;
  s.parameters = {row};
  CHECK_FALSE(has(validate(s), "interaction.focal_length_invalid", kLens + "/focal_length"));
}

TEST_CASE("interaction.object_distance_invalid", "[validate][ideal]") {
  for (const double d : {0.0, kInf}) {
    INFO("object_distance = " << d);
    CHECK(has(validate(with_thin(IdealLens{Param(50.0), Param(d)})),
              "interaction.object_distance_invalid", kLens + "/object_distance"));
    CHECK(has(validate(with_thin(IdealCylinderLens{Param(50.0), 0.0, Param(d)})),
              "interaction.object_distance_invalid", kLens + "/object_distance"));
  }
}

TEST_CASE("a non-finite axis_deg is value.not_finite", "[validate][ideal]") {
  const System s = with_thin(IdealCylinderLens{Param(50.0), std::nan(""), std::nullopt});
  CHECK(has(validate(s), "value.not_finite", kLens + "/axis_deg"));
}
