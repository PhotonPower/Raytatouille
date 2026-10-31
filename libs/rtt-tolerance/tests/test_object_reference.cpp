// object.reference in tolerance trials (#288; ADR 0038, point 5): the object stays at the
// nominal vertex of its reference surface. Neither a model offset of that surface (thickness,
// ADR 0036, point 1a) nor the overlay (decenter, point 1b) moves it; a tolerance on the distance
// itself does.

#include <bit>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

#include "rtt/coating/catalog.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
#include "rtt/model/tolerance.hpp"
#include "rtt/tolerance/trial.hpp"

using rtt::model::Param;
using rtt::model::System;
using rtt::model::ToleranceAnchor;
using rtt::model::ToleranceItem;
using rtt::model::ToleranceKind;
using rtt::tolerance::TrialContext;

namespace {

bool same(double a, double b) {
  return std::bit_cast<std::uint64_t>(a) == std::bit_cast<std::uint64_t>(b);
}

/// tests/reference/r3/object_reference (ideal lens at z = 20, object 75 mm in front of it, so
/// z_obj = -55) with one criterion, ray_y of field 1 on the image surface.
System referenced() {
  System s = rtt::io::load_system(std::filesystem::path(RTT_REFERENCE_DIR) / "r3" /
                                  "object_reference.rtt.json");
  rtt::model::RayOperand ray;
  ray.path = "main";
  ray.surface = rtt::model::SurfaceId("IMG");
  ray.field = 1;
  rtt::model::Criterion c;
  c.id = "image-y";
  c.operand = ray;
  c.max_change = 1.0;
  s.tolerances.criteria = {c};
  return s;
}

ToleranceItem item(const std::string& id,
                   ToleranceAnchor anchor,
                   const std::string& name,
                   ToleranceKind kind,
                   double limit) {
  ToleranceItem i;
  i.id = id;
  i.target.anchor = anchor;
  i.target.name = name;
  i.kind = kind;
  i.limit = limit;
  return i;
}

/// object_z of the trial with `values`, compiled as evaluate() does (model offsets and overlay).
double trial_object_z(const TrialContext& ctx, const std::vector<double>& values) {
  const rtt::tolerance::Applied a = rtt::tolerance::apply(ctx, values);
  const rtt::material::MaterialLibrary materials;
  const rtt::coating::CoatingLibrary coatings;
  return rtt::compile::compile(a.system, materials, coatings, 0, a.perturbation).object_z();
}

}  // namespace

TEST_CASE("object.reference: a trial keeps the object at the nominal reference",
          "[tolerance][trial][object-reference]") {
  const System s = referenced();
  const rtt::material::MaterialLibrary materials;
  const std::vector<ToleranceItem> items = {
      item("lens-dt", ToleranceAnchor::Node, "ideal lens", ToleranceKind::Thickness, 1.0),
      item("lens-dy", ToleranceAnchor::Node, "ideal lens", ToleranceKind::DecenterY, 1.0)};
  const TrialContext ctx = rtt::tolerance::prepare(s, materials, nullptr, items);
  REQUIRE(ctx.nominal().status == rtt::tolerance::TrialStatus::Valid);
  CHECK(same(trial_object_z(ctx, {0.0, 0.0}), -55.0));
  CHECK(same(trial_object_z(ctx, {0.5, 0.0}), -55.0));   // the lens 0.5 mm further in the model
  CHECK(same(trial_object_z(ctx, {0.0, 0.25}), -55.0));  // decentred in the overlay
  // Content guard: the lens did move, so the image ray changes.
  const rtt::tolerance::TrialResult r =
      rtt::tolerance::evaluate(ctx, std::vector<double>{0.5, 0.0});
  REQUIRE(r.status == rtt::tolerance::TrialStatus::Valid);
  REQUIRE(r.criteria.at(0).value.has_value());
  REQUIRE(ctx.nominal().criteria.at(0).value.has_value());
  CHECK(*r.criteria[0].value != *ctx.nominal().criteria[0].value);
}

TEST_CASE("object.reference: a tolerance on the distance moves the object",
          "[tolerance][trial][object-reference]") {
  // The distance bound to a row (75), a parameter item on the row: z_obj = 20 - (75 + v).
  System s = referenced();
  rtt::model::ParameterRow row;
  row.name = "OBJ";
  row.form = 75.0;
  s.parameters = {row};
  s.object.distance = Param::bound("OBJ");
  const rtt::material::MaterialLibrary materials;
  const std::vector<ToleranceItem> items = {
      item("obj", ToleranceAnchor::Parameter, "OBJ", ToleranceKind::Parameter, 1.0)};
  const TrialContext ctx = rtt::tolerance::prepare(s, materials, nullptr, items);
  CHECK(same(trial_object_z(ctx, {0.0}), -55.0));
  CHECK(same(trial_object_z(ctx, {0.5}), 20.0 - 75.5));
}

TEST_CASE("object.reference: trials give the same criteria as the system from the origin",
          "[tolerance][trial][object-reference]") {
  // The same system written once with object.reference (75 mm in front of IL at z = 20) and once
  // from the origin (55 mm): both put the object at z = -55, nominally and in every trial, so all
  // criteria agree bitwise. efl and bfl of a finite object go through the variant with the object
  // at infinity (addendum #256), which must compile with and without object.reference.
  // Field 0 only: the variant has one field, and validate checks every criterion against it.
  System referenced_s = referenced();
  rtt::model::RayOperand ray;
  ray.path = "main";
  ray.surface = rtt::model::SurfaceId("IMG");
  ray.py = 0.5;
  rtt::model::Criterion rc;
  rc.id = "image-y";
  rc.operand = ray;
  rc.max_change = 1.0;
  referenced_s.tolerances.criteria = {rc};
  for (const auto q : {rtt::model::FirstOrderQuantity::Efl, rtt::model::FirstOrderQuantity::Bfl}) {
    rtt::model::FirstOrderOperand f;
    f.quantity = q;
    f.path = "main";
    rtt::model::Criterion c;
    c.id = q == rtt::model::FirstOrderQuantity::Efl ? "efl" : "bfl";
    c.operand = f;
    c.max_change = 1.0;
    referenced_s.tolerances.criteria.push_back(c);
  }
  System origin_s = referenced_s;
  origin_s.object.reference.reset();
  origin_s.object.distance = Param(55.0);

  const rtt::material::MaterialLibrary materials;
  const std::vector<ToleranceItem> items = {
      item("lens-dt", ToleranceAnchor::Node, "ideal lens", ToleranceKind::Thickness, 1.0),
      item("lens-dy", ToleranceAnchor::Node, "ideal lens", ToleranceKind::DecenterY, 1.0)};
  const TrialContext with_ref = rtt::tolerance::prepare(referenced_s, materials, nullptr, items);
  const TrialContext from_origin = rtt::tolerance::prepare(origin_s, materials, nullptr, items);
  for (const std::vector<double>& values :
       {std::vector<double>{0.0, 0.0}, std::vector<double>{0.5, 0.0},
        std::vector<double>{0.0, 0.25}}) {
    const rtt::tolerance::TrialResult a = rtt::tolerance::evaluate(with_ref, values);
    const rtt::tolerance::TrialResult b = rtt::tolerance::evaluate(from_origin, values);
    REQUIRE(a.status == rtt::tolerance::TrialStatus::Valid);
    REQUIRE(b.status == rtt::tolerance::TrialStatus::Valid);
    REQUIRE(a.criteria.size() == 3);
    REQUIRE(b.criteria.size() == 3);
    for (std::size_t i = 0; i < a.criteria.size(); ++i) {
      REQUIRE(a.criteria[i].value.has_value());
      REQUIRE(b.criteria[i].value.has_value());
      CHECK(same(*a.criteria[i].value, *b.criteria[i].value));
    }
  }
  // Content guard: the thickness moves the lens, so the image ray changes; efl and bfl of the
  // single ideal lens stay f = 50 (the variant really compiled and gave values).
  const rtt::tolerance::TrialResult moved =
      rtt::tolerance::evaluate(with_ref, std::vector<double>{0.5, 0.0});
  const rtt::tolerance::TrialResult nominal = with_ref.nominal();
  REQUIRE(moved.criteria.at(0).value.has_value());
  REQUIRE(nominal.criteria.at(0).value.has_value());
  CHECK(*moved.criteria[0].value != *nominal.criteria[0].value);
  for (std::size_t i = 1; i < 3; ++i) {
    REQUIRE(moved.criteria.at(i).value.has_value());
    CHECK_THAT(*moved.criteria[i].value, Catch::Matchers::WithinRel(50.0, 1e-12));
  }
}
