// compile() of the ideal lens and the ideal cylinder lens (#178, ADR 0031, points 1 and 7): the
// resolved lens data of the compiled surface, the path rule paths.ideal_lens_event, a lens inside
// a crystal, and focal length and object distance bound to the parameter table (ADR 0029).

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <numbers>
#include <optional>
#include <string>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/compile/errors.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"

using rtt::compile::CompiledSystem;
using rtt::compile::CompileError;
using rtt::material::MaterialLibrary;
using rtt::model::Element;
using rtt::model::ElementKind;
using rtt::model::Event;
using rtt::model::EventKind;
using rtt::model::IdealCylinderLens;
using rtt::model::IdealLens;
using rtt::model::Param;
using rtt::model::Pose;
using rtt::model::Surface;
using rtt::model::SurfaceId;
using rtt::model::System;

namespace {

System load(const std::string& relative) {
  return rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/" + relative);
}

/// The surface `id` of `cs`.
const rtt::compile::CompiledSurface& surface(const CompiledSystem& cs, const std::string& id) {
  return cs.surfaces()[*cs.find_surface(SurfaceId{id})];
}

/// Valid system: one wavelength, EPD 10 mm, an empty root, the automatic path.
System bare_system() {
  System s;
  s.name = "ideal lens";
  s.wavelengths = {{0.5876, 1.0, true}};
  s.aperture = {rtt::model::SystemApertureType::EntrancePupilDiameter, Param(10.0)};
  s.fields = {rtt::model::FieldType::AngleDeg, {{0.0, 0.0, 1.0}}};
  s.root.name = "root";
  s.paths = {{"main", true, {}}};
  return s;
}

Surface plane(const std::string& id, double z_mm = 0.0) {
  Surface s;
  s.id = SurfaceId(id);
  s.pose = Pose::along_z(z_mm);
  return s;
}

/// Thin element `name` at z_mm whose only surface `name` carries `interaction`.
Element thin(const std::string& name, double z_mm, const rtt::model::Interaction& interaction) {
  Surface s = plane(name);
  s.aperture = rtt::model::CircularAperture{10.0, 0.0};
  s.interaction = interaction;
  return Element{name, ElementKind::ThinElement, Pose::along_z(z_mm), std::nullopt, {s}};
}

/// Stop at 0, the thin element "IL" at 10 with `interaction`, the detector at 60.
System lens_system(const rtt::model::Interaction& interaction) {
  System s = bare_system();
  Surface sto = plane("STO");
  sto.aperture = rtt::model::CircularAperture{5.0, 0.0};
  s.root.children.push_back({Element{"stop", ElementKind::Stop, {}, std::nullopt, {sto}}});
  s.root.children.push_back({thin("IL", 10.0, interaction)});
  s.root.children.push_back(
      {Element{"image", ElementKind::Detector, Pose::along_z(60.0), std::nullopt, {plane("IMG")}}});
  return s;
}

Event event(const std::string& surface, EventKind kind) {
  return {SurfaceId(surface), kind, 0};
}

/// The errors of compile(s, configuration); fails the test if it compiles.
std::vector<rtt::model::Diagnostic> errors(const System& s, std::size_t configuration = 0) {
  const MaterialLibrary lib;
  try {
    (void)rtt::compile::compile(s, lib, configuration);
  } catch (const CompileError& e) {
    return e.diagnostics();
  }
  FAIL("compile() did not throw");
  return {};
}

bool only(const std::vector<rtt::model::Diagnostic>& d,
          const std::string& code,
          const std::string& location) {
  INFO((d.empty() ? std::string("no diagnostics") : rtt::model::to_string(d.front())));
  return d.size() == 1 && d[0].code == code && d[0].location == location;
}

/// A table with the configurations a and b and the row `name` with the values `a`, `b`.
void table(System& s, const std::string& name, double a, double b) {
  s.configurations = {{"a"}, {"b"}};
  rtt::model::ParameterRow row;
  row.name = name;
  row.form = std::vector<double>{a, b};
  s.parameters.push_back(row);
}

}  // namespace

TEST_CASE("ideal lens: the compiled surface carries f, t_o and no power axis", "[ideal]") {
  const CompiledSystem cs =
      rtt::compile::compile(load("r2/ideal_lens.rtt.json"), MaterialLibrary{});
  const auto& il = surface(cs, "IL").ideal_lens;
  REQUIRE(il.has_value());
  CHECK(il->focal_length == 50.0);
  // t_o = 1/s with s = -object_distance (ADR 0031, point 3); the same expression as compile.
  CHECK(il->t_o == 1.0 / -75.0);
  CHECK_FALSE(il->power_axis.has_value());
  // A surface without an ideal lens has none.
  CHECK_FALSE(surface(cs, "STO").ideal_lens.has_value());
}

TEST_CASE("ideal lens: an object at infinity gives t_o = 0", "[ideal]") {
  const CompiledSystem cs =
      rtt::compile::compile(lens_system(IdealLens{Param(-40.0), std::nullopt}), MaterialLibrary{});
  const auto& il = surface(cs, "IL").ideal_lens;
  REQUIRE(il.has_value());
  CHECK(il->focal_length == -40.0);
  CHECK(il->t_o == 0.0);
}

TEST_CASE("ideal cylinder lens: the power axis b = (-sin psi, cos psi, 0)", "[ideal]") {
  SECTION("psi = 0: b along y (the reference file)") {
    const CompiledSystem cs =
        rtt::compile::compile(load("r2/ideal_lens.rtt.json"), MaterialLibrary{});
    const auto& cl = surface(cs, "CL").ideal_lens;
    REQUIRE(cl.has_value());
    REQUIRE(cl->power_axis.has_value());
    CHECK(cl->power_axis->x() == 0.0);
    CHECK(cl->power_axis->y() == 1.0);
    CHECK(cl->power_axis->z() == 0.0);
  }
  SECTION("psi = 30 deg, a unit vector") {
    const CompiledSystem cs = rtt::compile::compile(
        lens_system(IdealCylinderLens{Param(50.0), 30.0, std::nullopt}), MaterialLibrary{});
    const auto& b = surface(cs, "IL").ideal_lens->power_axis;
    REQUIRE(b.has_value());
    const double psi = 30.0 * std::numbers::pi / 180.0;
    CHECK(std::abs(b->x() + std::sin(psi)) <= 1e-16);
    CHECK(std::abs(b->y() - std::cos(psi)) <= 1e-16);
    CHECK(b->z() == 0.0);
    CHECK(std::abs(b->norm() - 1.0) <= 2e-16);
  }
}

TEST_CASE("paths.ideal_lens_event: only transmit at an ideal lens", "[ideal]") {
  for (const EventKind kind : {EventKind::Refract, EventKind::Reflect}) {
    for (const bool cylinder : {false, true}) {
      INFO("kind " << static_cast<int>(kind) << ", cylinder " << cylinder);
      System s = cylinder ? lens_system(IdealCylinderLens{Param(50.0), 0.0, std::nullopt})
                          : lens_system(IdealLens{Param(50.0), std::nullopt});
      s.paths = {{"x",
                  false,
                  {event("STO", EventKind::Transmit), event("IL", kind),
                   event("IMG", EventKind::Transmit)}}};
      CHECK(only(errors(s), "paths.ideal_lens_event", "/paths/0/events/1"));
    }
  }
  // transmit is fine, also in an explicit path.
  System ok = lens_system(IdealLens{Param(50.0), std::nullopt});
  ok.paths = {{"x",
               false,
               {event("STO", EventKind::Transmit), event("IL", EventKind::Transmit),
                event("IMG", EventKind::Transmit)}}};
  CHECK_NOTHROW(rtt::compile::compile(ok, MaterialLibrary{}));
}

TEST_CASE("an ideal lens inside a crystal is crystal.unsupported", "[ideal]") {
  System s = bare_system();
  Element plate{"P",
                ElementKind::Plate,
                Pose::along_z(10.0),
                std::nullopt,
                {plane("P.S1"), plane("P.S2", 10.0)}};
  plate.crystal = rtt::model::CrystalMaterial{"CONST:1.6584", "CONST:1.4864"};
  plate.optic_axis = std::array<double, 3>{0.0, 1.0, 1.0};
  s.root.children.push_back({plate});
  s.root.children.push_back({thin("IL", 15.0, IdealLens{Param(50.0), std::nullopt})});
  s.paths = {{"x",
              false,
              {event("P.S1", EventKind::Ordinary), event("IL", EventKind::Transmit),
               event("P.S2", EventKind::Refract)}}};
  CHECK(only(errors(s), "crystal.unsupported", "/paths/0/events/1"));
}

TEST_CASE("bound focal length and object distance per configuration (ADR 0029)", "[ideal]") {
  System s = lens_system(IdealLens{Param::bound("F"), Param::bound("OBJ")});
  table(s, "F", 50.0, 80.0);
  table(s, "OBJ", 75.0, 120.0);
  const MaterialLibrary lib;
  const CompiledSystem a = rtt::compile::compile(s, lib, std::size_t{0});
  const CompiledSystem b = rtt::compile::compile(s, lib, std::size_t{1});
  CHECK(surface(a, "IL").ideal_lens->focal_length == 50.0);
  CHECK(surface(b, "IL").ideal_lens->focal_length == 80.0);
  CHECK(surface(a, "IL").ideal_lens->t_o == 1.0 / -75.0);
  CHECK(surface(b, "IL").ideal_lens->t_o == 1.0 / -120.0);
}

TEST_CASE("a bound value 0 is an error at the Param after the resolution", "[ideal]") {
  const std::string at = "/root/children/1/surfaces/0/interaction";
  SECTION("focal length") {
    System s = lens_system(IdealLens{Param::bound("F"), std::nullopt});
    table(s, "F", 50.0, 0.0);
    CHECK_NOTHROW(rtt::compile::compile(s, MaterialLibrary{}, std::size_t{0}));
    CHECK(only(errors(s, 1), "interaction.focal_length_invalid", at + "/focal_length"));
  }
  SECTION("object distance") {
    System s = lens_system(IdealCylinderLens{Param(50.0), 0.0, Param::bound("OBJ")});
    table(s, "OBJ", 0.0, 75.0);
    CHECK(only(errors(s, 0), "interaction.object_distance_invalid", at + "/object_distance"));
    CHECK_NOTHROW(rtt::compile::compile(s, MaterialLibrary{}, std::size_t{1}));
  }
}
