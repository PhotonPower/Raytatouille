// Crystal elements in compile (#131, ADR 0026, points 2-4): one medium per crystal element with
// n_O, n_E and the global optic axis, the mode of each event inside a crystal, and the path
// rules of ADR 0026, point 4 (allowed cases compile, all others give a code at the event).

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cmath>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"

using Catch::Matchers::ContainsSubstring;
using rtt::compile::CompiledEvent;
using rtt::compile::CompiledSystem;
using rtt::compile::CompileError;
using rtt::compile::CrystalMode;
using rtt::compile::PathId;
using rtt::material::Material;
using rtt::material::MaterialLibrary;
using rtt::material::WavelengthRange;
using rtt::math::Complex;
using rtt::math::Vec3;
using rtt::model::CrystalMaterial;
using rtt::model::Element;
using rtt::model::ElementKind;
using rtt::model::Event;
using rtt::model::EventKind;
using rtt::model::Param;
using rtt::model::Pose;
using rtt::model::Surface;
using rtt::model::SurfaceId;
using rtt::model::System;

namespace {

/// Valid system: wavelengths 0.5876 um (reference) and 0.6563 um, EPD, empty root.
System bare_system() {
  System s;
  s.name = "crystal";
  s.wavelengths = {{0.5876, 1.0, true}, {0.6563, 1.0, false}};
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

/// Plate `name` of an n_O = 1.6584, n_E = 1.4864 crystal at z_mm, 5 mm thick, axis (0, 1, 1) in
/// element coordinates; surfaces <name>.S1 and <name>.S2.
Element crystal_plate(const std::string& name, double z_mm) {
  Element e{name,
            ElementKind::Plate,
            Pose::along_z(z_mm),
            std::nullopt,
            {plane(name + ".S1"), plane(name + ".S2", 5.0)}};
  e.crystal = CrystalMaterial{"CONST:1.6584", "CONST:1.4864"};
  e.optic_axis = std::array<double, 3>{0.0, 1.0, 1.0};
  return e;
}

Event event(const std::string& surface, EventKind kind, int order = 0) {
  return {SurfaceId(surface), kind, order};
}

/// The errors of compile(s); fails the test if it compiles.
std::vector<rtt::model::Diagnostic> errors(const System& s, const MaterialLibrary& lib) {
  try {
    (void)rtt::compile::compile(s, lib);
  } catch (const CompileError& e) {
    return e.diagnostics();
  }
  FAIL("compile() did not throw");
  return {};
}

/// True if `d` has exactly one diagnostic, with this code and location.
bool only(const std::vector<rtt::model::Diagnostic>& d,
          const std::string& code,
          const std::string& location) {
  INFO((d.empty() ? std::string("no diagnostics") : rtt::model::to_string(d.front())));
  return d.size() == 1 && d[0].code == code && d[0].location == location;
}

/// Constant index on a wavelength range (um).
class RangedMaterial final : public Material {
 public:
  RangedMaterial(double n, WavelengthRange range) : n_(n), range_(range) {}
  [[nodiscard]] Complex index(double /*wavelength_um*/,
                              double /*temperature_c*/,
                              double /*pressure_atm*/) const override {
    return {n_, 0.0};
  }
  [[nodiscard]] std::optional<WavelengthRange> wavelength_range_um() const override {
    return range_;
  }

 private:
  double n_;
  WavelengthRange range_;
};

}  // namespace

TEST_CASE("crystal plate: one medium with n_O, n_E and the global optic axis", "[crystal]") {
  System s = bare_system();
  Element plate = crystal_plate("P", 10.0);
  plate.pose.rotation_deg[0] = Param(90.0);  // axis (0, 1, 1) -> (0, -1, 1) globally (Rx(90))
  s.root.children.push_back({plate});
  s.paths = {{"e",
              false,
              {event("P.S1", EventKind::Extraordinary), event("P.S2", EventKind::Transmit),
               event("P.S2", EventKind::Refract)}}};
  const MaterialLibrary lib;
  const CompiledSystem cs = rtt::compile::compile(s, lib);

  const auto& events = cs.path(PathId{0}).events;
  REQUIRE(events.size() == 3);
  const std::uint32_t crystal = events[0].medium_after;
  const auto& m = cs.media()[crystal];
  REQUIRE(m.is_crystal());
  REQUIRE(m.reference == "CONST:1.6584");
  REQUIRE(m.reference_extraordinary == "CONST:1.4864");
  REQUIRE(m.index == std::vector<Complex>{{1.6584, 0.0}, {1.6584, 0.0}});
  REQUIRE(m.index_extraordinary == std::vector<double>{1.4864, 1.4864});
  REQUIRE((*m.optic_axis - Vec3(0.0, -1.0, 1.0).normalized()).norm() <= 1e-15);
  REQUIRE(std::abs(m.optic_axis->norm() - 1.0) <= 1e-15);

  // The environment stays isotropic and shared.
  const auto& air = cs.media()[cs.environment_medium()];
  REQUIRE_FALSE(air.is_crystal());
  REQUIRE(air.reference_extraordinary.empty());
  REQUIRE(air.index_extraordinary.empty());

  // Entry with a mode, the mode kept through a Transmit inside, exit with Refract.
  REQUIRE(events[0].medium_before == cs.environment_medium());
  REQUIRE(events[0].medium_beyond == crystal);
  REQUIRE_FALSE(events[0].from_inside);
  REQUIRE(events[0].crystal_mode == CrystalMode::None);
  REQUIRE(events[1].medium_before == crystal);
  REQUIRE(events[1].medium_after == crystal);
  REQUIRE(events[1].from_inside);
  REQUIRE(events[1].crystal_mode == CrystalMode::Extraordinary);
  REQUIRE(events[2].medium_before == crystal);
  REQUIRE(events[2].medium_after == cs.environment_medium());
  REQUIRE(events[2].crystal_mode == CrystalMode::Extraordinary);
}

TEST_CASE("crystal: the ordinary mode is recorded as well", "[crystal]") {
  System s = bare_system();
  s.root.children.push_back({crystal_plate("P", 10.0)});
  s.paths = {{"o", false, {event("P.S1", EventKind::Ordinary), event("P.S2", EventKind::Refract)}}};
  const MaterialLibrary lib;
  const CompiledSystem cs = rtt::compile::compile(s, lib);
  const auto& events = cs.path(PathId{0}).events;
  REQUIRE(events[0].crystal_mode == CrystalMode::None);
  REQUIRE(events[1].crystal_mode == CrystalMode::Ordinary);
}

TEST_CASE("crystal: two elements of the same crystal get two media with their own axes",
          "[crystal]") {
  System s = bare_system();
  Element a = crystal_plate("A", 10.0);
  Element b = crystal_plate("B", 20.0);
  b.optic_axis = std::array<double, 3>{1.0, 0.0, 0.0};
  s.root.children.push_back({a});
  s.root.children.push_back({b});
  s.paths = {{"both",
              false,
              {event("A.S1", EventKind::Ordinary), event("A.S2", EventKind::Refract),
               event("B.S1", EventKind::Extraordinary), event("B.S2", EventKind::Refract)}}};
  const MaterialLibrary lib;
  const CompiledSystem cs = rtt::compile::compile(s, lib);
  const auto& events = cs.path(PathId{0}).events;
  const auto& ma = cs.media()[events[0].medium_after];
  const auto& mb = cs.media()[events[2].medium_after];
  REQUIRE(events[0].medium_after != events[2].medium_after);
  REQUIRE((*ma.optic_axis - Vec3(0.0, 1.0, 1.0).normalized()).norm() <= 1e-15);
  REQUIRE((*mb.optic_axis - Vec3(1.0, 0.0, 0.0)).norm() <= 1e-15);
  REQUIRE(events[3].crystal_mode == CrystalMode::Extraordinary);
}

TEST_CASE("crystal: allowed events outside the crystal (ADR 0026, point 4)", "[crystal]") {
  const MaterialLibrary lib;
  System s = bare_system();
  Element plate = crystal_plate("P", 10.0);
  plate.surfaces[0].interaction = rtt::model::IdealAntiReflection{};
  s.root.children.push_back({plate});
  // Transmit at a crystal surface from outside: a dummy passage, the medium stays the
  // environment. Reflect from outside at an ideal_anti_reflection surface: r = 0.
  s.paths = {
      {"outside", false, {event("P.S1", EventKind::Transmit), event("P.S1", EventKind::Reflect)}}};
  const CompiledSystem cs = rtt::compile::compile(s, lib);
  for (const CompiledEvent& e : cs.path(PathId{0}).events) {
    REQUIRE(e.medium_before == cs.environment_medium());
    REQUIRE(e.medium_after == cs.environment_medium());
    REQUIRE(e.crystal_mode == CrystalMode::None);
  }
}

TEST_CASE("crystal: path rules give a code at the event (ADR 0026, point 4)", "[crystal]") {
  const MaterialLibrary lib;
  System s = bare_system();
  s.root.children.push_back({crystal_plate("P", 10.0)});
  Element glass{"G",
                ElementKind::Plate,
                Pose::along_z(30.0),
                "CONST:1.5",
                {plane("G.S1"), plane("G.S2", 5.0)}};
  s.root.children.push_back({glass});

  SECTION("refract into a crystal: a mode is required") {
    s.paths = {
        {"x", false, {event("P.S1", EventKind::Refract), event("P.S2", EventKind::Refract)}}};
    REQUIRE(only(errors(s, lib), "paths.crystal_mode_required", "/paths/0/events/0"));
  }
  SECTION("a mode at an element without crystal") {
    s.paths = {
        {"x", false, {event("G.S1", EventKind::Ordinary), event("G.S2", EventKind::Refract)}}};
    REQUIRE(only(errors(s, lib), "paths.mode_without_crystal", "/paths/0/events/0"));
  }
  SECTION("a mode at the exit of a crystal") {
    s.paths = {{"x",
                false,
                {event("P.S1", EventKind::Ordinary), event("P.S2", EventKind::Extraordinary)}}};
    REQUIRE(only(errors(s, lib), "paths.mode_without_crystal", "/paths/0/events/1"));
  }
  SECTION("reflection inside a crystal") {
    s.paths = {
        {"x", false, {event("P.S1", EventKind::Ordinary), event("P.S2", EventKind::Reflect)}}};
    REQUIRE(only(errors(s, lib), "crystal.unsupported", "/paths/0/events/1"));
  }
  SECTION("reflection from outside at a Fresnel crystal surface") {
    s.paths = {{"x", false, {event("P.S1", EventKind::Reflect)}}};
    REQUIRE(only(errors(s, lib), "crystal.unsupported", "/paths/0/events/0"));
  }
  SECTION("transmit inside a crystal with an order") {
    std::get<Element>(s.root.children[0].value)
        .surfaces[1]
        .phases.emplace_back(rtt::model::LinearGrating{Param(300.0), 0.0});
    s.paths = {{"x",
                false,
                {event("P.S1", EventKind::Ordinary), event("P.S2", EventKind::Transmit, 1),
                 event("P.S2", EventKind::Refract)}}};
    REQUIRE(only(errors(s, lib), "crystal.unsupported", "/paths/0/events/1"));
  }
  SECTION("the automatic path through a crystal element") {
    s.paths = {{"main", true, {}}};
    REQUIRE(only(errors(s, lib), "crystal.unsupported", "/paths/0/events"));
  }
  SECTION("crystal to crystal: the inner surface of a crystal lens") {
    Element lens{"L",
                 ElementKind::Lens,
                 Pose::along_z(50.0),
                 std::nullopt,
                 {plane("L.S1"), plane("L.S2", 2.0), plane("L.S3", 4.0)}};
    lens.crystal = CrystalMaterial{"CONST:1.6584", "CONST:1.4864"};
    lens.optic_axis = std::array<double, 3>{0.0, 0.0, 1.0};
    s.root.children.push_back({lens});
    s.paths = {{"x",
                false,
                {event("L.S1", EventKind::Ordinary), event("L.S2", EventKind::Refract),
                 event("L.S3", EventKind::Refract)}}};
    REQUIRE(only(errors(s, lib), "crystal.unsupported", "/paths/0/events/1"));
  }
}

TEST_CASE("crystal: interactions other than fresnel and ideal_anti_reflection", "[crystal]") {
  const MaterialLibrary lib;
  System s = bare_system();
  Element plate = crystal_plate("P", 10.0);
  plate.surfaces[1].interaction = rtt::model::IdealMirror{};
  s.root.children.push_back({plate});
  s.paths = {{"x", false, {event("P.S1", EventKind::Ordinary), event("P.S2", EventKind::Refract)}}};
  REQUIRE(only(errors(s, lib), "crystal.interaction_unsupported",
               "/root/children/0/surfaces/1/interaction"));
}

TEST_CASE("crystal: an absorbing part is an error at its reference", "[crystal]") {
  const MaterialLibrary lib;
  System s = bare_system();
  Element plate = crystal_plate("P", 10.0);
  s.paths = {{"x", false, {event("P.S1", EventKind::Ordinary), event("P.S2", EventKind::Refract)}}};

  SECTION("ordinary") {
    plate.crystal->ordinary = "CONST:1.6584,0.001";
    s.root.children.push_back({plate});
    const auto d = errors(s, lib);
    REQUIRE(only(d, "crystal.absorbing", "/root/children/0/material/ordinary"));
    REQUIRE_THAT(d[0].message, ContainsSubstring("0.5876"));
  }
  SECTION("extraordinary") {
    plate.crystal->extraordinary = "CONST:1.4864,0.002";
    s.root.children.push_back({plate});
    REQUIRE(only(errors(s, lib), "crystal.absorbing", "/root/children/0/material/extraordinary"));
  }
}

TEST_CASE("crystal: unknown parts are reported at their reference", "[crystal]") {
  const MaterialLibrary lib;
  System s = bare_system();
  Element plate = crystal_plate("P", 10.0);
  plate.crystal->extraordinary = "NOPE:X-E";
  s.root.children.push_back({plate});
  s.paths = {{"x", false, {event("P.S1", EventKind::Ordinary), event("P.S2", EventKind::Refract)}}};
  REQUIRE(only(errors(s, lib), "material.unknown", "/root/children/0/material/extraordinary"));
}

TEST_CASE("crystal: an unknown environment gives no follow-up errors at crystals", "[crystal]") {
  // With /environment/medium unresolved, the environment index is a placeholder that may alias
  // the crystal's medium; the path rules are skipped, so only the real error remains.
  const MaterialLibrary lib;
  System s = bare_system();
  s.environment.medium = "NOPE:AIR";
  s.root.children.push_back({crystal_plate("P", 10.0)});
  s.paths = {{"x", false, {event("P.S1", EventKind::Ordinary), event("P.S2", EventKind::Refract)}}};
  REQUIRE(only(errors(s, lib), "material.unknown", "/environment/medium"));
}

TEST_CASE("crystal: entering a second crystal from inside the first is crystal to crystal",
          "[crystal]") {
  // Elements do not nest: entering B while inside A leaves A (rules of #5), so both sides are
  // crystals; with a mode as well as with Refract.
  const MaterialLibrary lib;
  System s = bare_system();
  s.root.children.push_back({crystal_plate("A", 10.0)});
  s.root.children.push_back({crystal_plate("B", 20.0)});
  for (const EventKind kind : {EventKind::Extraordinary, EventKind::Refract}) {
    s.paths = {{"x",
                false,
                {event("A.S1", EventKind::Ordinary), event("B.S1", kind),
                 event("B.S2", EventKind::Refract)}}};
    REQUIRE(only(errors(s, lib), "crystal.unsupported", "/paths/0/events/1"));
  }
}

TEST_CASE("crystal: the optic axis turns with the assembly and element poses", "[crystal]") {
  // As the axes of the ideal elements (ADR 0021): assembly pose, then element pose; not the
  // surface pose. Assembly Rz(90): (0, -1, 1) of the element pose Rx(90) -> (1, 0, 1).
  System s = bare_system();
  rtt::model::Assembly group;
  group.name = "group";
  group.pose.rotation_deg[2] = Param(90.0);
  Element plate = crystal_plate("P", 10.0);
  plate.pose.rotation_deg[0] = Param(90.0);
  plate.surfaces[0].pose.rotation_deg[1] = Param(30.0);  // surface pose: no effect on the axis
  group.children.push_back({plate});
  s.root.children.push_back({group});
  s.paths = {{"o", false, {event("P.S1", EventKind::Ordinary), event("P.S2", EventKind::Refract)}}};
  const MaterialLibrary lib;
  const CompiledSystem cs = rtt::compile::compile(s, lib);
  const auto& m = cs.media()[cs.path(PathId{0}).events[0].medium_after];
  REQUIRE((*m.optic_axis - Vec3(1.0, 0.0, 1.0).normalized()).norm() <= 1e-15);
}

TEST_CASE("crystal: wavelengths outside the common range of both parts", "[crystal]") {
  MaterialLibrary lib;
  lib.add("TEST:O", std::make_shared<const RangedMaterial>(1.66, WavelengthRange{0.4, 1.0}));
  lib.add("TEST:E", std::make_shared<const RangedMaterial>(1.49, WavelengthRange{0.6, 2.0}));
  lib.add("TEST:LOW", std::make_shared<const RangedMaterial>(1.49, WavelengthRange{0.3, 0.35}));
  System s = bare_system();  // 0.5876 and 0.6563 um
  Element plate = crystal_plate("P", 10.0);
  s.paths = {{"x", false, {event("P.S1", EventKind::Ordinary), event("P.S2", EventKind::Refract)}}};

  SECTION("outside the intersection") {
    plate.crystal = CrystalMaterial{"TEST:O", "TEST:E"};
    s.root.children.push_back({plate});
    const auto d = errors(s, lib);
    REQUIRE(only(d, "material.wavelength_out_of_range", "/root/children/0/material"));
    REQUIRE_THAT(d[0].message, ContainsSubstring("[0.6, 1] um"));
  }
  SECTION("disjoint ranges: no common range") {
    plate.crystal = CrystalMaterial{"TEST:O", "TEST:LOW"};
    s.root.children.push_back({plate});
    const auto d = errors(s, lib);
    REQUIRE(only(d, "material.wavelength_out_of_range", "/root/children/0/material"));
    REQUIRE_THAT(d[0].message, ContainsSubstring("no common valid range"));
  }
}
