#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cmath>
#include <memory>
#include <optional>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/air.hpp"
#include "test_support.hpp"

using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinAbs;
using rtt::compile::compile;
using rtt::compile::CompiledSystem;
using rtt::compile::CompileError;
using rtt::material::MaterialLibrary;
using rtt::math::Complex;
using rtt::math::Vec3;
using rtt::model::Element;
using rtt::model::ElementKind;
using rtt::model::EventKind;
using rtt::model::Param;
using rtt::model::Pose;
using rtt::model::Surface;
using rtt::model::SurfaceId;
using rtt::model::System;

namespace {

constexpr double kTol = 1e-12;  // mm, issue #5

void require_near(const Vec3& actual, const Vec3& expected, double tol = kTol) {
  INFO("actual " << actual.transpose() << ", expected " << expected.transpose());
  REQUIRE((actual - expected).cwiseAbs().maxCoeff() <= tol);
}

System load(const std::string& relative) {
  return rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/" + relative);
}

Surface plane_surface(const std::string& id, double z_mm = 0.0) {
  Surface s;
  s.id = SurfaceId(id);
  s.pose = Pose::along_z(z_mm);
  s.aperture = rtt::model::CircularAperture{10.0, 0.0};
  return s;
}

/// Valid system with one wavelength, one field, EPD aperture, empty root and no path.
System bare_system() {
  System s;
  s.name = "test";
  s.wavelengths = {{0.5876, 1.0, true}};
  s.aperture = {rtt::model::SystemApertureType::EntrancePupilDiameter, Param(10.0)};
  s.fields = {rtt::model::FieldType::AngleDeg, {{0.0, 0.0, 1.0}}};
  s.root.name = "root";
  return s;
}

std::uint32_t surface_index(const CompiledSystem& cs, const std::string& id) {
  const auto i = cs.find_surface(SurfaceId(id));
  REQUIRE(i.has_value());
  return *i;
}

std::uint32_t medium_index(const CompiledSystem& cs, const std::string& reference) {
  const auto& m = cs.media();
  const auto it =
      std::find_if(m.begin(), m.end(), [&](const auto& x) { return x.reference == reference; });
  REQUIRE(it != m.end());
  return static_cast<std::uint32_t>(it - m.begin());
}

struct ExpectedEvent {
  std::string surface;
  EventKind kind;
  std::string before;  // medium reference
  std::string after;
};

void require_events(const CompiledSystem& cs,
                    const std::string& path_name,
                    const std::vector<ExpectedEvent>& expected) {
  const auto id = cs.find_path(path_name);
  REQUIRE(id.has_value());
  const auto& events = cs.path(*id).events;
  REQUIRE(events.size() == expected.size());
  for (std::size_t i = 0; i < expected.size(); ++i) {
    INFO("event " << i << " at " << expected[i].surface);
    REQUIRE(events[i].surface == surface_index(cs, expected[i].surface));
    REQUIRE(events[i].kind == expected[i].kind);
    REQUIRE(events[i].medium_before == medium_index(cs, expected[i].before));
    REQUIRE(events[i].medium_after == medium_index(cs, expected[i].after));
  }
}

/// Runs compile and returns the CompileError, failing the test if none is thrown.
CompileError compile_error(const System& s) {
  const MaterialLibrary lib;
  try {
    (void)compile(s, lib);
  } catch (const CompileError& e) {
    return e;
  }
  FAIL("compile did not throw CompileError");
  return CompileError({});
}

bool has_error_at(const CompileError& e, const std::string& location) {
  return rtt::model::test::has_error_at(e.diagnostics(), location);
}

}  // namespace

TEST_CASE("global surface pose chains root, assemblies, element and surface", "[compile]") {
  // Hand calculation with the pose convention p_parent = position + pivot + R (p_child - pivot),
  // R = Rx Ry Rz (docs/architecture.md, Konventionen):
  //   Rz(90 deg) (x, y, z) = (-y, x, z),  Rx(90 deg) (x, y, z) = (x, -z, y)
  //   surface:  pose along z by 2
  //             p_e = p_s + (0, 0, 2)
  //   element:  position (0, 0, 4), Rz(90 deg)
  //             p_a = (0, 0, 4) + Rz p_e
  //   assembly: position (1, 2, 3), Rx(90 deg), pivot (0, 0, 5)
  //             p_r = (1, 2, 3) + (0, 0, 5) + Rx (p_a - (0, 0, 5))
  //   root:     position (0, 0, 10)
  //             p_g = p_r + (0, 0, 10)
  // Local (1, 0, 0): p_e = (1, 0, 2), p_a = (0, 1, 6), p_r = (1, 2, 8) + (0, -1, 1) = (1, 1, 9),
  //                  p_g = (1, 1, 19).
  // Local origin:    p_e = (0, 0, 2), p_a = (0, 0, 6), p_r = (1, 2, 8) + (0, -1, 0) = (1, 1, 8),
  //                  p_g = (1, 1, 18).
  // Local +z axis:   Rx Rz (0, 0, 1) = Rx (0, 0, 1) = (0, -1, 0).
  System s = bare_system();
  s.root.pose.position[2] = Param(10.0);

  Element thin{"T", ElementKind::ThinElement, {}, std::nullopt, {plane_surface("T.S", 2.0)}};
  thin.pose.position[2] = Param(4.0);
  thin.pose.rotation_deg[2] = Param(90.0);

  rtt::model::Assembly group;
  group.name = "group";
  group.pose.position = {Param(1.0), Param(2.0), Param(3.0)};
  group.pose.rotation_deg[0] = Param(90.0);
  group.pose.pivot = {0.0, 0.0, 5.0};
  group.children.push_back({thin});
  s.root.children.push_back({group});
  s.paths = {{"main", true, {}}};

  const MaterialLibrary lib;
  const CompiledSystem cs = compile(s, lib);
  REQUIRE(cs.surfaces().size() == 1);
  const auto& surf = cs.surfaces()[0];
  REQUIRE(surf.id == SurfaceId("T.S"));
  REQUIRE(surf.element_name == "T");
  REQUIRE(surf.element_kind == ElementKind::ThinElement);
  require_near(surf.to_global.apply_point(Vec3(1.0, 0.0, 0.0)), Vec3(1.0, 1.0, 19.0));
  require_near(surf.to_global.apply_point(Vec3::Zero()), Vec3(1.0, 1.0, 18.0));
  require_near(surf.to_global.apply_vector(Vec3::UnitZ()), Vec3(0.0, -1.0, 0.0));
  // to_local is the inverse.
  require_near(surf.to_local.apply_point(Vec3(1.0, 1.0, 19.0)), Vec3(1.0, 0.0, 0.0));
  require_near(surf.to_local.apply_vector(Vec3(0.0, -1.0, 0.0)), Vec3::UnitZ());
}

TEST_CASE("Michelson reference system compiles with all media in air", "[compile]") {
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(load("m0/michelson.rtt.json"), lib);

  REQUIRE(cs.wavelengths_um() == std::vector<double>{0.6328});
  REQUIRE(cs.reference_wavelength() == 0);
  REQUIRE(cs.temperature_c() == 20.0);
  REQUIRE(cs.media().size() == 1);
  REQUIRE(cs.media()[cs.environment_medium()].reference == "AIR");
  // AIR is Ciddor air at 20 degC and 1 atm since #25 (was n = 1 in M1).
  REQUIRE(cs.media()[0].index ==
          std::vector<Complex>{Complex(rtt::material::ciddor_air_index(0.6328, 20.0, 1.0))});

  REQUIRE(cs.surfaces().size() == 4);
  REQUIRE(surface_index(cs, "BS") == 0);
  REQUIRE(surface_index(cs, "M1") == 1);
  REQUIRE(surface_index(cs, "M2") == 2);
  REQUIRE(surface_index(cs, "DET") == 3);

  // Beam splitter at z = 50 tilted by Rx(45 deg): normal Rx (0, 0, 1) = (0, -sin 45, cos 45).
  const auto& bs = cs.surfaces()[0];
  const double h = std::sqrt(0.5);
  require_near(bs.to_global.apply_point(Vec3::Zero()), Vec3(0.0, 0.0, 50.0));
  require_near(bs.to_global.apply_vector(Vec3::UnitZ()), Vec3(0.0, -h, h));
  REQUIRE(std::holds_alternative<rtt::model::IdealBeamSplitter>(bs.interaction));
  // Test mirror at (0, 50, 50) with Rx(-90 deg): normal (0, 1, 0).
  const auto& m2 = cs.surfaces()[2];
  require_near(m2.to_global.apply_point(Vec3::Zero()), Vec3(0.0, 50.0, 50.0));
  require_near(m2.to_global.apply_vector(Vec3::UnitZ()), Vec3(0.0, 1.0, 0.0));
  REQUIRE(std::holds_alternative<rtt::model::IdealMirror>(m2.interaction));
  REQUIRE(m2.aperture == rtt::model::Aperture{rtt::model::CircularAperture{12.5, 0.0}});
  // Detector at (0, -50, 50) with Rx(-90 deg): normal (0, 1, 0). Light arrives along -y, i.e.
  // from the local +z side; surfaces can be hit from both sides.
  const auto& det = cs.surfaces()[3];
  require_near(det.to_global.apply_point(Vec3::Zero()), Vec3(0.0, -50.0, 50.0));
  require_near(det.to_global.apply_vector(Vec3::UnitZ()), Vec3(0.0, 1.0, 0.0));
  REQUIRE(det.aperture == rtt::model::Aperture{rtt::model::RectangularAperture{4.0, 3.0}});

  REQUIRE(cs.paths().size() == 2);
  require_events(cs, "reference arm",
                 {{"BS", EventKind::Transmit, "AIR", "AIR"},
                  {"M1", EventKind::Reflect, "AIR", "AIR"},
                  {"BS", EventKind::Reflect, "AIR", "AIR"},
                  {"DET", EventKind::Transmit, "AIR", "AIR"}});
  require_events(cs, "test arm",
                 {{"BS", EventKind::Reflect, "AIR", "AIR"},
                  {"M2", EventKind::Reflect, "AIR", "AIR"},
                  {"BS", EventKind::Transmit, "AIR", "AIR"},
                  {"DET", EventKind::Transmit, "AIR", "AIR"}});
}

TEST_CASE("wavelength weights are copied unchanged from the model", "[compile]") {
  // Needed for polychromatic analyses (#28); normalisation is left to the consumer.
  System s = bare_system();
  s.wavelengths = {{0.4861, 1.0, false}, {0.5876, 3.0, true}, {0.6563, 0.5, false}};
  s.root.children.push_back(
      {Element{"D", ElementKind::Detector, {}, std::nullopt, {plane_surface("D.S")}}});
  s.paths = {{"main", true, {}}};
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(s, lib);
  REQUIRE(cs.wavelength_weights() == std::vector<double>{1.0, 3.0, 0.5});
  REQUIRE(cs.wavelength_weights().size() == cs.wavelengths_um().size());
}

TEST_CASE("singlet with CONST:1.5168 goes air -> glass -> air", "[compile]") {
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(load("m1/singlet_const.rtt.json"), lib);

  REQUIRE(cs.wavelengths_um() == std::vector<double>{0.4861, 0.5876, 0.6563});
  REQUIRE(cs.reference_wavelength() == 1);
  REQUIRE(cs.media().size() == 2);
  const auto glass = medium_index(cs, "CONST:1.5168");
  REQUIRE(cs.media()[glass].index == std::vector<Complex>(3, Complex(1.5168, 0.0)));
  // AIR is Ciddor air at 20 degC and 1 atm since #25 (was n = 1 in M1).
  for (std::size_t i = 0; i < 3; ++i) {
    REQUIRE(cs.media()[cs.environment_medium()].index[i] ==
            Complex(rtt::material::ciddor_air_index(cs.wavelengths_um()[i], 20.0, 1.0)));
  }

  // Vertices on the axis: stop at 0, lens at 5 with thickness 4, image at 106.363 mm.
  const std::vector<std::pair<std::string, double>> vertices{
      {"STO", 0.0}, {"L1.S1", 5.0}, {"L1.S2", 9.0}, {"IMG", 106.363}};
  REQUIRE(cs.surfaces().size() == vertices.size());
  for (std::size_t i = 0; i < vertices.size(); ++i) {
    const auto& surf = cs.surfaces()[i];
    INFO(vertices[i].first);
    REQUIRE(surf.id == SurfaceId(vertices[i].first));
    require_near(surf.to_global.apply_point(Vec3::Zero()), Vec3(0.0, 0.0, vertices[i].second));
    require_near(surf.to_global.apply_vector(Vec3::UnitZ()), Vec3::UnitZ());
  }

  // Shapes: L1.S1 conic with c = 1/R, R = 51.68 mm, k = 0; the others are planes.
  const auto* s1 = std::get_if<rtt::geom::Conic<double>>(&cs.surfaces()[1].shape);
  REQUIRE(s1 != nullptr);
  REQUIRE(s1->curvature() == 1.0 / 51.68);
  REQUIRE(s1->conic_constant() == 0.0);
  REQUIRE(std::holds_alternative<rtt::geom::Plane<double>>(cs.surfaces()[2].shape));
  REQUIRE(cs.surfaces()[1].aperture ==
          rtt::model::Aperture{rtt::model::CircularAperture{12.7, 0.0}});
  REQUIRE_FALSE(cs.surfaces()[3].aperture.has_value());

  require_events(cs, "main",
                 {{"STO", EventKind::Transmit, "AIR", "AIR"},
                  {"L1.S1", EventKind::Refract, "AIR", "CONST:1.5168"},
                  {"L1.S2", EventKind::Refract, "CONST:1.5168", "AIR"},
                  {"IMG", EventKind::Transmit, "AIR", "AIR"}});

  REQUIRE(cs.find_path("main") == rtt::compile::PathId{0});
  REQUIRE_FALSE(cs.find_path("missing").has_value());
  REQUIRE_FALSE(cs.find_surface(SurfaceId("missing")).has_value());
}

TEST_CASE("automatic path chooses the event from the element kind", "[compile]") {
  System s = bare_system();
  s.root.children.push_back(
      {Element{"STOP", ElementKind::Stop, {}, std::nullopt, {plane_surface("STO")}}});
  s.root.children.push_back({Element{"P",
                                     ElementKind::Plate,
                                     Pose::along_z(10.0),
                                     "CONST:1.5",
                                     {plane_surface("P.S1"), plane_surface("P.S2", 5.0)}}});
  s.root.children.push_back({Element{
      "M", ElementKind::Mirror, Pose::along_z(30.0), std::nullopt, {plane_surface("M.S")}}});
  s.root.children.push_back({Element{
      "T", ElementKind::ThinElement, Pose::along_z(20.0), std::nullopt, {plane_surface("T.S")}}});
  s.root.children.push_back({Element{
      "D", ElementKind::Detector, Pose::along_z(5.0), std::nullopt, {plane_surface("D.S")}}});
  s.paths = {{"auto", true, {}}};

  const MaterialLibrary lib;
  const CompiledSystem cs = compile(s, lib);
  require_events(cs, "auto",
                 {{"STO", EventKind::Transmit, "AIR", "AIR"},
                  {"P.S1", EventKind::Refract, "AIR", "CONST:1.5"},
                  {"P.S2", EventKind::Refract, "CONST:1.5", "AIR"},
                  {"M.S", EventKind::Reflect, "AIR", "AIR"},
                  {"T.S", EventKind::Transmit, "AIR", "AIR"},
                  {"D.S", EventKind::Transmit, "AIR", "AIR"}});
}

TEST_CASE("media along explicit paths follow the inside/outside rule", "[compile]") {
  // Rule decided for #5: start in the environment; Refract/Ordinary/Extraordinary at an element
  // with a material toggle inside <-> environment; Reflect, Transmit (with any order, ADR 0025)
  // and events at elements without material keep the medium.
  System s = bare_system();
  s.environment.medium = "CONST:1.333";
  // A phase layer at P.S2, so that the Transmit with order 1 below is valid (ADR 0025).
  Surface grating = plane_surface("P.S2", 5.0);
  const rtt::model::LinearGrating lines{Param(100.0), 0.0};
  grating.phases.emplace_back(lines);
  s.root.children.push_back({Element{"P",
                                     ElementKind::Plate,
                                     Pose::along_z(10.0),
                                     "CONST:1.5",
                                     {plane_surface("P.S1"), grating}}});
  // Mangin mirror: refracting front, reflecting back inside a substrate.
  Surface back = plane_surface("M.S2", 3.0);
  back.interaction = rtt::model::IdealMirror{};
  s.root.children.push_back({Element{
      "M", ElementKind::Mirror, Pose::along_z(30.0), "CONST:1.6", {plane_surface("M.S1"), back}}});
  s.root.children.push_back({Element{
      "T", ElementKind::ThinElement, Pose::along_z(50.0), std::nullopt, {plane_surface("T.S")}}});
  s.paths = {{"explicit",
              false,
              {{SurfaceId("P.S1"), EventKind::Refract, 0},
               {SurfaceId("P.S2"), EventKind::Transmit, 1},
               {SurfaceId("P.S2"), EventKind::Refract, 0},
               {SurfaceId("M.S1"), EventKind::Refract, 0},
               {SurfaceId("M.S2"), EventKind::Reflect, 0},
               {SurfaceId("M.S1"), EventKind::Refract, 0},
               {SurfaceId("T.S"), EventKind::Refract, 0},
               {SurfaceId("T.S"), EventKind::Transmit, 0},
               {SurfaceId("P.S1"), EventKind::Reflect, 0}}}};

  const MaterialLibrary lib;
  const CompiledSystem cs = compile(s, lib);
  REQUIRE(cs.media()[cs.environment_medium()].reference == "CONST:1.333");
  require_events(cs, "explicit",
                 {{"P.S1", EventKind::Refract, "CONST:1.333", "CONST:1.5"},
                  {"P.S2", EventKind::Transmit, "CONST:1.5", "CONST:1.5"},
                  {"P.S2", EventKind::Refract, "CONST:1.5", "CONST:1.333"},
                  {"M.S1", EventKind::Refract, "CONST:1.333", "CONST:1.6"},
                  {"M.S2", EventKind::Reflect, "CONST:1.6", "CONST:1.6"},
                  {"M.S1", EventKind::Refract, "CONST:1.6", "CONST:1.333"},
                  {"T.S", EventKind::Refract, "CONST:1.333", "CONST:1.333"},
                  {"T.S", EventKind::Transmit, "CONST:1.333", "CONST:1.333"},
                  {"P.S1", EventKind::Reflect, "CONST:1.333", "CONST:1.333"}});
  REQUIRE(cs.path(*cs.find_path("explicit")).events[1].order == 1);
}

TEST_CASE("elements with the same material share one medium entry", "[compile]") {
  System s = bare_system();
  s.root.children.push_back({Element{"P1",
                                     ElementKind::Plate,
                                     {},
                                     "CONST:1.5",
                                     {plane_surface("A.S1"), plane_surface("A.S2", 1.0)}}});
  s.root.children.push_back({Element{"P2",
                                     ElementKind::Plate,
                                     Pose::along_z(5.0),
                                     "CONST:1.5",
                                     {plane_surface("B.S1"), plane_surface("B.S2", 1.0)}}});
  s.paths = {{"main", true, {}}};
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(s, lib);
  REQUIRE(cs.media().size() == 2);
  require_events(cs, "main",
                 {{"A.S1", EventKind::Refract, "AIR", "CONST:1.5"},
                  {"A.S2", EventKind::Refract, "CONST:1.5", "AIR"},
                  {"B.S1", EventKind::Refract, "AIR", "CONST:1.5"},
                  {"B.S2", EventKind::Refract, "CONST:1.5", "AIR"}});
}

TEST_CASE("invalid models throw CompileError with the validation diagnostics", "[compile]") {
  System s = load("m1/singlet_const.rtt.json");
  s.wavelengths.clear();
  const CompileError e = compile_error(s);
  REQUIRE(has_error_at(e, "/wavelengths"));
  REQUIRE_THAT(e.what(), ContainsSubstring("/wavelengths"));
}

TEST_CASE("unknown materials throw CompileError pointing at the material", "[compile]") {
  SECTION("element material") {
    // The M0 singlet uses SCHOTT:N-BK7; no catalogue is loaded here.
    const CompileError e = compile_error(load("m0/singlet.rtt.json"));
    REQUIRE(has_error_at(e, "/root/children/1/material"));
    REQUIRE_THAT(e.what(), ContainsSubstring("SCHOTT:N-BK7"));
    REQUIRE_THAT(e.what(), ContainsSubstring("not loaded"));
  }
  SECTION("environment medium") {
    System s = load("m1/singlet_const.rtt.json");
    s.environment.medium = "WATER";
    const CompileError e = compile_error(s);
    REQUIRE(has_error_at(e, "/environment/medium"));
    REQUIRE_THAT(e.what(), ContainsSubstring("WATER"));
  }
  SECTION("malformed constant index") {
    System s = load("m1/singlet_const.rtt.json");
    std::get<Element>(s.root.children[1].value).material = "CONST:1.5,-0.1";
    REQUIRE(has_error_at(compile_error(s), "/root/children/1/material"));
  }
}

TEST_CASE("model features beyond M1 throw CompileError", "[compile]") {
  System s = load("m1/singlet_const.rtt.json");
  auto& lens = std::get<Element>(s.root.children[1].value);

  SECTION("Zernike sag term (M8)") {
    lens.surfaces[0].shape.terms.push_back(rtt::model::ZernikeSag{Param(10.0), {Param(1e-3)}});
    const CompileError e = compile_error(s);
    REQUIRE(has_error_at(e, "/root/children/1/surfaces/0/shape/terms/0"));
  }
}

TEST_CASE("even asphere compiles to geom::EvenAsphere with c = 1/R", "[compile]") {
  // Supported since #6; coefficients as in the model (coefficients[0] = A4).
  System s = load("m1/singlet_const.rtt.json");
  auto& lens = std::get<Element>(s.root.children[1].value);
  lens.surfaces[0].shape.base =
      rtt::model::EvenAsphere{Param(51.68), Param(-0.5), {Param(1e-6), Param(-2e-9)}};
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(s, lib);
  const auto& shape = cs.surfaces()[surface_index(cs, "L1.S1")].shape;
  const auto* asphere = std::get_if<rtt::geom::EvenAsphere<double>>(&shape);
  REQUIRE(asphere != nullptr);
  REQUIRE(asphere->base_conic() == std::pair{1.0 / 51.68, -0.5});
  REQUIRE(asphere->coefficients() == std::vector<double>{1e-6, -2e-9});
}

TEST_CASE("mirror with substrate material and several surfaces needs an explicit path",
          "[compile]") {
  // A Mangin mirror refracts at its front surface; the automatic path would reflect at every
  // surface of a Mirror, which is physically wrong (#6).
  System s = bare_system();
  Surface back = plane_surface("M.S2", 3.0);
  back.interaction = rtt::model::IdealMirror{};
  s.root.children.push_back({Element{
      "M", ElementKind::Mirror, Pose::along_z(30.0), "CONST:1.6", {plane_surface("M.S1"), back}}});
  s.paths = {{"auto", true, {}}};
  const CompileError e = compile_error(s);
  REQUIRE(has_error_at(e, "/root/children/0/surfaces"));
  REQUIRE_THAT(e.what(), ContainsSubstring("mirror with substrate material on the automatic path"));
  REQUIRE_THAT(e.what(), ContainsSubstring("explicit path"));

  SECTION("a single-surface mirror on a substrate stays on the automatic path") {
    System one = bare_system();
    one.root.children.push_back({Element{
        "M", ElementKind::Mirror, Pose::along_z(30.0), "CONST:1.6", {plane_surface("M.S1")}}});
    one.paths = {{"auto", true, {}}};
    const MaterialLibrary lib;
    const CompiledSystem cs = compile(one, lib);
    require_events(cs, "auto", {{"M.S1", EventKind::Reflect, "AIR", "AIR"}});
  }
}

TEST_CASE("plates with more than 2 surfaces are fine on explicit paths", "[compile]") {
  System s = bare_system();
  s.root.children.push_back(
      {Element{"cube",
               ElementKind::Plate,
               {},
               "CONST:1.5",
               {plane_surface("C.S1"), plane_surface("C.S2", 1.0), plane_surface("C.S3", 2.0)}}});
  s.paths = {{"through",
              false,
              {{SurfaceId("C.S1"), EventKind::Refract, 0},
               {SurfaceId("C.S2"), EventKind::Reflect, 0},
               {SurfaceId("C.S3"), EventKind::Refract, 0}}}};
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(s, lib);
  require_events(cs, "through",
                 {{"C.S1", EventKind::Refract, "AIR", "CONST:1.5"},
                  {"C.S2", EventKind::Reflect, "CONST:1.5", "CONST:1.5"},
                  {"C.S3", EventKind::Refract, "CONST:1.5", "AIR"}});
}

TEST_CASE("compiled system is independent of the model and const-only", "[compile]") {
  STATIC_REQUIRE_FALSE(std::is_default_constructible_v<CompiledSystem>);
  STATIC_REQUIRE(
      std::is_const_v<
          std::remove_reference_t<decltype(std::declval<const CompiledSystem&>().surfaces())>>);
  STATIC_REQUIRE(std::is_const_v<
                 std::remove_reference_t<decltype(std::declval<CompiledSystem&>().surfaces())>>);
  STATIC_REQUIRE(
      std::is_const_v<std::remove_reference_t<decltype(std::declval<CompiledSystem&>().paths())>>);
  STATIC_REQUIRE(
      std::is_const_v<std::remove_reference_t<decltype(std::declval<CompiledSystem&>().media())>>);

  std::optional<CompiledSystem> cs;
  {
    auto model = std::make_unique<System>(load("m1/singlet_const.rtt.json"));
    const auto lib = std::make_unique<MaterialLibrary>();
    cs.emplace(compile(*model, *lib));
    // Changing and then destroying the model and the library must not affect the result.
    std::get<Element>(model->root.children[1].value).pose.position[2] = Param(99.0);
    model->paths.clear();
  }
  require_near(cs->surfaces()[1].to_global.apply_point(Vec3::Zero()), Vec3(0.0, 0.0, 5.0));
  REQUIRE(cs->paths().size() == 1);
  REQUIRE(cs->media()[medium_index(*cs, "CONST:1.5168")].index[0] == Complex(1.5168, 0.0));
}

TEST_CASE("interim state until #163/#165: relative poses and bound Params are rejected",
          "[compile]") {
  // Schema 0.4 reads them (#162), compile evaluates them only with #163 (relative placement)
  // and #165 (parameter table). Until then compile rejects them instead of using the parent as
  // reference or the meaningless value 0 of a bound Param. These cases go with #163 and #165.
  const MaterialLibrary lib;
  System s = load("m1/singlet_const.rtt.json");
  std::get<Element>(s.root.children[1].value).pose.reference =
      rtt::model::PoseReference::RelativeToPreceding;
  REQUIRE(has_error_at(compile_error(s), "/root/children/1/pose/reference"));
  REQUIRE_THAT(compile_error(s).what(), ContainsSubstring("pose.reference_unsupported"));

  s = load("m1/singlet_const.rtt.json");
  rtt::model::ParameterRow d;
  d.name = "D";
  d.form = 5.0;
  s.parameters = {d};
  std::get<Element>(s.root.children[1].value).pose.position[2] = Param::bound("D");
  REQUIRE(has_error_at(compile_error(s), "/root/children/1/pose/position/2"));
  REQUIRE_THAT(compile_error(s).what(), ContainsSubstring("param.unresolved"));

  // rotate_first needs no evaluation beyond to_isometry: it compiles already.
  s = load("m1/singlet_const.rtt.json");
  auto& lens = std::get<Element>(s.root.children[1].value);
  lens.pose.order = rtt::model::PoseOrder::RotateFirst;
  const CompiledSystem cs = compile(s, lib);
  REQUIRE_THAT(cs.surfaces()[1].to_global.translation().z(), WithinAbs(5.0, kTol));
}

// ---------------------------------------------------------------------------------------------
// Cemented groups (ADR 0017, #27): segment i lies between surface i and surface i + 1.
// ---------------------------------------------------------------------------------------------

namespace {

/// Lens "L" at z = 10 with plane surfaces L.S1, L.S2, ... 2 mm apart and the given segment list.
Element segmented_lens(std::size_t surface_count, std::vector<std::string> segment_materials) {
  Element lens{"L", ElementKind::Lens, Pose::along_z(10.0), std::nullopt, {}};
  for (std::size_t i = 0; i < surface_count; ++i) {
    lens.surfaces.push_back(
        plane_surface("L.S" + std::to_string(i + 1), 2.0 * static_cast<double>(i)));
  }
  lens.segment_materials = std::move(segment_materials);
  return lens;
}

rtt::model::Event refract(const std::string& id) {
  return {SurfaceId(id), EventKind::Refract, 0};
}

}  // namespace

TEST_CASE("cemented doublet on the automatic path goes air -> glass A -> glass B -> air",
          "[compile][cemented]") {
  System s = bare_system();
  s.wavelengths = {{0.4861, 1.0, false}, {0.5876, 1.0, true}};
  s.root.children.push_back({segmented_lens(3, {"CONST:1.5168", "CONST:1.62"})});
  s.root.children.push_back({Element{
      "D", ElementKind::Detector, Pose::along_z(50.0), std::nullopt, {plane_surface("D.S")}}});
  s.paths = {{"auto", true, {}}};

  const MaterialLibrary lib;
  const CompiledSystem cs = compile(s, lib);
  REQUIRE(cs.media().size() == 3);
  // Every segment material is evaluated once per system wavelength (ADR 0014).
  REQUIRE(cs.media()[medium_index(cs, "CONST:1.5168")].index ==
          std::vector<Complex>(2, Complex(1.5168, 0.0)));
  REQUIRE(cs.media()[medium_index(cs, "CONST:1.62")].index ==
          std::vector<Complex>(2, Complex(1.62, 0.0)));
  require_events(cs, "auto",
                 {{"L.S1", EventKind::Refract, "AIR", "CONST:1.5168"},
                  {"L.S2", EventKind::Refract, "CONST:1.5168", "CONST:1.62"},
                  {"L.S3", EventKind::Refract, "CONST:1.62", "AIR"},
                  {"D.S", EventKind::Transmit, "AIR", "AIR"}});
}

TEST_CASE("cemented doublet in reverse and in double pass", "[compile][cemented]") {
  System s = bare_system();
  s.root.children.push_back({segmented_lens(3, {"CONST:1.5", "CONST:1.7"})});
  s.root.children.push_back({Element{
      "M", ElementKind::Mirror, Pose::along_z(30.0), std::nullopt, {plane_surface("M.S")}}});
  s.paths = {{"reverse", false, {refract("L.S3"), refract("L.S2"), refract("L.S1")}},
             {"double pass",
              false,
              {refract("L.S1"),
               refract("L.S2"),
               refract("L.S3"),
               {SurfaceId("M.S"), EventKind::Reflect, 0},
               refract("L.S3"),
               refract("L.S2"),
               refract("L.S1")}}};

  const MaterialLibrary lib;
  const CompiledSystem cs = compile(s, lib);
  // Entering through the last surface goes into the last segment; Refract at surface i from
  // segment i goes back into segment i - 1, at the first surface into the environment.
  require_events(cs, "reverse",
                 {{"L.S3", EventKind::Refract, "AIR", "CONST:1.7"},
                  {"L.S2", EventKind::Refract, "CONST:1.7", "CONST:1.5"},
                  {"L.S1", EventKind::Refract, "CONST:1.5", "AIR"}});
  require_events(cs, "double pass",
                 {{"L.S1", EventKind::Refract, "AIR", "CONST:1.5"},
                  {"L.S2", EventKind::Refract, "CONST:1.5", "CONST:1.7"},
                  {"L.S3", EventKind::Refract, "CONST:1.7", "AIR"},
                  {"M.S", EventKind::Reflect, "AIR", "AIR"},
                  {"L.S3", EventKind::Refract, "AIR", "CONST:1.7"},
                  {"L.S2", EventKind::Refract, "CONST:1.7", "CONST:1.5"},
                  {"L.S1", EventKind::Refract, "CONST:1.5", "AIR"}});
}

TEST_CASE("cemented group of two equal glasses has the media of a single lens",
          "[compile][cemented]") {
  const MaterialLibrary lib;
  System single = bare_system();
  single.root.children.push_back({Element{"L",
                                          ElementKind::Lens,
                                          Pose::along_z(10.0),
                                          "CONST:1.5",
                                          {plane_surface("L.S1"), plane_surface("L.S3", 4.0)}}});
  single.paths = {{"auto", true, {}}};
  const CompiledSystem one = compile(single, lib);

  System list = bare_system();
  list.root.children.push_back({segmented_lens(3, {"CONST:1.5", "CONST:1.5"})});
  list.paths = {{"auto", true, {}}};
  const CompiledSystem two = compile(list, lib);

  System shorthand = bare_system();
  Element lens = segmented_lens(3, {});
  lens.material = "CONST:1.5";  // the shorthand holds for all segments
  shorthand.root.children.push_back({lens});
  shorthand.paths = {{"auto", true, {}}};
  const CompiledSystem three = compile(shorthand, lib);

  require_events(one, "auto",
                 {{"L.S1", EventKind::Refract, "AIR", "CONST:1.5"},
                  {"L.S3", EventKind::Refract, "CONST:1.5", "AIR"}});
  REQUIRE(one.media().size() == 2);
  for (const CompiledSystem* cs : {&two, &three}) {
    REQUIRE(cs->media().size() == one.media().size());
    for (std::size_t m = 0; m < one.media().size(); ++m) {
      REQUIRE(cs->media()[m].reference == one.media()[m].reference);
      REQUIRE(cs->media()[m].index == one.media()[m].index);
    }
    require_events(*cs, "auto",
                   {{"L.S1", EventKind::Refract, "AIR", "CONST:1.5"},
                    {"L.S2", EventKind::Refract, "CONST:1.5", "CONST:1.5"},
                    {"L.S3", EventKind::Refract, "CONST:1.5", "AIR"}});
  }
}

TEST_CASE("two-surface lens with a one-entry material list has that material",
          "[compile][cemented]") {
  // Closes the gap left by #26: rtt-compile used to read only the shorthand.
  System s = bare_system();
  s.root.children.push_back({segmented_lens(2, {"CONST:1.6"})});
  s.paths = {{"auto", true, {}}};
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(s, lib);
  REQUIRE(cs.media().size() == 2);
  require_events(cs, "auto",
                 {{"L.S1", EventKind::Refract, "AIR", "CONST:1.6"},
                  {"L.S2", EventKind::Refract, "CONST:1.6", "AIR"}});
}

TEST_CASE("plate with different segment materials follows its segments on any path",
          "[compile][cemented]") {
  System s = bare_system();
  s.environment.medium = "CONST:1.333";
  Element plate{"P",
                ElementKind::Plate,
                Pose::along_z(5.0),
                std::nullopt,
                {plane_surface("P.S1"), plane_surface("P.S2", 1.0), plane_surface("P.S3", 2.0),
                 plane_surface("P.S4", 3.0)}};
  plate.segment_materials = {"CONST:1.5", "CONST:1.6", "CONST:1.5"};
  s.root.children.push_back({plate});
  s.paths = {
      {"auto", true, {}},
      {"explicit", false, {refract("P.S4"), refract("P.S3"), refract("P.S2"), refract("P.S1")}}};
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(s, lib);
  REQUIRE(cs.media().size() == 3);
  require_events(cs, "auto",
                 {{"P.S1", EventKind::Refract, "CONST:1.333", "CONST:1.5"},
                  {"P.S2", EventKind::Refract, "CONST:1.5", "CONST:1.6"},
                  {"P.S3", EventKind::Refract, "CONST:1.6", "CONST:1.5"},
                  {"P.S4", EventKind::Refract, "CONST:1.5", "CONST:1.333"}});
  require_events(cs, "explicit",
                 {{"P.S4", EventKind::Refract, "CONST:1.333", "CONST:1.5"},
                  {"P.S3", EventKind::Refract, "CONST:1.5", "CONST:1.6"},
                  {"P.S2", EventKind::Refract, "CONST:1.6", "CONST:1.5"},
                  {"P.S1", EventKind::Refract, "CONST:1.5", "CONST:1.333"}});
}

TEST_CASE("plate of one material with more than 2 surfaces needs an explicit path",
          "[compile][cemented]") {
  // Decided for #27: a prism or cube of one glass toggles inside <-> environment at every face,
  // so the tree order of its faces is no meaningful automatic path (compare the Mangin case).
  const MaterialLibrary lib;
  System s = bare_system();
  s.root.children.push_back(
      {Element{"cube",
               ElementKind::Plate,
               Pose::along_z(20.0),
               "CONST:1.7",
               {plane_surface("C.S1"), plane_surface("C.S2", 1.0), plane_surface("C.S3", 2.0)}}});
  s.paths = {{"auto", true, {}}, {"again", true, {}}};
  const CompileError e = compile_error(s);
  REQUIRE(e.diagnostics().size() == 1);  // reported once per element
  REQUIRE(has_error_at(e, "/root/children/0/surfaces"));
  REQUIRE_THAT(e.what(),
               ContainsSubstring("plate with more than 2 surfaces and uniform material needs an "
                                 "explicit path"));

  SECTION("also with a list of equal entries") {
    std::get<Element>(s.root.children[0].value).material.reset();
    std::get<Element>(s.root.children[0].value).segment_materials = {"CONST:1.7", "CONST:1.7"};
    REQUIRE(has_error_at(compile_error(s), "/root/children/0/surfaces"));
  }
  SECTION("a plate of one material with 2 surfaces stays on the automatic path") {
    std::get<Element>(s.root.children[0].value).surfaces.pop_back();
    const CompiledSystem cs = compile(s, lib);
    require_events(cs, "auto",
                   {{"C.S1", EventKind::Refract, "AIR", "CONST:1.7"},
                    {"C.S2", EventKind::Refract, "CONST:1.7", "AIR"}});
  }
}

TEST_CASE("achromat reference system goes air -> crown -> flint -> air", "[compile][cemented]") {
  // Catalogs follow with #24; the file's SCHOTT glasses are replaced by constants.
  System s = load("m2/achromat.rtt.json");
  auto& lens = std::get<Element>(s.root.children[1].value);
  REQUIRE(lens.segment_materials == std::vector<std::string>{"SCHOTT:N-BK7", "SCHOTT:F2"});
  lens.segment_materials = {"CONST:1.5168", "CONST:1.62"};
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(s, lib);
  REQUIRE(cs.media().size() == 3);
  require_events(cs, "main",
                 {{"STO", EventKind::Transmit, "AIR", "AIR"},
                  {"L1.S1", EventKind::Refract, "AIR", "CONST:1.5168"},
                  {"L1.S2", EventKind::Refract, "CONST:1.5168", "CONST:1.62"},
                  {"L1.S3", EventKind::Refract, "CONST:1.62", "AIR"},
                  {"IMG", EventKind::Transmit, "AIR", "AIR"}});
}

TEST_CASE("unknown segment materials point at the list entry", "[compile][cemented]") {
  SECTION("list: one error per unknown entry") {
    System s = bare_system();
    s.root.children.push_back({segmented_lens(4, {"CONST:1.5", "NOPE:A", "NOPE:B"})});
    s.paths = {{"auto", true, {}}};
    const CompileError e = compile_error(s);
    REQUIRE(e.diagnostics().size() == 2);
    REQUIRE(has_error_at(e, "/root/children/0/material/1"));
    REQUIRE(has_error_at(e, "/root/children/0/material/2"));
    REQUIRE_THAT(e.what(), ContainsSubstring("NOPE:A"));
    REQUIRE_THAT(e.what(), ContainsSubstring("NOPE:B"));
  }
  SECTION("shorthand: one error at the material for all segments") {
    System s = bare_system();
    Element lens = segmented_lens(3, {});
    lens.material = "NOPE:A";
    s.root.children.push_back({lens});
    s.paths = {{"auto", true, {}}};
    const CompileError e = compile_error(s);
    REQUIRE(e.diagnostics().size() == 1);
    REQUIRE(has_error_at(e, "/root/children/0/material"));
  }
}

TEST_CASE("entering element B while inside element A goes from the medium of A to that of B",
          "[compile][cemented]") {
  // Decided for #27: elements do not nest. A Refract at B while the ray is in A leaves A and
  // enters B (medium of A before, of B after); leaving B goes to the environment, not back to
  // A. Cemented groups are one element with one material per segment (ADR 0017).
  System s = bare_system();
  s.root.children.push_back({Element{"A",
                                     ElementKind::Plate,
                                     Pose::along_z(10.0),
                                     "CONST:1.5",
                                     {plane_surface("A.S1"), plane_surface("A.S2", 2.0)}}});
  s.root.children.push_back({segmented_lens(3, {"CONST:1.6", "CONST:1.7"})});
  s.paths = {
      {"A then B", false, {refract("A.S1"), refract("L.S3"), refract("L.S2"), refract("L.S1")}},
      // Exit of A and entry of B directly one after the other: through the environment.
      {"touching", false, {refract("A.S1"), refract("A.S2"), refract("L.S1")}}};
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(s, lib);
  require_events(cs, "A then B",
                 {{"A.S1", EventKind::Refract, "AIR", "CONST:1.5"},
                  {"L.S3", EventKind::Refract, "CONST:1.5", "CONST:1.7"},
                  {"L.S2", EventKind::Refract, "CONST:1.7", "CONST:1.6"},
                  {"L.S1", EventKind::Refract, "CONST:1.6", "AIR"}});
  require_events(cs, "touching",
                 {{"A.S1", EventKind::Refract, "AIR", "CONST:1.5"},
                  {"A.S2", EventKind::Refract, "CONST:1.5", "AIR"},
                  {"L.S1", EventKind::Refract, "AIR", "CONST:1.6"}});
}

TEST_CASE("Refract at a surface that does not bound the current segment", "[compile][cemented]") {
  System s = bare_system();
  s.root.children.push_back({segmented_lens(4, {"CONST:1.5", "CONST:1.6", "CONST:1.7"})});

  SECTION("from inside the element: the ray leaves it into the environment") {
    // Segment rule: through a surface, only the segments on its two sides are reached.
    s.paths = {{"skip",
                false,
                {refract("L.S1"),
                 {SurfaceId("L.S2"), EventKind::Reflect, 0},
                 refract("L.S3"),
                 refract("L.S1")}}};
    const MaterialLibrary lib;
    const CompiledSystem cs = compile(s, lib);
    require_events(cs, "skip",
                   {{"L.S1", EventKind::Refract, "AIR", "CONST:1.5"},
                    {"L.S2", EventKind::Reflect, "CONST:1.5", "CONST:1.5"},
                    {"L.S3", EventKind::Refract, "CONST:1.5", "AIR"},
                    {"L.S1", EventKind::Refract, "AIR", "CONST:1.5"}});
  }
  SECTION("from outside through an inner surface between different glasses: CompileError") {
    s.paths = {{"inner",
                false,
                {refract("L.S1"), refract("L.S1"), {SurfaceId("L.S3"), EventKind::Refract, 0}}}};
    const CompileError e = compile_error(s);
    REQUIRE(e.diagnostics().size() == 1);
    REQUIRE(has_error_at(e, "/paths/0/events/2"));
    REQUIRE_THAT(e.what(), ContainsSubstring("inner surface"));
    REQUIRE_THAT(e.what(), ContainsSubstring("ambiguous"));
  }
}

TEST_CASE("entering a cemented doublet through the cemented surface is ambiguous",
          "[compile][cemented]") {
  System s = bare_system();
  s.root.children.push_back({segmented_lens(3, {"CONST:1.5", "CONST:1.7"})});
  s.paths = {{"auto", true, {}}, {"cemented surface", false, {refract("L.S2")}}};
  const CompileError e = compile_error(s);
  REQUIRE(e.diagnostics().size() == 1);
  REQUIRE(has_error_at(e, "/paths/1/events/0"));
  REQUIRE_THAT(e.what(), ContainsSubstring("L.S2"));
  REQUIRE_THAT(e.what(), ContainsSubstring("CONST:1.5"));
  REQUIRE_THAT(e.what(), ContainsSubstring("CONST:1.7"));
}

TEST_CASE("a prism of one glass can be entered through an inner surface", "[compile][cemented]") {
  // Decided for #27: a plate of one material toggles inside <-> environment at every face, so
  // it can be entered through any of them.
  const MaterialLibrary lib;
  SECTION("plate with the shorthand, 3 surfaces") {
    System s = bare_system();
    s.root.children.push_back(
        {Element{"prism",
                 ElementKind::Plate,
                 Pose::along_z(10.0),
                 "CONST:1.5",
                 {plane_surface("P.S1"), plane_surface("P.S2", 1.0), plane_surface("P.S3", 2.0)}}});
    s.paths = {{"side",
                false,
                {refract("P.S2"), {SurfaceId("P.S1"), EventKind::Reflect, 0}, refract("P.S3")}}};
    const CompiledSystem cs = compile(s, lib);
    require_events(cs, "side",
                   {{"P.S2", EventKind::Refract, "AIR", "CONST:1.5"},
                    {"P.S1", EventKind::Reflect, "CONST:1.5", "CONST:1.5"},
                    {"P.S3", EventKind::Refract, "CONST:1.5", "AIR"}});
  }
  SECTION("plate with a list of equal entries, 4 surfaces") {
    System s = bare_system();
    Element prism{"prism",
                  ElementKind::Plate,
                  Pose::along_z(10.0),
                  std::nullopt,
                  {plane_surface("P.S1"), plane_surface("P.S2", 1.0), plane_surface("P.S3", 2.0),
                   plane_surface("P.S4", 3.0)}};
    prism.segment_materials = {"CONST:1.6", "CONST:1.6", "CONST:1.6"};
    s.root.children.push_back({prism});
    s.paths = {{"side", false, {{SurfaceId("P.S3"), EventKind::Refract, 0}, refract("P.S1")}}};
    const CompiledSystem cs = compile(s, lib);
    require_events(cs, "side",
                   {{"P.S3", EventKind::Refract, "AIR", "CONST:1.6"},
                    {"P.S1", EventKind::Refract, "CONST:1.6", "AIR"}});
  }
}

TEST_CASE("prism of one glass is left through a side face into the environment",
          "[compile][cemented]") {
  // Decided for #27: a Plate of one material (shorthand, or a list of equal entries) keeps the
  // #5 rule, inside <-> environment at every surface, as before #27.
  System s = bare_system();
  s.root.children.push_back(
      {Element{"prism",
               ElementKind::Plate,
               Pose::along_z(10.0),
               "CONST:1.5",
               {plane_surface("P.S1"), plane_surface("P.S2", 1.0), plane_surface("P.S3", 2.0)}}});
  s.paths = {{"through a side face", false, {refract("P.S1"), refract("P.S2")}}};
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(s, lib);
  require_events(cs, "through a side face",
                 {{"P.S1", EventKind::Refract, "AIR", "CONST:1.5"},
                  {"P.S2", EventKind::Refract, "CONST:1.5", "AIR"}});
}

TEST_CASE("explicit path through a cemented doublet follows the segments", "[compile][cemented]") {
  System s = bare_system();
  s.root.children.push_back({segmented_lens(3, {"CONST:1.5", "CONST:1.7"})});
  s.paths = {{"forward", false, {refract("L.S1"), refract("L.S2"), refract("L.S3")}}};
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(s, lib);
  require_events(cs, "forward",
                 {{"L.S1", EventKind::Refract, "AIR", "CONST:1.5"},
                  {"L.S2", EventKind::Refract, "CONST:1.5", "CONST:1.7"},
                  {"L.S3", EventKind::Refract, "CONST:1.7", "AIR"}});
}

TEST_CASE("lens of one glass with 3 surfaces follows the segments on an explicit path too",
          "[compile][cemented]") {
  // Decided for #27: the rule depends on the element, not on the kind of path. An inner surface
  // of a lens is no exit, so the same events give the same media on automatic and explicit
  // paths.
  System s = bare_system();
  Element lens = segmented_lens(3, {});
  lens.material = "CONST:1.5";
  s.root.children.push_back({lens});
  s.paths = {{"auto", true, {}},
             {"explicit", false, {refract("L.S1"), refract("L.S2"), refract("L.S3")}}};
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(s, lib);
  for (const char* path : {"auto", "explicit"}) {
    require_events(cs, path,
                   {{"L.S1", EventKind::Refract, "AIR", "CONST:1.5"},
                    {"L.S2", EventKind::Refract, "CONST:1.5", "CONST:1.5"},
                    {"L.S3", EventKind::Refract, "CONST:1.5", "AIR"}});
  }
}

TEST_CASE("achromat reference system with the SCHOTT test catalogue", "[compile][cemented]") {
  // Recommended in the review of #27: the file's own glasses through the AGF test catalogue
  // (#24). rtt-paraxial may not be used here (ADR 0016: rtt-compile < rtt-paraxial), so the
  // focal length is checked with a hand y-nu trace on the compiled data.
  MaterialLibrary lib;
  lib.add_catalog(std::string(RTT_CATALOG_DIR) + "/schott.agf");
  const CompiledSystem cs = compile(load("m2/achromat.rtt.json"), lib);
  require_events(cs, "main",
                 {{"STO", EventKind::Transmit, "AIR", "AIR"},
                  {"L1.S1", EventKind::Refract, "AIR", "SCHOTT:N-BK7"},
                  {"L1.S2", EventKind::Refract, "SCHOTT:N-BK7", "SCHOTT:F2"},
                  {"L1.S3", EventKind::Refract, "SCHOTT:F2", "AIR"},
                  {"IMG", EventKind::Transmit, "AIR", "AIR"}});

  // Paraxial marginal ray at the reference wavelength (d line), Greivenkamp, OPTI-201/202
  // Sec. 9, p. 9-2 (docs/quellen.md): n'u' = nu - y phi with phi = (n' - n) C, transfer
  // y' = y + u' t'; p. 9-12: phi = -n'u'_k / y_1, f_E = 1/phi, BFD = -y_k / u'_k.
  // Since #25 the indices are absolute and the surrounding AIR has n_air = 1.00027 (Ciddor), so
  // the trace starts in n_air. The catalogue data are relative to air: n_glass = n_rel * n_air,
  // hence every phi and nu is n_air times its value in a relative trace (air = 1) and
  // f' = n_air / phi equals the relative design focal length; BFD is a physical distance and
  // stays as well. Evaluating the catalogue at the air wavelength 0.5876 / n_air (0.16 nm
  // shorter) raises n_rel by 6.7e-6 (N-BK7) and 1.4e-5 (F2); for the achromat this changes f'
  // and BFD by only about 1e-6 relative. The design values are rounded to 1 um (relative trace
  // at 0.5876 um: f' = 101.01333 mm); the tolerance 1e-4 is unchanged.
  const std::size_t wl = cs.reference_wavelength();
  const auto& events = cs.path(*cs.find_path("main")).events;
  const double n_air = cs.media()[cs.environment_medium()].index[wl].real();
  REQUIRE(n_air == rtt::material::ciddor_air_index(cs.wavelengths_um()[wl], 20.0, 1.0));
  double y = 1.0;
  double nu = 0.0;
  double n = n_air;
  double z = 0.0;
  bool first = true;
  for (const auto& event : events) {
    if (event.kind != EventKind::Refract) continue;
    const auto& surface = cs.surfaces()[event.surface];
    const double vertex_z = surface.to_global.translation().z();
    if (!first) y += nu / n * (vertex_z - z);
    first = false;
    z = vertex_z;
    const auto* conic = std::get_if<rtt::geom::Conic<double>>(&surface.shape);
    REQUIRE(conic != nullptr);
    const double n_after = cs.media()[event.medium_after].index[wl].real();
    REQUIRE(cs.media()[event.medium_before].index[wl].real() == n);
    nu -= y * (n_after - n) * conic->curvature();
    n = n_after;
  }
  REQUIRE(n == n_air);
  const double phi = -nu;                    // y_1 = 1
  const double rear_focal_length = n / phi;  // f' = n' / phi
  const double bfd = -y * n / nu;            // -y_k / u'_k with u'_k = nu / n'
  // Design values (relative to air): EFL 101.013 mm; F' at the detector, z = 110.614 mm (L1.S3
  // at z = 14 mm).
  REQUIRE_THAT(rear_focal_length, Catch::Matchers::WithinRel(101.013, 1e-4));
  REQUIRE_THAT(z + bfd, Catch::Matchers::WithinRel(110.614, 1e-4));
}
