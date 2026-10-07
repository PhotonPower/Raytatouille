// Every registered diagnostic code (ADR 0022) with one case that produces it at the expected
// JSON pointer, from model::validate or from compile().

#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <functional>
#include <limits>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <type_traits>
#include <vector>

#include "rtt/coating/catalog.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/diagnostics/codes.hpp"
#include "rtt/material/material.hpp"
#include "test_support.hpp"

using rtt::coating::CoatingCatalog;
using rtt::coating::CoatingDesign;
using rtt::coating::CoatingLibrary;
using rtt::coating::LayerSpec;
using rtt::coating::PhysicalThickness;
using rtt::coating::QuarterWaves;
using rtt::compile::CompileError;
using rtt::material::MaterialLibrary;
using rtt::model::Diagnostic;
using rtt::model::Element;
using rtt::model::ElementKind;
using rtt::model::EventKind;
using rtt::model::Param;
using rtt::model::Pose;
using rtt::model::Severity;
using rtt::model::SurfaceId;
using rtt::model::System;
using rtt::model::test::element;

// The model's Severity is the registry's (rtt-diagnostics, layer Basis).
static_assert(std::is_same_v<rtt::model::Severity, rtt::diagnostics::Severity>);

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

/// Constant index on a limited wavelength range, in um.
class RangedMaterial final : public rtt::material::Material {
 public:
  RangedMaterial(double n, rtt::material::WavelengthRange range) : n_(n), range_(range) {}
  [[nodiscard]] rtt::math::Complex index(double /*wavelength_um*/,
                                         double /*temperature_c*/,
                                         double /*pressure_atm*/) const override {
    return {n_, 0.0};
  }
  [[nodiscard]] std::optional<rtt::material::WavelengthRange> wavelength_range_um() const override {
    return range_;
  }

 private:
  double n_;
  rtt::material::WavelengthRange range_;
};

/// Materials of the compile cases: TEST:NARROW (n = 1.6 on [0.6, 0.7] um, outside the system
/// wavelength 0.5876 um), TEST:WIDE (n = 1.6 on [0.5, 0.7] um) and TEST:NEG (n = -0.5).
const MaterialLibrary& materials() {
  static const auto library = [] {
    auto m = std::make_unique<MaterialLibrary>();  // holds a mutex, cannot be moved
    m->add("TEST:NARROW",
           std::make_shared<const RangedMaterial>(1.6, rtt::material::WavelengthRange{0.6, 0.7}));
    m->add("TEST:WIDE",
           std::make_shared<const RangedMaterial>(1.6, rtt::material::WavelengthRange{0.5, 0.7}));
    m->add("TEST:NEG",
           std::make_shared<const RangedMaterial>(-0.5, rtt::material::WavelengthRange{0.3, 2.0}));
    return m;
  }();
  return *library;
}

/// Coatings T:OK, T:NARROW (layer outside the system wavelength), T:BAD_LAYER (unknown layer
/// material), T:QW_OUT (design wavelength outside the layer range), T:QW_NEG (QWOT with n < 0).
const CoatingLibrary& coatings() {
  static const auto library = [] {
    auto c = std::make_unique<CoatingLibrary>();
    c->add(CoatingCatalog{
        "T",
        {CoatingDesign{"OK", "", std::nullopt, {LayerSpec{"CONST:1.38", PhysicalThickness{0.1}}}},
         CoatingDesign{
             "NARROW", "", std::nullopt, {LayerSpec{"TEST:NARROW", PhysicalThickness{0.1}}}},
         CoatingDesign{
             "BAD_LAYER", "", std::nullopt, {LayerSpec{"NOPE:X", PhysicalThickness{0.1}}}},
         CoatingDesign{"QW_OUT", "", 0.9, {LayerSpec{"TEST:WIDE", QuarterWaves{1.0, 0.9}}}},
         CoatingDesign{"QW_NEG", "", 0.55, {LayerSpec{"TEST:NEG", QuarterWaves{1.0, 0.55}}}}}});
    return c;
  }();
  return *library;
}

/// Reference singlet (stop /root/children/0, lens L1 /root/children/1 with L1.S1 and L1.S2,
/// detector /root/children/2) in CONST:1.5168, so that it compiles without catalogue.
System singlet() {
  System s = rtt::model::test::make_singlet();
  element(s, 1).material = "CONST:1.5168";
  return s;
}

/// Lens L1 as a cemented doublet with three plane surfaces in CONST:1.5 and CONST:1.7.
void make_doublet(System& s) {
  Element& lens = element(s, 1);
  for (auto& surface : lens.surfaces) surface.shape.base = rtt::model::Plane{};
  rtt::model::Surface back = lens.surfaces[1];
  back.id = SurfaceId("L1.S3");
  back.pose = Pose::along_z(6.0);
  lens.surfaces.push_back(back);
  lens.material.reset();
  lens.segment_materials = {"CONST:1.5", "CONST:1.7"};
}

/// All diagnostics of `s`: validate() and, if it found no error, the errors of compile().
std::vector<Diagnostic> diagnose(const System& s) {
  std::vector<Diagnostic> d = rtt::model::validate(s);
  if (rtt::model::has_errors(d)) return d;
  try {
    (void)rtt::compile::compile(s, materials(), coatings());
  } catch (const CompileError& e) {
    d.insert(d.end(), e.diagnostics().begin(), e.diagnostics().end());
  }
  return d;
}

struct Case {
  const char* code;
  const char* location;
  std::function<void(System&)> mutate;
};

std::vector<Case> cases() {
  using namespace rtt::model;
  return {
      // --- validate()
      {"wavelengths.empty", "/wavelengths", [](System& s) { s.wavelengths.clear(); }},
      {"wavelengths.value_invalid", "/wavelengths/0/um",
       [](System& s) { s.wavelengths[0].um = -1.0; }},
      {"wavelengths.weight_invalid", "/wavelengths/0/weight",
       [](System& s) { s.wavelengths[0].weight = -1.0; }},
      {"wavelengths.reference_count", "/wavelengths",
       [](System& s) { s.wavelengths.push_back({0.4861, 1.0, true}); }},
      {"aperture.value_invalid", "/aperture/value",
       [](System& s) { s.aperture.value = Param(0.0); }},
      {"aperture.na_not_physical", "/aperture/value",
       [](System& s) {
         s.aperture = {SystemApertureType::ObjectSpaceNA, Param(1.2)};
         s.object = {false, Param(100.0)};
       }},
      {"aperture.stop_missing", "/aperture/type",
       [](System& s) {
         s.aperture.type = SystemApertureType::StopSize;
         element(s, 0).kind = ElementKind::ThinElement;
       }},
      {"fields.empty", "/fields/points", [](System& s) { s.fields.points.clear(); }},
      {"fields.coordinate_invalid", "/fields/points/0",
       [](System& s) { s.fields.points[0].x = kNaN; }},
      {"fields.weight_invalid", "/fields/points/0/weight",
       [](System& s) { s.fields.points[0].weight = -1.0; }},
      {"object.distance_invalid", "/object/distance",
       [](System& s) { s.object = {false, Param(-5.0)}; }},
      {"environment.temperature_invalid", "/environment/temperature_c",
       [](System& s) { s.environment.temperature_c = -300.0; }},
      {"environment.pressure_invalid", "/environment/pressure_atm",
       [](System& s) { s.environment.pressure_atm = -1.0; }},
      {"environment.medium_empty", "/environment/medium",
       [](System& s) { s.environment.medium.clear(); }},
      {"node.name_empty", "/root/children/1/name", [](System& s) { element(s, 1).name.clear(); }},
      {"node.name_duplicate", "/root/children/2/name",
       [](System& s) { element(s, 2).name = "L1"; }},
      {"element.surface_count", "/root/children/1/surfaces",
       [](System& s) { element(s, 1).surfaces.pop_back(); }},
      {"element.plate_surface_not_plane", "/root/children/1/surfaces/0/shape",
       [](System& s) { element(s, 1).kind = ElementKind::Plate; }},
      {"element.material_missing", "/root/children/1/material",
       [](System& s) { element(s, 1).material.reset(); }},
      {"element.material_not_allowed", "/root/children/2/material",
       [](System& s) { element(s, 2).material = "AIR"; }},
      {"element.material_list_not_allowed", "/root/children/1/material",
       [](System& s) {
         element(s, 1).kind = ElementKind::Mirror;
         element(s, 1).material.reset();
         element(s, 1).segment_materials = {"CONST:1.5"};
       }},
      {"element.material_both", "/root/children/1/material",
       [](System& s) { element(s, 1).segment_materials = {"CONST:1.5"}; }},
      {"element.material_empty", "/root/children/1/material",
       [](System& s) { element(s, 1).material = ""; }},
      {"element.segment_count", "/root/children/1/material",
       [](System& s) {
         element(s, 1).material.reset();
         element(s, 1).segment_materials = {"CONST:1.5", "CONST:1.6"};
       }},
      {"stop.aperture_missing", "/root/children/0/surfaces/0/aperture",
       [](System& s) { element(s, 0).surfaces[0].aperture.reset(); }},
      {"stop.multiple", "/root/children/2",
       [](System& s) {
         element(s, 2).kind = ElementKind::Stop;
         element(s, 2).surfaces[0].aperture = CircularAperture{10.0, 0.0};
       }},
      {"surface.id_empty", "/root/children/1/surfaces/0/id",
       [](System& s) { element(s, 1).surfaces[0].id = SurfaceId(""); }},
      {"surface.id_duplicate", "/root/children/1/surfaces/1/id",
       [](System& s) { element(s, 1).surfaces[1].id = SurfaceId("L1.S1"); }},
      {"shape.radius_invalid", "/root/children/1/surfaces/0/shape/base/radius",
       [](System& s) { element(s, 1).surfaces[0].shape.base = Conic{Param(0.0), Param(0.0)}; }},
      {"shape.asphere_without_coefficients", "/root/children/1/surfaces/0/shape/base/coefficients",
       [](System& s) {
         element(s, 1).surfaces[0].shape.base = EvenAsphere{Param(51.68), Param(0.0), {}};
       }},
      {"shape.zernike_radius_invalid",
       "/root/children/1/surfaces/0/shape/terms/0/normalization_radius",
       [](System& s) {
         element(s, 1).surfaces[0].shape.terms.push_back(ZernikeSag{Param(0.0), {Param(1e-3)}});
       }},
      {"shape.zernike_without_coefficients",
       "/root/children/1/surfaces/0/shape/terms/0/coefficients",
       [](System& s) {
         element(s, 1).surfaces[0].shape.terms.push_back(ZernikeSag{Param(10.0), {}});
       }},
      {"surface_aperture.radius_invalid", "/root/children/1/surfaces/0/aperture/radius",
       [](System& s) { element(s, 1).surfaces[0].aperture = CircularAperture{0.0, 0.0}; }},
      {"surface_aperture.inner_radius_invalid", "/root/children/1/surfaces/0/aperture/inner_radius",
       [](System& s) { element(s, 1).surfaces[0].aperture = CircularAperture{12.7, 13.0}; }},
      {"surface_aperture.half_width_invalid", "/root/children/1/surfaces/0/aperture",
       [](System& s) { element(s, 1).surfaces[0].aperture = RectangularAperture{0.0, 1.0}; }},
      {"surface_aperture.semi_axis_invalid", "/root/children/1/surfaces/0/aperture",
       [](System& s) { element(s, 1).surfaces[0].aperture = EllipticalAperture{1.0, 0.0}; }},
      {"phase.lines_per_mm_invalid", "/root/children/1/surfaces/0/phases/0/lines_per_mm",
       [](System& s) {
         element(s, 1).surfaces[0].phases.push_back(LinearGrating{Param(0.0), 0.0});
       }},
      {"phase.radius_invalid", "/root/children/1/surfaces/0/phases/0/normalization_radius",
       [](System& s) { element(s, 1).surfaces[0].phases.push_back(RadialPhase{Param(0.0), {}}); }},
      {"interaction.reflectance_invalid", "/root/children/1/surfaces/0/interaction",
       [](System& s) { element(s, 1).surfaces[0].interaction = IdealBeamSplitter{1.5, 0.5}; }},
      {"interaction.coating_name_empty", "/root/children/1/surfaces/0/interaction/name",
       [](System& s) { element(s, 1).surfaces[0].interaction = CoatingRef{""}; }},
      {"interaction.axis_invalid", "/root/children/1/surfaces/0/interaction/transmission_axis",
       [](System& s) {
         element(s, 1).surfaces[0].interaction = IdealPolarizer{{0.0, 0.0, 0.0}, 0.0};
       }},
      {"interaction.extinction_ratio_invalid",
       "/root/children/1/surfaces/0/interaction/extinction_ratio",
       [](System& s) {
         element(s, 1).surfaces[0].interaction = IdealPolarizer{{1.0, 0.0, 0.0}, 2.0};
       }},
      {"interaction.retardance_invalid", "/root/children/1/surfaces/0/interaction/retardance_waves",
       [](System& s) {
         element(s, 1).surfaces[0].interaction = IdealRetarder{{1.0, 0.0, 0.0}, kNaN};
       }},
      {"paths.empty", "/paths", [](System& s) { s.paths.clear(); }},
      {"paths.name_empty", "/paths/0/name", [](System& s) { s.paths[0].name.clear(); }},
      {"paths.name_duplicate", "/paths/1/name", [](System& s) { s.paths.push_back(s.paths[0]); }},
      {"paths.events_empty", "/paths/0/events", [](System& s) { s.paths = {{"x", false, {}}}; }},
      {"paths.unknown_surface", "/paths/0/events/0/surface",
       [](System& s) { s.paths = {{"x", false, {{SurfaceId("NOPE"), EventKind::Transmit, 0}}}}; }},
      {"paths.order_not_allowed", "/paths/0/events/0/order",
       [](System& s) { s.paths = {{"x", false, {{SurfaceId("STO"), EventKind::Transmit, 1}}}}; }},
      // --- compile()
      {"wavelengths.too_many", "/wavelengths",
       [](System& s) {
         s.wavelengths.assign(65536, {0.5876, 1.0, false});
         s.wavelengths[0].reference = true;
       }},
      {"material.unknown", "/root/children/1/material",
       [](System& s) { element(s, 1).material = "NOPE:X"; }},
      {"material.wavelength_out_of_range", "/root/children/1/material",
       [](System& s) { element(s, 1).material = "TEST:NARROW"; }},
      {"shape.zernike_unsupported", "/root/children/1/surfaces/0/shape/terms/0",
       [](System& s) {
         element(s, 1).surfaces[0].shape.terms.push_back(ZernikeSag{Param(10.0), {Param(1e-3)}});
       }},
      {"coating.unknown", "/root/children/1/surfaces/0/interaction/name",
       [](System& s) { element(s, 1).surfaces[0].interaction = CoatingRef{"T:NOPE"}; }},
      {"coating.no_substrate", "/root/children/1/surfaces/0/interaction",
       [](System& s) {
         element(s, 1).kind = ElementKind::Mirror;
         element(s, 1).material.reset();
         element(s, 1).surfaces[0].interaction = CoatingRef{"T:OK"};
       }},
      {"coating.not_allowed", "/root/children/0/surfaces/0/interaction",
       [](System& s) { element(s, 0).surfaces[0].interaction = CoatingRef{"T:OK"}; }},
      {"coating.substrate_ambiguous", "/root/children/1/surfaces/1/interaction",
       [](System& s) {
         make_doublet(s);
         element(s, 1).surfaces[1].interaction = CoatingRef{"T:OK"};
       }},
      {"coating.layer_material_unknown", "/root/children/1/surfaces/0/interaction/name",
       [](System& s) { element(s, 1).surfaces[0].interaction = CoatingRef{"T:BAD_LAYER"}; }},
      {"coating.design_wavelength_out_of_range", "/root/children/1/surfaces/0/interaction/name",
       [](System& s) { element(s, 1).surfaces[0].interaction = CoatingRef{"T:QW_OUT"}; }},
      {"coating.thickness_invalid", "/root/children/1/surfaces/0/interaction/name",
       [](System& s) { element(s, 1).surfaces[0].interaction = CoatingRef{"T:QW_NEG"}; }},
      {"coating.wavelength_out_of_range", "/root/children/1/surfaces/0/interaction/name",
       [](System& s) { element(s, 1).surfaces[0].interaction = CoatingRef{"T:NARROW"}; }},
      {"paths.mangin_mirror_automatic", "/root/children/1/surfaces",
       [](System& s) { element(s, 1).kind = ElementKind::Mirror; }},
      {"paths.uniform_plate_automatic", "/root/children/1/surfaces",
       [](System& s) {
         make_doublet(s);
         element(s, 1).kind = ElementKind::Plate;
         element(s, 1).segment_materials.clear();
         element(s, 1).material = "CONST:1.5";
       }},
      {"paths.inner_surface_ambiguous", "/paths/0/events/0",
       [](System& s) {
         make_doublet(s);
         s.paths = {{"x", false, {{SurfaceId("L1.S2"), EventKind::Refract, 0}}}};
       }},
  };
}

}  // namespace

TEST_CASE("the unchanged singlet has no diagnostics", "[diagnostics]") {
  REQUIRE(diagnose(singlet()).empty());
}

TEST_CASE("every registered code has a case that produces it at its location", "[diagnostics]") {
  std::set<std::string> covered;
  for (const Case& c : cases()) {
    INFO(c.code);
    const rtt::diagnostics::CodeInfo* info = rtt::diagnostics::find_code(c.code);
    REQUIRE(info != nullptr);
    System s = singlet();
    c.mutate(s);
    const std::vector<Diagnostic> d = diagnose(s);
    bool found = false;
    for (const Diagnostic& x : d) {
      INFO(rtt::model::to_string(x));
      // Every diagnostic has a registered code and the severity of its registry entry.
      const rtt::diagnostics::CodeInfo* registered = rtt::diagnostics::find_code(x.code);
      REQUIRE(registered != nullptr);
      REQUIRE(x.severity == registered->severity);
      found = found || (x.code == c.code && x.location == c.location);
    }
    REQUIRE(found);
    covered.insert(c.code);
  }
  for (const rtt::diagnostics::CodeInfo& info : rtt::diagnostics::kCodes) {
    INFO("no case for " << info.code);
    REQUIRE(covered.contains(std::string(info.code)));
  }
}

TEST_CASE("to_string puts the code in brackets", "[diagnostics]") {
  const Diagnostic error{Severity::Error, "/root/children/1/material", "material 'X' not found",
                         "material.unknown"};
  REQUIRE(rtt::model::to_string(error) ==
          "error [material.unknown] /root/children/1/material: material 'X' not found");
  const Diagnostic warning{Severity::Warning, "/aperture/value",
                           "object-space NA >= 1 in air is not physical",
                           "aperture.na_not_physical"};
  REQUIRE(rtt::model::to_string(warning) ==
          "warning [aperture.na_not_physical] /aperture/value: object-space NA >= 1 in air is "
          "not physical");
}

TEST_CASE("a Diagnostic without code prints as before", "[diagnostics]") {
  const Diagnostic d{Severity::Error, "/paths", "at least one path is required", ""};
  REQUIRE(d.code.empty());
  REQUIRE(rtt::model::to_string(d) == "error /paths: at least one path is required");
}

TEST_CASE("CompileError lists the codes in its message", "[diagnostics]") {
  System s = singlet();
  element(s, 1).material = "NOPE:X";
  try {
    (void)rtt::compile::compile(s, materials(), coatings());
    FAIL("no CompileError");
  } catch (const CompileError& e) {
    REQUIRE(e.diagnostics().size() == 1);
    REQUIRE(e.diagnostics()[0].code == "material.unknown");
    REQUIRE(std::string(e.what()).find("error [material.unknown] /root/children/1/material: ") !=
            std::string::npos);
  }
}
