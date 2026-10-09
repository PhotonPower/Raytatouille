// Every registered diagnostic code (ADR 0022) with one case that produces it at the expected
// JSON pointer, from model::validate or from compile().

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <functional>
#include <limits>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <type_traits>
#include <utility>
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
using rtt::model::CrystalMaterial;
using rtt::model::Diagnostic;
using rtt::model::DiffractionEfficiency;
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

/// Row of the parameter table (ADR 0029).
rtt::model::ParameterRow row(std::string name, rtt::model::ParameterForm form) {
  rtt::model::ParameterRow r;
  r.name = std::move(name);
  r.form = std::move(form);
  return r;
}

/// RMS spot operand on `path`, field 0, target 0 (ADR 0030).
rtt::model::SpotRmsOperand spot_rms(const std::string& path) {
  rtt::model::SpotRmsOperand o;
  o.path = path;
  return o;
}

/// Ray operand y on path "main" at `surface`, field 0, chief ray (ADR 0030).
rtt::model::RayOperand ray_y(const std::string& surface) {
  rtt::model::RayOperand o;
  o.path = "main";
  o.surface = SurfaceId(surface);
  return o;
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

/// All diagnostics of `s`: validate() and, if it found no error, the errors of compile() or,
/// if it succeeds, its warnings (CompiledSystem::diagnostics(), validate's warnings included).
std::vector<Diagnostic> diagnose(const System& s, std::size_t configuration = 0) {
  std::vector<Diagnostic> d = rtt::model::validate(s);
  if (rtt::model::has_errors(d)) return d;
  try {
    return rtt::compile::compile(s, materials(), coatings(), configuration).diagnostics();
  } catch (const CompileError& e) {
    d.insert(d.end(), e.diagnostics().begin(), e.diagnostics().end());
  }
  return d;
}

struct Case {
  const char* code;
  const char* location;
  std::function<void(System&)> mutate;
  std::size_t configuration = 0;  ///< compiled configuration (#165)
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
      {"surface.efficiency_invalid", "/root/children/1/surfaces/0/diffraction_efficiency",
       [](System& s) {
         element(s, 1).surfaces[0].diffraction_efficiency =
             std::vector<DiffractionEfficiency>{{1, 0.5}};
       }},
      {"crystal.kind_not_allowed", "/root/children/2/material",
       [](System& s) {
         element(s, 2).crystal = CrystalMaterial{"CONST:1.66", "CONST:1.49"};
         element(s, 2).optic_axis = std::array<double, 3>{0.0, 0.0, 1.0};
       }},
      {"crystal.material_conflict", "/root/children/1/material",
       [](System& s) {
         element(s, 1).crystal = CrystalMaterial{"CONST:1.66", "CONST:1.49"};
         element(s, 1).optic_axis = std::array<double, 3>{0.0, 0.0, 1.0};
       }},
      {"crystal.optic_axis_missing", "/root/children/1/optic_axis",
       [](System& s) {
         element(s, 1).material.reset();
         element(s, 1).crystal = CrystalMaterial{"CONST:1.66", "CONST:1.49"};
       }},
      {"crystal.optic_axis_invalid", "/root/children/1/optic_axis",
       [](System& s) {
         element(s, 1).material.reset();
         element(s, 1).crystal = CrystalMaterial{"CONST:1.66", "CONST:1.49"};
         element(s, 1).optic_axis = std::array<double, 3>{0.0, 0.0, 0.0};
       }},
      {"crystal.optic_axis_not_allowed", "/root/children/1/optic_axis",
       [](System& s) { element(s, 1).optic_axis = std::array<double, 3>{0.0, 0.0, 1.0}; }},
      // --- compile()
      {"crystal.unsupported", "/paths/0/events",
       [](System& s) {
         element(s, 1).material.reset();
         element(s, 1).crystal = CrystalMaterial{"CONST:1.66", "CONST:1.49"};
         element(s, 1).optic_axis = std::array<double, 3>{0.0, 1.0, 1.0};
       }},
      {"crystal.absorbing", "/root/children/1/material/ordinary",
       [](System& s) {
         element(s, 1).material.reset();
         element(s, 1).crystal = CrystalMaterial{"CONST:1.66,0.01", "CONST:1.49"};
         element(s, 1).optic_axis = std::array<double, 3>{0.0, 1.0, 1.0};
       }},
      {"crystal.interaction_unsupported", "/root/children/1/surfaces/0/interaction",
       [](System& s) {
         element(s, 1).material.reset();
         element(s, 1).crystal = CrystalMaterial{"CONST:1.66", "CONST:1.49"};
         element(s, 1).optic_axis = std::array<double, 3>{0.0, 1.0, 1.0};
         element(s, 1).surfaces[0].interaction = IdealMirror{};
       }},
      {"paths.crystal_mode_required", "/paths/0/events/1",
       [](System& s) {
         element(s, 1).material.reset();
         element(s, 1).crystal = CrystalMaterial{"CONST:1.66", "CONST:1.49"};
         element(s, 1).optic_axis = std::array<double, 3>{0.0, 1.0, 1.0};
         s.paths = {{"x",
                     false,
                     {{SurfaceId("STO"), EventKind::Transmit, 0},
                      {SurfaceId("L1.S1"), EventKind::Refract, 0},
                      {SurfaceId("L1.S2"), EventKind::Refract, 0}}}};
       }},
      {"paths.mode_without_crystal", "/paths/0/events/1",
       [](System& s) {
         s.paths = {{"x",
                     false,
                     {{SurfaceId("STO"), EventKind::Transmit, 0},
                      {SurfaceId("L1.S1"), EventKind::Ordinary, 0},
                      {SurfaceId("L1.S2"), EventKind::Refract, 0}}}};
       }},
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
      // Second producing sites of codes with several (second review of #98, H1).
      {"element.surface_count", "/root/children/1/surfaces",
       [](System& s) {
         element(s, 1).kind = ElementKind::Plate;
         element(s, 1).surfaces.pop_back();
       }},
      {"element.surface_count", "/root/children/1/surfaces",
       [](System& s) {
         element(s, 1).kind = ElementKind::Mirror;
         element(s, 1).surfaces.clear();
       }},
      {"element.surface_count", "/root/children/0/surfaces",
       [](System& s) { element(s, 0).surfaces.push_back(element(s, 0).surfaces[0]); }},
      {"element.material_empty", "/root/children/1/material/0",
       [](System& s) {
         element(s, 1).material.reset();
         element(s, 1).segment_materials = {""};
       }},
      {"interaction.axis_invalid", "/root/children/1/surfaces/0/interaction/fast_axis",
       [](System& s) {
         element(s, 1).surfaces[0].interaction = IdealRetarder{{0.0, 0.0, 0.0}, 0.25};
       }},
      {"node.name_empty", "/root/name", [](System& s) { s.root.name.clear(); }},
      // coating.no_substrate for a lens or plate is not reachable after validate(): both need a
      // material (element.material_missing), and an unresolved one keeps a placeholder medium.
      {"stop.not_on_path", "/paths/1",
       [](System& s) {
         // The stop inside an assembly is found as well (review of #86 B).
         rtt::model::Assembly group;
         group.name = "front";
         group.children.push_back(s.root.children[0]);
         s.root.children[0] = {group};
         s.paths.push_back({"lens only",
                            false,
                            {{SurfaceId("L1.S1"), EventKind::Refract, 0},
                             {SurfaceId("L1.S2"), EventKind::Refract, 0}}});
       }},
      {"stop.not_on_path", "/paths/1",
       [](System& s) {
         s.paths.push_back({"lens only",
                            false,
                            {{SurfaceId("L1.S1"), EventKind::Refract, 0},
                             {SurfaceId("L1.S2"), EventKind::Refract, 0}}});
       }},
      // --- schema 0.4 (#162): relative placement (ADR 0028), parameter table (ADR 0029)
      {"pose.relative_first_surface", "/root/children/1/surfaces/0/pose/reference",
       [](System& s) {
         element(s, 1).surfaces[0].pose.reference = PoseReference::RelativeToSibling;
       }},
      {"pose.no_preceding", "/root/children/0/pose/reference",
       [](System& s) { element(s, 0).pose.reference = PoseReference::RelativeToPreceding; }},
      {"pose.no_sibling", "/root/children/0/pose/reference",
       [](System& s) { element(s, 0).pose.reference = PoseReference::RelativeToSibling; }},
      {"parameters.name_invalid", "/parameters/0/name",
       [](System& s) { s.parameters = {row("2D", 1.0)}; }},
      {"parameters.name_duplicate", "/parameters/1/name",
       [](System& s) { s.parameters = {row("D", 1.0), row("D", 2.0)}; }},
      {"parameters.values_count", "/parameters/0/values",
       [](System& s) { s.parameters = {row("D", std::vector<double>{1.0, 2.0})}; }},
      {"parameters.variable_expression", "/parameters/1/variable",
       [](System& s) {
         s.parameters = {row("D", 1.0), row("E", ParameterExpression{"2 * D"})};
         s.parameters[1].variable = true;
       }},
      // evaluating codes of the table (#164)
      {"parameters.expression_syntax", "/parameters/0/expression",
       [](System& s) { s.parameters = {row("D", ParameterExpression{"2 +"})}; }},
      {"parameters.unknown_name", "/parameters/0/expression",
       [](System& s) { s.parameters = {row("D", ParameterExpression{"2 * E"})}; }},
      {"parameters.forward_reference", "/parameters/0/expression",
       [](System& s) { s.parameters = {row("D", ParameterExpression{"2 * E"}), row("E", 1.0)}; }},
      {"parameters.not_finite", "/parameters/1/expression",
       [](System& s) {
         s.configurations = {{"a"}, {"b"}};
         s.parameters = {row("D", std::vector<double>{1.0, 0.0}),
                         row("E", ParameterExpression{"1 / D"})};
       }},
      {"param.unknown_parameter", "/root/children/2/pose/position/2/param",
       [](System& s) { element(s, 2).pose.position[2] = Param::bound("D"); }},
      {"param.bound_conflict", "/root/children/2/pose/position/2",
       [](System& s) {
         s.parameters = {row("D", 106.363)};
         Param p = Param::bound("D");
         p.variable = true;
         element(s, 2).pose.position[2] = p;
       }},
      {"bounds.invalid", "/parameters/0/min",
       [](System& s) {
         s.parameters = {row("D", 1.0)};
         s.parameters[0].min = 2.0;
         s.parameters[0].max = 0.0;
       }},
      {"bounds.value_outside", "/parameters/0/value",
       [](System& s) {
         s.parameters = {row("D", 1.0)};
         s.parameters[0].max = 0.5;
       }},
      {"configurations.name_invalid", "/configurations/0/name",
       [](System& s) { s.configurations = {{""}}; }},
      {"configurations.name_duplicate", "/configurations/1/name",
       [](System& s) { s.configurations = {{"a"}, {"a"}}; }},
      {"value.not_finite", "/root/children/1/pose/pivot/2",
       [](System& s) { element(s, 1).pose.pivot[2] = kNaN; }},
      // validate, merit function (ADR 0030, #162 part B)
      {"merit.unknown_path", "/optimization/operands/0/path",
       [](System& s) { s.optimization.operands = {spot_rms("other")}; }},
      {"merit.unknown_configuration", "/optimization/operands/0/configuration",
       [](System& s) {
         rtt::model::SpotRmsOperand o = spot_rms("main");
         o.common.configuration = "far";
         s.optimization.operands = {o};
       }},
      {"merit.unknown_parameter", "/optimization/operands/0/parameter",
       [](System& s) { s.optimization.operands = {rtt::model::ParamValueOperand{{}, "D"}}; }},
      {"merit.surface_not_on_path", "/optimization/operands/0/surface",
       [](System& s) { s.optimization.operands = {ray_y("NONE")}; }},
      {"merit.surface_ambiguous", "/optimization/operands/0/surface",
       [](System& s) {
         s.paths = {{"main",
                     false,
                     {{SurfaceId("STO"), EventKind::Transmit, 0},
                      {SurfaceId("L1.S1"), EventKind::Refract, 0},
                      {SurfaceId("L1.S2"), EventKind::Reflect, 0},
                      {SurfaceId("L1.S1"), EventKind::Refract, 0}}}};
         s.optimization.operands = {ray_y("L1.S1")};
       }},
      {"merit.index_out_of_range", "/optimization/operands/0/field",
       [](System& s) {
         rtt::model::SpotRmsOperand o = spot_rms("main");
         o.field = 9;
         s.optimization.operands = {o};
       }},
      {"merit.weight_invalid", "/optimization/operands/0/weight",
       [](System& s) {
         rtt::model::SpotRmsOperand o = spot_rms("main");
         o.common.weight = -1.0;
         s.optimization.operands = {o};
       }},
      {"merit.sampling_invalid", "/optimization/operands/0/rings",
       [](System& s) {
         rtt::model::SpotRmsOperand o = spot_rms("main");
         o.rings = 0;
         s.optimization.operands = {o};
       }},
      {"merit.selection_empty", "/optimization/generators/0/fields",
       [](System& s) {
         rtt::model::SpotGenerator g;
         g.path = "main";
         g.fields = std::vector<std::uint16_t>{};
         s.optimization.generators = {g};
       }},
      {"merit.polychromatic_wavelength", "/optimization/operands/0/wavelength",
       [](System& s) {
         rtt::model::SpotRmsOperand o = spot_rms("main");
         o.polychromatic = true;
         o.wavelength = 0;
         s.optimization.operands = {o};
       }},
      // compile: a configuration index beyond the configurations (#165)
      {"config.unknown", "/configurations", [](System& s) { s.configurations = {{"a"}}; }, 1},
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
    const std::vector<Diagnostic> d = diagnose(s, c.configuration);
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
    if (info.producer != "validate" && info.producer != "compile") continue;  // see rtt-analysis
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
