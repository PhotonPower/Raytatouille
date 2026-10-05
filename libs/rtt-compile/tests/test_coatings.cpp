#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cmath>
#include <complex>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "rtt/coating/catalog.hpp"
#include "rtt/coating/transfer_matrix.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "test_support.hpp"

using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinAbs;
using rtt::coating::CoatingCatalog;
using rtt::coating::CoatingDesign;
using rtt::coating::CoatingLibrary;
using rtt::coating::LayerSpec;
using rtt::coating::PhysicalThickness;
using rtt::coating::QuarterWaves;
using rtt::compile::CompiledSystem;
using rtt::compile::CompileError;
using rtt::material::MaterialLibrary;
using rtt::model::CoatingRef;
using rtt::model::Element;
using rtt::model::ElementKind;
using rtt::model::System;
using rtt::model::test::element;
using rtt::model::test::has_error_at;

namespace {

const std::string kDemo = std::string(RTT_CATALOG_DIR) + "/coatings/demo.json";

/// Singlet (stop, lens L1 = /root/children/1 with surfaces L1.S1, L1.S2, detector) in glass
/// n = 1.52 with `coating` on both lens surfaces; wavelengths 0.45, 0.55 (reference), 0.65 um.
System coated_singlet(const std::string& coating) {
  System s = rtt::model::test::make_singlet();
  s.wavelengths = {{0.45, 1.0, false}, {0.55, 1.0, true}, {0.65, 1.0, false}};
  element(s, 1).material = "CONST:1.52";
  for (auto& surface : element(s, 1).surfaces) surface.interaction = CoatingRef{coating};
  return s;
}

/// The demo catalogue (CoatingLibrary holds a mutex and cannot be moved).
std::unique_ptr<CoatingLibrary> demo() {
  auto library = std::make_unique<CoatingLibrary>();
  library->add_catalog(kDemo);
  return library;
}

CompileError compile_error(const System& s,
                           const MaterialLibrary& materials,
                           const CoatingLibrary& coatings) {
  try {
    (void)rtt::compile::compile(s, materials, coatings);
  } catch (const CompileError& e) {
    return e;
  }
  FAIL("compile did not throw CompileError");
  return CompileError({});
}

std::uint32_t medium_index(const CompiledSystem& cs, const std::string& reference) {
  for (std::size_t i = 0; i < cs.media().size(); ++i) {
    if (cs.media()[i].reference == reference) return static_cast<std::uint32_t>(i);
  }
  FAIL("medium " << reference << " not found");
  return 0;
}

/// Constant index 1.6 on a limited wavelength range.
class RangedMaterial final : public rtt::material::Material {
 public:
  explicit RangedMaterial(rtt::material::WavelengthRange range) : range_(range) {}
  [[nodiscard]] rtt::math::Complex index(double /*wavelength_um*/,
                                         double /*temperature_c*/,
                                         double /*pressure_atm*/) const override {
    return {1.6, 0.0};
  }
  [[nodiscard]] std::optional<rtt::material::WavelengthRange> wavelength_range_um() const override {
    return range_;
  }

 private:
  rtt::material::WavelengthRange range_;
};

}  // namespace

TEST_CASE("a CoatingRef compiles into layers per wavelength with its substrate",
          "[compile][coating]") {
  const MaterialLibrary materials;
  const auto demo_library = demo();
  const CoatingLibrary& coatings = *demo_library;
  const CompiledSystem cs =
      rtt::compile::compile(coated_singlet("DEMO:AR_MGF2"), materials, coatings);
  // One compiled coating, shared by both surfaces.
  REQUIRE(cs.coatings().size() == 1);
  const rtt::compile::CompiledCoating& ar = cs.coatings()[0];
  REQUIRE(ar.reference == "DEMO:AR_MGF2");
  REQUIRE(ar.layers.size() == 3);  // one stack per system wavelength
  for (const auto& stack : ar.layers) {
    REQUIRE(stack.size() == 1);
    REQUIRE(stack[0].index == std::complex<double>(1.38));
    // QWOT at the design wavelength 0.55 um: d = lambda0 / (4 n) (thickness.hpp).
    REQUIRE_THAT(stack[0].thickness_um, WithinAbs(0.55 / (4.0 * 1.38), 1e-16));
  }
  // Substrate: the lens glass, on its front and its back surface (ADR 0019).
  const std::uint32_t glass = medium_index(cs, "CONST:1.52");
  for (const char* id : {"L1.S1", "L1.S2"}) {
    const auto surface = cs.find_surface(rtt::model::SurfaceId(id));
    REQUIRE(surface.has_value());
    const auto& coating = cs.surfaces()[*surface].coating;
    REQUIRE(coating.has_value());
    REQUIRE(coating->coating == 0);
    REQUIRE(coating->substrate_medium == glass);
  }
  REQUIRE_FALSE(cs.surfaces()[*cs.find_surface(rtt::model::SurfaceId("STO"))].coating);
}

TEST_CASE("the compiled AR coating gives the quarter-wave reflectance in AIR",
          "[compile][coating]") {
  // Reference case of #58 with the compiled data: in AIR (n0 = n_air, Ciddor) at the design
  // wavelength r = (n0 n2 - n1^2) / (n0 n2 + n1^2) (derivation in test_transfer_matrix.cpp).
  const MaterialLibrary materials;
  const auto demo_library = demo();
  const CoatingLibrary& coatings = *demo_library;
  const CompiledSystem cs =
      rtt::compile::compile(coated_singlet("DEMO:AR_MGF2"), materials, coatings);
  const std::size_t wl = cs.reference_wavelength();
  REQUIRE(cs.wavelengths_um()[wl] == 0.55);
  const std::complex<double> n0 = cs.media()[cs.environment_medium()].index[wl];
  const std::complex<double> n2 = cs.media()[medium_index(cs, "CONST:1.52")].index[wl];
  const auto a = rtt::coating::stack_amplitudes<double>(n0, cs.coatings()[0].layers[wl], n2, 0.0,
                                                        cs.wavelengths_um()[wl]);
  const double r = (n0.real() * 1.52 - 1.38 * 1.38) / (n0.real() * 1.52 + 1.38 * 1.38);
  REQUIRE_THAT(std::norm(a.rs), WithinAbs(r * r, 1e-10));
}

TEST_CASE("layer materials resolve through the MaterialLibrary", "[compile][coating]") {
  const auto demo_library = demo();
  const CoatingLibrary& coatings = *demo_library;
  const System s = coated_singlet("DEMO:GLASS_SPACER");
  // Without the SCHOTT catalogue the layer material is unknown.
  const CompileError e = compile_error(s, MaterialLibrary{}, coatings);
  REQUIRE(has_error_at(e.diagnostics(), "/root/children/1/surfaces/0/interaction/name"));
  REQUIRE_THAT(e.what(), ContainsSubstring("layer 0 of coating 'DEMO:GLASS_SPACER'"));
  REQUIRE_THAT(e.what(), ContainsSubstring("SCHOTT"));
  REQUIRE(e.diagnostics().size() == 1);  // reported once, not per surface
  // With it, the layer index is the catalogue glass at the environment temperature and pressure
  // (here 30 degC and 0.9 atm, away from the catalogue's 20 degC / 1 atm).
  MaterialLibrary materials;
  materials.add_catalog(std::string(RTT_CATALOG_DIR) + "/schott.agf");
  System warm = s;
  warm.environment.temperature_c = 30.0;
  warm.environment.pressure_atm = 0.9;
  const CompiledSystem cs = rtt::compile::compile(warm, materials, coatings);
  const auto glass = materials.resolve("SCHOTT:N-BK7");
  for (std::size_t w = 0; w < cs.wavelengths_um().size(); ++w) {
    const double wl = cs.wavelengths_um()[w];
    const std::complex<double> n = cs.coatings()[0].layers[w][0].index;
    REQUIRE_THAT(n.real(), WithinAbs(glass->index(wl, 30.0, 0.9).real(), 1e-15));
    REQUIRE(std::abs(n.real() - glass->index(wl, 20.0, 1.0).real()) > 1e-6);
    REQUIRE(cs.coatings()[0].layers[w][0].thickness_um == 0.2);
  }
}

TEST_CASE("QWOT uses the index at the design wavelength and the environment conditions",
          "[compile][coating]") {
  // Dispersive test material n(lambda, T, p) = 1.3 + 0.1 lambda + 0.001 (T - 20) + 0.01 (p - 1)
  // and a design wavelength 0.6 um that is no system wavelength: d = lambda0 / (4 n(lambda0))
  // at the environment (30 degC, 0.9 atm), the same at every system wavelength, while the
  // layer index follows n(lambda) (thickness.hpp, Byrnes Eq. (8) with (2)).
  class Dispersive final : public rtt::material::Material {
   public:
    [[nodiscard]] rtt::math::Complex index(double wavelength_um,
                                           double temperature_c,
                                           double pressure_atm) const override {
      return {
          1.3 + 0.1 * wavelength_um + 0.001 * (temperature_c - 20.0) + 0.01 * (pressure_atm - 1.0),
          0.0};
    }
  };
  MaterialLibrary materials;
  materials.add("TEST:DISP", std::make_shared<const Dispersive>());
  CoatingLibrary coatings;
  coatings.add(CoatingCatalog{
      "T", {CoatingDesign{"QW", "", 0.6, {LayerSpec{"TEST:DISP", QuarterWaves{1.0, 0.6}}}}}});
  System s = coated_singlet("T:QW");
  s.environment.temperature_c = 30.0;
  s.environment.pressure_atm = 0.9;
  const CompiledSystem cs = rtt::compile::compile(s, materials, coatings);
  const auto n = [](double wl) { return 1.3 + 0.1 * wl + 0.001 * 10.0 + 0.01 * -0.1; };
  const double d = 0.6 / (4.0 * n(0.6));
  for (std::size_t w = 0; w < cs.wavelengths_um().size(); ++w) {
    const double wl = cs.wavelengths_um()[w];
    INFO("wavelength " << wl);
    REQUIRE_THAT(cs.coatings()[0].layers[w][0].thickness_um, WithinAbs(d, 1e-15));
    REQUIRE_THAT(cs.coatings()[0].layers[w][0].index.real(), WithinAbs(n(wl), 1e-15));
  }
  // The design index differs from the index at every system wavelength.
  for (const double wl : cs.wavelengths_um()) REQUIRE(std::abs(n(wl) - n(0.6)) > 1e-3);
}

TEST_CASE("QWOT with Re n <= 0 at the design wavelength is a CompileError", "[compile][coating]") {
  class Negative final : public rtt::material::Material {
   public:
    [[nodiscard]] rtt::math::Complex index(double /*wavelength_um*/,
                                           double /*temperature_c*/,
                                           double /*pressure_atm*/) const override {
      return {-0.5, 0.0};
    }
  };
  MaterialLibrary materials;
  materials.add("TEST:NEG", std::make_shared<const Negative>());
  CoatingLibrary coatings;
  coatings.add(CoatingCatalog{
      "T", {CoatingDesign{"QW", "", 0.55, {LayerSpec{"TEST:NEG", QuarterWaves{1.0, 0.55}}}}}});
  const CompileError e = compile_error(coated_singlet("T:QW"), materials, coatings);
  REQUIRE(has_error_at(e.diagnostics(), "/root/children/1/surfaces/0/interaction/name"));
  REQUIRE_THAT(e.what(), ContainsSubstring("QWOT needs a finite index > 0"));
}

TEST_CASE("cemented lens: first and last surface have their own segment as substrate",
          "[compile][coating]") {
  const MaterialLibrary materials;
  const auto demo_library = demo();
  System s = coated_singlet("DEMO:AR_MGF2");
  Element& lens = element(s, 1);
  rtt::model::Surface back = lens.surfaces[1];
  back.id = rtt::model::SurfaceId("L1.S3");
  back.pose = rtt::model::Pose::along_z(6.0);
  lens.surfaces[1].interaction = rtt::model::Fresnel{};  // the cemented surface stays uncoated
  lens.surfaces.push_back(back);
  lens.material.reset();
  lens.segment_materials = {"CONST:1.52", "CONST:1.62"};
  const CompiledSystem cs = rtt::compile::compile(s, materials, *demo_library);
  const auto front = cs.surfaces()[*cs.find_surface(rtt::model::SurfaceId("L1.S1"))].coating;
  const auto rear = cs.surfaces()[*cs.find_surface(rtt::model::SurfaceId("L1.S3"))].coating;
  REQUIRE(front->substrate_medium == medium_index(cs, "CONST:1.52"));
  REQUIRE(rear->substrate_medium == medium_index(cs, "CONST:1.62"));
  REQUIRE_FALSE(cs.surfaces()[*cs.find_surface(rtt::model::SurfaceId("L1.S2"))].coating);
}

TEST_CASE("unknown coatings are CompileErrors at the interaction name", "[compile][coating]") {
  const MaterialLibrary materials;
  SECTION("no catalogue loaded (compile without CoatingLibrary)") {
    try {
      (void)rtt::compile::compile(coated_singlet("DEMO:AR_MGF2"), materials);
      FAIL("no CompileError");
    } catch (const CompileError& e) {
      REQUIRE(has_error_at(e.diagnostics(), "/root/children/1/surfaces/0/interaction/name"));
      REQUIRE_THAT(e.what(), ContainsSubstring("catalog DEMO not loaded"));
    }
  }
  SECTION("unknown name in a loaded catalogue") {
    const CompileError e = compile_error(coated_singlet("DEMO:NOPE"), materials, *demo());
    REQUIRE_THAT(e.what(), ContainsSubstring("not found in catalog DEMO"));
    REQUIRE(e.diagnostics().size() == 1);
  }
}

TEST_CASE("the substrate side must be unambiguous (ADR 0019)", "[compile][coating]") {
  const MaterialLibrary materials;
  const auto demo_library = demo();
  const CoatingLibrary& coatings = *demo_library;
  SECTION("stop surface") {
    System s = coated_singlet("DEMO:AR_MGF2");
    element(s, 0).surfaces[0].interaction = CoatingRef{"DEMO:AR_MGF2"};
    const CompileError e = compile_error(s, materials, coatings);
    REQUIRE(has_error_at(e.diagnostics(), "/root/children/0/surfaces/0/interaction"));
    REQUIRE_THAT(e.what(), ContainsSubstring("only lens, plate and mirror surfaces"));
  }
  SECTION("detector surface") {
    System s = coated_singlet("DEMO:AR_MGF2");
    element(s, 2).surfaces[0].interaction = CoatingRef{"DEMO:AR_MGF2"};
    REQUIRE(has_error_at(compile_error(s, materials, coatings).diagnostics(),
                         "/root/children/2/surfaces/0/interaction"));
  }
  SECTION("inner surface of a cemented lens") {
    System s = coated_singlet("DEMO:AR_MGF2");
    Element& lens = element(s, 1);
    rtt::model::Surface back = lens.surfaces[1];
    back.id = rtt::model::SurfaceId("L1.S3");
    back.pose = rtt::model::Pose::along_z(6.0);
    lens.surfaces.push_back(back);
    lens.material.reset();
    lens.segment_materials = {"CONST:1.52", "CONST:1.62"};
    const CompileError e = compile_error(s, materials, coatings);
    REQUIRE(has_error_at(e.diagnostics(), "/root/children/1/surfaces/1/interaction"));
    REQUIRE_THAT(e.what(), ContainsSubstring("inner surface between two segments"));
    REQUIRE_FALSE(has_error_at(e.diagnostics(), "/root/children/1/surfaces/0/interaction"));
    REQUIRE_FALSE(has_error_at(e.diagnostics(), "/root/children/1/surfaces/2/interaction"));
  }
  SECTION("mirror without substrate material") {
    System s = coated_singlet("DEMO:AR_MGF2");
    for (auto& surface : element(s, 1).surfaces) surface.interaction = rtt::model::Fresnel{};
    rtt::model::Surface face;
    face.id = rtt::model::SurfaceId("M");
    face.interaction = CoatingRef{"DEMO:AR_MGF2"};
    s.root.children.push_back(
        {Element{"M", ElementKind::Mirror, rtt::model::Pose::along_z(50.0), std::nullopt, {face}}});
    s.paths = {{"main",
                false,
                {{rtt::model::SurfaceId("STO"), rtt::model::EventKind::Transmit, 0},
                 {rtt::model::SurfaceId("M"), rtt::model::EventKind::Reflect, 0}}}};
    const CompileError e = compile_error(s, materials, coatings);
    REQUIRE(has_error_at(e.diagnostics(), "/root/children/3/surfaces/0/interaction"));
    REQUIRE_THAT(e.what(), ContainsSubstring("mirror without substrate"));
  }
}

TEST_CASE("plate of one material and mirror with substrate have their inside as substrate",
          "[compile][coating]") {
  const MaterialLibrary materials;
  const auto demo_library = demo();
  const CoatingLibrary& coatings = *demo_library;
  System s = coated_singlet("DEMO:AR_MGF2");
  element(s, 1).kind = ElementKind::Plate;
  element(s, 1).surfaces[0].shape = element(s, 1).surfaces[1].shape;  // plates are plane
  rtt::model::Surface face;
  face.id = rtt::model::SurfaceId("M");
  face.interaction = CoatingRef{"DEMO:V_AR"};
  s.root.children.push_back(
      {Element{"M", ElementKind::Mirror, rtt::model::Pose::along_z(50.0), "CONST:1.7", {face}}});
  s.paths = {{"main",
              false,
              {{rtt::model::SurfaceId("STO"), rtt::model::EventKind::Transmit, 0},
               {rtt::model::SurfaceId("L1.S1"), rtt::model::EventKind::Refract, 0},
               {rtt::model::SurfaceId("L1.S2"), rtt::model::EventKind::Refract, 0},
               {rtt::model::SurfaceId("M"), rtt::model::EventKind::Reflect, 0}}}};
  const CompiledSystem cs = rtt::compile::compile(s, materials, coatings);
  REQUIRE(cs.coatings().size() == 2);
  const auto plate = cs.surfaces()[*cs.find_surface(rtt::model::SurfaceId("L1.S2"))].coating;
  const auto mirror = cs.surfaces()[*cs.find_surface(rtt::model::SurfaceId("M"))].coating;
  REQUIRE(plate->substrate_medium == medium_index(cs, "CONST:1.52"));
  REQUIRE(mirror->substrate_medium == medium_index(cs, "CONST:1.7"));
  // V_AR: physical and QWOT layer, from the ambient side to the substrate.
  const auto& v = cs.coatings()[mirror->coating].layers[1];
  REQUIRE(v.size() == 2);
  REQUIRE_THAT(v[0].thickness_um, WithinAbs(0.55 / (4.0 * 1.38), 1e-16));
  REQUIRE(v[1].thickness_um == 0.0655);
  REQUIRE(v[1].index == std::complex<double>(2.1));
}

TEST_CASE("wavelength ranges of layer materials are checked for coatings on a path",
          "[compile][coating]") {
  MaterialLibrary materials;
  materials.add("TEST:NARROW",
                std::make_shared<const RangedMaterial>(rtt::material::WavelengthRange{0.5, 0.6}));
  CoatingLibrary coatings;
  coatings.add(CoatingCatalog{
      "T",
      {CoatingDesign{"PHYS", "", std::nullopt, {LayerSpec{"TEST:NARROW", PhysicalThickness{0.1}}}},
       CoatingDesign{"QW_OUT", "", 0.7, {LayerSpec{"TEST:NARROW", QuarterWaves{1.0, 0.7}}}}}});
  SECTION("system wavelength 0.45 um outside [0.5, 0.6]") {
    const CompileError e = compile_error(coated_singlet("T:PHYS"), materials, coatings);
    REQUIRE(has_error_at(e.diagnostics(), "/root/children/1/surfaces/0/interaction/name"));
    REQUIRE_THAT(e.what(), ContainsSubstring("0.45"));
    REQUIRE_THAT(e.what(), ContainsSubstring("layer 0 of coating 'T:PHYS'"));
  }
  SECTION("a coating that no path meets is not checked") {
    System s = coated_singlet("T:PHYS");
    s.paths = {{"skip lens",
                false,
                {{rtt::model::SurfaceId("STO"), rtt::model::EventKind::Transmit, 0},
                 {rtt::model::SurfaceId("IMG"), rtt::model::EventKind::Transmit, 0}}}};
    REQUIRE_NOTHROW(rtt::compile::compile(s, materials, coatings));
  }
  SECTION("QWOT design wavelength outside the material range") {
    const CompileError e = compile_error(coated_singlet("T:QW_OUT"), materials, coatings);
    REQUIRE_THAT(e.what(), ContainsSubstring("design wavelength 0.7 um is outside"));
  }
}

TEST_CASE("the coated reference system validates, round-trips and compiles with its catalogue",
          "[compile][coating]") {
  // tests/reference/m3/ar_singlet.rtt.json: validate and the file round trip need no
  // catalogue (only compile() resolves coatings, ADR 0019).
  const System s = rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m3/ar_singlet.rtt.json");
  REQUIRE_FALSE(rtt::model::has_errors(rtt::model::validate(s)));
  const MaterialLibrary materials;
  REQUIRE_THROWS_AS(rtt::compile::compile(s, materials), CompileError);
  const CompiledSystem cs = rtt::compile::compile(s, materials, *demo());
  REQUIRE(cs.coatings().size() == 1);
  REQUIRE(cs.coatings()[0].reference == "DEMO:AR_MGF2");
}
