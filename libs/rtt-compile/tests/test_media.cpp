#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cmath>
#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "test_support.hpp"

using Catch::Matchers::ContainsSubstring;
using rtt::compile::CompiledSystem;
using rtt::compile::CompileError;
using rtt::material::Material;
using rtt::material::MaterialLibrary;
using rtt::material::WavelengthRange;
using rtt::math::Complex;
using rtt::model::Element;
using rtt::model::System;

namespace {

/// Returns n + i*kappa = pressure + i*temperature, to check what compile() passes on.
class EchoMaterial final : public Material {
 public:
  [[nodiscard]] Complex index(double /*wavelength_um*/,
                              double temperature_c,
                              double pressure_atm) const override {
    return {pressure_atm, temperature_c};
  }
};

/// Constant index 1.5 on a limited wavelength range.
class RangedMaterial final : public Material {
 public:
  explicit RangedMaterial(WavelengthRange range) : range_(range) {}
  [[nodiscard]] Complex index(double /*wavelength_um*/,
                              double /*temperature_c*/,
                              double /*pressure_atm*/) const override {
    return {1.5, 0.0};
  }
  [[nodiscard]] std::optional<WavelengthRange> wavelength_range_um() const override {
    return range_;
  }

 private:
  WavelengthRange range_;
};

/// Singlet (stop, lens L1 at /root/children/1, detector) with the lens made of `material`
/// and the wavelengths 0.4861, 0.5876 (reference), 0.6563 um.
System singlet(const std::string& material) {
  System s = rtt::model::test::make_singlet();
  s.wavelengths = {{0.4861, 1.0, false}, {0.5876, 1.0, true}, {0.6563, 1.0, false}};
  rtt::model::test::element(s, 1).material = material;
  return s;
}

std::uint32_t medium_index(const CompiledSystem& cs, const std::string& reference) {
  for (std::size_t i = 0; i < cs.media().size(); ++i) {
    if (cs.media()[i].reference == reference) return static_cast<std::uint32_t>(i);
  }
  FAIL("medium " << reference << " not found");
  return 0;
}

CompileError compile_error(const System& s, const MaterialLibrary& lib) {
  try {
    (void)rtt::compile::compile(s, lib);
  } catch (const CompileError& e) {
    return e;
  }
  FAIL("compile did not throw CompileError");
  return CompileError({});
}

}  // namespace

TEST_CASE("SCHOTT:N-BK7 in a system file resolves through an AGF test catalogue",
          "[compile][media]") {
  // Acceptance criterion of #24: tests/reference/m0/singlet.rtt.json uses SCHOTT:N-BK7.
  MaterialLibrary lib;
  lib.add_catalog(std::string(RTT_CATALOG_DIR) + "/schott.agf");
  const System s = rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m0/singlet.rtt.json");
  const CompiledSystem cs = rtt::compile::compile(s, lib);
  const auto& glass = cs.media()[medium_index(cs, "SCHOTT:N-BK7")];
  // Same value as the catalogue material itself; n_d = 1.5168 in the SCHOTT data sheet.
  const auto expected = lib.resolve("SCHOTT:N-BK7")
                            ->index(cs.wavelengths_um()[cs.reference_wavelength()],
                                    s.environment.temperature_c, s.environment.pressure_atm);
  REQUIRE(glass.index[cs.reference_wavelength()] == expected);
  REQUIRE(std::abs(expected.real() - 1.5168) < 5e-6);
}

TEST_CASE("compile passes temperature and pressure of the environment to the materials",
          "[compile][media]") {
  MaterialLibrary lib;
  lib.add("TEST:ECHO", std::make_shared<const EchoMaterial>());
  System s = singlet("TEST:ECHO");
  s.environment.temperature_c = 25.0;
  s.environment.pressure_atm = 0.8;
  s.environment.medium = "TEST:ECHO";
  const CompiledSystem cs = rtt::compile::compile(s, lib);
  const auto& index = cs.media()[medium_index(cs, "TEST:ECHO")].index;
  REQUIRE(index == std::vector<Complex>(3, Complex(0.8, 25.0)));
  REQUIRE(cs.environment_medium() == medium_index(cs, "TEST:ECHO"));
}

TEST_CASE("system wavelengths outside a material's range are a CompileError", "[compile][media]") {
  MaterialLibrary lib;
  lib.add("TEST:NARROW", std::make_shared<const RangedMaterial>(WavelengthRange{0.55, 0.6}));
  lib.add("TEST:WIDE", std::make_shared<const RangedMaterial>(WavelengthRange{0.4861, 0.6563}));

  SECTION("element material") {
    const CompileError e = compile_error(singlet("TEST:NARROW"), lib);
    REQUIRE(rtt::model::test::has_error_at(e.diagnostics(), "/root/children/1/material"));
    REQUIRE_THAT(e.what(), ContainsSubstring("TEST:NARROW"));
    REQUIRE_THAT(e.what(), ContainsSubstring("0.4861"));
    REQUIRE_THAT(e.what(), ContainsSubstring("0.55"));
    REQUIRE_THAT(e.what(), ContainsSubstring("0.6"));
  }
  SECTION("environment medium") {
    System s = singlet("CONST:1.5");
    s.environment.medium = "TEST:NARROW";
    const CompileError e = compile_error(s, lib);
    REQUIRE(rtt::model::test::has_error_at(e.diagnostics(), "/environment/medium"));
  }
  SECTION("the range is closed: its end points are valid") {
    REQUIRE_NOTHROW(rtt::compile::compile(singlet("TEST:WIDE"), lib));
  }
  SECTION("unknown environment and unknown element material give a clean CompileError") {
    // Regression (review of #23): with no resolvable medium at all, media is empty while the
    // paths still refer to medium 0; the range check must not index past the end.
    System s = singlet("UNKNOWN:LENS");
    s.environment.medium = "UNKNOWN:ENV";
    const CompileError e = compile_error(s, lib);
    REQUIRE(rtt::model::test::has_error_at(e.diagnostics(), "/environment/medium"));
    REQUIRE(rtt::model::test::has_error_at(e.diagnostics(), "/root/children/1/material"));
  }
  SECTION("the error points at an element on the path, not at an unused one") {
    // Two lenses with the same material; only the second one is on the explicit path.
    System s = singlet("TEST:NARROW");
    Element second = rtt::model::test::element(s, 1);
    second.name = "L2";
    for (auto& surface : second.surfaces) {
      surface.id = rtt::model::SurfaceId("L2." + surface.id.str().substr(3));
    }
    s.root.children.push_back({second});
    s.paths = {{"second lens only",
                false,
                {{rtt::model::SurfaceId("L2.S1"), rtt::model::EventKind::Refract, 0},
                 {rtt::model::SurfaceId("L2.S2"), rtt::model::EventKind::Refract, 0}}}};
    const CompileError e = compile_error(s, lib);
    REQUIRE(rtt::model::test::has_error_at(e.diagnostics(), "/root/children/3/material"));
    REQUIRE_FALSE(rtt::model::test::has_error_at(e.diagnostics(), "/root/children/1/material"));
  }
  SECTION("media that no path uses are not checked") {
    // Explicit path that skips the lens: its material never meets a ray.
    System s = singlet("TEST:NARROW");
    s.paths = {{"skip lens",
                false,
                {{rtt::model::SurfaceId("STO"), rtt::model::EventKind::Transmit, 0},
                 {rtt::model::SurfaceId("IMG"), rtt::model::EventKind::Transmit, 0}}}};
    REQUIRE_NOTHROW(rtt::compile::compile(s, lib));
  }
}
