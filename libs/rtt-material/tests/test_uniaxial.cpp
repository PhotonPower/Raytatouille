// Uniaxial crystal material (#129; ADR 0014, point 4; ADR 0026, point 1): the principal indices
// n_O and n_E come from two ordinary materials. The tests check the mechanism (resolution,
// ranges, errors); no crystal data are asserted.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <filesystem>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>

#include "rtt/material/material.hpp"
#include "rtt/material/uniaxial.hpp"

namespace fs = std::filesystem;
using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::StartsWith;
using rtt::material::Material;
using rtt::material::MaterialLibrary;
using rtt::material::UniaxialMaterial;
using rtt::material::UnknownMaterial;
using rtt::material::WavelengthRange;
using rtt::math::Complex;

namespace {

const fs::path kCatalogDir = RTT_CATALOG_DIR;

/// Constant index on a given wavelength range (um).
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

TEST_CASE("uniaxial material from two CONST references", "[uniaxial]") {
  const MaterialLibrary lib;
  const UniaxialMaterial m = lib.resolve_uniaxial("CONST:1.6584", "CONST:1.4864");
  // The parts are the objects of resolve() (same cache), in the given order.
  REQUIRE(&m.ordinary() == lib.resolve("CONST:1.6584").get());
  REQUIRE(&m.extraordinary() == lib.resolve("CONST:1.4864").get());
  REQUIRE(m.n_ordinary(0.5876, 20.0, 1.0) == Complex(1.6584, 0.0));
  REQUIRE(m.n_extraordinary(0.5876, 20.0, 1.0) == Complex(1.4864, 0.0));
  REQUIRE_FALSE(m.wavelength_range_um().has_value());

  // A value type: a copy shares the two materials.
  const UniaxialMaterial copy = m;
  REQUIRE(&copy.ordinary() == &m.ordinary());
  REQUIRE(&copy.extraordinary() == &m.extraordinary());
}

TEST_CASE("uniaxial material passes wavelength, temperature and pressure on", "[uniaxial]") {
  // AIR (Ciddor) depends on all three, CONST on none: each part must get exactly the arguments
  // of the call, in this order.
  const MaterialLibrary lib;
  const auto air = lib.resolve("AIR");
  const UniaxialMaterial o_air = lib.resolve_uniaxial("AIR", "CONST:1.5");
  const UniaxialMaterial e_air = lib.resolve_uniaxial("CONST:1.5", "AIR");
  REQUIRE(o_air.n_ordinary(0.55, 35.0, 0.8) == air->index(0.55, 35.0, 0.8));
  REQUIRE(e_air.n_extraordinary(0.55, 35.0, 0.8) == air->index(0.55, 35.0, 0.8));
  REQUIRE(o_air.n_ordinary(0.55, 35.0, 0.8) != air->index(0.55, 20.0, 1.0));
  REQUIRE(o_air.n_ordinary(0.55, 35.0, 0.8) != air->index(0.55, 0.8, 35.0));  // not swapped
  REQUIRE(o_air.n_ordinary(0.55, 35.0, 0.8) != air->index(0.65, 35.0, 0.8));
  REQUIRE(e_air.n_extraordinary(0.55, 35.0, 0.8) != air->index(0.55, 20.0, 1.0));
}

TEST_CASE("uniaxial material passes the absorption through", "[uniaxial]") {
  // ADR 0026, point 1: kappa != 0 is an error of compile (crystal.absorbing per wavelength, #131),
  // so the material neither checks nor drops it.
  const MaterialLibrary lib;
  const UniaxialMaterial m = lib.resolve_uniaxial("CONST:1.6,0.01", "CONST:1.5,0.002");
  REQUIRE(m.n_ordinary(0.5876, 20.0, 1.0) == Complex(1.6, 0.01));
  REQUIRE(m.n_extraordinary(0.5876, 20.0, 1.0) == Complex(1.5, 0.002));
}

TEST_CASE("uniaxial wavelength range is the intersection of both parts", "[uniaxial]") {
  MaterialLibrary lib;
  lib.add("TEST:A", std::make_shared<const RangedMaterial>(1.6, WavelengthRange{0.4, 1.0}));
  lib.add("TEST:B", std::make_shared<const RangedMaterial>(1.5, WavelengthRange{0.5, 2.0}));
  lib.add("TEST:C", std::make_shared<const RangedMaterial>(1.5, WavelengthRange{0.6, 0.7}));
  lib.add("TEST:D", std::make_shared<const RangedMaterial>(1.5, WavelengthRange{0.3, 0.5}));

  REQUIRE(lib.resolve_uniaxial("TEST:A", "TEST:B").wavelength_range_um() ==
          WavelengthRange{0.5, 1.0});
  REQUIRE(lib.resolve_uniaxial("TEST:B", "TEST:A").wavelength_range_um() ==
          WavelengthRange{0.5, 1.0});
  // One unbounded part: the range of the other one.
  REQUIRE(lib.resolve_uniaxial("CONST:1.6", "TEST:B").wavelength_range_um() ==
          WavelengthRange{0.5, 2.0});
  REQUIRE(lib.resolve_uniaxial("TEST:A", "CONST:1.5").wavelength_range_um() ==
          WavelengthRange{0.4, 1.0});
  // Disjoint ranges: an empty range (min > max), so that no wavelength is contained.
  const auto empty = lib.resolve_uniaxial("TEST:C", "TEST:D").wavelength_range_um();
  REQUIRE(empty == WavelengthRange{0.6, 0.5});
  for (const double wavelength : {0.3, 0.5, 0.55, 0.6, 0.7}) {
    REQUIRE_FALSE(empty->contains(wavelength));
  }
  // Touching ranges [0.3, 0.5] and [0.5, 2.0]: exactly one wavelength.
  const auto point = lib.resolve_uniaxial("TEST:D", "TEST:B").wavelength_range_um();
  REQUIRE(point == WavelengthRange{0.5, 0.5});
  REQUIRE(point->contains(0.5));
}

TEST_CASE("uniaxial resolution errors name the failing part", "[uniaxial]") {
  const MaterialLibrary lib;
  REQUIRE_THROWS_AS((void)lib.resolve_uniaxial("NOPE:X", "CONST:1.5"), UnknownMaterial);
  REQUIRE_THROWS_WITH((void)lib.resolve_uniaxial("NOPE:X", "CONST:1.5"),
                      StartsWith("ordinary: ") && ContainsSubstring("NOPE"));
  REQUIRE_THROWS_AS((void)lib.resolve_uniaxial("CONST:1.6", "CONST:-1"), UnknownMaterial);
  REQUIRE_THROWS_WITH((void)lib.resolve_uniaxial("CONST:1.6", "CONST:-1"),
                      StartsWith("extraordinary: ") && ContainsSubstring("CONST:-1"));
  // Both parts invalid: the ordinary one is reported (documented order).
  REQUIRE_THROWS_WITH((void)lib.resolve_uniaxial("NOPE:X", "NOPE:Y"),
                      StartsWith("ordinary: ") && ContainsSubstring("NOPE:X"));
  // The message after the prefix is the one of resolve().
  std::string plain;
  try {
    (void)lib.resolve("NOPE:X");
  } catch (const UnknownMaterial& e) {
    plain = e.what();
  }
  REQUIRE_FALSE(plain.empty());
  REQUIRE_THROWS_WITH((void)lib.resolve_uniaxial("NOPE:X", "CONST:1.5"), "ordinary: " + plain);
}

TEST_CASE("uniaxial material rejects null parts", "[uniaxial]") {
  const auto glass = std::make_shared<const RangedMaterial>(1.5, WavelengthRange{0.4, 1.0});
  REQUIRE_THROWS_AS((void)UniaxialMaterial(nullptr, glass), std::invalid_argument);
  REQUIRE_THROWS_AS((void)UniaxialMaterial(glass, nullptr), std::invalid_argument);
}

TEST_CASE("uniaxial crystal as the pair X / X-E of one AGF catalogue", "[uniaxial][agf]") {
  // tests/catalogs/crystals/demo_crystals.agf: synthetic DEMO (n_O) and DEMO-E (n_E), the
  // convention of the Zemax birefringent.agf (ADR 0026, "Kataloge").
  MaterialLibrary lib;
  lib.add_catalog(kCatalogDir / "crystals" / "demo_crystals.agf");
  const UniaxialMaterial m = lib.resolve_uniaxial("DEMO_CRYSTALS:DEMO", "DEMO_CRYSTALS:DEMO-E");
  const auto o = lib.resolve("DEMO_CRYSTALS:DEMO");
  const auto e = lib.resolve("DEMO_CRYSTALS:DEMO-E");
  REQUIRE(&m.ordinary() == o.get());
  REQUIRE(&m.extraordinary() == e.get());
  for (const double wavelength : {0.45, 0.5876, 1.0}) {
    INFO(wavelength);
    // Bitwise the indices of the two catalogue glasses (same objects, same call).
    REQUIRE(m.n_ordinary(wavelength, 25.0, 1.0) == o->index(wavelength, 25.0, 1.0));
    REQUIRE(m.n_extraordinary(wavelength, 25.0, 1.0) == e->index(wavelength, 25.0, 1.0));
    REQUIRE(m.n_ordinary(wavelength, 25.0, 1.0).real() >
            m.n_extraordinary(wavelength, 25.0, 1.0).real());
  }
  // LD records: DEMO 0.4 ... 2.0 um, DEMO-E 0.35 ... 1.5 um.
  REQUIRE(m.wavelength_range_um() == WavelengthRange{0.4, 1.5});
}
