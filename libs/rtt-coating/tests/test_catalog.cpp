#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <complex>
#include <filesystem>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#include "rtt/coating/catalog.hpp"
#include "rtt/coating/thickness.hpp"
#include "rtt/coating/transfer_matrix.hpp"

using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinAbs;
using rtt::coating::CoatingCatalogError;
using rtt::coating::CoatingLibrary;
using rtt::coating::QuarterWaves;

namespace fs = std::filesystem;

namespace {

const fs::path kDemo = fs::path(RTT_CATALOG_DIR) / "coatings" / "demo.json";

/// Fresh empty directory under the temporary directory.
fs::path temp_dir(const std::string& name) {
  const fs::path dir = fs::temp_directory_path() / ("rtt_coating_test_" + name);
  fs::remove_all(dir);
  fs::create_directories(dir);
  return dir;
}

void write(const fs::path& file, const std::string& text) {
  std::ofstream(file, std::ios::binary) << text;
}

std::string catalog_text(const std::string& name, const std::string& coatings) {
  return R"({"format": "raytatouille-coatings", "schema_version": "0.1.0", "catalog": ")" + name +
         R"(", "coatings": [)" + coatings + "]}";
}

const std::string kOneLayer =
    R"({"name": "C", "layers": [{"material": "CONST:1.5", "thickness_um": 0.1}]})";

CoatingCatalogError parse_error(const std::string& text) {
  try {
    (void)rtt::coating::parse_coating_catalog(text, "test.json");
  } catch (const CoatingCatalogError& e) {
    return e;
  }
  FAIL("no CoatingCatalogError for " << text);
  return CoatingCatalogError("", "", "");
}

}  // namespace

TEST_CASE("the demo catalogue loads with its designs", "[coating][catalog]") {
  const rtt::coating::CoatingCatalog catalog = rtt::coating::load_coating_catalog(kDemo);
  REQUIRE(catalog.name == "DEMO");
  REQUIRE(catalog.coatings.size() == 3);
  const rtt::coating::CoatingDesign& ar = catalog.coatings[0];
  REQUIRE(ar.name == "AR_MGF2");
  REQUIRE(ar.design_wavelength_um == 0.55);
  REQUIRE(ar.layers.size() == 1);
  REQUIRE(ar.layers[0].material == "CONST:1.38");
  REQUIRE(ar.layers[0].thickness == rtt::coating::Thickness{QuarterWaves{1.0, 0.55}});
  const rtt::coating::CoatingDesign& v = catalog.coatings[1];
  REQUIRE(std::get<rtt::coating::PhysicalThickness>(v.layers[1].thickness).um == 0.0655);
  REQUIRE_FALSE(catalog.coatings[2].design_wavelength_um.has_value());
}

TEST_CASE("the AR design of the catalogue gives R = 1.26 % on n = 1.52", "[coating][catalog]") {
  // End to end with the reference case of #58: QWOT -> thickness (thickness.hpp) -> stack.
  CoatingLibrary library;
  library.add_catalog(kDemo);
  const auto ar = library.resolve("DEMO:AR_MGF2");
  const double d = rtt::coating::physical_thickness_um(ar->layers[0].thickness, 1.38);
  const std::vector<rtt::coating::Layer<double>> layers = {{1.38, d}};
  const auto a =
      rtt::coating::stack_amplitudes<double>(1.0, layers, std::complex<double>(1.52), 0.0, 0.55);
  const double r = (1.52 - 1.38 * 1.38) / (1.52 + 1.38 * 1.38);
  REQUIRE_THAT(std::norm(a.rs), WithinAbs(r * r, 1e-10));
}

TEST_CASE("CoatingLibrary resolves CATALOG:NAME", "[coating][catalog]") {
  CoatingLibrary library;
  REQUIRE(library.find("DEMO:AR_MGF2") == nullptr);
  library.add_catalog(kDemo);
  const auto ar = library.resolve("DEMO:AR_MGF2");
  REQUIRE(ar != nullptr);
  REQUIRE(library.resolve("DEMO:AR_MGF2") == ar);  // same object
  REQUIRE(library.find("DEMO:V_AR") != nullptr);
  REQUIRE(library.find("DEMO:NOPE") == nullptr);
  REQUIRE_THROWS_WITH(library.resolve("DEMO:NOPE"), ContainsSubstring("not found in catalog DEMO"));
  REQUIRE_THROWS_WITH(library.resolve("OTHER:AR"), ContainsSubstring("catalog OTHER not loaded"));
  REQUIRE_THROWS_WITH(library.resolve("AR_MGF2"), ContainsSubstring("expected CATALOG:NAME"));
  REQUIRE_THROWS_AS(library.resolve("X:Y"), rtt::coating::UnknownCoating);
  // A catalogue name can be loaded once.
  REQUIRE_THROWS_AS(library.add_catalog(kDemo), std::invalid_argument);
}

TEST_CASE("CoatingLibrary loads a directory, all or nothing", "[coating][catalog]") {
  const fs::path dir = temp_dir("directory");
  write(dir / "a.json", catalog_text("A", kOneLayer));
  write(dir / "b.JSON", catalog_text("B", kOneLayer));
  write(dir / "notes.txt", "not a catalogue");
  CoatingLibrary library;
  library.add_catalog(dir);
  REQUIRE(library.find("A:C") != nullptr);
  REQUIRE(library.find("B:C") != nullptr);

  const fs::path twice = temp_dir("twice");
  write(twice / "1.json", catalog_text("X", kOneLayer));
  write(twice / "2.json", catalog_text("X", kOneLayer));
  CoatingLibrary other;
  REQUIRE_THROWS_AS(other.add_catalog(twice), std::invalid_argument);
  REQUIRE(other.find("X:C") == nullptr);  // nothing registered

  const fs::path broken = temp_dir("broken");
  write(broken / "1.json", catalog_text("OK", kOneLayer));
  write(broken / "2.json", "{");
  REQUIRE_THROWS_AS(other.add_catalog(broken), CoatingCatalogError);
  REQUIRE(other.find("OK:C") == nullptr);

  REQUIRE_THROWS_AS(other.add_catalog(temp_dir("empty")), std::invalid_argument);
  fs::remove_all(dir);
  fs::remove_all(twice);
  fs::remove_all(broken);
}

TEST_CASE("catalogue errors name the file and the JSON pointer", "[coating][catalog]") {
  struct Case {
    std::string text;
    std::string pointer;
    std::string message;
  };
  const std::string prefix = R"({"format": "raytatouille-coatings", "schema_version": "0.1.0", )";
  const std::vector<Case> cases = {
      {"{", "", "JSON syntax error"},
      {"[]", "", "expected an object"},
      {prefix + R"("catalog": "D", "coatings": [], "extra": 1})", "/extra", "unknown key"},
      {R"({"format": "x", "schema_version": "0.1.0", "catalog": "D", "coatings": []})", "/format",
       "raytatouille-coatings"},
      {R"({"format": "raytatouille-coatings", "schema_version": "9.0.0", "catalog": "D",
           "coatings": []})",
       "/schema_version", "unsupported version"},
      {prefix + R"("coatings": []})", "", "missing key 'catalog'"},
      {prefix + R"("catalog": "A:B", "coatings": []})", "/catalog", "contains ':'"},
      {prefix + R"("catalog": "", "coatings": []})", "/catalog", "empty name"},
      {prefix + R"("catalog": "D", "coatings": {}})", "/coatings", "expected an array"},
      {catalog_text("D", R"({"layers": []})"), "/coatings/0", "missing key 'name'"},
      {catalog_text("D", kOneLayer + "," + kOneLayer), "/coatings/1/name", "duplicate"},
      {catalog_text("D", R"({"name": "C", "layers": []})"), "/coatings/0/layers",
       "at least one layer"},
      {catalog_text("D", R"({"name": "C", "layers": [{"material": "M", "qwot": 1,
                                                      "thickness_um": 0.1}]})"),
       "/coatings/0/layers/0", "exactly one"},
      {catalog_text("D", R"({"name": "C", "layers": [{"material": "M"}]})"), "/coatings/0/layers/0",
       "exactly one"},
      {catalog_text("D", R"({"name": "C", "layers": [{"material": "M", "thickness_um": -0.1}]})"),
       "/coatings/0/layers/0/thickness_um", ">= 0"},
      {catalog_text("D", R"({"name": "C", "layers": [{"material": "M", "thickness_um": "1"}]})"),
       "/coatings/0/layers/0/thickness_um", "expected a number"},
      {catalog_text("D", R"({"name": "C", "layers": [{"material": "M", "qwot": 1}]})"),
       "/coatings/0/layers/0/qwot", "design_wavelength_um"},
      {catalog_text("D", R"({"name": "C", "design_wavelength_um": 0,
                             "layers": [{"material": "M", "qwot": 1}]})"),
       "/coatings/0/design_wavelength_um", "> 0"},
      {catalog_text("D", R"({"name": "C", "layers": [{"material": "", "thickness_um": 0.1}]})"),
       "/coatings/0/layers/0/material", "empty material"},
      {catalog_text("D", R"({"name": "C", "layers": [{"material": "M", "thickness_um": 0.1,
                                                      "n": 1.5}]})"),
       "/coatings/0/layers/0/n", "unknown key"},
      {catalog_text("D", R"({"name": "C:D", "layers": [{"material": "M", "thickness_um": 0.1}]})"),
       "/coatings/0/name", "contains ':'"},
  };
  for (const Case& c : cases) {
    INFO(c.text);
    const CoatingCatalogError e = parse_error(c.text);
    REQUIRE(e.file() == "test.json");
    REQUIRE(e.pointer() == c.pointer);
    REQUIRE_THAT(e.what(), ContainsSubstring(c.message));
    REQUIRE_THAT(e.what(), ContainsSubstring("test.json"));
  }
}

TEST_CASE("a missing catalogue file is a CoatingCatalogError", "[coating][catalog]") {
  const fs::path missing = temp_dir("missing") / "nope.json";
  try {
    (void)rtt::coating::load_coating_catalog(missing);
    FAIL("no error");
  } catch (const CoatingCatalogError& e) {
    REQUIRE(e.pointer().empty());
    REQUIRE_THAT(e.what(), ContainsSubstring("cannot open"));
  }
}
