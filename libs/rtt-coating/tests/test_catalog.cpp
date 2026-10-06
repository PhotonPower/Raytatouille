#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <system_error>
#include <utility>
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

/// Fresh empty directory under the temporary directory, removed again at the end of the test.
/// The name is unique (clock ticks), so tests in parallel worktrees do not collide.
class TempDir {
 public:
  explicit TempDir(const std::string& name)
      : path_(fs::temp_directory_path() /
              ("rtt_coating_test_" + name + "_" +
               std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()))) {
    fs::create_directories(path_);
  }
  TempDir(const TempDir&) = delete;
  TempDir& operator=(const TempDir&) = delete;
  TempDir(TempDir&&) = delete;
  TempDir& operator=(TempDir&&) = delete;
  ~TempDir() {
    std::error_code ignored;
    fs::remove_all(path_, ignored);
  }
  [[nodiscard]] const fs::path& path() const noexcept { return path_; }

 private:
  fs::path path_;
};

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
  const TempDir directory("directory");
  const fs::path& dir = directory.path();
  write(dir / "a.json", catalog_text("A", kOneLayer));
  write(dir / "b.JSON", catalog_text("B", kOneLayer));
  write(dir / "notes.txt", "not a catalogue");
  CoatingLibrary library;
  library.add_catalog(dir);
  REQUIRE(library.find("A:C") != nullptr);
  REQUIRE(library.find("B:C") != nullptr);

  const TempDir twice_dir("twice");
  const fs::path& twice = twice_dir.path();
  write(twice / "1.json", catalog_text("X", kOneLayer));
  write(twice / "2.json", catalog_text("X", kOneLayer));
  CoatingLibrary other;
  REQUIRE_THROWS_AS(other.add_catalog(twice), std::invalid_argument);
  REQUIRE(other.find("X:C") == nullptr);  // nothing registered

  const TempDir broken_dir("broken");
  const fs::path& broken = broken_dir.path();
  write(broken / "1.json", catalog_text("OK", kOneLayer));
  write(broken / "2.json", "{");
  REQUIRE_THROWS_AS(other.add_catalog(broken), CoatingCatalogError);
  REQUIRE(other.find("OK:C") == nullptr);

  const TempDir empty("empty");
  REQUIRE_THROWS_AS(other.add_catalog(empty.path()), std::invalid_argument);

  // A directory with a catalogue that the library already has: nothing of it is added.
  const TempDir again_dir("again");
  write(again_dir.path() / "a.json", catalog_text("A", kOneLayer));
  write(again_dir.path() / "new.json", catalog_text("NEW", kOneLayer));
  REQUIRE_THROWS_AS(library.add_catalog(again_dir.path()), std::invalid_argument);
  REQUIRE(library.find("NEW:C") == nullptr);
}

TEST_CASE("CoatingLibrary::add checks a catalogue built in code like a file",
          "[coating][catalog]") {
  using rtt::coating::CoatingCatalog;
  using rtt::coating::CoatingDesign;
  using rtt::coating::LayerSpec;
  using rtt::coating::PhysicalThickness;
  const LayerSpec ok{"CONST:1.5", PhysicalThickness{0.1}};
  const auto design = [&](std::string name, std::vector<LayerSpec> layers,
                          std::optional<double> l0 = std::nullopt) {
    return CoatingDesign{std::move(name), "", l0, std::move(layers)};
  };
  const std::vector<CoatingCatalog> bad = {
      {"", {design("C", {ok})}},
      {"A:B", {design("C", {ok})}},
      {"K", {design("", {ok})}},
      {"K", {design("C:D", {ok})}},
      {"K", {design("C", {ok}), design("C", {ok})}},
      {"K", {design("C", {})}},
      {"K", {design("C", {LayerSpec{"", PhysicalThickness{0.1}}})}},
      {"K", {design("C", {LayerSpec{"M", PhysicalThickness{-0.1}}})}},
      {"K", {design("C", {LayerSpec{"M", PhysicalThickness{std::nan("")}}})}},
      {"K", {design("C", {LayerSpec{"M", QuarterWaves{-1.0, 0.55}}})}},
      {"K", {design("C", {LayerSpec{"M", QuarterWaves{1.0, 0.0}}})}},
      {"K", {design("C", {ok}, 0.0)}},
  };
  for (std::size_t i = 0; i < bad.size(); ++i) {
    INFO("case " << i);
    CoatingLibrary library;
    REQUIRE_THROWS_AS(library.add(bad[i]), std::invalid_argument);
    REQUIRE(library.find("K:C") == nullptr);
  }
  CoatingLibrary library;
  REQUIRE_NOTHROW(library.add(CoatingCatalog{"K", {design("C", {ok})}}));
  REQUIRE(library.find("K:C") != nullptr);
}

TEST_CASE("catalogue errors name the file and the JSON pointer", "[coating][catalog]") {
  struct Case {
    std::string text;
    std::string pointer;
    std::string message;
  };
  const std::string prefix = R"({"format": "raytatouille-coatings", "schema_version": "0.1.0", )";
  const std::vector<Case> cases = {
      {"{", "", "invalid JSON"},
      // A number that overflows double: the only non-finite JSON number (nlohmann out_of_range).
      {catalog_text("D", R"({"name": "C", "layers": [{"material": "M", "qwot": 1e400}]})"), "",
       "invalid JSON"},
      {"[]", "", "expected an object"},
      // Duplicate keys (#68): the pointer names the second occurrence.
      {prefix + R"("catalog": "D", "catalog": "E", "coatings": []})", "/catalog",
       "duplicate key 'catalog'"},
      {catalog_text("D", R"({"name": "C", "layers": [{"material": "M", "thickness_um": -1,
                                                      "thickness_um": 0.1}]})"),
       "/coatings/0/layers/0/thickness_um", "duplicate key 'thickness_um'"},
      {catalog_text("D", R"({"name": "C", "a/b~c": 1, "a/b~c": 2, "layers": []})"),
       "/coatings/0/a~1b~0c", "duplicate key 'a/b~c'"},
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
      // RFC 6901: '/' and '~' in a key are escaped as ~1 and ~0 in the pointer.
      {catalog_text("D", R"({"name": "C", "layers": [{"material": "M", "thickness_um": 0.1,
                                                      "a/b~c": 1}]})"),
       "/coatings/0/layers/0/a~1b~0c", "unknown key 'a/b~c'"},
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
  const TempDir dir("missing");
  const fs::path missing = dir.path() / "nope.json";
  try {
    (void)rtt::coating::load_coating_catalog(missing);
    FAIL("no error");
  } catch (const CoatingCatalogError& e) {
    REQUIRE(e.pointer().empty());
    REQUIRE_THAT(e.what(), ContainsSubstring("cannot open"));
  }
}
