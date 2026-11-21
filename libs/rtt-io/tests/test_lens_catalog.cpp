// Lens catalogues (ADR 0053, points 2, 3 and 6.5): the demo catalogue and every error with the
// file and the JSON pointer.

#include <catch2/catch_test_macros.hpp>
#include <optional>
#include <string>
#include <vector>

#include "rtt/io/lens_catalog.hpp"

namespace {

using rtt::io::LensCatalog;
using rtt::io::LensCatalogError;
using rtt::io::LensEntry;
using rtt::io::LensKind;

const std::string kDemo = std::string(RTT_CATALOG_DIR) + "/lenses/demo.json";

/// A catalogue "D" with the given lens objects.
std::string catalog_text(const std::string& lenses) {
  return R"({"format": "raytatouille-lenses", "schema_version": "0.1.0", "catalog": "D",
             "lenses": [)" +
         lenses + "]}";
}

const std::string kSinglet =
    R"({"name": "S", "kind": "singlet", "radii": [50.0, null], "thicknesses": [4.0],
        "materials": ["CONST:1.5"]})";

}  // namespace

TEST_CASE("the demo lens catalogue loads with its three kinds", "[io][lens_catalog]") {
  const LensCatalog catalog = rtt::io::load_lens_catalog(kDemo);
  CHECK(catalog.name == "DEMO");
  REQUIRE(catalog.lenses.size() == 3);
  const LensEntry& pcx = catalog.lenses[0];
  CHECK(pcx.name == "PCX-25-100");
  CHECK(pcx.kind == LensKind::Singlet);
  CHECK(pcx.radii == std::vector<std::optional<double>>{51.68, std::nullopt});
  CHECK(pcx.thicknesses == std::vector<double>{4.0});
  CHECK(pcx.gaps.empty());
  CHECK(pcx.materials == std::vector<std::string>{"CONST:1.5168"});
  CHECK(pcx.diameter == 25.4);
  CHECK_FALSE(pcx.description.empty());
  const LensEntry& ach = catalog.lenses[1];
  CHECK(ach.kind == LensKind::Doublet);
  CHECK(ach.radii == std::vector<std::optional<double>>{61.5, -35.2, -136.9});
  CHECK(ach.thicknesses == std::vector<double>{6.0, 3.0});
  CHECK(ach.materials == std::vector<std::string>{"SCHOTT:N-BK7", "SCHOTT:F2"});
  const LensEntry& trp = catalog.lenses[2];
  CHECK(trp.kind == LensKind::Triplet);
  CHECK(trp.radii ==
        std::vector<std::optional<double>>{23.713, 7331.288, -24.456, 21.896, 86.759, -20.4942});
  CHECK(trp.materials ==
        std::vector<std::string>{"SCHOTT:N-LAK9", "SCHOTT:N-SF5", "SCHOTT:N-LAK9"});
  CHECK(ach.diameter == 25.4);
  CHECK(trp.thicknesses == std::vector<double>{4.831, 0.975, 3.127});
  CHECK(trp.gaps == std::vector<double>{5.86, 4.822});
  CHECK_FALSE(trp.diameter.has_value());
  // Parsing the same text again gives the same catalogue.
  CHECK(rtt::io::parse_lens_catalog(catalog_text(kSinglet)) ==
        rtt::io::parse_lens_catalog(catalog_text(kSinglet)));
}

TEST_CASE("the names of the lens kinds", "[io][lens_catalog]") {
  CHECK(rtt::io::lens_kind_name(LensKind::Singlet) == "singlet");
  CHECK(rtt::io::lens_kind_name(LensKind::Doublet) == "doublet");
  CHECK(rtt::io::lens_kind_name(LensKind::Triplet) == "triplet");
}

TEST_CASE("lens catalogue errors name the file and the JSON pointer", "[io][lens_catalog]") {
  struct Case {
    std::string text;
    std::string pointer;
    std::string message;
  };
  const std::string prefix = R"({"format": "raytatouille-lenses", "schema_version": "0.1.0", )";
  const auto lens = [](const std::string& fields) {
    return catalog_text(R"({"name": "L", )" + fields + "}");
  };
  const std::vector<Case> cases = {
      {"{", "", "invalid JSON"},
      {"[]", "", "expected an object"},
      {prefix + R"("catalog": "D", "catalog": "E", "lenses": []})", "/catalog", "duplicate key"},
      {prefix + R"("catalog": "D", "lenses": [], "extra": 1})", "/extra", "unknown key"},
      {R"({"format": "x", "schema_version": "0.1.0", "catalog": "D", "lenses": []})", "/format",
       "raytatouille-lenses"},
      {R"({"format": "raytatouille-lenses", "schema_version": "0.2.0", "catalog": "D",
           "lenses": []})",
       "/schema_version", "unsupported version"},
      {prefix + R"("lenses": []})", "", "missing key 'catalog'"},
      {prefix + R"("catalog": "A:B", "lenses": []})", "/catalog", "contains ':'"},
      {prefix + R"("catalog": "", "lenses": []})", "/catalog", "empty name"},
      {prefix + R"("catalog": "D", "lenses": {}})", "/lenses", "expected an array"},
      {catalog_text(R"({"kind": "singlet"})"), "/lenses/0", "missing key 'name'"},
      {catalog_text(kSinglet + "," + kSinglet), "/lenses/1/name", "duplicate lens name"},
      {lens(R"("kind": "quartet", "radii": [], "thicknesses": [], "materials": [])"),
       "/lenses/0/kind", "unknown kind"},
      {lens(R"("kind": "singlet", "radii": [1.0, 2.0], "thicknesses": [1.0],
               "materials": ["M"], "colour": "red")"),
       "/lenses/0/colour", "unknown key"},
      {lens(R"("kind": "singlet", "radii": [1.0], "thicknesses": [1.0], "materials": ["M"])"),
       "/lenses/0/radii", "expected 2 entries, got 1"},
      {lens(R"("kind": "doublet", "radii": [1.0, 2.0, 3.0], "thicknesses": [1.0],
               "materials": ["M", "N"])"),
       "/lenses/0/thicknesses", "expected 2 entries"},
      {lens(R"("kind": "doublet", "radii": [1.0, 2.0, 3.0], "thicknesses": [1.0, 1.0],
               "materials": ["M"])"),
       "/lenses/0/materials", "expected 2 entries"},
      {lens(R"("kind": "triplet", "radii": [1, 2, 3, 4, 5, 6], "thicknesses": [1, 1, 1],
               "materials": ["M", "N", "M"])"),
       "/lenses/0", "missing key 'gaps'"},
      {lens(R"("kind": "triplet", "radii": [1, 2, 3, 4, 5, 6], "thicknesses": [1, 1, 1],
               "gaps": [1], "materials": ["M", "N", "M"])"),
       "/lenses/0/gaps", "expected 2 entries"},
      {lens(R"("kind": "singlet", "radii": [1.0, 2.0], "thicknesses": [1.0], "gaps": [],
               "materials": ["M"])"),
       "/lenses/0/gaps", "only a triplet"},
      {lens(R"("kind": "doublet", "radii": [1.0, 2.0, 3.0], "thicknesses": [1.0, 1.0],
               "gaps": [1.0], "materials": ["M", "N"])"),
       "/lenses/0/gaps", "only a triplet"},
      {lens(R"("kind": "singlet", "radii": [-0.0, 2.0], "thicknesses": [1.0], "materials": ["M"])"),
       "/lenses/0/radii/0", "radius != 0"},
      {lens(R"("kind": "singlet", "radii": [1.0, 2.0], "thicknesses": [-1.0], "materials": ["M"])"),
       "/lenses/0/thicknesses/0", "> 0"},
      {lens(R"("kind": "singlet", "radii": [0.0, 2.0], "thicknesses": [1.0], "materials": ["M"])"),
       "/lenses/0/radii/0", "radius != 0"},
      {lens(R"("kind": "singlet", "radii": [1.0, "x"], "thicknesses": [1.0], "materials": ["M"])"),
       "/lenses/0/radii/1", "number or null"},
      {lens(R"("kind": "singlet", "radii": [1.0, 2.0], "thicknesses": [0.0], "materials": ["M"])"),
       "/lenses/0/thicknesses/0", "> 0"},
      {lens(R"("kind": "triplet", "radii": [1, 2, 3, 4, 5, 6], "thicknesses": [1, 1, 1],
               "gaps": [1, -1], "materials": ["M", "N", "M"])"),
       "/lenses/0/gaps/1", ">= 0"},
      {lens(R"("kind": "singlet", "radii": [1.0, 2.0], "thicknesses": [1.0], "materials": [""])"),
       "/lenses/0/materials/0", "empty material"},
      {lens(R"("kind": "singlet", "radii": [1.0, 2.0], "thicknesses": [1.0], "materials": ["M"],
               "diameter": 0)"),
       "/lenses/0/diameter", "> 0"},
      {catalog_text(R"({"name": "A:B", "kind": "singlet", "radii": [1.0, 2.0],
                        "thicknesses": [1.0], "materials": ["M"]})"),
       "/lenses/0/name", "contains ':'"},
  };
  for (const Case& c : cases) {
    INFO(c.text);
    try {
      static_cast<void>(rtt::io::parse_lens_catalog(c.text, "cat.json"));
      FAIL("no error");
    } catch (const LensCatalogError& e) {
      CHECK(e.file() == "cat.json");
      CHECK(e.pointer() == c.pointer);
      CHECK(std::string(e.what()).find(c.message) != std::string::npos);
      CHECK(std::string(e.what()).rfind("cat.json: ", 0) == 0);
    }
  }
}

TEST_CASE("a missing lens catalogue file is a LensCatalogError", "[io][lens_catalog]") {
  try {
    static_cast<void>(rtt::io::load_lens_catalog(std::string(RTT_CATALOG_DIR) + "/lenses/none"));
    FAIL("no error");
  } catch (const LensCatalogError& e) {
    CHECK(e.pointer().empty());
    CHECK(std::string(e.what()).find("cannot open") != std::string::npos);
  }
}
