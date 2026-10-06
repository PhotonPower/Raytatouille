// Material library for a GUI (#85, G9): listing of AGF catalogues, catalogue alias, loading from
// memory, vectorized index. AGF records after the Ansys OpticStudio User Guide 2025 R1, "The AGF &
// BGF File Formats"; placeholders and value ranges as written by the manufacturers' files; units
// checked against the SCHOTT N-BK7 data sheet (docs/quellen.md).

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "rtt/material/agf.hpp"
#include "rtt/material/material.hpp"

namespace fs = std::filesystem;
using Catch::Matchers::ContainsSubstring;
using rtt::material::AgfCatalog;
using rtt::material::AgfClassRange;
using rtt::material::AgfError;
using rtt::material::AgfGlass;
using rtt::material::AgfOtherData;
using rtt::material::MaterialLibrary;
using rtt::material::UnknownMaterial;
using rtt::material::WavelengthRange;

namespace {

const fs::path kCatalogDir = RTT_CATALOG_DIR;

std::string read_bytes(const fs::path& file) {
  std::ifstream in(file, std::ios::binary);
  std::ostringstream text;
  text << in.rdbuf();
  return text.str();
}

/// The AgfError thrown by `f`; fails the test if there is none.
template <class F>
AgfError agf_error(F&& f) {
  try {
    f();
  } catch (const AgfError& e) {
    return e;
  }
  FAIL("no AgfError");
  return AgfError("", 0, "");
}

const AgfGlass& glass_named(const AgfCatalog& cat, const std::string& name) {
  for (const AgfGlass& g : cat.glasses) {
    if (g.name == name) return g;
  }
  FAIL("no glass " << name);
  return cat.glasses.front();
}

/// Invented catalogue with every record that the listing exposes, including the placeholders of
/// the manufacturers' files: "_" in MD, OD and the NM extras, "-" in the NM extras, -1 and 0 as
/// "melt frequency not given".
constexpr std::string_view kAll =
    "CC listing test\n"
    "NM ALL 2 517642.251 1.5168 64.17 1 3 5\n"
    "GC   free text comment  \n"
    "ED 7.1 8.3 2.51 -0.0009 0\n"
    "CD 1.039612120E+00 6.000698670E-03 2.317923440E-01 2.001791440E-02 1.010469450E+00 "
    "1.035606530E+02 0 0\n"
    "TD 1.86e-6 1.31e-8 -1.37e-11 4.34e-7 6.27e-10 0.17 20\n"
    "MD 82.00 0.21 610 _ 1.11\n"
    "OD 1.0 _ 0.0 1.0 2.3 -1\n"
    "LD 0.3 2.5\n"
    "IT 0.3 0.05 25\n"
    "IT 0.31 0.25 25\n"
    "BD _ _ _ _\n"
    "NM DASH 2 1 1.5 60 0 2 -\n"
    "CD 1 0.01 0 0 0 0\n"
    "NM MINUS 2 1 1.5 60 _ 1 -1\n"
    "CD 1 0.01 0 0 0 0\n"
    "NM ZERO 2 1 1.5 60 0 4 0\n"
    "CD 1 0.01 0 0 0 0\n"
    "NM SHORT 2 1 1.5 60\n"
    "CD 1 0.01 0 0 0 0\n";

}  // namespace

// ----------------------------------------------------------------- AGF records -----

TEST_CASE("AGF listing: NM extras, GC, MD, OD and IT as in the file", "[library]") {
  const AgfCatalog cat = rtt::material::parse_agf(kAll, "ALLCAT", "all.agf");
  REQUIRE(cat.glasses.size() == 5);
  const AgfGlass& g = cat.glasses[0];
  REQUIRE(g.exclude_substitution == 1);
  REQUIRE(g.status == 3);
  REQUIRE(g.melt_frequency == 5);
  REQUIRE(g.comment == "free text comment");
  REQUIRE(g.extra == std::vector<double>{7.1, 8.3, 2.51, -0.0009, 0.0});
  REQUIRE(g.thermal ==
          std::vector<double>{1.86e-6, 1.31e-8, -1.37e-11, 4.34e-7, 6.27e-10, 0.17, 20});
  REQUIRE(g.mechanical ==
          std::vector<std::optional<double>>{82.0, 0.21, 610.0, std::nullopt, 1.11});
  // OD keeps -1 as written (the format's "not available"); "_" is not available as well.
  REQUIRE(g.other == AgfOtherData{1.0, std::nullopt, AgfClassRange{0.0, 0.0},
                                  AgfClassRange{1.0, 1.0}, AgfClassRange{2.3, 2.3},
                                  AgfClassRange{-1.0, -1.0}});
  REQUIRE(g.range == WavelengthRange{0.3, 2.5});
  REQUIRE(g.transmission.size() == 2);
  REQUIRE(g.transmission[0].wavelength_um == 0.3);
  REQUIRE(g.transmission[0].transmittance == 0.05);
  REQUIRE(g.transmission[0].thickness_mm == 25.0);
  REQUIRE(g.transmission[1].wavelength_um == 0.31);
  REQUIRE(g.transmission[1].transmittance == 0.25);

  // Melt frequency "-", -1 and 0: not given; "_" in the NM extras: not given.
  REQUIRE(cat.glasses[1].status == 2);
  REQUIRE_FALSE(cat.glasses[1].melt_frequency.has_value());
  REQUIRE_FALSE(cat.glasses[2].exclude_substitution.has_value());
  REQUIRE(cat.glasses[2].status == 1);
  REQUIRE_FALSE(cat.glasses[2].melt_frequency.has_value());
  REQUIRE(cat.glasses[3].status == 4);
  REQUIRE_FALSE(cat.glasses[3].melt_frequency.has_value());
  // Without extras and without the optional records.
  const AgfGlass& s = cat.glasses[4];
  REQUIRE_FALSE(s.exclude_substitution.has_value());
  REQUIRE_FALSE(s.status.has_value());
  REQUIRE_FALSE(s.melt_frequency.has_value());
  REQUIRE(s.comment.empty());
  REQUIRE_FALSE(s.mechanical.has_value());
  REQUIRE_FALSE(s.other.has_value());
  REQUIRE(s.transmission.empty());
}

TEST_CASE("AGF listing of a manufacturer excerpt: NICF-V of NIKON-HIKARI", "[library]") {
  // tests/catalogs/nikon/nikon-hikari.agf (byte copy, #42), lines of NICF-V:
  //   NM NICF-V 6 1 1.433837 95.260792 0 0 0
  //   GC TCE value is available for 0 to 25 degrees Celsius.
  //   ED 1.840000000E+001 0.000000000E+000 3.180000000E+000 5.530000000E-002 0 0
  //   OD -1.00000 1.00000 0.00000 4.50000 2.30000 1.00000
  //   LD 1.84890000E-001 2.32542000E+000
  //   IT 1.57000E-001 9.95000E-001 1.00000E+001 ... IT 2.50000E+000 9.98000E-001 1.00000E+001
  //   (36 IT lines)
  MaterialLibrary lib;
  lib.add_catalog(kCatalogDir / "nikon" / "nikon-hikari.agf");
  REQUIRE(lib.catalogs() == std::vector<std::string>{"NIKON-HIKARI"});
  const auto cat = lib.catalog("NIKON-HIKARI");
  REQUIRE(cat->glasses.size() == 6);
  const AgfGlass g = glass_named(*cat, "NICF-V");  // copy: no reference into a temporary
  REQUIRE(g.formula == 6);
  REQUIRE(g.nd == 1.433837);
  REQUIRE(g.vd == 95.260792);
  REQUIRE(g.exclude_substitution == 0);
  REQUIRE(g.status == 0);
  REQUIRE_FALSE(g.melt_frequency.has_value());  // 0: not given
  REQUIRE(g.comment == "TCE value is available for 0 to 25 degrees Celsius.");
  REQUIRE(g.extra == std::vector<double>{18.4, 0.0, 3.18, 0.0553, 0.0, 0.0});
  REQUIRE(g.other == AgfOtherData{-1.0, AgfClassRange{1.0, 1.0}, AgfClassRange{0.0, 0.0},
                                  AgfClassRange{4.5, 4.5}, AgfClassRange{2.3, 2.3},
                                  AgfClassRange{1.0, 1.0}});
  REQUIRE_FALSE(g.mechanical.has_value());
  REQUIRE(g.range == WavelengthRange{0.18489, 2.32542});
  REQUIRE(g.transmission.size() == 36);
  REQUIRE(g.transmission.front().wavelength_um == 0.157);
  REQUIRE(g.transmission.front().transmittance == 0.995);
  REQUIRE(g.transmission.front().thickness_mm == 10.0);
  REQUIRE(g.transmission.back().wavelength_um == 2.5);
  REQUIRE(g.transmission.back().transmittance == 0.998);
  // A GC record with only tabs gives an empty comment.
  REQUIRE(glass_named(*cat, "E-LAKH1").comment.empty());
  REQUIRE_THROWS_AS((void)lib.catalog("NIKON"), UnknownMaterial);
}

TEST_CASE("AGF listing: resistance classes as ranges, as SCHOTT writes them", "[library]") {
  // OD records of the current SCHOTT catalogue ("schott glasses preferred and special
  // June-2025-B.AGF", schott.com), lines 318 (LAK9G15) and 404: CR as a class range "a-b".
  const std::string text =
      "CC c\n"
      "NM LAK9G15 2 691548.353 1.69064 54.76 0 1\n"
      "CD 1 2\n"
      "OD 32.6 1-2 2.0 53.0 1.3 4.3\n"
      "NM SECOND 2 1 1.5 60\n"
      "CD 1 2\n"
      "OD 3.1 2-3 2.0 3.4 2.2 3.0\n"
      // Line 6427 (KZFSN4, "old, was replaced by N-KZFS4"): "-" as relative cost.
      "NM KZFSN4 2 613443.320 1.6134 44.29 0 2 -\n"
      "CD 1 2\n"
      "OD - 3.0000 2.0000 52.3000 4.3000 4.3000\n";
  const AgfCatalog cat = rtt::material::parse_agf(text, "S", "s.agf");
  REQUIRE(cat.glasses[0].other == AgfOtherData{32.6, AgfClassRange{1.0, 2.0},
                                               AgfClassRange{2.0, 2.0}, AgfClassRange{53.0, 53.0},
                                               AgfClassRange{1.3, 1.3}, AgfClassRange{4.3, 4.3}});
  REQUIRE(cat.glasses[1].other->cr == AgfClassRange{2.0, 3.0});
  const AgfGlass& old = cat.glasses[2];
  REQUIRE(old.status == 2);
  REQUIRE_FALSE(old.melt_frequency.has_value());
  REQUIRE(old.other == AgfOtherData{std::nullopt, AgfClassRange{3.0, 3.0}, AgfClassRange{2.0, 2.0},
                                    AgfClassRange{52.3, 52.3}, AgfClassRange{4.3, 4.3},
                                    AgfClassRange{4.3, 4.3}});
  // A single number with an exponent is one class, not a range; -1 is kept as written.
  const AgfCatalog exp =
      rtt::material::parse_agf("CC c\nNM A 2 1 1.5 60\nCD 1 2\nOD 1 1e-5 -1 3 4 5\n", "S", "s.agf");
  REQUIRE(exp.glasses[0].other->cr == AgfClassRange{1e-5, 1e-5});
  REQUIRE(exp.glasses[0].other->fr == AgfClassRange{-1.0, -1.0});
  // "-" also stands for a missing class in OD, but not in MD or IT (error cases below).
  REQUIRE_FALSE(
      rtt::material::parse_agf("CC c\nNM A 2 1 1.5 60\nCD 1 2\nOD 1 - 3 4 5 6\n", "S", "s.agf")
          .glasses[0]
          .other->cr.has_value());
}

TEST_CASE("AGF listing: malformed or out-of-range extra data are errors with the line",
          "[library]") {
  struct Case {
    std::string text;
    std::size_t line;
    std::string message;
  };
  const std::string glass = "CC c\nNM A 2 1 1.5 60\nCD 1 2\n";
  const std::vector<Case> cases{
      {glass + "OD 1 2 3 4 5\n", 4, "OD"},         // OD needs exactly 6 values
      {glass + "OD 1 2 3 4 5 6 7\n", 4, "OD"},     // ... not more
      {glass + "OD 1 x 3 4 5 6\n", 4, "'x'"},      // not a number or "_"
      {glass + "OD x 1 3 4 5 6\n", 4, "'x'"},      // relative cost: a number or "_"
      {glass + "OD 1-2 1 3 4 5 6\n", 4, "'1-2'"},  // ... not a range
      {glass + "OD 1 2-1 3 4 5 6\n", 4, "'2-1'"},  // range with a > b
      {glass + "OD 1 1-x 3 4 5 6\n", 4, "'1-x'"},  // range of two numbers
      {glass + "OD 1 1- 3 4 5 6\n", 4, "'1-'"},
      {glass + "OD 1 -1-2 3 4 5 6\n", 4, "'-1-2'"},      // classes are >= 0
      {glass + "OD 1 1e-5-2 3 4 5 6\n", 4, "'1e-5-2'"},  // no exponent in a range
      {glass + "GC first\nGC second\n", 5, "second GC"},
      {glass + "MD 1 2 3 4\n", 4, "MD"},      // MD needs exactly 5 values
      {glass + "MD 1 2 3 4 - \n", 4, "'-'"},  // "-" only in the NM extras
      {glass + "IT 0.5 0.9\n", 4, "IT"},      // IT needs exactly 3 values
      {glass + "IT 0.5 _ 10\n", 4, "'_'"},    // no placeholder in IT
      {glass + "OD 1 2 3 4 5 6\nOD 1 2 3 4 5 6\n", 5, "second OD"},
      {glass + "MD 1 2 3 4 5\nMD 1 2 3 4 5\n", 5, "second MD"},
      {"CC c\nNM A 2 1 1.5 60 0 1 3 9\nCD 1 2\n", 2, "NM"},  // more than 3 extras
      {"CC c\nNM A 2 1 1.5 60 2\nCD 1 2\n", 2, "exclude"},   // exclude sub not 0/1
      {"CC c\nNM A 2 1 1.5 60 0 5\nCD 1 2\n", 2, "status"},  // status not 0..4
      {"CC c\nNM A 2 1 1.5 60 0 1 6\nCD 1 2\n", 2, "melt"},  // melt freq not -1..5
      {"CC c\nNM A 2 1 1.5 60 0 1 -2\nCD 1 2\n", 2, "melt"},
      {"CC c\nNM A 2 1 1.5 60 0 1.5\nCD 1 2\n", 2, "status"},  // not an integer
      {"CC c\nNM A 2 1 1.5 60 0 1 y\nCD 1 2\n", 2, "'y'"},
  };
  for (const auto& c : cases) {
    INFO(c.text);
    const AgfError e =
        agf_error([&] { (void)rtt::material::parse_agf(c.text, "TEST", "bad.agf"); });
    REQUIRE(e.line() == c.line);
    REQUIRE_THAT(e.what(), ContainsSubstring("bad.agf:" + std::to_string(c.line)));
    REQUIRE_THAT(e.what(), ContainsSubstring(c.message));
  }
}

// --------------------------------------------------------- alias and catalogues -----

TEST_CASE("catalogue alias: two files with the same stem load side by side", "[library]") {
  // tests/catalogs/schott.agf (N-BK7, F2) and tests/catalogs/m2/schott.agf (N-LAK9, N-SF5) are
  // both named SCHOTT; without an alias the second one fails as before (#85).
  MaterialLibrary lib;
  lib.add_catalog(kCatalogDir / "schott.agf");
  REQUIRE_THROWS_AS(lib.add_catalog(kCatalogDir / "m2" / "schott.agf"), std::invalid_argument);
  lib.add_catalog(kCatalogDir / "m2" / "schott.agf", "SCHOTT_M2");
  REQUIRE(lib.catalogs() == std::vector<std::string>{"SCHOTT", "SCHOTT_M2"});
  REQUIRE(lib.resolve("SCHOTT:N-BK7") != nullptr);
  REQUIRE(lib.resolve("SCHOTT_M2:N-LAK9") != nullptr);
  REQUIRE_THROWS_AS((void)lib.resolve("SCHOTT:N-LAK9"), UnknownMaterial);
  REQUIRE(lib.catalog("SCHOTT_M2")->name == "SCHOTT_M2");
  // The alias only renames the catalogue: same glass data and bit-identical indices as the
  // catalogue loaded under its own name.
  MaterialLibrary plain;
  plain.add_catalog(kCatalogDir / "m2" / "schott.agf");
  for (const double wl : {0.4, 0.5876, 1.0}) {
    REQUIRE(lib.resolve("SCHOTT_M2:N-LAK9")->index(wl, 25.0, 1.0) ==
            plain.resolve("SCHOTT:N-LAK9")->index(wl, 25.0, 1.0));
  }
}

TEST_CASE("catalogue names are case-sensitive like all references", "[library]") {
  // MaterialLibrary::resolve matches references exactly (material.hpp), so "schott" and
  // "SCHOTT" are different catalogues, not a collision.
  MaterialLibrary lib;
  lib.add_catalog(kCatalogDir / "schott.agf");
  lib.add_catalog(kCatalogDir / "schott.agf", "schott");
  REQUIRE(lib.catalogs() == std::vector<std::string>{"SCHOTT", "schott"});
  REQUIRE(lib.resolve("schott:N-BK7") != lib.resolve("SCHOTT:N-BK7"));
  REQUIRE(lib.resolve("schott:N-BK7")->index(0.5876, 20.0, 1.0) ==
          lib.resolve("SCHOTT:N-BK7")->index(0.5876, 20.0, 1.0));
}

TEST_CASE("catalogue alias: invalid names and collisions register nothing", "[library]") {
  MaterialLibrary lib;
  lib.add_catalog(kCatalogDir / "schott.agf");
  for (const std::string name : {"", "A:B", "A B", "A\tB", "CONST", "SCHOTT"}) {
    INFO("'" << name << "'");
    REQUIRE_THROWS_AS(lib.add_catalog(kCatalogDir / "m2" / "schott.agf", name),
                      std::invalid_argument);
  }
  // A directory with a name is ambiguous.
  REQUIRE_THROWS_AS(lib.add_catalog(kCatalogDir / "m2", std::string("M2")), std::invalid_argument);
  REQUIRE(lib.catalogs() == std::vector<std::string>{"SCHOTT"});
  REQUIRE_THROWS_AS((void)lib.resolve("M2:N-LAK9"), UnknownMaterial);
  // A glass reference already in use: nothing of the catalogue is registered.
  lib.add("A2:N-LAK9", lib.resolve("CONST:1.5"));
  REQUIRE_THROWS_AS(lib.add_catalog(kCatalogDir / "m2" / "schott.agf", std::string("A2")),
                    std::invalid_argument);
  REQUIRE(lib.catalogs() == std::vector<std::string>{"SCHOTT"});
  REQUIRE_THROWS_AS((void)lib.resolve("A2:N-SF5"), UnknownMaterial);
  // add_catalog_text checks its name like an alias.
  for (const std::string name : {"", "A:B", "A B", "CONST", "SCHOTT"}) {
    INFO("'" << name << "'");
    REQUIRE_THROWS_AS(lib.add_catalog_text("CC c\nNM A 2 1 1.5 60\nCD 1 2\n", name),
                      std::invalid_argument);
  }
  REQUIRE(lib.catalogs() == std::vector<std::string>{"SCHOTT"});
}

TEST_CASE("catalogue from memory: same listing and indices as from the file", "[library]") {
  for (const fs::path& file : {kCatalogDir / "schott.agf", kCatalogDir / "utf16" / "schott.agf",
                               kCatalogDir / "nikon" / "nikon-hikari.agf"}) {
    INFO(file.string());
    MaterialLibrary from_file;
    from_file.add_catalog(file, "CAT");
    MaterialLibrary from_memory;
    from_memory.add_catalog_text(read_bytes(file), "CAT", "memory.agf");
    const auto a = from_file.catalog("CAT");
    const auto b = from_memory.catalog("CAT");
    REQUIRE(b->file == "memory.agf");
    REQUIRE(a->glasses.size() == b->glasses.size());
    for (std::size_t i = 0; i < a->glasses.size(); ++i) {
      const AgfGlass& ga = a->glasses[i];
      const AgfGlass& gb = b->glasses[i];
      REQUIRE(ga.name == gb.name);
      REQUIRE(ga.coefficients == gb.coefficients);
      REQUIRE(ga.thermal == gb.thermal);
      REQUIRE(ga.other == gb.other);
      REQUIRE(ga.transmission.size() == gb.transmission.size());
      const std::string ref = "CAT:" + ga.name;
      REQUIRE(from_file.resolve(ref)->index(0.55, 30.0, 0.9) ==
              from_memory.resolve(ref)->index(0.55, 30.0, 0.9));
    }
  }
  // Errors name the given source; the name is required and checked like an alias.
  MaterialLibrary lib;
  const AgfError e = agf_error([&] { lib.add_catalog_text("CC c\nXX 1\n", "MEM", "pasted text"); });
  REQUIRE_THAT(e.what(), ContainsSubstring("pasted text:2"));
  REQUIRE_THROWS_AS(lib.add_catalog_text("CC c\n", "", "x"), std::invalid_argument);
  REQUIRE(lib.catalogs().empty());
}

// ------------------------------------------------------------- vectorized index -----

TEST_CASE("index_many is bit-identical with the scalar index", "[library]") {
  MaterialLibrary lib;
  lib.add_catalog(kCatalogDir / "schott.agf");
  lib.add_catalog(kCatalogDir / "nikon" / "nikon-hikari.agf");
  std::vector<double> wavelengths;
  for (int i = 0; i < 1000; ++i) wavelengths.push_back(0.4 + 0.3 * i / 999.0);
  for (const std::string ref :
       {"SCHOTT:N-BK7", "SCHOTT:F2", "NIKON-HIKARI:NICF-V", "NIKON-HIKARI:E-LAKH1",
        "NIKON-HIKARI:J-SFH1", "AIR", "VACUUM", "CONST:1.5,1e-3"}) {
    INFO(ref);
    const auto material = lib.resolve(ref);
    const std::vector<rtt::math::Complex> many =
        rtt::material::index_many(*material, wavelengths, 23.5, 0.95);
    REQUIRE(many.size() == wavelengths.size());
    for (std::size_t i = 0; i < wavelengths.size(); ++i) {
      const rtt::math::Complex one = material->index(wavelengths[i], 23.5, 0.95);
      REQUIRE(std::isfinite(one.real()));
      REQUIRE(many[i] == one);
    }
  }
  REQUIRE(rtt::material::index_many(*lib.resolve("AIR"), {}, 20.0, 1.0).empty());
}
