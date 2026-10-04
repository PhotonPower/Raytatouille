#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <variant>
#include <vector>

#include "rtt/material/agf.hpp"
#include "rtt/material/air.hpp"
#include "rtt/material/dispersion.hpp"
#include "rtt/material/material.hpp"

namespace fs = std::filesystem;
using Catch::Matchers::ContainsSubstring;
using rtt::material::AgfCatalog;
using rtt::material::AgfError;
using rtt::material::AgfGlass;
using rtt::material::CatalogMaterial;
using rtt::material::decode_agf_text;
using rtt::material::DispersionMaterial;
using rtt::material::load_agf;
using rtt::material::MaterialLibrary;
using rtt::material::parse_agf;
using rtt::material::SchottCoefficients;
using rtt::material::Sellmeier1Coefficients;
using rtt::material::UnknownMaterial;
using rtt::material::WavelengthRange;

namespace {

const fs::path kCatalogDir = RTT_CATALOG_DIR;
const fs::path kSchottFile = kCatalogDir / "schott.agf";
const fs::path kSchottUtf16 = kCatalogDir / "utf16" / "schott.agf";

// N-BK7 as in tests/catalogs/schott.agf (free SCHOTT catalogue, docs/quellen.md).
const Sellmeier1Coefficients kNbk7{{1.039612120, 2.317923440e-1, 1.010469450},
                                   {6.000698670e-3, 2.001791440e-2, 1.035606530e2}};

/// Small catalogue with every record type of the format description. GLASS-A is invented. B270
/// is copied (NM and CD only) from the free SCHOTT catalogue "schott glasses preferred and
/// special June-2025-B.AGF" (schott.com, downloaded 2026-10-04; docs/quellen.md, #24).
constexpr std::string_view kText =
    "CC test catalogue\n"
    "! optional comment line\n"
    "NM GLASS-A 2 500600.000 1.5 60.0 0 1\n"
    "GC some comment with spaces\n"
    "ED 7.1 8.3 2.51 -0.0009 0\n"
    "CD 1.0 0.01 0.2 0.02 1.0 100.0 0 0\n"
    "TD 1e-6 1e-8 -1e-11 4e-7 6e-10 0.17 20\n"
    "MD 82.00 0.21 610 858.000 1.11\n"
    "OD 1.0 1.0 0.0 1.0 2.3 2.3\n"
    "LD 0.3 2.5\n"
    "IT 0.3 0.05 25\n"
    "IT 0.31 0.25 25\n"
    "BD 0.5 3.0 -0.5 -3.5\n"
    "\n"
    "NM B270 1 1 1.523080 58.571369 0 3 0\n"
    "CD 2.286575000E+000 -8.733458200E-003 1.174288400E-002 2.904175600E-004 -1.250669500E-005 "
    "9.264625300E-007 0.000000000E+000 0.000000000E+000 0.000000000E+000 0.000000000E+000\n";

/// Writes a temporary catalogue file and removes it at the end of the scope.
class TempCatalog {
 public:
  TempCatalog(const std::string& directory, const std::string& file, const std::string& text)
      : dir_(fs::temp_directory_path() /
             (directory + "_" +
              std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()))) {
    fs::remove_all(dir_);
    fs::create_directories(dir_);
    std::ofstream(dir_ / file, std::ios::binary) << text;
  }
  TempCatalog(const TempCatalog&) = delete;
  TempCatalog& operator=(const TempCatalog&) = delete;
  ~TempCatalog() {
    std::error_code ec;
    fs::remove_all(dir_, ec);
  }
  [[nodiscard]] const fs::path& dir() const { return dir_; }

 private:
  fs::path dir_;
};

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

}  // namespace

// ------------------------------------------------------------------- encoding -----

TEST_CASE("AGF text: UTF-16LE with BOM, UTF-8 with BOM and ANSI decode to the same text", "[agf]") {
  const std::string text = "CC caf\xC3\xA9\r\nNM A 2 1 1.5 60\r\n";  // "café" in UTF-8
  // UTF-16LE with BOM: 'C' 'C' ' ' 'c' 'a' 'f' U+00E9 ...
  std::string utf16 = "\xFF\xFE";
  for (const char32_t c : std::u32string(U"CC café\r\nNM A 2 1 1.5 60\r\n")) {
    utf16 += static_cast<char>(c & 0xFF);
    utf16 += static_cast<char>((c >> 8) & 0xFF);
  }
  REQUIRE(decode_agf_text(utf16, "f.agf") == text);
  REQUIRE(decode_agf_text("\xEF\xBB\xBF" + text, "f.agf") == text);
  REQUIRE(decode_agf_text(text, "f.agf") == text);
}

TEST_CASE("AGF text: UTF-16 surrogate pairs and errors", "[agf]") {
  // U+1F600 = D83D DE00 -> F0 9F 98 80 in UTF-8.
  const std::string pair = std::string("\xFF\xFE", 2) + std::string("\x3D\xD8\x00\xDE", 4);
  REQUIRE(decode_agf_text(pair, "f.agf") == "\xF0\x9F\x98\x80");
  // Odd number of bytes, lone surrogate, UTF-16BE.
  REQUIRE_THROWS_AS(decode_agf_text(std::string("\xFF\xFE\x41", 3), "f.agf"), AgfError);
  REQUIRE_THROWS_AS(decode_agf_text(std::string("\xFF\xFE\x00\xDE", 4), "f.agf"), AgfError);
  REQUIRE_THROWS_AS(decode_agf_text(std::string("\xFE\xFF\x00\x41", 4), "f.agf"), AgfError);
}

// -------------------------------------------------------------------- parsing -----

TEST_CASE("AGF records NM, CD, LD, TD and ED are read; others are skipped", "[agf]") {
  for (const bool crlf : {false, true}) {
    std::string text(kText);
    if (crlf) {
      std::string converted;
      for (const char c : text) converted += c == '\n' ? std::string("\r\n") : std::string(1, c);
      text = converted;
    }
    INFO("CRLF: " << crlf);
    const AgfCatalog cat = parse_agf(text, "TEST", "test.agf");
    REQUIRE(cat.name == "TEST");
    REQUIRE(cat.comment == "test catalogue");
    REQUIRE(cat.glasses.size() == 2);

    const AgfGlass& a = cat.glasses[0];
    REQUIRE(a.name == "GLASS-A");
    REQUIRE(a.formula == 2);
    REQUIRE(a.nd == 1.5);
    REQUIRE(a.vd == 60.0);
    REQUIRE(a.line == 3);
    REQUIRE(a.coefficients == std::vector<double>{1.0, 0.01, 0.2, 0.02, 1.0, 100.0, 0.0, 0.0});
    REQUIRE(a.range == WavelengthRange{0.3, 2.5});
    REQUIRE(a.thermal == std::vector<double>{1e-6, 1e-8, -1e-11, 4e-7, 6e-10, 0.17, 20.0});
    REQUIRE(a.extra == std::vector<double>{7.1, 8.3, 2.51, -0.0009, 0.0});

    const AgfGlass& b = cat.glasses[1];
    REQUIRE(b.name == "B270");
    REQUIRE(b.formula == 1);
    REQUIRE(b.line == 15);
    REQUIRE(b.coefficients.size() == 10);
    REQUIRE_FALSE(b.range.has_value());
    REQUIRE_FALSE(b.thermal.has_value());
  }
}

TEST_CASE("malformed AGF lines are errors with file and line", "[agf]") {
  struct Case {
    std::string text;
    std::size_t line;
  };
  const std::vector<Case> cases{
      {"CC c\nCD 1 2 3\n", 2},                                         // data before the first NM
      {"CC c\nNM A 2 1 1.5 60\nXX 1\nCD 1 2 3 4 5 6\n", 3},            // unknown mnemonic
      {"CC c\nNM A 2 1 1.5 60\nCD 1 x 3 4 5 6\n", 3},                  // not a number
      {"CC c\nNM A 2 1 1.5\nCD 1 2 3 4 5 6\n", 2},                     // NM without V(d)
      {"CC c\nNM A 2.5 1 1.5 60\nCD 1 2 3 4 5 6\n", 2},                // formula not an integer
      {"CC c\nNM A 2 1 1.5 60\nCD 1 2 3 4 5 6 7 8 9 10 11\n", 3},      // more than 10 coefficients
      {"CC c\nNM A 2 1 1.5 60\nCD 1 2 3 4 5 6\nLD 0.3\n", 4},          // LD with one value
      {"CC c\nNM A 2 1 1.5 60\nCD 1 2 3 4 5 6\nLD 2.5 0.3\n", 4},      // LD min >= max
      {"CC c\nNM A 2 1 1.5 60\nCD 1 2 3 4 5 6\nCD 1 2 3 4 5 6\n", 4},  // second CD
      {"CC c\nNM A 2 1 1.5 60\nLD 0.3 2.5\n", 2},                      // glass without CD
      {"CC c\nNM A 2 1 1.5 60\nCD 1 2 3 4 5 6\nNM A 2 1 1.5 60\nCD 1 2 3 4 5 6\n", 4},  // twice
      {"CC c\nNM A 2 1 1.5 60\nCD 1 2 3 4 5 6\nNM\n", 4},               // NM without name
      {"CC c\nNM A 2 1 1.5 60\nCD\n", 3},                               // CD without coefficients
      {"CC c\nNM A 2 1 1.5 60\nCD 1 2\nLD 0.3 2.5\nLD 0.4 2.0\n", 5},   // second LD
      {"CC c\nNM A 2 1 1.5 60\nCD 1 2\nTD 1 2 3 4 5 6 20\nTD 1\n", 5},  // second TD
      {"CC c\nNM A 2 1 1.5 60\nCD 1 2\nED 1 2 3 4 0\nED 1\n", 5},       // second ED
  };
  for (const auto& c : cases) {
    INFO(c.text);
    const AgfError e = agf_error([&] { (void)parse_agf(c.text, "TEST", "bad.agf"); });
    REQUIRE(e.file() == "bad.agf");
    REQUIRE(e.line() == c.line);
    REQUIRE_THAT(e.what(), ContainsSubstring("bad.agf:" + std::to_string(c.line)));
  }
}

// ------------------------------------------------------------ formula mapping -----

TEST_CASE("AGF formula 1 is Schott a0..a5, formula 2 is Sellmeier 1 K1 L1 K2 L2 K3 L3", "[agf]") {
  // Order verified against N(d) and V(d) of the free SCHOTT catalogue (docs/quellen.md).
  AgfGlass g;
  g.name = "G";
  g.formula = 1;
  g.coefficients = {1, 2, 3, 4, 5, 6, 0, 0, 0, 0};
  REQUIRE(rtt::material::agf_formula(g, "w") ==
          rtt::material::DispersionFormula{SchottCoefficients{{1, 2, 3, 4, 5, 6}}});
  g.formula = 2;
  REQUIRE(rtt::material::agf_formula(g, "w") ==
          rtt::material::DispersionFormula{Sellmeier1Coefficients{{1, 3, 5}, {2, 4, 6}}});
  // Missing trailing coefficients count as 0 (catalogues often write fewer than 10).
  g.coefficients = {1, 2, 3, 4};
  REQUIRE(rtt::material::agf_formula(g, "w") ==
          rtt::material::DispersionFormula{Sellmeier1Coefficients{{1, 3, 0}, {2, 4, 0}}});
}

TEST_CASE("CD order of formulas 1 and 2 reproduces N(d) and V(d) of the NM record", "[agf]") {
  // Reference for the coefficient order (docs/quellen.md, #24): with the order a0..a5
  // (formula 1, B270) and K1 L1 K2 L2 K3 L3 (formula 2, N-BK7 and F2), the catalogue's own N(d)
  // and V(d) come out. Abbe number V(d) = (n_d - 1) / (n_F - n_C) with n_d at 587.56 nm, n_F at
  // 486.13 nm and n_C at 656.27 nm (SCHOTT TIE-29, April 2005, p. 1 and Eq. (2.1-1)).
  // Tolerances: half a unit of the last digit in the NM record (n_d: 5e-6 covers 1.5168;
  // V(d): 0.005). A wrong order is far off (B270 with a1/a2 swapped: V(d) = -96; N-BK7 with
  // K1 K2 K3 L1 L2 L3: V(d) = 26.8).
  constexpr double kD = 0.58756;
  constexpr double kF = 0.48613;
  constexpr double kC = 0.65627;
  const auto check = [&](const AgfGlass& glass, const std::string& catalog) {
    INFO(glass.name);
    const CatalogMaterial m(glass, catalog);
    // The formula itself (relative index over relative wavelength, as in the NM record); since
    // #25 index() is absolute and takes the vacuum wavelength, see the test below.
    const auto n = [&](double wl) { return rtt::material::refractive_index(m.formula(), wl); };
    REQUIRE(std::abs(n(kD) - glass.nd) <= 5e-6);
    REQUIRE(std::abs((n(kD) - 1.0) / (n(kF) - n(kC)) - glass.vd) <= 0.005);
  };
  const AgfCatalog text = parse_agf(kText, "TEST", "test.agf");
  REQUIRE(text.glasses[1].formula == 1);
  check(text.glasses[1], "TEST");  // B270
  const AgfCatalog schott = load_agf(kSchottFile);
  for (const AgfGlass& glass : schott.glasses) {
    REQUIRE(glass.formula == 2);
    check(glass, schott.name);  // N-BK7, F2
  }
}

TEST_CASE("unsupported AGF formula numbers are a clear error", "[agf]") {
  AgfGlass g;
  g.name = "IRG";
  g.coefficients = {1, 2, 3, 4, 5, 6};
  for (const int formula : {0, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14}) {
    INFO(formula);
    g.formula = formula;
    try {
      (void)rtt::material::agf_formula(g, "cat.agf:7");
      FAIL("no exception");
    } catch (const std::invalid_argument& e) {
      REQUIRE_THAT(e.what(), ContainsSubstring("formula " + std::to_string(formula)));
      REQUIRE_THAT(e.what(), ContainsSubstring("IRG"));
      REQUIRE_THAT(e.what(), ContainsSubstring("cat.agf:7"));
      REQUIRE_THAT(e.what(), ContainsSubstring("#42"));
    }
  }
}

// ------------------------------------------------------- files and materials -----

TEST_CASE("test catalogue: name from the file, glasses N-BK7 and F2", "[agf]") {
  const AgfCatalog cat = load_agf(kSchottFile);
  REQUIRE(cat.name == "SCHOTT");
  REQUIRE(cat.glasses.size() == 2);
  REQUIRE(cat.glasses[0].name == "N-BK7");
  REQUIRE(cat.glasses[1].name == "F2");
  REQUIRE(cat.glasses[0].range == WavelengthRange{0.3, 2.5});
  REQUIRE(cat.glasses[1].range == WavelengthRange{0.32, 2.5});
}

TEST_CASE("UTF-16 and ANSI catalogues with the same content give the same glasses", "[agf]") {
  const AgfCatalog ansi = load_agf(kSchottFile);
  const AgfCatalog utf16 = load_agf(kSchottUtf16);
  REQUIRE(utf16.name == "SCHOTT");
  REQUIRE(utf16.comment == ansi.comment);
  REQUIRE(utf16.glasses.size() == ansi.glasses.size());
  for (std::size_t i = 0; i < ansi.glasses.size(); ++i) {
    const AgfGlass& a = ansi.glasses[i];
    const AgfGlass& u = utf16.glasses[i];
    INFO(a.name);
    REQUIRE(u.name == a.name);
    REQUIRE(u.formula == a.formula);
    REQUIRE(u.coefficients == a.coefficients);
    REQUIRE(u.range == a.range);
    REQUIRE(u.thermal == a.thermal);
    REQUIRE(u.extra == a.extra);
  }
}

TEST_CASE("catalogue glass: absolute index consistent with the relative catalogue data",
          "[agf][air]") {
  // #25: the catalogue formula gives the index relative to air at T_ref and 1 atm as a function
  // of the wavelength in that air. So at T = T_ref: n_abs / n_air(lambda_vac, T_ref, 1 atm) =
  // n_rel(lambda_vac / n_air) with the formula of #23 (relative 1e-14), for any pressure of the
  // medium (a solid ignores it).
  const AgfCatalog cat = load_agf(kSchottFile);
  const CatalogMaterial nbk7(cat.glasses[0], cat.name);
  REQUIRE(nbk7.formula() == rtt::material::DispersionFormula{kNbk7});
  const double t_ref = 20.0;  // TD record of N-BK7
  for (const double wl : {0.3, 0.4861, 0.5876, 0.6563, 1.064, 2.5}) {
    INFO(wl);
    const double n_air = rtt::material::ciddor_air_index(wl, t_ref, 101325.0);
    const double n_rel = rtt::material::refractive_index(kNbk7, wl / n_air);
    for (const double p : {0.0, 1.0, 2.0}) {
      const auto n = nbk7.index(wl, t_ref, p);
      REQUIRE_THAT(n.real() / n_air, Catch::Matchers::WithinRel(n_rel, 1e-14));
      REQUIRE(n.imag() == 0.0);
    }
  }
  // n_d of the SCHOTT data sheet (1.51680 at the d line, 587.56 nm in air; SCHOTT TIE-29 p. 1):
  // the vacuum wavelength is the air wavelength times n_air.
  const double lambda_d = 0.58756 * rtt::material::ciddor_air_index(0.58756, t_ref, 101325.0);
  const double n_d = nbk7.index(lambda_d, t_ref, 1.0).real() /
                     rtt::material::ciddor_air_index(lambda_d, t_ref, 101325.0);
  REQUIRE(std::abs(n_d - 1.51680) <= 5e-6);
  REQUIRE(nbk7.wavelength_range_um() == WavelengthRange{0.3, 2.5});
  REQUIRE(nbk7.glass().name == "N-BK7");
  REQUIRE(nbk7.glass().thermal.has_value());  // stored for #25

  // The constructor checks the LD range of glasses built by hand.
  AgfGlass bad = cat.glasses[0];
  bad.range = WavelengthRange{2.5, 0.3};
  REQUIRE_THROWS_AS(CatalogMaterial(bad, "SCHOTT"), std::invalid_argument);
  bad.range = WavelengthRange{0.0, 2.5};
  REQUIRE_THROWS_AS(CatalogMaterial(bad, "SCHOTT"), std::invalid_argument);
}

TEST_CASE("catalogue glass: dn/dT matches the SCHOTT data sheet", "[agf][air]") {
  // SCHOTT N-BK7 data sheet (as of 01-Dec-2023, p. 13): "Temperature Coefficients of the
  // Refractive Index", Delta n_abs / Delta T in 1e-6/K for +20/+40 degC: 1.1 at 1060.0 nm,
  // 1.6 at the e line (546.1 nm), 2.1 at the g line (435.8 nm); rounded to 0.1e-6/K, so the
  // tolerance is 0.05e-6/K. Model: TIE-19 Eq. (3)/(4) with the TD record of the catalogue.
  const AgfCatalog cat = load_agf(kSchottFile);
  const CatalogMaterial nbk7(cat.glasses[0], cat.name);
  for (const auto& [wl, sheet] : {std::pair{1.0600, 1.1}, {0.5461, 1.6}, {0.4358, 2.1}}) {
    INFO(wl);
    const double slope =
        (nbk7.index(wl, 40.0, 1.0).real() - nbk7.index(wl, 20.0, 1.0).real()) / 20.0;
    REQUIRE(std::abs(slope * 1e6 - sheet) <= 0.05);
  }
}

TEST_CASE("catalogue glass without TD record has no temperature dependence", "[agf][air]") {
  // Without the six coefficients there is no thermal model (#25, P4): T_ref = 20 degC.
  AgfGlass glass = load_agf(kSchottFile).glasses[0];
  glass.thermal.reset();
  const CatalogMaterial m(glass, "SCHOTT");
  REQUIRE(m.index(0.5876, -20.0, 1.0) == m.index(0.5876, 60.0, 1.0));
  const double n_air = rtt::material::ciddor_air_index(0.5876, 20.0, 101325.0);
  REQUIRE_THAT(
      m.index(0.5876, 60.0, 1.0).real() / n_air,
      Catch::Matchers::WithinRel(rtt::material::refractive_index(kNbk7, 0.5876 / n_air), 1e-14));
  // A TD record must have the seven values D0 D1 D2 E0 E1 Ltk Temp.
  glass.thermal = std::vector<double>{1e-6, 1e-8};
  REQUIRE_THROWS_AS(CatalogMaterial(glass, "SCHOTT"), std::invalid_argument);
}

TEST_CASE("add_catalog: resolve KATALOG:NAME, same object, thread-safe", "[agf]") {
  MaterialLibrary lib;
  lib.add_catalog(kSchottFile);
  const auto nbk7 = lib.resolve("SCHOTT:N-BK7");
  REQUIRE(lib.resolve("SCHOTT:N-BK7") == nbk7);
  const auto* catalog_glass = dynamic_cast<const CatalogMaterial*>(nbk7.get());
  REQUIRE(catalog_glass != nullptr);
  REQUIRE(catalog_glass->formula() == rtt::material::DispersionFormula{kNbk7});
  REQUIRE(lib.resolve("SCHOTT:F2") != nbk7);

  std::vector<const rtt::material::Material*> seen(8, nullptr);
  {
    std::vector<std::jthread> threads;
    for (std::size_t t = 0; t < seen.size(); ++t) {
      threads.emplace_back([&lib, &seen, t] { seen[t] = lib.resolve("SCHOTT:F2").get(); });
    }
  }
  for (const auto* m : seen) REQUIRE(m == lib.resolve("SCHOTT:F2").get());
}

TEST_CASE("add_catalog loads the .agf files of a directory, not recursively", "[agf]") {
  MaterialLibrary lib;
  lib.add_catalog(kCatalogDir);  // schott.agf; utf16/schott.agf is in a subdirectory
  REQUIRE(lib.resolve("SCHOTT:N-BK7") != nullptr);
  MaterialLibrary lib16;
  lib16.add_catalog(kCatalogDir / "utf16");
  REQUIRE(lib16.resolve("SCHOTT:F2")->index(0.5876, 20.0, 1.0) ==
          lib.resolve("SCHOTT:F2")->index(0.5876, 20.0, 1.0));
}

TEST_CASE("unknown glasses and catalogues name the catalogue", "[agf]") {
  MaterialLibrary lib;
  lib.add_catalog(kSchottFile);
  try {
    (void)lib.resolve("SCHOTT:N-SF6");
    FAIL("no exception");
  } catch (const UnknownMaterial& e) {
    REQUIRE_THAT(e.what(), ContainsSubstring("N-SF6"));
    REQUIRE_THAT(e.what(), ContainsSubstring("catalog SCHOTT"));
    REQUIRE_THAT(e.what(), ContainsSubstring("schott.agf"));
  }
  try {
    (void)lib.resolve("OHARA:S-BSL7");
    FAIL("no exception");
  } catch (const UnknownMaterial& e) {
    REQUIRE_THAT(e.what(), ContainsSubstring("catalog OHARA"));
    REQUIRE_THAT(e.what(), ContainsSubstring("not loaded"));
  }
}

TEST_CASE("glass with an unsupported formula fails on resolve with file and line", "[agf]") {
  const TempCatalog tmp("rtt_agf_unsupported", "ir.AGF",
                        "CC c\nNM GOOD 2 1 1.5 60\nCD 1 0.01 0 0 0 0\n"
                        "NM IRG 3 1 2.5 20\nCD 1 2 3 4 5 6\n");
  MaterialLibrary lib;
  lib.add_catalog(tmp.dir());  // upper-case extension counts
  REQUIRE(lib.resolve("IR:GOOD") != nullptr);
  try {
    (void)lib.resolve("IR:IRG");
    FAIL("no exception");
  } catch (const UnknownMaterial& e) {
    REQUIRE_THAT(e.what(), ContainsSubstring("formula 3"));
    REQUIRE_THAT(e.what(), ContainsSubstring("IRG"));
    REQUIRE_THAT(e.what(), ContainsSubstring("ir.AGF:4"));
    REQUIRE_THAT(e.what(), ContainsSubstring("#42"));
  }
}

TEST_CASE("add_catalog rejects conflicts and registers nothing then", "[agf]") {
  MaterialLibrary lib;
  lib.add_catalog(kSchottFile);
  REQUIRE_THROWS_AS(lib.add_catalog(kSchottFile), std::invalid_argument);  // SCHOTT loaded
  REQUIRE_THROWS_AS(lib.add_catalog(kCatalogDir / "utf16"), std::invalid_argument);

  MaterialLibrary other;
  const auto mine = std::make_shared<const DispersionMaterial>(kNbk7, WavelengthRange{0.3, 2.5});
  other.add("SCHOTT:N-BK7", mine);
  REQUIRE_THROWS_AS(other.add_catalog(kSchottFile), std::invalid_argument);
  REQUIRE(other.resolve("SCHOTT:N-BK7") == mine);
  REQUIRE_THROWS_AS(other.resolve("SCHOTT:F2"), UnknownMaterial);  // nothing registered

  const TempCatalog reserved("rtt_agf_reserved", "const.agf", "CC c\nNM A 2 1 1.5 60\nCD 1\n");
  REQUIRE_THROWS_AS(lib.add_catalog(reserved.dir() / "const.agf"), std::invalid_argument);
  const TempCatalog empty("rtt_agf_empty", "readme.txt", "no catalogue here\n");
  REQUIRE_THROWS_AS(lib.add_catalog(empty.dir()), std::invalid_argument);
  REQUIRE_THROWS_AS(lib.add_catalog(empty.dir() / "missing.agf"), AgfError);
}
