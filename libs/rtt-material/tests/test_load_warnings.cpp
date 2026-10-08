// Narrow, documented exceptions to the strict AGF reading (ADR 0008 addendum, #71) and the
// LoadWarning channel of rtt-material (ADR 0022). The catalogues are excerpts: single lines
// quoted from the manufacturers' files (docs/quellen.md), the rest invented.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cstddef>
#include <set>
#include <string>
#include <vector>

#include "rtt/diagnostics/codes.hpp"
#include "rtt/material/agf.hpp"
#include "rtt/material/load_warning.hpp"
#include "rtt/material/material.hpp"

using Catch::Matchers::ContainsSubstring;
using rtt::material::AgfCatalog;
using rtt::material::AgfError;
using rtt::material::LoadWarning;
using rtt::material::MaterialLibrary;
using rtt::material::parse_agf;
using rtt::material::UnknownMaterial;

namespace {

AgfError agf_error(const std::string& text) {
  try {
    (void)parse_agf(text, "TEST", "bad.agf");
  } catch (const AgfError& e) {
    return e;
  }
  FAIL("no AgfError for\n" << text);
  return AgfError("", 0, "");
}

bool has_warning(const std::vector<LoadWarning>& w, const std::string& code, std::size_t line) {
  for (const LoadWarning& x : w) {
    if (x.code == code && x.line == line) return true;
  }
  return false;
}

/// A complete glass block (Schott formula, valid coefficients).
std::string block(const std::string& name,
                  const std::string& vd = "64.17",
                  const std::string& ed = "") {
  return "NM " + name + " 2 517642 1.5168 " + vd +
         " 0 1 1\nCD 1.03961212 0.00600069867 0.231792344 0.0200179144 1.01046945 103.560653\n" +
         (ed.empty() ? "" : "ED " + ed + "\n") + "LD 0.3 2.5\n";
}

}  // namespace

TEST_CASE("R1: identical duplicate glass blocks are merged with a warning", "[agf][load]") {
  // NIKON-HIKARI_201911.AGF repeats E-BAK1, E-BAK2 and E-BAK4 with identical blocks
  // (lines 7228/7300, 7259/7331, 7372/7410). Line ends and blanks do not matter.
  const std::string first = block("E-BAK1");
  std::string second = block("E-BAK1");
  for (std::size_t i = second.find(' '); i != std::string::npos; i = second.find(' ', i + 3)) {
    second.replace(i, 1, "  ");  // more blanks, CRLF line ends
  }
  for (std::size_t i = second.find('\n'); i != std::string::npos; i = second.find('\n', i + 2)) {
    second.replace(i, 1, "\r\n");
  }
  const AgfCatalog cat = parse_agf("CC excerpt\n" + first + block("B") + second, "N", "n.agf");
  REQUIRE(cat.glasses.size() == 2);  // E-BAK1 once, B
  REQUIRE(cat.glasses[0].name == "E-BAK1");
  REQUIRE(cat.glasses[0].line == 2);
  REQUIRE(cat.warnings.size() == 1);
  const LoadWarning& w = cat.warnings[0];
  REQUIRE(w.code == "agf.duplicate_glass");
  REQUIRE(w.file == "n.agf");
  REQUIRE(w.line == 8);  // the repetition
  REQUIRE_THAT(w.message, ContainsSubstring("E-BAK1") && ContainsSubstring("line 2"));

  MaterialLibrary lib;
  lib.add_catalog_text("CC excerpt\n" + first + block("B") + second, "N", "n.agf");
  REQUIRE(lib.resolve("N:E-BAK1") != nullptr);
  REQUIRE(lib.load_warnings() == cat.warnings);
}

TEST_CASE("the former error case 'glass appears twice' under R1 and R2", "[agf][load]") {
  // Moved from the error table of test_agf.cpp (approved with #71): two identical blocks are
  // one glass and a warning at line 4; with another V(d) they are a conflict at line 4.
  const AgfCatalog same = parse_agf(
      "CC c\nNM A 2 1 1.5 60\nCD 1 2 3 4 5 6\nNM A 2 1 1.5 60\nCD 1 2 3 4 5 6\n", "T", "t.agf");
  REQUIRE(same.glasses.size() == 1);
  REQUIRE(same.warnings.size() == 1);
  REQUIRE(same.warnings[0].code == "agf.duplicate_glass");
  REQUIRE(same.warnings[0].line == 4);
  const AgfCatalog other = parse_agf(
      "CC c\nNM A 2 1 1.5 60\nCD 1 2 3 4 5 6\nNM A 2 1 1.5 61\nCD 1 2 3 4 5 6\n", "T", "t.agf");
  REQUIRE(other.glasses.size() == 2);
  REQUIRE(other.warnings.size() == 1);
  REQUIRE(other.warnings[0].code == "agf.duplicate_glass_conflict");
  REQUIRE(other.warnings[0].line == 4);
}

TEST_CASE("R2: differing duplicates make the glass ambiguous", "[agf][load]") {
  // NIKON E-F2 (lines 8160/8191): V(d) 36.258938 against 36.258932, further digits in CD.
  const std::string text =
      "CC excerpt\n" + block("E-F2", "36.258938") + block("B") + block("E-F2", "36.258932");
  const AgfCatalog cat = parse_agf(text, "N", "n.agf");
  REQUIRE(cat.glasses.size() == 3);  // both blocks stay in the listing
  REQUIRE(cat.warnings.size() == 1);
  REQUIRE(cat.warnings[0].code == "agf.duplicate_glass_conflict");
  REQUIRE(cat.warnings[0].line == 8);
  REQUIRE_THAT(cat.warnings[0].message, ContainsSubstring("line 2") && ContainsSubstring("NM"));

  MaterialLibrary lib;
  lib.add_catalog_text(text, "N", "n.agf");
  REQUIRE(lib.resolve("N:B") != nullptr);  // the rest of the catalogue loads
  try {
    (void)lib.resolve("N:E-F2");
    FAIL("no UnknownMaterial");
  } catch (const UnknownMaterial& e) {
    const std::string message = e.what();
    REQUIRE_THAT(message, ContainsSubstring("ambiguous"));
    REQUIRE_THAT(message, ContainsSubstring("lines 2 and 8"));
    REQUIRE_THAT(message, ContainsSubstring("alias"));  // how to use it anyway
  }

  SECTION("complementary data (old HOYA FD1: ED only in one block) is a conflict as well") {
    const AgfCatalog fd = parse_agf("CC c\n" + block("FD1", "29.500280", "7.9 89 4.46 0.0034 0") +
                                        block("FD1", "29.500288", "0 0 0 0 0"),
                                    "H", "h.agf");
    REQUIRE(fd.warnings.size() == 1);
    REQUIRE(fd.warnings[0].code == "agf.duplicate_glass_conflict");
    REQUIRE_THAT(fd.warnings[0].message, ContainsSubstring("ED"));
  }
}

TEST_CASE("warnings are in file order, duplicates included", "[agf][load]") {
  // Duplicates are found after the last line; the warnings still follow the lines.
  const AgfCatalog cat =
      parse_agf("CC c\n" + block("A") + block("A") + "E\n" + block("B"), "T", "t.agf");
  REQUIRE(cat.warnings.size() == 2);
  REQUIRE((cat.warnings[0].code == "agf.duplicate_glass" && cat.warnings[0].line == 5));
  REQUIRE((cat.warnings[1].code == "agf.stray_line" && cat.warnings[1].line == 8));
}

TEST_CASE("more than two blocks of one name", "[agf][load]") {
  // Lines 2, 5, 8, 11: A, A with another V(d), a copy of the first A, A with a third V(d).
  const std::string text = "CC c\n" + block("A", "60") + block("A", "61") + block("A", "60") +
                           block("A", "62") + block("B");
  const AgfCatalog cat = parse_agf(text, "T", "t.agf");
  REQUIRE(cat.glasses.size() == 4);  // A at 2, 5 and 11; B
  REQUIRE(cat.warnings.size() == 3);
  REQUIRE((cat.warnings[0].code == "agf.duplicate_glass_conflict" && cat.warnings[0].line == 5));
  REQUIRE((cat.warnings[1].code == "agf.duplicate_glass" && cat.warnings[1].line == 8));
  REQUIRE((cat.warnings[2].code == "agf.duplicate_glass_conflict" && cat.warnings[2].line == 11));
  MaterialLibrary lib;
  lib.add_catalog_text(text, "T", "t.agf");
  REQUIRE_THROWS_WITH(lib.resolve("T:A"), ContainsSubstring("lines 2, 5 and 11"));
  REQUIRE(lib.resolve("T:B") != nullptr);
}

TEST_CASE("R3: lines before the first record are skipped with a warning", "[agf][load]") {
  // Line 1 of nine files of the ZemaxGlass distribution (lightpath.agf, schott.agf, ...).
  const std::string header =
      "Reproduced here by permission of RadiantZemax (www.radiantzemax.com).";
  const AgfCatalog cat =
      parse_agf(header + "\n\n" + header + "\nCC c\n" + block("A"), "Z", "z.agf");
  REQUIRE(cat.glasses.size() == 1);
  REQUIRE(cat.warnings.size() == 2);
  REQUIRE(has_warning(cat.warnings, "agf.preamble_skipped", 1));
  REQUIRE(has_warning(cat.warnings, "agf.preamble_skipped", 3));
  // Only before the first record: after CC (or NM) the strict reading holds.
  REQUIRE(agf_error("CC c\n" + header + "\n" + block("A")).line() == 2);
  REQUIRE(agf_error(block("A") + header + "\n").line() == 4);
}

TEST_CASE("R3: a file without any CC or NM record is not an AGF catalogue", "[agf][load]") {
  // R3 only skips lines before a first record. Prose and a renamed JSON file have none and stay
  // an error at the first skipped line, as before #71.
  const AgfError prose = agf_error("This is not a glass catalogue.\nJust prose.\n");
  REQUIRE(prose.line() == 1);
  REQUIRE_THAT(std::string(prose.what()), ContainsSubstring("no CC or NM record"));
  const std::string json = "\n{\n  \"name\": \"N-BK7\",\n  \"nd\": \"1.5168\"\n}\n";
  REQUIRE(agf_error(json).line() == 2);
  MaterialLibrary lib;
  REQUIRE_THROWS_AS(lib.add_catalog_text(json, "J", "j.agf"), AgfError);
  REQUIRE(lib.load_warnings().empty());
  // Unchanged: a file without any text line is an empty catalogue.
  const AgfCatalog empty = parse_agf("", "T", "t.agf");
  REQUIRE((empty.glasses.empty() && empty.warnings.empty()));
  const AgfCatalog comments = parse_agf("! comment\n\n", "T", "t.agf");
  REQUIRE((comments.glasses.empty() && comments.warnings.empty()));
}

TEST_CASE("R4: a single word before the next NM is skipped with a warning", "[agf][load]") {
  // NIKON-HIKARI_201911.AGF lines 10101-10104 verbatim: empty line, "E", empty line,
  // "NM E-LAKH1"; only empty lines may come between the word and the NM.
  const AgfCatalog cat =
      parse_agf("CC c\n" + block("A") + "\r\nE\r\n\r\n" + block("E-LAKH1"), "N", "n.agf");
  REQUIRE(cat.glasses.size() == 2);
  REQUIRE(cat.glasses[1].name == "E-LAKH1");
  REQUIRE(cat.warnings.size() == 1);
  REQUIRE(cat.warnings[0].code == "agf.stray_line");
  REQUIRE(cat.warnings[0].line == 6);
  SECTION("anything else stays an error at the stray line") {
    REQUIRE(agf_error("CC c\n" + block("A") + "E X\n" + block("B")).line() == 5);     // two words
    REQUIRE(agf_error("CC c\n" + block("A") + "E1\n" + block("B")).line() == 5);      // not a word
    REQUIRE(agf_error("CC c\n" + block("A") + "E\nIT 0.3 0.9 10\n").line() == 5);     // no NM
    REQUIRE(agf_error("CC c\n" + block("A") + "E\n! c\n" + block("B")).line() == 5);  // comment
    REQUIRE(agf_error("CC c\n" + block("A") + "E\n\nIT 0.3 0.9 10\n").line() == 5);   // no NM
    REQUIRE(agf_error("CC c\n" + block("A") + "E\n").line() == 5);                    // at the end
  }
}

TEST_CASE("R5: data errors of old catalogue copies stay errors", "[agf][load]") {
  // Old SCHOTT of the ZemaxGlass distribution, line 1609: IT without transmittance.
  REQUIRE(agf_error("CC c\n" + block("A") + "IT 2.90000E-01  2.50000E+01      \n").line() == 5);
  // zeon.agf line 32: resistance class -11.000.
  REQUIRE(agf_error("CC c\n" + block("A") + "OD  -1.000 -1.000 -1.000 -1.000 -1.000 -11.000\n")
              .line() == 5);
}

TEST_CASE("a failed catalogue adds no warnings", "[agf][load]") {
  MaterialLibrary lib;
  const std::string header =
      "Reproduced here by permission of RadiantZemax (www.radiantzemax.com).";
  REQUIRE_THROWS_AS(lib.add_catalog_text(header + "\nCC c\n" + block("A") + "XX 1\n", "Z", "z.agf"),
                    AgfError);
  REQUIRE(lib.load_warnings().empty());
  lib.add_catalog_text(header + "\nCC c\n" + block("A"), "Z", "z.agf");
  lib.add_catalog_text("CC c\n" + block("A") + block("A"), "Y", "y.agf");
  // In loading order.
  const std::vector<LoadWarning> w = lib.load_warnings();
  REQUIRE(w.size() == 2);
  REQUIRE((w[0].file == "z.agf" && w[0].code == "agf.preamble_skipped"));
  REQUIRE((w[1].file == "y.agf" && w[1].code == "agf.duplicate_glass"));
}

TEST_CASE("every agf code of the registry has a case in this file", "[agf][load]") {
  const std::set<std::string> covered = {"agf.duplicate_glass", "agf.duplicate_glass_conflict",
                                         "agf.preamble_skipped", "agf.stray_line"};
  for (const rtt::diagnostics::CodeInfo& info : rtt::diagnostics::kCodes) {
    if (info.producer != "agf") continue;
    INFO("no case for " << info.code);
    REQUIRE(covered.contains(std::string(info.code)));
    REQUIRE(info.severity == rtt::diagnostics::Severity::Warning);
  }
  for (const std::string& code : covered) REQUIRE(rtt::diagnostics::find_code(code) != nullptr);
}
