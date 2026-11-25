// ZMX import of several configurations (#361; ADR 0044, addendum #361). The importer builds the
// system once per configuration with the operand values of that configuration and binds every
// Param that differs to a row of the parameter table (ADR 0029). Acceptance: the system of each
// configuration k, resolve_parameters(import, k), equals bitwise the import of a single
// configuration file with the values of k. The files are the cases of tests/zmx with single
// lines changed and a configuration block appended.

#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <fstream>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

#include "rtt/analysis/spot.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/zmx.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/parameters.hpp"
#include "rtt/model/system.hpp"

namespace {

using rtt::io::import_zmx;
using rtt::io::ZmxError;
using rtt::io::ZmxImport;
using rtt::model::System;

const std::string kZmxDir = std::string(RTT_ZMX_DIR) + "/";

std::string read_file(const std::string& path) {
  const std::ifstream in(path, std::ios::binary);
  REQUIRE(in);
  std::ostringstream buffer;
  buffer << in.rdbuf();
  return buffer.str();
}

const rtt::material::MaterialLibrary& schott() {
  static rtt::material::MaterialLibrary lib;
  static const bool loaded = [] {
    lib.add_catalog(std::string(RTT_CATALOG_DIR) + "/schott.agf");
    return true;
  }();
  (void)loaded;
  return lib;
}

/// `text` with the first line starting with `prefix` at or after `from` replaced by
/// `replacement` (without its line end).
std::string replace_line(const std::string& text,
                         std::size_t from,
                         std::string_view prefix,
                         std::string_view replacement) {
  const std::size_t at = text.find(prefix, from);
  REQUIRE(at != std::string::npos);
  const std::size_t end = text.find("\r\n", at);
  REQUIRE(end != std::string::npos);
  return text.substr(0, at) + std::string(replacement) + text.substr(end);
}

/// Surface block `k` of `text`: its line starting with `prefix` (two blanks of indentation).
std::string in_surface(const std::string& text,
                       int k,
                       std::string_view prefix,
                       std::string_view replacement) {
  const std::size_t at = text.find("SURF " + std::to_string(k) + "\r\n");
  REQUIRE(at != std::string::npos);
  return replace_line(text, at, prefix, replacement);
}

/// `text` with everything from the line MNUM on replaced by `lines` (CRLF line ends).
std::string with_block(const std::string& text, const std::vector<std::string>& lines) {
  const std::size_t at = text.find("\r\nMNUM");
  REQUIRE(at != std::string::npos);
  std::string out = text.substr(0, at + 2);
  for (const std::string& l : lines) out += l + "\r\n";
  return out;
}

/// One operand row as OpticStudio writes it: type, Int1, configuration, value, Int2 (field 5).
std::string op(std::string_view type, int int1, int config, std::string_view value, int int2 = 0) {
  return std::string(type) + " " + std::to_string(int1) + " " + std::to_string(config) + " " +
         std::string(value) + " 0 " + std::to_string(int2) + " 0 1 1 1.0E+000 0.0E+000 0";
}

/// The singlet of make_zmx_cases.py with an air dummy surface without aperture between its back
/// surface (SURF 3, DISZ 2) and the image (now SURF 5): the import drops the dummy and adds its
/// thickness `t4` to the 2 mm of SURF 3 (ADR 0044, point 5). Front curvature `c2`.
std::string zoom(std::string_view t4, std::string_view c2) {
  std::string t = read_file(kZmxDir + "singlet.zmx");
  t = in_surface(t, 3, "  DISZ", "  DISZ 2.0E+000");
  t = in_surface(t, 2, "  CURV", "  CURV " + std::string(c2) + " 0 0 0 0");
  const std::size_t image = t.find("SURF 4\r\n");
  REQUIRE(image != std::string::npos);
  const std::string dummy = "SURF 4\r\n  TYPE STANDARD\r\n  CURV 0 0 0 0 0\r\n  DISZ " +
                            std::string(t4) + "\r\n  DIAM 0 0 0 0 1\r\n";
  return t.substr(0, image) + dummy + "SURF 5\r\n" + t.substr(image + 8);
}

/// The system of configuration k without the table, as a single configuration import gives it.
System configuration(const System& s, std::size_t k) {
  System r = rtt::model::resolve_parameters(s, k);
  r.configurations.clear();
  r.parameters.clear();
  return r;
}

std::vector<double> row_values(const System& s, std::string_view name) {
  const auto row = rtt::model::find_parameter(s, name);
  REQUIRE(row);
  const auto* v = std::get_if<std::vector<double>>(&s.parameters[*row].form);
  REQUIRE(v != nullptr);
  return *v;
}

double rms(const System& s, std::size_t k) {
  const auto d = rtt::analysis::spot(rtt::compile::compile(s, schott(), k), {0}, 1, std::nullopt);
  REQUIRE(d.rays_arrived > 0);
  return d.stats.rms_centroid;
}

void expect_error(const std::string& text, const std::string& code, std::string_view needle) {
  INFO(code << " / " << needle);
  try {
    (void)import_zmx(text, schott());
    FAIL("no ZmxError");
  } catch (const ZmxError& e) {
    INFO(e.what());
    REQUIRE(e.code() == code);
    if (!needle.empty()) {
      // the reported line holds `needle`
      std::istringstream lines(text);
      std::string line;
      for (std::size_t n = 1; n <= e.line() && std::getline(lines, line); ++n) {
      }
      INFO("line " << e.line() << ": " << line);
      REQUIRE(line.find(needle) != std::string::npos);
    }
  }
}

}  // namespace

TEST_CASE("a zoom over a dropped dummy surface becomes rows of the parameter table",
          "[zmx][config][361]") {
  // THIC of the dummy surface 4 and CRVT of the front: three configurations, the second active
  // (MNUM 3 2), so the surface data hold the values of C2.
  const std::vector<std::string> t4{"4.3E+001", "3.0E+001", "5.0E+001"};
  const std::vector<std::string> c2{"2.0E-002", "2.5E-002", "1.8E-002"};
  std::vector<std::string> block{"MNUM 3 2", "MOFF 0 1 \"\" 0 0 0 1 1 0 0.0 "};
  for (int k = 1; k <= 3; ++k) block.push_back(op("THIC", 4, k, t4[k - 1]));
  for (int k = 1; k <= 3; ++k) block.push_back(op("CRVT", 2, k, c2[k - 1]));
  block.push_back(op("THIC", 3, 1, "2.0E+000"));  // the same in every configuration: no row
  block.push_back(op("THIC", 3, 2, "2.0E+000"));
  block.push_back(op("THIC", 3, 3, "2.0E+000"));
  block.push_back(op("MCOM", 0, 1, "\"wide\""));
  const ZmxImport multi = import_zmx(with_block(zoom(t4[1], c2[1]), block), schott());
  const System& s = multi.system;
  REQUIRE(s.configurations.size() == 3);
  REQUIRE(s.configurations[0].name == "C1");
  REQUIRE(s.configurations[2].name == "C3");
  // two rows: the position of the image (2 + THIC 4) and the radius of the front (1 / c)
  REQUIRE(s.parameters.size() == 2);
  for (std::size_t k = 0; k < 3; ++k) {
    INFO("configuration " << k + 1);
    const ZmxImport single = import_zmx(with_block(zoom(t4[k], c2[k]), {"MNUM 1"}), schott());
    REQUIRE(configuration(s, k) == single.system);
    REQUIRE(multi.warnings.size() == single.warnings.size());
  }
  // content guard: the rows hold the model values and the spot changes with the configuration
  const std::vector<double> image_z = row_values(s, "image_position_z");
  REQUIRE(image_z == std::vector<double>{2.0 + 43.0, 2.0 + 30.0, 2.0 + 50.0});
  const std::vector<double> radius = row_values(s, "S2_radius");
  REQUIRE(radius == std::vector<double>{1.0 / 2.0e-2, 1.0 / 2.5e-2, 1.0 / 1.8e-2});
  REQUIRE(rms(s, 0) != rms(s, 2));
}

TEST_CASE("asphere coefficients and the system aperture change per configuration",
          "[zmx][config][361]") {
  // asphere.zmx (FLOA): PRAM 3 3 (alpha_3 of the EVENASPH surface 3) and its obsolete form PAR2;
  // MNUM 2 without the active configuration: the surface data hold C1.
  const std::string base = read_file(kZmxDir + "asphere.zmx");
  const std::vector<std::string> a3{"3.000000000000E-010", "5.0E-010"};
  const std::vector<std::string> a2{"-1.500000000000E-006", "-1.0E-006"};
  const auto file = [&](std::size_t k) {
    std::string t = in_surface(base, 3, "  PARM 3", "  PARM 3 " + a3[k]);
    return in_surface(t, 3, "  PARM 2", "  PARM 2 " + a2[k]);
  };
  const std::vector<std::string> block{"MNUM 2", op("PRAM", 3, 1, a3[0], 3),
                                       op("PRAM", 3, 2, a3[1], 3), op("PAR2", 3, 1, a2[0]),
                                       op("PAR2", 3, 2, a2[1])};
  const ZmxImport multi = import_zmx(with_block(file(0), block), schott());
  REQUIRE(multi.system.configurations.size() == 2);
  REQUIRE(multi.system.parameters.size() == 2);
  for (std::size_t k = 0; k < 2; ++k) {
    INFO("configuration " << k + 1);
    const ZmxImport single = import_zmx(with_block(file(k), {"MNUM 1"}), schott());
    REQUIRE(configuration(multi.system, k) == single.system);
  }
  REQUIRE(rms(multi.system, 0) != rms(multi.system, 1));

  // APER on the singlet (ENPD): the entrance pupil diameter per configuration
  const std::string singlet = read_file(kZmxDir + "singlet.zmx");
  const auto pupil = [&](std::string_view d) {
    return replace_line(singlet, 0, "ENPD", "ENPD " + std::string(d));
  };
  const ZmxImport aper = import_zmx(
      with_block(pupil("1.0E+001"),
                 {"MNUM 2 1", op("APER", 0, 1, "1.0E+001"), op("APER", 0, 2, "6.0E+000")}),
      schott());
  REQUIRE(row_values(aper.system, "aperture_value") == std::vector<double>{10.0, 6.0});
  REQUIRE(configuration(aper.system, 0) ==
          import_zmx(with_block(pupil("1.0E+001"), {"MNUM 1"}), schott()).system);
  REQUIRE(configuration(aper.system, 1) ==
          import_zmx(with_block(pupil("6.0E+000"), {"MNUM 1"}), schott()).system);
  REQUIRE(rms(aper.system, 0) != rms(aper.system, 1));

  // CRVT of an ODDASPHE surface: the model holds the curvature itself (ADR 0043, point 1)
  const std::string odd = read_file(kZmxDir + "oddasph.zmx");
  const auto curved = [&](std::string_view c) {
    return in_surface(odd, 2, "  CURV", "  CURV " + std::string(c) + " 0 0 0 0");
  };
  const std::vector<std::string> c{"2.222222222222222307E-002", "2.0E-002"};
  const ZmxImport oc = import_zmx(
      with_block(curved(c[0]), {"MNUM 2", op("CRVT", 2, 1, c[0]), op("CRVT", 2, 2, c[1])}),
      schott());
  REQUIRE(row_values(oc.system, "S2_curvature").size() == 2);
  for (std::size_t k = 0; k < 2; ++k) {
    INFO("configuration " << k + 1);
    REQUIRE(configuration(oc.system, k) ==
            import_zmx(with_block(curved(c[k]), {"MNUM 1"}), schott()).system);
  }
  REQUIRE(rms(oc.system, 0) != rms(oc.system, 1));
}

TEST_CASE("the focal length of a paraxial lens before a fold changes per configuration",
          "[zmx][config][361]") {
  // fold.zmx: PRAM 1 1 (focal length of the PARAXIAL surface 1) through the folded path
  const std::string base = read_file(kZmxDir + "fold.zmx");
  const auto file = [&](std::string_view f) {
    return in_surface(base, 1, "  PARM 1", "  PARM 1 " + std::string(f));
  };
  const ZmxImport multi = import_zmx(
      with_block(file("5.0E+001"),
                 {"MNUM 2", op("PRAM", 1, 1, "5.0E+001", 1), op("PRAM", 1, 2, "6.0E+001", 1)}),
      schott());
  REQUIRE(multi.system.parameters.size() == 1);
  REQUIRE(multi.system.parameters[0].form ==
          rtt::model::ParameterForm{std::vector<double>{50.0, 60.0}});
  REQUIRE(configuration(multi.system, 0) ==
          import_zmx(with_block(file("5.0E+001"), {"MNUM 1"}), schott()).system);
  REQUIRE(configuration(multi.system, 1) ==
          import_zmx(with_block(file("6.0E+001"), {"MNUM 1"}), schott()).system);
}

TEST_CASE("configurations that change more than parameter values are rejected",
          "[zmx][config][361]") {
  const std::string z = zoom("4.3E+001", "2.0E-002");
  const auto two = [&](const std::string& text, std::vector<std::string> rows) {
    rows.insert(rows.begin(), "MNUM 2");
    return with_block(text, rows);
  };
  // a glass, a wavelength and a field of their own per configuration (ADR 0029: not
  // configurable)
  expect_error(two(z, {op("GLSS", 2, 1, "N-BK7"), op("GLSS", 2, 2, "N-SF6")}),
               "zmx.config_unsupported", "GLSS 2 1");
  expect_error(two(z, {op("WAVE", 1, 1, "4.861327E-001"), op("WAVE", 1, 2, "5.0E-001")}),
               "zmx.config_unsupported", "WAVE 1 1");
  expect_error(two(z, {op("YFIE", 2, 1, "5.0E+000"), op("YFIE", 2, 2, "6.0E+000")}),
               "zmx.config_unsupported", "YFIE 2 1");
  expect_error(two(z, {op("IGNR", 4, 1, "0"), op("IGNR", 4, 2, "1")}), "zmx.config_unsupported",
               "IGNR 4 1");
  // a tilt of a coordinate break (decision 6 on #361, follow-up #440)
  const std::string fold = read_file(kZmxDir + "fold.zmx");
  expect_error(two(fold, {op("PRAM", 2, 1, "4.5E+001", 3), op("PRAM", 2, 2, "4.0E+001", 3)}),
               "zmx.config_unsupported", "PRAM 2 1");
  // the thickness of the object at infinity
  expect_error(two(z, {op("THIC", 0, 1, "1.0E+010"), op("THIC", 0, 2, "1.0E+003")}),
               "zmx.config_unsupported", "THIC 0 1");
  // a curvature 0 in one configuration turns the sphere into a plane: another surface type
  expect_error(two(z, {op("CRVT", 2, 1, "2.0E-002"), op("CRVT", 2, 2, "0")}),
               "zmx.config_unsupported", "CRVT 2 1");
  // a parameter of a tilted plane (decision 6 on #361)
  const std::string plate = read_file(kZmxDir + "tilt_plate.zmx");
  expect_error(two(plate, {op("PRAM", 2, 1, "5.0E-002", 1), op("PRAM", 2, 2, "1.0E-001", 1)}),
               "zmx.config_unsupported", "PRAM 2 1");
  // APER with the system aperture floating by stop size (asphere.zmx has FLOA)
  const std::string floating = read_file(kZmxDir + "asphere.zmx");
  expect_error(two(floating, {op("APER", 0, 1, "5.0E+000"), op("APER", 0, 2, "4.0E+000")}),
               "zmx.config_unsupported", "APER 0 1");
  // the last asphere coefficient 0 in one configuration: the import cuts trailing zeros, the
  // coefficient list is shorter (decision on #361, a deviation from plan v1)
  expect_error(
      two(floating, {op("PRAM", 3, 1, "-2.000000000000E-013", 4), op("PRAM", 3, 2, "0", 4)}),
      "zmx.config_unsupported", "PRAM 3 1");
  // the same in every configuration: accepted without a row
  REQUIRE(import_zmx(two(z, {op("GLSS", 2, 1, "N-BK7"), op("GLSS", 2, 2, "N-BK7")}), schott())
              .system.parameters.empty());
}

TEST_CASE("a configuration block that contradicts the surface data or is incomplete",
          "[zmx][config][361]") {
  const std::string z = zoom("4.3E+001", "2.0E-002");
  // the active configuration (MNUM 2 2) differs from the surface data
  expect_error(
      with_block(z, {"MNUM 2 2", op("CRVT", 2, 1, "2.0E-002"), op("CRVT", 2, 2, "2.1E-002")}),
      "zmx.config_mismatch", "CRVT 2 2");
  // without the active configuration (MNUM 2) one configuration must agree with the surface
  // data, here none; in the OISC corpus the surface data hold C1, C2, C3 or C4
  expect_error(
      with_block(z, {"MNUM 2", op("THIC", 4, 1, "4.0E+001"), op("THIC", 4, 2, "3.0E+001")}),
      "zmx.config_mismatch", "MNUM 2");
  const ZmxImport second = import_zmx(
      with_block(z, {"MNUM 2", op("THIC", 4, 1, "3.0E+001"), op("THIC", 4, 2, "4.3E+001")}),
      schott());
  REQUIRE(row_values(second.system, "image_position_z") == std::vector<double>{32.0, 45.0});
  // a missing cell, a cell twice with different values, a configuration or active configuration
  // outside 1 ... n
  expect_error(with_block(z, {"MNUM 2", op("THIC", 4, 1, "4.3E+001")}), "zmx.multi_config",
               "MNUM 2");
  expect_error(with_block(z, {"MNUM 2", op("THIC", 4, 1, "4.3E+001"), op("THIC", 4, 2, "3.0E+001"),
                              op("THIC", 4, 2, "3.5E+001")}),
               "zmx.multi_config", "THIC 4 2 3.5");
  expect_error(
      with_block(z, {"MNUM 2", op("THIC", 4, 1, "4.3E+001"), op("THIC", 4, 3, "3.0E+001")}),
      "zmx.multi_config", "THIC 4 3");
  expect_error(with_block(z, {"MNUM 2 3"}), "zmx.multi_config", "MNUM 2 3");
  // a surface that does not exist, an operand line without a value
  expect_error(
      with_block(z, {"MNUM 2", op("THIC", 9, 1, "1.0E+000"), op("THIC", 9, 2, "2.0E+000")}),
      "zmx.config_mismatch", "THIC 9 1");
  expect_error(with_block(z, {"MNUM 2", "THIC 4 1"}), "zmx.parse", "THIC 4 1");
  expect_error(with_block(z, {"MNUM 161"}), "zmx.multi_config", "MNUM 161");
}
