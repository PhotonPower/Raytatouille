// ZMX import (ADR 0044). The cases under libs/rtt-io/tests/zmx come from make_zmx_cases.py,
// which writes each .zmx file together with the expected system (built in Python from the
// prescription, not with this code) and the first order from an own y-nu trace (acceptance 1
// and 3). The error and warning cases change single lines of those files.

#include <algorithm>
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <map>
#include <nlohmann/json.hpp>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "rtt/analysis/spot.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/diagnostics/codes.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/io/zmx.hpp"
#include "rtt/material/material.hpp"
#include "rtt/math/types.hpp"
#include "rtt/model/system.hpp"
#include "rtt/paraxial/paraxial.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/sequential.hpp"

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

/// A library with one AGF catalogue, loaded once (MaterialLibrary is neither copyable nor
/// movable).
const rtt::material::MaterialLibrary& catalogue(const std::string& file) {
  static std::map<std::string, rtt::material::MaterialLibrary> libs;
  const auto [it, inserted] = libs.try_emplace(file);
  if (inserted) it->second.add_catalog(std::string(RTT_CATALOG_DIR) + "/" + file);
  return it->second;
}

const rtt::material::MaterialLibrary& schott() {
  return catalogue("schott.agf");
}

const rtt::material::MaterialLibrary& schott_m2() {
  return catalogue("m2/schott.agf");
}

const rtt::material::MaterialLibrary& library_for(const std::string& name) {
  return name == "sasian" ? schott_m2() : schott();
}

double rel(double a, double b) {
  return std::abs(a - b) / std::abs(b);
}

/// Compares an import with the expected system: wavelengths to 1e-14 relative, F/# and NA to
/// 1e-15 relative (both hang on n_air from two Ciddor implementations), everything else
/// bitwise (ADR 0044, acceptance 1).
void require_matches(System imported, const System& expected) {
  INFO("imported:\n"
       << rtt::io::to_json(imported) << "\nexpected:\n"
       << rtt::io::to_json(expected));
  REQUIRE(imported.wavelengths.size() == expected.wavelengths.size());
  for (std::size_t i = 0; i < expected.wavelengths.size(); ++i) {
    REQUIRE(rel(imported.wavelengths[i].um, expected.wavelengths[i].um) <= 1e-14);
    imported.wavelengths[i].um = expected.wavelengths[i].um;
  }
  using T = rtt::model::SystemApertureType;
  if (expected.aperture.type == T::ImageSpaceFNumber ||
      expected.aperture.type == T::ObjectSpaceNA) {
    REQUIRE(imported.aperture.type == expected.aperture.type);
    REQUIRE(rel(imported.aperture.value.value, expected.aperture.value.value) <= 1e-15);
    imported.aperture.value.value = expected.aperture.value.value;
  }
  REQUIRE(imported == expected);
}

ZmxImport import_case(const std::string& name) {
  return import_zmx(read_file(kZmxDir + name + ".zmx"), library_for(name));
}

System expected_case(const std::string& name) {
  return rtt::io::load_system(kZmxDir + name + ".rtt.json");
}

std::vector<std::string> codes_of(const std::vector<rtt::model::Diagnostic>& d) {
  std::vector<std::string> out;
  for (const auto& x : d) out.push_back(x.code);
  return out;
}

const std::vector<std::string> kCases{"singlet", "doublet", "mirror",   "finite",
                                      "asphere", "sasian",  "immersion"};

}  // namespace

TEST_CASE("each ZMX case imports to the independently built system", "[zmx]") {
  const std::map<std::string, std::vector<std::string>> expected_warnings{
      {"doublet", {"zmx.aperture_dropped"}},
  };
  for (const std::string& name : kCases) {
    INFO(name);
    const ZmxImport r = import_case(name);
    const auto w = expected_warnings.find(name);
    REQUIRE(codes_of(r.warnings) ==
            (w == expected_warnings.end() ? std::vector<std::string>{} : w->second));
    require_matches(r.system, expected_case(name));
  }
}

TEST_CASE("a file in inches, centimetres or metres gives the system of its expected file",
          "[zmx]") {
  for (const char* name : {"singlet_in", "singlet_m", "asphere_in", "asphere_cm"}) {
    INFO(name);
    const ZmxImport r = import_case(name);
    REQUIRE(r.warnings.empty());
    require_matches(r.system, expected_case(name));
  }
}

TEST_CASE("inches and metres agree with the millimetre file to the printed digits", "[zmx]") {
  // The files print 13 significant digits; a value in inches times 25.4 meets the mm value to
  // about 5e-14 relative, a curvature c/s and R = 1/(c/s) likewise. 1e-12 leaves room.
  const System mm = import_case("singlet").system;
  for (const char* name : {"singlet_in", "singlet_m"}) {
    INFO(name);
    const System other = import_case(name).system;
    const auto& a = std::get<rtt::model::Element>(mm.root.children[1].value);
    const auto& b = std::get<rtt::model::Element>(other.root.children[1].value);
    const auto radius = [](const rtt::model::Element& e, std::size_t i) {
      return std::get<rtt::model::Conic>(e.surfaces[i].shape.base).radius.value;
    };
    REQUIRE(rel(radius(b, 0), radius(a, 0)) <= 1e-12);
    REQUIRE(rel(radius(b, 1), radius(a, 1)) <= 1e-12);
    REQUIRE(rel(b.surfaces[1].pose.position[2].value, a.surfaces[1].pose.position[2].value) <=
            1e-12);
    REQUIRE(rel(other.aperture.value.value, mm.aperture.value.value) <= 1e-12);
  }
}

TEST_CASE("asphere coefficients in inches and centimetres agree with the mm file", "[zmx]") {
  // alpha_i in L^(1 - 2i) scaled by (2i - 1) multiplications with 1/s (ADR 0044, point 3); the
  // files print 13 significant digits, so 1e-12 relative leaves room.
  const auto coefficients = [](const System& s, std::size_t element) {
    const auto& e = std::get<rtt::model::Element>(s.root.children[element].value);
    std::vector<double> out;
    for (const rtt::model::Surface& surface : e.surfaces) {
      if (const auto* odd = std::get_if<rtt::model::OddAsphere>(&surface.shape.base)) {
        for (const rtt::model::Param& p : odd->coefficients) out.push_back(p.value);
      }
      if (const auto* even = std::get_if<rtt::model::EvenAsphere>(&surface.shape.base)) {
        for (const rtt::model::Param& p : even->coefficients) out.push_back(p.value);
      }
    }
    return out;
  };
  const System mm = import_case("asphere").system;
  for (const char* name : {"asphere_in", "asphere_cm"}) {
    INFO(name);
    const System other = import_case(name).system;
    for (const std::size_t element : {std::size_t{0}, std::size_t{1}}) {
      const std::vector<double> a = coefficients(mm, element);
      const std::vector<double> b = coefficients(other, element);
      REQUIRE(a.size() == b.size());
      REQUIRE(!a.empty());
      for (std::size_t i = 0; i < a.size(); ++i) {
        if (a[i] == 0.0) {
          REQUIRE(b[i] == 0.0);
        } else {
          REQUIRE(rel(b[i], a[i]) <= 1e-12);
        }
      }
    }
  }
}

TEST_CASE("UTF-16 LE with BOM imports like the ASCII file", "[zmx]") {
  const std::string ascii = read_file(kZmxDir + "singlet.zmx");
  std::string utf16 = "\xFF\xFE";
  for (const char ch : ascii) {
    utf16.push_back(ch);
    utf16.push_back('\0');
  }
  const ZmxImport a = import_zmx(ascii, schott());
  const ZmxImport b = import_zmx(utf16, schott());
  REQUIRE(a.system == b.system);
}

TEST_CASE("first order of each case agrees with the independent y-nu trace", "[zmx]") {
  // make_zmx_cases.py: own Sellmeier and Ciddor, absolute and relative to air. Bound 1e-12
  // relative times the condition sum|y_i phi_i| / |y_1 Phi| (ADR 0044, acceptance 3).
  const nlohmann::json fo = nlohmann::json::parse(read_file(kZmxDir + "first_order.json"));
  for (const char* name : {"singlet", "singlet_in", "singlet_m", "doublet", "mirror", "finite",
                           "asphere", "asphere_in", "asphere_cm", "sasian", "immersion"}) {
    INFO(name);
    const auto& e = fo.at(name);
    const System s = import_case(name).system;
    const auto cs = rtt::compile::compile(s, library_for(name));
    std::uint16_t ref = 0;
    for (std::size_t i = 0; i < s.wavelengths.size(); ++i) {
      if (s.wavelengths[i].reference) ref = static_cast<std::uint16_t>(i);
    }
    const rtt::paraxial::FirstOrder f = rtt::paraxial::first_order(cs, {0}, ref);
    REQUIRE(f.efl);
    REQUIRE(f.bfl);
    REQUIRE(f.rear_focal_length);
    const double bound = 1e-12 * std::max(1.0, e.at("condition").get<double>());
    INFO("efl " << *f.efl << " bfl " << *f.bfl << " f'_R " << *f.rear_focal_length);
    REQUIRE(rel(*f.efl, e.at("efl").get<double>()) <= bound);
    REQUIRE(rel(*f.bfl, e.at("bfl").get<double>()) <= bound);
    REQUIRE(rel(*f.rear_focal_length, e.at("rear_focal_length").get<double>()) <= bound);
    // Relative to air: rear focal length = EFL of OpticStudio, efl n_air likewise, BFL equal.
    const double n_air = e.at("n_air").get<double>();
    REQUIRE(rel(*f.rear_focal_length, e.at("efl_relative").get<double>()) <= bound);
    REQUIRE(rel(*f.efl * n_air, e.at("efl_relative").get<double>()) <= bound);
    REQUIRE(rel(*f.bfl, e.at("bfl_relative").get<double>()) <= bound);
    // System apertures (acceptance 3): the EPD from image_fnumber is EFL_rel / F# of the file,
    // and NA / n of object space (absolute here, relative there) is the same slope, because
    // n_air cancels. Both against the independent script, not against the expected file.
    if (e.contains("epd_relative")) {
      REQUIRE(f.entrance_pupil);
      REQUIRE(f.entrance_pupil->diameter);
      REQUIRE(rel(*f.entrance_pupil->diameter, e.at("epd_relative").get<double>()) <= bound);
    }
    if (e.contains("na_over_n_relative")) {
      const auto& path = cs.path(rtt::compile::PathId{0});
      const double n_obj = cs.media()[path.events.front().medium_before].index[ref].real();
      REQUIRE(rel(s.aperture.value.value / n_obj, e.at("na_over_n_relative").get<double>()) <=
              1e-14);
    }
  }
}

TEST_CASE("the Sasian triplet as ZMX meets the ray-optics first order and the reference", "[zmx]") {
  // Acceptance 4: ray-optics 0.9.10, relative to air (docs/quellen.md, "Vergleichsrechnung
  // Cooke-Triplet"; values as in test_m2_acceptance.cpp), 1e-9 mm.
  const System s = import_case("sasian").system;
  const auto cs = rtt::compile::compile(s, schott_m2());
  const rtt::paraxial::FirstOrder f = rtt::paraxial::first_order(cs, {0}, 1);
  REQUIRE(f.rear_focal_length);
  REQUIRE(f.bfl);
  REQUIRE(f.ffl);
  REQUIRE(std::abs(*f.rear_focal_length - 50.00098928993641) <= 1e-9);
  REQUIRE(std::abs(*f.bfl - 41.23764895362289) <= 1e-9);
  REQUIRE(std::abs(*f.ffl - 37.102124967857094) <= 1e-9);

  // Against tests/reference/m2/cooke_triplet_surface_stop.rtt.json: element names, kinds and
  // materials exactly, surface ids through a fixed map, vertices after compile to 1e-12 mm (the
  // reference has absolute positions; 4.831 + 5.86 need not be 10.691 in double).
  const System ref = rtt::io::load_system(std::string(RTT_REFERENCE_DIR) +
                                          "/m2/cooke_triplet_surface_stop.rtt.json");
  REQUIRE(s.root.children.size() == ref.root.children.size());
  for (std::size_t i = 0; i < ref.root.children.size(); ++i) {
    const auto& a = std::get<rtt::model::Element>(s.root.children[i].value);
    const auto& b = std::get<rtt::model::Element>(ref.root.children[i].value);
    INFO(b.name);
    REQUIRE(a.name == b.name);
    REQUIRE(a.kind == b.kind);
    REQUIRE(a.material == b.material);
    REQUIRE(a.surfaces.size() == b.surfaces.size());
  }
  const std::map<std::string, std::string> ids{{"S1", "L1.S1"}, {"S2", "L1.S2"}, {"S3", "L2.S1"},
                                               {"S4", "L2.S2"}, {"S5", "L3.S1"}, {"S6", "L3.S2"},
                                               {"S7", "IMG"}};
  REQUIRE(s.aperture.stop);
  REQUIRE(ref.aperture.stop);
  REQUIRE(ids.at(s.aperture.stop->str()) == ref.aperture.stop->str());
  // Wavelengths: lambda_vac of the fixed point against the file to 1e-12 relative; radii
  // R = 1/c of the printed curvature to a few ulp (1e-15 relative).
  REQUIRE(s.wavelengths.size() == ref.wavelengths.size());
  for (std::size_t i = 0; i < ref.wavelengths.size(); ++i) {
    REQUIRE(rel(s.wavelengths[i].um, ref.wavelengths[i].um) <= 1e-12);
  }
  for (std::size_t i = 0; i + 1 < ref.root.children.size(); ++i) {
    const auto& a = std::get<rtt::model::Element>(s.root.children[i].value);
    const auto& b = std::get<rtt::model::Element>(ref.root.children[i].value);
    for (std::size_t k = 0; k < b.surfaces.size(); ++k) {
      const double ra = std::get<rtt::model::Conic>(a.surfaces[k].shape.base).radius.value;
      const double rb = std::get<rtt::model::Conic>(b.surfaces[k].shape.base).radius.value;
      REQUIRE(rel(ra, rb) <= 1e-15);
    }
  }
  const auto cr = rtt::compile::compile(ref, schott_m2());
  REQUIRE(cs.surfaces().size() == cr.surfaces().size());
  for (std::size_t i = 0; i < cr.surfaces().size(); ++i) {
    REQUIRE(ids.at(cs.surfaces()[i].id.str()) == cr.surfaces()[i].id.str());
    const auto pa = cs.surfaces()[i].to_global.apply_point({0.0, 0.0, 0.0});
    const auto pb = cr.surfaces()[i].to_global.apply_point({0.0, 0.0, 0.0});
    REQUIRE((pa - pb).norm() <= 1e-12);
  }
}

TEST_CASE("content guard: every case traces, all rays arrive, the fields differ", "[zmx]") {
  for (const std::string& name : kCases) {
    INFO(name);
    const System s = import_case(name).system;
    const auto cs = rtt::compile::compile(s, library_for(name));
    std::set<double> rms;
    for (std::uint16_t field = 0; field < s.fields.points.size(); ++field) {
      const rtt::analysis::SpotDiagram d = rtt::analysis::spot(cs, {0}, field, std::nullopt);
      INFO("field " << field);
      REQUIRE(d.rays_launched > 0);
      REQUIRE(d.rays_arrived == d.rays_launched);
      rms.insert(d.stats.rms_centroid);
    }
    REQUIRE(rms.size() == s.fields.points.size());
  }
}

namespace {

/// Replaces the first line of `text` that starts with `prefix` (after any indentation) by
/// `replacement` (several lines with \n allowed; empty removes the line).
std::string with(std::string text, std::string_view prefix, std::string_view replacement) {
  std::size_t pos = 0;
  while (pos < text.size()) {
    const std::size_t end = text.find('\n', pos);
    const std::size_t stop = end == std::string::npos ? text.size() : end + 1;
    std::size_t first = pos;
    while (first < stop && text[first] == ' ') ++first;
    if (text.compare(first, prefix.size(), prefix) == 0) {
      std::string rep(replacement);
      if (!rep.empty()) rep += "\r\n";
      text.replace(pos, stop - pos, rep);
      return text;
    }
    pos = stop;
  }
  // A throw, not FAIL plus return: MSVC flags the return as unreachable (C4702).
  throw std::logic_error("prefix not found: " + std::string(prefix));
}

/// 1-based number of the first line containing `needle`.
std::size_t line_of(const std::string& text, std::string_view needle) {
  std::istringstream in(text);
  std::string line;
  std::size_t n = 0;
  while (std::getline(in, line)) {
    ++n;
    if (line.find(needle) != std::string::npos) return n;
  }
  throw std::logic_error("needle not found");
}

std::string singlet() {
  return read_file(kZmxDir + "singlet.zmx");
}

/// doublet.zmx with a curvature on its first surface, a dummy without aperture that is dropped.
std::string doublet_curved_dummy() {
  const std::string t = read_file(kZmxDir + "doublet.zmx");
  const std::size_t at = t.find("SURF 1\r\n");
  REQUIRE(at != std::string::npos);
  return t.substr(0, at) + with(t.substr(at), "CURV", "  CURV 1.0E-002 0 0 0 0");
}

/// Surface block `k` of `text`: replaces its line starting with `prefix`.
std::string in_surface(const std::string& text,
                       int k,
                       std::string_view prefix,
                       std::string_view replacement) {
  const std::string& t = text;
  const std::string head = "SURF " + std::to_string(k) + "\r\n";
  const std::size_t at = t.find(head);
  REQUIRE(at != std::string::npos);
  return t.substr(0, at) + with(t.substr(at), prefix, replacement);
}

/// Surface block `k` of the singlet.
std::string in_surface(int k, std::string_view prefix, std::string_view replacement) {
  return in_surface(singlet(), k, prefix, replacement);
}

struct ErrorCase {
  std::string code;
  std::string text;
  std::string needle;  ///< text of the offending line; empty: line 0 (whole file)
  bool m2 = false;
};

std::vector<ErrorCase> error_cases() {
  const std::string s = singlet();
  const std::string asph = read_file(kZmxDir + "asphere.zmx");
  const std::string fold = read_file(kZmxDir + "fold.zmx");
  const std::string par = read_file(kZmxDir + "paraxial.zmx");
  const std::string plate = read_file(kZmxDir + "tilt_plate.zmx");
  return {
      {"zmx.parse", in_surface(2, "DISZ", "  DISZ six"), "DISZ six"},
      {"zmx.encoding_unsupported", "\xFE\xFF" + s, ""},
      {"zmx.mode_unsupported", with(s, "MODE", "MODE NSC"), "MODE NSC"},
      {"zmx.unit_unsupported", with(s, "UNIT", "UNIT FT X W X CM MR"), "UNIT FT"},
      {"zmx.environment_unsupported", with(s, "ENVD", "ENVD 2.5E+001 1.0E+000 0"), "ENVD 2.5"},
      {"zmx.aperture_unsupported", with(s, "ENPD", "FNUM 5.0E+000 1"), "FNUM 5.0E+000 1"},
      {"zmx.aperture_unsupported", with(s, "ENPD", "OBNA 1.0E-001 1"), "OBNA 1.0E-001 1"},
      {"zmx.aperture_unsupported", with(s, "FTYP", "FTYP 0 1 3 3 0 0 0"), "FTYP 0 1"},
      {"zmx.aperture_unsupported", in_surface(2, "FLAP", "  FLAP 1.0E+000 8.0E+000 0"),
       "FLAP 1.0E+000"},
      {"zmx.aperture_unsupported", in_surface(2, "FLAP", "  SQOB 4.0E+000 4.0E+000 0"), "SQOB"},
      {"zmx.field_unsupported", with(s, "FTYP", "FTYP 3 0 3 3 0 0 0"), "FTYP 3"},
      {"zmx.vignetting_unsupported", with(s, "VDYN", "VDYN 0 1.0E-001 0 0 0 0 0 0 0 0 0 0"),
       "VDYN 0 1.0E-001"},
      {"zmx.multi_config", with(s, "MNUM", "MNUM 2"), "MNUM 2"},
      {"zmx.config_mismatch", with(s, "THIC", "THIC 3 1 4.4E+001 0 0 0 1 1 0"), "THIC 3 1 4.4"},
      {"zmx.surface_unsupported", in_surface(2, "TYPE", "  TYPE TOROIDAL"), "TYPE TOROIDAL"},
      {"zmx.asphere_r2", with(asph, "PARM 1", "  PARM 1 1.0E-003"), "PARM 1 1.0E-003"},
      {"zmx.object_medium_infinite", in_surface(0, "DISZ", "  DISZ INFINITY\r\n  GLAS N-BK7 0 0"),
       "GLAS N-BK7 0 0"},
      {"zmx.object_medium_unsupported",
       in_surface(0, "DISZ", "  DISZ 1.0E+002\r\n  GLAS N-BK7 0 0"), "GLAS N-BK7 0 0"},
      {"zmx.mirror_in_glass", in_surface(3, "DISZ", "  DISZ 4.5E+001\r\n  GLAS MIRROR 0 0"),
       "GLAS MIRROR"},
      {"zmx.image_in_glass", in_surface(3, "DISZ", "  DISZ 4.5E+001\r\n  GLAS F2 0 0"), "GLAS F2"},
      {"zmx.glass_unknown", in_surface(2, "GLAS", "  GLAS BK7 0 0 1.5 40"), "GLAS BK7"},
      {"zmx.model_glass", in_surface(2, "GLAS", "  GLAS N-BK7 1 0 1.5168 64.17"), "GLAS N-BK7 1"},
      {"zmx.model_glass", in_surface(2, "GLAS", "  GLAS ___BLANK 1 0 1.5168 64.17"), "___BLANK"},
      {"zmx.model_glass", in_surface(2, "GLAS", "  GLAS 517642 0 0 1.517 64.2"), "517642"},
      {"zmx.model_glass", in_surface(2, "GLAS", "  GLAS 517.642 0 0 1.517 64.2"), "517.642"},
      {"zmx.aperture_conflict", in_surface(2, "FLAP", "  FLAP 0 8.0E+000 0\r\n  CLAP 0 7.0E+000 0"),
       "CLAP 0 7"},
      {"zmx.aperture_invalid", in_surface(2, "FLAP", "  CLAP 0 0 0"), "CLAP 0 0 0"},
      {"zmx.aperture_invalid", in_surface(2, "FLAP", "  FLAP 0 0 0"), "FLAP 0 0 0"},
      {"zmx.aperture_invalid", in_surface(2, "FLAP", "  CLAP 5.0E+000 4.0E+000 0"), "CLAP 5.0"},
      {"zmx.aperture_decenter", in_surface(2, "FLAP", "  FLAP 0 8.0E+000 0\r\n  OBDC 1.0E-001 0"),
       "OBDC"},
      {"zmx.stop_radius_missing", in_surface(1, "DIAM", "  DIAM 0 0 0 0 1"), "DIAM 0 0 0 0 1"},
      {"zmx.stop_radius_conflict",
       with(in_surface(1, "DIAM", "  DIAM 5.0E+000 1 0 0 1\r\n  CLAP 0 4.0E+000 0"), "ENPD",
            "FLOA"),
       "CLAP 0 4"},
      {"zmx.invalid_result", with(s, "WAVM 1 ", "WAVM 1 4.861327E-001 -1.0E+000"), ""},
      // Review of PR B: further rules that reject a file.
      {"zmx.parse", in_surface(1, "STOP", ""), ""},
      {"zmx.parse", in_surface(in_surface(1, "STOP", ""), 4, "TYPE", "  STOP\r\n  TYPE STANDARD"),
       "STOP"},
      {"zmx.parse", with(s, "ENPD", "ENPD 1.0E+001\r\nFNUM 5.0E+000 0"), "FNUM 5.0E+000 0"},
      {"zmx.parse", with(s, "SURF 2", "XXXX 1\r\nSURF 2"), "SURF 2"},
      {"zmx.parse", in_surface(2, "CURV", "  CURV nan 0 0 0 0"), "CURV nan"},
      {"zmx.parse", in_surface(2, "DISZ", "  DISZ Infinity"), "DISZ Infinity"},
      {"zmx.surface_unsupported", in_surface(0, "CURV", "  CURV 1.0E-002 0 0 0 0"),
       "CURV 1.0E-002"},
      {"zmx.aperture_invalid", in_surface(2, "FLAP", "  CLAP -1.0E+000 4.0E+000 0"), "CLAP -1.0"},
      {"zmx.aperture_unsupported", in_surface(2, "FLAP", "  ELOB 4.0E+000 3.0E+000 0"), "ELOB"},
      {"zmx.aperture_unsupported", in_surface(2, "FLAP", "  OBSC 0 1.0E+000 0"), "OBSC"},
      {"zmx.config_mismatch", with(s, "THIC", "CRVT 2 1 2.1E-002 0 0 0 1 1 0"), "CRVT 2 1"},
      // PR C: coordinate breaks, tilted planes, the paraxial lens, apertures, surface keywords.
      {"zmx.coordbrk_surface",
       in_surface(in_surface(fold, 1, "STOP", ""), 2, "TYPE", "  STOP\r\n  TYPE COORDBRK"), "STOP"},
      {"zmx.coordbrk_surface", in_surface(fold, 2, "HIDE", "  CLAP 0 5.0E+000 0"),
       "CLAP 0 5.0E+000"},
      {"zmx.parse", in_surface(fold, 2, "HIDE", "  GLAS F2 0 0 1.5 40"), "GLAS F2"},
      {"zmx.parse", in_surface(fold, 2, "HIDE", "  GLAS MIRROR 0 0"), "GLAS MIRROR"},
      {"zmx.surface_unsupported", in_surface(4, "TYPE", "  TYPE COORDBRK"), "TYPE COORDBRK"},
      {"zmx.surface_unsupported", in_surface(0, "TYPE", "  TYPE COORDBRK"), "TYPE COORDBRK"},
      {"zmx.paraxial_medium", in_surface(par, 1, "HIDE", "  GLAS N-BK7 0 0 1.5 40"), "GLAS N-BK7"},
      {"zmx.surface_unsupported", in_surface(par, 1, "HIDE", "  GLAS MIRROR 0 0"), "GLAS MIRROR"},
      {"zmx.surface_unsupported", in_surface(par, 1, "CURV", "  CURV 1.0E-002 0 0 0 0"),
       "CURV 1.0E-002"},
      {"zmx.surface_unsupported", in_surface(4, "TYPE", "  TYPE PARAXIAL"), "TYPE PARAXIAL"},
      {"zmx.aperture_unsupported", in_surface(plate, 2, "HIDE", "  FLAP 0 5.0E+000 0"),
       "FLAP 0 5.0E+000"},
      {"zmx.aperture_unsupported",
       in_surface(in_surface(plate, 1, "STOP", ""), 2, "TYPE", "  STOP\r\n  TYPE TILTSURF"),
       "STOP"},
      {"zmx.aperture_unsupported", in_surface(1, "HIDE", "  SQAP 4.0E+000 4.0E+000 0"), "SQAP"},
      {"zmx.aperture_invalid", in_surface(2, "FLAP", "  SQAP 0 4.0E+000 0"), "SQAP 0"},
      {"zmx.aperture_invalid", in_surface(2, "FLAP", "  ELAP 4.0E+000 -1.0E+000 0"), "ELAP 4"},
      {"zmx.aperture_conflict",
       in_surface(2, "FLAP", "  CLAP 0 7.0E+000 0\r\n  SQAP 4.0E+000 4.0E+000 0"), "SQAP 4"},
      {"zmx.keyword_unsupported", in_surface(2, "HIDE", "  ZZZZ 1"), "ZZZZ 1"},
      {"zmx.keyword_unsupported", in_surface(2, "HIDE", "  VDSZ 1.0E+000 0 0"), "VDSZ 1"},
  };
}

}  // namespace

TEST_CASE("every rule that rejects a file throws its code at its line", "[zmx]") {
  std::set<std::string> covered;
  for (const ErrorCase& c : error_cases()) {
    INFO(c.code << " / " << c.needle);
    const std::size_t line = c.needle.empty() ? 0 : line_of(c.text, c.needle);
    try {
      (void)import_zmx(c.text, c.m2 ? schott_m2() : schott());
      FAIL("no ZmxError");
    } catch (const ZmxError& e) {
      INFO(e.what());
      REQUIRE(e.code() == c.code);
      REQUIRE(e.line() == line);
    }
    covered.insert(c.code);
  }
  for (const rtt::diagnostics::CodeInfo& info : rtt::diagnostics::kCodes) {
    if (info.producer != "zmx" || info.severity != rtt::diagnostics::Severity::Error) continue;
    INFO("no case for " << info.code);
    REQUIRE(covered.contains(std::string(info.code)));
  }
}

TEST_CASE("every warning of the import has a case at its line", "[zmx]") {
  const std::string s = singlet();
  // A glass with T_ref = 25 degC: N-BK7 of the test catalogue with another TD temperature.
  std::string agf = read_file(std::string(RTT_CATALOG_DIR) + "/schott.agf");
  const std::size_t td = agf.find("2.000000E+01", agf.find("NM N-BK7"));
  REQUIRE(td != std::string::npos);
  agf.replace(td, 12, "2.500000E+01");
  rtt::material::MaterialLibrary warm;
  warm.add_catalog_text(agf, "SCHOTT");

  struct WarningCase {
    std::string code;
    std::string text;
    std::string needle;
    const rtt::material::MaterialLibrary* lib;
  };
  const std::vector<WarningCase> cases{
      {"zmx.keyword_ignored", with(s, "ROPD", "ROPD 2\r\nZZZZ 1\r\nZZZZ 2"), "ZZZZ 1", &schott()},
      {"zmx.surface_dropped", doublet_curved_dummy(), "CURV 1.0E-002", &schott()},
      {"zmx.aperture_dropped",
       in_surface(in_surface(3, "DIAM", "  DIAM 8.0E+000 0 0 0 1"), 3, "FLAP",
                  "  FLAP 0 8.00E+000 0"),
       "FLAP 0 8.00E+000", &schott()},
      {"zmx.glass_reference_temperature", s, "GLAS N-BK7", &warm},
      {"zmx.afocal_image_space", with(s, "FTYP", "FTYP 0 0 3 3 0 0 1"), "FTYP 0 0 3 3 0 0 1",
       &schott()},
      {"zmx.stop_radius_small", with(s, "ENPD", "ENPD 1.2E+001"), "DIAM 5.000000000000E+000",
       &schott()},
      // PR C.
      {"zmx.paraxial_unavailable", read_file(kZmxDir + "fold.zmx"), "", &schott()},
      {"zmx.paraxial_opd",
       in_surface(read_file(kZmxDir + "paraxial.zmx"), 1, "PARM 2", "  PARM 2 3.0E+000"),
       "PARM 2 3.0E+000", &schott()},
      {"zmx.tilted_image", read_file(kZmxDir + "tilt_image.zmx"), "TYPE TILTSURF", &schott()},
  };
  std::set<std::string> covered;
  for (const WarningCase& c : cases) {
    INFO(c.code);
    const ZmxImport r = import_zmx(c.text, *c.lib);
    bool found = false;
    for (const auto& w : r.warnings) {
      INFO(w.code << " " << w.location << " " << w.message);
      if (w.code != c.code) continue;
      found = true;
      if (!c.needle.empty())
        REQUIRE(w.location == "line " + std::to_string(line_of(c.text, c.needle)));
    }
    REQUIRE(found);
    covered.insert(c.code);
  }
  for (const rtt::diagnostics::CodeInfo& info : rtt::diagnostics::kCodes) {
    if (info.producer != "zmx" || info.severity != rtt::diagnostics::Severity::Warning) continue;
    INFO("no case for " << info.code);
    REQUIRE(covered.contains(std::string(info.code)));
  }
}

TEST_CASE("the merit, tolerance and configuration blocks give no warnings", "[zmx]") {
  // singlet.zmx has a merit block (EFFL, BLNK, DMFS, TRAC, OPDX), TOL lines and MNUM 1 with
  // MOFF and a matching THIC.
  REQUIRE(import_case("singlet").warnings.empty());
}

TEST_CASE("keyword_ignored is reported once per keyword with its count", "[zmx]") {
  const std::string t = with(singlet(), "ROPD", "ROPD 2\r\nZZZZ 1\r\nZZZZ 2\r\nYYYY 3");
  const ZmxImport r = import_zmx(t, schott());
  REQUIRE(codes_of(r.warnings) ==
          std::vector<std::string>{"zmx.keyword_ignored", "zmx.keyword_ignored"});
  REQUIRE(r.warnings[0].message.find("ZZZZ") != std::string::npos);
  REQUIRE(r.warnings[0].message.find("2 times") != std::string::npos);
}

TEST_CASE("variants that import like the case they come from", "[zmx]") {
  const System base = import_case("singlet").system;
  // GLAS with field 2 = 3 (as at catalogue glasses in OpticStudio files) is a catalogue glass.
  REQUIRE(import_zmx(in_surface(2, "GLAS", "  GLAS N-BK7 3 0 1.62040997 60.32364916"), schott())
              .system == base);
  // EVENASPH with all alpha = 0 is the STANDARD surface.
  REQUIRE(import_zmx(in_surface(2, "TYPE", "  TYPE EVENASPH\r\n  PARM 2 0"), schott()).system ==
          base);
}

TEST_CASE("a curved stop between air and air becomes a stop element with its shape", "[zmx]") {
  const ZmxImport r = import_zmx(in_surface(1, "CURV", "  CURV 1.0E-002 0 0 0 0"), schott());
  REQUIRE(r.warnings.empty());
  const auto& stop = std::get<rtt::model::Element>(r.system.root.children[0].value);
  REQUIRE(stop.kind == rtt::model::ElementKind::Stop);
  const auto* conic = std::get_if<rtt::model::Conic>(&stop.surfaces[0].shape.base);
  REQUIRE(conic != nullptr);
  REQUIRE(conic->radius.value == 1.0 / 1.0E-002);
  const auto cs = rtt::compile::compile(r.system, schott());
  for (std::uint16_t field = 0; field < 3; ++field) {
    const auto d = rtt::analysis::spot(cs, {0}, field, std::nullopt);
    REQUIRE(d.rays_arrived == d.rays_launched);
  }
}

// ------------------------------------------------------------------------------------- PR C
// Coordinate breaks, tilted planes, fold mirrors, the paraxial lens and rectangular and
// elliptical apertures (ADR 0044, addendum "PR C"). The folded cases have no expected file:
// make_zmx_cases.py writes the frame of OpticStudio for every surface into frames.json (UG
// "Coordinate Break", "Tilted", "Thickness"), and the analytic anchors are computed here.

namespace {

using rtt::compile::CompiledSystem;
using rtt::math::Vec3;
using rtt::trace::RayBatch;
using rtt::trace::RayStatus;

const nlohmann::json& frames_json() {
  static const nlohmann::json j = nlohmann::json::parse(read_file(kZmxDir + "frames.json"));
  return j;
}

const rtt::compile::CompiledSurface& surface_by_id(const CompiledSystem& cs,
                                                   const std::string& id) {
  for (const auto& s : cs.surfaces()) {
    if (s.id.str() == id) return s;
  }
  throw std::logic_error("no surface " + id);
}

/// The element that holds surface `id` (a pointer: a reference result bound to a reference with
/// a temporary string argument trips GCC's -Wdangling-reference, AGENTS.md).
const rtt::model::Element* element_of(const System& s, const std::string& id) {
  for (const rtt::model::Node& n : s.root.children) {
    const auto& e = std::get<rtt::model::Element>(n.value);
    for (const rtt::model::Surface& surface : e.surfaces) {
      if (surface.id.str() == id) return &e;
    }
  }
  throw std::logic_error("no element with " + id);
}

Vec3 vec(const nlohmann::json& a) {
  return {a.at(0).get<double>(), a.at(1).get<double>(), a.at(2).get<double>()};
}

/// Traces start rays (global, unit directions) along the path; returns the batch after it.
RayBatch trace_rays(const CompiledSystem& cs, const std::vector<std::pair<Vec3, Vec3>>& rays) {
  RayBatch b(rays.size());
  for (std::size_t i = 0; i < rays.size(); ++i) {
    const auto& [p, d] = rays[i];
    b.pos_x()[i] = p.x();
    b.pos_y()[i] = p.y();
    b.pos_z()[i] = p.z();
    b.dir_x()[i] = d.x();
    b.dir_y()[i] = d.y();
    b.dir_z()[i] = d.z();
  }
  (void)rtt::trace::SequentialTracer().trace(cs, rtt::compile::PathId{0}, b);
  return b;
}

Vec3 pos(const RayBatch& b, std::size_t i) {
  return {b.pos_x()[i], b.pos_y()[i], b.pos_z()[i]};
}

/// Parallel rays from z = -1 (in front of surface 1, AGENTS.md) on a grid inside radius r, at a
/// field angle about x.
std::vector<std::pair<Vec3, Vec3>> parallel_bundle(double r, double angle_deg = 0.0) {
  const Vec3 d = Vec3(0.0, std::tan(angle_deg * 3.14159265358979323846 / 180.0), 1.0).normalized();
  std::vector<std::pair<Vec3, Vec3>> out;
  for (int iy = -4; iy <= 4; ++iy) {
    for (int ix = -4; ix <= 4; ++ix) {
      const double x = r * ix / 4.0;
      const double y = r * iy / 4.0;
      if (x * x + y * y <= r * r) out.emplace_back(Vec3(x, y, -1.0), d);
    }
  }
  return out;
}

// The last three: breaks in a row, multiplied out and decomposed (PR C, point 3, rule 3), at
// b = 90 deg as a composition (R02 = 1 +- 1 ulp) and with angles from a fixed seed; the frame of
// the thin element behind them checks the decomposition to 1e-14.
const std::vector<std::string> kFolded{"fold",       "oap",        "pivot",
                                       "tilt_image", "tilt_plate", "tilt_prism",
                                       "gimbal",     "gimbal_x",   "random_breaks"};

}  // namespace

TEST_CASE("the coaxial cases of PR C import to their expected systems", "[zmx]") {
  for (const char* name : {"paraxial", "apertures"}) {
    INFO(name);
    const ZmxImport r = import_case(name);
    REQUIRE(r.warnings.empty());
    require_matches(r.system, expected_case(name));
  }
}

TEST_CASE("first order of the PR C coaxial cases agrees with the y-nu trace", "[zmx]") {
  // As acceptance 3 of PR B; the paraxial lens enters the script with phi = |n| / f.
  const nlohmann::json fo = nlohmann::json::parse(read_file(kZmxDir + "first_order.json"));
  for (const char* name : {"paraxial", "apertures"}) {
    INFO(name);
    const auto& e = fo.at(name);
    const System s = import_case(name).system;
    const auto cs = rtt::compile::compile(s, schott());
    const rtt::paraxial::FirstOrder f = rtt::paraxial::first_order(cs, {0}, 0);
    REQUIRE(f.efl);
    REQUIRE(f.bfl);
    REQUIRE(f.rear_focal_length);
    const double bound = 1e-12 * std::max(1.0, e.at("condition").get<double>());
    REQUIRE(rel(*f.efl, e.at("efl").get<double>()) <= bound);
    REQUIRE(rel(*f.bfl, e.at("bfl").get<double>()) <= bound);
    REQUIRE(rel(*f.rear_focal_length, e.at("efl_relative").get<double>()) <= bound);
  }
}

TEST_CASE("content guard of the PR C coaxial cases", "[zmx]") {
  // apertures: all rays pass both apertures, the fields differ (as PR B). paraxial: the ideal
  // lens images every field to a point (spot < 1e-9 mm), at different image heights.
  {
    const System s = import_case("apertures").system;
    const auto cs = rtt::compile::compile(s, schott());
    std::set<double> rms;
    for (std::uint16_t field = 0; field < 3; ++field) {
      const auto d = rtt::analysis::spot(cs, {0}, field, std::nullopt);
      REQUIRE(d.rays_launched > 0);
      REQUIRE(d.rays_arrived == d.rays_launched);
      rms.insert(d.stats.rms_centroid);
    }
    REQUIRE(rms.size() == 3);
  }
  const System s = import_case("paraxial").system;
  const auto cs = rtt::compile::compile(s, schott());
  std::set<double> heights;
  for (std::uint16_t field = 0; field < 3; ++field) {
    const auto d = rtt::analysis::spot(cs, {0}, field, std::nullopt);
    REQUIRE(d.rays_launched > 0);
    REQUIRE(d.rays_arrived == d.rays_launched);
    REQUIRE(d.stats.rms_centroid < 1e-9);
    heights.insert(d.points.front().y);
  }
  REQUIRE(heights.size() == 3);
}

TEST_CASE("SQAP and ELAP clip at their half widths, the corner only by the ellipse", "[zmx]") {
  // apertures.zmx: SQAP 7 x 6 mm (printed in inches) on surface 1, ELAP 7 x 6 on surface 2.
  const System s = import_case("apertures").system;
  const auto& sq =
      std::get<rtt::model::RectangularAperture>(*element_of(s, "S1")->surfaces[0].aperture);
  const auto& el =
      std::get<rtt::model::EllipticalAperture>(*element_of(s, "S2")->surfaces[0].aperture);
  const double a = sq.half_width_x;
  const double b = sq.half_width_y;
  REQUIRE(el.semi_axis_x == a);
  REQUIRE(el.semi_axis_y == b);
  const auto cs = rtt::compile::compile(s, schott());
  const Vec3 z(0.0, 0.0, 1.0);
  const std::vector<std::pair<Vec3, Vec3>> rays{{Vec3(0.999 * a, 0.0, -1.0), z},
                                                {Vec3(1.001 * a, 0.0, -1.0), z},
                                                {Vec3(0.0, 0.999 * b, -1.0), z},
                                                {Vec3(0.0, 1.001 * b, -1.0), z},
                                                {Vec3(0.9 * a, 0.9 * b, -1.0), z}};
  const RayBatch out = trace_rays(cs, rays);
  CHECK(out.status()[0] == RayStatus::Alive);
  CHECK(out.status()[1] == RayStatus::Vignetted);
  CHECK(out.status()[2] == RayStatus::Alive);
  CHECK(out.status()[3] == RayStatus::Vignetted);
  // (0.9 a, 0.9 b) is inside the rectangle and outside the ellipse (0.81 + 0.81 > 1).
  CHECK(out.status()[4] == RayStatus::Vignetted);
  CHECK(out.last_surface()[4] == 1);
}

TEST_CASE("every surface of a folded case lies in the frame of OpticStudio", "[zmx]") {
  // Vertices to 1e-12 mm, axes (or the normal of a tilted plane) to 1e-14; the poses of a
  // decomposed transform round, so the frames are compared, not the poses.
  const std::map<std::string, std::vector<std::string>> expected_warnings{
      {"fold", {"zmx.paraxial_unavailable"}},
      {"oap", {"zmx.paraxial_unavailable"}},
      {"pivot", {"zmx.paraxial_unavailable"}},
      {"tilt_image", {"zmx.tilted_image", "zmx.paraxial_unavailable"}},
      {"tilt_plate", {"zmx.paraxial_unavailable"}},
      {"tilt_prism", {"zmx.paraxial_unavailable"}},
      {"gimbal", {"zmx.paraxial_unavailable"}},
      {"gimbal_x", {"zmx.paraxial_unavailable"}},
      {"random_breaks", {"zmx.paraxial_unavailable"}},
  };
  for (const std::string& name : kFolded) {
    INFO(name);
    const ZmxImport r = import_case(name);
    REQUIRE(codes_of(r.warnings) == expected_warnings.at(name));
    const auto cs = rtt::compile::compile(r.system, schott());
    std::size_t checked = 0;
    for (const auto& [id, f] : frames_json().at(name).at("frames").items()) {
      bool present = false;
      for (const auto& s : cs.surfaces()) present = present || s.id.str() == id;
      if (!present) continue;  // a dropped dummy
      INFO(id);
      const auto& g = surface_by_id(cs, id).to_global;
      REQUIRE((g.translation() - vec(f.at("origin"))).norm() <= 1e-12);
      if (f.contains("normal")) {
        REQUIRE((g.apply_vector({0.0, 0.0, 1.0}) - vec(f.at("normal"))).norm() <= 1e-14);
      } else {
        for (int i = 0; i < 3; ++i) {
          for (int j = 0; j < 3; ++j) {
            REQUIRE(std::abs(g.rotation()(i, j) - f.at("axes").at(i).at(j).get<double>()) <= 1e-14);
          }
        }
      }
      ++checked;
    }
    REQUIRE(checked == cs.surfaces().size());
  }
}

TEST_CASE("a single coordinate break gives an exact rotate_first pose", "[zmx]") {
  // PR C, point 3, rule 2: flag 0 -> pivot (dx, dy, t1), position (dx, dy, t1 + t2), the tilts
  // as in the file.
  using rtt::model::PoseOrder;
  using rtt::model::PoseReference;
  struct Expected {
    const char* file;
    const char* surface;
    PoseReference reference;
    std::array<double, 3> rotation, pivot, position;
  };
  const std::vector<Expected> cases{
      {"fold", "S3", PoseReference::RelativeToPreceding, {45, 0, 0}, {0, 0, 20}, {0, 0, 20}},
      {"fold", "S5", PoseReference::RelativeToPreceding, {45, 0, 0}, {0, 0, 0}, {0, 0, -32}},
      {"pivot", "S2", PoseReference::Absolute, {5, 0, 0}, {0, 0, 0}, {0, 0, 0}},
      {"pivot", "S6", PoseReference::RelativeToPreceding, {-5, 0, 0}, {0, 0, -9}, {0, 0, 141}},
      {"oap", "S3", PoseReference::RelativeToPreceding, {0, 0, 0}, {0, -10, 100}, {0, -10, 100}},
  };
  for (const Expected& c : cases) {
    INFO(c.file << " " << c.surface);
    const System s = import_case(c.file).system;
    const rtt::model::Pose p = element_of(s, c.surface)->pose;
    REQUIRE(p.reference == c.reference);
    REQUIRE(p.order == PoseOrder::RotateFirst);
    for (std::size_t i = 0; i < 3; ++i) {
      REQUIRE(p.rotation_deg[i].value == c.rotation[i]);
      REQUIRE(p.pivot[i] == c.pivot[i]);
      REQUIRE(p.position[i].value == c.position[i]);
    }
  }
  // The exit face of the tilted plate takes the tilt out: rotate_first about its vertex, the
  // thickness stays one number (PR C, point 4).
  const System plate = import_case("tilt_plate").system;
  const rtt::model::Element* lens = element_of(plate, "S3");
  REQUIRE(lens->surfaces[1].pose.order == PoseOrder::RotateFirst);
  REQUIRE(lens->surfaces[1].pose.position[2].value == 5.0);
}

TEST_CASE("the fold as ZMX meets tests/reference/m6/fold45_focus and focuses there", "[zmx]") {
  const System s = import_case("fold").system;
  const auto cs = rtt::compile::compile(s, schott());
  const System ref =
      rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m6/fold45_focus.rtt.json");
  const auto cr = rtt::compile::compile(ref, schott());
  // Vertices and normals; the normals up to the sign (the reference turns its detector by
  // -90 deg, the frame of OpticStudio by +90 deg).
  for (const auto& [mine, theirs] :
       std::vector<std::pair<std::string, std::string>>{{"S1", "IL"}, {"S3", "M"}, {"S5", "IMG"}}) {
    INFO(mine);
    const auto& a = surface_by_id(cs, mine).to_global;
    const auto& b = surface_by_id(cr, theirs).to_global;
    REQUIRE((a.translation() - b.translation()).norm() <= 1e-12);
    const Vec3 na = a.apply_vector({0.0, 0.0, 1.0});
    const Vec3 nb = b.apply_vector({0.0, 0.0, 1.0});
    REQUIRE(std::min((na - nb).norm(), (na + nb).norm()) <= 1e-14);
  }
  // The detector is 2 mm behind the focus (0, 30, 20): a ray at height h hits it 2 h / 50 off
  // (0, 32, 20); with DISZ -30 every parallel ray meets in (0, 30, 20).
  const RayBatch at32 = trace_rays(cs, {{Vec3(0.0, 4.0, -1.0), Vec3(0.0, 0.0, 1.0)},
                                        {Vec3(0.0, 0.0, -1.0), Vec3(0.0, 0.0, 1.0)}});
  REQUIRE(at32.status()[0] == RayStatus::Alive);
  REQUIRE((pos(at32, 1) - Vec3(0.0, 32.0, 20.0)).norm() <= 1e-12);
  REQUIRE(std::abs((pos(at32, 0) - Vec3(0.0, 32.0, 20.0)).norm() - 2.0 * 4.0 / 50.0) <= 1e-12);
  const std::string text = read_file(kZmxDir + "fold.zmx");
  const System s30 =
      import_zmx(in_surface(text, 4, "DISZ", "  DISZ -3.000000000000E+001"), schott()).system;
  const auto c30 = rtt::compile::compile(s30, schott());
  std::set<double> centres;
  for (const double angle : {0.0, 1.0}) {
    const RayBatch b = trace_rays(c30, parallel_bundle(4.0, angle));
    const Vec3 first = pos(b, 0);
    for (std::size_t i = 0; i < b.size(); ++i) {
      REQUIRE(b.status()[i] == RayStatus::Alive);
      REQUIRE((pos(b, i) - first).norm() < 1e-9);
    }
    if (angle == 0.0) REQUIRE((first - Vec3(0.0, 30.0, 20.0)).norm() <= 1e-12);
    centres.insert(first.z());  // content guard: the field moves the image point
  }
  REQUIRE(centres.size() == 2);
}

TEST_CASE("the off-axis parabola focuses in its focus, R/2 off the decentred vertex", "[zmx]") {
  // Vertex at (0, -10, 100) (decenter -10 in y), R = -200, k = -1: every ray parallel to the
  // axis meets in (0, -10, 100 - 100) = (0, -10, 0), on the image plane z = 0. With the wrong
  // sign of the decenter the focus would be at y = +10.
  const System s = import_case("oap").system;
  const auto cs = rtt::compile::compile(s, schott());
  const RayBatch b = trace_rays(cs, parallel_bundle(2.5));
  for (std::size_t i = 0; i < b.size(); ++i) {
    INFO(i);
    REQUIRE(b.status()[i] == RayStatus::Alive);
    REQUIRE((pos(b, i) - Vec3(0.0, -10.0, 0.0)).norm() <= 1e-9);
  }
}

TEST_CASE("a ray hits the tilted image plane on z = x tan_x + y tan_y", "[zmx]") {
  // UG "Tilted": the plane keeps the frame of the image (vertex (0, 0, 56), axes as global).
  const System s = import_case("tilt_image").system;
  const auto cs = rtt::compile::compile(s, schott());
  const RayBatch b = trace_rays(cs, {{Vec3(1.0, 3.0, -1.0), Vec3(0.0, 0.0, 1.0)}});
  REQUIRE(b.status()[0] == RayStatus::Alive);
  const Vec3 p = pos(b, 0) - Vec3(0.0, 0.0, 56.0);
  REQUIRE(std::abs(p.y()) > 0.1);
  REQUIRE(std::abs(p.z() - (0.1 * p.x() - 0.3 * p.y())) <= 1e-12);
}

TEST_CASE("a ray through the tilted plate and prism meets the hand-traced Snell ray", "[zmx]") {
  // make_zmx_cases.py traces (0, 1, -1) along +z with its own vector Snell's law, own Sellmeier
  // and Ciddor; 1e-11 mm leaves room for the rounding of two index computations.
  for (const char* name : {"tilt_plate", "tilt_prism"}) {
    INFO(name);
    const auto& ray = frames_json().at(name).at("ray");
    const System s = import_case(name).system;
    const auto cs = rtt::compile::compile(s, schott());
    const RayBatch b = trace_rays(cs, {{vec(ray.at("start")), Vec3(0.0, 0.0, 1.0)}});
    REQUIRE(b.status()[0] == RayStatus::Alive);
    REQUIRE((pos(b, 0) - vec(ray.at("hit"))).norm() <= 1e-11);
  }
}

TEST_CASE("variants of PR C that import like the case they come from", "[zmx]") {
  // Pickup and variable notes, TCED, SSID, FIMP and VDSZ with zeros in a surface are silent.
  const System base = import_case("singlet").system;
  const ZmxImport notes = import_zmx(
      in_surface(2, "HIDE",
                 "  HIDE 0 0 0 0 0 0 0 0\r\n  PPAR 1 2 1.0E+000 0 0\r\n  PZUP 1 -1 0 0\r\n"
                 "  VPAR 2\r\n  TCED 2.36E+001 0 0\r\n  SSID 7\r\n  FIMP\r\n  VDSZ 0 0 0"),
      schott());
  REQUIRE(notes.warnings.empty());
  REQUIRE(notes.system == base);
  // A GLAS line on a coordinate break that names the medium before it.
  const std::string prism = read_file(kZmxDir + "tilt_prism.zmx");
  REQUIRE(import_zmx(
              in_surface(prism, 3, "DISZ", "  DISZ 0.000000000000E+000\r\n  GLAS N-BK7 0 0 1.5 40"),
              schott())
              .system == import_case("tilt_prism").system);
  // OPD mode 0 of a paraxial lens: the same system, with the warning zmx.paraxial_opd.
  const std::string par = read_file(kZmxDir + "paraxial.zmx");
  const ZmxImport mode0 = import_zmx(in_surface(par, 1, "PARM 2", "  PARM 2 0"), schott());
  REQUIRE(codes_of(mode0.warnings) == std::vector<std::string>{"zmx.paraxial_opd"});
  REQUIRE(mode0.system == import_case("paraxial").system);
}
