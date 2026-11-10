// ZMX export of the odd asphere (ADR 0045, addendum #368): ODDASPHE against values written by
// hand, the choice between EVENASPH and ODDASPHE, the surfaces neither form holds, the odd
// asphere cases of tests/zmx through import, export and import, and the mirrored surface of a
// turned lens.

#include <catch2/catch_test_macros.hpp>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <optional>
#include <sstream>
#include <string>
#include <system_error>
#include <utility>
#include <variant>
#include <vector>

#include "rtt/analysis/spot.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/io/surface_table.hpp"
#include "rtt/io/zmx.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/system.hpp"

namespace {

using rtt::io::export_zmx;
using rtt::io::ExportError;
using rtt::io::import_zmx;
using rtt::io::ZmxExport;
using rtt::model::Element;
using rtt::model::OddAsphere;
using rtt::model::Param;
using rtt::model::Surface;
using rtt::model::System;

const std::string kZmxDir = std::string(RTT_ZMX_DIR) + "/";

std::string read_file(const std::string& path) {
  const std::ifstream in(path, std::ios::binary);
  REQUIRE(in);
  std::ostringstream buffer;
  buffer << in.rdbuf();
  return buffer.str();
}

/// SCHOTT, loaded once (MaterialLibrary is neither copyable nor movable).
const rtt::material::MaterialLibrary& schott() {
  static rtt::material::MaterialLibrary lib;
  static const bool loaded = [] {
    lib.add_catalog(std::string(RTT_CATALOG_DIR) + "/schott.agf");
    return true;
  }();
  REQUIRE(loaded);
  return lib;
}

double rel(double a, double b) {
  return std::abs(a - b) / std::abs(b);
}

Element& element(System& s, std::size_t i) {
  return std::get<Element>(s.root.children.at(i).value);
}

std::vector<std::string> codes_of(const std::vector<rtt::model::Diagnostic>& d) {
  std::vector<std::string> out;
  for (const auto& x : d) out.push_back(x.code);
  return out;
}

std::vector<std::string> lines_of(const std::string& text) {
  std::vector<std::string> out;
  std::size_t pos = 0;
  while (pos < text.size()) {
    const std::size_t end = text.find("\r\n", pos);
    REQUIRE(end != std::string::npos);
    out.push_back(text.substr(pos, end - pos));
    pos = end + 2;
  }
  return out;
}

std::vector<std::string> words_of(const std::string& line) {
  std::vector<std::string> out;
  std::istringstream in(line);
  for (std::string w; in >> w;) out.push_back(w);
  return out;
}

std::optional<double> as_number(const std::string& w) {
  double v = 0.0;
  const auto [ptr, ec] = std::from_chars(w.data(), w.data() + w.size(), v);
  if (ec != std::errc() || ptr != w.data() + w.size()) return std::nullopt;
  return v;
}

/// The surfaces of a flat tree in order.
std::vector<Surface*> surfaces_of(System& s) {
  std::vector<Surface*> out;
  for (rtt::model::Node& n : s.root.children) {
    for (Surface& surface : std::get<Element>(n.value).surfaces) out.push_back(&surface);
  }
  return out;
}

/// Radius of a conic or even asphere, none otherwise.
Param* radius_of(Surface& s) {
  if (auto* c = std::get_if<rtt::model::Conic>(&s.shape.base)) return &c->radius;
  if (auto* a = std::get_if<rtt::model::EvenAsphere>(&s.shape.base)) return &a->radius;
  return nullptr;
}

/// The export of `s` must throw ExportError with `code`.
void require_export_error(const System& s, const std::string& code) {
  try {
    (void)export_zmx(s, schott());
    FAIL("no ExportError, expected " << code);
  } catch (const ExportError& e) {
    INFO(e.what());
    REQUIRE(e.code() == code);
  }
}

/// The words of the indented lines of `SURF n` in an exported file.
std::vector<std::vector<std::string>> surface_block(const std::string& text, int n) {
  const std::vector<std::string> lines = lines_of(text);
  const std::string head = "SURF " + std::to_string(n);
  std::vector<std::vector<std::string>> out;
  std::size_t i = 0;
  while (i < lines.size() && lines[i] != head) ++i;
  REQUIRE(i < lines.size());
  for (++i; i < lines.size() && lines[i].starts_with("  "); ++i) out.push_back(words_of(lines[i]));
  return out;
}

/// The value of `key` (and, for PARM, of parameter `index`) in a block; none if absent.
std::optional<double> block_value(const std::vector<std::vector<std::string>>& block,
                                  const std::string& key,
                                  const std::string& index = "") {
  for (const auto& w : block) {
    if (w.empty() || w[0] != key) continue;
    if (index.empty() && w.size() >= 2) return as_number(w[1]);
    if (!index.empty() && w.size() >= 3 && w[1] == index) return as_number(w[2]);
  }
  return std::nullopt;
}

std::string block_type(const std::vector<std::vector<std::string>>& block) {
  for (const auto& w : block) {
    if (w.size() == 2 && w[0] == "TYPE") return w[1];
  }
  return "";
}

OddAsphere odd(double c, double k, const std::vector<double>& a) {
  OddAsphere o{Param(c), Param(k), {}};
  for (const double x : a) o.coefficients.emplace_back(x);
  return o;
}

/// The singlet of make_export_cases.py with `o` as the back surface L1.S2 (SURF 3 in the file).
System singlet_with(const OddAsphere& o) {
  System s = rtt::io::load_system(kZmxDir + "export_singlet.rtt.json");
  element(s, 1).surfaces[1].shape.base = o;
  return s;
}

}  // namespace

TEST_CASE("an odd asphere is written as ODDASPHE with CURV, CONI and PARM 1 to 8",
          "[zmx][export][368]") {
  // UG "Odd Asphere": z = c r^2 / (1 + sqrt(1 - (1 + k) c^2 r^2)) + beta_1 r + ... + beta_8 r^8,
  // parameters 1 to 8 = beta_1 ... beta_8; in mm (UNIT MM) beta_n = A_n. Expected values written
  // by hand, every coefficient different, a zero in the middle.
  const std::vector<double> a{1e-3, -2e-4, 3e-6, -4e-7, 5e-9, 0.0, -7e-12, 8e-14};
  const System s = singlet_with(odd(-0.02, -0.5, a));
  const ZmxExport out = export_zmx(s, schott());
  INFO(out.text);
  REQUIRE(out.warnings.empty());
  const auto block = surface_block(out.text, 3);
  REQUIRE(block_type(block) == "ODDASPHE");
  REQUIRE(block_value(block, "CURV") == -0.02);
  REQUIRE(block_value(block, "CONI") == -0.5);
  for (std::size_t n = 1; n <= 8; ++n) {
    INFO("PARM " << n);
    REQUIRE(block_value(block, "PARM", std::to_string(n)) == a[n - 1]);
  }
  REQUIRE_FALSE(block_value(block, "PARM", "0"));
  REQUIRE_FALSE(block_value(block, "PARM", "9"));
  REQUIRE_FALSE(block_value(block, "XDAT", "1"));

  // back through the import: the same odd asphere, bitwise; A1 != 0 is a cone at the vertex
  const auto back = import_zmx(out.text, schott());
  REQUIRE(codes_of(back.warnings) == std::vector<std::string>{"zmx.paraxial_unavailable"});
  const auto& lens = std::get<Element>(back.system.root.children.at(1).value);
  REQUIRE(lens.surfaces.size() == 2);
  REQUIRE(std::get<OddAsphere>(lens.surfaces[1].shape.base) == odd(-0.02, -0.5, a));
}

TEST_CASE("EVENASPH comes first; ODDASPHE takes what EVENASPH cannot hold", "[zmx][export][368]") {
  // Decision of the maintainer (2026-11-10): a surface both forms hold stays EVENASPH (c = 0,
  // only even terms); ODDASPHE for c != 0 (EVENASPH with c != 0 comes back as even_asphere), A1
  // to A3 or odd terms. Either way the import gives the same odd asphere; zeros at the end of the
  // list are dropped.
  struct Case {
    const char* what;
    OddAsphere shape;
    const char* type;
    OddAsphere back;
  };
  const std::vector<Case> cases{
      {"c = 0, A4 and A6: both forms, EVENASPH", odd(0.0, 0.0, {0, 0, 0, 1e-6, 0, -2e-9}),
       "EVENASPH", odd(0.0, 0.0, {0, 0, 0, 1e-6, 0, -2e-9})},
      {"c = 0, A4 and A10: EVENASPH only", odd(0.0, 0.0, {0, 0, 0, 1e-6, 0, 0, 0, 0, 0, 3e-15}),
       "EVENASPH", odd(0.0, 0.0, {0, 0, 0, 1e-6, 0, 0, 0, 0, 0, 3e-15})},
      {"c = 0, A3: ODDASPHE", odd(0.0, 0.0, {0, 0, 2e-6}), "ODDASPHE", odd(0.0, 0.0, {0, 0, 2e-6})},
      {"c = 0, A4 and A5: ODDASPHE", odd(0.0, 0.0, {0, 0, 0, 1e-6, -3e-8}), "ODDASPHE",
       odd(0.0, 0.0, {0, 0, 0, 1e-6, -3e-8})},
      {"c != 0, A4 only: ODDASPHE", odd(-0.01, -1.0, {0, 0, 0, 1e-6}), "ODDASPHE",
       odd(-0.01, -1.0, {0, 0, 0, 1e-6})},
      {"c != 0, A2: ODDASPHE", odd(-0.01, 0.25, {0, 1e-5}), "ODDASPHE",
       odd(-0.01, 0.25, {0, 1e-5})},
      {"c != 0, zeros beyond A8: ODDASPHE, the zeros dropped",
       odd(-0.01, 0.0, {0, 0, 0, 1e-6, 0, 0, 0, 0, 0, 0}), "ODDASPHE",
       odd(-0.01, 0.0, {0, 0, 0, 1e-6})},
  };
  for (const Case& c : cases) {
    INFO(c.what);
    const ZmxExport out = export_zmx(singlet_with(c.shape), schott());
    INFO(out.text);
    REQUIRE(out.warnings.empty());
    REQUIRE(block_type(surface_block(out.text, 3)) == c.type);
    const auto back = import_zmx(out.text, schott());
    REQUIRE(back.warnings.empty());
    const auto& lens = std::get<Element>(back.system.root.children.at(1).value);
    REQUIRE(std::get<OddAsphere>(lens.surfaces[1].shape.base) == c.back);
  }
}

TEST_CASE("an odd asphere neither form holds is not exported", "[zmx][export][368]") {
  for (const auto& [what, shape] : std::vector<std::pair<const char*, OddAsphere>>{
           {"c != 0 and A10", odd(-0.01, 0.0, {0, 0, 0, 1e-6, 0, 0, 0, 0, 0, 1e-20})},
           {"c = 0, A3 and A10", odd(0.0, 0.0, {0, 0, 2e-6, 0, 0, 0, 0, 0, 0, 1e-20})},
           {"c = 0, A9", odd(0.0, 0.0, {0, 0, 0, 1e-6, 0, 0, 0, 0, 1e-19})},
           {"c = 0, A18",
            odd(0.0, 0.0, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1e-30})},
       }) {
    INFO(what);
    require_export_error(singlet_with(shape), "export.surface_unsupported");
  }
}

TEST_CASE("the odd asphere cases survive import, export and import", "[zmx][export][368]") {
  // Acceptance of #368: oddasph, oddasph_in and oddasph_cone through import, export and import,
  // bitwise; R of the STANDARD surfaces (1/R twice) and the wavelengths (n_air there and back) to
  // 1e-15 relative, as in acceptance 2.
  for (const char* name : {"oddasph", "oddasph_in", "oddasph_cone"}) {
    INFO(name);
    const auto first = import_zmx(read_file(kZmxDir + name + ".zmx"), schott());
    System a = first.system;
    const ZmxExport out = export_zmx(a, schott());
    INFO(out.text);
    REQUIRE(out.warnings.empty());
    REQUIRE(out.text.find("  TYPE ODDASPHE\r\n") != std::string::npos);
    const auto again = import_zmx(out.text, schott());
    REQUIRE(codes_of(again.warnings) == codes_of(first.warnings));
    System b = again.system;
    REQUIRE(element(a, 0).pose.position[2].value == 0.0);  // canonical, thicknesses bitwise

    REQUIRE(a.wavelengths.size() == b.wavelengths.size());
    for (std::size_t i = 0; i < a.wavelengths.size(); ++i) {
      REQUIRE(rel(b.wavelengths[i].um, a.wavelengths[i].um) <= 1e-15);
      b.wavelengths[i].um = a.wavelengths[i].um;
    }
    const auto sa = surfaces_of(a);
    const auto sb = surfaces_of(b);
    REQUIRE(sa.size() == sb.size());
    std::size_t odd_surfaces = 0;
    for (std::size_t i = 0; i < sa.size(); ++i) {
      if (b.aperture.stop && *b.aperture.stop == sb[i]->id) b.aperture.stop = sa[i]->id;
      sb[i]->id = sa[i]->id;
      Param* ra = radius_of(*sa[i]);
      Param* rb = radius_of(*sb[i]);
      if (ra != nullptr && rb != nullptr) {
        REQUIRE(rel(rb->value, ra->value) <= 1e-15);
        rb->value = ra->value;
      }
      // content guard: an odd asphere with an odd term, compared bitwise below
      if (const auto* o = std::get_if<OddAsphere>(&sb[i]->shape.base)) {
        REQUIRE(o->coefficients.size() >= 3);
        bool odd_term = false;
        for (std::size_t n = 0; n < o->coefficients.size(); n += 2) {
          odd_term = odd_term || o->coefficients[n].value != 0.0;  // A1, A3, ...
        }
        REQUIRE(odd_term);
        ++odd_surfaces;
      }
    }
    REQUIRE(odd_surfaces == 1);
    INFO("a:\n" << rtt::io::to_json(a) << "\nb:\n" << rtt::io::to_json(b));
    REQUIRE(a == b);

    // content guard: the rays arrive in both systems, the same number; the spot aims
    // paraxially, which the cone at the vertex of oddasph_cone does not have
    if (first.warnings.empty()) {
      const auto ca = rtt::compile::compile(a, schott());
      const auto cb = rtt::compile::compile(b, schott());
      const auto da = rtt::analysis::spot(ca, {0}, 1, std::nullopt);
      const auto db = rtt::analysis::spot(cb, {0}, 1, std::nullopt);
      REQUIRE(da.rays_arrived > 0);
      REQUIRE(db.rays_arrived == da.rays_arrived);
    }

    // the surface table stays without odd_asphere (decision G1 on #318)
    try {
      (void)rtt::io::write_surface_table(a, schott());
      FAIL("no ExportError from the table");
    } catch (const ExportError& e) {
      REQUIRE(e.code() == "export.surface_unsupported");
    }
  }
}

TEST_CASE("a turned lens writes its odd asphere mirrored", "[zmx][export][368]") {
  // ADR 0045, point 2: with the axis -e_z, c -> -c and A_n -> -A_n, the conic stays. The turned
  // lens has the same global surfaces as the plain one, so the file is the same.
  const std::vector<double> a{1e-3, -2e-4, 3e-6, -4e-7};
  const System plain = singlet_with(odd(-0.02, -0.5, a));
  System s = plain;
  Element& lens = element(s, 1);
  lens.pose.rotation_deg = {Param(180.0), Param(0.0), Param(0.0)};
  auto& front = std::get<rtt::model::EvenAsphere>(lens.surfaces[0].shape.base);
  front.radius = Param(-front.radius.value);
  for (Param& x : front.coefficients) x.value = -x.value;
  std::vector<double> minus;
  for (const double x : a) minus.push_back(-x);
  lens.surfaces[1].shape.base = odd(0.02, -0.5, minus);
  lens.surfaces[1].pose.position[2] = Param(-lens.surfaces[1].pose.position[2].value);
  element(s, 2).pose.position[2] = Param(-element(s, 2).pose.position[2].value);
  const std::string text = export_zmx(plain, schott()).text;
  REQUIRE(block_type(surface_block(text, 3)) == "ODDASPHE");
  REQUIRE(export_zmx(s, schott()).text == text);
}
