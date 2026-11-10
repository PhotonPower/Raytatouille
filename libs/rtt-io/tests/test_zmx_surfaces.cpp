// ZMX import of the surface types of S1 and S2 (#360, PR A): ODDASPHE and TOROIDAL. The cases
// come from make_zmx_cases.py, which writes each .zmx file together with the expected system
// (built in Python from the prescription) and the first order from an own y-nu trace, per
// section for the toroids (ADR 0044, addendum #360).

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
#include <variant>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
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
using rtt::math::Vec3;
using rtt::model::Element;
using rtt::model::System;
using rtt::trace::RayBatch;
using rtt::trace::RayStatus;

const std::string kZmxDir = std::string(RTT_ZMX_DIR) + "/";

std::string read_file(const std::string& path) {
  const std::ifstream in(path, std::ios::binary);
  REQUIRE(in);
  std::ostringstream buffer;
  buffer << in.rdbuf();
  return buffer.str();
}

/// The test catalogue, loaded once (MaterialLibrary is neither copyable nor movable).
const rtt::material::MaterialLibrary& schott() {
  static std::map<std::string, rtt::material::MaterialLibrary> libs;
  const auto [it, inserted] = libs.try_emplace("schott");
  if (inserted) it->second.add_catalog(std::string(RTT_CATALOG_DIR) + "/schott.agf");
  return it->second;
}

/// Parallel rays from z = -1 (in front of surface 1, AGENTS.md) on a grid inside radius r, in
/// the direction (tan ax, tan ay, 1).
std::vector<std::pair<Vec3, Vec3>> bundle(double r, double ax_deg, double ay_deg) {
  const double k = 3.14159265358979323846 / 180.0;
  const Vec3 d = Vec3(std::tan(ax_deg * k), std::tan(ay_deg * k), 1.0).normalized();
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

RayBatch trace_rays(const rtt::compile::CompiledSystem& cs,
                    const std::vector<std::pair<Vec3, Vec3>>& rays) {
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

double rel(double a, double b) {
  return std::abs(a - b) / std::abs(b);
}

ZmxImport import_case(const std::string& name) {
  return import_zmx(read_file(kZmxDir + name + ".zmx"), schott());
}

std::vector<std::string> codes_of(const std::vector<rtt::model::Diagnostic>& d) {
  std::vector<std::string> out;
  for (const auto& x : d) out.push_back(x.code);
  return out;
}

const Element& element(const System& s, std::size_t i) {
  return std::get<Element>(s.root.children.at(i).value);
}

const nlohmann::json& first_order_data() {
  static const nlohmann::json fo = nlohmann::json::parse(read_file(kZmxDir + "first_order.json"));
  return fo;
}

/// Replaces the line of surface block `k` that starts with `prefix` (after the indentation).
std::string in_surface(const std::string& text,
                       int k,
                       std::string_view prefix,
                       const std::string& replacement) {
  const std::size_t at = text.find("SURF " + std::to_string(k) + "\r\n");
  if (at == std::string::npos) throw std::logic_error("no surface");
  std::size_t pos = at;
  while (pos < text.size()) {
    const std::size_t end = text.find('\n', pos) + 1;
    std::size_t first = pos;
    while (text[first] == ' ') ++first;
    if (text.compare(first, prefix.size(), prefix) == 0) {
      return text.substr(0, pos) + replacement + "\r\n" + text.substr(end);
    }
    pos = end;
  }
  throw std::logic_error("prefix not found");
}

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

}  // namespace

TEST_CASE("the ODDASPHE and TOROIDAL cases import to their expected systems", "[zmx][360]") {
  // Bitwise as acceptance 1 of ADR 0044; wavelengths to 1e-14 (two Ciddor implementations).
  const std::map<std::string, std::vector<std::string>> warnings{
      {"oddasph", {}},
      {"oddasph_in", {}},
      {"oddasph_cone", {"zmx.paraxial_unavailable"}},
      {"toroid", {"zmx.paraxial_unavailable"}},
      {"toroid_in", {"zmx.paraxial_unavailable"}},
  };
  for (const auto& [name, codes] : warnings) {
    INFO(name);
    const ZmxImport r = import_case(name);
    REQUIRE(codes_of(r.warnings) == codes);
    System imported = r.system;
    const System expected = rtt::io::load_system(kZmxDir + name + ".rtt.json");
    REQUIRE(imported.wavelengths.size() == expected.wavelengths.size());
    for (std::size_t i = 0; i < expected.wavelengths.size(); ++i) {
      REQUIRE(rel(imported.wavelengths[i].um, expected.wavelengths[i].um) <= 1e-14);
      imported.wavelengths[i].um = expected.wavelengths[i].um;
    }
    INFO("imported:\n"
         << rtt::io::to_json(imported) << "\nexpected:\n"
         << rtt::io::to_json(expected));
    REQUIRE(imported == expected);
  }
}

TEST_CASE("the odd asphere meets the y-nu trace with the vertex curvature c + 2 beta2",
          "[zmx][360]") {
  const auto& fo = first_order_data();
  for (const char* name : {"oddasph", "oddasph_in"}) {
    INFO(name);
    const auto& e = fo.at(name);
    const auto cs = rtt::compile::compile(import_case(name).system, schott());
    const auto f = rtt::paraxial::first_order(cs, {0}, 1);
    REQUIRE(f.efl);
    REQUIRE(f.bfl);
    const double bound = 1e-12 * std::max(1.0, e.at("condition").get<double>());
    REQUIRE(rel(*f.efl, e.at("efl").get<double>()) <= bound);
    REQUIRE(rel(*f.bfl, e.at("bfl").get<double>()) <= bound);
    REQUIRE(rel(*f.rear_focal_length, e.at("efl_relative").get<double>()) <= bound);
  }
}

TEST_CASE("the toroids meet the y-nu trace in both sections", "[zmx][360]") {
  // x-z: the circle of the radius of rotation (0 for the cylinder); y-z: the profile. The
  // cylinder turned by 90 deg (decision H2) has its power in x.
  const auto& fo = first_order_data();
  for (const char* name : {"toroid", "toroid_in", "xcylinder"}) {
    INFO(name);
    const auto cs = rtt::compile::compile(import_case(name).system, schott());
    for (const auto& [key, section] :
         {std::pair{"x", rtt::paraxial::Section::x}, std::pair{"y", rtt::paraxial::Section::y}}) {
      INFO("section " << key);
      const auto& e = fo.at(name).at(key);
      const auto f = rtt::paraxial::first_order_section(cs, {0}, 0, section);
      REQUIRE(f.efl);
      REQUIRE(f.bfl);
      const double bound = 1e-12 * std::max(1.0, e.at("condition").get<double>());
      REQUIRE(rel(*f.efl, e.at("efl").get<double>()) <= bound);
      REQUIRE(rel(*f.bfl, e.at("bfl").get<double>()) <= bound);
    }
  }
}

TEST_CASE("the cases in inches agree with the millimetre files to the printed digits",
          "[zmx][360]") {
  // 13 printed digits: the radius of rotation and the coefficients to 1e-12 relative.
  const auto coefficients = [](const rtt::model::Surface& s) {
    std::vector<double> out;
    if (const auto* o = std::get_if<rtt::model::OddAsphere>(&s.shape.base)) {
      for (const auto& p : o->coefficients) out.push_back(p.value);
    } else if (const auto* t = std::get_if<rtt::model::Toroid>(&s.shape.base)) {
      for (const auto& p : t->coefficients) out.push_back(p.value);
      out.push_back(t->radius_of_rotation.value);
    } else if (const auto* c = std::get_if<rtt::model::Cylinder>(&s.shape.base)) {
      for (const auto& p : c->coefficients) out.push_back(p.value);
    }
    return out;
  };
  for (const auto& [mm, in] :
       {std::pair{"oddasph", "oddasph_in"}, std::pair{"toroid", "toroid_in"}}) {
    INFO(in);
    const System a = import_case(mm).system;
    const System b = import_case(in).system;
    std::size_t compared = 0;
    for (std::size_t i = 0; i < a.root.children.size(); ++i) {
      const Element& ea = element(a, i);
      const Element& eb = element(b, i);
      for (std::size_t k = 0; k < ea.surfaces.size(); ++k) {
        const auto ca = coefficients(ea.surfaces[k]);
        const auto cb = coefficients(eb.surfaces[k]);
        REQUIRE(ca.size() == cb.size());
        for (std::size_t j = 0; j < ca.size(); ++j) {
          if (ca[j] == 0.0) {
            REQUIRE(cb[j] == 0.0);
          } else {
            REQUIRE(rel(cb[j], ca[j]) <= 1e-12);
            ++compared;
          }
        }
      }
    }
    REQUIRE(compared >= 3);
  }
}

TEST_CASE("a cylinder with power in x is turned by 90 degrees about z", "[zmx][360]") {
  // Decision H2: TOROIDAL with an infinite y-z radius. The surface frame is turned, its
  // elliptical aperture swaps the half widths (OpticStudio tests it in the unturned frame), and
  // the next surface lies untilted on the axis again.
  const ZmxImport r = import_case("xcylinder");
  REQUIRE(codes_of(r.warnings) == std::vector<std::string>{"zmx.paraxial_unavailable"});
  const Element& lens = element(r.system, 1);
  const auto& front = lens.surfaces.at(0);
  const auto* cyl = std::get_if<rtt::model::Cylinder>(&front.shape.base);
  REQUIRE(cyl != nullptr);
  REQUIRE(cyl->radius.value == 70.0);
  REQUIRE(front.pose.rotation_deg[2].value == 90.0);
  const auto* ap = std::get_if<rtt::model::EllipticalAperture>(&front.aperture.value());
  REQUIRE(ap != nullptr);
  REQUIRE(ap->semi_axis_x == 5.0);
  REQUIRE(ap->semi_axis_y == 7.0);
  const auto cs = rtt::compile::compile(r.system, schott());
  REQUIRE(cs.surfaces().size() == 4);
  const std::array<double, 4> z{0.0, 3.0, 8.0, 98.0};
  for (std::size_t i = 0; i < 4; ++i) {
    INFO(cs.surfaces()[i].id.str());
    const auto& g = cs.surfaces()[i].to_global;
    REQUIRE(std::abs(g.translation().x()) <= 1e-12);
    REQUIRE(std::abs(g.translation().y()) <= 1e-12);
    REQUIRE(std::abs(g.translation().z() - z[i]) <= 1e-12);
    // every frame but the cylinder's is the global one again
    const rtt::math::Mat3 expected =
        i == 1 ? rtt::math::Mat3{{0.0, -1.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 0.0, 1.0}}
               : rtt::math::Mat3::Identity();
    REQUIRE((g.rotation() - expected).norm() <= 1e-12);
  }
}

TEST_CASE("a toroid as OpticStudio writes it, with extra data but no Zernike terms", "[zmx][360]") {
  // As in the ray-optics test file HoO-V2C18Ex46.zmx: PARM 0 (extrapolate), XDAT 1 = 0 terms,
  // XDAT 2 = normalization radius. It imports as the toroid without them.
  const std::string tor = read_file(kZmxDir + "toroid.zmx");
  const std::string real =
      in_surface(tor, 2, "PARM 1",
                 "  PARM 0 0\r\n  PARM 1 8.000000000000E+001\r\n"
                 "  XDAT 1 0.000000000000E+00 0 0 1.000000000000E+00 0.000000000000E+00 0 \"\"\r\n"
                 "  XDAT 2 1.000000000000E+02 0 0 1.000000000000E+00 0.000000000000E+00 0 \"\"");
  REQUIRE(real != tor);
  REQUIRE(import_zmx(real, schott()).system == import_case("toroid").system);
}

TEST_CASE("the turned cylinder swaps a rectangular aperture too", "[zmx][360]") {
  const std::string text = in_surface(read_file(kZmxDir + "xcylinder.zmx"), 2, "ELAP",
                                      "  SQAP 7.000000000000E+000 5.000000000000E+000 0");
  const System s = import_zmx(text, schott()).system;
  const auto& front = element(s, 1).surfaces.at(0);
  const auto* r = std::get_if<rtt::model::RectangularAperture>(&front.aperture.value());
  REQUIRE(r != nullptr);
  REQUIRE(r->half_width_x == 5.0);
  REQUIRE(r->half_width_y == 7.0);
}

TEST_CASE("content guard: the new cases trace, all rays arrive, the fields differ", "[zmx][360]") {
  // Given start rays: the toroids and the cone have no paraxial model for the aiming. Three
  // fields (on axis, about x, about y) give three different mean image points.
  for (const char* name : {"oddasph", "oddasph_cone", "toroid", "xcylinder"}) {
    INFO(name);
    const auto cs = rtt::compile::compile(import_case(name).system, schott());
    std::set<std::pair<double, double>> centres;
    for (const auto& [ax, ay] : {std::pair{0.0, 0.0}, std::pair{2.0, 0.0}, std::pair{0.0, 2.0}}) {
      const RayBatch b = trace_rays(cs, bundle(3.5, ax, ay));
      REQUIRE(b.size() > 40);
      Vec3 mean = Vec3::Zero();
      for (std::size_t i = 0; i < b.size(); ++i) {
        REQUIRE(b.status()[i] == RayStatus::Alive);
        mean += Vec3(b.pos_x()[i], b.pos_y()[i], b.pos_z()[i]);
      }
      mean /= static_cast<double>(b.size());
      centres.insert({mean.x(), mean.y()});
    }
    REQUIRE(centres.size() == 3);
  }
}

TEST_CASE("the rules of ODDASPHE and TOROIDAL reject what the model cannot hold", "[zmx][360]") {
  const std::string odd = read_file(kZmxDir + "oddasph.zmx");
  const std::string tor = read_file(kZmxDir + "toroid.zmx");
  struct Case {
    std::string code;
    std::string text;
    std::string needle;
  };
  const std::vector<Case> cases{
      // alpha1 (y^2) of the toroid: the model has no such term (ADR 0042, F2)
      {"zmx.asphere_r2", in_surface(tor, 2, "PARM 2", "  PARM 2 1.0E-003"), "PARM 2 1.0E-003"},
      // a parameter beyond the table of the type
      {"zmx.parse", in_surface(odd, 2, "PARM 8", "  PARM 9 1.0E-015"), "PARM 9 1.0E-015"},
      {"zmx.parse", in_surface(tor, 2, "PARM 4", "  PARM 9 1.0E-015"), "PARM 9 1.0E-015"},
      // a toroid with an infinite y-z radius and terms in y: neither a cylinder in x nor a toroid
      {"zmx.surface_unsupported",
       in_surface(tor, 2, "CURV", "  CURV 0.000000000000000000E+000 0 0 0 0"), "TYPE TOROIDAL"},
      // the same with radius of rotation 0: a profile flat in y with terms in y
      {"zmx.surface_unsupported",
       in_surface(in_surface(tor, 2, "CURV", "  CURV 0.000000000000000000E+000 0 0 0 0"), 2,
                  "PARM 1", "  PARM 1 0.000000000000E+000"),
       "TYPE TOROIDAL"},
      // Zernike terms on the toroid (extra data, UG parameter 13)
      {"zmx.surface_unsupported",
       in_surface(tor, 2, "PARM 1",
                  "  PARM 1 8.000000000000E+001\r\n  XDAT 1 3.000000000000E+00 0 0 1 0 0"),
       "XDAT 1 3"},
      // extra data on another type
      {"zmx.keyword_unsupported",
       in_surface(odd, 2, "PARM 1", "  PARM 1 0.0E+000\r\n  XDAT 1 0.0E+000 0 0 1 0 0"),
       "XDAT 1 0.0E+000"},
      // the turned cylinder as the image surface would turn the image coordinates
      {"zmx.surface_unsupported",
       in_surface(read_file(kZmxDir + "xcylinder.zmx"), 4, "TYPE",
                  "  TYPE TOROIDAL\r\n  PARM 1 5.000000000000E+001"),
       "PARM 1 5.000000000000E+001"},
  };
  for (const Case& c : cases) {
    INFO(c.code << " / " << c.needle);
    try {
      (void)import_zmx(c.text, schott());
      FAIL("no ZmxError");
    } catch (const ZmxError& e) {
      INFO(e.what());
      REQUIRE(e.code() == c.code);
      REQUIRE(e.line() == line_of(c.text, c.needle));
    }
  }
}

// ------------------------------------------------------------------- PR B: BICONICX, QED_TYPE

TEST_CASE("the BICONICX and QED_TYPE cases import to their expected systems", "[zmx][360]") {
  const std::map<std::string, std::vector<std::string>> warnings{
      {"biconic", {"zmx.paraxial_unavailable"}},
      {"biconic_in", {"zmx.paraxial_unavailable"}},
      {"qcon", {}},
      {"qcon_in", {}},
      {"qbfs", {}},
  };
  for (const auto& [name, codes] : warnings) {
    INFO(name);
    const ZmxImport r = import_case(name);
    REQUIRE(codes_of(r.warnings) == codes);
    System imported = r.system;
    const System expected = rtt::io::load_system(kZmxDir + name + ".rtt.json");
    REQUIRE(imported.wavelengths.size() == expected.wavelengths.size());
    for (std::size_t i = 0; i < expected.wavelengths.size(); ++i) {
      REQUIRE(rel(imported.wavelengths[i].um, expected.wavelengths[i].um) <= 1e-14);
      imported.wavelengths[i].um = expected.wavelengths[i].um;
    }
    INFO("imported:\n"
         << rtt::io::to_json(imported) << "\nexpected:\n"
         << rtt::io::to_json(expected));
    REQUIRE(imported == expected);
  }
}

TEST_CASE("the biconics meet the y-nu trace in both sections, the Q-con on axis", "[zmx][360]") {
  // Biconic x-z: 1/Rx (0 for Rx = 0), y-z: c; the biconic with an infinite Ry is the cylinder
  // turned by 90 deg (H2). Q-con: the vertex curvature is c (ADR 0043).
  const auto& fo = first_order_data();
  for (const char* name : {"biconic", "biconic_in", "xbiconic"}) {
    INFO(name);
    const auto cs = rtt::compile::compile(import_case(name).system, schott());
    for (const auto& [key, section] :
         {std::pair{"x", rtt::paraxial::Section::x}, std::pair{"y", rtt::paraxial::Section::y}}) {
      INFO("section " << key);
      const auto& e = fo.at(name).at(key);
      const auto f = rtt::paraxial::first_order_section(cs, {0}, 0, section);
      REQUIRE(f.efl);
      REQUIRE(f.bfl);
      const double bound = 1e-12 * std::max(1.0, e.at("condition").get<double>());
      REQUIRE(rel(*f.efl, e.at("efl").get<double>()) <= bound);
      REQUIRE(rel(*f.bfl, e.at("bfl").get<double>()) <= bound);
    }
  }
  for (const char* name : {"qcon", "qcon_in"}) {
    INFO(name);
    const auto& e = fo.at(name);
    const auto cs = rtt::compile::compile(import_case(name).system, schott());
    const auto f = rtt::paraxial::first_order(cs, {0}, 0);
    REQUIRE(f.efl);
    REQUIRE(f.bfl);
    const double bound = 1e-12 * std::max(1.0, e.at("condition").get<double>());
    REQUIRE(rel(*f.efl, e.at("efl").get<double>()) <= bound);
    REQUIRE(rel(*f.bfl, e.at("bfl").get<double>()) <= bound);
  }
}

TEST_CASE("the biconic and the Q-con in inches agree with the millimetre files", "[zmx][360]") {
  // 13 printed digits: radius in x and the Q values to 1e-12 relative.
  {
    const System sa = import_case("biconic").system;
    const System sb = import_case("biconic_in").system;
    const auto& a = std::get<rtt::model::Biconic>(element(sa, 1).surfaces[0].shape.base);
    const auto& b = std::get<rtt::model::Biconic>(element(sb, 1).surfaces[0].shape.base);
    REQUIRE(rel(b.radius_x.value, a.radius_x.value) <= 1e-12);
    REQUIRE(rel(b.radius_y.value, a.radius_y.value) <= 1e-12);
    REQUIRE(b.conic_x.value == a.conic_x.value);
  }
  const System sa = import_case("qcon").system;
  const System sb = import_case("qcon_in").system;
  const auto& a = std::get<rtt::model::QCon>(element(sa, 1).surfaces[0].shape.base);
  const auto& b = std::get<rtt::model::QCon>(element(sb, 1).surfaces[0].shape.base);
  REQUIRE(rel(b.normalization_radius.value, a.normalization_radius.value) <= 1e-12);
  REQUIRE(a.coefficients.size() == 3);
  REQUIRE(b.coefficients.size() == 3);
  for (std::size_t i = 0; i < 3; ++i) {
    REQUIRE(rel(b.coefficients[i].value, a.coefficients[i].value) <= 1e-12);
  }
}

TEST_CASE("a biconic with an infinite Ry is the cylinder in x, turned", "[zmx][360]") {
  const ZmxImport r = import_case("xbiconic");
  REQUIRE(codes_of(r.warnings) == std::vector<std::string>{"zmx.paraxial_unavailable"});
  const auto& front = element(r.system, 1).surfaces.at(0);
  const auto* cyl = std::get_if<rtt::model::Cylinder>(&front.shape.base);
  REQUIRE(cyl != nullptr);
  REQUIRE(cyl->radius.value == 60.0);
  REQUIRE(cyl->conic.value == -0.4);
  REQUIRE(front.pose.rotation_deg[2].value == 90.0);
  const auto* ap = std::get_if<rtt::model::RectangularAperture>(&front.aperture.value());
  REQUIRE(ap != nullptr);
  REQUIRE(ap->half_width_x == 5.0);
  REQUIRE(ap->half_width_y == 7.0);
}

TEST_CASE("a Q coefficient beyond the maximum term number is ignored, as in OpticStudio",
          "[zmx][360]") {
  const std::string q = read_file(kZmxDir + "qcon.zmx");
  const std::string more = in_surface(
      q, 2, "XDAT 6",
      "  XDAT 6 5.000000000000E-005 0 0 1 0 0\r\n  XDAT 7 9.000000000000E-001 0 0 1 0 0");
  REQUIRE(more != q);
  REQUIRE(import_zmx(more, schott()).system == import_case("qcon").system);
}

TEST_CASE("content guard: the PR B cases trace, all rays arrive, the fields differ", "[zmx][360]") {
  for (const char* name : {"biconic", "qcon", "qbfs", "xbiconic"}) {
    INFO(name);
    const auto cs = rtt::compile::compile(import_case(name).system, schott());
    std::set<std::pair<double, double>> centres;
    for (const auto& [ax, ay] : {std::pair{0.0, 0.0}, std::pair{2.0, 0.0}, std::pair{0.0, 2.0}}) {
      const RayBatch b = trace_rays(cs, bundle(3.5, ax, ay));
      REQUIRE(b.size() > 40);
      Vec3 mean = Vec3::Zero();
      for (std::size_t i = 0; i < b.size(); ++i) {
        REQUIRE(b.status()[i] == RayStatus::Alive);
        mean += Vec3(b.pos_x()[i], b.pos_y()[i], b.pos_z()[i]);
      }
      mean /= static_cast<double>(b.size());
      centres.insert({mean.x(), mean.y()});
    }
    REQUIRE(centres.size() == 3);
  }
}

TEST_CASE("the rules of BICONICX and QED_TYPE reject what the model cannot hold", "[zmx][360]") {
  const std::string bic = read_file(kZmxDir + "biconic.zmx");
  const std::string q = read_file(kZmxDir + "qcon.zmx");
  struct Case {
    std::string code;
    std::string text;
    std::string needle;
  };
  const std::vector<Case> cases{
      // a type flag other than 0 (Qbfs) and 1 (Qcon)
      {"zmx.parse", in_surface(q, 2, "XDAT 1", "  XDAT 1 2.000000000000E+000 0 0 1 0 0"),
       "XDAT 1 2.0"},
      // more than 100 terms (UG "Q-Type Asphere")
      {"zmx.parse", in_surface(q, 2, "XDAT 2", "  XDAT 2 1.010000000000E+002 0 0 1 0 0"),
       "XDAT 2 1.01"},
      // parameters: the Q type has none, the biconic two
      {"zmx.parse",
       in_surface(q, 2, "XDAT 1", "  PARM 3 1.0E-003\r\n  XDAT 1 1.000000000000E+000 0 0 1 0 0"),
       "PARM 3 1.0E-003"},
      {"zmx.parse", in_surface(bic, 2, "PARM 2", "  PARM 2 3.0E-001\r\n  PARM 3 1.0E-003"),
       "PARM 3 1.0E-003"},
      // extra data on a biconic
      {"zmx.keyword_unsupported",
       in_surface(bic, 2, "PARM 2", "  PARM 2 3.0E-001\r\n  XDAT 1 0.5 0 0 1 0 0"), "XDAT 1 0.5"},
      // the turned cylinder of a biconic as the image surface
      {"zmx.surface_unsupported",
       in_surface(read_file(kZmxDir + "xbiconic.zmx"), 4, "TYPE",
                  "  TYPE BICONICX\r\n  PARM 1 5.000000000000E+001"),
       "PARM 1 5.000000000000E+001"},
  };
  for (const Case& c : cases) {
    INFO(c.code << " / " << c.needle);
    try {
      (void)import_zmx(c.text, schott());
      FAIL("no ZmxError");
    } catch (const ZmxError& e) {
      INFO(e.what());
      REQUIRE(e.code() == c.code);
      REQUIRE(e.line() == line_of(c.text, c.needle));
    }
  }
}
