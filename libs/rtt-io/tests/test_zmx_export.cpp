// ZMX export (ADR 0045, PR B): the singlet against the file of make_export_cases.py, written
// independently of this code (acceptance 4), the Cooke triplet through export and import
// (acceptance 1 and 8), the import cases of tests/zmx through import, export and import
// (acceptance 2), every reference system (acceptance 3), and one case per code (acceptance 6).

#include <algorithm>
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <variant>
#include <vector>

#include "rtt/analysis/spot.hpp"
#include "rtt/coating/catalog.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/diagnostics/codes.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/io/zmx.hpp"
#include "rtt/material/material.hpp"
#include "rtt/math/types.hpp"
#include "rtt/model/system.hpp"
#include "rtt/paraxial/paraxial.hpp"

namespace {

namespace fs = std::filesystem;
using rtt::io::export_zmx;
using rtt::io::ExportError;
using rtt::io::ExportOptions;
using rtt::io::import_zmx;
using rtt::io::ZmxExport;
using rtt::model::Element;
using rtt::model::ElementKind;
using rtt::model::EventKind;
using rtt::model::Node;
using rtt::model::Param;
using rtt::model::Surface;
using rtt::model::SurfaceId;
using rtt::model::System;

const std::string kZmxDir = std::string(RTT_ZMX_DIR) + "/";
const std::string kRefDir = std::string(RTT_REFERENCE_DIR) + "/";

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

/// The coatings of the reference systems (CoatingLibrary is neither copyable nor movable).
const rtt::coating::CoatingLibrary& coatings() {
  static rtt::coating::CoatingLibrary lib;
  static const bool loaded = [] {
    lib.add_catalog(std::string(RTT_CATALOG_DIR) + "/coatings/demo.json");
    lib.add_catalog(std::string(RTT_CATALOG_DIR) + "/coatings/m3.json");
    return true;
  }();
  REQUIRE(loaded);
  return lib;
}

const rtt::material::MaterialLibrary& library_for(const std::string& name) {
  return name == "sasian" ? schott_m2() : schott();
}

double rel(double a, double b) {
  return std::abs(a - b) / std::abs(b);
}

System singlet() {
  return rtt::io::load_system(kZmxDir + "export_singlet.rtt.json");
}

Element& element(System& s, std::size_t i) {
  return std::get<Element>(s.root.children.at(i).value);
}

std::uint16_t reference_of(const System& s) {
  std::uint16_t ref = 0;
  for (std::size_t i = 0; i < s.wavelengths.size(); ++i) {
    if (s.wavelengths[i].reference) ref = static_cast<std::uint16_t>(i);
  }
  return ref;
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
    REQUIRE(end != std::string::npos);  // every line ends with CRLF, also the last one
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

/// The surfaces of a flat tree in order (the systems of the import and the singlet are flat).
std::vector<Surface*> surfaces_of(System& s) {
  std::vector<Surface*> out;
  for (Node& n : s.root.children) {
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

/// Global vertices of all surfaces relative to the first one, in mm.
std::vector<rtt::math::Vec3> vertices(const rtt::compile::CompiledSystem& cs) {
  std::vector<rtt::math::Vec3> out;
  const rtt::math::Vec3 first = cs.surfaces().front().to_global.apply_point({0.0, 0.0, 0.0});
  for (const auto& s : cs.surfaces()) {
    out.emplace_back(s.to_global.apply_point({0.0, 0.0, 0.0}) - first);
  }
  return out;
}

/// The export of `s` must throw ExportError with `code`.
void require_export_error(const System& s,
                          const rtt::material::MaterialLibrary& lib,
                          const std::string& code,
                          const ExportOptions& options = {}) {
  try {
    (void)export_zmx(s, lib, options);
    FAIL("no ExportError, expected " << code);
  } catch (const ExportError& e) {
    INFO(e.what());
    REQUIRE(e.code() == code);
  }
}

const std::vector<std::string> kCoaxialCases{
    "singlet",    "singlet_in", "singlet_m", "doublet",   "mirror",   "finite",   "asphere",
    "asphere_in", "asphere_cm", "sasian",    "immersion", "paraxial", "apertures"};

}  // namespace

TEST_CASE("the singlet is written as the independently built file", "[zmx][export]") {
  // Acceptance 4: make_export_cases.py builds export_singlet.zmx from the prescription with its
  // own Ciddor. Word for word; numbers as doubles, bitwise except the air wavelengths and FNUM.
  const ZmxExport out = export_zmx(singlet(), schott());
  REQUIRE(out.warnings.empty());
  REQUIRE(out.text.ends_with("\r\n"));
  const std::vector<std::string> got = lines_of(out.text);
  const std::vector<std::string> want = lines_of(read_file(kZmxDir + "export_singlet.zmx"));
  REQUIRE(got.size() == want.size());
  std::size_t numbers = 0;
  for (std::size_t i = 0; i < want.size(); ++i) {
    INFO("line " << (i + 1) << ": '" << got[i] << "' vs '" << want[i] << "'");
    // indentation: the lines of a surface block start with blanks, the others do not
    REQUIRE(got[i].starts_with(" ") == want[i].starts_with(" "));
    const std::vector<std::string> g = words_of(got[i]);
    const std::vector<std::string> w = words_of(want[i]);
    REQUIRE(g.size() == w.size());
    for (std::size_t k = 0; k < w.size(); ++k) {
      const std::optional<double> a = as_number(g[k]);
      const std::optional<double> b = as_number(w[k]);
      if (!b) {
        REQUIRE(g[k] == w[k]);
        continue;
      }
      REQUIRE(a);
      ++numbers;
      const bool air = (w[0] == "WAVM" && k == 2) || (w[0] == "FNUM" && k == 1);
      if (air) {
        REQUIRE(rel(*a, *b) <= 1e-15);
      } else {
        REQUIRE(*a == *b);
      }
    }
  }
  REQUIRE(numbers > 50);
}

TEST_CASE("the Cooke triplet through ZMX keeps first order, spots and geometry", "[zmx][export]") {
  // Acceptance 1 and 8 (ADR 0045, point 9).
  const System s = rtt::io::load_system(kRefDir + "m2/cooke_triplet_surface_stop.rtt.json");
  const ZmxExport out = export_zmx(s, schott_m2());
  REQUIRE(out.warnings.empty());
  const auto back = import_zmx(out.text, schott_m2());
  REQUIRE(back.warnings.empty());

  const auto ca = rtt::compile::compile(s, schott_m2());
  const auto cb = rtt::compile::compile(back.system, schott_m2());
  const std::uint16_t ref = reference_of(s);
  const auto fa = rtt::paraxial::first_order(ca, {0}, ref);
  const auto fb = rtt::paraxial::first_order(cb, {0}, ref);
  REQUIRE(fa.efl);
  REQUIRE(fb.efl);
  REQUIRE(fa.bfl);
  REQUIRE(fb.bfl);
  REQUIRE(rel(*fb.efl, *fa.efl) <= 1e-12);
  REQUIRE(rel(*fb.bfl, *fa.bfl) <= 1e-12);

  std::set<double> rms;
  for (std::uint16_t field = 0; field < s.fields.points.size(); ++field) {
    INFO("field " << field);
    const auto da = rtt::analysis::spot(ca, {0}, field, std::nullopt);
    const auto db = rtt::analysis::spot(cb, {0}, field, std::nullopt);
    REQUIRE(da.rays_launched > 0);
    REQUIRE(da.rays_arrived == da.rays_launched);
    REQUIRE(db.rays_arrived == db.rays_launched);
    REQUIRE(std::abs(db.stats.rms_centroid - da.stats.rms_centroid) <= 1e-9);
    rms.insert(da.stats.rms_centroid);
  }
  REQUIRE(rms.size() == s.fields.points.size());

  const auto va = vertices(ca);
  const auto vb = vertices(cb);
  REQUIRE(va.size() == vb.size());
  for (std::size_t i = 0; i < va.size(); ++i) REQUIRE((va[i] - vb[i]).norm() <= 1e-12);
  REQUIRE(s.wavelengths.size() == back.system.wavelengths.size());
  for (std::size_t i = 0; i < s.wavelengths.size(); ++i) {
    REQUIRE(rel(back.system.wavelengths[i].um, s.wavelengths[i].um) <= 1e-15);
  }
  System a = s;
  System b = back.system;
  const auto sa = surfaces_of(a);
  const auto sb = surfaces_of(b);
  REQUIRE(sa.size() == sb.size());
  std::size_t radii = 0;
  for (std::size_t i = 0; i < sa.size(); ++i) {
    const Param* ra = radius_of(*sa[i]);
    const Param* rb = radius_of(*sb[i]);
    REQUIRE((ra == nullptr) == (rb == nullptr));
    if (ra != nullptr) {
      REQUIRE(rel(rb->value, ra->value) <= 1e-15);
      ++radii;
    }
  }
  REQUIRE(radii == 6);
}

TEST_CASE("the Cooke triplet with the stop element in L2 is rejected", "[zmx][export]") {
  const System s = rtt::io::load_system(kRefDir + "m2/cooke_triplet.rtt.json");
  try {
    (void)export_zmx(s, schott_m2());
    FAIL("no ExportError");
  } catch (const ExportError& e) {
    REQUIRE(e.code() == "export.element_in_glass");
    REQUIRE_FALSE(e.surface().empty());
  }
}

TEST_CASE("every coaxial ZMX case survives import, export and import", "[zmx][export]") {
  // Acceptance 2: through a map of the surfaces in tree order. Vertices and object distance
  // relative to the first vertex to 1e-12 mm; R, wavelengths, F/# and NA to 1e-15 relative;
  // everything else bitwise, the thicknesses too where no dummy surface was dropped (the first
  // element at z = 0).
  for (const std::string& name : kCoaxialCases) {
    INFO(name);
    const auto& lib = library_for(name);
    System a = import_zmx(read_file(kZmxDir + name + ".zmx"), lib).system;
    const ZmxExport out = export_zmx(a, lib);
    INFO(out.text);
    REQUIRE(out.warnings.empty());
    const auto again = import_zmx(out.text, lib);
    REQUIRE(again.warnings.empty());
    System b = again.system;

    REQUIRE(a.wavelengths.size() == b.wavelengths.size());
    for (std::size_t i = 0; i < a.wavelengths.size(); ++i) {
      REQUIRE(rel(b.wavelengths[i].um, a.wavelengths[i].um) <= 1e-15);
      b.wavelengths[i].um = a.wavelengths[i].um;
    }
    using T = rtt::model::SystemApertureType;
    if (a.aperture.type == T::ImageSpaceFNumber || a.aperture.type == T::ObjectSpaceNA) {
      REQUIRE(rel(b.aperture.value.value, a.aperture.value.value) <= 1e-15);
      b.aperture.value.value = a.aperture.value.value;
    }
    const auto sa = surfaces_of(a);
    const auto sb = surfaces_of(b);
    REQUIRE(sa.size() == sb.size());
    for (std::size_t i = 0; i < sa.size(); ++i) {
      if (b.aperture.stop && *b.aperture.stop == sb[i]->id) b.aperture.stop = sa[i]->id;
      sb[i]->id = sa[i]->id;
      Param* ra = radius_of(*sa[i]);
      Param* rb = radius_of(*sb[i]);
      if (ra != nullptr && rb != nullptr) {
        REQUIRE(rel(rb->value, ra->value) <= 1e-15);
        rb->value = ra->value;
      }
    }
    const bool canonical = element(a, 0).pose.position[2].value == 0.0;
    if (!canonical) {
      const auto ca = rtt::compile::compile(a, lib);
      const auto cb = rtt::compile::compile(b, lib);
      const auto va = vertices(ca);
      const auto vb = vertices(cb);
      for (std::size_t i = 0; i < va.size(); ++i) REQUIRE((va[i] - vb[i]).norm() <= 1e-12);
      const double z1a = ca.surfaces().front().to_global.translation().z();
      const double z1b = cb.surfaces().front().to_global.translation().z();
      if (!a.object.at_infinity) {
        REQUIRE(std::abs((z1b - cb.object_z()) - (z1a - ca.object_z())) <= 1e-12);
      }
      b.object.distance = a.object.distance;
      for (std::size_t i = 0; i < b.root.children.size(); ++i) {
        element(b, i).pose = element(a, i).pose;
      }
      for (std::size_t i = 0; i < sa.size(); ++i) sb[i]->pose = sa[i]->pose;
    }
    INFO("a:\n" << rtt::io::to_json(a) << "\nb:\n" << rtt::io::to_json(b));
    REQUIRE(a == b);
  }
}

TEST_CASE("the folded ZMX cases are not coaxial", "[zmx][export]") {
  for (const char* name : {"fold", "oap", "pivot", "tilt_plate", "gimbal"}) {
    INFO(name);
    const System s = import_zmx(read_file(kZmxDir + name + ".zmx"), schott()).system;
    require_export_error(s, schott(), "export.not_coaxial");
  }
}

TEST_CASE("every reference system exports and imports again, or fails with its export code",
          "[zmx][export]") {
  // Acceptance 3. A system that does not compile with either test catalogue is listed, so that
  // the list stays visible; the export never gives a file the import rejects.
  std::vector<fs::path> files;
  for (const auto& entry : fs::recursive_directory_iterator(RTT_REFERENCE_DIR)) {
    if (entry.is_regular_file() && entry.path().filename().string().ends_with(".rtt.json")) {
      files.push_back(entry.path());
    }
  }
  std::sort(files.begin(), files.end());
  std::map<std::string, std::string> outcome;  // relative path -> "ok" or the code
  for (const fs::path& file : files) {
    const std::string name = fs::relative(file, RTT_REFERENCE_DIR).generic_string();
    INFO(name);
    const System s = rtt::io::load_system(file.string());
    std::string result = "not compiled";
    for (const rtt::material::MaterialLibrary* lib : {&schott(), &schott_m2()}) {
      ZmxExport out;
      try {
        out = export_zmx(s, *lib);
      } catch (const ExportError& e) {
        REQUIRE(e.code().starts_with("export."));
        result = e.code();
        break;
      } catch (const rtt::compile::CompileError&) {
        continue;
      }
      result = "ok";
      const auto back = import_zmx(out.text, *lib);  // never throws for an exported file
      const auto ca = rtt::compile::compile(s, *lib, coatings());
      const auto cb = rtt::compile::compile(back.system, *lib);
      const std::uint16_t ref = reference_of(s);
      std::optional<rtt::paraxial::FirstOrder> fa;
      try {
        fa = rtt::paraxial::first_order(ca, {0}, ref);
      } catch (const std::exception&) {
      }
      if (fa) {
        const auto fb = rtt::paraxial::first_order(cb, {0}, ref);
        REQUIRE(fa->efl.has_value() == fb.efl.has_value());
        if (fa->efl) REQUIRE(rel(*fb.efl, *fa->efl) <= 1e-12);
        REQUIRE(fa->bfl.has_value() == fb.bfl.has_value());
        if (fa->bfl) REQUIRE(rel(*fb.bfl, *fa->bfl) <= 1e-12);
        for (std::uint16_t field = 0; field < s.fields.points.size(); ++field) {
          const auto da = rtt::analysis::spot(ca, {0}, field, std::nullopt);
          const auto db = rtt::analysis::spot(cb, {0}, field, std::nullopt);
          REQUIRE(db.rays_arrived == da.rays_arrived);
          REQUIRE(std::abs(db.stats.rms_centroid - da.stats.rms_centroid) <= 1e-9);
        }
      }
      break;
    }
    outcome[name] = result;
  }
  std::size_t ok = 0;
  for (const auto& [name, result] : outcome) {
    std::cout << "export " << name << ": " << result << "\n";
    if (result == "ok") ++ok;
  }
  // Content guard: the exportable systems of M0 to M6 in air go through.
  for (const char* name : {"m0/singlet.rtt.json", "m2/achromat.rtt.json",
                           "m2/cooke_triplet_surface_stop.rtt.json", "m2/paraboloid_stop.rtt.json",
                           "m6/paraboloid_f10.rtt.json", "m7/achromat_tolerances.rtt.json"}) {
    INFO(name);
    REQUIRE(outcome.at(name) == "ok");
  }
  REQUIRE(outcome.at("m1/paraboloid_mirror.rtt.json") == "export.stop_missing");
  REQUIRE(outcome.at("m2/cooke_triplet.rtt.json") == "export.element_in_glass");
  REQUIRE(outcome.at("m3/ar_mgf2.rtt.json") == "export.environment_unsupported");
  // Its explicit path refracts at the stop and the image (the default event kind), which keeps
  // the medium as transmit (ADR 0045, addendum PR B).
  REQUIRE(outcome.at("r2/ideal_lens.rtt.json") == "ok");
  REQUIRE(ok >= 7);
}

TEST_CASE("a lens turned by 180 degrees about x is written with mirrored shapes", "[zmx][export]") {
  // ADR 0045, point 2: with the axis -e_z, R -> -R and A_i -> -A_i, the conic stays. The turned
  // lens has the same global surfaces as the singlet, so the file is the same.
  System s = singlet();
  Element& lens = element(s, 1);
  lens.pose.rotation_deg = {Param(180.0), Param(0.0), Param(0.0)};
  auto& front = std::get<rtt::model::EvenAsphere>(lens.surfaces[0].shape.base);
  front.radius = Param(-front.radius.value);
  for (Param& a : front.coefficients) a = Param(-a.value);
  auto& back = std::get<rtt::model::Conic>(lens.surfaces[1].shape.base);
  back.radius = Param(-back.radius.value);
  lens.surfaces[1].pose.position[2] = Param(-lens.surfaces[1].pose.position[2].value);
  // the image is placed relative to the turned frame of L1.S2
  element(s, 2).pose.position[2] = Param(-element(s, 2).pose.position[2].value);
  REQUIRE(export_zmx(s, schott()).text == export_zmx(singlet(), schott()).text);
}

TEST_CASE("a decentre up to 1e-12 mm is dropped silently, more is not coaxial", "[zmx][export]") {
  System s = singlet();
  element(s, 1).pose.position[0] = Param(1e-13);
  const ZmxExport out = export_zmx(s, schott());
  REQUIRE(out.warnings.empty());
  REQUIRE(out.text == export_zmx(singlet(), schott()).text);
  element(s, 1).pose.position[0] = Param(1e-9);
  require_export_error(s, schott(), "export.not_coaxial");
  System t = singlet();
  element(t, 1).pose.rotation_deg[1] = Param(0.01);
  require_export_error(t, schott(), "export.not_coaxial");
}

TEST_CASE("a configuration is written with its resolved parameters", "[zmx][export]") {
  System s = singlet();
  s.configurations = {{"short"}, {"long"}};
  s.parameters.push_back({"t", std::vector<double>{6.0, 7.5}, false, std::nullopt, std::nullopt});
  element(s, 1).surfaces[1].pose.position[2] = Param::bound("t");
  require_export_error(s, schott(), "export.multi_config");
  const std::string first = export_zmx(s, schott(), {std::nullopt, 0}).text;
  const std::string second = export_zmx(s, schott(), {std::nullopt, 1}).text;
  REQUIRE(first == export_zmx(singlet(), schott()).text);
  REQUIRE(second.find("  DISZ 7.5\r\n") != std::string::npos);
  REQUIRE(second.find("  DISZ 6\r\n") == std::string::npos);
  REQUIRE_THROWS_AS(export_zmx(s, schott(), {std::nullopt, 2}), std::invalid_argument);
}

TEST_CASE("a path is chosen by name; the others are dropped with a warning", "[zmx][export]") {
  System s = singlet();
  rtt::model::Path p{"no stop", false, {}};
  for (const char* id : {"L1.S1", "L1.S2", "IMG"}) p.events.push_back({SurfaceId(id), {}, 0});
  p.events[2].kind = EventKind::Transmit;
  s.paths.push_back(p);
  const ZmxExport main_path = export_zmx(s, schott());
  REQUIRE(codes_of(main_path.warnings) == std::vector<std::string>{"export.paths_dropped"});
  try {
    (void)export_zmx(s, schott(), {std::string("no stop"), std::nullopt});
    FAIL("no ExportError");
  } catch (const ExportError& e) {
    REQUIRE(e.code() == "export.stop_missing");
  }
  REQUIRE_THROWS_AS(export_zmx(s, schott(), {std::string("none"), std::nullopt}),
                    std::invalid_argument);
}

TEST_CASE("refract at the stop and the detector is written as transmit", "[zmx][export]") {
  // ADR 0045, addendum PR B: refract (the default kind of an event) keeps the medium at an
  // element without material, so the explicit path gives the file of the automatic one.
  System s = singlet();
  s.paths = {{"main",
              false,
              {{SurfaceId("STO"), EventKind::Refract, 0},
               {SurfaceId("L1.S1"), EventKind::Refract, 0},
               {SurfaceId("L1.S2"), EventKind::Refract, 0},
               {SurfaceId("IMG"), EventKind::Refract, 0}}}};
  const ZmxExport out = export_zmx(s, schott());
  REQUIRE(out.warnings.empty());
  REQUIRE(out.text == export_zmx(singlet(), schott()).text);
  // a reflection at the stop stays unsupported
  s.paths.front().events[0].kind = EventKind::Reflect;
  require_export_error(s, schott(), "export.event_unsupported");
}

TEST_CASE("GCAT is ordered so that each glass is found in its own catalogue", "[zmx][export]") {
  // CZ holds N-BK7, CA holds N-BK7 and F2: L1 of CZ:N-BK7 needs CZ before CA, against the order
  // of the names.
  const std::string agf = read_file(std::string(RTT_CATALOG_DIR) + "/schott.agf");
  const std::size_t f2 = agf.find("NM F2 ");
  REQUIRE(f2 != std::string::npos);
  rtt::material::MaterialLibrary lib;
  lib.add_catalog_text(agf.substr(0, f2), "CZ");
  lib.add_catalog_text(agf, "CA");
  System s = singlet();
  Element l2 = element(s, 1);
  l2.name = "L2";
  l2.surfaces[0].id = SurfaceId("L2.S1");
  l2.surfaces[1].id = SurfaceId("L2.S2");
  l2.pose.position[2] = Param(2.0);
  l2.material = "CA:F2";
  element(s, 1).material = "CZ:N-BK7";
  s.root.children.insert(s.root.children.begin() + 2, Node{l2});
  const ZmxExport out = export_zmx(s, lib);
  REQUIRE(out.text.find("\r\nGCAT CZ CA\r\n") != std::string::npos);
  const System back = import_zmx(out.text, lib).system;
  REQUIRE(std::get<Element>(back.root.children[1].value).material == "CZ:N-BK7");
  REQUIRE(std::get<Element>(back.root.children[2].value).material == "CA:F2");
  // Both catalogues hold both glasses: no order finds both.
  rtt::material::MaterialLibrary both;
  both.add_catalog_text(agf, "CZ");
  both.add_catalog_text(agf, "CA");
  require_export_error(s, both, "export.glass_ambiguous");
}

TEST_CASE("every export error has a case", "[zmx][export]") {
  // Acceptance 6: one case per code; the registry has no export error without a case.
  std::set<std::string> seen;
  const auto check = [&](const System& s, const rtt::material::MaterialLibrary& lib,
                         const std::string& code, const ExportOptions& options = {}) {
    require_export_error(s, lib, code, options);
    seen.insert(code);
  };

  {  // not coaxial: a decentred lens
    System s = singlet();
    element(s, 1).pose.position[1] = Param(0.5);
    check(s, schott(), "export.not_coaxial");
  }
  {  // two configurations without a choice
    System s = singlet();
    s.configurations = {{"a"}, {"b"}};
    check(s, schott(), "export.multi_config");
  }
  {  // a reflection at a lens surface
    System s = singlet();
    s.paths = {{"p",
                false,
                {{SurfaceId("STO"), EventKind::Transmit, 0},
                 {SurfaceId("L1.S1"), EventKind::Reflect, 0},
                 {SurfaceId("IMG"), EventKind::Transmit, 0}}}};
    check(s, schott(), "export.event_unsupported");
  }
  {  // the stop twice (double pass)
    System s = singlet();
    s.paths = {{"p",
                false,
                {{SurfaceId("STO"), EventKind::Transmit, 0},
                 {SurfaceId("STO"), EventKind::Transmit, 0},
                 {SurfaceId("L1.S1"), EventKind::Refract, 0},
                 {SurfaceId("L1.S2"), EventKind::Refract, 0},
                 {SurfaceId("IMG"), EventKind::Transmit, 0}}}};
    check(s, schott(), "export.event_unsupported");
  }
  {  // a detector before the last event
    System s = singlet();
    Element second = element(s, 2);
    second.name = "image 2";
    second.surfaces[0].id = SurfaceId("IMG2");
    second.pose.position[2] = Param(1.0);
    s.root.children.push_back(Node{second});
    check(s, schott(), "export.event_unsupported");
  }
  {  // a mirror behind the first lens surface, inside the glass
    System s = singlet();
    Element m{"M1", ElementKind::Mirror, {}, std::nullopt, {}};
    Surface ms;
    ms.id = SurfaceId("M1.S1");
    ms.pose.position[2] = Param(8.0);  // between L1.S1 at z = 5 and L1.S2 at z = 11
    ms.interaction = rtt::model::IdealMirror{};
    m.surfaces.push_back(ms);
    s.root.children.insert(s.root.children.begin() + 2, Node{m});
    s.paths = {{"p",
                false,
                {{SurfaceId("STO"), EventKind::Transmit, 0},
                 {SurfaceId("L1.S1"), EventKind::Refract, 0},
                 {SurfaceId("M1.S1"), EventKind::Reflect, 0},
                 {SurfaceId("IMG"), EventKind::Transmit, 0}}}};
    check(s, schott(), "export.mirror_in_glass");
  }
  check(rtt::io::load_system(kRefDir + "m2/cooke_triplet.rtt.json"), schott_m2(),
        "export.element_in_glass");
  {  // glass in object space in front of the stop
    System s = singlet();
    s.object.at_infinity = false;
    s.object.distance = Param(100.0);
    s.object.medium = "SCHOTT:N-BK7";
    check(s, schott(), "export.object_medium_unsupported");
  }
  check(rtt::io::load_system(kRefDir + "m1/paraboloid_mirror.rtt.json"), schott(),
        "export.stop_missing");
  {  // no detector at the end
    System s = singlet();
    s.root.children.pop_back();
    check(s, schott(), "export.image_missing");
  }
  {  // a rectangular stop
    System s = singlet();
    element(s, 0).surfaces[0].aperture = rtt::model::RectangularAperture{5.0, 4.0};
    check(s, schott(), "export.aperture_unsupported");
  }
  {  // object_na read paraxially
    System s = singlet();
    s.object.at_infinity = false;
    s.object.distance = Param(100.0);
    s.aperture.type = rtt::model::SystemApertureType::ObjectSpaceNA;
    s.aperture.value = Param(0.05);
    check(s, schott(), "export.aperture_unsupported");
  }
  {  // environment not at 20 degC
    System s = singlet();
    s.environment.temperature_c = 25.0;
    check(s, schott(), "export.environment_unsupported");
  }
  {  // 13 fields
    System s = singlet();
    s.fields.points.assign(13, rtt::model::Field{0.0, 1.0, 1.0});
    check(s, schott(), "export.limit_exceeded");
  }
  {  // 25 wavelengths
    System s = singlet();
    s.wavelengths.assign(25, rtt::model::Wavelength{0.55, 1.0, false});
    s.wavelengths[0].reference = true;
    check(s, schott(), "export.limit_exceeded");
  }
  {  // a model glass, vacuum and a MIL number as a glass name
    System s = singlet();
    element(s, 1).material = "CONST:1.5";
    check(s, schott(), "export.material_unsupported");
    element(s, 1).material = "VACUUM";
    check(s, schott(), "export.material_unsupported");
    std::string agf = read_file(std::string(RTT_CATALOG_DIR) + "/schott.agf");
    const std::size_t at = agf.find("NM N-BK7 ");
    REQUIRE(at != std::string::npos);
    agf.replace(at, 9, "NM 517642 ");
    rtt::material::MaterialLibrary lib;
    lib.add_catalog_text(agf, "MIL");
    element(s, 1).material = "MIL:517642";
    check(s, lib, "export.material_unsupported");
  }
  {  // N-BK7 and F2 in both catalogues
    const std::string agf = read_file(std::string(RTT_CATALOG_DIR) + "/schott.agf");
    rtt::material::MaterialLibrary lib;
    lib.add_catalog_text(agf, "CA");
    lib.add_catalog_text(agf, "CB");
    System s = singlet();
    element(s, 1).material = "CA:N-BK7";
    Element l0 = element(s, 1);
    l0.name = "L0";
    l0.surfaces[0].id = SurfaceId("L0.S1");
    l0.surfaces[1].id = SurfaceId("L0.S2");
    l0.material = "CB:F2";
    l0.pose = {};
    element(s, 0).pose.reference = rtt::model::PoseReference::RelativeToPreceding;
    element(s, 0).pose.position[2] = Param(3.0);
    s.root.children.insert(s.root.children.begin(), Node{l0});
    check(s, lib, "export.glass_ambiguous");
  }
  {  // shapes without a mapping: biconic, A18, odd asphere with c != 0
    System s = singlet();
    element(s, 1).surfaces[1].shape.base =
        rtt::model::Biconic{Param(-50.0), Param(-60.0), Param(0.0), Param(0.0)};
    check(s, schott(), "export.surface_unsupported");
    System t = singlet();
    auto& a = std::get<rtt::model::EvenAsphere>(element(t, 1).surfaces[0].shape.base);
    a.coefficients.resize(8, Param(0.0));
    a.coefficients[7] = Param(1e-20);
    check(t, schott(), "export.surface_unsupported");
    System u = singlet();
    element(u, 1).surfaces[1].shape.base = rtt::model::OddAsphere{
        Param(-0.02), Param(0.0), {Param(0.0), Param(0.0), Param(0.0), Param(1e-6)}};
    check(u, schott(), "export.surface_unsupported");
  }
  {  // a phase layer, an ideal mirror at a lens surface, a polarizer
    System s = singlet();
    element(s, 1).surfaces[1].phases.emplace_back(
        rtt::model::RadialPhase{Param(8.0), {Param(1.0)}});
    check(s, schott(), "export.interaction_unsupported");
    System t = singlet();
    element(t, 1).surfaces[0].interaction = rtt::model::IdealMirror{};
    check(t, schott(), "export.interaction_unsupported");
    System u = singlet();
    Element pol{"P", ElementKind::ThinElement, {}, std::nullopt, {}};
    Surface ps;
    ps.id = SurfaceId("P.S1");
    ps.interaction = rtt::model::IdealPolarizer{};
    pol.surfaces.push_back(ps);
    pol.pose.reference = rtt::model::PoseReference::RelativeToPreceding;
    pol.pose.position[2] = Param(2.0);
    u.root.children.insert(u.root.children.begin() + 1, Node{pol});
    element(u, 2).pose.position[2] = Param(3.0);
    check(u, schott(), "export.interaction_unsupported");
  }

  std::set<std::string> registered;
  for (const auto& info : rtt::diagnostics::kCodes) {
    if (info.producer == "export" && info.severity == rtt::diagnostics::Severity::Error) {
      registered.insert(std::string(info.code));
    }
  }
  REQUIRE(seen == registered);
}

TEST_CASE("every export warning has a case", "[zmx][export]") {
  std::set<std::string> seen;
  const auto check = [&](const System& s, const std::vector<std::string>& codes) {
    const ZmxExport out = export_zmx(s, schott());
    REQUIRE(codes_of(out.warnings) == codes);
    for (const auto& d : out.warnings) {
      REQUIRE(d.severity == rtt::model::Severity::Warning);
      REQUIRE_FALSE(d.location.empty());
    }
    (void)import_zmx(out.text, schott());  // the file imports
    seen.insert(codes.begin(), codes.end());
  };

  {  // a second path
    System s = singlet();
    s.paths.push_back(s.paths.front());
    s.paths.back().name = "copy";
    check(s, {"export.paths_dropped"});
  }
  {  // a thin element that the path skips
    System s = singlet();
    Element t{"T", ElementKind::ThinElement, {}, std::nullopt, {}};
    Surface ts;
    ts.id = SurfaceId("T.S1");
    ts.aperture = rtt::model::CircularAperture{9.0, 0.0};
    t.surfaces.push_back(ts);
    t.pose.reference = rtt::model::PoseReference::RelativeToPreceding;
    t.pose.position[2] = Param(2.0);
    s.root.children.insert(s.root.children.begin() + 1, Node{t});
    element(s, 2).pose.position[2] = Param(3.0);
    s.paths = {{"main",
                false,
                {{SurfaceId("STO"), EventKind::Transmit, 0},
                 {SurfaceId("L1.S1"), EventKind::Refract, 0},
                 {SurfaceId("L1.S2"), EventKind::Refract, 0},
                 {SurfaceId("IMG"), EventKind::Transmit, 0}}}};
    check(s, {"export.surface_off_path"});
  }
  {  // an ideal lens with a finite design conjugate
    System s = singlet();
    Element t{"T", ElementKind::ThinElement, {}, std::nullopt, {}};
    Surface ts;
    ts.id = SurfaceId("T.S1");
    ts.interaction = rtt::model::IdealLens{Param(50.0), Param(1000.0)};
    t.surfaces.push_back(ts);
    t.pose = element(s, 1).pose;
    s.root.children[1] = Node{t};
    element(s, 2).pose.position[2] = Param(50.0);
    check(s, {"export.paraxial_opd"});
  }
  {  // a mirror with a substrate
    System s = singlet();
    s.root.children.erase(s.root.children.begin() + 1);
    Element m{"M1", ElementKind::Mirror, {}, std::string("SCHOTT:N-BK7"), {}};
    Surface ms;
    ms.id = SurfaceId("M1.S1");
    ms.shape.base = rtt::model::Conic{Param(-200.0), Param(-1.0)};
    m.surfaces.push_back(ms);
    m.pose.reference = rtt::model::PoseReference::RelativeToPreceding;
    m.pose.position[2] = Param(20.0);
    s.root.children.insert(s.root.children.begin() + 1, Node{m});
    element(s, 2).pose.position[2] = Param(-100.0);
    check(s, {"export.mirror_material_dropped"});
  }
  {  // a coating and an ideal anti-reflection layer: no coating library needed
    System s = singlet();
    element(s, 1).surfaces[0].interaction = rtt::model::CoatingRef{"DEMO:AR_MGF2"};
    element(s, 1).surfaces[1].interaction = rtt::model::IdealAntiReflection{};
    check(s, {"export.coating_dropped", "export.coating_dropped"});
  }
  {  // a name with a line break and a tab
    System s = singlet();
    s.name = "two\nlines\tand a tab";
    check(s, {"export.name_changed"});
    REQUIRE(export_zmx(s, schott()).text.find("\r\nNAME two lines and a tab\r\n") !=
            std::string::npos);
  }

  std::set<std::string> registered;
  for (const auto& info : rtt::diagnostics::kCodes) {
    if (info.producer == "export" && info.severity == rtt::diagnostics::Severity::Warning) {
      registered.insert(std::string(info.code));
    }
  }
  REQUIRE(seen == registered);
}

TEST_CASE("an ideal lens is a PARAXIAL surface with its OPD mode", "[zmx][export]") {
  System s = singlet();
  Element t{"T1", ElementKind::ThinElement, {}, std::nullopt, {}};
  Surface ts;
  ts.id = SurfaceId("T1.S1");
  ts.interaction = rtt::model::IdealLens{Param(50.0), std::nullopt};
  t.surfaces.push_back(ts);
  t.pose = element(s, 1).pose;
  s.root.children[1] = Node{t};
  element(s, 2).pose.position[2] = Param(50.0);
  const ZmxExport out = export_zmx(s, schott());
  REQUIRE(out.warnings.empty());
  REQUIRE(out.text.find("  TYPE PARAXIAL\r\n") != std::string::npos);
  REQUIRE(out.text.find("  PARM 1 50\r\n") != std::string::npos);
  REQUIRE(out.text.find("  PARM 2 2\r\n") != std::string::npos);
  const System back = import_zmx(out.text, schott()).system;
  const auto& lens = std::get<Element>(back.root.children[1].value);
  REQUIRE(lens.kind == ElementKind::ThinElement);
  const auto& il = std::get<rtt::model::IdealLens>(lens.surfaces[0].interaction);
  REQUIRE(il.focal_length.value == 50.0);
  REQUIRE_FALSE(il.object_distance);
}

TEST_CASE("an object in glass and the mirror keep their media", "[zmx][export]") {
  // The object medium is a GLAS line at surface 0; after a mirror the thickness changes sign.
  const System immersion = import_zmx(read_file(kZmxDir + "immersion.zmx"), schott()).system;
  REQUIRE(immersion.object.medium);
  const std::string text = export_zmx(immersion, schott()).text;
  const std::size_t s0 = text.find("SURF 0\r\n");
  const std::size_t s1 = text.find("SURF 1\r\n");
  REQUIRE(s0 != std::string::npos);
  REQUIRE(text.substr(s0, s1 - s0).find("  GLAS ") != std::string::npos);
  const System mirror = import_zmx(read_file(kZmxDir + "mirror.zmx"), schott()).system;
  const std::string m = export_zmx(mirror, schott()).text;
  REQUIRE(m.find("  GLAS MIRROR 0\r\n") != std::string::npos);
  REQUIRE(m.find("  DISZ -") != std::string::npos);
}
