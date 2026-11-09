// CSV surface table (ADR 0045, point 7; PR C): the singlet against the table of
// make_export_cases.py, written independently of this code (acceptance 4), the Cooke triplet
// through the table (acceptance 1), the round trips and the idempotence (acceptance 5), the
// reading cases of RFC 4180, one case per table code with its position (acceptance 6) and every
// reference system (acceptance 3).

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <utility>
#include <variant>
#include <vector>

#include "rtt/analysis/spot.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/diagnostics/codes.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/io/surface_table.hpp"
#include "rtt/io/zmx.hpp"
#include "rtt/material/material.hpp"
#include "rtt/math/types.hpp"
#include "rtt/model/system.hpp"
#include "rtt/paraxial/paraxial.hpp"

namespace {

namespace fs = std::filesystem;
using rtt::io::ExportError;
using rtt::io::read_surface_table;
using rtt::io::TableError;
using rtt::io::TableExport;
using rtt::io::TableImport;
using rtt::io::write_surface_table;
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

double rel(double a, double b) {
  return std::abs(a - b) / std::abs(b);
}

System singlet() {
  return rtt::io::load_system(kZmxDir + "export_singlet.rtt.json");
}

Element& element(System& s, std::size_t i) {
  return std::get<Element>(s.root.children.at(i).value);
}

std::vector<std::string> codes_of(const std::vector<rtt::model::Diagnostic>& d) {
  std::vector<std::string> out;
  for (const auto& x : d) out.push_back(x.code);
  return out;
}

/// Records of a CSV text without quotes (the written tables quote nothing here).
std::vector<std::vector<std::string>> records(const std::string& text) {
  std::vector<std::vector<std::string>> out;
  std::size_t pos = 0;
  while (pos < text.size()) {
    const std::size_t end = text.find("\r\n", pos);
    REQUIRE(end != std::string::npos);  // CRLF after every record, also the last one
    std::vector<std::string> fields;
    std::string line = text.substr(pos, end - pos);
    std::size_t from = 0;
    while (true) {
      const std::size_t comma = line.find(',', from);
      fields.push_back(
          line.substr(from, comma == std::string::npos ? std::string::npos : comma - from));
      if (comma == std::string::npos) break;
      from = comma + 1;
    }
    out.push_back(fields);
    pos = end + 2;
  }
  return out;
}

std::optional<double> as_number(const std::string& w) {
  double v = 0.0;
  const auto [ptr, ec] = std::from_chars(w.data(), w.data() + w.size(), v);
  if (w.empty() || ec != std::errc() || ptr != w.data() + w.size()) return std::nullopt;
  return v;
}

std::vector<rtt::math::Vec3> vertices(const rtt::compile::CompiledSystem& cs) {
  std::vector<rtt::math::Vec3> out;
  const rtt::math::Vec3 first = cs.surfaces().front().to_global.apply_point({0.0, 0.0, 0.0});
  for (const auto& s : cs.surfaces()) {
    out.emplace_back(s.to_global.apply_point({0.0, 0.0, 0.0}) - first);
  }
  return out;
}

/// The table of the singlet with `field` of the record that starts with `first` replaced.
std::string edited(const std::string& csv,
                   const std::string& first,
                   std::size_t field,
                   const std::string& value) {
  std::string out;
  std::size_t pos = 0;
  while (pos < csv.size()) {
    const std::size_t end = csv.find("\r\n", pos) + 2;
    std::string line = csv.substr(pos, end - pos - 2);
    if (line.starts_with(first)) {
      std::vector<std::string> f;
      std::stringstream in(line);
      for (std::string x; std::getline(in, x, ',');) f.push_back(x);
      while (f.size() < 22) f.emplace_back();
      f.at(field - 1) = value;
      line.clear();
      for (std::size_t i = 0; i < f.size(); ++i) line += (i == 0 ? "" : ",") + f[i];
    }
    out += line + "\r\n";
    pos = end;
  }
  return out;
}

/// Reading `csv` must throw TableError with `code` at `line` and `column`.
void require_table_error(const std::string& csv,
                         const std::string& code,
                         std::size_t line,
                         std::size_t column,
                         const System& base = singlet()) {
  try {
    (void)read_surface_table(csv, base, schott());
    FAIL("no TableError, expected " << code);
  } catch (const TableError& e) {
    INFO(e.what());
    REQUIRE(e.code() == code);
    REQUIRE(e.line() == line);
    REQUIRE(e.column() == column);
  }
}

// The coaxial import cases without an odd asphere (the asphere cases have one at c = 0, which
// the table does not hold, G1).
const std::vector<std::string> kCoaxialCases{"singlet",  "singlet_in", "singlet_m", "doublet",
                                             "mirror",   "finite",     "sasian",    "immersion",
                                             "paraxial", "apertures"};

}  // namespace

TEST_CASE("the singlet is written as the independently built table", "[table]") {
  // Acceptance 4: make_export_cases.py builds export_singlet.csv from the prescription. Field by
  // field; numbers as doubles, bitwise (the table holds no air quantities).
  const TableExport out = write_surface_table(singlet(), schott());
  REQUIRE(out.warnings.empty());
  const auto got = records(out.text);
  const auto want = records(read_file(kZmxDir + "export_singlet.csv"));
  REQUIRE(got.size() == want.size());
  std::size_t numbers = 0;
  for (std::size_t i = 0; i < want.size(); ++i) {
    INFO("record " << (i + 1));
    REQUIRE(got[i].size() == want[i].size());
    for (std::size_t k = 0; k < want[i].size(); ++k) {
      INFO("field " << (k + 1) << ": '" << got[i][k] << "' vs '" << want[i][k] << "'");
      const std::optional<double> b = as_number(want[i][k]);
      if (b) {
        const std::optional<double> a = as_number(got[i][k]);
        REQUIRE(a);
        REQUIRE(*a == *b);
        ++numbers;
      } else {
        REQUIRE(got[i][k] == want[i][k]);
      }
    }
  }
  REQUIRE(numbers > 15);
}

TEST_CASE("the singlet comes back bitwise and the table is idempotent", "[table]") {
  // Acceptance 5: the singlet is canonical (ADR 0045, point 6) and writes without a warning.
  const System s = singlet();
  const TableExport out = write_surface_table(s, schott());
  const TableImport back = read_surface_table(out.text, s, schott());
  REQUIRE(back.warnings.empty());
  INFO(rtt::io::to_json(back.system));
  REQUIRE(back.system == s);
  REQUIRE(write_surface_table(back.system, schott()).text == out.text);
}

TEST_CASE("the Cooke triplet through the table keeps first order, spots and geometry", "[table]") {
  // Acceptance 1 and 8 over the table; the reference has absolute positions (not canonical).
  const System s = rtt::io::load_system(kRefDir + "m2/cooke_triplet_surface_stop.rtt.json");
  const TableExport out = write_surface_table(s, schott_m2());
  REQUIRE(out.warnings.empty());
  const TableImport back = read_surface_table(out.text, s, schott_m2());
  REQUIRE(back.warnings.empty());
  const auto ca = rtt::compile::compile(s, schott_m2());
  const auto cb = rtt::compile::compile(back.system, schott_m2());
  const auto fa = rtt::paraxial::first_order(ca, {0}, 1);
  const auto fb = rtt::paraxial::first_order(cb, {0}, 1);
  REQUIRE(fa.efl);
  REQUIRE(fb.efl);
  REQUIRE(rel(*fb.efl, *fa.efl) <= 1e-12);
  REQUIRE(rel(*fb.bfl, *fa.bfl) <= 1e-12);
  std::set<double> rms;
  for (std::uint16_t field = 0; field < s.fields.points.size(); ++field) {
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
  // the wavelengths stay in the JSON: bitwise
  REQUIRE(back.system.wavelengths == s.wavelengths);
  REQUIRE(write_surface_table(back.system, schott_m2()).text == out.text);
}

TEST_CASE("the coaxial ZMX cases come back through the table", "[table]") {
  // The import gives canonical systems (with z1 != 0 if dummies were dropped): bitwise where the
  // first element is at z = 0, otherwise vertices to 1e-12 mm; idempotent in both cases.
  for (const std::string& name : kCoaxialCases) {
    INFO(name);
    const auto& lib = name == "sasian" ? schott_m2() : schott();
    const System a = rtt::io::import_zmx(read_file(kZmxDir + name + ".zmx"), lib).system;
    const TableExport out = write_surface_table(a, lib);
    INFO(out.text);
    REQUIRE(out.warnings.empty());
    const TableImport back = read_surface_table(out.text, a, lib);
    REQUIRE(back.warnings.empty());
    REQUIRE(write_surface_table(back.system, lib).text == out.text);
    const bool canonical =
        std::get<Element>(a.root.children.front().value).pose.position[2].value == 0.0;
    if (canonical) {
      INFO("a:\n" << rtt::io::to_json(a) << "\nb:\n" << rtt::io::to_json(back.system));
      REQUIRE(back.system == a);
    } else {
      const auto va = vertices(rtt::compile::compile(a, lib));
      const auto vb = vertices(rtt::compile::compile(back.system, lib));
      REQUIRE(va.size() == vb.size());
      for (std::size_t i = 0; i < va.size(); ++i) REQUIRE((va[i] - vb[i]).norm() <= 1e-12);
    }
  }
}

TEST_CASE("the table reads quotes, LF, a byte order mark and columns in any order", "[table]") {
  const System s = singlet();
  const std::string csv = write_surface_table(s, schott()).text;
  const System expected = read_surface_table(csv, s, schott()).system;

  // an empty line at the end, as text editors add one
  REQUIRE(read_surface_table(csv + "\r\n", s, schott()).system == expected);

  // LF only, and no line end after the last record (RFC 4180, rule 2)
  std::string lf;
  for (const char c : csv) {
    if (c != '\r') lf += c;
  }
  lf.pop_back();
  REQUIRE(read_surface_table(lf, s, schott()).system == expected);

  // UTF-8 byte order mark
  REQUIRE(read_surface_table("\xEF\xBB\xBF" + csv, s, schott()).system == expected);

  // every field quoted, a quote doubled inside a name, a comma and a line break in a name
  std::string quoted = edited(csv, "surface,L1.S1", 3, "\"L1, \"\"front\"\"\r\nlens\"");
  quoted = edited(quoted, "surface,L1.S2", 3, "\"L1, \"\"front\"\"\r\nlens\"");
  const System q = read_surface_table(quoted, s, schott()).system;
  REQUIRE(std::get<Element>(q.root.children[1].value).name == "L1, \"front\"\r\nlens");
  // written back, the name is quoted again (RFC 4180, rules 6 and 7)
  const std::string again = write_surface_table(q, schott()).text;
  REQUIRE(again.find(",\"L1, \"\"front\"\"\r\nlens\",") != std::string::npos);
  REQUIRE(read_surface_table(again, s, schott()).system == q);

  // the columns in another order: swap "radius" and "thickness" in every record
  std::string swapped;
  std::size_t pos = 0;
  while (pos < csv.size()) {
    const std::size_t end = csv.find("\r\n", pos);
    std::vector<std::string> f;
    std::stringstream in(csv.substr(pos, end - pos));
    for (std::string x; std::getline(in, x, ',');) f.push_back(x);
    while (f.size() < 22) f.emplace_back();
    std::swap(f[5], f[14]);
    for (std::size_t i = 0; i < f.size(); ++i) swapped += (i == 0 ? "" : ",") + f[i];
    swapped += "\r\n";
    pos = end + 2;
  }
  REQUIRE(swapped != csv);
  REQUIRE(read_surface_table(swapped, s, schott()).system == expected);
}

TEST_CASE("the table keeps what the ZMX file cannot hold", "[table]") {
  // CONST and VACUUM, any environment, a design conjugate and a thin element without aperture.
  System s = singlet();
  s.environment.medium = "VACUUM";
  element(s, 1).material = "CONST:1.5";
  Element t{"T", ElementKind::ThinElement, {}, std::nullopt, {}};
  Surface ts;
  ts.id = SurfaceId("T.S1");
  ts.interaction = rtt::model::IdealLens{Param(500.0), Param(1000.0)};
  t.surfaces.push_back(ts);
  t.pose.reference = rtt::model::PoseReference::RelativeToPreceding;
  t.pose.position[2] = Param(2.0);
  s.root.children.insert(s.root.children.begin() + 1, Node{t});
  element(s, 2).pose.position[2] = Param(3.0);
  const TableExport out = write_surface_table(s, rtt::material::MaterialLibrary{});
  REQUIRE(out.warnings.empty());
  REQUIRE(out.text.find(",CONST:1.5,") != std::string::npos);
  const TableImport back = read_surface_table(out.text, s, rtt::material::MaterialLibrary{});
  REQUIRE(back.system == s);
}

TEST_CASE("an air gap inside a cemented lens is written and read back", "[table]") {
  // The segment of the environment medium needs its name in the table (review of PR C).
  System s = singlet();
  Element& lens = element(s, 1);
  Surface third = lens.surfaces[1];
  third.id = SurfaceId("L1.S3");
  Surface fourth = lens.surfaces[1];
  fourth.id = SurfaceId("L1.S4");
  lens.surfaces[1].shape.base = rtt::model::Plane{};
  third.shape.base = rtt::model::Plane{};
  third.pose.position[2] = Param(0.5);
  fourth.pose.position[2] = Param(4.0);
  lens.surfaces.push_back(third);
  lens.surfaces.push_back(fourth);
  lens.material.reset();
  lens.segment_materials = {"SCHOTT:N-BK7", "AIR", "SCHOTT:N-BK7"};
  const TableExport out = write_surface_table(s, schott());
  REQUIRE(out.warnings.empty());
  REQUIRE(out.text.find(",AIR,") != std::string::npos);
  const TableImport back = read_surface_table(out.text, s, schott());
  REQUIRE(back.system == s);
}

TEST_CASE("the table refuses what it cannot hold", "[table]") {
  const auto refuse = [](const System& s, const std::string& code) {
    try {
      (void)write_surface_table(s, schott());
      FAIL("no ExportError, expected " << code);
    } catch (const ExportError& e) {
      INFO(e.what());
      REQUIRE(e.code() == code);
    }
  };
  {  // G1: an odd asphere, even with c = 0
    System s = singlet();
    element(s, 1).surfaces[1].shape.base = rtt::model::OddAsphere{
        Param(0.0), Param(0.0), {Param(0.0), Param(0.0), Param(0.0), Param(1e-6)}};
    refuse(s, "export.surface_unsupported");
  }
  for (const char* name : {"asphere", "asphere_in", "asphere_cm"}) {
    INFO(name);
    refuse(rtt::io::import_zmx(read_file(kZmxDir + name + ".zmx"), schott()).system,
           "export.surface_unsupported");
  }
  {  // A18
    System s = singlet();
    auto& a = std::get<rtt::model::EvenAsphere>(element(s, 1).surfaces[0].shape.base);
    a.coefficients.resize(8, Param(0.0));
    a.coefficients[7] = Param(1e-20);
    refuse(s, "export.surface_unsupported");
  }
  {  // G2: the path leaves L1 and comes back to it; L1 of AIR, so no other rule fires first
    System s = singlet();
    element(s, 1).material = "AIR";
    Element l2 = element(s, 1);
    l2.name = "L2";
    l2.surfaces[0].id = SurfaceId("L2.S1");
    l2.surfaces[1].id = SurfaceId("L2.S2");
    l2.pose.position[2] = Param(2.0);
    s.root.children.insert(s.root.children.begin() + 2, Node{l2});
    element(s, 3).pose.position[2] = Param(40.0);
    s.paths = {{"p",
                false,
                {{SurfaceId("STO"), EventKind::Transmit, 0},
                 {SurfaceId("L1.S1"), EventKind::Refract, 0},
                 {SurfaceId("L2.S1"), EventKind::Refract, 0},
                 {SurfaceId("L2.S2"), EventKind::Refract, 0},
                 {SurfaceId("L1.S2"), EventKind::Refract, 0},
                 {SurfaceId("IMG"), EventKind::Transmit, 0}}}};
    refuse(s, "export.event_unsupported");
  }
  {  // the shared rules of points 2 and 3
    System s = singlet();
    element(s, 1).pose.position[0] = Param(0.5);
    refuse(s, "export.not_coaxial");
  }
}

TEST_CASE("every table error has a case at its line and column", "[table]") {
  // Acceptance 6 with G3: line of the file where the field starts (header = 1), field number.
  const std::string csv = write_surface_table(singlet(), schott()).text;
  std::set<std::string> seen;
  const auto check = [&](const std::string& text, const std::string& code, std::size_t line,
                         std::size_t column) {
    require_table_error(text, code, line, column);
    seen.insert(code);
  };

  // table.parse
  std::string semicolons = csv;
  std::replace(semicolons.begin(), semicolons.end(), ',', ';');
  check(semicolons, "table.parse", 1, 1);                                // Excel, German settings
  check(edited(csv, "surface,L1.S1", 3, "\"L1"), "table.parse", 4, 3);   // quote not closed
  check(edited(csv, "surface,L1.S1", 3, "L\"1"), "table.parse", 4, 3);   // quote inside a field
  check(edited(csv, "surface,L1.S1", 3, "L1,x"), "table.parse", 4, 23);  // 23 fields
  // after a field with a line break, the extra field starts on the next line (G3)
  check(edited(csv, "surface,L1.S1", 3, "\"L\r\n1\",x"), "table.parse", 5, 23);
  // and a later error counts the lines of that field
  check(edited(edited(csv, "surface,L1.S1", 3, "\"L\r\n1\""), "surface,L1.S2", 6, "fifty"),
        "table.value", 6, 6);
  check(edited(csv, "surface,L1.S1", 3, "L\xFF"), "table.parse", 4, 3);  // not UTF-8

  // table.column
  check(edited(csv, "row,", 6, "curvature"), "table.column", 1, 6);   // unknown
  check(edited(csv, "row,", 6, "thickness"), "table.column", 1, 15);  // twice
  {
    std::string missing;  // drop the column "stop" from every record
    std::size_t pos = 0;
    while (pos < csv.size()) {
      const std::size_t end = csv.find("\r\n", pos);
      std::vector<std::string> f;
      std::stringstream in(csv.substr(pos, end - pos));
      for (std::string x; std::getline(in, x, ',');) f.push_back(x);
      while (f.size() < 22) f.emplace_back();
      f.erase(f.begin() + 19);
      for (std::size_t i = 0; i < f.size(); ++i) missing += (i == 0 ? "" : ",") + f[i];
      missing += "\r\n";
      pos = end + 2;
    }
    check(missing, "table.column", 0, 0);
  }

  // table.value
  check(edited(csv, "surface,L1.S1", 6, "fifty"), "table.value", 4, 6);
  check(edited(csv, "surface,L1.S1", 6, " 50"), "table.value", 4, 6);  // a blank is part of it
  check(edited(csv, "surface,L1.S1", 6, "Infinity"), "table.value", 4, 6);
  check(edited(csv, "surface,L1.S1", 4, "doublet"), "table.value", 4, 4);
  check(edited(csv, "surface,L1.S1", 5, "toroid"), "table.value", 4, 5);
  check(edited(csv, "surface,STO", 6, "20"), "table.value", 3, 6);        // a radius at a plane
  check(edited(csv, "surface,STO", 20, "yes"), "table.value", 3, 20);     // G5
  check(edited(csv, "surface,L1.S2", 20, "true"), "table.value", 5, 20);  // a second stop
  check(edited(csv, "object", 3, "L0"), "table.value", 2, 3);       // G4: only thickness, material
  check(edited(csv, "object", 1, "surface"), "table.value", 2, 1);  // the object row first
  check(edited(csv, "surface,STO", 1, "object"), "table.value", 3, 1);
  check(edited(csv, "surface,L1.S2", 16, "SCHOTT:N-BK7"), "table.value", 5, 16);  // lens end
  check(edited(csv, "surface,L1.S1", 21, "100"), "table.value", 4, 21);           // f at a lens
  check(edited(csv, "surface,L1.S1", 17, "square"), "table.value", 4, 17);

  // table.duplicate
  check(edited(csv, "surface,IMG", 3, "stop"), "table.duplicate", 6, 3);  // element again
  check(edited(csv, "surface,L1.S2", 2, "L1.S1"), "table.duplicate", 5, 2);

  // table.invalid_result: the reading is fine, validate or compile are not
  check(edited(csv, "surface,L1.S1", 16, "SCHOTT:NOPE"), "table.invalid_result", 0, 0);

  std::set<std::string> registered;
  for (const auto& info : rtt::diagnostics::kCodes) {
    if (info.producer == "table" && info.severity == rtt::diagnostics::Severity::Error) {
      registered.insert(std::string(info.code));
    }
  }
  REQUIRE(seen == registered);
}

TEST_CASE("the blocks of the base system that refer to its tree are dropped", "[table]") {
  System base = singlet();
  const std::string csv = write_surface_table(base, schott()).text;
  base.paths.push_back(base.paths.front());
  base.paths.back().name = "second";
  base.configurations = {{"a"}, {"b"}};
  base.parameters.push_back(
      {"epd", std::vector<double>{10.0, 12.0}, false, std::nullopt, std::nullopt});
  base.aperture.type = rtt::model::SystemApertureType::EntrancePupilDiameter;
  base.aperture.value = Param::bound("epd");
  rtt::model::FirstOrderOperand efl;
  efl.common.target = 50.0;
  efl.path = "main";
  base.optimization.operands.emplace_back(efl);
  base.tolerances.tolerance_class = rtt::model::ToleranceClass::Normal;
  const TableImport r = read_surface_table(csv, base, schott());
  REQUIRE(codes_of(r.warnings) == std::vector<std::string>(5, "table.base_dropped"));
  std::set<std::string> where;
  for (const auto& d : r.warnings) where.insert(d.location);
  REQUIRE(where == std::set<std::string>{"/paths", "/optimization", "/tolerances", "/parameters",
                                         "/configurations"});
  REQUIRE(r.system.paths.size() == 1);
  REQUIRE(r.system.paths.front().automatic);
  REQUIRE(r.system.configurations.empty());
  REQUIRE(r.system.parameters.empty());
  // the bound EPD keeps its value of configuration 0
  REQUIRE(r.system.aperture.value == Param(10.0));
  bool registered = false;
  for (const auto& info : rtt::diagnostics::kCodes) {
    if (info.code == "table.base_dropped") registered = info.producer == "table";
  }
  REQUIRE(registered);
}

TEST_CASE("every reference system writes a table that reads back, or fails with its code",
          "[table]") {
  // Acceptance 3 for the table: an export.* code, or a table that reads back without an error
  // and writes the same text again. First order and spots of a path through the export are
  // compared in test_zmx_export.cpp; the table has the same geometry.
  std::size_t ok = 0;
  for (const auto& entry : fs::recursive_directory_iterator(RTT_REFERENCE_DIR)) {
    if (!entry.is_regular_file() || !entry.path().filename().string().ends_with(".rtt.json")) {
      continue;
    }
    INFO(entry.path().string());
    const System s = rtt::io::load_system(entry.path().string());
    for (const rtt::material::MaterialLibrary* lib : {&schott(), &schott_m2()}) {
      TableExport out;
      try {
        out = write_surface_table(s, *lib);
      } catch (const ExportError& e) {
        REQUIRE(e.code().starts_with("export."));
        break;
      } catch (const rtt::compile::CompileError&) {
        continue;
      }
      ++ok;
      const TableImport back = read_surface_table(out.text, s, *lib);
      REQUIRE(write_surface_table(back.system, *lib).text == out.text);
      break;
    }
  }
  REQUIRE(ok >= 8);
}
