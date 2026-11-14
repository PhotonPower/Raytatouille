// The ZMX/CSV export with inactive nodes (#380, ADR 0049, point 6): inactive elements are left
// out with export.inactive_skipped; thicknesses and the object distance run to the next active
// vertex, so the files equal those of the system without the element.

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <sstream>
#include <string>
#include <vector>

#include "rtt/io/json_io.hpp"
#include "rtt/io/surface_table.hpp"
#include "rtt/io/zmx.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/system.hpp"
#include "test_support.hpp"

using rtt::model::CircularAperture;
using rtt::model::Diagnostic;
using rtt::model::Element;
using rtt::model::ElementKind;
using rtt::model::Node;
using rtt::model::Pose;
using rtt::model::Surface;
using rtt::model::SurfaceId;
using rtt::model::System;

namespace {

const rtt::material::MaterialLibrary& library() {
  static rtt::material::MaterialLibrary lib;
  static const bool loaded = [] {
    lib.add_catalog(std::string(RTT_CATALOG_DIR) + "/schott.agf");
    return true;
  }();
  static_cast<void>(loaded);
  return lib;
}

/// An inactive plane window "W" of N-BK7, 2 mm thick, its first surface at global z = `z`.
Element inactive_window(double z) {
  Surface a;
  a.id = SurfaceId("W.S1");
  a.aperture = CircularAperture{12.7, 0.0};
  Surface b = a;
  b.id = SurfaceId("W.S2");
  b.pose = Pose::along_z(2.0);
  Element w{"W", ElementKind::Plate, Pose::along_z(z), "SCHOTT:N-BK7", {a, b}};
  w.active = false;
  return w;
}

std::size_t count(const std::vector<Diagnostic>& d,
                  const std::string& code,
                  const std::string& location = {}) {
  return static_cast<std::size_t>(std::count_if(d.begin(), d.end(), [&](const Diagnostic& x) {
    return x.code == code && (location.empty() || x.location == location);
  }));
}

/// The values of the DISZ lines of a ZMX text, in order (the thickness after each surface,
/// the object distance first).
std::vector<double> disz(const std::string& text) {
  std::vector<double> out;
  std::istringstream in(text);
  std::string line;
  while (std::getline(in, line)) {
    const std::size_t at = line.find("DISZ ");
    if (at == std::string::npos) continue;
    const std::string value = line.substr(at + 5);
    out.push_back(value.starts_with("INFINITY") ? -1.0 : std::stod(value));
  }
  return out;
}

/// The canonical export singlet (ADR 0045, point 6): stop at z = 0, L1 relative (5, then 6),
/// image relative 45; object at infinity.
System canonical_singlet() {
  return rtt::io::load_system(std::string(RTT_ZMX_DIR) + "/export_singlet.rtt.json");
}

}  // namespace

TEST_CASE(
    "export: a canonical system with an inactive element sums the t of the skipped "
    "surfaces",
    "[io][0.10][export]") {
  // Review of PR 5: row k is no longer tree surface k. With L1 inactive the stop is followed by
  // the image at 5 + 6 + 45 = 56 mm (exact in binary).
  System s = canonical_singlet();
  const std::vector<double> on = disz(rtt::io::export_zmx(s, library()).text);
  REQUIRE(on.size() == 5);  // object, STO, L1.S1, L1.S2, IMG
  CHECK(on[1] == 5.0);      // content guard: the active system as in its file
  rtt::model::test::element(s, 1).active = false;
  const rtt::io::ZmxExport off = rtt::io::export_zmx(s, library());
  const std::vector<double> d = disz(off.text);
  REQUIRE(d.size() == 3);  // object, STO, IMG
  CHECK(d[1] == 56.0);
  CHECK(count(off.warnings, "export.inactive_skipped", "/root/children/1") == 1);

  // Thicknesses that are not exact in binary stay the exact t of the model where nothing is
  // skipped, and are the sum of the t where something is: L1.S2 at 0.1 after L1.S1, an inactive
  // thin element 0.2 after L1.S2, the image 0.7 after it.
  System r = canonical_singlet();
  rtt::model::test::element(r, 1).surfaces[1].pose.position[2] = rtt::model::Param(0.1);
  Surface ts;
  ts.id = SurfaceId("T.S1");
  Element t{"T", ElementKind::ThinElement, Pose::along_z(0.2), std::nullopt, {ts}};
  t.pose.reference = rtt::model::PoseReference::RelativeToPreceding;
  t.active = false;
  r.root.children.insert(r.root.children.begin() + 2, Node{t});
  rtt::model::test::element(r, 3).pose.position[2] = rtt::model::Param(0.7);
  const std::vector<double> e = disz(rtt::io::export_zmx(r, library()).text);
  REQUIRE(e.size() == 5);  // object, STO, L1.S1, L1.S2, IMG
  volatile double a = 0.2;
  volatile double b = 0.7;
  CHECK(e[2] == 0.1);
  CHECK(e[3] == a + b);
}

TEST_CASE(
    "export: a canonical system with an inactive first element adds its t to the object "
    "distance",
    "[io][0.10][export]") {
  // An inactive thin element at z = 0 in front of the stop, which follows 3 mm after it; the
  // object 100 mm before the first element. The exported object distance runs to the stop.
  System s = canonical_singlet();
  s.object.at_infinity = false;
  s.object.distance = rtt::model::Param(100.0);
  s.fields = {rtt::model::FieldType::ObjectHeight, {{0.0, 0.0, 1.0}, {0.0, 2.0, 1.0}}};
  Surface ws;
  ws.id = SurfaceId("W.S1");
  Element w{"W", ElementKind::ThinElement, Pose::along_z(0.0), std::nullopt, {ws}};
  w.active = false;
  rtt::model::test::element(s, 0).pose = Pose::along_z(3.0);
  rtt::model::test::element(s, 0).pose.reference = rtt::model::PoseReference::RelativeToPreceding;
  s.root.children.insert(s.root.children.begin(), Node{w});
  const rtt::io::ZmxExport out = rtt::io::export_zmx(s, library());
  const std::vector<double> d = disz(out.text);
  REQUIRE(d.size() == 5);  // object, STO, L1.S1, L1.S2, IMG
  CHECK(d[0] == 103.0);
  CHECK(d[1] == 5.0);
  CHECK(count(out.warnings, "export.inactive_skipped", "/root/children/0") == 1);
}

TEST_CASE("export: an inactive element is left out, the files equal those without it",
          "[io][0.10][export]") {
  // make_singlet places every node absolutely, so the thicknesses of both systems are
  // differences of the compiled z and agree bit for bit.
  const System without = rtt::model::test::make_singlet();
  System with = without;
  with.root.children.insert(with.root.children.begin() + 2, Node{inactive_window(50.0)});

  const rtt::io::ZmxExport a = rtt::io::export_zmx(with, library());
  const rtt::io::ZmxExport b = rtt::io::export_zmx(without, library());
  CHECK(a.text == b.text);
  CHECK(a.text.find("W.S1") == std::string::npos);
  CHECK(count(a.warnings, "export.inactive_skipped", "/root/children/2") == 1);
  CHECK(count(a.warnings, "export.inactive_skipped") == 1);
  CHECK(count(a.warnings, "export.surface_off_path") == 0);
  CHECK(count(b.warnings, "export.inactive_skipped") == 0);

  const rtt::io::TableExport ta = rtt::io::write_surface_table(with, library());
  const rtt::io::TableExport tb = rtt::io::write_surface_table(without, library());
  CHECK(ta.text == tb.text);
  CHECK(count(ta.warnings, "export.inactive_skipped", "/root/children/2") == 1);
}

TEST_CASE("export: an inactive first element, the object distance runs to the first active one",
          "[io][0.10][export]") {
  System without = rtt::model::test::make_singlet();
  without.object.at_infinity = false;
  without.object.distance = rtt::model::Param(100.0);
  without.fields = {rtt::model::FieldType::ObjectHeight, {{0.0, 0.0, 1.0}, {0.0, 2.0, 1.0}}};
  System with = without;
  with.root.children.insert(with.root.children.begin(), Node{inactive_window(-10.0)});

  const rtt::io::ZmxExport a = rtt::io::export_zmx(with, library());
  const rtt::io::ZmxExport b = rtt::io::export_zmx(without, library());
  CHECK(a.text == b.text);
  CHECK(count(a.warnings, "export.inactive_skipped", "/root/children/0") == 1);
  // Content guard: the object distance of 100 mm is in the file.
  CHECK(b.text.find("DISZ 100") != std::string::npos);
}
