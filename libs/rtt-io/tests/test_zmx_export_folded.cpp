// ZMX export of folded and tilted paths (ADR 0045, addendum #366): the non-coaxial import cases
// through import, export and import (acceptance 1) with fixed rays (acceptance 2), the sign rule
// of the thicknesses on the text (acceptance 3), the fold mirror against the file of
// make_export_cases.py (acceptance 4), folded reference systems in air (acceptance 6), random
// folded systems (acceptance 7) and the new error cases (acceptance 8).

#include <algorithm>
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <system_error>
#include <utility>
#include <variant>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/io/zmx.hpp"
#include "rtt/material/material.hpp"
#include "rtt/math/types.hpp"
#include "rtt/model/system.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/sequential.hpp"

namespace {

using rtt::compile::CompiledSystem;
using rtt::io::export_zmx;
using rtt::io::ExportError;
using rtt::io::import_zmx;
using rtt::io::ZmxExport;
using rtt::math::Mat3;
using rtt::math::Vec3;
using rtt::model::Element;
using rtt::model::ElementKind;
using rtt::model::Node;
using rtt::model::Param;
using rtt::model::Surface;
using rtt::model::SurfaceId;
using rtt::model::System;
using rtt::trace::RayBatch;
using rtt::trace::RayStatus;

const std::string kZmxDir = std::string(RTT_ZMX_DIR) + "/";
const std::string kRefDir = std::string(RTT_REFERENCE_DIR) + "/";

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
std::vector<const Surface*> surfaces_of(const System& s) {
  std::vector<const Surface*> out;
  for (const Node& n : s.root.children) {
    for (const Surface& surface : std::get<Element>(n.value).surfaces) out.push_back(&surface);
  }
  return out;
}

/// The surface as seen from a frame turned by 180 deg about x or y: c -> -c, A -> -A, the conic
/// stays (ADR 0045, point 2), for the shapes the export writes.
rtt::model::ShapeStack mirrored(rtt::model::ShapeStack shape) {
  const auto negate = [](Param& p) { p.value = -p.value + 0.0; };
  if (auto* c = std::get_if<rtt::model::Conic>(&shape.base)) {
    negate(c->radius);
  } else if (auto* a = std::get_if<rtt::model::EvenAsphere>(&shape.base)) {
    negate(a->radius);
    for (Param& x : a->coefficients) negate(x);
  } else if (auto* o = std::get_if<rtt::model::OddAsphere>(&shape.base)) {
    negate(o->curvature);
    for (Param& x : o->coefficients) negate(x);
  }
  return shape;
}

/// Radius of a conic or even asphere, none otherwise.
Param* radius_of(rtt::model::ShapeStack& s) {
  if (auto* c = std::get_if<rtt::model::Conic>(&s.base)) return &c->radius;
  if (auto* a = std::get_if<rtt::model::EvenAsphere>(&s.base)) return &a->radius;
  return nullptr;
}

/// Acceptance 1 of the addendum: the same surfaces in global coordinates. Vertices relative to
/// the first one to `tol` mm; the frames equal up to a turn by 180 deg about x, y or z (ADR 0045,
/// addendum #366, point 2) to `tol`, the shape mirrored for x and y; R to 1e-15 relative, the
/// rest of each surface (aperture, interaction, phases) bitwise. Returns the number of surfaces
/// whose frame is turned.
std::size_t require_same_surfaces(const System& a,
                                  const System& b,
                                  const rtt::material::MaterialLibrary& lib,
                                  double tol = 1e-12) {
  const CompiledSystem ca = rtt::compile::compile(a, lib);
  const CompiledSystem cb = rtt::compile::compile(b, lib);
  const auto sa = surfaces_of(a);
  const auto sb = surfaces_of(b);
  REQUIRE(sa.size() == sb.size());
  REQUIRE(ca.surfaces().size() == sa.size());
  REQUIRE(cb.surfaces().size() == sb.size());
  const Vec3 a0 = ca.surfaces().front().to_global.translation();
  const Vec3 b0 = cb.surfaces().front().to_global.translation();
  std::size_t turned = 0;
  for (std::size_t i = 0; i < sa.size(); ++i) {
    INFO("surface " << i << ": " << sa[i]->id.str() << " / " << sb[i]->id.str());
    const auto& ta = ca.surfaces()[i].to_global;
    const auto& tb = cb.surfaces()[i].to_global;
    REQUIRE(((tb.translation() - b0) - (ta.translation() - a0)).norm() <= tol);
    const Mat3 s = ta.rotation().transpose() * tb.rotation();
    std::optional<Vec3> which;
    for (const Vec3& d : {Vec3(1.0, 1.0, 1.0), Vec3(-1.0, -1.0, 1.0), Vec3(1.0, -1.0, -1.0),
                          Vec3(-1.0, 1.0, -1.0)}) {
      if ((s - Mat3(d.asDiagonal())).cwiseAbs().maxCoeff() <= tol) which = d;
    }
    INFO("R_a^T R_b =\n" << s);
    REQUIRE(which.has_value());
    const bool flip = which.value().z() < 0.0;
    if (which.value() != Vec3(1.0, 1.0, 1.0)) ++turned;
    Surface want = *sa[i];
    if (flip) want.shape = mirrored(want.shape);
    Surface got = *sb[i];
    Param* rw = radius_of(want.shape);
    Param* rg = radius_of(got.shape);
    if (rw != nullptr && rg != nullptr) {
      REQUIRE(rel(rg->value, rw->value) <= 1e-15);
      rg->value = rw->value;
    }
    got.id = want.id;
    got.pose = want.pose;
    REQUIRE(got == want);
  }
  return turned;
}

/// Traces start rays (global, unit directions) along path 0; returns the batch after it.
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

/// Acceptance 2: parallel rays along +z on a grid of radius r, 1 mm in front of the first vertex
/// (not in its plane, #151), through both systems; the same status per ray and the same hit
/// points relative to the first vertex to 1e-9 mm. Returns the hit points of the alive rays.
std::vector<Vec3> require_same_rays(const System& a,
                                    const System& b,
                                    const rtt::material::MaterialLibrary& lib,
                                    double r) {
  const CompiledSystem ca = rtt::compile::compile(a, lib);
  const CompiledSystem cb = rtt::compile::compile(b, lib);
  const Vec3 a0 = ca.surfaces().front().to_global.translation();
  const Vec3 b0 = cb.surfaces().front().to_global.translation();
  std::vector<std::pair<Vec3, Vec3>> ra;
  std::vector<std::pair<Vec3, Vec3>> rb;
  for (int iy = -2; iy <= 2; ++iy) {
    for (int ix = -2; ix <= 2; ++ix) {
      const Vec3 p(r * ix / 2.0, r * iy / 2.0, -1.0);
      ra.emplace_back(a0 + p, Vec3(0.0, 0.0, 1.0));
      rb.emplace_back(b0 + p, Vec3(0.0, 0.0, 1.0));
    }
  }
  const RayBatch ta = trace_rays(ca, ra);
  const RayBatch tb = trace_rays(cb, rb);
  std::vector<Vec3> hits;
  for (std::size_t i = 0; i < ta.size(); ++i) {
    INFO("ray " << i);
    REQUIRE(ta.status()[i] == tb.status()[i]);
    if (ta.status()[i] != RayStatus::Alive) continue;
    const Vec3 pa = Vec3(ta.pos_x()[i], ta.pos_y()[i], ta.pos_z()[i]) - a0;
    const Vec3 pb = Vec3(tb.pos_x()[i], tb.pos_y()[i], tb.pos_z()[i]) - b0;
    REQUIRE((pa - pb).norm() <= 1e-9);
    hits.push_back(pa);
  }
  return hits;
}

/// Acceptance 3, on the text alone: every thickness != 0 (DISZ of a surface or of a coordinate
/// break) has the sign (-1)^m, m the number of MIRROR surfaces up to it (UG "Thickness";
/// KA-01588). Returns the number of negative thicknesses.
std::size_t require_sign_rule(const std::string& text) {
  int mirrors = 0;
  std::size_t negative = 0;
  std::vector<std::vector<std::string>> block;
  const auto close = [&]() {
    for (const auto& w : block) {
      if (w.size() >= 2 && w[0] == "GLAS" && w[1] == "MIRROR") ++mirrors;
    }
    for (const auto& w : block) {
      if (w.size() < 2 || w[0] != "DISZ" || w[1] == "INFINITY") continue;
      const std::optional<double> t = as_number(w[1]);
      REQUIRE(t.has_value());
      if (t.value() == 0.0) continue;
      INFO("DISZ " << w[1] << " after " << mirrors << " mirrors");
      REQUIRE((t.value() > 0.0) == (mirrors % 2 == 0));
      if (t.value() < 0.0) ++negative;
    }
    block.clear();
  };
  for (const std::string& line : lines_of(text)) {
    if (line.starts_with("SURF ") || line == "MNUM 1") close();
    if (line.starts_with("  ")) block.push_back(words_of(line));
  }
  close();
  return negative;
}

struct RoundTrip {
  System a;          ///< the import of the case
  System b;          ///< the import of its export
  std::string text;  ///< the export
};

/// Import, export and import of a ZMX case; the export without warnings.
RoundTrip round_trip(const std::string& name) {
  const System a = import_zmx(read_file(kZmxDir + name + ".zmx"), schott()).system;
  const ZmxExport out = export_zmx(a, schott());
  INFO(out.text);
  REQUIRE(out.warnings.empty());
  require_sign_rule(out.text);
  REQUIRE(out.text.find("  TYPE COORDBRK\r\n") != std::string::npos);
  return {a, import_zmx(out.text, schott()).system, out.text};
}

/// The number of coordinate breaks with a GLAS line (a break inside glass, ADR 0045, addendum
/// #366, point 6). OpticStudio reads a break without GLAS as air; the import of this repository
/// takes the medium before it either way, so only the text shows it.
std::size_t breaks_in_glass(const std::string& text) {
  std::size_t n = 0;
  bool in_break = false;
  for (const std::string& line : lines_of(text)) {
    if (line.starts_with("SURF ")) in_break = false;
    if (line == "  TYPE COORDBRK") in_break = true;
    if (in_break && line.starts_with("  GLAS ")) ++n;
  }
  return n;
}

void require_export_error(const System& s, const std::string& code) {
  try {
    (void)export_zmx(s, schott());
    FAIL("no ExportError, expected " << code);
  } catch (const ExportError& e) {
    INFO(e.what());
    REQUIRE(e.code() == code);
  }
}

Element& element(System& s, std::size_t i) {
  return std::get<Element>(s.root.children.at(i).value);
}

/// Deterministic numbers in [0, 1) for every platform (std:: distributions are not, rule 7).
class Numbers {
 public:
  explicit Numbers(std::uint64_t seed) : state_(seed) {}
  double next() {
    state_ = state_ * 6364136223846793005ULL + 1442695040888963407ULL;  // Knuth's MMIX LCG
    return static_cast<double>(state_ >> 11U) * (1.0 / 9007199254740992.0);
  }
  double in(double lo, double hi) { return lo + (hi - lo) * next(); }

 private:
  std::uint64_t state_;
};

}  // namespace

TEST_CASE("the folded ZMX cases survive import, export and import", "[zmx][export][366]") {
  // Acceptances 1 and 2 of the addendum; every case of tests/zmx that is not coaxial. In some
  // cases (gimbal) the grid of parallel rays misses the tilted apertures; the status is compared.
  std::size_t turned = 0;
  std::size_t alive = 0;  // rays that arrive, over all cases
  for (const char* name : {"fold", "oap", "pivot", "tilt_image", "tilt_plate", "tilt_prism",
                           "gimbal", "gimbal_x", "random_breaks"}) {
    INFO(name);
    const RoundTrip r = round_trip(name);
    turned += require_same_surfaces(r.a, r.b, schott());
    const std::vector<Vec3> hits = require_same_rays(r.a, r.b, schott(), 0.5);
    alive += hits.size();
    // content guard: in the folds and the tilted glass every ray arrives, at different points
    const std::string n = name;
    if (n == "fold" || n == "tilt_plate" || n == "tilt_prism") {
      REQUIRE(hits.size() == 25);
      std::set<std::pair<double, double>> distinct;
      for (const Vec3& h : hits) distinct.emplace(h.x(), h.y());
      REQUIRE(distinct.size() > 1);
    }
    // the faces of the prism are tilted against each other: two breaks inside N-BK7
    if (n == "tilt_prism") REQUIRE(breaks_in_glass(r.text) == 2);
  }
  // the imported files already have the frames of OpticStudio, so nothing is turned here; the
  // turn is checked with the fold of make_export_cases.py and the random systems
  REQUIRE(turned == 0);
  REQUIRE(alive > 100);
}

TEST_CASE("the fold mirror is written as the independently built file", "[zmx][export][366]") {
  // Acceptance 4: make_export_cases.py writes export_fold.zmx after the pattern of "Add Fold
  // Mirror" (KA-01573): a break with the tilt, MIRROR, a break with the same tilt and the
  // negative thickness. Numbers bitwise, except the air wavelengths and FNUM (two Ciddor) and,
  // in the coordinate breaks, tilts to 1e-12 deg and thicknesses to 1e-12 mm (the frames come
  // from compile through sin and cos).
  const System s = rtt::io::load_system(kZmxDir + "export_fold.rtt.json");
  const ZmxExport out = export_zmx(s, schott());
  INFO(out.text);
  REQUIRE(out.warnings.empty());
  const std::vector<std::string> got = lines_of(out.text);
  const std::vector<std::string> want = lines_of(read_file(kZmxDir + "export_fold.zmx"));
  REQUIRE(got.size() == want.size());
  bool in_break = false;
  std::size_t breaks = 0;
  for (std::size_t i = 0; i < want.size(); ++i) {
    INFO("line " << (i + 1) << ": '" << got[i] << "' vs '" << want[i] << "'");
    REQUIRE(got[i].starts_with(" ") == want[i].starts_with(" "));
    const std::vector<std::string> g = words_of(got[i]);
    const std::vector<std::string> w = words_of(want[i]);
    REQUIRE(g.size() == w.size());
    if (w[0] == "SURF") in_break = false;
    if (w[0] == "TYPE" && w[1] == "COORDBRK") {
      in_break = true;
      ++breaks;
    }
    for (std::size_t k = 0; k < w.size(); ++k) {
      const std::optional<double> x = as_number(g[k]);
      const std::optional<double> y = as_number(w[k]);
      if (!y) {
        REQUIRE(g[k] == w[k]);
        continue;
      }
      REQUIRE(x);
      if ((w[0] == "WAVM" && k == 2) || (w[0] == "FNUM" && k == 1)) {
        REQUIRE(rel(*x, *y) <= 1e-15);
      } else if (in_break && k > 0 && (w[0] == "PARM" || w[0] == "DISZ")) {
        REQUIRE(std::abs(*x - *y) <= 1e-12);
      } else {
        REQUIRE(*x == *y);
      }
    }
  }
  REQUIRE(breaks == 2);
  REQUIRE(require_sign_rule(out.text) == 1);
  // and back: the same surfaces, the image frame turned by 180 deg about x
  const System back = import_zmx(out.text, schott()).system;
  REQUIRE(require_same_surfaces(s, back, schott()) == 1);
}

TEST_CASE("folded reference systems in air go through export and import", "[zmx][export][366]") {
  // Acceptance 6: m6/fold45_focus and m3/metal_periscope use VACUUM and coated mirrors; here in
  // AIR with ideal mirrors. The fold of 45 deg and the periscope of two mirrors (its path
  // "periscope" only).
  for (const char* name : {"m6/fold45_focus.rtt.json", "m3/metal_periscope.rtt.json"}) {
    INFO(name);
    System s = rtt::io::load_system(kRefDir + name);
    s.environment.medium = "AIR";
    // the periscope has a second path, single, that ends at the first mirror: keep periscope
    std::erase_if(s.paths, [](const rtt::model::Path& p) { return p.name == "single"; });
    for (Node& n : s.root.children) {
      Element& e = std::get<Element>(n.value);
      if (e.kind == ElementKind::Mirror) {
        e.material.reset();
        for (Surface& f : e.surfaces) f.interaction = rtt::model::IdealMirror{};
      }
    }
    if (!s.aperture.stop) {
      // a stop on the first surface, as the export needs one (export.stop_missing)
      s.aperture.stop = surfaces_of(s).front()->id;
      auto& first = std::get<Element>(s.root.children.front().value).surfaces.front();
      if (!first.aperture) first.aperture = rtt::model::CircularAperture{5.0, 0.0};
    }
    const ZmxExport out = export_zmx(s, schott());
    INFO(out.text);
    REQUIRE(out.text.find("TYPE COORDBRK") != std::string::npos);
    REQUIRE(require_sign_rule(out.text) >= 1);
    const System back = import_zmx(out.text, schott()).system;
    require_same_surfaces(s, back, schott());
    const std::vector<Vec3> hits = require_same_rays(s, back, schott(), 0.5);
    REQUIRE(hits.size() == 25);  // content guard: every ray arrives
  }
}

TEST_CASE("random folded systems survive export and import", "[zmx][export][366]") {
  // Acceptance 7: 200 systems with a fixed seed: a stop, 1 to 5 mirrors (plane or curved) in
  // random frames, the light folded by random angles, a detector in a random frame.
  Numbers rng(366);
  std::size_t turned = 0;
  for (int k = 0; k < 200; ++k) {
    INFO("system " << k);
    System s = rtt::io::load_system(kZmxDir + "export_singlet.rtt.json");
    s.root.children.clear();
    const auto place = [](Element& e, const Vec3& p, const Vec3& deg) {
      e.pose.position = {Param(p.x()), Param(p.y()), Param(p.z())};
      e.pose.rotation_deg = {Param(deg.x()), Param(deg.y()), Param(deg.z())};
    };
    Element stop{"stop", ElementKind::Stop, {}, std::nullopt, {}};
    Surface st;
    st.id = SurfaceId("STO");
    st.aperture = rtt::model::CircularAperture{5.0, 0.0};
    stop.surfaces.push_back(st);
    s.root.children.push_back(Node{stop});
    Vec3 v(0.0, 0.0, 0.0);
    Vec3 dir(0.0, 0.0, 1.0);
    const int mirrors = 1 + static_cast<int>(rng.next() * 5.0);
    for (int m = 0; m < mirrors; ++m) {
      v += dir * rng.in(5.0, 40.0);
      Element e{"M" + std::to_string(m + 1), ElementKind::Mirror, {}, std::nullopt, {}};
      Surface f;
      f.id = SurfaceId("M" + std::to_string(m + 1) + ".S");
      if (rng.next() < 0.5)
        f.shape.base = rtt::model::Conic{Param(rng.in(-300.0, 300.0)), Param(rng.in(-1.0, 0.0))};
      f.aperture = rtt::model::CircularAperture{10.0, 0.0};
      f.interaction = rtt::model::IdealMirror{};
      e.surfaces.push_back(f);
      place(e, v, Vec3(rng.in(-180.0, 180.0), rng.in(-80.0, 80.0), rng.in(-180.0, 180.0)));
      s.root.children.push_back(Node{e});
      // a new direction of the light, at least 10 deg away from the plane of nothing in
      // particular: any unit vector
      const double theta = rng.in(0.2, 2.9);
      const double phi = rng.in(0.0, 6.28);
      dir = Vec3(std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta));
    }
    v += dir * rng.in(5.0, 40.0);
    Element image{"image", ElementKind::Detector, {}, std::nullopt, {}};
    Surface im;
    im.id = SurfaceId("IMG");
    image.surfaces.push_back(im);
    place(image, v, Vec3(rng.in(-180.0, 180.0), rng.in(-80.0, 80.0), rng.in(-180.0, 180.0)));
    s.root.children.push_back(Node{image});
    s.aperture.stop.reset();
    ZmxExport out;
    try {
      out = export_zmx(s, schott());
    } catch (const ExportError& e) {
      // a segment in the plane of a surface or a detector seen edge-on is possible by chance
      INFO(e.what());
      REQUIRE(e.code() == "export.frame_ambiguous");
      continue;
    }
    INFO(out.text);
    REQUIRE(out.warnings.empty());
    require_sign_rule(out.text);
    const System back = import_zmx(out.text, schott()).system;
    turned += require_same_surfaces(s, back, schott(), 1e-11);
  }
  REQUIRE(turned > 0);
}

TEST_CASE("a folded path that the export cannot write", "[zmx][export][366]") {
  // Acceptance 8: export.frame_ambiguous for light in the plane of a surface; an ideal lens
  // after a mirror is export.surface_unsupported (ADR 0045, addendum #366, point 5).
  {
    System s = rtt::io::load_system(kZmxDir + "export_singlet.rtt.json");
    Element thin{"T", ElementKind::ThinElement, {}, std::nullopt, {}};
    Surface t;
    t.id = SurfaceId("T.S1");
    t.aperture = rtt::model::CircularAperture{3.0, 0.0};
    thin.surfaces.push_back(t);
    thin.pose.reference = rtt::model::PoseReference::RelativeToPreceding;
    thin.pose.position[2] = Param(10.0);
    thin.pose.rotation_deg[0] = Param(90.0);  // seen edge-on: the light runs in its plane
    s.root.children.insert(s.root.children.begin() + 2, Node{thin});
    element(s, 3).pose.position[2] = Param(35.0);
    require_export_error(s, "export.frame_ambiguous");
  }
  {
    System s = rtt::io::load_system(kZmxDir + "export_fold.rtt.json");
    Element thin{"T", ElementKind::ThinElement, {}, std::nullopt, {}};
    Surface t;
    t.id = SurfaceId("T.S1");
    t.interaction = rtt::model::IdealLens{Param(50.0), std::nullopt};
    thin.surfaces.push_back(t);
    thin.pose.reference = rtt::model::PoseReference::RelativeToPreceding;
    thin.pose.order = rtt::model::PoseOrder::RotateFirst;
    thin.pose.rotation_deg[0] = Param(-135.0);
    thin.pose.position[2] = Param(10.0);
    s.root.children.insert(s.root.children.begin() + 3, Node{thin});
    Element& image = element(s, 4);
    image.pose.order = rtt::model::PoseOrder::TranslateFirst;
    image.pose.rotation_deg[0] = Param(0.0);
    image.pose.position[2] = Param(15.0);
    require_export_error(s, "export.surface_unsupported");
  }
}
