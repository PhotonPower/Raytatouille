// Nodes taken out of the computation (#380, ADR 0049, points 2 and 3): the geometry stays, the
// paths skip inactive elements, an inactive stop counts as no stop, and the consequences of
// switching a part off are warnings.

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <variant>
#include <vector>

#include "rtt/coating/catalog.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/compile/errors.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/validate.hpp"
#include "test_support.hpp"

using rtt::compile::CompiledSystem;
using rtt::compile::CompileError;
using rtt::compile::PathId;
using rtt::material::MaterialLibrary;
using rtt::model::Assembly;
using rtt::model::CircularAperture;
using rtt::model::CoatingRef;
using rtt::model::Diagnostic;
using rtt::model::Element;
using rtt::model::ElementKind;
using rtt::model::EventKind;
using rtt::model::Node;
using rtt::model::Pose;
using rtt::model::Surface;
using rtt::model::SurfaceId;
using rtt::model::System;
using rtt::model::test::element;

namespace {

/// The reference singlet (stop /root/children/0, lens L1 /root/children/1, detector IMG
/// /root/children/2) in CONST:1.5168, so that it compiles without catalogue.
System singlet() {
  System s = rtt::model::test::make_singlet();
  element(s, 1).material = "CONST:1.5168";
  return s;
}

CompiledSystem compiled(const System& s) {
  const MaterialLibrary lib;
  return rtt::compile::compile(s, lib);
}

/// Surface ids of the events of path `p`.
std::vector<std::string> event_ids(const CompiledSystem& cs, std::uint32_t p = 0) {
  std::vector<std::string> out;
  for (const auto& e : cs.path(PathId{p}).events) out.push_back(cs.surfaces()[e.surface].id.str());
  return out;
}

/// Number of diagnostics with `code` at `location`.
std::size_t count(const std::vector<Diagnostic>& d,
                  const std::string& code,
                  const std::string& location) {
  return static_cast<std::size_t>(std::count_if(d.begin(), d.end(), [&](const Diagnostic& x) {
    return x.code == code && x.location == location;
  }));
}

std::size_t count_code(const std::vector<Diagnostic>& d, const std::string& code) {
  return static_cast<std::size_t>(
      std::count_if(d.begin(), d.end(), [&](const Diagnostic& x) { return x.code == code; }));
}

/// A flat plate (window) of two planes 2 mm apart at z = 1 mm, absolute, in CONST:1.5.
Element window() {
  Surface a;
  a.id = SurfaceId("W.S1");
  a.aperture = CircularAperture{12.7, 0.0};
  Surface b = a;
  b.id = SurfaceId("W.S2");
  b.pose = Pose::along_z(2.0);
  return Element{"window", ElementKind::Plate, Pose::along_z(1.0), "CONST:1.5", {a, b}};
}

/// A second stop element "stop2" with surface STO2 at z = 2 mm.
Element second_stop() {
  Surface s;
  s.id = SurfaceId("STO2");
  s.aperture = CircularAperture{8.0, 0.0};
  return Element{"stop2", ElementKind::Stop, Pose::along_z(2.0), std::nullopt, {s}};
}

}  // namespace

TEST_CASE("active: an inactive element keeps its surfaces, frames and media", "[active]") {
  const System on = singlet();
  System off = on;
  element(off, 1).active = false;
  const CompiledSystem a = compiled(on);
  const CompiledSystem b = compiled(off);
  REQUIRE(b.surfaces().size() == a.surfaces().size());
  REQUIRE(b.elements().size() == a.elements().size());
  REQUIRE(b.node_frames().size() == a.node_frames().size());
  for (std::size_t i = 0; i < a.surfaces().size(); ++i) {
    INFO(a.surfaces()[i].id.str());
    CHECK(b.surfaces()[i].id == a.surfaces()[i].id);
    // Geometry stays (ADR 0049, point 2): the same frames bit for bit.
    CHECK(b.surfaces()[i].to_global.rotation() == a.surfaces()[i].to_global.rotation());
    CHECK(b.surfaces()[i].to_global.translation() == a.surfaces()[i].to_global.translation());
    CHECK(b.surfaces()[i].medium_front == a.surfaces()[i].medium_front);
    CHECK(b.surfaces()[i].medium_back == a.surfaces()[i].medium_back);
    const bool in_lens = a.surfaces()[i].element == 1;
    CHECK(a.surfaces()[i].active);
    CHECK(b.surfaces()[i].active == !in_lens);
  }
  CHECK(a.elements()[1].active);
  CHECK_FALSE(b.elements()[1].active);
  CHECK(b.elements()[0].active);
  CHECK(b.elements()[1].media == a.elements()[1].media);
  // The automatic path leaves the lens out.
  CHECK(event_ids(a) == std::vector<std::string>{"STO", "L1.S1", "L1.S2", "IMG"});
  CHECK(event_ids(b) == std::vector<std::string>{"STO", "IMG"});
}

TEST_CASE("active: a relatively placed follower keeps its place", "[active]") {
  // ADR 0049, point 2: the image placed relative to the preceding surface (L1.S2) does not move
  // when L1 is switched off.
  System on = singlet();
  element(on, 2).pose = Pose::along_z(106.363 - 9.0);
  element(on, 2).pose.reference = rtt::model::PoseReference::RelativeToPreceding;
  System off = on;
  element(off, 1).active = false;
  const CompiledSystem a = compiled(on);
  const CompiledSystem b = compiled(off);
  const std::uint32_t img = a.find_surface(SurfaceId("IMG")).value();
  CHECK(b.surfaces()[img].to_global.translation() == a.surfaces()[img].to_global.translation());
  CHECK(b.surfaces()[img].to_global.rotation() == a.surfaces()[img].to_global.rotation());
}

TEST_CASE("active: an inactive assembly takes its whole subtree out (reference case 2)",
          "[active]") {
  // #380, reference case 2: a branch behind a beam splitter. The splitter is a thin element,
  // the branch an assembly with a fold mirror and a second detector. Path 0 is the automatic
  // path of the singlet, path 1 an explicit path that transmits at the splitter and reflects at
  // the fold mirror. Switching the branch off takes its surfaces out of the automatic path and
  // out of the explicit branch path.
  System s = singlet();
  Surface bs;
  bs.id = SurfaceId("BS");
  bs.aperture = CircularAperture{20.0, 0.0};
  const Element splitter{
      "splitter", ElementKind::ThinElement, Pose::along_z(60.0), std::nullopt, {bs}};
  Surface fold;
  fold.id = SurfaceId("FOLD");
  fold.aperture = CircularAperture{20.0, 0.0};
  Surface det;
  det.id = SurfaceId("IMG2");
  Assembly branch;
  branch.name = "branch";
  branch.children.push_back(
      {Element{"fold", ElementKind::Mirror, Pose::along_z(80.0), std::nullopt, {fold}}});
  branch.children.push_back(
      {Element{"image2", ElementKind::Detector, Pose::along_z(90.0), std::nullopt, {det}}});
  s.root.children.insert(s.root.children.begin() + 2, Node{splitter});
  s.root.children.insert(s.root.children.begin() + 3, Node{branch});
  s.paths.push_back({"branch",
                     false,
                     {{SurfaceId("STO"), EventKind::Transmit, 0},
                      {SurfaceId("L1.S1"), EventKind::Refract, 0},
                      {SurfaceId("L1.S2"), EventKind::Refract, 0},
                      {SurfaceId("BS"), EventKind::Transmit, 0},
                      {SurfaceId("FOLD"), EventKind::Reflect, 0},
                      {SurfaceId("IMG2"), EventKind::Transmit, 0}}});
  const CompiledSystem on = compiled(s);
  CHECK(event_ids(on, 0) ==
        std::vector<std::string>{"STO", "L1.S1", "L1.S2", "BS", "FOLD", "IMG2", "IMG"});

  std::get<Assembly>(s.root.children[3].value).active = false;
  const std::vector<Diagnostic> v = rtt::model::validate(s);
  CHECK_FALSE(rtt::model::has_errors(v));
  CHECK(count(v, "paths.inactive_surface", "/paths/1/events/4") == 1);
  CHECK(count(v, "paths.inactive_surface", "/paths/1/events/5") == 1);
  const CompiledSystem off = compiled(s);
  CHECK(event_ids(off, 0) == std::vector<std::string>{"STO", "L1.S1", "L1.S2", "BS", "IMG"});
  CHECK(event_ids(off, 1) == std::vector<std::string>{"STO", "L1.S1", "L1.S2", "BS"});
  // The fold mirror inside the inactive assembly is inactive too.
  CHECK(count(off.diagnostics(), "paths.inactive_mirror", "/root/children/3/children/0") == 1);
  const std::uint32_t fold_index = off.find_surface(SurfaceId("FOLD")).value();
  CHECK_FALSE(off.surfaces()[fold_index].active);
  CHECK_FALSE(off.elements()[off.surfaces()[fold_index].element].active);
}

TEST_CASE("active: explicit path events at inactive surfaces are skipped", "[active]") {
  System s = singlet();
  s.root.children.insert(s.root.children.begin() + 1, Node{window()});
  s.paths = {{"main",
              false,
              {{SurfaceId("STO"), EventKind::Transmit, 0},
               {SurfaceId("W.S1"), EventKind::Refract, 0},
               {SurfaceId("W.S2"), EventKind::Refract, 0},
               {SurfaceId("L1.S1"), EventKind::Refract, 0},
               {SurfaceId("L1.S2"), EventKind::Refract, 0},
               {SurfaceId("IMG"), EventKind::Transmit, 0}}}};
  element(s, 1).active = false;
  const CompiledSystem cs = compiled(s);
  CHECK(event_ids(cs) == std::vector<std::string>{"STO", "L1.S1", "L1.S2", "IMG"});
  CHECK(count(cs.diagnostics(), "paths.inactive_surface", "/paths/0/events/1") == 1);
  CHECK(count(cs.diagnostics(), "paths.inactive_surface", "/paths/0/events/2") == 1);
  CHECK(count_code(cs.diagnostics(), "paths.inactive_surface") == 2);
  // Media as without the window: environment up to the lens, glass inside, environment after.
  const auto& events = cs.path(PathId{0}).events;
  const std::uint32_t env = cs.environment_medium();
  CHECK(events[0].medium_after == env);
  CHECK(events[1].medium_before == env);
  CHECK(cs.media()[events[1].medium_after].reference == "CONST:1.5168");
  CHECK(events[2].medium_after == env);
  CHECK(events[3].medium_before == env);
}

TEST_CASE("active: inactive mirror and inactive last element give warnings", "[active]") {
  System s = singlet();
  element(s, 2).active = false;
  const CompiledSystem cs = compiled(s);
  CHECK(event_ids(cs) == std::vector<std::string>{"STO", "L1.S1", "L1.S2"});
  CHECK(count(cs.diagnostics(), "paths.image_changed", "/paths/0") == 1);
  // An inactive element in the middle changes nothing at the image.
  System mid = singlet();
  element(mid, 1).active = false;
  CHECK(count_code(compiled(mid).diagnostics(), "paths.image_changed") == 0);
  // An inactive mirror: one warning per mirror, also with two automatic paths.
  System m = singlet();
  element(m, 2).kind = ElementKind::Mirror;
  element(m, 2).active = false;
  m.paths.push_back({"second", true, {}});
  const CompiledSystem cm = compiled(m);
  CHECK(count(cm.diagnostics(), "paths.inactive_mirror", "/root/children/2") == 1);
  CHECK(count(cm.diagnostics(), "paths.image_changed", "/paths/0") == 1);
  CHECK(count(cm.diagnostics(), "paths.image_changed", "/paths/1") == 1);
}

TEST_CASE("active: automatic path rules hold only for active elements", "[active]") {
  // A Mangin mirror on an automatic path is an error (paths.mangin_mirror_automatic); switched
  // off it is skipped and the error is gone.
  System s = singlet();
  element(s, 1).kind = ElementKind::Mirror;
  CHECK_THROWS_AS(compiled(s), CompileError);
  element(s, 1).active = false;
  const CompiledSystem cs = compiled(s);
  CHECK(count_code(cs.diagnostics(), "paths.mangin_mirror_automatic") == 0);
  CHECK(event_ids(cs) == std::vector<std::string>{"STO", "IMG"});
}

TEST_CASE("active: materials of inactive elements are still resolved", "[active]") {
  // ADR 0049, point 2: an unknown material stays an error, so that switching on never brings a
  // new error in the middle of an edit.
  System s = singlet();
  element(s, 1).material = "NOPE:GLASS";
  element(s, 1).active = false;
  try {
    static_cast<void>(compiled(s));
    FAIL("no CompileError");
  } catch (const CompileError& e) {
    CHECK(count(e.diagnostics(), "material.unknown", "/root/children/1/material") == 1);
  }
}

TEST_CASE("active: path-dependent wavelength checks see only active elements", "[active]") {
  // A material valid only on [0.6, 0.7] um (system wavelength 0.5876 um): an error on a path,
  // none if its element is inactive.
  class Narrow final : public rtt::material::Material {
   public:
    [[nodiscard]] rtt::math::Complex index(double /*wavelength_um*/,
                                           double /*temperature_c*/,
                                           double /*pressure_atm*/) const override {
      return {1.6, 0.0};
    }
    [[nodiscard]] std::optional<rtt::material::WavelengthRange> wavelength_range_um()
        const override {
      return rtt::material::WavelengthRange{0.6, 0.7};
    }
  };
  MaterialLibrary lib;
  lib.add("TEST:NARROW", std::make_shared<const Narrow>());
  System s = singlet();
  element(s, 1).material = "TEST:NARROW";
  CHECK_THROWS_AS(rtt::compile::compile(s, lib), CompileError);
  element(s, 1).active = false;
  CHECK_NOTHROW(rtt::compile::compile(s, lib));
}

TEST_CASE("active: a path without events gives paths.all_inactive", "[active]") {
  System s = singlet();
  for (Node& child : s.root.children) std::get<Element>(child.value).active = false;
  const CompiledSystem cs = compiled(s);
  CHECK(cs.path(PathId{0}).events.empty());
  CHECK(count(cs.diagnostics(), "paths.all_inactive", "/paths/0") == 1);
  CHECK(count_code(cs.diagnostics(), "paths.image_changed") == 0);
  // The same for an explicit path whose surfaces are all inactive.
  System e = singlet();
  e.paths.push_back(
      {"lens",
       false,
       {{SurfaceId("L1.S1"), EventKind::Refract, 0}, {SurfaceId("L1.S2"), EventKind::Refract, 0}}});
  element(e, 1).active = false;
  const CompiledSystem ce = compiled(e);
  CHECK(ce.path(PathId{1}).events.empty());
  CHECK(count(ce.diagnostics(), "paths.all_inactive", "/paths/1") == 1);
  CHECK(count_code(ce.diagnostics(), "paths.all_inactive") == 1);
}

TEST_CASE("active: an inactive stop counts as no stop", "[active][stop]") {
  System s = singlet();
  element(s, 0).active = false;
  const CompiledSystem cs = compiled(s);
  CHECK_FALSE(cs.stop_surface().has_value());
  CHECK(count(cs.diagnostics(), "stop.inactive", "/root/children/0") == 1);
  // stop.not_on_path is not reported for an inactive stop.
  CHECK(count_code(cs.diagnostics(), "stop.not_on_path") == 0);
  CHECK_THROWS_AS(rtt::compile::require_stop(cs, PathId{0}), rtt::compile::NoStopError);
  // The structure of an inactive stop is still checked.
  element(s, 0).surfaces[0].aperture.reset();
  CHECK(rtt::model::test::has_error_at(rtt::model::validate(s),
                                       "/root/children/0/surfaces/0/aperture"));
}

TEST_CASE("active: stop.multiple counts only active stops; switching between stops", "[active]") {
  System s = singlet();
  s.root.children.insert(s.root.children.begin() + 1, Node{second_stop()});
  // Two active stops: an error.
  CHECK(rtt::model::test::has_error_at(rtt::model::validate(s), "/root/children/1"));
  // One of them inactive: the other is the stop, no warning.
  element(s, 1).active = false;
  const CompiledSystem a = compiled(s);
  CHECK(a.surfaces()[a.stop_surface().value()].id == SurfaceId("STO"));
  CHECK(count_code(a.diagnostics(), "stop.inactive") == 0);
  // Switched over: the second stop is the stop.
  element(s, 1).active = true;
  element(s, 0).active = false;
  const CompiledSystem b = compiled(s);
  CHECK(b.surfaces()[b.stop_surface().value()].id == SurfaceId("STO2"));
  CHECK(count_code(b.diagnostics(), "stop.inactive") == 0);
  // Both inactive: one warning at the first stop.
  element(s, 1).active = false;
  const CompiledSystem c = compiled(s);
  CHECK_FALSE(c.stop_surface().has_value());
  CHECK(count(c.diagnostics(), "stop.inactive", "/root/children/0") == 1);
  CHECK(count_code(c.diagnostics(), "stop.inactive") == 1);
}

TEST_CASE("active: aperture.stop on an inactive surface falls back to an active stop element",
          "[active]") {
  // ADR 0049, point 3: the surface of aperture.stop if its element is active, else the first
  // active Stop element. With both active they are stop.multiple (ADR 0039).
  System s = singlet();
  element(s, 1).surfaces[0].aperture = CircularAperture{10.0, 0.0};
  s.aperture.stop = SurfaceId("L1.S1");
  CHECK(rtt::model::test::has_error_at(rtt::model::validate(s), "/aperture/stop"));
  element(s, 0).active = false;  // only aperture.stop is active
  const CompiledSystem a = compiled(s);
  CHECK(a.surfaces()[a.stop_surface().value()].id == SurfaceId("L1.S1"));
  element(s, 0).active = true;
  element(s, 1).active = false;  // the stop surface is inactive, the Stop element active
  const CompiledSystem b = compiled(s);
  CHECK(b.surfaces()[b.stop_surface().value()].id == SurfaceId("STO"));
  CHECK(count_code(b.diagnostics(), "stop.inactive") == 0);
  element(s, 0).active = false;  // neither
  const CompiledSystem c = compiled(s);
  CHECK_FALSE(c.stop_surface().has_value());
  CHECK(count(c.diagnostics(), "stop.inactive", "/aperture/stop") == 1);
}

TEST_CASE("active: switching off never adds an error", "[active]") {
  // ADR 0049, point 8: every element of the singlet and of the window system can be switched
  // off alone, and all together, without an error.
  System base = singlet();
  base.root.children.insert(base.root.children.begin() + 1, Node{window()});
  base.paths.push_back({"explicit",
                        false,
                        {{SurfaceId("STO"), EventKind::Transmit, 0},
                         {SurfaceId("W.S1"), EventKind::Refract, 0},
                         {SurfaceId("W.S2"), EventKind::Refract, 0},
                         {SurfaceId("L1.S1"), EventKind::Refract, 0},
                         {SurfaceId("L1.S2"), EventKind::Refract, 0},
                         {SurfaceId("IMG"), EventKind::Transmit, 0}}});
  // Ray operands at every surface of the explicit path, and the system aperture from the stop.
  for (const char* id : {"STO", "W.S1", "L1.S2", "IMG"}) {
    rtt::model::RayOperand op;
    op.common = {0.0, 1.0, std::nullopt};
    op.path = "explicit";
    op.surface = SurfaceId(id);
    op.py = 1.0;
    base.optimization.operands.emplace_back(op);
  }
  base.aperture = {rtt::model::SystemApertureType::StopSize, rtt::model::Param(10.0)};
  REQUIRE_FALSE(rtt::model::has_errors(rtt::model::validate(base)));
  for (std::size_t i = 0; i < base.root.children.size(); ++i) {
    INFO("element " << i);
    System s = base;
    element(s, i).active = false;
    CHECK_FALSE(rtt::model::has_errors(rtt::model::validate(s)));
    CHECK_NOTHROW(compiled(s));
  }
  System all = base;
  for (Node& child : all.root.children) std::get<Element>(child.value).active = false;
  CHECK_FALSE(rtt::model::has_errors(rtt::model::validate(all)));
  CHECK_NOTHROW(compiled(all));
}

TEST_CASE("active: object.reference may name an inactive surface", "[active]") {
  System s = singlet();
  s.object.at_infinity = false;
  s.object.distance = rtt::model::Param(100.0);
  s.object.reference = SurfaceId("L1.S1");
  s.aperture = {rtt::model::SystemApertureType::ObjectSpaceNA, rtt::model::Param(0.05)};
  s.fields = {rtt::model::FieldType::ObjectHeight, {{0.0, 0.0, 1.0}}};
  const CompiledSystem a = compiled(s);
  element(s, 1).active = false;
  const CompiledSystem b = compiled(s);
  CHECK(b.object_z() == a.object_z());
  CHECK(b.object_z() == 5.0 - 100.0);
}

TEST_CASE("active: errors at events behind skipped events name the event of the model",
          "[active]") {
  // Review of PR 3: after skipped events the compiled index is not the index in the model.
  // Ordinary at L1.S1 without a crystal is paths.mode_without_crystal at the event of L1.S1,
  // index 3 in the model, index 1 in the compiled path.
  System s = singlet();
  s.root.children.insert(s.root.children.begin() + 1, Node{window()});
  element(s, 1).active = false;
  s.paths = {{"main",
              false,
              {{SurfaceId("STO"), EventKind::Transmit, 0},
               {SurfaceId("W.S1"), EventKind::Refract, 0},
               {SurfaceId("W.S2"), EventKind::Refract, 0},
               {SurfaceId("L1.S1"), EventKind::Ordinary, 0},
               {SurfaceId("L1.S2"), EventKind::Refract, 0},
               {SurfaceId("IMG"), EventKind::Transmit, 0}}}};
  try {
    static_cast<void>(compiled(s));
    FAIL("no CompileError");
  } catch (const CompileError& e) {
    CHECK(count(e.diagnostics(), "paths.mode_without_crystal", "/paths/0/events/3") == 1);
    CHECK(count_code(e.diagnostics(), "paths.mode_without_crystal") == 1);
  }
}

TEST_CASE("active: refract at an inactive ideal lens is still an error", "[active]") {
  // ADR 0049, Nachtrag PR 3: paths.ideal_lens_event is a static path rule, checked also at
  // inactive surfaces, so that switching on brings no new error.
  System s = singlet();
  element(s, 0).kind = ElementKind::ThinElement;
  element(s, 0).surfaces[0].interaction =
      rtt::model::IdealLens{rtt::model::Param(50.0), std::nullopt};
  element(s, 0).active = false;
  s.paths = {{"x",
              false,
              {{SurfaceId("STO"), EventKind::Refract, 0},
               {SurfaceId("L1.S1"), EventKind::Refract, 0},
               {SurfaceId("L1.S2"), EventKind::Refract, 0}}}};
  try {
    static_cast<void>(compiled(s));
    FAIL("no CompileError");
  } catch (const CompileError& e) {
    CHECK(count(e.diagnostics(), "paths.ideal_lens_event", "/paths/0/events/0") == 1);
  }
}

TEST_CASE("active: an active assembly inside an inactive one stays inactive", "[active]") {
  // ADR 0049, point 1: inactive if the node or an ancestor has active false. The lens sits in an
  // active inner assembly inside an inactive outer one; an explicit event behind the outer
  // assembly is not reported (the state is restored after the subtree).
  System s = singlet();
  Assembly inner;
  inner.name = "inner";
  inner.children.push_back(s.root.children[1]);
  Assembly outer;
  outer.name = "outer";
  outer.active = false;
  outer.children.push_back(Node{inner});
  s.root.children[1] = Node{outer};
  s.paths.push_back({"explicit",
                     false,
                     {{SurfaceId("STO"), EventKind::Transmit, 0},
                      {SurfaceId("L1.S1"), EventKind::Refract, 0},
                      {SurfaceId("L1.S2"), EventKind::Refract, 0},
                      {SurfaceId("IMG"), EventKind::Transmit, 0}}});
  const std::vector<Diagnostic> v = rtt::model::validate(s);
  CHECK(count(v, "paths.inactive_surface", "/paths/1/events/1") == 1);
  CHECK(count(v, "paths.inactive_surface", "/paths/1/events/2") == 1);
  CHECK(count_code(v, "paths.inactive_surface") == 2);  // not at STO or IMG
  const CompiledSystem cs = compiled(s);
  CHECK(event_ids(cs, 0) == std::vector<std::string>{"STO", "IMG"});
  CHECK(event_ids(cs, 1) == std::vector<std::string>{"STO", "IMG"});
  const std::uint32_t l1 = cs.find_surface(SurfaceId("L1.S1")).value();
  CHECK_FALSE(cs.surfaces()[l1].active);
  CHECK(cs.surfaces()[cs.find_surface(SurfaceId("IMG")).value()].active);
}

TEST_CASE("active: coating wavelength checks see only active surfaces", "[active]") {
  // A coating whose layer material is valid only on [0.6, 0.7] um (system wavelength
  // 0.5876 um): coating.wavelength_out_of_range on a path, none at an inactive surface.
  class Narrow final : public rtt::material::Material {
   public:
    [[nodiscard]] rtt::math::Complex index(double /*wavelength_um*/,
                                           double /*temperature_c*/,
                                           double /*pressure_atm*/) const override {
      return {1.38, 0.0};
    }
    [[nodiscard]] std::optional<rtt::material::WavelengthRange> wavelength_range_um()
        const override {
      return rtt::material::WavelengthRange{0.6, 0.7};
    }
  };
  MaterialLibrary lib;
  lib.add("TEST:NARROW", std::make_shared<const Narrow>());
  rtt::coating::CoatingLibrary coatings;
  coatings.add(rtt::coating::CoatingCatalog{
      "T",
      {rtt::coating::CoatingDesign{
          "NARROW",
          "",
          std::nullopt,
          {rtt::coating::LayerSpec{"TEST:NARROW", rtt::coating::PhysicalThickness{0.1}}}}}});
  System s = singlet();
  element(s, 1).surfaces[0].interaction = CoatingRef{"T:NARROW"};
  CHECK_THROWS_AS(rtt::compile::compile(s, lib, coatings), CompileError);
  element(s, 1).active = false;
  CHECK_NOTHROW(rtt::compile::compile(s, lib, coatings));
}
