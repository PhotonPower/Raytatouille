// Schema 0.3 (ADR 0025, ADR 0026): orders at every event without "diffract", diffraction
// efficiency, crystal material and optic axis; migration of 0.1 and 0.2 files.

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#include "rtt/io/edit.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/model/validate.hpp"
#include "test_support.hpp"

using Catch::Matchers::ContainsSubstring;
using rtt::model::CrystalMaterial;
using rtt::model::DiffractionEfficiency;
using rtt::model::Element;
using rtt::model::EventKind;
using rtt::model::Param;
using rtt::model::SurfaceId;
using rtt::model::System;
using rtt::model::test::element;

namespace {

/// `text` with the only occurrence of `from` replaced by `to`; fails if there is none or more.
std::string replace_once(const std::string& text, const std::string& from, const std::string& to) {
  const std::size_t at = text.find(from);
  REQUIRE(at != std::string::npos);
  REQUIRE(text.find(from, at + 1) == std::string::npos);
  std::string out = text;
  out.replace(at, from.size(), to);
  return out;
}

std::string error_pointer(const std::string& text) {
  try {
    (void)rtt::io::parse_system(text);
  } catch (const rtt::io::ParseError& e) {
    return e.pointer();
  }
  return "<no error>";
}

/// Singlet with a grating on L1.S1 and an explicit path with a transmitted first order there.
System grating_singlet() {
  System s = rtt::model::test::make_singlet();
  element(s, 1).surfaces[0].phases.emplace_back(rtt::model::LinearGrating{Param(300.0), 0.0});
  s.paths = {{"first order",
              false,
              {{SurfaceId("STO"), EventKind::Transmit, 0},
               {SurfaceId("L1.S1"), EventKind::Transmit, 1},
               {SurfaceId("L1.S2"), EventKind::Refract, 0},
               {SurfaceId("IMG"), EventKind::Transmit, 0}}}};
  return s;
}

/// Singlet whose lens L1 is a calcite crystal with its optic axis at 45 degree in the y-z plane.
System crystal_singlet() {
  System s = rtt::model::test::make_singlet();
  Element& l = element(s, 1);
  l.material.reset();
  l.crystal = CrystalMaterial{"BIREFRINGENT:CALCITE", "BIREFRINGENT:CALCITE-E"};
  l.optic_axis = std::array<double, 3>{0.0, 1.0, 1.0};
  return s;
}

constexpr const char* kTransmitFirst = R"("surface": "L1.S1", "kind": "transmit", "order": 1)";
constexpr const char* kDiffractFirst = R"("surface": "L1.S1", "kind": "diffract", "order": 1)";

}  // namespace

TEST_CASE("schema 0.3: diffract of a 0.1 or 0.2 file is read as transmit with its order",
          "[io][format03]") {
  // ADR 0025, point 4: both kept the medium, so the migration does not change the path.
  const System s = grating_singlet();
  const std::string current = rtt::io::to_json(s);
  REQUIRE_THAT(current, ContainsSubstring(R"("schema_version": "0.5.0")"));
  for (const std::string version : {"0.2.0", "0.1.0"}) {
    INFO(version);
    const std::string old = replace_once(replace_once(current, R"("schema_version": "0.5.0")",
                                                      R"("schema_version": ")" + version + R"(")"),
                                         kTransmitFirst, kDiffractFirst);
    const System back = rtt::io::parse_system(old);
    REQUIRE(back == s);
    REQUIRE(back.schema_version == "0.5.0");
    REQUIRE(rtt::io::to_json(back) == current);
  }
}

TEST_CASE("schema 0.3: diffract is no event kind of a 0.3 file", "[io][format03]") {
  const std::string text =
      replace_once(rtt::io::to_json(grating_singlet()), kTransmitFirst, kDiffractFirst);
  REQUIRE(error_pointer(text) == "/paths/0/events/1/kind");
}

TEST_CASE("schema 0.3: a migrated order at a surface without phase layer is a validate error",
          "[io][format03]") {
  // ADR 0025, point 4: the migration does not set the order to 0; such paths were never
  // traceable (EventImpossible), and validate now reports them.
  System s = grating_singlet();
  element(s, 1).surfaces[0].phases.clear();
  const std::string old =
      replace_once(replace_once(rtt::io::to_json(s), R"("schema_version": "0.5.0")",
                                R"("schema_version": "0.2.0")"),
                   kTransmitFirst, kDiffractFirst);
  const System back = rtt::io::parse_system(old);
  REQUIRE(back.paths[0].events[1].kind == EventKind::Transmit);
  REQUIRE(back.paths[0].events[1].order == 1);
  const auto diagnostics = rtt::model::validate(back);
  REQUIRE(diagnostics.size() == 1);
  REQUIRE(diagnostics[0].code == "paths.order_not_allowed");
  REQUIRE(diagnostics[0].location == "/paths/0/events/1/order");
}

TEST_CASE("schema 0.3: crystal material and optic axis survive file -> model -> file",
          "[io][format03][roundtrip]") {
  const System s = crystal_singlet();
  const std::string text = rtt::io::to_json(s);
  REQUIRE_THAT(text, ContainsSubstring(R"("ordinary": "BIREFRINGENT:CALCITE")"));
  REQUIRE_THAT(text, ContainsSubstring(R"("extraordinary": "BIREFRINGENT:CALCITE-E")"));
  REQUIRE_THAT(text, ContainsSubstring(R"("optic_axis": [0.0, 1.0, 1.0])"));
  const System back = rtt::io::parse_system(text);
  REQUIRE(back == s);
  REQUIRE_FALSE(std::get<Element>(back.root.children[1].value).material.has_value());
  REQUIRE(rtt::io::to_json(back) == text);
  REQUIRE(rtt::model::validate(back).empty());
}

TEST_CASE("schema 0.3: diffraction efficiencies survive file -> model -> file",
          "[io][format03][roundtrip]") {
  System s = grating_singlet();
  element(s, 1).surfaces[0].diffraction_efficiency =
      std::vector<DiffractionEfficiency>{{1, 0.8}, {0, 0.05}, {-1, 0.0}};
  const std::string text = rtt::io::to_json(s);
  REQUIRE_THAT(text, ContainsSubstring(R"("diffraction_efficiency")"));
  const System back = rtt::io::parse_system(text);
  REQUIRE(back == s);
  REQUIRE(rtt::io::to_json(back) == text);
  REQUIRE(rtt::model::validate(back).empty());
}

TEST_CASE("schema 0.3: a crystal object in a 0.1 or 0.2 file is an error", "[io][format03]") {
  // ADR 0026, point 6: like a material list in a 0.1 file (ADR 0017), a newer form never
  // silently upgrades the version.
  const std::string current = rtt::io::to_json(crystal_singlet());
  for (const std::string version : {"0.2.0", "0.1.0"}) {
    INFO(version);
    const std::string old = replace_once(current, R"("schema_version": "0.5.0")",
                                         R"("schema_version": ")" + version + R"(")");
    REQUIRE(error_pointer(old) == "/root/children/1/material");
    REQUIRE_THROWS_WITH(rtt::io::parse_system(old),
                        ContainsSubstring("crystal materials need schema_version 0.3"));
  }
}

TEST_CASE("schema 0.3: structural errors of the new keys carry a JSON pointer", "[io][format03]") {
  const std::string crystal = rtt::io::to_json(crystal_singlet());
  REQUIRE(error_pointer(replace_once(crystal, R"("extraordinary": "BIREFRINGENT:CALCITE-E")",
                                     R"("e": "BIREFRINGENT:CALCITE-E")")) ==
          "/root/children/1/material/e");
  REQUIRE(error_pointer(replace_once(crystal, R"("extraordinary": "BIREFRINGENT:CALCITE-E")",
                                     R"("extraordinary": 1)")) ==
          "/root/children/1/material/extraordinary");
  REQUIRE(error_pointer(replace_once(crystal, R"("optic_axis": [0.0, 1.0, 1.0])",
                                     R"("optic_axis": [0.0, 1.0])")) ==
          "/root/children/1/optic_axis");

  System g = grating_singlet();
  element(g, 1).surfaces[0].diffraction_efficiency = std::vector<DiffractionEfficiency>{{1, 0.8}};
  const std::string grating = rtt::io::to_json(g);
  const std::string entry = R"("order": 1, "efficiency": 0.8)";
  REQUIRE(error_pointer(replace_once(grating, entry, R"("efficiency": 0.8)")) ==
          "/root/children/1/surfaces/0/diffraction_efficiency/0");
  REQUIRE(error_pointer(replace_once(grating, entry, R"("order": 1.5, "efficiency": 0.8)")) ==
          "/root/children/1/surfaces/0/diffraction_efficiency/0/order");
  REQUIRE(error_pointer(replace_once(grating, entry, R"("order": 1, "efficiency": "high")")) ==
          "/root/children/1/surfaces/0/diffraction_efficiency/0/efficiency");
}

TEST_CASE("schema 0.3: the edit form carries crystal, optic axis and efficiencies",
          "[io][format03][edit]") {
  // ADR 0026, point 6: the edit form (ADR 0024) writes the material object as read and the
  // optic axis if set; patches reach the new keys like any other.
  System s = crystal_singlet();
  element(s, 1).surfaces[0].phases.emplace_back(rtt::model::LinearGrating{Param(300.0), 0.0});
  const std::string edit = rtt::io::to_edit_json(s);
  REQUIRE_THAT(edit, ContainsSubstring(R"("optic_axis": [0.0, 1.0, 1.0])"));
  REQUIRE_THAT(edit, ContainsSubstring(R"("extraordinary": "BIREFRINGENT:CALCITE-E")"));
  REQUIRE_THAT(edit, !ContainsSubstring("diffraction_efficiency"));  // absent stays absent

  const rtt::io::PatchResult r = rtt::io::apply_patch(s, R"([
    {"op": "replace", "path": "/root/children/1/optic_axis", "value": [1.0, 0.0, 0.0]},
    {"op": "replace", "path": "/root/children/1/material/ordinary", "value": "CONST:1.66"},
    {"op": "add", "path": "/root/children/1/surfaces/0/diffraction_efficiency",
     "value": [{"order": 1, "efficiency": 0.5}]}
  ])");
  const Element& patched = std::get<Element>(r.system.root.children[1].value);
  REQUIRE(patched.optic_axis == std::array<double, 3>{1.0, 0.0, 0.0});
  REQUIRE(patched.crystal->ordinary == "CONST:1.66");
  REQUIRE(patched.surfaces[0].diffraction_efficiency ==
          std::vector<DiffractionEfficiency>{{1, 0.5}});

  const rtt::io::PatchResult back =
      rtt::io::apply_patch(r.system, r.inverse, rtt::io::EditCheck::StructureOnly);
  REQUIRE(back.system == s);
}

TEST_CASE("schema 0.3: crystal and isotropic material together cannot be written",
          "[io][format03]") {
  System s = crystal_singlet();
  element(s, 1).material = "SCHOTT:N-BK7";
  REQUIRE_THROWS_AS(rtt::io::to_json(s), std::invalid_argument);
}
