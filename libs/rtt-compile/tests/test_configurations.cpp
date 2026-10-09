// Configurations in compile (ADR 0029, point 5; #165): compile(system, materials[, coatings],
// configuration) evaluates the parameter table for one column and sets the bound Params before
// the poses are composed. Without bound Params nothing changes (bitwise); column k is bitwise
// the system with the values of column k entered by hand.

#include <algorithm>
#include <bit>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <tuple>
#include <variant>
#include <vector>

#include "rtt/coating/catalog.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/compile/ghosts.hpp"
#include "rtt/compile/layout.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"

namespace fs = std::filesystem;
using Catch::Matchers::ContainsSubstring;
using rtt::compile::compile;
using rtt::compile::CompiledSystem;
using rtt::compile::CompileError;
using rtt::material::MaterialLibrary;
using rtt::math::Isometry3;
using rtt::math::Vec3;
using rtt::model::CircularAperture;
using rtt::model::Conic;
using rtt::model::Element;
using rtt::model::ElementKind;
using rtt::model::Param;
using rtt::model::ParameterExpression;
using rtt::model::ParameterRow;
using rtt::model::Pose;
using rtt::model::PoseReference;
using rtt::model::Surface;
using rtt::model::SurfaceId;
using rtt::model::System;

namespace {

bool same(double a, double b) {
  return std::bit_cast<std::uint64_t>(a) == std::bit_cast<std::uint64_t>(b);
}

bool same(const Isometry3& a, const Isometry3& b) {
  for (int r = 0; r < 3; ++r) {
    if (!same(a.translation()(r), b.translation()(r))) return false;
    for (int c = 0; c < 3; ++c) {
      if (!same(a.rotation()(r, c), b.rotation()(r, c))) return false;
    }
  }
  return true;
}

/// Asserts that two compiled systems have bitwise the same geometry, media, paths and system
/// values: surface transforms, apertures, sags at a few points, media indices and path events.
void require_same(const CompiledSystem& a, const CompiledSystem& b) {
  REQUIRE(a.surfaces().size() == b.surfaces().size());
  for (std::uint32_t i = 0; i < a.surfaces().size(); ++i) {
    const auto& sa = a.surfaces()[i];
    const auto& sb = b.surfaces()[i];
    INFO("surface " << sa.id.str());
    REQUIRE(sa.id == sb.id);
    CHECK(same(sa.to_global, sb.to_global));
    CHECK(same(sa.to_local, sb.to_local));
    CHECK(sa.aperture == sb.aperture);
    const double xs[] = {0.0, 0.5, 1.0};
    const double ys[] = {1.0, 0.3, -0.7};
    double za[3] = {};
    double zb[3] = {};
    rtt::compile::surface_sag(a, i, xs, ys, za);
    rtt::compile::surface_sag(b, i, xs, ys, zb);
    for (int k = 0; k < 3; ++k) CHECK(same(za[k], zb[k]));
  }
  REQUIRE(a.media().size() == b.media().size());
  for (std::size_t m = 0; m < a.media().size(); ++m) {
    CHECK(a.media()[m].reference == b.media()[m].reference);
    for (std::size_t w = 0; w < a.media()[m].index.size(); ++w) {
      CHECK(same(a.media()[m].index[w].real(), b.media()[m].index[w].real()));
      CHECK(same(a.media()[m].index[w].imag(), b.media()[m].index[w].imag()));
    }
  }
  REQUIRE(a.paths().size() == b.paths().size());
  for (std::size_t p = 0; p < a.paths().size(); ++p) {
    REQUIRE(a.paths()[p].events.size() == b.paths()[p].events.size());
    for (std::size_t e = 0; e < a.paths()[p].events.size(); ++e) {
      const auto& ea = a.paths()[p].events[e];
      const auto& eb = b.paths()[p].events[e];
      CHECK(ea.surface == eb.surface);
      CHECK(ea.kind == eb.kind);
      CHECK(ea.medium_before == eb.medium_before);
      CHECK(ea.medium_after == eb.medium_after);
    }
  }
  CHECK(a.aperture() == b.aperture());
  CHECK(a.object() == b.object());
  CHECK(a.wavelengths_um() == b.wavelengths_um());
}

Surface plane(const std::string& id, double z, double radius = 5.0) {
  Surface s;
  s.id = SurfaceId(id);
  s.pose = Pose::along_z(z);
  s.aperture = CircularAperture{radius, 0.0};
  return s;
}

Surface sphere(const std::string& id, double z, double r) {
  Surface s = plane(id, z);
  s.shape.base = Conic{Param(r), Param(0.0)};
  return s;
}

ParameterRow row(std::string name, rtt::model::ParameterForm form) {
  ParameterRow r;
  r.name = std::move(name);
  r.form = std::move(form);
  return r;
}

/// A two-lens zoom (ADR 0029, point 8) in vacuum: stop at z = 0, L1 (CONST:1.5, R 40 / plane,
/// 3 mm) at z = 10; L2 (R 30 / plane, 2 mm) placed relative to L1.S2 at the air gap G, the image
/// relative to L2.S2 at B = TOTAL - G. Configurations "wide" (G = 20) and "tele" (G = 5).
/// With `gap` and `back` set, the same system with these numbers entered by hand and no table.
System zoom(std::optional<double> gap = std::nullopt, std::optional<double> back = std::nullopt) {
  System s;
  s.name = "zoom";
  s.environment.medium = "VACUUM";
  s.wavelengths = {{0.5876, 1.0, true}};
  s.aperture = {rtt::model::SystemApertureType::EntrancePupilDiameter, Param(2.0)};
  s.fields = {rtt::model::FieldType::AngleDeg, {{0.0, 0.0, 1.0}}};
  s.root.name = "root";
  s.root.children.push_back(
      {Element{"stop", ElementKind::Stop, {}, std::nullopt, {plane("STO", 0.0)}}});
  s.root.children.push_back({Element{"L1",
                                     ElementKind::Lens,
                                     Pose::along_z(10.0),
                                     "CONST:1.5",
                                     {sphere("L1.S1", 0.0, 40.0), plane("L1.S2", 3.0)}}});
  Element l2{
      "L2", ElementKind::Lens, {}, "CONST:1.5", {sphere("L2.S1", 0.0, 30.0), plane("L2.S2", 2.0)}};
  l2.pose.reference = PoseReference::RelativeToPreceding;
  Element image{"image", ElementKind::Detector, {}, std::nullopt, {plane("IMG", 0.0, 20.0)}};
  image.pose.reference = PoseReference::RelativeToPreceding;
  if (gap) {
    l2.pose.position[2] = Param(*gap);
    image.pose.position[2] = Param(*back);
  } else {
    s.configurations = {{"wide"}, {"tele"}};
    ParameterRow g = row("G", std::vector<double>{20.0, 5.0});
    g.variable = true;
    s.parameters = {row("TOTAL", 40.0), g, row("B", ParameterExpression{"TOTAL - G"})};
    l2.pose.position[2] = Param::bound("G");
    image.pose.position[2] = Param::bound("B");
  }
  s.root.children.push_back({l2});
  s.root.children.push_back({image});
  s.paths = {{"main", true, {}}};
  return s;
}

std::vector<fs::path> reference_files() {
  std::vector<fs::path> files;
  for (const auto& entry : fs::recursive_directory_iterator(RTT_REFERENCE_DIR)) {
    if (entry.is_regular_file() && entry.path().filename().string().ends_with(".rtt.json")) {
      files.push_back(entry.path());
    }
  }
  std::sort(files.begin(), files.end());
  return files;
}

}  // namespace

TEST_CASE("configurations: the overloads are unambiguous", "[compile][configurations]") {
  const MaterialLibrary lib;
  const rtt::coating::CoatingLibrary none;
  const System s = zoom();
  const std::size_t k = 1;
  CHECK(compile(s, lib).configuration() == 0);
  CHECK(compile(s, lib, 0).configuration() == 0);
  CHECK(compile(s, lib, 1).configuration() == 1);
  CHECK(compile(s, lib, k).configuration() == 1);
  CHECK(compile(s, lib, none).configuration() == 0);
  CHECK(compile(s, lib, none, 0).configuration() == 0);
  CHECK(compile(s, lib, none, k).configuration() == 1);
}

TEST_CASE("configurations: index and name of the compiled configuration",
          "[compile][configurations]") {
  const MaterialLibrary lib;
  const CompiledSystem tele = compile(zoom(), lib, 1);
  CHECK(tele.configuration() == 1);
  CHECK(tele.configuration_name() == "tele");
  CHECK(compile(zoom(), lib, 0).configuration_name() == "wide");
  const CompiledSystem nominal = compile(zoom(20.0, 20.0), lib);
  CHECK(nominal.configuration() == 0);
  CHECK(nominal.configuration_name().empty());
}

TEST_CASE("configurations: column k is bitwise the system with hand values",
          "[compile][configurations]") {
  // Also column 0 with a table and bindings (review of the #165 plan).
  const MaterialLibrary lib;
  const System s = zoom();
  for (const auto& [k, g, b] : {std::tuple{0U, 20.0, 20.0}, std::tuple{1U, 5.0, 35.0}}) {
    INFO("configuration " << k);
    require_same(compile(s, lib, k), compile(zoom(g, b), lib));
  }
  // The image stays at 10 + 3 + G + 2 + B = 55 in both configurations; L2 moves with G.
  for (const std::size_t k : {0U, 1U}) {
    const CompiledSystem cs = compile(s, lib, k);
    const auto z = [&](const char* id) {
      return cs.surfaces()[*cs.find_surface(SurfaceId(id))].to_global.translation().z();
    };
    CHECK(z("IMG") == 55.0);
    CHECK(z("L2.S1") == (k == 0 ? 33.0 : 18.0));
  }
}

TEST_CASE("configurations: a larger air gap moves the following element by exactly Delta D",
          "[compile][configurations]") {
  // #165 acceptance: a variable air gap D through {"param": "D"} in front of a relatively
  // placed element. Column 1 has D + 0.25; L2 and the image (relative to L2, fixed distance 50)
  // move by exactly 0.25 along z: the sums 13 + 20 = 33 and 13 + 20.25 = 33.25 are exact in
  // binary, so are the later ones.
  System s = zoom(20.0, 50.0);
  s.configurations = {{"D"}, {"D + 1/4"}};
  s.parameters = {row("D", std::vector<double>{20.0, 20.25})};
  std::get<Element>(s.root.children[2].value).pose.position[2] = Param::bound("D");
  const MaterialLibrary lib;
  const CompiledSystem a = compile(s, lib, 0);
  const CompiledSystem b = compile(s, lib, 1);
  for (const char* id : {"L2.S1", "L2.S2", "IMG"}) {
    INFO(id);
    const auto i = *a.find_surface(SurfaceId(id));
    const Vec3 pa = a.surfaces()[i].to_global.translation();
    const Vec3 pb = b.surfaces()[i].to_global.translation();
    CHECK(pb.z() - pa.z() == 0.25);
    CHECK(same(pb.x(), pa.x()));
    CHECK(same(pb.y(), pa.y()));
  }
  // Surfaces before the gap do not move.
  for (const char* id : {"STO", "L1.S1", "L1.S2"}) {
    const auto i = *a.find_surface(SurfaceId(id));
    CHECK(same(a.surfaces()[i].to_global, b.surfaces()[i].to_global));
  }
}

TEST_CASE("configurations: without bound Params compile is unchanged",
          "[compile][configurations]") {
  // Every reference system that compiles: explicit configuration 0 and an added table without
  // bindings give bitwise the same compiled system as the plain call.
  MaterialLibrary m1;
  m1.add_catalog(std::string(RTT_CATALOG_DIR) + "/schott.agf");
  MaterialLibrary m2;
  m2.add_catalog(std::string(RTT_CATALOG_DIR) + "/m2/schott.agf");
  rtt::coating::CoatingLibrary coatings;
  coatings.add_catalog(std::string(RTT_CATALOG_DIR) + "/coatings/demo.json");
  coatings.add_catalog(std::string(RTT_CATALOG_DIR) + "/coatings/m3.json");
  std::size_t compiled = 0;
  for (const fs::path& file : reference_files()) {
    INFO(file.string());
    const System s = rtt::io::load_system(file);
    for (const MaterialLibrary* lib : {&m1, &m2}) {
      std::optional<CompiledSystem> plain;
      try {
        plain = compile(s, *lib, coatings);
      } catch (const CompileError&) {
        continue;
      }
      require_same(*plain, compile(s, *lib, coatings, 0));
      System with_table = s;
      with_table.parameters.push_back(row("UNUSED", 1.0));
      with_table.parameters.push_back(row("ALSO", ParameterExpression{"2 * UNUSED"}));
      require_same(*plain, compile(with_table, *lib, coatings, 0));
      ++compiled;
      break;
    }
  }
  CHECK(compiled >= 20U);
}

TEST_CASE("configurations: errors come as CompileError", "[compile][configurations]") {
  const MaterialLibrary lib;
  // An index beyond the configurations, with and without a configurations section.
  try {
    static_cast<void>(compile(zoom(), lib, 2));
    FAIL("no CompileError");
  } catch (const CompileError& e) {
    REQUIRE(e.diagnostics().size() == 1);
    CHECK(e.diagnostics()[0].code == "config.unknown");
    CHECK(e.diagnostics()[0].location == "/configurations");
  }
  try {
    static_cast<void>(compile(zoom(20.0, 20.0), lib, 1));
    FAIL("no CompileError");
  } catch (const CompileError& e) {
    REQUIRE(e.diagnostics().size() == 1);
    CHECK(e.diagnostics()[0].code == "config.unknown");
    CHECK(e.diagnostics()[0].location.empty());
  }
  // A table error in any configuration: validate reports it, compile never throws
  // std::invalid_argument from resolve_parameters.
  System s = zoom();
  s.parameters.push_back(row("Q", ParameterExpression{"1 / (G - 5)"}));  // not finite in "tele"
  CHECK_THROWS_AS(compile(s, lib, 0), CompileError);
  try {
    static_cast<void>(compile(s, lib, 0));
  } catch (const CompileError& e) {
    CHECK_THAT(e.what(), ContainsSubstring("parameters.not_finite"));
  }
  // A binding to an unknown row.
  s = zoom();
  std::get<Element>(s.root.children[2].value).pose.position[2] = Param::bound("nobody");
  CHECK_THROWS_AS(compile(s, lib, 0), CompileError);
}

TEST_CASE("configurations: ghosts of one configuration", "[compile][configurations]") {
  const MaterialLibrary lib;
  const rtt::coating::CoatingLibrary none;
  const System s = zoom();
  const rtt::compile::GhostSystem g =
      rtt::compile::compile_with_ghosts(s, "main", lib, none, {}, 1);
  const CompiledSystem tele = compile(s, lib, none, 1);
  CHECK(g.system.configuration() == 1);
  REQUIRE(g.system.surfaces().size() == tele.surfaces().size());
  for (std::uint32_t i = 0; i < tele.surfaces().size(); ++i) {
    CHECK(same(g.system.surfaces()[i].to_global, tele.surfaces()[i].to_global));
  }
  CHECK(!g.ghosts.empty());
}
