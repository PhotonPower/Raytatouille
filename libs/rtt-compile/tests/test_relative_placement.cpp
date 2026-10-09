// Relative placement in compile (ADR 0028, #163): global(X) = global(reference(X)) * pose(X), with
// the reference the parent (absolute), the last surface before X in tree order
// (relative_to_preceding) or the preceding sibling (relative_to_sibling).

#include <algorithm>
#include <bit>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <optional>
#include <set>
#include <string>
#include <variant>
#include <vector>

#include "relative_placement.hpp"
#include "rtt/coating/catalog.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"

using rtt::compile::compile;
using rtt::compile::CompiledSystem;
using rtt::material::MaterialLibrary;
using rtt::math::Isometry3;
using rtt::math::Vec3;
using rtt::model::Assembly;
using rtt::model::Element;
using rtt::model::ElementKind;
using rtt::model::Param;
using rtt::model::Pose;
using rtt::model::PoseOrder;
using rtt::model::PoseReference;
using rtt::model::Surface;
using rtt::model::SurfaceId;
using rtt::model::System;

namespace fs = std::filesystem;

namespace {

/// Tolerances agreed for #163 (coordinator): 1e-10 mm and 1e-12 rad. Each test that uses them
/// first checks that its a priori error bound lies below them.
constexpr double kTolMm = 1e-10;
constexpr double kTolRad = 1e-12;
constexpr double kEps = std::numeric_limits<double>::epsilon();

std::vector<fs::path> reference_files() {
  std::vector<fs::path> files;
  for (const auto& entry : fs::recursive_directory_iterator(RTT_REFERENCE_DIR)) {
    const std::string name = entry.path().filename().string();
    if (entry.is_regular_file() && name.ends_with(".rtt.json")) files.push_back(entry.path());
  }
  std::sort(files.begin(), files.end());
  return files;
}

/// The libraries the reference systems need: N-BK7 of the M0/M1 files or the M2 SCHOTT
/// catalogue, and the coating catalogues of M3.
struct Libraries {
  MaterialLibrary m1;
  MaterialLibrary m2;
  rtt::coating::CoatingLibrary coatings;
  Libraries() {
    m1.add_catalog(std::string(RTT_CATALOG_DIR) + "/schott.agf");
    m2.add_catalog(std::string(RTT_CATALOG_DIR) + "/m2/schott.agf");
    coatings.add_catalog(std::string(RTT_CATALOG_DIR) + "/coatings/demo.json");
    coatings.add_catalog(std::string(RTT_CATALOG_DIR) + "/coatings/m3.json");
  }
  /// Compiles with the first material library that resolves every material.
  [[nodiscard]] std::optional<CompiledSystem> compile_any(const System& s) const {
    for (const MaterialLibrary* lib : {&m1, &m2}) {
      try {
        return compile(s, *lib, coatings);
      } catch (const rtt::compile::CompileError&) {
      }
    }
    return std::nullopt;
  }
};

/// The feature tour as it compiles (the same changes as test_m4_acceptance.cpp): without the
/// Zernike term of A.S1 (M8) and with a bare Fresnel A.S1 instead of the coating "AR_VIS".
System compilable_tour() {
  System tour = rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m0/feature_tour.rtt.json");
  auto& a_s1 = std::get<Element>(tour.root.children[2].value).surfaces[0];
  a_s1.shape.terms.clear();
  a_s1.interaction = rtt::model::Fresnel{};
  return tour;
}

/// The reference files that use relative poses themselves (#163, #169, #167).
const std::set<std::string> kRelativeFiles = {"feature_tour.rtt.json", "zoom.rtt.json",
                                              "two_lens_gap.rtt.json"};

/// True if a pose of the system uses a field of schema 0.4 (reference or order).
bool uses_relative_fields(const System& s) {
  bool found = false;
  const auto check = [&](const Pose& p) {
    found = found || p.reference != PoseReference::Absolute || p.order != PoseOrder::TranslateFirst;
  };
  struct Walk {
    decltype(check)& f;
    void assembly(const Assembly& a) {
      f(a.pose);
      for (const auto& child : a.children) {
        if (const auto* sub = std::get_if<Assembly>(&child.value)) {
          assembly(*sub);
        } else {
          const auto& e = std::get<Element>(child.value);
          f(e.pose);
          for (const Surface& surface : e.surfaces) f(surface.pose);
        }
      }
    }
  };
  Walk{check}.assembly(s.root);
  return found;
}

bool same_bits(const Isometry3& a, const Isometry3& b) {
  for (int r = 0; r < 3; ++r) {
    if (std::bit_cast<std::uint64_t>(a.translation()(r)) !=
        std::bit_cast<std::uint64_t>(b.translation()(r))) {
      return false;
    }
    for (int c = 0; c < 3; ++c) {
      if (std::bit_cast<std::uint64_t>(a.rotation()(r, c)) !=
          std::bit_cast<std::uint64_t>(b.rotation()(r, c))) {
        return false;
      }
    }
  }
  return true;
}

double rotation_error(const Isometry3& a, const Isometry3& b) {
  return (a.rotation() - b.rotation()).cwiseAbs().maxCoeff();
}

double translation_error(const Isometry3& a, const Isometry3& b) {
  return (a.translation() - b.translation()).cwiseAbs().maxCoeff();
}

void require_near(const Vec3& actual, const Vec3& expected, double tol) {
  INFO("actual " << actual.transpose() << ", expected " << expected.transpose());
  REQUIRE((actual - expected).cwiseAbs().maxCoeff() <= tol);
}

const rtt::compile::CompiledSurface& surface(const CompiledSystem& cs, const char* id) {
  const auto i = cs.find_surface(SurfaceId(id));
  REQUIRE(i.has_value());
  return cs.surfaces()[*i];
}

Surface plane(const std::string& id, double z) {
  Surface s;
  s.id = SurfaceId(id);
  s.pose = Pose::along_z(z);
  return s;
}

/// Vacuum, 0.5876 um, EPD 2, one field, automatic path "main"; the root is filled by the test.
System bare_system() {
  System s;
  s.name = "relative placement";
  s.environment.medium = "VACUUM";
  s.wavelengths = {{0.5876, 1.0, true}};
  s.aperture = {rtt::model::SystemApertureType::EntrancePupilDiameter, Param(2.0)};
  s.fields = {rtt::model::FieldType::AngleDeg, {{0.0, 0.0, 1.0}}};
  s.root.name = "root";
  s.paths = {{"main", true, {}}};
  return s;
}

/// Element with one surface (validate: a thin element has exactly one).
Element thin(const std::string& name, const Pose& pose, Surface s) {
  return Element{name, ElementKind::ThinElement, pose, std::nullopt, {std::move(s)}};
}

/// Plate of CONST:1.5 with plane surfaces.
Element plate(const std::string& name, const Pose& pose, std::vector<Surface> surfaces) {
  return Element{name, ElementKind::Plate, pose, "CONST:1.5", std::move(surfaces)};
}

}  // namespace

TEST_CASE("relative placement: absolute poses are bitwise the formula before schema 0.4",
          "[compile][relative]") {
  // ADR 0028, point 6: for absolute translate_first poses compile keeps the expression
  // global(parent) * to_isometry(pose). The expected transforms are computed here with the old
  // formula parent_global * Isometry3::from_pose(position, rotation, pivot), node by node
  // (relative_placement.hpp), not with to_isometry, so the test also pins the absolute branch
  // of to_isometry. Every reference system without the fields of schema 0.4.
  const Libraries libs;
  std::size_t compared = 0;
  for (const fs::path& file : reference_files()) {
    INFO(file.filename().string());
    const System s =
        file.filename() == "feature_tour.rtt.json" ? compilable_tour() : rtt::io::load_system(file);
    if (uses_relative_fields(s)) {
      // Only these files show the fields of schema 0.4: the feature tour (#163: the image placed
      // relative to the dump mirror; its geometry is checked in a test of its own below) and
      // the zoom of #169 (relative placement with the parameter table, test_configurations).
      REQUIRE(kRelativeFiles.contains(file.filename().string()));
      continue;
    }
    const std::optional<CompiledSystem> cs = libs.compile_any(s);
    REQUIRE(cs.has_value());
    const std::vector<Isometry3> expected = rtt::model::test::old_surface_globals(s);
    REQUIRE(cs->surfaces().size() == expected.size());
    for (std::size_t i = 0; i < expected.size(); ++i) {
      INFO("surface " << cs->surfaces()[i].id.str());
      REQUIRE(same_bits(cs->surfaces()[i].to_global, expected[i]));
      REQUIRE(same_bits(cs->surfaces()[i].to_local, expected[i].inverse()));
    }
    ++compared;
  }
  REQUIRE(compared + kRelativeFiles.size() == reference_files().size());
}

TEST_CASE(
    "relative placement: every reference system rewritten to relative poses keeps its "
    "geometry",
    "[compile][relative]") {
  // ADR 0028, acceptance: each reference system with relative poses (thicknesses instead of
  // places, relative_placement.hpp: to_relative) has the same global transforms as the
  // original, to an a priori tolerance.
  //
  // Bound: a rewritten pose T = F^-1 X carries the rounding of the product (three-term sums),
  // of atan2 and of the conversion to degree and back, of sin, cos and Rx Ry Rz in to_isometry
  // and of the composition in compile: a few ulp each, together below 64 eps as an angle and
  // 64 eps |t| in the translation. Along a chain of N relative poses the angle errors add up,
  // at most 64 N eps, and each moves every later origin by at most that angle times the lever
  // arm, the sum of the |t| of the chain. With N the number of rewritten poses and S the sum of
  // their |t| (both bound any chain) the errors are at most 64 N eps (rotation entries) and
  // 64 eps (N S + S) (mm). The test requires these bounds to lie below the agreed tolerances
  // 1e-12 and 1e-10 mm, then the errors below the tolerances.
  const Libraries libs;
  std::size_t preceding = 0;
  std::size_t sibling = 0;
  std::size_t rotate_first = 0;
  for (const fs::path& file : reference_files()) {
    if (kRelativeFiles.contains(file.filename().string())) continue;  // relative already
    INFO(file.filename().string());
    const System s = rtt::io::load_system(file);
    const rtt::model::test::RelativeSystem r = rtt::model::test::to_relative(s);
    const auto n = static_cast<double>(r.preceding + r.sibling);
    const double bound_rad = 64.0 * n * kEps;
    const double bound_mm = 64.0 * kEps * (n + 1.0) * r.translation_sum;
    INFO("N = " << n << ", S = " << r.translation_sum << " mm, bounds " << bound_rad << " rad, "
                << bound_mm << " mm");
    REQUIRE(bound_rad <= kTolRad);
    REQUIRE(bound_mm <= kTolMm);
    const std::optional<CompiledSystem> absolute = libs.compile_any(s);
    const std::optional<CompiledSystem> relative = libs.compile_any(r.system);
    REQUIRE(absolute.has_value());
    REQUIRE(relative.has_value());
    REQUIRE(relative->surfaces().size() == absolute->surfaces().size());
    for (std::size_t i = 0; i < absolute->surfaces().size(); ++i) {
      INFO("surface " << absolute->surfaces()[i].id.str());
      const Isometry3& a = absolute->surfaces()[i].to_global;
      const Isometry3& b = relative->surfaces()[i].to_global;
      CHECK(rotation_error(a, b) <= kTolRad);
      CHECK(translation_error(a, b) <= kTolMm);
    }
    preceding += r.preceding;
    sibling += r.sibling;
    rotate_first += r.rotate_first;
  }
  // Both references and both orders occur (the rewrite is not empty).
  REQUIRE(preceding > 10);
  REQUIRE(sibling > 10);
  REQUIRE(rotate_first > 10);
}

TEST_CASE("relative placement: preceding is the last surface, sibling the sibling's own frame",
          "[compile][relative]") {
  // ADR 0028, point 2. Plate E1 at z = 10, rotated by 90 deg about z, with E1.S1 at 0 and E1.S2
  // at z = 3 (element coordinates). Rz(90 deg) (x, y, z) = (-y, x, z), Rx(90 deg) (x, y, z) =
  // (x, -z, y) (docs/architecture.md).
  // - E2 relative_to_sibling, position (1, 0, 0): the frame of E1, not its last surface E1.S2:
  //   (0, 0, 10) + Rz (1, 0, 0) = (0, 1, 10).
  // - Plate E3 relative_to_preceding, position (1, 0, 0): the last surface before it is E2.S at
  //   the origin of E2 (rotation of E1): (0, 1, 10) + Rz (1, 0, 0) = (0, 2, 10).
  // - In E3, E3.S1 is tilted by 90 deg about x; E3.S2 relative_to_sibling with position
  //   (0, 0, 2) follows it: (0, 2, 10) + Rz Rx (0, 0, 2) = (0, 2, 10) + Rz (0, -2, 0) =
  //   (2, 2, 10), z axis Rz Rx z = (1, 0, 0).
  // - E4 relative_to_preceding, position (0, 0, 4): the last surface before it is E3.S2:
  //   (2, 2, 10) + 4 (1, 0, 0) = (6, 2, 10).
  // Hand calculation; exact in binary up to the rounding of cos(90 deg) ~ 6e-17: 1e-12.
  System s = bare_system();
  Element e1 = plate("E1", Pose::along_z(10.0), {plane("E1.S1", 0.0), plane("E1.S2", 3.0)});
  e1.pose.rotation_deg[2] = Param(90.0);
  Element e2 = thin("E2", {}, plane("E2.S", 0.0));
  e2.pose.reference = PoseReference::RelativeToSibling;
  e2.pose.position[0] = Param(1.0);
  Element e3 = plate("E3", {}, {plane("E3.S1", 0.0), plane("E3.S2", 0.0)});
  e3.pose.reference = PoseReference::RelativeToPreceding;
  e3.pose.position[0] = Param(1.0);
  e3.surfaces[0].pose.rotation_deg[0] = Param(90.0);
  e3.surfaces[1].pose.reference = PoseReference::RelativeToSibling;
  e3.surfaces[1].pose.position[2] = Param(2.0);
  Element e4 = thin("E4", {}, plane("E4.S", 0.0));
  e4.pose.reference = PoseReference::RelativeToPreceding;
  e4.pose.position[2] = Param(4.0);
  s.root.children = {{e1}, {e2}, {e3}, {e4}};
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(s, lib);
  require_near(surface(cs, "E1.S2").to_global.apply_point(Vec3::Zero()), Vec3(0.0, 0.0, 13.0),
               1e-12);
  require_near(surface(cs, "E2.S").to_global.apply_point(Vec3::Zero()), Vec3(0.0, 1.0, 10.0),
               1e-12);
  require_near(surface(cs, "E3.S1").to_global.apply_point(Vec3::Zero()), Vec3(0.0, 2.0, 10.0),
               1e-12);
  require_near(surface(cs, "E3.S2").to_global.apply_point(Vec3::Zero()), Vec3(2.0, 2.0, 10.0),
               1e-12);
  // E3.S2 keeps the tilt of E3.S1: its z axis is Rz Rx z = Rz (0, -1, 0) = (1, 0, 0).
  require_near(surface(cs, "E3.S2").to_global.apply_vector(Vec3::UnitZ()), Vec3(1.0, 0.0, 0.0),
               1e-12);
  // E4 refers to E3.S2 (its own frame, tilted): origin (2, 2, 10) + 4 (1, 0, 0).
  require_near(surface(cs, "E4.S").to_global.apply_point(Vec3::Zero()), Vec3(6.0, 2.0, 10.0),
               1e-12);

  SECTION("the children of a relative assembly follow it, siblings of assemblies too") {
    // G1 at z = 20 with element A (A.S at 0); G2 relative_to_sibling (G1), position (0, 0, 5):
    // G2 at z = 25, its element B (absolute in G2, z = 1) at z = 26.
    System t = bare_system();
    Assembly g1;
    g1.name = "G1";
    g1.pose = Pose::along_z(20.0);
    g1.children = {{thin("A", {}, plane("A.S", 0.0))}};
    Assembly g2;
    g2.name = "G2";
    g2.pose.reference = PoseReference::RelativeToSibling;
    g2.pose.position[2] = Param(5.0);
    g2.children = {{thin("B", Pose::along_z(1.0), plane("B.S", 0.0))}};
    t.root.children = {{g1}, {g2}};
    const CompiledSystem ct = compile(t, lib);
    require_near(surface(ct, "B.S").to_global.apply_point(Vec3::Zero()), Vec3(0.0, 0.0, 26.0),
                 1e-12);
  }
}

TEST_CASE("relative placement: the fold mirror of ADR 0028, point 4", "[compile][relative]") {
  // ADR 0028, point 4 (numbers in the ADR): L relative_to_preceding (M.S) with rotate_first,
  // rotation (-135 deg, 0, 0), position (0, 0, 30): origin (0, 30, 50), z axis +y. L.S2 at
  // (0, 35, 50); the image d behind it at (0, 35 + d, 50). The absolute counterpart has the same
  // transforms: rotations of +-45 and 90 deg and translations below 100 mm, a few operations
  // each: below 1e-12 (13 measured 1e-14 in z in #171).
  const MaterialLibrary lib;
  for (const double d : {20.0, 25.0}) {
    INFO("d = " << d);
    const CompiledSystem rel = compile(rtt::model::test::fold_mirror(true, d), lib);
    const CompiledSystem abs = compile(rtt::model::test::fold_mirror(false, d), lib);
    require_near(surface(rel, "L.S1").to_global.apply_point(Vec3::Zero()), Vec3(0.0, 30.0, 50.0),
                 1e-12);
    require_near(surface(rel, "L.S1").to_global.apply_vector(Vec3::UnitZ()), Vec3::UnitY(), 1e-12);
    require_near(surface(rel, "L.S2").to_global.apply_point(Vec3::Zero()), Vec3(0.0, 35.0, 50.0),
                 1e-12);
    require_near(surface(rel, "IMG").to_global.apply_point(Vec3::Zero()), Vec3(0.0, 35.0 + d, 50.0),
                 1e-12);
    for (std::size_t i = 0; i < abs.surfaces().size(); ++i) {
      INFO("surface " << abs.surfaces()[i].id.str());
      CHECK(rotation_error(rel.surfaces()[i].to_global, abs.surfaces()[i].to_global) <= 1e-12);
      CHECK(translation_error(rel.surfaces()[i].to_global, abs.surfaces()[i].to_global) <= 1e-12);
    }
  }
}

TEST_CASE("relative placement: an element relative to a surface outside its assembly stays there",
          "[compile][relative]") {
  // ADR 0028, point 8 (limit): E in assembly G refers to the last surface before it, A.S
  // outside G. Moving G moves its absolute child F but not E.
  const auto make = [](double g_z) {
    System s = bare_system();
    Assembly g;
    g.name = "G";
    g.pose = Pose::along_z(g_z);
    Element e = thin("E", {}, plane("E.S", 0.0));
    e.pose.reference = PoseReference::RelativeToPreceding;
    e.pose.position[2] = Param(5.0);
    g.children = {{e}, {thin("F", Pose::along_z(30.0), plane("F.S", 0.0))}};
    s.root.children = {{thin("A", Pose::along_z(10.0), plane("A.S", 0.0))}, {g}};
    return s;
  };
  const MaterialLibrary lib;
  for (const double g_z : {20.0, 100.0}) {
    INFO("G at z = " << g_z);
    const CompiledSystem cs = compile(make(g_z), lib);
    require_near(surface(cs, "E.S").to_global.apply_point(Vec3::Zero()), Vec3(0.0, 0.0, 15.0),
                 1e-12);
    require_near(surface(cs, "F.S").to_global.apply_point(Vec3::Zero()), Vec3(0.0, 0.0, g_z + 30.0),
                 1e-12);
  }
}

TEST_CASE("relative placement: the image of the feature tour follows the dump mirror",
          "[compile][relative]") {
  // tests/reference/m0/feature_tour.rtt.json (#163): the image is placed relative to DUMP
  // (mirror at (0, 30, 70), rotated by (90 deg, 0, 0)) with rotate_first, rotation
  // (-90 deg, 0, 0) and position (0, -30, 50): p = Rx(90) Rx(-90) (t + p_image) + (0, 30, 70),
  // so the image lies at (0, 0, 120) with z axis +z, where it was absolute before. A few
  // operations on numbers below 130 mm: 1e-12.
  const System tour = compilable_tour();
  const auto& image = std::get<Element>(tour.root.children.back().value);
  REQUIRE(image.pose.reference == PoseReference::RelativeToPreceding);
  REQUIRE(image.pose.order == PoseOrder::RotateFirst);
  MaterialLibrary lib;
  lib.add_catalog(std::string(RTT_CATALOG_DIR) + "/schott.agf");
  const CompiledSystem cs = compile(tour, lib);
  require_near(surface(cs, "IMG").to_global.apply_point(Vec3::Zero()), Vec3(0.0, 0.0, 120.0),
               1e-12);
  require_near(surface(cs, "IMG").to_global.apply_vector(Vec3::UnitZ()), Vec3::UnitZ(), 1e-12);
}
