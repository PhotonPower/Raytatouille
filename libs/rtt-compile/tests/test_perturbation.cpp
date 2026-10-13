// Overlay of tolerance perturbations on compile (ADR 0036, points 1b and 1c; #254): rigid
// motions of nodes and surfaces in their nominal frames, G'_s = [prod over X in ancestors(s) and
// s, outer to inner, of G_X T_X G_X^-1] * G_s with T = Tr(t) Tr(P) R Tr(-P), and material
// changes that keep every relative partial dispersion. Without overlay compile is bitwise
// unchanged.
//
// Expected values: positions after a tilt by theta about x of a point at z in the frame of the
// pivot are (0, -z sin theta, z cos theta) (R = Rx, positive theta turns +y towards +z,
// docs/architecture.md, Konventionen); the other expected values are the definitions of ADR 0036
// written out by hand. Tolerance 1e-12 mm (coordinates <= 50 mm, a few roundings of 1.1e-16
// relative each, so the error stays below 1e-13).

#include <algorithm>
#include <array>
#include <bit>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <functional>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "rtt/coating/catalog.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/compile/perturbation.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/air.hpp"
#include "rtt/material/material.hpp"
#include "rtt/math/isometry.hpp"
#include "rtt/model/model.hpp"

namespace fs = std::filesystem;
using rtt::compile::compile;
using rtt::compile::CompiledSystem;
using rtt::compile::CompileError;
using rtt::compile::MaterialChange;
using rtt::compile::Perturbation;
using rtt::compile::RigidMotion;
using rtt::material::MaterialLibrary;
using rtt::math::Isometry3;
using rtt::math::Vec3;
using rtt::model::Assembly;
using rtt::model::CircularAperture;
using rtt::model::Conic;
using rtt::model::Element;
using rtt::model::ElementKind;
using rtt::model::Event;
using rtt::model::EventKind;
using rtt::model::Param;
using rtt::model::Pose;
using rtt::model::PoseOrder;
using rtt::model::PoseReference;
using rtt::model::Surface;
using rtt::model::SurfaceId;
using rtt::model::System;

namespace {

constexpr double kTol = 1e-12;
constexpr double kTheta = 0.01;  // rad

bool same(double a, double b) {
  return std::bit_cast<std::uint64_t>(a) == std::bit_cast<std::uint64_t>(b);
}

bool same(const Vec3& a, const Vec3& b) {
  return same(a.x(), b.x()) && same(a.y(), b.y()) && same(a.z(), b.z());
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

bool near(const Vec3& a, const Vec3& b, double tol = kTol) {
  return (a - b).cwiseAbs().maxCoeff() <= tol;
}

/// Asserts that two compiled systems are bitwise the same in everything the overlay could touch:
/// surface frames and axes, node frames, media, elements and path events.
void require_same(const CompiledSystem& a, const CompiledSystem& b) {
  REQUIRE(a.surfaces().size() == b.surfaces().size());
  for (std::size_t i = 0; i < a.surfaces().size(); ++i) {
    const auto& sa = a.surfaces()[i];
    const auto& sb = b.surfaces()[i];
    INFO("surface " << sa.id.str());
    REQUIRE(sa.id == sb.id);
    CHECK(same(sa.to_global, sb.to_global));
    CHECK(same(sa.to_local, sb.to_local));
    CHECK(sa.ideal_axis.has_value() == sb.ideal_axis.has_value());
    if (sa.ideal_axis && sb.ideal_axis) CHECK(same(*sa.ideal_axis, *sb.ideal_axis));
    CHECK(sa.medium_front == sb.medium_front);
    CHECK(sa.medium_back == sb.medium_back);
    CHECK(sa.coating.has_value() == sb.coating.has_value());
    if (sa.coating && sb.coating)
      CHECK(sa.coating->substrate_medium == sb.coating->substrate_medium);
  }
  REQUIRE(a.node_frames().size() == b.node_frames().size());
  for (std::size_t i = 0; i < a.node_frames().size(); ++i) {
    INFO("node " << a.node_frames()[i].location);
    CHECK(a.node_frames()[i].location == b.node_frames()[i].location);
    CHECK(same(a.node_frames()[i].reference, b.node_frames()[i].reference));
    CHECK(same(a.node_frames()[i].to_global, b.node_frames()[i].to_global));
  }
  REQUIRE(a.media().size() == b.media().size());
  for (std::size_t m = 0; m < a.media().size(); ++m) {
    const auto& ma = a.media()[m];
    const auto& mb = b.media()[m];
    CHECK(ma.reference == mb.reference);
    REQUIRE(ma.index.size() == mb.index.size());
    for (std::size_t w = 0; w < ma.index.size(); ++w) {
      CHECK(same(ma.index[w].real(), mb.index[w].real()));
      CHECK(same(ma.index[w].imag(), mb.index[w].imag()));
    }
    CHECK(ma.optic_axis.has_value() == mb.optic_axis.has_value());
    if (ma.optic_axis && mb.optic_axis) CHECK(same(*ma.optic_axis, *mb.optic_axis));
  }
  REQUIRE(a.elements().size() == b.elements().size());
  for (std::size_t e = 0; e < a.elements().size(); ++e) {
    CHECK(a.elements()[e].media == b.elements()[e].media);
    CHECK(a.elements()[e].segmented == b.elements()[e].segmented);
  }
  REQUIRE(a.paths().size() == b.paths().size());
  for (std::size_t p = 0; p < a.paths().size(); ++p) {
    REQUIRE(a.paths()[p].events.size() == b.paths()[p].events.size());
    for (std::size_t k = 0; k < a.paths()[p].events.size(); ++k) {
      const auto& ea = a.paths()[p].events[k];
      const auto& eb = b.paths()[p].events[k];
      CHECK(ea.surface == eb.surface);
      CHECK(ea.medium_before == eb.medium_before);
      CHECK(ea.medium_after == eb.medium_after);
      CHECK(ea.medium_beyond == eb.medium_beyond);
      CHECK(ea.from_inside == eb.from_inside);
    }
  }
}

std::vector<fs::path> reference_files() {
  std::vector<fs::path> files;
  for (const auto& entry : fs::recursive_directory_iterator(RTT_REFERENCE_DIR)) {
    const std::string name = entry.path().filename().string();
    if (entry.is_regular_file() && name.ends_with(".rtt.json")) files.push_back(entry.path());
  }
  std::sort(files.begin(), files.end());
  return files;
}

/// Plane surface `id` with its pose `z` mm along the axis of its reference.
Surface plane(const std::string& id, double z, PoseReference reference = PoseReference::Absolute) {
  Surface s;
  s.id = SurfaceId(id);
  s.pose = Pose::along_z(z);
  s.pose.reference = reference;
  s.aperture = CircularAperture{5.0, 0.0};
  return s;
}

/// Lens `name` in CONST:1.5 at `pose`: plane <name>.S1 at 0, plane <name>.S2 `thickness` mm
/// behind it (relative_to_preceding).
Element lens(const std::string& name, Pose pose, double thickness) {
  return Element{name,
                 ElementKind::Lens,
                 std::move(pose),
                 "CONST:1.5",
                 {plane(name + ".S1", 0.0),
                  plane(name + ".S2", thickness, PoseReference::RelativeToPreceding)}};
}

Pose relative_z(double z) {
  Pose p = Pose::along_z(z);
  p.reference = PoseReference::RelativeToPreceding;
  return p;
}

/// Valid system in vacuum without elements: wavelength 0.5876 um, EPD 2 mm, root "system".
System empty_system() {
  System s;
  s.name = "perturbation";
  s.environment.medium = "VACUUM";
  s.wavelengths = {{0.5876, 1.0, true}};
  s.aperture = {rtt::model::SystemApertureType::EntrancePupilDiameter, Param(2.0)};
  s.fields = {rtt::model::FieldType::AngleDeg, {{0.0, 0.0, 1.0}}};
  s.root.name = "system";
  s.paths = {{"main", true, {}}};
  return s;
}

/// One lens L at z = 10: L.S1 at 0, L.S2 5 mm behind it.
System single_lens(Pose pose = Pose::along_z(10.0)) {
  System s = empty_system();
  s.root.children.push_back({lens("L", std::move(pose), 5.0)});
  return s;
}

/// Assembly A at z = 20 with L1 (at A's origin, 4 mm thick) and L2 (6 mm behind L1.S2, 3 mm
/// thick, relative_to_preceding); outside A the lens L3 10 mm behind L2.S2 (relative_to_preceding,
/// 2 mm thick). In A coordinates the vertices are at z = 0, 4, 10, 13.
/// Locations: A /root/children/0, L1 .../children/0/children/0, L2 .../children/0/children/1,
/// L3 /root/children/1.
System assembly_system(Pose a_pose = Pose::along_z(20.0)) {
  System s = empty_system();
  Assembly a{"A", std::move(a_pose), {}};
  a.children.push_back({lens("L1", Pose{}, 4.0)});
  a.children.push_back({lens("L2", relative_z(6.0), 3.0)});
  s.root.children.push_back({std::move(a)});
  s.root.children.push_back({lens("L3", relative_z(10.0), 2.0)});
  return s;
}

const MaterialLibrary& plain_materials() {
  static const MaterialLibrary library;
  return library;
}

const rtt::coating::CoatingLibrary& no_coatings() {
  static const rtt::coating::CoatingLibrary library;
  return library;
}

CompiledSystem compile_with(const System& s, const Perturbation& p) {
  return compile(s, plain_materials(), no_coatings(), 0, p);
}

CompiledSystem compile_plain(const System& s) {
  return compile(s, plain_materials(), no_coatings());
}

/// Motion with translation t and angles in rad.
RigidMotion motion(const Vec3& t, const Vec3& angles, std::optional<Vec3> pivot = std::nullopt) {
  RigidMotion m;
  m.translation_mm = t;
  m.rotation_rad = angles;
  m.pivot_mm = pivot;
  return m;
}

RigidMotion tilt_x(double theta) {
  return motion(Vec3::Zero(), Vec3(theta, 0.0, 0.0));
}

const rtt::compile::CompiledSurface& surface(const CompiledSystem& cs, const std::string& id) {
  const auto i = cs.find_surface(SurfaceId(id));
  REQUIRE(i.has_value());
  const std::uint32_t index = i.value_or(0);  // bound once after the check (clang-tidy)
  return cs.surfaces()[index];
}

Vec3 vertex(const CompiledSystem& cs, const std::string& id) {
  return surface(cs, id).to_global.translation();
}

const rtt::compile::NodeFrame& frame(const CompiledSystem& cs, const std::string& location) {
  const auto& frames = cs.node_frames();
  const auto it = std::find_if(frames.begin(), frames.end(),
                               [&](const auto& f) { return f.location == location; });
  REQUIRE(it != frames.end());
  return *it;
}

/// T = Tr(t) Tr(P) R Tr(-P) (ADR 0036, point 1b), written out for the expected values.
Isometry3 motion_transform(const Vec3& t, const Vec3& angles, const Vec3& pivot) {
  const rtt::math::Mat3 r = rtt::math::rotation_xyz_intrinsic(angles);
  return Isometry3(r, t + pivot - r * pivot);
}

/// Calls `f` with every assembly (except the root) and every element of the tree.
void for_each_node(const Assembly& a,
                   const std::function<void(const Assembly*, const Element*)>& f) {
  for (const auto& child : a.children) {
    if (const auto* sub = std::get_if<Assembly>(&child.value)) {
      f(sub, nullptr);
      for_each_node(*sub, f);
    } else {
      f(nullptr, &std::get<Element>(child.value));
    }
  }
}

/// An overlay with an entry for every node, surface and isotropic element material of `s`, all
/// identity motions (with a pivot, which does not matter then) and zero material changes.
Perturbation zero_overlay(const System& s) {
  Perturbation p;
  for_each_node(s.root, [&](const Assembly* a, const Element* e) {
    const RigidMotion zero = motion(Vec3::Zero(), Vec3::Zero(), Vec3(1.0, 2.0, 3.0));
    if (a != nullptr) {
      p.nodes[a->name] = zero;
      return;
    }
    p.nodes[e->name] = zero;
    for (const Surface& surface : e->surfaces) p.surfaces[surface.id] = zero;
    if (e->crystal) return;
    if (!e->segment_materials.empty()) {
      // A list of equal entries is one material: segment 0 only (ADR 0035, point 2).
      const auto& list = e->segment_materials;
      const bool one =
          std::adjacent_find(list.begin(), list.end(), std::not_equal_to<>()) == list.end();
      for (std::size_t k = 0; k < (one ? std::size_t{1} : list.size()); ++k) {
        p.materials[{e->name, k}] = MaterialChange{};
      }
    } else if (e->material) {
      p.materials[{e->name, std::nullopt}] = MaterialChange{};
    }
  });
  return p;
}

}  // namespace

TEST_CASE("perturbation: no overlay and identity entries are bitwise compile without overlay",
          "[compile][perturbation]") {
  // Every reference system that compiles (#254 acceptance "Byte-Dump ohne Überlagerung").
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
      require_same(*plain, compile(s, *lib, coatings, 0, Perturbation{}));
      require_same(*plain, compile(s, *lib, coatings, 0, zero_overlay(s)));
      ++compiled;
      break;
    }
  }
  CHECK(compiled >= 20U);
}

TEST_CASE("perturbation: sign of tilts against a hand calculation", "[compile][perturbation]") {
  // Acceptance 11: L at z = 10, L.S2 5 mm behind L.S1; pivot of the pose (0) at L.S1.
  const System s = single_lens();
  const double sn = std::sin(kTheta);
  const double cs = std::cos(kTheta);

  SECTION("tilt_x turns +y towards +z") {
    Perturbation p;
    p.nodes["L"] = tilt_x(kTheta);
    const CompiledSystem c = compile_with(s, p);
    CHECK(near(vertex(c, "L.S1"), Vec3(0.0, 0.0, 10.0)));
    CHECK(near(vertex(c, "L.S2"), Vec3(0.0, -5.0 * sn, 10.0 + 5.0 * cs)));
    CHECK(near(surface(c, "L.S2").to_global.rotation().col(2), Vec3(0.0, -sn, cs)));
    // The hand values (Python, math.sin/cos of 0.01): -0.049999166670833324, 14.999750002083326.
    CHECK(std::abs(vertex(c, "L.S2").y() - -0.049999166670833324) <= kTol);
    CHECK(std::abs(vertex(c, "L.S2").z() - 14.999750002083326) <= kTol);
  }
  SECTION("tilt_y turns +z towards +x") {
    Perturbation p;
    p.nodes["L"] = motion(Vec3::Zero(), Vec3(0.0, kTheta, 0.0));
    const CompiledSystem c = compile_with(s, p);
    CHECK(near(vertex(c, "L.S2"), Vec3(5.0 * sn, 0.0, 10.0 + 5.0 * cs)));
  }
  SECTION("tilt_x of a surface turns it about its vertex (its pose pivot)") {
    Perturbation p;
    p.surfaces[SurfaceId("L.S2")] = tilt_x(kTheta);
    const CompiledSystem c = compile_with(s, p);
    CHECK(near(vertex(c, "L.S2"), Vec3(0.0, 0.0, 15.0)));
    CHECK(near(surface(c, "L.S2").to_global.rotation().col(2), Vec3(0.0, -sn, cs)));
    CHECK(near(surface(c, "L.S2").to_global.rotation().col(1), Vec3(0.0, cs, sn)));
  }
  SECTION("tilt_z turns +x towards +y") {
    Perturbation p;
    p.nodes["L"] = motion(Vec3::Zero(), Vec3(0.0, 0.0, kTheta));
    const CompiledSystem c = compile_with(s, p);
    CHECK(near(surface(c, "L.S2").to_global.rotation().col(0), Vec3(cs, sn, 0.0)));
    CHECK(near(vertex(c, "L.S2"), Vec3(0.0, 0.0, 15.0)));
  }
}

TEST_CASE("perturbation: decenter acts along the axes of the nominal frame",
          "[compile][perturbation]") {
  // L turned by 90 degrees about z: its local x is the global y.
  Pose pose = Pose::along_z(10.0);
  pose.rotation_deg[2] = Param(90.0);
  const System s = single_lens(pose);
  const CompiledSystem nominal = compile_plain(s);
  Perturbation p;
  p.nodes["L"] = motion(Vec3(0.1, 0.0, 0.0), Vec3::Zero());
  const CompiledSystem c = compile_with(s, p);
  for (const char* id : {"L.S1", "L.S2"}) {
    INFO(id);
    CHECK(near(vertex(c, id) - vertex(nominal, id), Vec3(0.0, 0.1, 0.0)));
    CHECK(near(surface(c, id).to_global.rotation().col(2),
               surface(nominal, id).to_global.rotation().col(2)));
  }
}

TEST_CASE("perturbation: default and explicit pivots", "[compile][perturbation]") {
  const Vec3 angles(kTheta, -0.02, 0.005);
  SECTION("translate_first: the pose pivot P stays fixed") {
    Pose pose = Pose::along_z(10.0);
    pose.rotation_deg[1] = Param(2.0);
    pose.pivot = {0.0, 1.0, 2.0};
    const System s = single_lens(pose);
    const Isometry3 nominal = frame(compile_plain(s), "/root/children/0").to_global;
    Perturbation p;
    p.nodes["L"] = motion(Vec3::Zero(), angles);
    const Isometry3 moved = frame(compile_with(s, p), "/root/children/0").to_global;
    const Vec3 pivot(0.0, 1.0, 2.0);
    CHECK(near(moved.apply_point(pivot), nominal.apply_point(pivot)));
    CHECK(!near(moved.apply_point(Vec3::Zero()), nominal.apply_point(Vec3::Zero()), 1e-6));
  }
  SECTION("rotate_first: the node point Q - t, which the pose maps onto Q, stays fixed") {
    Pose pose = Pose::along_z(10.0);
    pose.order = PoseOrder::RotateFirst;
    pose.rotation_deg[0] = Param(3.0);
    pose.pivot = {0.0, 1.0, 2.0};  // Q, in the reference frame (the root, identity)
    const System s = single_lens(pose);
    const Isometry3 nominal = frame(compile_plain(s), "/root/children/0").to_global;
    const Vec3 node_point = Vec3(0.0, 1.0, 2.0) - Vec3(0.0, 0.0, 10.0);
    REQUIRE(near(nominal.apply_point(node_point), Vec3(0.0, 1.0, 2.0)));
    Perturbation p;
    p.nodes["L"] = motion(Vec3::Zero(), angles);
    const Isometry3 moved = frame(compile_with(s, p), "/root/children/0").to_global;
    CHECK(near(moved.apply_point(node_point), Vec3(0.0, 1.0, 2.0)));
    CHECK(!near(moved.apply_point(Vec3::Zero()), nominal.apply_point(Vec3::Zero()), 1e-6));
  }
  SECTION("an explicit pivot replaces the default") {
    Pose pose = Pose::along_z(10.0);
    pose.pivot = {0.0, 1.0, 2.0};
    const System s = single_lens(pose);
    const Isometry3 nominal = frame(compile_plain(s), "/root/children/0").to_global;
    const Vec3 pivot(0.5, -1.0, 3.0);
    Perturbation p;
    p.nodes["L"] = motion(Vec3::Zero(), angles, pivot);
    const Isometry3 moved = frame(compile_with(s, p), "/root/children/0").to_global;
    CHECK(near(moved.apply_point(pivot), nominal.apply_point(pivot)));
    CHECK(!near(moved.apply_point(Vec3(0.0, 1.0, 2.0)), nominal.apply_point(Vec3(0.0, 1.0, 2.0)),
                1e-6));
  }
}

TEST_CASE("perturbation: an assembly motion moves all children, not the nodes behind it",
          "[compile][perturbation]") {
  // Acceptance: tilt_x 0.01 rad of A about its origin (pose pivot 0). Vertices at z_A = 0, 4,
  // 10, 13 go to (0, -z sin theta, 20 + z cos theta); L2 is placed relative to L1.S2 and moves
  // with A; L3 outside A stays at its nominal place although it is placed relative to L2.S2.
  const System s = assembly_system();
  const CompiledSystem nominal = compile_plain(s);
  Perturbation p;
  p.nodes["A"] = tilt_x(kTheta);
  const CompiledSystem c = compile_with(s, p);
  const double sn = std::sin(kTheta);
  const double cs = std::cos(kTheta);
  const std::array<std::pair<const char*, double>, 4> vertices{
      {{"L1.S1", 0.0}, {"L1.S2", 4.0}, {"L2.S1", 10.0}, {"L2.S2", 13.0}}};
  for (const auto& [id, z] : vertices) {
    INFO(id);
    CHECK(near(vertex(c, id), Vec3(0.0, -z * sn, 20.0 + z * cs)));
    CHECK(near(surface(c, id).to_global.rotation().col(2), Vec3(0.0, -sn, cs)));
  }
  // Hand values of L2.S2 (Python): -0.12999783334416665, 32.999350005416645.
  CHECK(std::abs(vertex(c, "L2.S2").y() - -0.12999783334416665) <= kTol);
  CHECK(std::abs(vertex(c, "L2.S2").z() - 32.999350005416645) <= kTol);
  for (const char* id : {"L3.S1", "L3.S2"}) {
    INFO(id);
    CHECK(same(surface(c, id).to_global, surface(nominal, id).to_global));
    CHECK(same(surface(c, id).to_local, surface(nominal, id).to_local));
  }
  CHECK(same(frame(c, "/root/children/1").to_global, frame(nominal, "/root/children/1").to_global));
  CHECK(same(frame(c, "/root/children/1").reference, frame(nominal, "/root/children/1").reference));
}

TEST_CASE("perturbation: thickness moves the following nodes, shift_z does not",
          "[compile][perturbation]") {
  // Acceptance 12 (ADR 0036, point 1b, "Asymmetrie, gewollt").
  const System s = assembly_system();
  const CompiledSystem nominal = compile_plain(s);
  const Vec3 dz(0.0, 0.0, 0.05);

  SECTION("thickness: an offset on position[2] of L1.S2 (a Param, ADR 0036, point 1a)") {
    System thicker = s;
    auto& a = std::get<Assembly>(thicker.root.children[0].value);
    std::get<Element>(a.children[0].value).surfaces[1].pose.position[2] = Param(4.05);
    const CompiledSystem c = compile_plain(thicker);
    CHECK(same(surface(c, "L1.S1").to_global, surface(nominal, "L1.S1").to_global));
    for (const char* id : {"L1.S2", "L2.S1", "L2.S2", "L3.S1", "L3.S2"}) {
      INFO(id);
      CHECK(near(vertex(c, id), vertex(nominal, id) + dz));
    }
  }
  SECTION("shift_z of L1 moves L1 only") {
    Perturbation p;
    p.nodes["L1"] = motion(dz, Vec3::Zero());
    const CompiledSystem c = compile_with(s, p);
    for (const char* id : {"L1.S1", "L1.S2"}) {
      INFO(id);
      CHECK(near(vertex(c, id), vertex(nominal, id) + dz));
    }
    for (const char* id : {"L2.S1", "L2.S2", "L3.S1", "L3.S2"}) {
      INFO(id);
      CHECK(same(surface(c, id).to_global, surface(nominal, id).to_global));
    }
  }
  SECTION("decenter of the surface L1.S1 moves only that surface") {
    Perturbation p;
    p.surfaces[SurfaceId("L1.S1")] = motion(Vec3(0.1, 0.0, 0.0), Vec3::Zero());
    const CompiledSystem c = compile_with(s, p);
    CHECK(near(vertex(c, "L1.S1"), vertex(nominal, "L1.S1") + Vec3(0.1, 0.0, 0.0)));
    for (const char* id : {"L1.S2", "L2.S1", "L2.S2", "L3.S1", "L3.S2"}) {
      INFO(id);
      CHECK(same(surface(c, id).to_global, surface(nominal, id).to_global));
    }
  }
}

TEST_CASE("perturbation: a child placed relative to a surface outside the assembly moves with it",
          "[compile][perturbation]") {
  // ADR 0028, addendum #252 ("Punkt 8 für Toleranzen erledigt"): L0 at z = 0 (2 mm thick) in
  // front of the assembly A at z = 20; the first child L1 of A is placed relative_to_preceding
  // 5 mm behind L0.S2, i.e. at global z = 7, z_A = -13 in A coordinates (L1.S2 at z_A = -9).
  // A tilt of A about its origin moves L1 rigidly with A: (0, -z sin theta, 20 + z cos theta);
  // L0 outside A stays bitwise.
  System s = empty_system();
  s.root.children.push_back({lens("L0", Pose{}, 2.0)});
  Assembly a{"A", Pose::along_z(20.0), {}};
  a.children.push_back({lens("L1", relative_z(5.0), 4.0)});
  s.root.children.push_back({std::move(a)});
  const CompiledSystem nominal = compile_plain(s);
  REQUIRE(near(vertex(nominal, "L1.S1"), Vec3(0.0, 0.0, 7.0)));
  Perturbation p;
  p.nodes["A"] = tilt_x(kTheta);
  const CompiledSystem c = compile_with(s, p);
  const double sn = std::sin(kTheta);
  const double cs = std::cos(kTheta);
  for (const auto& [id, z] :
       std::array<std::pair<const char*, double>, 2>{{{"L1.S1", -13.0}, {"L1.S2", -9.0}}}) {
    INFO(id);
    CHECK(near(vertex(c, id), Vec3(0.0, -z * sn, 20.0 + z * cs)));
    CHECK(near(surface(c, id).to_global.rotation().col(2), Vec3(0.0, -sn, cs)));
  }
  for (const char* id : {"L0.S1", "L0.S2"}) {
    INFO(id);
    CHECK(same(surface(c, id).to_global, surface(nominal, id).to_global));
  }
}

TEST_CASE("perturbation: a perturbed sibling does not move the next one",
          "[compile][perturbation]") {
  // M2 is placed relative_to_sibling 10 mm behind M1 (ADR 0028, point 2); a motion of M1 leaves
  // M2 at its nominal place (placements use the nominal frames, ADR 0036, point 1b).
  System s = empty_system();
  s.root.children.push_back({lens("M1", Pose::along_z(40.0), 2.0)});
  Pose sibling = Pose::along_z(10.0);
  sibling.reference = PoseReference::RelativeToSibling;
  s.root.children.push_back({lens("M2", sibling, 2.0)});
  const CompiledSystem nominal = compile_plain(s);
  REQUIRE(near(vertex(nominal, "M2.S1"), Vec3(0.0, 0.0, 50.0)));
  Perturbation p;
  p.nodes["M1"] = motion(Vec3(0.1, 0.0, 0.0), Vec3(kTheta, 0.0, 0.0));
  const CompiledSystem c = compile_with(s, p);
  CHECK(!near(vertex(c, "M1.S2"), vertex(nominal, "M1.S2"), 1e-3));
  for (const char* id : {"M2.S1", "M2.S2"}) {
    INFO(id);
    CHECK(same(surface(c, id).to_global, surface(nominal, id).to_global));
  }
}

TEST_CASE("perturbation: a motion of the assembly equals the motion of its children",
          "[compile][perturbation]") {
  // Acceptance: A with a base rotation; children whose frames in A are pure translations by p_c
  // (L1: 0, L2: (0, 0, 10)). Then L_c^-1 T_A L_c = Tr(t) Tr(P_A - p_c) R Tr(-(P_A - p_c)) for
  // L_c = Tr(p_c): the same t and angles about the pivot P_A - p_c.
  Pose a_pose = Pose::along_z(20.0);
  a_pose.rotation_deg[0] = Param(5.0);
  const System s = assembly_system(a_pose);
  const Vec3 t(0.01, 0.02, -0.03);
  const Vec3 angles(0.002, -0.001, 0.0005);
  const Vec3 pivot_a(1.0, 2.0, 3.0);
  Perturbation assembly;
  assembly.nodes["A"] = motion(t, angles, pivot_a);
  Perturbation children;
  children.nodes["L1"] = motion(t, angles, pivot_a);
  children.nodes["L2"] = motion(t, angles, pivot_a - Vec3(0.0, 0.0, 10.0));
  const CompiledSystem ca = compile_with(s, assembly);
  const CompiledSystem cc = compile_with(s, children);
  for (const char* id : {"L1.S1", "L1.S2", "L2.S1", "L2.S2", "L3.S1", "L3.S2"}) {
    INFO(id);
    CHECK(surface(ca, id).to_global.is_approx(surface(cc, id).to_global, kTol));
    CHECK(surface(ca, id).to_local.is_approx(surface(cc, id).to_local, kTol));
  }
  CHECK(!near(vertex(ca, "L2.S2"), vertex(compile_plain(s), "L2.S2"), 1e-3));
}

TEST_CASE("perturbation: nested motions compose outer to inner in nominal frames",
          "[compile][perturbation]") {
  // G'_s = (G_A T_A G_A^-1) (G_L1 T_L1 G_L1^-1) G_s with the nominal frames; L2 only sees A.
  Pose a_pose = Pose::along_z(20.0);
  a_pose.rotation_deg[1] = Param(4.0);
  const System s = assembly_system(a_pose);
  const CompiledSystem nominal = compile_plain(s);
  const Vec3 t_a(0.02, -0.01, 0.0);
  const Vec3 r_a(0.003, 0.0, -0.002);
  const Vec3 t_l(0.0, 0.03, 0.01);
  const Vec3 r_l(-0.004, 0.006, 0.0);
  Perturbation p;
  p.nodes["A"] = motion(t_a, r_a);
  p.nodes["L1"] = motion(t_l, r_l, Vec3(0.0, 0.0, 2.0));
  const CompiledSystem c = compile_with(s, p);

  const Isometry3 g_a = frame(nominal, "/root/children/0").to_global;
  const Isometry3 g_l = frame(nominal, "/root/children/0/children/0").to_global;
  const Isometry3 o_a = g_a * motion_transform(t_a, r_a, Vec3::Zero()) * g_a.inverse();
  const Isometry3 o_l = g_l * motion_transform(t_l, r_l, Vec3(0.0, 0.0, 2.0)) * g_l.inverse();
  for (const char* id : {"L1.S1", "L1.S2"}) {
    INFO(id);
    CHECK(surface(c, id).to_global.is_approx(o_a * o_l * surface(nominal, id).to_global, kTol));
  }
  for (const char* id : {"L2.S1", "L2.S2"}) {
    INFO(id);
    CHECK(surface(c, id).to_global.is_approx(o_a * surface(nominal, id).to_global, kTol));
  }
}

TEST_CASE("perturbation: node frames are the perturbed frames", "[compile][perturbation]") {
  // Decision E1 (#254): reference' = O * reference, to_global' = O G T, so to_global' =
  // reference' * to_isometry(pose) * T, and the frame of a surface is bitwise its to_global.
  const System s = assembly_system();
  const Vec3 t(0.0, 0.01, 0.0);
  const Vec3 angles(0.002, 0.0, 0.001);
  Perturbation p;
  p.nodes["A"] = tilt_x(kTheta);
  p.nodes["L2"] = motion(t, angles);
  const CompiledSystem c = compile_with(s, p);
  for (const auto& surface : c.surfaces()) {
    INFO(surface.id.str());
    CHECK(same(frame(c, surface.location).to_global, surface.to_global));
  }
  const auto& a = std::get<Assembly>(s.root.children[0].value);
  const auto& l1 = std::get<Element>(a.children[0].value);
  const auto& l2 = std::get<Element>(a.children[1].value);
  const rtt::compile::NodeFrame f1 = frame(c, "/root/children/0/children/0");
  CHECK(f1.to_global.is_approx(f1.reference * rtt::model::to_isometry(l1.pose), kTol));
  const rtt::compile::NodeFrame f2 = frame(c, "/root/children/0/children/1");
  CHECK(f2.to_global.is_approx(
      f2.reference * rtt::model::to_isometry(l2.pose) * motion_transform(t, angles, Vec3::Zero()),
      kTol));
}

TEST_CASE("perturbation: crystal and polarizer axes turn with the element",
          "[compile][perturbation]") {
  System s = empty_system();
  Element crystal{"C",
                  ElementKind::Plate,
                  Pose::along_z(10.0),
                  std::nullopt,
                  {plane("C.S1", 0.0), plane("C.S2", 5.0, PoseReference::RelativeToPreceding)}};
  crystal.crystal = rtt::model::CrystalMaterial{"CONST:1.6584", "CONST:1.4864"};
  crystal.optic_axis = std::array<double, 3>{0.0, 0.0, 1.0};
  Surface polarizer = plane("P.S1", 0.0);
  polarizer.interaction = rtt::model::IdealPolarizer{{1.0, 0.0, 0.0}, 0.0};
  s.root.children.push_back({crystal});
  s.root.children.push_back(
      {Element{"P", ElementKind::ThinElement, relative_z(5.0), std::nullopt, {polarizer}}});
  s.paths = {{"main",
              false,
              {Event{SurfaceId("C.S1"), EventKind::Ordinary, 0},
               Event{SurfaceId("C.S2"), EventKind::Refract, 0},
               Event{SurfaceId("P.S1"), EventKind::Transmit, 0}}}};
  Perturbation p;
  p.nodes["C"] = tilt_x(kTheta);
  p.nodes["P"] = motion(Vec3::Zero(), Vec3(0.0, kTheta, 0.0));
  const CompiledSystem c = compile_with(s, p);
  const auto crystal_medium = c.paths()[0].events[0].medium_after;
  const auto& optic_axis = c.media()[crystal_medium].optic_axis;
  REQUIRE(optic_axis.has_value());
  const Vec3 crystal_axis = optic_axis.value_or(Vec3::Zero());  // bound once (clang-tidy)
  CHECK(near(crystal_axis, Vec3(0.0, -std::sin(kTheta), std::cos(kTheta))));
  // By value: surface() gets a temporary string (-Wdangling-reference, AGENTS.md).
  const auto ideal_axis = surface(c, "P.S1").ideal_axis;
  REQUIRE(ideal_axis.has_value());
  const Vec3 polarizer_axis = ideal_axis.value_or(Vec3::Zero());  // bound once (clang-tidy)
  CHECK(near(polarizer_axis, Vec3(std::cos(kTheta), 0.0, -std::sin(kTheta))));
}

namespace {

/// Lens `name` of `glass` with radii R1 = 50 and R2 = -50 mm (K = c1 - c2 = 0.04 / mm).
Element biconvex(const std::string& name, Pose pose, const std::string& glass) {
  Element e = lens(name, std::move(pose), 4.0);
  e.material = glass;
  e.surfaces[0].shape.base = Conic{Param(50.0), Param(0.0)};
  e.surfaces[1].shape.base = Conic{Param(-50.0), Param(0.0)};
  return e;
}

/// N-BK7 lenses L (z = 10) and M (z = 30) in air; system wavelengths d, F, C (vacuum
/// wavelengths of the lines, as compile uses them) and 0.5, 0.7 um.
System glass_system() {
  System s = empty_system();
  s.environment.medium = "AIR";
  using rtt::compile::line_vacuum_wavelength_um;
  s.wavelengths = {{line_vacuum_wavelength_um(rtt::compile::kLineDAirUm), 1.0, true},
                   {line_vacuum_wavelength_um(rtt::compile::kLineFAirUm), 1.0, false},
                   {line_vacuum_wavelength_um(rtt::compile::kLineCAirUm), 1.0, false},
                   {0.5, 1.0, false},
                   {0.7, 1.0, false}};
  s.root.children.push_back({biconvex("L", Pose::along_z(10.0), "SCHOTT:N-BK7")});
  s.root.children.push_back({biconvex("M", Pose::along_z(30.0), "SCHOTT:N-BK7")});
  return s;
}

const MaterialLibrary& schott() {
  static const auto library = [] {
    auto m = std::make_unique<MaterialLibrary>();  // holds a mutex, cannot be moved
    m->add_catalog(std::string(RTT_CATALOG_DIR) + "/schott.agf");
    return m;
  }();
  return *library;
}

/// Real index of the medium of element `e` at system wavelength `w`.
double n_of(const CompiledSystem& c, std::size_t e, std::size_t w) {
  return c.media()[c.elements()[e].media[0]].index[w].real();
}

}  // namespace

TEST_CASE("perturbation: the lines are vacuum wavelengths of standard air",
          "[compile][perturbation]") {
  // lambda_air = lambda_vac / n_air(lambda_vac, 20 degC, 1 atm) (Ciddor, air.hpp).
  for (const double air :
       {rtt::compile::kLineDAirUm, rtt::compile::kLineFAirUm, rtt::compile::kLineCAirUm}) {
    const double vac = rtt::compile::line_vacuum_wavelength_um(air);
    CHECK(vac > air);
    CHECK(std::abs(vac / rtt::material::ciddor_air_index(vac, 20.0, 1.0) - air) <= 1e-15);
  }
}

TEST_CASE("perturbation: material changes keep the partial dispersions",
          "[compile][perturbation]") {
  // ADR 0036, point 1c: n'_d = n_d + dn exactly (one rounding), nu'_d = nu_d + dnu, every
  // P_x,y = (n_x - n_y) / (n_F - n_C) unchanged, kappa unchanged; the other N-BK7 lens keeps
  // the nominal medium.
  const System s = glass_system();
  const CompiledSystem nominal = compile(s, schott(), no_coatings());
  constexpr double kDn = 1e-3;
  constexpr double kDv = -0.5;
  Perturbation p;
  p.materials[{"L", std::nullopt}] = {kDn, kDv};
  const CompiledSystem c = compile(s, schott(), no_coatings(), 0, p);

  const double nd = n_of(nominal, 0, 0);
  const double nf = n_of(nominal, 0, 1);
  const double nc = n_of(nominal, 0, 2);
  const double nu = (nd - 1.0) / (nf - nc);
  CHECK(same(n_of(c, 0, 0), nd + kDn));
  const double nu_new = (n_of(c, 0, 0) - 1.0) / (n_of(c, 0, 1) - n_of(c, 0, 2));
  CHECK(std::abs(nu_new / (nu + kDv) - 1.0) <= 1e-12);
  for (const std::size_t w : {3U, 4U}) {
    INFO("wavelength " << w);
    const double p_nominal = (n_of(nominal, 0, w) - nd) / (nf - nc);
    const double p_new = (n_of(c, 0, w) - n_of(c, 0, 0)) / (n_of(c, 0, 1) - n_of(c, 0, 2));
    CHECK(std::abs(p_new - p_nominal) <= 1e-12);
  }
  for (std::size_t w = 0; w < 5; ++w) {
    CHECK(same(c.media()[c.elements()[0].media[0]].index[w].imag(),
               nominal.media()[nominal.elements()[0].media[0]].index[w].imag()));
  }
  // M keeps the nominal (shared) medium bitwise; L has a medium of its own.
  CHECK(c.elements()[1].media == nominal.elements()[1].media);
  CHECK(c.elements()[0].media[0] != c.elements()[1].media[0]);
  for (std::size_t w = 0; w < 5; ++w) CHECK(same(n_of(c, 1, w), n_of(nominal, 1, w)));
  CHECK(c.media()[c.elements()[0].media[0]].reference == "SCHOTT:N-BK7");
}

TEST_CASE("perturbation: thin-lens focal length and colour follow dn and dnu",
          "[compile][perturbation]") {
  // Thin lens in air: f(lambda) = 1 / ((n - 1) K), K = c1 - c2 = 0.04 / mm, from the compiled
  // indices (the paraxial trace of a perturbed system belongs to #262).
  const System s = glass_system();
  const CompiledSystem nominal = compile(s, schott(), no_coatings());
  constexpr double kK = 1.0 / 50.0 + 1.0 / 50.0;
  const auto f = [&](const CompiledSystem& c, std::size_t w) {
    return 1.0 / ((n_of(c, 0, w) - 1.0) * kK);
  };
  const double nd = n_of(nominal, 0, 0);
  const double nf = n_of(nominal, 0, 1);
  const double nc = n_of(nominal, 0, 2);

  SECTION("dn: dEFL/EFL = -dn/(n_d - 1) to first order") {
    constexpr double kDn = 2e-3;
    Perturbation p;
    p.materials[{"L", std::nullopt}] = {kDn, 0.0};
    const CompiledSystem c = compile(s, schott(), no_coatings(), 0, p);
    // Exact: f'/f = (n_d - 1)/(n_d - 1 + dn); first order -x with x = dn/(n_d - 1), error
    // x^2 / (1 + x) <= x^2.
    CHECK(std::abs(f(c, 0) / f(nominal, 0) - (nd - 1.0) / (nd - 1.0 + kDn)) <= 1e-12);
    const double x = kDn / (nd - 1.0);
    CHECK(std::abs((f(c, 0) - f(nominal, 0)) / f(nominal, 0) + x) <= x * x);
  }
  SECTION("dnu: d(f_F - f_C) = f dnu / nu^2 to first order") {
    constexpr double kDv = -0.5;
    Perturbation p;
    p.materials[{"L", std::nullopt}] = {0.0, kDv};
    const CompiledSystem c = compile(s, schott(), no_coatings(), 0, p);
    const double nu = (nd - 1.0) / (nf - nc);
    const double predicted = f(nominal, 0) * kDv / (nu * nu);
    const double actual = (f(c, 1) - f(c, 2)) - (f(nominal, 1) - f(nominal, 2));
    // f_F - f_C = -D / (K (n_F - 1)(n_C - 1)). With s = D'/D = nu / (nu + dnu) the change is
    // f dnu/nu^2 times ((n_d - 1)^2 / ((n_F - 1)(n_C - 1))) (nu / (nu + dnu)) (1 - w) up to
    // O((dnu/nu)^2), w = (n_F - n_d)/(n_F - 1) + (n_C - n_d)/(n_C - 1). Each factor differs from
    // 1 by at most the bound below (factor 1.5 and 2 for the neglected higher orders); for
    // N-BK7 and dnu = -0.5 it is about 3.5 %, the actual deviation about 0.4 %.
    const double a = (nf - 1.0) * (nc - 1.0);
    const double w = (nf - nd) / (nf - 1.0) + (nc - nd) / (nc - 1.0);
    const double bound = 1.5 * std::abs(kDv) / nu +
                         2.0 * std::abs((nd - 1.0) * (nd - 1.0) - a) / a + 2.0 * std::abs(w);
    REQUIRE(bound < 0.05);
    CHECK(std::abs(actual / predicted - 1.0) <= bound);
    // Sharper: against the leading term itself (review of #266). With eps = s - 1 =
    // -dnu/(nu + dnu), alpha = (n_F - n_d)/(n_F - 1), beta = (n_C - n_d)/(n_C - 1) and
    // w = alpha + beta the exact change is f dnu/nu^2 times (n_d - 1)^2/a * nu/(nu + dnu) *
    // [(1 - w) - eps alpha beta] / ((1 + eps alpha)(1 + eps beta)); the denominator is
    // 1 + eps w + O(eps^2), so the relative rest of the leading term is about -eps w (4.5e-5 for
    // N-BK7). Bound: twice that.
    const double leading = (nd - 1.0) * (nd - 1.0) / a * nu / (nu + kDv) * (1.0 - w);
    const double eps = -kDv / (nu + kDv);
    CHECK(std::abs(actual / predicted - leading) <= 2.0 * std::abs(eps * w));
  }
}

TEST_CASE("perturbation: CONST, segments and the structure of an element",
          "[compile][perturbation]") {
  SECTION("CONST: n' = n + dn at every wavelength") {
    System s = single_lens();
    s.wavelengths = {{0.5876, 1.0, true}, {0.4, 1.0, false}};
    Perturbation p;
    p.materials[{"L", std::nullopt}] = {2e-3, 0.0};
    const CompiledSystem c = compile_with(s, p);
    for (std::size_t w = 0; w < 2; ++w) CHECK(same(n_of(c, 0, w), 1.5 + 2e-3));
  }
  SECTION("cemented doublet: only the perturbed segment changes") {
    // At the d line the correction (D'/D - 1)(n - n_d) is exactly 0, so n'_d = n_d + dn.
    System s = empty_system();
    s.environment.medium = "AIR";
    s.wavelengths = {
        {rtt::compile::line_vacuum_wavelength_um(rtt::compile::kLineDAirUm), 1.0, true}};
    Element doublet = biconvex("D", Pose::along_z(10.0), "SCHOTT:N-BK7");
    doublet.material.reset();
    doublet.surfaces.push_back(plane("D.S3", 2.0, PoseReference::RelativeToPreceding));
    doublet.segment_materials = {"SCHOTT:N-BK7", "SCHOTT:F2"};
    s.root.children.push_back({doublet});
    const CompiledSystem nominal = compile(s, schott(), no_coatings());
    Perturbation p;
    p.materials[{"D", std::size_t{1}}] = {1e-3, 0.0};
    const CompiledSystem c = compile(s, schott(), no_coatings(), 0, p);
    CHECK(c.elements()[0].media[0] == nominal.elements()[0].media[0]);
    CHECK(c.elements()[0].media[1] != nominal.elements()[0].media[1]);
    const auto& f2 = c.media()[c.elements()[0].media[1]];
    const auto& f2_nominal = nominal.media()[nominal.elements()[0].media[1]];
    CHECK(f2.reference == "SCHOTT:F2");
    CHECK(same(f2.index[0].real(), f2_nominal.index[0].real() + 1e-3));
  }
  SECTION("the structure comes from the nominal media") {
    // Segments 0 and 1 are the same material, so the inner surface T.S2 may be entered from
    // outside (ADR 0017). After a change of segment 0 alone it still may: the path rules come
    // from the nominal media, and the ray enters the unchanged segment 1.
    System s = empty_system();
    Element lens4{"T",
                  ElementKind::Lens,
                  Pose::along_z(10.0),
                  std::nullopt,
                  {plane("T.S1", 0.0), plane("T.S2", 2.0, PoseReference::RelativeToPreceding),
                   plane("T.S3", 2.0, PoseReference::RelativeToPreceding),
                   plane("T.S4", 2.0, PoseReference::RelativeToPreceding)}};
    lens4.segment_materials = {"CONST:1.5", "CONST:1.5", "CONST:1.7"};
    s.root.children.push_back({lens4});
    s.paths = {{"main",
                false,
                {Event{SurfaceId("T.S2"), EventKind::Refract, 0},
                 Event{SurfaceId("T.S3"), EventKind::Refract, 0},
                 Event{SurfaceId("T.S4"), EventKind::Refract, 0}}}};
    const CompiledSystem nominal = compile_plain(s);
    Perturbation p;
    p.materials[{"T", std::size_t{0}}] = {1e-3, 0.0};
    const CompiledSystem c = compile_with(s, p);
    CHECK(c.elements()[0].media[0] != c.elements()[0].media[1]);
    CHECK(same(n_of(c, 0, 0), 1.5 + 1e-3));
    CHECK(c.paths()[0].events[0].medium_after == c.elements()[0].media[1]);
    CHECK(c.paths()[0].events[0].medium_after == nominal.paths()[0].events[0].medium_after);
  }
  SECTION("a list of equal materials is one material: segment 0 changes all segments") {
    // Decision for #254 (ADR 0035, point 2, addendum): like the shorthand.
    for (const ElementKind kind : {ElementKind::Lens, ElementKind::Plate}) {
      System s = empty_system();
      Element body{"B",
                   kind,
                   Pose::along_z(10.0),
                   std::nullopt,
                   {plane("B.S1", 0.0), plane("B.S2", 2.0, PoseReference::RelativeToPreceding),
                    plane("B.S3", 2.0, PoseReference::RelativeToPreceding)}};
      body.segment_materials = {"CONST:1.5", "CONST:1.5"};
      s.root.children.push_back({body});
      s.paths = {{"main",
                  false,
                  {Event{SurfaceId("B.S1"), EventKind::Refract, 0},
                   Event{SurfaceId("B.S3"), EventKind::Refract, 0}}}};
      const CompiledSystem nominal = compile_plain(s);
      Perturbation p;
      p.materials[{"B", std::size_t{0}}] = {1e-3, 0.0};
      const CompiledSystem c = compile_with(s, p);
      CHECK(c.elements()[0].media[0] == c.elements()[0].media[1]);
      CHECK(c.elements()[0].segmented == nominal.elements()[0].segmented);
      CHECK(same(n_of(c, 0, 0), 1.5 + 1e-3));
      CHECK(c.paths()[0].events[0].medium_after == c.elements()[0].media[0]);
      Perturbation second;
      second.materials[{"B", std::size_t{1}}] = {1e-3, 0.0};
      CHECK_THROWS_AS((void)compile_with(s, second), std::invalid_argument);
    }
  }
  SECTION("a plate of one material stays one body") {
    System s = empty_system();
    Element plate{"P",
                  ElementKind::Plate,
                  Pose::along_z(10.0),
                  "CONST:1.5",
                  {plane("P.S1", 0.0), plane("P.S2", 2.0, PoseReference::RelativeToPreceding)}};
    s.root.children.push_back({plate});
    Perturbation p;
    p.materials[{"P", std::nullopt}] = {1e-3, 0.0};
    const CompiledSystem c = compile_with(s, p);
    CHECK_FALSE(c.elements()[0].segmented);
    CHECK(same(n_of(c, 0, 0), 1.5 + 1e-3));
  }
}

namespace {

/// The diagnostics of a CompileError of compile(s, materials, overlay `p`); fails if it compiles.
std::vector<rtt::model::Diagnostic> errors(const System& s,
                                           const Perturbation& p,
                                           const MaterialLibrary& lib = plain_materials()) {
  try {
    (void)compile(s, lib, no_coatings(), 0, p);
  } catch (const CompileError& e) {
    return e.diagnostics();
  }
  FAIL("compile() did not throw");
  return {};
}

bool only(const std::vector<rtt::model::Diagnostic>& d,
          const std::string& code,
          const std::string& location) {
  INFO((d.empty() ? std::string("no diagnostics") : rtt::model::to_string(d.front())));
  return d.size() == 1 && d[0].code == code && d[0].location == location;
}

Perturbation material(const std::string& element,
                      std::optional<std::size_t> segment,
                      double dn,
                      double dv) {
  Perturbation p;
  p.materials[{element, segment}] = {dn, dv};
  return p;
}

}  // namespace

TEST_CASE("perturbation: material errors come as CompileError at the material",
          "[compile][perturbation]") {
  SECTION("Abbe number at CONST") {
    CHECK(only(errors(single_lens(), material("L", std::nullopt, 0.0, 1.0)),
               "tolerances.abbe_on_const", "/root/children/0/material"));
  }
  SECTION("n_d + dn <= 1 and nu_d + dnu <= 0") {
    CHECK(only(errors(single_lens(), material("L", std::nullopt, -0.5, 0.0)),
               "tolerances.limit_invalid", "/root/children/0/material"));
    CHECK(only(errors(glass_system(), material("L", std::nullopt, 0.0, -70.0), schott()),
               "tolerances.limit_invalid", "/root/children/0/material"));
  }
  SECTION("complex index, crystal, no material") {
    System metal = single_lens();
    std::get<Element>(metal.root.children[0].value).material = "CONST:1.5,0.1";
    CHECK(only(errors(metal, material("L", std::nullopt, 1e-3, 0.0)),
               "tolerances.material_unsupported", "/root/children/0/material"));
    System crystal = single_lens();
    auto& c = std::get<Element>(crystal.root.children[0].value);
    c.kind = ElementKind::Plate;
    c.material.reset();
    c.crystal = rtt::model::CrystalMaterial{"CONST:1.6584", "CONST:1.4864"};
    c.optic_axis = std::array<double, 3>{0.0, 0.0, 1.0};
    crystal.paths = {{"main",
                      false,
                      {Event{SurfaceId("L.S1"), EventKind::Ordinary, 0},
                       Event{SurfaceId("L.S2"), EventKind::Refract, 0}}}};
    CHECK(only(errors(crystal, material("L", std::nullopt, 1e-3, 0.0)),
               "tolerances.material_unsupported", "/root/children/0/material"));
    System mirror = empty_system();  // a mirror without substrate: no material at all
    mirror.root.children.push_back({Element{
        "M", ElementKind::Mirror, Pose::along_z(10.0), std::nullopt, {plane("M.S1", 0.0)}}});
    CHECK(only(errors(mirror, material("M", std::nullopt, 1e-3, 0.0)),
               "tolerances.material_unsupported", "/root/children/0"));
  }
  SECTION("a segment of a doublet") {
    System s = single_lens();
    auto& e = std::get<Element>(s.root.children[0].value);
    e.material.reset();
    e.surfaces.push_back(plane("L.S3", 2.0, PoseReference::RelativeToPreceding));
    e.segment_materials = {"CONST:1.5", "CONST:1.7"};
    CHECK(only(errors(s, material("L", 1, 0.0, 1.0)), "tolerances.abbe_on_const",
               "/root/children/0/material/1"));
  }
}

TEST_CASE("perturbation: entries without target are invalid arguments", "[compile][perturbation]") {
  const System s = assembly_system();
  const auto throws = [&](const Perturbation& p) {
    CHECK_THROWS_AS((void)compile_with(s, p), std::invalid_argument);
  };
  Perturbation p;
  p.nodes["NONE"] = tilt_x(kTheta);
  throws(p);
  p = {};
  p.nodes["system"] = tilt_x(kTheta);  // the root
  throws(p);
  p = {};
  p.surfaces[SurfaceId("NONE")] = tilt_x(kTheta);
  throws(p);
  throws(material("NONE", std::nullopt, 1e-3, 0.0));
  throws(material("L1", 0, 1e-3, 0.0));  // one material: no segment
  p = {};
  p.nodes["L1"] = motion(Vec3(std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0), Vec3::Zero());
  throws(p);
  p = {};
  p.nodes["L1"] =
      motion(Vec3::Zero(), Vec3::Zero(), Vec3(0.0, std::numeric_limits<double>::infinity(), 0.0));
  throws(p);
  throws(material("L1", std::nullopt, std::numeric_limits<double>::quiet_NaN(), 0.0));

  System doublet = single_lens();
  auto& e = std::get<Element>(doublet.root.children[0].value);
  e.material.reset();
  e.surfaces.push_back(plane("L.S3", 2.0, PoseReference::RelativeToPreceding));
  e.segment_materials = {"CONST:1.5", "CONST:1.7"};
  CHECK_THROWS_AS((void)compile_with(doublet, material("L", std::nullopt, 1e-3, 0.0)),
                  std::invalid_argument);
  CHECK_THROWS_AS((void)compile_with(doublet, material("L", 2, 1e-3, 0.0)), std::invalid_argument);
}

namespace {

/// Constant index 1.6 on [0.55, 2] um: the line F (0.486 um) is outside.
class RedOnly final : public rtt::material::Material {
 public:
  [[nodiscard]] rtt::math::Complex index(double, double, double) const override {
    return {1.6, 0.0};
  }
  [[nodiscard]] std::optional<rtt::material::WavelengthRange> wavelength_range_um() const override {
    return rtt::material::WavelengthRange{0.55, 2.0};
  }
};

}  // namespace

TEST_CASE("perturbation: indices at the Fraunhofer lines and nu_d", "[compile][perturbation]") {
  // Shared by the material overlay and the class limits of abbe in rtt-tolerance (#260).
  using rtt::compile::fraunhofer_indices;
  using rtt::compile::FraunhoferIndices;
  using rtt::compile::line_vacuum_wavelength_um;
  MaterialLibrary lib;
  lib.add_catalog(std::string(RTT_CATALOG_DIR) + "/schott.agf");
  const double t = 25.0;
  const double p = 1.0;

  SECTION("a catalog glass: the material at the vacuum wavelengths of the lines") {
    const auto glass = lib.resolve("SCHOTT:N-BK7");
    const FraunhoferIndices f = fraunhofer_indices(*glass, t, p);
    REQUIRE(f.status == FraunhoferIndices::Status::Ok);
    const auto at = [&](double air_um) {
      return glass->index(line_vacuum_wavelength_um(air_um), t, p).real();
    };
    const double nd = at(rtt::compile::kLineDAirUm);
    const double nf = at(rtt::compile::kLineFAirUm);
    const double nc = at(rtt::compile::kLineCAirUm);
    CHECK(same(f.nd, nd));
    CHECK(same(f.nf, nf));
    CHECK(same(f.nc, nc));
    const std::optional<double> abbe = f.abbe_number();
    REQUIRE(abbe.has_value());
    const double vd = abbe.value_or(0.0);  // bound once after the check (clang-tidy)
    CHECK(same(vd, (nd - 1.0) / (nf - nc)));
    // Content guard: V(d) = 64.17 of the NM line of N-BK7 in tests/catalogs/schott.agf (relative
    // to air at the catalogue temperature; here absolute at 25 degC, which differs by < 0.1).
    CHECK(std::abs(vd - 64.17) < 0.1);
  }
  SECTION("a real CONST material has no nu_d") {
    const FraunhoferIndices f = fraunhofer_indices(*lib.resolve("CONST:1.5"), t, p);
    REQUIRE(f.status == FraunhoferIndices::Status::Ok);
    CHECK(f.nd == 1.5);
    CHECK(f.dispersion() == 0.0);
    CHECK_FALSE(f.abbe_number().has_value());
  }
  SECTION("a metal: complex at the line d") {
    const FraunhoferIndices f = fraunhofer_indices(*lib.resolve("CONST:0.2,3.0"), t, p);
    CHECK(f.status == FraunhoferIndices::Status::Complex);
    CHECK(f.line == 0U);
  }
  SECTION("a line outside the range: F") {
    const RedOnly red;
    const FraunhoferIndices f = fraunhofer_indices(red, t, p);
    CHECK(f.status == FraunhoferIndices::Status::OutOfRange);
    CHECK(f.line == 1U);
    CHECK(f.wavelength_um == line_vacuum_wavelength_um(rtt::compile::kLineFAirUm));
  }
}
