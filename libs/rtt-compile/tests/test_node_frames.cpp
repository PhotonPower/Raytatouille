// Frames of the model nodes in a compiled system (#169; ADR 0028, points 2 and 6): for every
// assembly, element and surface the global frame its pose is given in and its own global
// transform, as a GUI needs them to move a node or to turn an absolute pose into a relative one.

#include <algorithm>
#include <bit>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <map>
#include <optional>
#include <string>
#include <variant>
#include <vector>

#include "rtt/coating/catalog.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"

namespace fs = std::filesystem;
using rtt::compile::compile;
using rtt::compile::CompiledSystem;
using rtt::compile::CompileError;
using rtt::compile::NodeFrame;
using rtt::material::MaterialLibrary;
using rtt::math::Isometry3;
using rtt::math::Vec3;
using rtt::model::Assembly;
using rtt::model::Element;
using rtt::model::ElementKind;
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

std::map<std::string, NodeFrame> by_location(const CompiledSystem& cs) {
  std::map<std::string, NodeFrame> out;
  for (const NodeFrame& f : cs.node_frames()) {
    REQUIRE(out.emplace(f.location, f).second);  // each node once
  }
  return out;
}

/// Pointers and poses of every node of the model in pre-order (an element before its surfaces).
struct ModelNode {
  std::string location;
  Pose pose;
};

void collect(const Assembly& a, const std::string& loc, std::vector<ModelNode>& out) {
  out.push_back({loc, a.pose});
  for (std::size_t i = 0; i < a.children.size(); ++i) {
    const std::string child = loc + "/children/" + std::to_string(i);
    if (const auto* sub = std::get_if<Assembly>(&a.children[i].value)) {
      collect(*sub, child, out);
    } else {
      const auto& e = std::get<Element>(a.children[i].value);
      out.push_back({child, e.pose});
      for (std::size_t j = 0; j < e.surfaces.size(); ++j) {
        out.push_back({child + "/surfaces/" + std::to_string(j), e.surfaces[j].pose});
      }
    }
  }
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

Surface plane(const std::string& id, double z) {
  Surface s;
  s.id = SurfaceId(id);
  s.pose = Pose::along_z(z);
  s.aperture = rtt::model::CircularAperture{5.0, 0.0};
  return s;
}

Element thin(const std::string& name, Pose pose, Surface s) {
  return Element{name, ElementKind::ThinElement, pose, std::nullopt, {std::move(s)}};
}

}  // namespace

TEST_CASE("node frames: every node of every reference system, in pre-order", "[compile][frames]") {
  // to_global = reference * to_isometry(pose) bit for bit (the expression of compile), the
  // surfaces' frames equal CompiledSurface::to_global, and the root has the identity as
  // reference.
  MaterialLibrary m1;
  m1.add_catalog(std::string(RTT_CATALOG_DIR) + "/schott.agf");
  MaterialLibrary m2;
  m2.add_catalog(std::string(RTT_CATALOG_DIR) + "/m2/schott.agf");
  rtt::coating::CoatingLibrary coatings;
  coatings.add_catalog(std::string(RTT_CATALOG_DIR) + "/coatings/demo.json");
  coatings.add_catalog(std::string(RTT_CATALOG_DIR) + "/coatings/m3.json");
  std::size_t checked = 0;
  for (const fs::path& file : reference_files()) {
    INFO(file.string());
    const System s = rtt::io::load_system(file);
    std::optional<CompiledSystem> cs;
    for (const MaterialLibrary* lib : {&m1, &m2}) {
      try {
        cs = compile(s, *lib, coatings);
        break;
      } catch (const CompileError&) {
      }
    }
    if (!cs) continue;
    // Bound Params (tests/reference/m5/zoom.rtt.json) carry their value only after
    // resolve_parameters, as compile uses them.
    const System resolved = rtt::model::resolve_parameters(s, cs->configuration());
    std::vector<ModelNode> nodes;
    collect(resolved.root, "/root", nodes);
    REQUIRE(cs->node_frames().size() == nodes.size());
    for (std::size_t i = 0; i < nodes.size(); ++i) {
      const NodeFrame& f = cs->node_frames()[i];
      INFO(nodes[i].location);
      CHECK(f.location == nodes[i].location);
      CHECK(same(f.to_global, f.reference * rtt::model::to_isometry(nodes[i].pose)));
    }
    CHECK(same(cs->node_frames()[0].reference, Isometry3::identity()));
    const auto frames = by_location(*cs);
    for (const auto& surface : cs->surfaces()) {
      INFO(surface.id.str());
      CHECK(same(frames.at(surface.location).to_global, surface.to_global));
    }
    ++checked;
  }
  CHECK(checked >= 20U);
}

TEST_CASE("node frames: the reference of each kind of pose", "[compile][frames]") {
  // A at z = 10 (absolute); B relative_to_preceding (A.S) at 5; group G absolute at z = 40 with
  // E (absolute in G, z = 1) and F relative_to_sibling (E) at 2; within F a second surface
  // relative to the first (surface j > 0).
  System s;
  s.name = "frames";
  s.environment.medium = "VACUUM";
  s.wavelengths = {{0.5876, 1.0, true}};
  s.aperture = {rtt::model::SystemApertureType::EntrancePupilDiameter, rtt::model::Param(2.0)};
  s.fields = {rtt::model::FieldType::AngleDeg, {{0.0, 0.0, 1.0}}};
  s.root.name = "root";
  Element b = thin("B", Pose::along_z(5.0), plane("B.S", 0.0));
  b.pose.reference = PoseReference::RelativeToPreceding;
  Assembly g;
  g.name = "G";
  g.pose = Pose::along_z(40.0);
  Element f{"F",
            ElementKind::Plate,
            Pose::along_z(2.0),
            "CONST:1.5",
            {plane("F.S1", 0.0), plane("F.S2", 3.0)}};
  f.pose.reference = PoseReference::RelativeToSibling;
  f.surfaces[1].pose.reference = PoseReference::RelativeToPreceding;
  g.children = {{thin("E", Pose::along_z(1.0), plane("E.S", 0.0))}, {f}};
  s.root.children = {{thin("A", Pose::along_z(10.0), plane("A.S", 0.0))}, {b}, {g}};
  s.paths = {{"main", true, {}}};
  const CompiledSystem cs = compile(s, MaterialLibrary{});
  const auto frames = by_location(cs);
  const auto origin = [&](const char* loc, bool reference) {
    const NodeFrame& fr = frames.at(loc);
    return (reference ? fr.reference : fr.to_global).apply_point(Vec3::Zero()).z();
  };
  CHECK(origin("/root/children/0", true) == 0.0);   // A: the root
  CHECK(origin("/root/children/1", true) == 10.0);  // B: A.S
  CHECK(origin("/root/children/1", false) == 15.0);
  CHECK(origin("/root/children/2", true) == 0.0);              // G: the root
  CHECK(origin("/root/children/2/children/0", true) == 40.0);  // E: G
  CHECK(origin("/root/children/2/children/1", true) == 41.0);  // F: its sibling E
  CHECK(origin("/root/children/2/children/1", false) == 43.0);
  CHECK(origin("/root/children/2/children/1/surfaces/0", true) == 43.0);  // F.S1: F
  CHECK(origin("/root/children/2/children/1/surfaces/1", true) == 43.0);  // F.S2: F.S1
  CHECK(origin("/root/children/2/children/1/surfaces/1", false) == 46.0);
  CHECK(frames.size() == 1 + 2 + 2 + 1 + 2 + 3);  // root, A + A.S, B + B.S, G, E + E.S, F + 2
}
