#pragma once

/// @file relative_placement.hpp
/// Test helpers for relative placement (ADR 0028, #163), shared with the tests of the libraries
/// above rtt-model:
/// - the global transforms of a system computed with the formula used before schema 0.4,
///   parent_global * Isometry3::from_pose(position, rotation, pivot), node by node;
/// - the rewrite of an absolutely placed system to relative poses with the same geometry.

#include <cmath>
#include <cstddef>
#include <optional>
#include <stdexcept>
#include <variant>
#include <vector>

#include "rtt/math/isometry.hpp"
#include "rtt/math/units.hpp"
#include "rtt/model/model.hpp"

namespace rtt::model::test {

/// The pose formula before schema 0.4: math::Isometry3::from_pose(position, rotation, pivot),
/// written out here on purpose instead of calling to_isometry(). Throws std::invalid_argument for
/// a pose with a reference or an order other than the defaults (no such pose before 0.4).
inline math::Isometry3 old_pose(const Pose& p) {
  if (p.reference != PoseReference::Absolute || p.order != PoseOrder::TranslateFirst) {
    throw std::invalid_argument("old_pose: only absolute translate_first poses existed before 0.4");
  }
  return math::Isometry3::from_pose(
      math::Vec3(p.position[0].value, p.position[1].value, p.position[2].value),
      math::Vec3(p.rotation_deg[0].value, p.rotation_deg[1].value, p.rotation_deg[2].value),
      math::Vec3(p.pivot[0], p.pivot[1], p.pivot[2]));
}

/// Global transforms of all surfaces in tree order (pre-order), the order of
/// CompiledSystem::surfaces(), with the formula before schema 0.4: the root is old_pose(root),
/// every other node parent_global * old_pose(node).
inline std::vector<math::Isometry3> old_surface_globals(const System& s) {
  std::vector<math::Isometry3> out;
  struct Walk {
    std::vector<math::Isometry3>& out;
    void assembly(const Assembly& a, const math::Isometry3& global) {
      for (const Node& child : a.children) {
        if (const auto* sub = std::get_if<Assembly>(&child.value)) {
          assembly(*sub, global * old_pose(sub->pose));
        } else {
          const auto& e = std::get<Element>(child.value);
          const math::Isometry3 element_global = global * old_pose(e.pose);
          for (const Surface& surface : e.surfaces) {
            out.push_back(element_global * old_pose(surface.pose));
          }
        }
      }
    }
  };
  Walk{out}.assembly(s.root, old_pose(s.root.pose));
  return out;
}

/// The pose of an assembly or element node (without std::visit and a temporary visitor, whose
/// returned reference GCC 13 reports as possibly dangling).
inline const Pose& pose_of(const Node& n) {
  if (const auto* a = std::get_if<Assembly>(&n.value)) return a->pose;
  return std::get<Element>(n.value).pose;
}

inline Pose& pose_of(Node& n) {
  if (auto* a = std::get_if<Assembly>(&n.value)) return a->pose;
  return std::get<Element>(n.value).pose;
}

/// An absolutely placed system rewritten to relative poses (ADR 0028), and what was rewritten.
struct RelativeSystem {
  System system;
  std::size_t preceding = 0;     ///< poses rewritten to relative_to_preceding
  std::size_t sibling = 0;       ///< poses rewritten to relative_to_sibling
  std::size_t rotate_first = 0;  ///< of these, written with order rotate_first
  /// Sum of |t| over the rewritten poses, mm (the lever arms of the a priori error bound).
  double translation_sum = 0.0;
};

/// Rewrites every pose that can be relative (ADR 0028, point 5) to a relative pose with the same
/// global transform, computed from the transforms of old_surface_globals():
/// - surface j > 0 of an element: relative_to_preceding for odd j, relative_to_sibling for even
///   j (both refer to surface j - 1);
/// - an element or assembly with a preceding sibling at an odd index: relative_to_sibling;
///   otherwise relative_to_preceding if a surface precedes it in tree order; otherwise it stays
///   absolute (the root, surface 0 of every element, the first nodes before any surface).
/// The relative pose of node X with reference F is T = F^-1 X, written with pivot 0 and the
/// rotation as intrinsic X -> Y -> Z angles (R = Rx(a) Ry(b) Rz(c), docs/architecture.md):
/// b = atan2(R02, hypot(R00, R01)), a = atan2(-R12, R22), c = atan2(-R01, R00). Every third
/// rewritten pose uses rotate_first with Q = 0, so p_ref = R (t + p) and t = R^T T.t; the others
/// translate_first, p_ref = t + R p. Throws std::invalid_argument for a pose that is not absolute
/// translate_first in the input, and for a rotation near the gimbal lock (|cos b| < 0.1), where
/// the angles would lose accuracy.
inline RelativeSystem to_relative(const System& s) {
  struct Rewrite {
    RelativeSystem r;
    std::optional<math::Isometry3> last_surface;  // global transform, input geometry

    Pose relative(const math::Isometry3& reference, const math::Isometry3& global) {
      const math::Isometry3 t = reference.inverse() * global;
      const math::Mat3& m = t.rotation();
      const double cos_b = std::hypot(m(0, 0), m(0, 1));
      if (cos_b < 0.1) throw std::invalid_argument("to_relative: rotation near the gimbal lock");
      const double a = std::atan2(-m(1, 2), m(2, 2));
      const double b = std::atan2(m(0, 2), cos_b);
      const double c = std::atan2(-m(0, 1), m(0, 0));
      Pose p;
      p.rotation_deg = {Param(math::rad_to_deg(a)), Param(math::rad_to_deg(b)),
                        Param(math::rad_to_deg(c))};
      math::Vec3 position = t.translation();
      if ((r.preceding + r.sibling) % 3 == 2) {
        p.order = PoseOrder::RotateFirst;
        position = m.transpose() * t.translation();
        ++r.rotate_first;
      }
      p.position = {Param(position.x()), Param(position.y()), Param(position.z())};
      r.translation_sum += position.norm();
      return p;
    }

    void assembly(const Assembly& in, Assembly& out, const math::Isometry3& global) {
      std::optional<math::Isometry3> previous;  // global transform of the preceding sibling
      for (std::size_t i = 0; i < in.children.size(); ++i) {
        const Node& child = in.children[i];
        Node& rewritten = out.children[i];
        const math::Isometry3 child_global = global * old_pose(pose_of(child));
        Pose& target = pose_of(rewritten);
        if (previous && i % 2 == 1) {
          target = relative(*previous, child_global);
          target.reference = PoseReference::RelativeToSibling;
          ++r.sibling;
        } else if (last_surface) {
          target = relative(*last_surface, child_global);
          target.reference = PoseReference::RelativeToPreceding;
          ++r.preceding;
        }
        if (const auto* sub = std::get_if<Assembly>(&child.value)) {
          assembly(*sub, std::get<Assembly>(rewritten.value), child_global);
        } else {
          element(std::get<Element>(child.value), std::get<Element>(rewritten.value), child_global);
        }
        previous = child_global;
      }
    }

    void element(const Element& in, Element& out, const math::Isometry3& global) {
      for (std::size_t j = 0; j < in.surfaces.size(); ++j) {
        const math::Isometry3 surface_global = global * old_pose(in.surfaces[j].pose);
        if (j > 0) {
          Pose& target = out.surfaces[j].pose;
          target = relative(*last_surface, surface_global);
          target.reference =
              j % 2 == 1 ? PoseReference::RelativeToPreceding : PoseReference::RelativeToSibling;
          ++(j % 2 == 1 ? r.preceding : r.sibling);
        }
        last_surface = surface_global;
      }
    }
  };
  Rewrite w{{s, 0, 0, 0, 0.0}, std::nullopt};
  w.assembly(s.root, w.r.system.root, old_pose(s.root.pose));
  return w.r;
}

/// The fold mirror of ADR 0028, point 4, in vacuum at 0.5876 um, EPD 4, one field on axis:
/// mirror M at z = 50 tilted by (45 deg, 0, 0) with the plane surface M.S, the lens L
/// (CONST:1.5, L.S1 a sphere R = 40 mm, L.S2 plane 5 mm behind it) and the detector "image" with
/// IMG d mm behind L.S2. Path "main": M.S reflect, L.S1, L.S2, IMG transmit.
/// - relative: L relative_to_preceding (M.S) with rotate_first, rotation (-135 deg, 0, 0),
///   position (0, 0, 30); image relative_to_preceding (L.S2), position (0, 0, d).
/// - absolute counterpart: L at (0, 30, 50) and image at (0, 35 + d, 50), both rotated by
///   (-90 deg, 0, 0), so that their z axis is the reflected axis +y (ADR 0028, point 4).
inline System fold_mirror(bool relative, double d) {
  System s;
  s.name = relative ? "fold mirror, relative" : "fold mirror, absolute";
  s.environment.medium = "VACUUM";
  s.wavelengths = {{0.5876, 1.0, true}};
  s.aperture = {SystemApertureType::EntrancePupilDiameter, Param(4.0)};
  s.fields = {FieldType::AngleDeg, {{0.0, 0.0, 1.0}}};
  s.root.name = "root";
  const auto surface = [](const char* id, double z) {
    Surface x;
    x.id = SurfaceId(id);
    x.pose = Pose::along_z(z);
    return x;
  };
  Element mirror{"M", ElementKind::Mirror, Pose::along_z(50.0), std::nullopt, {surface("M.S", 0)}};
  mirror.pose.rotation_deg[0] = Param(45.0);
  Surface s1 = surface("L.S1", 0.0);
  s1.shape.base = Conic{Param(40.0), Param(0.0)};
  Element lens{"L", ElementKind::Lens, {}, "CONST:1.5", {s1, surface("L.S2", 5.0)}};
  Element image{"image", ElementKind::Detector, {}, std::nullopt, {surface("IMG", 0.0)}};
  if (relative) {
    lens.pose.reference = PoseReference::RelativeToPreceding;
    lens.pose.order = PoseOrder::RotateFirst;
    lens.pose.rotation_deg[0] = Param(-135.0);
    lens.pose.position[2] = Param(30.0);
    image.pose.reference = PoseReference::RelativeToPreceding;
    image.pose.position[2] = Param(d);
  } else {
    lens.pose.position = {Param(0.0), Param(30.0), Param(50.0)};
    lens.pose.rotation_deg[0] = Param(-90.0);
    image.pose.position = {Param(0.0), Param(35.0 + d), Param(50.0)};
    image.pose.rotation_deg[0] = Param(-90.0);
  }
  s.root.children = {{mirror}, {lens}, {image}};
  s.paths = {{"main",
              false,
              {{SurfaceId("M.S"), EventKind::Reflect, 0},
               {SurfaceId("L.S1"), EventKind::Refract, 0},
               {SurfaceId("L.S2"), EventKind::Refract, 0},
               {SurfaceId("IMG"), EventKind::Transmit, 0}}}};
  return s;
}

}  // namespace rtt::model::test
