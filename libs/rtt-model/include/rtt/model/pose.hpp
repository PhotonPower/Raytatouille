#pragma once

/// @file pose.hpp
/// Placement of a node (assembly, element, surface) relative to its reference (ADR 0028).

#include <array>
#include <cstdint>

#include "rtt/math/isometry.hpp"
#include "rtt/model/param.hpp"

namespace rtt::model {

/// What a pose is relative to (ADR 0028, point 2): global(X) = global(reference) * pose(X).
enum class PoseReference : std::uint8_t {
  Absolute,             ///< the parent node (the only kind before schema 0.4)
  RelativeToPreceding,  ///< the last surface before the node in tree order (pre-order)
  RelativeToSibling,    ///< the preceding sibling in the same assembly or element
};

/// Order of translation and rotation (ADR 0028, point 3).
enum class PoseOrder : std::uint8_t {
  /// p_ref = t + P + R (p - P): shift in the reference frame, then rotate about the pivot P,
  /// given in the coordinates of the node (the only order before schema 0.4).
  TranslateFirst,
  /// p_ref = Q + R (t + p - Q): rotate the reference frame about the pivot Q, given in the
  /// reference frame, then shift along the rotated axes (coordinate break).
  RotateFirst,
};

/// Placement of a node: translation t = `position` in mm, intrinsic X->Y->Z rotation R from
/// `rotation_deg` in degree (R = Rx Ry Rz), and a pivot in mm; `order` decides the formula and
/// the frame of the pivot (see PoseOrder), `reference` the frame the pose is given in.
struct Pose {
  std::array<Param, 3> position{};      ///< mm, along the axes of PoseOrder
  std::array<Param, 3> rotation_deg{};  ///< degree, intrinsic X -> Y -> Z
  /// mm; in the node's coordinates for TranslateFirst, in the reference frame for RotateFirst.
  /// Tolerances tilt about it.
  std::array<double, 3> pivot{};
  PoseReference reference = PoseReference::Absolute;
  PoseOrder order = PoseOrder::TranslateFirst;

  /// Convenience constructor for the common case of a pure shift along the local z axis.
  [[nodiscard]] static Pose along_z(double z_mm);

  [[nodiscard]] bool is_identity() const;

  bool operator==(const Pose&) const = default;
};

/// Evaluates the current parameter values into a rigid transform from the node's coordinates to
/// those of its reference. TranslateFirst calls math::Isometry3::from_pose unchanged (bitwise as
/// before schema 0.4); RotateFirst builds Isometry3(R, Q + R (t - Q)). Bound Params are read with
/// their current `value` (resolved before, ADR 0029).
[[nodiscard]] math::Isometry3 to_isometry(const Pose& pose);

}  // namespace rtt::model
