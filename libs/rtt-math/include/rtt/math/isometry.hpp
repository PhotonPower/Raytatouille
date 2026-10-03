#pragma once

/// @file isometry.hpp
/// Rigid transformations between coordinate systems.
///
/// Convention (docs/architecture.md): a node is placed in its parent by
/// translation, then rotation about a pivot. Rotations are intrinsic in the
/// order X -> Y -> Z, i.e. R = Rx(a) * Ry(b) * Rz(c). Coordinates are right-handed,
/// the optical axis is +z and y is the meridional direction.

#include "rtt/math/types.hpp"

namespace rtt::math {

/// Rotation matrix for intrinsic rotations about X, then the new Y, then the new Z.
/// @param angles_rad rotation angles (a, b, c) in radian.
[[nodiscard]] Mat3 rotation_xyz_intrinsic(const Vec3& angles_rad);

/// Rigid transform p_parent = R * p_child + t. Always orthonormal with det(R) = +1.
class Isometry3 {
 public:
  /// Identity transform.
  Isometry3() = default;

  /// Builds a transform from a rotation (must be orthonormal) and a translation in mm.
  Isometry3(const Mat3& rotation, const Vec3& translation);

  /// Identity transform.
  [[nodiscard]] static Isometry3 identity() { return {}; }

  /// Placement of a child node in its parent.
  /// p_parent = position + pivot + R * (p_child - pivot), R = rotation_xyz_intrinsic(rotation).
  /// @param position     origin offset of the child in parent coordinates, mm
  /// @param rotation_deg intrinsic X->Y->Z rotation angles in degree
  /// @param pivot        rotation centre in child coordinates, mm
  [[nodiscard]] static Isometry3 from_pose(const Vec3& position,
                                           const Vec3& rotation_deg,
                                           const Vec3& pivot);

  [[nodiscard]] const Mat3& rotation() const noexcept { return rotation_; }
  [[nodiscard]] const Vec3& translation() const noexcept { return translation_; }

  /// Transforms a point (rotation and translation).
  [[nodiscard]] Vec3 apply_point(const Vec3& p) const { return rotation_ * p + translation_; }

  /// Transforms a direction or normal (rotation only).
  [[nodiscard]] Vec3 apply_vector(const Vec3& v) const { return rotation_ * v; }

  /// Inverse transform.
  [[nodiscard]] Isometry3 inverse() const;

  /// Element-wise comparison with an absolute tolerance.
  [[nodiscard]] bool is_approx(const Isometry3& other, double tol = 1e-12) const;

  /// Composition: (a * b).apply_point(p) == a.apply_point(b.apply_point(p)).
  friend Isometry3 operator*(const Isometry3& a, const Isometry3& b);

 private:
  Mat3 rotation_ = Mat3::Identity();
  Vec3 translation_ = Vec3::Zero();
};

}  // namespace rtt::math
