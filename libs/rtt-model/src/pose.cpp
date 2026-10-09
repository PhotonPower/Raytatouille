#include "rtt/model/pose.hpp"

#include "rtt/math/units.hpp"

namespace rtt::model {

Pose Pose::along_z(double z_mm) {
  Pose p;
  p.position[2] = Param(z_mm);
  return p;
}

bool Pose::is_identity() const {
  return *this == Pose{};
}

math::Isometry3 to_isometry(const Pose& pose) {
  const math::Vec3 position(pose.position[0].value, pose.position[1].value, pose.position[2].value);
  const math::Vec3 rotation(pose.rotation_deg[0].value, pose.rotation_deg[1].value,
                            pose.rotation_deg[2].value);
  const math::Vec3 pivot(pose.pivot[0], pose.pivot[1], pose.pivot[2]);
  if (pose.order == PoseOrder::TranslateFirst) {
    // Unchanged since before schema 0.4, so that such poses stay bitwise (ADR 0028, Folgen).
    return math::Isometry3::from_pose(position, rotation, pivot);
  }
  // RotateFirst (ADR 0028, point 3): p_ref = Q + R (t + p - Q) with the pivot Q in the reference
  // frame: rotation R and translation Q + R (t - Q).
  const math::Vec3 angles(math::deg_to_rad(rotation.x()), math::deg_to_rad(rotation.y()),
                          math::deg_to_rad(rotation.z()));
  const math::Mat3 r = math::rotation_xyz_intrinsic(angles);
  return {r, math::Vec3(pivot + r * (position - pivot))};
}

}  // namespace rtt::model
