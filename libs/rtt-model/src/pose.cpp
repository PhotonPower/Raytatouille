#include "rtt/model/pose.hpp"

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
  return math::Isometry3::from_pose(position, rotation, pivot);
}

}  // namespace rtt::model
