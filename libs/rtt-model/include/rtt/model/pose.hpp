#pragma once

/// @file pose.hpp
/// Placement of a node (assembly, element, surface) in its parent coordinate system.

#include <array>

#include "rtt/math/isometry.hpp"
#include "rtt/model/param.hpp"

namespace rtt::model {

/// Translation, then intrinsic X->Y->Z rotation about the pivot (see rtt/math/isometry.hpp).
struct Pose {
  std::array<Param, 3> position{};      ///< mm, in parent coordinates
  std::array<Param, 3> rotation_deg{};  ///< degree, intrinsic X -> Y -> Z
  std::array<double, 3> pivot{};        ///< mm, in local coordinates; tolerances tilt about it

  /// Convenience constructor for the common case of a pure shift along the local z axis.
  [[nodiscard]] static Pose along_z(double z_mm);

  [[nodiscard]] bool is_identity() const;

  bool operator==(const Pose&) const = default;
};

/// Evaluates the current parameter values into a rigid transform child -> parent.
[[nodiscard]] math::Isometry3 to_isometry(const Pose& pose);

}  // namespace rtt::model
