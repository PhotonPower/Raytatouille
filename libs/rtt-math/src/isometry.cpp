#include "rtt/math/isometry.hpp"

#include <cmath>

#include "rtt/math/units.hpp"

namespace rtt::math {

namespace {

Mat3 rot_x(double a) {
  const double c = std::cos(a);
  const double s = std::sin(a);
  Mat3 m;
  m << 1.0, 0.0, 0.0,  //
      0.0, c, -s,      //
      0.0, s, c;
  return m;
}

Mat3 rot_y(double a) {
  const double c = std::cos(a);
  const double s = std::sin(a);
  Mat3 m;
  m << c, 0.0, s,     //
      0.0, 1.0, 0.0,  //
      -s, 0.0, c;
  return m;
}

Mat3 rot_z(double a) {
  const double c = std::cos(a);
  const double s = std::sin(a);
  Mat3 m;
  m << c, -s, 0.0,  //
      s, c, 0.0,    //
      0.0, 0.0, 1.0;
  return m;
}

}  // namespace

Mat3 rotation_xyz_intrinsic(const Vec3& angles_rad) {
  return rot_x(angles_rad.x()) * rot_y(angles_rad.y()) * rot_z(angles_rad.z());
}

Isometry3::Isometry3(const Mat3& rotation, const Vec3& translation)
    : rotation_(rotation), translation_(translation) {}

Isometry3 Isometry3::from_pose(const Vec3& position, const Vec3& rotation_deg, const Vec3& pivot) {
  const Vec3 angles(deg_to_rad(rotation_deg.x()), deg_to_rad(rotation_deg.y()),
                    deg_to_rad(rotation_deg.z()));
  const Mat3 r = rotation_xyz_intrinsic(angles);
  return {r, position + pivot - r * pivot};
}

Isometry3 Isometry3::inverse() const {
  const Mat3 rt = rotation_.transpose();
  return {rt, -(rt * translation_)};
}

bool Isometry3::is_approx(const Isometry3& other, double tol) const {
  return (rotation_ - other.rotation_).cwiseAbs().maxCoeff() <= tol &&
         (translation_ - other.translation_).cwiseAbs().maxCoeff() <= tol;
}

Isometry3 operator*(const Isometry3& a, const Isometry3& b) {
  return {a.rotation_ * b.rotation_, a.rotation_ * b.translation_ + a.translation_};
}

}  // namespace rtt::math
