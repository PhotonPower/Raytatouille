#include <Eigen/LU>
#include <catch2/catch_test_macros.hpp>
#include <cmath>

#include "rtt/math/isometry.hpp"
#include "rtt/math/units.hpp"

using rtt::math::Isometry3;
using rtt::math::Mat3;
using rtt::math::Vec3;

namespace {
constexpr double kTol = 1e-12;

bool near(const Vec3& a, const Vec3& b, double tol = kTol) {
  return (a - b).cwiseAbs().maxCoeff() <= tol;
}
}  // namespace

TEST_CASE("rotations are right-handed", "[isometry]") {
  // +90 deg about X turns +y into +z.
  const Isometry3 rx = Isometry3::from_pose(Vec3::Zero(), Vec3(90, 0, 0), Vec3::Zero());
  REQUIRE(near(rx.apply_vector(Vec3::UnitY()), Vec3::UnitZ()));
  // +90 deg about Y turns +z into +x.
  const Isometry3 ry = Isometry3::from_pose(Vec3::Zero(), Vec3(0, 90, 0), Vec3::Zero());
  REQUIRE(near(ry.apply_vector(Vec3::UnitZ()), Vec3::UnitX()));
  // +90 deg about Z turns +x into +y.
  const Isometry3 rz = Isometry3::from_pose(Vec3::Zero(), Vec3(0, 0, 90), Vec3::Zero());
  REQUIRE(near(rz.apply_vector(Vec3::UnitX()), Vec3::UnitY()));
}

TEST_CASE("rotation order is intrinsic X then Y then Z", "[isometry]") {
  const Vec3 a(10.0, 20.0, 30.0);
  const Mat3 r = rtt::math::rotation_xyz_intrinsic(Vec3(
      rtt::math::deg_to_rad(a.x()), rtt::math::deg_to_rad(a.y()), rtt::math::deg_to_rad(a.z())));
  const Isometry3 x = Isometry3::from_pose(Vec3::Zero(), Vec3(a.x(), 0, 0), Vec3::Zero());
  const Isometry3 y = Isometry3::from_pose(Vec3::Zero(), Vec3(0, a.y(), 0), Vec3::Zero());
  const Isometry3 z = Isometry3::from_pose(Vec3::Zero(), Vec3(0, 0, a.z()), Vec3::Zero());
  REQUIRE((x * y * z).is_approx(Isometry3(r, Vec3::Zero())));
  // Rotation matrix stays orthonormal with determinant +1.
  REQUIRE((r * r.transpose() - Mat3::Identity()).cwiseAbs().maxCoeff() <= kTol);
  REQUIRE(std::abs(r.determinant() - 1.0) <= kTol);
}

TEST_CASE("pose: translation, then rotation about the pivot", "[isometry]") {
  const Vec3 position(1.0, 2.0, 3.0);
  const Vec3 pivot(0.0, 0.0, 5.0);
  const Isometry3 t = Isometry3::from_pose(position, Vec3(0, 0, 90), pivot);
  // The pivot point does not move by the rotation, only by the translation.
  REQUIRE(near(t.apply_point(pivot), position + pivot));
  // The child origin rotates about the pivot (rotation about z leaves an on-axis point fixed).
  REQUIRE(near(t.apply_point(Vec3::Zero()), position));
  // Off-axis point (1,0,5) rotates to (0,1,5) relative to the pivot.
  REQUIRE(near(t.apply_point(Vec3(1.0, 0.0, 5.0)), position + Vec3(0.0, 1.0, 5.0)));
}

TEST_CASE("inverse and composition", "[isometry]") {
  const Isometry3 a = Isometry3::from_pose(Vec3(1, -2, 3), Vec3(12, -7, 33), Vec3(0.5, 0, 2));
  const Isometry3 b = Isometry3::from_pose(Vec3(-4, 0, 9), Vec3(-45, 3, 0), Vec3(0, 1, 0));
  const Vec3 p(0.3, -1.7, 42.0);

  REQUIRE((a * a.inverse()).is_approx(Isometry3::identity()));
  REQUIRE((a.inverse() * a).is_approx(Isometry3::identity()));
  REQUIRE(near((a * b).apply_point(p), a.apply_point(b.apply_point(p))));
  REQUIRE(near(a.inverse().apply_point(a.apply_point(p)), p));
  // Directions are not translated.
  REQUIRE(near(a.apply_vector(Vec3::Zero()), Vec3::Zero()));
}
