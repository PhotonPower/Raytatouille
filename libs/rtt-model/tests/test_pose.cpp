#include <catch2/catch_test_macros.hpp>

#include "rtt/model/pose.hpp"

using rtt::math::Vec3;
using rtt::model::Param;
using rtt::model::Pose;

TEST_CASE("default pose is the identity", "[pose]") {
  const Pose p;
  REQUIRE(p.is_identity());
  REQUIRE(rtt::model::to_isometry(p).is_approx(rtt::math::Isometry3::identity()));
}

TEST_CASE("along_z shifts the origin", "[pose]") {
  const Pose p = Pose::along_z(4.0);
  REQUIRE_FALSE(p.is_identity());
  const Vec3 o = rtt::model::to_isometry(p).apply_point(Vec3::Zero());
  REQUIRE(o.isApprox(Vec3(0.0, 0.0, 4.0)));
}

TEST_CASE("variable flags do not change the transform", "[pose]") {
  Pose a = Pose::along_z(4.0);
  Pose b = a;
  b.position[2].variable = true;
  REQUIRE_FALSE(a == b);
  REQUIRE(rtt::model::to_isometry(a).is_approx(rtt::model::to_isometry(b)));
}

TEST_CASE("rotation about the pivot", "[pose]") {
  Pose p;
  p.rotation_deg[0] = Param(90.0);
  p.pivot = {0.0, 0.0, 10.0};
  const auto t = rtt::model::to_isometry(p);
  REQUIRE(t.apply_point(Vec3(0.0, 0.0, 10.0)).isApprox(Vec3(0.0, 0.0, 10.0)));
  // Origin sits 10 mm below the pivot along -z; +90 deg about x moves it to +y.
  REQUIRE((t.apply_point(Vec3::Zero()) - Vec3(0.0, 10.0, 10.0)).norm() < 1e-12);
}

// ---------------------------------------------------------- ADR 0028, Pose.order -----

TEST_CASE("translate_first is the default and keeps the transform of from_pose", "[pose]") {
  Pose p;
  REQUIRE(p.order == rtt::model::PoseOrder::TranslateFirst);
  REQUIRE(p.reference == rtt::model::PoseReference::Absolute);
  p.position = {Param(1.0), Param(-2.0), Param(3.0)};
  p.rotation_deg = {Param(10.0), Param(20.0), Param(30.0)};
  p.pivot = {0.5, 0.25, -4.0};
  const auto expected = rtt::math::Isometry3::from_pose(
      Vec3(1.0, -2.0, 3.0), Vec3(10.0, 20.0, 30.0), Vec3(0.5, 0.25, -4.0));
  const auto t = rtt::model::to_isometry(p);
  // Bitwise: translate_first calls Isometry3::from_pose unchanged (ADR 0028, Folgen).
  for (const Vec3& q : {Vec3(0.0, 0.0, 0.0), Vec3(1.0, 2.0, 3.0), Vec3(-7.0, 0.5, 11.0)}) {
    REQUIRE(t.apply_point(q) == expected.apply_point(q));
  }
}

TEST_CASE("rotate_first rotates the reference frame about the pivot, then shifts in it", "[pose]") {
  // ADR 0028, point 3: p_ref = Q + R (t + p - Q), Q in the reference frame.
  Pose p;
  p.order = rtt::model::PoseOrder::RotateFirst;
  p.rotation_deg[0] = Param(90.0);
  p.position = {Param(0.0), Param(0.0), Param(5.0)};
  SECTION("pivot 0: t runs along the rotated axes") {
    // Rx(90) z = (0, -1, 0): 5 mm along it.
    const auto t = rtt::model::to_isometry(p);
    REQUIRE((t.apply_point(Vec3::Zero()) - Vec3(0.0, -5.0, 0.0)).norm() < 1e-12);
    REQUIRE((t.apply_vector(Vec3::UnitZ()) - Vec3(0.0, -1.0, 0.0)).norm() < 1e-12);
  }
  SECTION("pivot in the reference frame") {
    // Q = (0, 0, 10): Q + R (t - Q) = (0, 0, 10) + Rx(90) (0, 0, -5) = (0, 5, 10).
    p.pivot = {0.0, 0.0, 10.0};
    const auto t = rtt::model::to_isometry(p);
    REQUIRE((t.apply_point(Vec3::Zero()) - Vec3(0.0, 5.0, 10.0)).norm() < 1e-12);
  }
  SECTION("without rotation both orders agree") {
    p.rotation_deg[0] = Param(0.0);
    p.pivot = {1.0, 2.0, 3.0};
    Pose q = p;
    q.order = rtt::model::PoseOrder::TranslateFirst;
    REQUIRE(rtt::model::to_isometry(p).is_approx(rtt::model::to_isometry(q)));
  }
}

TEST_CASE("rotate_first puts a lens on the axis behind a 45 degree fold mirror (ADR 0028)",
          "[pose]") {
  // Mirror at z = 50, rotation (45, 0, 0); the lens relative to it with rotate_first,
  // rotation (-135, 0, 0) and position (0, 0, 30): origin (0, 30, 50), z axis +y.
  Pose mirror = Pose::along_z(50.0);
  mirror.rotation_deg[0] = Param(45.0);
  Pose lens;
  lens.order = rtt::model::PoseOrder::RotateFirst;
  lens.rotation_deg[0] = Param(-135.0);
  lens.position[2] = Param(30.0);
  const auto global = rtt::model::to_isometry(mirror) * rtt::model::to_isometry(lens);
  REQUIRE((global.apply_point(Vec3::Zero()) - Vec3(0.0, 30.0, 50.0)).norm() < 1e-12);
  REQUIRE((global.apply_vector(Vec3::UnitZ()) - Vec3(0.0, 1.0, 0.0)).norm() < 1e-12);
}
