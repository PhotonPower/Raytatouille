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
