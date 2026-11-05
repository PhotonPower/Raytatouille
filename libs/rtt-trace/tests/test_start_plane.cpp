#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <limits>
#include <optional>

#include "rtt/math/isometry.hpp"
#include "rtt/math/types.hpp"
#include "start_plane.hpp"

// ADR 0041, point 6: the start plane of the aiming begins at the vertex (0, 0, sag(0, 0)), not
// at the origin of the surface. The test shape exists only here, not in the model.

using rtt::math::Isometry3;
using rtt::math::Mat3;
using rtt::math::Vec3;
using rtt::trace::detail::lowest_z;
using rtt::trace::detail::vertex_global_z;

namespace {

/// z = value + a r^2 with the vertex shifted by `value` along the local z axis.
struct ShiftedBowl {
  double value = 0.0;
  double a = 0.0;
  [[nodiscard]] double sag(double x, double y) const noexcept {
    return value + a * (x * x + y * y);
  }
};

constexpr double kInf = std::numeric_limits<double>::infinity();

}  // namespace

TEST_CASE("the start plane begins at the shifted vertex", "[sources][vertex]") {
  // Local z along global +z, origin at z = 5, vertex at 4.7; the surface rises from the vertex,
  // so its lowest point is the vertex itself, which the sampling from r = rho/64 on misses.
  const Isometry3 along(Mat3::Identity(), Vec3(0.0, 0.0, 5.0));
  const ShiftedBowl bowl{-0.3, 0.01};
  const std::optional<double> low = lowest_z(bowl, along, 10.0, 0.0, 5.0, 0.0);
  REQUIRE(low.has_value());
  CHECK(*low == 5.0 + -0.3);
  // Content guard: every sampled point lies above the vertex.
  CHECK(*low < 5.0 + bowl.sag(0.0, 5.0 / 64.0));

  // Local z along global -z (rotated by pi about x): the vertex lies at 5.3, and z = value - a r^2
  // rises in global z away from it.
  const Isometry3 against(Vec3(1.0, -1.0, -1.0).asDiagonal().toDenseMatrix(), Vec3(0.0, 0.0, 5.0));
  const ShiftedBowl flipped{-0.3, -0.01};
  const std::optional<double> low_flipped = lowest_z(flipped, against, 10.0, 0.0, 5.0, 0.0);
  REQUIRE(low_flipped.has_value());
  CHECK(*low_flipped == 5.0 - -0.3);
}

TEST_CASE("without a constant part the vertex is the origin bit for bit", "[sources][vertex]") {
  for (const double z0 : {7.25, 0.0, -0.0}) {
    const Isometry3 along(Mat3::Identity(), Vec3(0.0, 0.0, z0));
    for (const double sag0 : {0.0, -0.0}) {
      const double z = vertex_global_z(along, sag0);
      CHECK(z == z0);
      CHECK(std::signbit(z) == std::signbit(z0));
    }
  }
  const Isometry3 along(Mat3::Identity(), Vec3(0.0, 0.0, 2.0));
  CHECK(vertex_global_z(along, 0.5) == 2.5);
  CHECK(lowest_z(ShiftedBowl{0.0, 0.01}, along, kInf, 0.0, 5.0, 0.0) == 2.0);
}
