#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>

#include "../src/thickness.hpp"
#include "rtt/math/isometry.hpp"
#include "rtt/math/types.hpp"

// ADR 0041, point 6: the centre thickness of the dimension report is the distance of the
// vertices (0, 0, sag(0, 0)), not of the origins of the two surfaces.

using Catch::Matchers::WithinAbs;
using rtt::analysis::detail::centre_thickness;
using rtt::math::Isometry3;
using rtt::math::Mat3;
using rtt::math::Vec3;

TEST_CASE("the centre thickness follows shifted vertices", "[reports][vertex]") {
  const Isometry3 first;  // at the global origin, local z along +z
  const Isometry3 second(Mat3::Identity(), Vec3(0.0, 0.0, 8.0));
  // Vertices at z = 0.2 and z = 8 - 0.1.
  CHECK_THAT(centre_thickness(first.inverse(), 0.2, second, -0.1), WithinAbs(7.7, 1e-14));
  // Second surface with its local z along -z (a flipped element): its vertex lies at 8 + 0.1.
  const Isometry3 flipped(Vec3(1.0, -1.0, -1.0).asDiagonal().toDenseMatrix(), Vec3(0.0, 0.0, 8.0));
  CHECK_THAT(centre_thickness(first.inverse(), 0.2, flipped, -0.1), WithinAbs(7.9, 1e-14));
  // Content guard: the shift changes the result.
  CHECK(centre_thickness(first.inverse(), 0.2, second, -0.1) !=
        centre_thickness(first.inverse(), 0.0, second, 0.0));
}

TEST_CASE("without a constant part the centre thickness is the distance of the origins",
          "[reports][vertex]") {
  const Isometry3 first(Mat3::Identity(), Vec3(0.0, 0.0, 3.0));
  for (const double z2 : {11.5, 3.0}) {
    const Isometry3 second(Mat3::Identity(), Vec3(0.0, 0.0, z2));
    const double origins = first.inverse().apply_point(second.apply_point(Vec3::Zero())).z();
    for (const double sag0 : {0.0, -0.0}) {
      const double t = centre_thickness(first.inverse(), sag0, second, sag0);
      CHECK(t == origins);
      CHECK(std::signbit(t) == std::signbit(origins));
    }
  }
}
