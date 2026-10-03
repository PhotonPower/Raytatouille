#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "rtt/math/real.hpp"
#include "rtt/math/units.hpp"

using rtt::math::deg_to_rad;
using rtt::math::kPi;
using rtt::math::rad_to_deg;

TEST_CASE("degree/radian conversion round-trips", "[units]") {
  REQUIRE(deg_to_rad(180.0) == Catch::Approx(kPi).epsilon(1e-15));
  REQUIRE(rad_to_deg(deg_to_rad(37.5)) == Catch::Approx(37.5).epsilon(1e-15));
}

TEST_CASE("wavelength unit conversion", "[units]") {
  REQUIRE(rtt::math::um_to_mm(0.5876) == Catch::Approx(0.5876e-3).epsilon(1e-15));
  REQUIRE(rtt::math::mm_to_um(rtt::math::um_to_mm(1.55)) == Catch::Approx(1.55).epsilon(1e-15));
}

TEST_CASE("Real concept accepts floating point types only", "[units]") {
  STATIC_REQUIRE(rtt::math::Real<double>);
  STATIC_REQUIRE(rtt::math::Real<float>);
  STATIC_REQUIRE_FALSE(rtt::math::Real<int>);
}
