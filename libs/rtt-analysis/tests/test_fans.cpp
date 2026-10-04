#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#include "rtt/analysis/spot.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
#include "rtt/paraxial/paraxial.hpp"

using rtt::compile::compile;
using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::material::MaterialLibrary;
using rtt::model::Element;
using rtt::model::Pose;
using rtt::model::System;
using rtt::trace::RayStatus;

// Ray fans: transverse aberration of the tangential (px = 0) and sagittal (py = 0) fans
// relative to the chief ray on the image surface (decided for #28). Reference cases of #28.

namespace {

System load(const std::string& relative) {
  return rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/" + relative);
}

Element& element(System& s, std::size_t i) {
  return std::get<Element>(s.root.children[i].value);
}

}  // namespace

TEST_CASE("ideal system: ray fans vanish", "[fans]") {
  // Paraboloid on axis is free of aberrations: issue #28, < 1e-12 mm.
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(load("m2/paraboloid_stop.rtt.json"), lib);
  const auto fan = rtt::analysis::ray_fan(cs, PathId{0}, 0, 0);
  REQUIRE(fan.tangential.size() == 21);
  REQUIRE(fan.sagittal.size() == 21);
  for (const auto* points : {&fan.tangential, &fan.sagittal}) {
    for (const auto& p : *points) {
      REQUIRE(p.status == RayStatus::Alive);
      REQUIRE(std::abs(p.ex) < 1e-12);
      REQUIRE(std::abs(p.ey) < 1e-12);
    }
  }
  REQUIRE(fan.tangential.front().p == -1.0);
  REQUIRE(fan.tangential.back().p == 1.0);
}

TEST_CASE("plano-convex lens at small NA: tangential aberration scales with py^3", "[fans]") {
  // Third-order transverse spherical aberration is proportional to the cube of the pupil
  // coordinate, provided the image surface lies in the paraxial focus (otherwise a defocus term
  // proportional to py dominates). Issue #28: fit exponent 3 +- 0.05.
  // The independent reference value below assumes a surrounding index of 1, so the file's AIR
  // (Ciddor since #25) is replaced by VACUUM.
  System s = load("m1/singlet_const.rtt.json");
  s.environment.medium = "VACUUM";
  const MaterialLibrary lib;
  const CompiledSystem first = compile(s, lib);
  const double z_focus = *rtt::paraxial::first_order(first, PathId{0}, 1).rear_focal_z;
  element(s, 2).pose = Pose::along_z(z_focus);
  const CompiledSystem cs = compile(s, lib);

  rtt::analysis::FanOptions options;
  options.points = 201;  // step 0.01 in py
  const auto fan = rtt::analysis::ray_fan(cs, PathId{0}, 0, 1, options);
  const std::vector<std::size_t> index{102, 104, 108, 116};  // py = 0.02, 0.04, 0.08, 0.16
  std::vector<double> lp;
  std::vector<double> le;
  for (const std::size_t i : index) {
    const auto& p = fan.tangential[i];
    REQUIRE(p.status == RayStatus::Alive);
    // On axis the fan is antisymmetric: epsilon_y(-py) = -epsilon_y(py).
    const auto& m = fan.tangential[200 - i];
    REQUIRE(std::abs(m.ey + p.ey) <= 1e-12);
    lp.push_back(std::log(p.p));
    le.push_back(std::log(std::abs(p.ey)));
  }
  REQUIRE(std::abs(fan.tangential[102].p - 0.02) <= 1e-15);
  // Sign and size: the plano-convex lens is undercorrected, the marginal-zone rays cross the
  // axis before the paraxial focus, so epsilon_y < 0 for py > 0. Independent check (single-ray
  // trace through the sphere R = 51.68, n = 1.5168, plane at z = 9, image at the paraxial
  // focus, reviews of #28): epsilon_y(0.16) = -4.46348e-4 mm.
  REQUIRE(fan.tangential[116].ey < 0.0);
  REQUIRE(std::abs(fan.tangential[116].ey + 4.46348e-4) <= 2e-5 * 4.46348e-4);
  double mx = 0.0, my = 0.0;
  for (std::size_t i = 0; i < lp.size(); ++i) {
    mx += lp[i];
    my += le[i];
  }
  mx /= static_cast<double>(lp.size());
  my /= static_cast<double>(lp.size());
  double sxy = 0.0, sxx = 0.0;
  for (std::size_t i = 0; i < lp.size(); ++i) {
    sxy += (lp[i] - mx) * (le[i] - my);
    sxx += (lp[i] - mx) * (lp[i] - mx);
  }
  const double exponent = sxy / sxx;
  INFO("fit exponent " << exponent);
  REQUIRE(std::abs(exponent - 3.0) <= 0.05);
}

TEST_CASE("fans are relative to the chief ray: zero at the pupil centre", "[fans]") {
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(load("m1/singlet_const.rtt.json"), lib);
  const auto fan = rtt::analysis::ray_fan(cs, PathId{0}, 2, 1);
  // The centre point of each fan is the chief ray itself (reference wavelength).
  REQUIRE(fan.tangential[10].p == 0.0);
  REQUIRE(fan.tangential[10].ex == 0.0);
  REQUIRE(fan.tangential[10].ey == 0.0);
  REQUIRE(fan.sagittal[10].ex == 0.0);
  // Off axis the chief ray lands away from the axis (field 5 deg).
  REQUIRE(fan.chief.y > 1.0);
}

TEST_CASE("vignetted fan points keep their status and zero aberration", "[fans]") {
  System s = load("m2/paraboloid_stop.rtt.json");
  element(s, 1).surfaces[0].aperture = rtt::model::CircularAperture{21.0, 0.0};
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(s, lib);
  const auto fan = rtt::analysis::ray_fan(cs, PathId{0}, 0, 0);
  // py = -1 (radius 30 mm) is outside the mirror aperture of 21 mm.
  REQUIRE(fan.tangential.front().status == RayStatus::Vignetted);
  REQUIRE(fan.tangential.front().ex == 0.0);
  REQUIRE(fan.tangential.front().ey == 0.0);
  REQUIRE(fan.tangential[10].status == RayStatus::Alive);
}

TEST_CASE("invalid fan input throws std::invalid_argument", "[fans]") {
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(load("m1/singlet_const.rtt.json"), lib);
  rtt::analysis::FanOptions none;
  none.points = 0;
  REQUIRE_THROWS_AS(rtt::analysis::ray_fan(cs, PathId{0}, 0, 1, none), std::invalid_argument);
  REQUIRE_THROWS_AS(rtt::analysis::ray_fan(cs, PathId{0}, 0, 5), std::invalid_argument);
  REQUIRE_THROWS_AS(rtt::analysis::ray_fan(cs, PathId{0}, 9, 1), std::invalid_argument);
}
