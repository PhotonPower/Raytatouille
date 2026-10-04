#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#include "rtt/analysis/opd.hpp"
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

// OPD convention and defocus formulas: J. C. Wyant, K. Creath, "Basic Wavefront Aberration
// Theory for Optical Metrology", Applied Optics and Optical Engineering XI (1992):
// Sec. I: W > 0 if the wavefront leads the reference (curves in more).
// Eq. (17): a shift of the image plane by eps_z along the propagation direction needs
// Delta W = -eps_z (x^2 + y^2) / (2 R^2) to move the focus there; Eq. (18): at the pupil edge
// Delta W = -1/2 eps_z sin^2 U. A wavefront that still focuses on the old point therefore has
// W = +1/2 eps_z sin^2 U against the reference sphere about the shifted point.

namespace {

System load(const std::string& relative) {
  return rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/" + relative);
}

Element& element(System& s, std::size_t i) {
  return std::get<Element>(s.root.children[i].value);
}

constexpr double kLambdaMm = 0.5876e-3;  // reference wavelength of the test systems, mm

/// Least-squares slope of log|w| over log p.
double fit_exponent(const std::vector<double>& p, const std::vector<double>& w) {
  double mx = 0.0, my = 0.0;
  for (std::size_t i = 0; i < p.size(); ++i) {
    mx += std::log(p[i]);
    my += std::log(std::abs(w[i]));
  }
  mx /= static_cast<double>(p.size());
  my /= static_cast<double>(p.size());
  double sxy = 0.0, sxx = 0.0;
  for (std::size_t i = 0; i < p.size(); ++i) {
    sxy += (std::log(p[i]) - mx) * (std::log(std::abs(w[i])) - my);
    sxx += (std::log(p[i]) - mx) * (std::log(p[i]) - mx);
  }
  return sxy / sxx;
}

}  // namespace

TEST_CASE("paraboloid on axis: OPD below 1e-6 waves", "[opd]") {
  // Architecture, Validierung: paraboloid, object at infinity, on axis, OPD < 1e-6 lambda.
  // Mirror case: the image space runs towards -z.
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(load("m2/paraboloid_stop.rtt.json"), lib);
  const auto map = rtt::analysis::opd_map(cs, PathId{0}, 0, 0);
  REQUIRE(map.arrived == map.points.size());
  REQUIRE(map.vignetted == 0);
  for (const auto& p : map.points) {
    REQUIRE(p.status == RayStatus::Alive);
    REQUIRE(std::abs(p.w) < 1e-6);
  }
  REQUIRE(map.rms < 1e-6);
  REQUIRE(map.pv < 1e-6);
  // Exit pupil: the stop 50 mm before the mirror (f = 100 mm) is imaged to z = +100, the image
  // lies at z = -100, so R = 200 mm.
  REQUIRE(std::abs(map.sphere.radius - 200.0) < 1e-9);
  REQUIRE(std::abs(map.sphere.centre.z() + 100.0) < 1e-9);
}

TEST_CASE("defocus of an ideal image: W = eps_z sin^2(U) / 2 with both signs", "[opd]") {
  // Paraboloid on axis with the image surface moved by dz along global z. The image space runs
  // towards -z, so the shift along the propagation direction is eps_z = -dz. The exit pupil lies
  // at z = +100: dz > 0 moves the detector towards the pupil (eps_z < 0, W < 0).
  // The source is first order in eps_z and uses sin^2 U instead of 2 (1 - cos U): relative
  // deviation O(U^2) ~ 0.3 % at p <= 0.3 (U < 5 deg) and O(eps_z / R) ~ 5e-5, hence 1e-2.
  for (const double dz : {0.01, -0.01}) {
    System s = load("m2/paraboloid_stop.rtt.json");
    element(s, 2).pose = Pose::along_z(-100.0 + dz);
    const MaterialLibrary lib;
    const CompiledSystem cs = compile(s, lib);
    rtt::analysis::OpdOptions options;
    options.fan_points = 21;  // p = -1, -0.9, ..., 1
    const auto fan = rtt::analysis::opd_fan(cs, PathId{0}, 0, 0, options);
    const double eps_z = -dz;
    int checked = 0;
    for (const auto& p : fan.tangential) {
      if (std::abs(p.py) > 0.3 + 1e-12 || p.py == 0.0) continue;
      REQUIRE(p.status == RayStatus::Alive);
      // sin U of the ray in image space: height on the 30 mm pupil over the distance to the
      // focus of the paraboloid (z_mirror = -h^2 / 400, focus at z = -100).
      const double h = 30.0 * std::abs(p.py);
      const double z_m = -h * h / 400.0;
      const double sin_u = h / std::hypot(h, z_m + 100.0);
      const double expected = 0.5 * eps_z * sin_u * sin_u / kLambdaMm;
      INFO("dz " << dz << ", py " << p.py << ", W " << p.w << ", expected " << expected);
      REQUIRE((p.w > 0.0) == (eps_z > 0.0));
      REQUIRE(std::abs(p.w - expected) <= 1e-2 * std::abs(expected));
      ++checked;
    }
    REQUIRE(checked == 6);  // py = +-0.1, +-0.2, +-0.3
  }
}

TEST_CASE("plano-convex lens at small NA: OPD scales with p^4, positive at the rim", "[opd]") {
  // Third-order spherical aberration W040 p^4 with the reference sphere in the paraxial focus;
  // the undercorrected lens focuses the rim in front of the paraxial focus, i.e. the rim of
  // the wavefront curves in more: W > 0 (Wyant & Creath, Sec. I). Exponent 4 +- 0.05.
  System s = load("m1/singlet_const.rtt.json");
  const MaterialLibrary lib;
  const CompiledSystem first = compile(s, lib);
  element(s, 2).pose = Pose::along_z(*rtt::paraxial::first_order(first, PathId{0}, 1).rear_focal_z);
  const CompiledSystem cs = compile(s, lib);
  rtt::analysis::OpdOptions options;
  options.fan_points = 81;  // step 0.025
  const auto fan = rtt::analysis::opd_fan(cs, PathId{0}, 0, 1, options);
  std::vector<double> p;
  std::vector<double> w;
  for (const std::size_t i : {42, 44, 48, 56}) {  // p = 0.05, 0.1, 0.2, 0.4
    REQUIRE(fan.tangential[i].status == RayStatus::Alive);
    p.push_back(fan.tangential[i].py);
    w.push_back(fan.tangential[i].w);
    // Rotational symmetry on axis: W(-p) = W(p), sagittal equals tangential.
    REQUIRE(std::abs(fan.tangential[80 - i].w - fan.tangential[i].w) <= 1e-9);
    REQUIRE(std::abs(fan.sagittal[i].w - fan.tangential[i].w) <= 1e-9);
  }
  REQUIRE(std::abs(p[0] - 0.05) <= 1e-15);
  const double exponent = fit_exponent(p, w);
  INFO("fit exponent " << exponent);
  REQUIRE(std::abs(exponent - 4.0) <= 0.05);
  REQUIRE(fan.tangential.back().w > 0.0);
}

TEST_CASE("OPD map: chief ray zero, RMS as standard deviation, PV", "[opd]") {
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(load("m1/singlet_const.rtt.json"), lib);
  for (const std::uint16_t wl : {0, 1, 2}) {
    rtt::analysis::OpdOptions options;
    options.grid = 9;
    const auto map = rtt::analysis::opd_map(cs, PathId{0}, 2, wl, options);
    double sum = 0.0;
    double sum2 = 0.0;
    double lo = 1e300;
    double hi = -1e300;
    std::size_t n = 0;
    for (const auto& q : map.points) {
      if (q.px == 0.0 && q.py == 0.0) REQUIRE(q.w == 0.0);  // chief ray of the same wavelength
      if (q.status != RayStatus::Alive) {
        REQUIRE(q.w == 0.0);
        continue;
      }
      sum += q.w;
      sum2 += q.w * q.w;
      lo = std::min(lo, q.w);
      hi = std::max(hi, q.w);
      ++n;
    }
    REQUIRE(n == map.arrived);
    REQUIRE(map.arrived + map.vignetted == map.points.size());
    const double mean = sum / static_cast<double>(n);
    REQUIRE(std::abs(map.rms - std::sqrt(sum2 / static_cast<double>(n) - mean * mean)) <= 1e-9);
    REQUIRE(map.pv == hi - lo);
  }
}

TEST_CASE("vignetted OPD points are marked and excluded", "[opd]") {
  System s = load("m2/paraboloid_stop.rtt.json");
  element(s, 1).surfaces[0].aperture = rtt::model::CircularAperture{21.0, 0.0};
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(s, lib);
  const auto map = rtt::analysis::opd_map(cs, PathId{0}, 0, 0);
  REQUIRE(map.vignetted > 0);
  for (const auto& q : map.points) {
    if (q.px * q.px + q.py * q.py > 0.75 * 0.75) {  // radius > 22.5 mm on the 30 mm pupil
      REQUIRE(q.status == RayStatus::Vignetted);
      REQUIRE(q.w == 0.0);
    }
  }
  REQUIRE(map.pv < 1e-6);
}

TEST_CASE("invalid OPD input", "[opd]") {
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(load("m1/singlet_const.rtt.json"), lib);
  rtt::analysis::OpdOptions bad;
  bad.grid = 0;
  REQUIRE_THROWS_AS(rtt::analysis::opd_map(cs, PathId{0}, 0, 1, bad), std::invalid_argument);
  bad = {};
  bad.fan_points = 0;
  REQUIRE_THROWS_AS(rtt::analysis::opd_fan(cs, PathId{0}, 0, 1, bad), std::invalid_argument);
  REQUIRE_THROWS_AS(rtt::analysis::opd_map(cs, PathId{0}, 0, 7), std::invalid_argument);
  REQUIRE_THROWS_AS(rtt::analysis::opd_map(cs, PathId{0}, 9, 1), std::invalid_argument);
}
