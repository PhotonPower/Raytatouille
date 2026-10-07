#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstdint>
#include <optional>
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
  // deviation tan^2(U/2) ~ 0.2 % at p <= 0.3 (U ~ 5.2 deg) and O(eps_z / R) ~ 5e-5, hence 1e-2.
  for (const double dz : {0.01, -0.01}) {
    System s = load("m2/paraboloid_stop.rtt.json");
    s.environment.medium = "VACUUM";  // the source formula assumes n' = 1
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
  for (const std::size_t i : std::vector<std::size_t>{42, 44, 48, 56}) {  // p = 0.05, 0.1, 0.2, 0.4
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
  for (const std::uint16_t wl : std::vector<std::uint16_t>{0, 1, 2}) {
    rtt::analysis::OpdOptions options;
    options.grid = 9;
    const auto map = rtt::analysis::opd_map(cs, PathId{0}, 2, wl, options);
    double sum = 0.0;
    double sum2 = 0.0;
    double lo = 1e300;
    double hi = -1e300;
    std::size_t n = 0;
    int centre = 0;
    for (const auto& q : map.points) {
      if (q.px == 0.0 && q.py == 0.0) {
        REQUIRE(q.w == 0.0);  // chief ray of the same wavelength
        ++centre;
      }
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
    REQUIRE(centre == 1);
  }
}

TEST_CASE("OPD of every wavelength is in waves at the reference wavelength", "[opd]") {
  // Vacuum outside and a constant glass index: all wavelengths trace identical rays, so their
  // OPD in mm is identical. Normalised to the reference wavelength (decided for #29) the waves
  // must agree too; normalising to each ray's own wavelength would scale them by
  // 0.5876 / 0.4861.
  System s = load("m1/singlet_const.rtt.json");
  s.environment.medium = "VACUUM";
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(s, lib);
  REQUIRE(cs.reference_wavelength() == 1);
  rtt::analysis::OpdOptions options;
  options.grid = 9;
  const auto ref = rtt::analysis::opd_map(cs, PathId{0}, 2, 1, options);
  const auto blue = rtt::analysis::opd_map(cs, PathId{0}, 2, 0, options);
  REQUIRE(ref.points.size() == blue.points.size());
  double largest = 0.0;
  for (std::size_t k = 0; k < ref.points.size(); ++k) {
    REQUIRE(blue.points[k].w == ref.points[k].w);
    largest = std::max(largest, std::abs(ref.points[k].w));
  }
  REQUIRE(largest > 0.1);  // the field is aberrated enough to see a wrong normalisation
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

namespace {

// Image-space telecentric systems (#102). Sources: the OPD definition of Wyant & Creath (Sec. I,
// see the top of this file) with the reference sphere about the chief-ray image point C and the
// radius R = |C - XP| -> infinity: the OPL up to the sphere is OPL - n' s with
// s = +-R + b +- (b^2 - |p - C|^2) / (2 R) + ..., b = d . (p - C); the common +-n' R cancels in
// W = OPL_chief - OPL_ray, so W_inf = OPL_chief - OPL_ray + n' d . (p - C) (derivation in
// libs/rtt-analysis/src/opd.cpp and docs/quellen.md). Descartes lens: J. Sasian, OPTI 518
// lecture notes L14 "Aspheric surfaces and stop shifting", p. 5: a plano-hyperbolic lens with
// K = -n^2 focuses a collimated beam entering the plane side without aberration; the reverse
// path collimates a point at its focus. Focal distance from the hyperbola vertex |R| / (n - 1)
// (single refracting surface, y-nu equations, Greivenkamp OPTI-201/202, Sec. 9, p. 9-2).

/// Largest |W| over the arrived points of a map, waves.
double max_abs_w(const rtt::analysis::OpdMap& map) {
  double m = 0.0;
  for (const auto& q : map.points) {
    if (q.status == RayStatus::Alive) m = std::max(m, std::abs(q.w));
  }
  return m;
}

/// Largest |W_a(px, py) - W_b(label px, label py)| over the grid, waves; every point of a must
/// have its counterpart in b (the grid is symmetric).
double max_abs_dw(const rtt::analysis::OpdMap& a, const rtt::analysis::OpdMap& b, double label) {
  double m = 0.0;
  for (const auto& q : a.points) {
    REQUIRE(q.status == RayStatus::Alive);
    const auto it = std::find_if(b.points.begin(), b.points.end(), [&](const auto& r) {
      return r.px == label * q.px && r.py == label * q.py;
    });
    REQUIRE(it != b.points.end());
    REQUIRE(it->status == RayStatus::Alive);
    m = std::max(m, std::abs(q.w - it->w));
  }
  return m;
}

}  // namespace

TEST_CASE("double-telecentric 4f relay: OPD against the reference at infinity is 0 on axis (#102)",
          "[opd][telecentric]") {
  // tests/reference/m2/telecentric_4f.rtt.json: two Descartes lenses (n = 1.5, R = +-64,
  // K = -2.25, d = 3, f = 128), object at z = -128 (front focus of L1), stop at z = 129, image
  // at z = 386. Front matrix to the stop: y = 1 - 2/128 = 63/64 at L1.S2 and 63/64 - 126/128 = 0
  // at the stop (a = 0), back matrix: nu' = 1 - 128/128 = 0 (d = 0), both exact in binary:
  // entrance and exit pupil at infinity. On axis the relay images without aberration (Sasian
  // L14, p. 5), so W = 0 over the whole pupil.
  // Tolerance 1e-7 waves: the OPL over about 520 mm and six segments carries a few roundings of
  // 1e-13 mm each, about 1e-12 mm = 2e-9 waves; the term n' d . (p - C) is 0 up to the same
  // rounding (every ray meets C); the aiming residual (< 1e-9 mm at the stop) enters W only in
  // second order (Fermat). Margin about 50.
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(load("m2/telecentric_4f.rtt.json"), lib);
  const auto fo = rtt::paraxial::first_order(cs, PathId{0}, 0);
  REQUIRE(fo.entrance_pupil);
  REQUIRE_FALSE(fo.entrance_pupil->z);
  REQUIRE(fo.exit_pupil);
  REQUIRE_FALSE(fo.exit_pupil->z);

  const auto map = rtt::analysis::opd_map(cs, PathId{0}, 0, 0);
  REQUIRE(std::isinf(map.sphere.radius));
  REQUIRE(map.vignetted == 0);
  REQUIRE(map.arrived == map.points.size());
  REQUIRE(max_abs_w(map) <= 1e-7);

  // Off axis the relay is not perfect, but fans and map run and W(0, 0) = 0.
  const auto fan = rtt::analysis::opd_fan(cs, PathId{0}, 1, 0);
  for (const auto& q : fan.tangential) {
    REQUIRE(q.status == RayStatus::Alive);
    if (q.py == 0.0) REQUIRE(q.w == 0.0);
  }
  const auto off = rtt::analysis::opd_map(cs, PathId{0}, 1, 0);
  REQUIRE(off.vignetted == 0);
  REQUIRE(std::isfinite(off.rms));
}

TEST_CASE("image-space telecentric singlet: defocus W = eps (1 - cos U) exactly (#102)",
          "[opd][telecentric]") {
  // Descartes lens of the relay (plane side first, hyperbola R = -64, K = -2.25 at z = 3),
  // object at infinity, stop (EPD 20) in the front focal plane z = -126 (H = V1 + 2, f = 128):
  // the exit pupil is at infinity (back matrix: y = 126, then 128 at L.S2, nu' = 1 - 1 = 0).
  // Focus F at z = 3 + 128 = 131. With the image plane at z_F + eps every ray passes F with the
  // same OPL (perfect focus) and reaches the plane after eps / cos U more; C = F + eps z_hat and
  // p - C = eps (d / cos U - z_hat), so d . (p - C) = eps (1 / cos U - cos U) and
  // W_inf = (OPL_F + eps) - (OPL_F + eps / cos U) + eps (1 / cos U - cos U) = eps (1 - cos U).
  // To first order this is 1/2 eps sin^2 U, Wyant & Creath Eq. (18), with the sign of Sec. I
  // (eps > 0 along the propagation: W > 0). U from the ray at height h on the hyperbola
  // (collimated, normal incidence on the plane side): z_h = 3 + sag(h), cos U = (131 - z_h) /
  // |(h, 131 - z_h)|. Tolerance 1e-7 waves (W about 2.6 waves; roundings as above).
  for (const double eps : {0.5, -0.5}) {
    System s;
    s.name = "image-space telecentric Descartes lens";
    s.environment.medium = "VACUUM";
    s.wavelengths = {{0.5876, 1.0, true}};
    s.aperture = {rtt::model::SystemApertureType::EntrancePupilDiameter, rtt::model::Param(20.0)};
    s.fields = {rtt::model::FieldType::AngleDeg, {{0.0, 0.0, 1.0}}};
    s.root.name = "root";
    rtt::model::Surface sto;
    sto.id = rtt::model::SurfaceId("STO");
    sto.aperture = rtt::model::CircularAperture{10.0, 0.0};
    rtt::model::Surface s1;
    s1.id = rtt::model::SurfaceId("L.S1");
    rtt::model::Surface s2;
    s2.id = rtt::model::SurfaceId("L.S2");
    s2.pose = Pose::along_z(3.0);
    s2.shape.base = rtt::model::Conic{rtt::model::Param(-64.0), rtt::model::Param(-2.25)};
    rtt::model::Surface img;
    img.id = rtt::model::SurfaceId("IMG");
    s.root.children = {
        {rtt::model::Element{
            "stop", rtt::model::ElementKind::Stop, Pose::along_z(-126.0), std::nullopt, {sto}}},
        {rtt::model::Element{
            "L", rtt::model::ElementKind::Lens, Pose::along_z(0.0), "CONST:1.5", {s1, s2}}},
        {rtt::model::Element{"image",
                             rtt::model::ElementKind::Detector,
                             Pose::along_z(131.0 + eps),
                             std::nullopt,
                             {img}}}};
    s.paths = {{"main", true, {}}};
    const MaterialLibrary lib;
    const CompiledSystem cs = compile(s, lib);
    const auto fo = rtt::paraxial::first_order(cs, PathId{0}, 0);
    REQUIRE(fo.exit_pupil);
    REQUIRE_FALSE(fo.exit_pupil->z);

    rtt::analysis::OpdOptions options;
    options.fan_points = 21;
    const auto fan = rtt::analysis::opd_fan(cs, PathId{0}, 0, 0, options);
    REQUIRE(std::isinf(fan.sphere.radius));
    int checked = 0;
    for (const auto& q : fan.tangential) {
      REQUIRE(q.status == RayStatus::Alive);
      const double h = 10.0 * std::abs(q.py);
      const double c = -1.0 / 64.0;
      const double sag = c * h * h / (1.0 + std::sqrt(1.0 - (1.0 - 2.25) * c * c * h * h));
      const double dz = 131.0 - (3.0 + sag);
      const double cos_u = dz / std::hypot(h, dz);
      const double expected = eps * (1.0 - cos_u) / kLambdaMm;
      INFO("eps " << eps << ", py " << q.py << ", W " << q.w << ", expected " << expected);
      REQUIRE(std::abs(q.w - expected) <= 1e-7);
      ++checked;
    }
    REQUIRE(checked == 21);
  }
}

TEST_CASE("finite, distant exit pupil converges to the reference at infinity (#102)",
          "[opd][telecentric]") {
  // The 4f relay with the stop at 129 + delta: the exit pupil is finite but far away (about
  // 128^2 / delta), the entrance pupil too. W(delta) - W_inf is O(delta): the bundle changes
  // smoothly with the stop (labels mirrored, (px, py) -> (-px, -py), when the EP lies on the
  // far side of the object, #96) and the reference term n' (b^2 - |p - C|^2) / (2 R) is
  // O(1 / R) = O(delta). Checked off axis (field 1, aberrated) as linear convergence:
  // max |dW(1e-4)| <= 2e-2 max |dW(1e-2)| (ratio 1e-2, quadratic remainder <= 1e-2 relative)
  // and, for delta = 1e-10 (R about 1.6e14 mm), max |dW| <= 2e-8 max |dW(1e-2)| + 1e-7 waves:
  // there the OPL to the sphere must not be formed as OPL - n' s with s about R: one ulp of R
  // (ulp(1e14) = 0.016 mm) is up to about 30 waves; with that old form the test measured
  // max |dW| = 0.88 waves at delta = 1e-10 and 6.4e-5 waves already at delta = 1e-4.
  const MaterialLibrary lib;
  const CompiledSystem tele = compile(load("m2/telecentric_4f.rtt.json"), lib);
  const auto map_t = rtt::analysis::opd_map(tele, PathId{0}, 1, 0);
  for (const double sign : {1.0, -1.0}) {
    double dw_coarse = 0.0;
    for (const double magnitude : {1e-2, 1e-4, 1e-10}) {
      const double delta = sign * magnitude;
      INFO("delta = " << delta);
      System s = load("m2/telecentric_4f.rtt.json");
      element(s, 1).pose = Pose::along_z(129.0 + delta);
      const CompiledSystem near = compile(s, lib);
      const auto fo = rtt::paraxial::first_order(near, PathId{0}, 0);
      REQUIRE(fo.exit_pupil);
      REQUIRE(fo.exit_pupil->z);
      REQUIRE(fo.entrance_pupil);
      REQUIRE(fo.entrance_pupil->z);
      const double label = *fo.entrance_pupil->z < -128.0 ? -1.0 : 1.0;
      const auto map_n = rtt::analysis::opd_map(near, PathId{0}, 1, 0);
      REQUIRE(std::isfinite(map_n.sphere.radius));
      const double dw = max_abs_dw(map_n, map_t, label);
      INFO("max |dW| = " << dw << " waves");
      if (magnitude == 1e-2) {
        dw_coarse = dw;
        REQUIRE(dw > 0.0);
      } else if (magnitude == 1e-4) {
        REQUIRE(dw <= 2e-2 * dw_coarse);
      } else {
        REQUIRE(dw <= 2e-8 * dw_coarse + 1e-7);
      }
    }
  }
}
