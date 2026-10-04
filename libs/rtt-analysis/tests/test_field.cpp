#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstdint>
#include <numbers>
#include <optional>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#include "rtt/analysis/field.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
#include "rtt/paraxial/paraxial.hpp"
#include "rtt/paraxial/seidel.hpp"

using rtt::compile::compile;
using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::material::MaterialLibrary;
using rtt::model::Element;
using rtt::model::ElementKind;
using rtt::model::Field;
using rtt::model::FieldType;
using rtt::model::Param;
using rtt::model::Pose;
using rtt::model::Surface;
using rtt::model::SurfaceId;
using rtt::model::System;

namespace {

constexpr double kDeg = std::numbers::pi / 180.0;

System load(const std::string& relative) {
  return rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/" + relative);
}

Element& element(System& s, std::size_t i) {
  return std::get<Element>(s.root.children[i].value);
}

Surface surface(const std::string& id, double z_mm = 0.0, std::optional<double> radius = {}) {
  Surface s;
  s.id = SurfaceId(id);
  s.pose = Pose::along_z(z_mm);
  if (radius) s.shape.base = rtt::model::Conic{Param(*radius), Param(0.0)};
  return s;
}

System base_system(const std::string& name) {
  System s;
  s.name = name;
  s.environment.medium = "VACUUM";
  s.wavelengths = {{0.5876, 1.0, true}};
  s.root.name = "root";
  s.paths = {{"main", true, {}}};
  return s;
}

/// Stop element at z with a circular aperture.
Element stop_at(double z, double radius) {
  Surface st = surface("STO");
  st.aperture = rtt::model::CircularAperture{radius, 0.0};
  return Element{"stop", ElementKind::Stop, Pose::along_z(z), std::nullopt, {st}};
}

/// Biconvex lens R = +-50 mm, 4 mm thick, constant index, front vertex at z.
Element biconvex(const std::string& name, double z) {
  return Element{name,
                 ElementKind::Lens,
                 Pose::along_z(z),
                 "CONST:1.5168",
                 {surface(name + ".S1", 0.0, 50.0), surface(name + ".S2", 4.0, -50.0)}};
}

/// Least-squares slope of log|v| over log x.
double fit_exponent(const std::vector<double>& x, const std::vector<double>& v) {
  double mx = 0.0, my = 0.0;
  for (std::size_t i = 0; i < x.size(); ++i) {
    mx += std::log(x[i]);
    my += std::log(std::abs(v[i]));
  }
  mx /= static_cast<double>(x.size());
  my /= static_cast<double>(x.size());
  double sxy = 0.0, sxx = 0.0;
  for (std::size_t i = 0; i < x.size(); ++i) {
    sxy += (std::log(x[i]) - mx) * (std::log(std::abs(v[i])) - my);
    sxx += (std::log(x[i]) - mx) * (std::log(x[i]) - mx);
  }
  return sxy / sxx;
}

}  // namespace

TEST_CASE("pinhole camera: distortion exactly zero", "[distortion]") {
  // Stop at z = 0 and image plane at z = 50 in vacuum: the chief ray of field angle theta is the
  // straight line through the stop centre, h_real = 50 tan(theta). The paraxial chief ray with
  // slope u = tan(theta) (field convention of #8) gives h_par = 50 tan(theta) in the same plane.
  System s = base_system("pinhole");
  s.aperture = {rtt::model::SystemApertureType::EntrancePupilDiameter, Param(2.0)};
  s.fields = {FieldType::AngleDeg, {{0.0, 0.0, 1.0}, {0.0, 10.0, 1.0}, {0.0, 20.0, 1.0}}};
  s.root.children = {
      {stop_at(0.0, 1.0)},
      {Element{"D", ElementKind::Detector, Pose::along_z(50.0), std::nullopt, {surface("IMG")}}}};
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(s, lib);
  rtt::analysis::FieldSweepOptions options;
  options.samples = 5;
  const auto d = rtt::analysis::distortion(cs, PathId{0}, 0, options);
  REQUIRE(d.size() == 5);
  for (const auto& p : d) {
    INFO("fraction " << p.fraction);
    REQUIRE(std::abs(p.paraxial_height - 50.0 * std::tan(p.field.y * kDeg)) <= 1e-12);
    REQUIRE(std::abs(p.percent) <= 1e-12);
  }
  REQUIRE(d.front().percent == 0.0);
  REQUIRE(std::abs(d.back().field.y - 20.0) <= 1e-15);
  REQUIRE(std::abs(d[2].fraction - 0.5) <= 1e-15);
}

TEST_CASE("symmetric 1:1 system: distortion zero", "[distortion]") {
  // Two identical biconvex lenses mirror-symmetric about the stop at z = 0, object in the front
  // focal plane of the first lens, image in the mirrored plane: the system is symmetric with
  // m = -1, so every chief ray leaves symmetric to its entry and h_real = -h = h_par.
  System half = base_system("front lens");
  half.aperture = {rtt::model::SystemApertureType::EntrancePupilDiameter, Param(4.0)};
  half.fields = {FieldType::AngleDeg, {{0.0, 0.0, 1.0}}};
  half.root.children = {{biconvex("L1", -20.0)}, {stop_at(0.0, 5.0)}};
  const MaterialLibrary lib;
  const double z_front =
      *rtt::paraxial::first_order(compile(half, lib), PathId{0}, 0).front_focal_z;

  System s = base_system("symmetric 1:1");
  s.object.at_infinity = false;
  s.object.distance = Param(-z_front);
  s.aperture = {rtt::model::SystemApertureType::StopSize, Param(0.0)};
  s.fields = {FieldType::ObjectHeight, {{0.0, 0.0, 1.0}, {0.0, 3.0, 1.0}, {0.0, 6.0, 1.0}}};
  s.root.children = {
      {biconvex("L1", -20.0)},
      {stop_at(0.0, 2.0)},
      {biconvex("L2", 16.0)},
      {Element{
          "D", ElementKind::Detector, Pose::along_z(-z_front), std::nullopt, {surface("IMG")}}}};
  const CompiledSystem cs = compile(s, lib);
  const auto fo = rtt::paraxial::first_order(cs, PathId{0}, 0);
  REQUIRE(std::abs(*fo.lateral_magnification + 1.0) <= 1e-12);
  for (std::uint16_t f = 1; f < 3; ++f) {
    const auto p = rtt::analysis::distortion_at(cs, PathId{0}, cs.fields().points[f], 0);
    INFO("field " << f << ", D = " << p.percent << " %");
    REQUIRE(std::abs(p.paraxial_height + cs.fields().points[f].y) <= 1e-9);
    REQUIRE(std::abs(p.percent) <= 1e-9);
  }
}

TEST_CASE("field curvature on axis: tangential equals sagittal", "[field-curvature]") {
  // Rotational symmetry: the tangential pair (0, +-delta) and the sagittal pair (+-delta, 0) are
  // the same rays rotated by 90 deg. With the image surface in the paraxial focus both foci
  // differ from it only by the spherical aberration of the neighbour rays, O(delta^2).
  System s = load("m1/singlet_const.rtt.json");
  const MaterialLibrary lib;
  const CompiledSystem first = compile(s, lib);
  element(s, 2).pose = Pose::along_z(*rtt::paraxial::first_order(first, PathId{0}, 1).rear_focal_z);
  const CompiledSystem cs = compile(s, lib);
  const auto p = rtt::analysis::field_curvature_at(cs, PathId{0}, Field{0.0, 0.0, 1.0}, 1);
  REQUIRE(std::abs(p.tangential - p.sagittal) <= 1e-9);
  REQUIRE(std::abs(p.astigmatism) <= 1e-9);
  REQUIRE(std::abs(p.tangential) <= 1e-5);
}

TEST_CASE("field curvature at small fields agrees with the Seidel sums", "[field-curvature]") {
  // Third-order wavefront at field H (Sasian, OPTI 517 L4 p. 23, see docs/quellen.md):
  //   sagittal defocus W020_S = W220 H^2, tangential W020_T = (W220 + W222) H^2,
  //   W220 = (S_IV + S_III) / 4, W222 = S_III / 2,
  // with H relative to the maximum field of rtt::paraxial::seidel (here tan(theta) / tan(5 deg)).
  // A defocus coefficient W020 moves the focus by
  //   eps_z = -2 W020 / (n' u'^2)
  // along z: Wyant & Creath, Eq. (18), Delta W = -1/2 eps_z sin^2 U at the pupil rim (n' = 1),
  // with the optical path n' times the geometric one and sin U -> u' of the paraxial marginal ray
  // of seidel() (the ray that defines rho = 1 in Sasian's normalisation). The same follows from
  // the transverse relation of seidel.hpp, eps = (dW/drho) / (n' u'), checked against real rays
  // in #30: a ray at rho with image-space slope rho u' and eps = 2 W020 rho / (n' u') meets the
  // chief ray at Delta z = -eps / (rho u') = -2 W020 / (n' u'^2). n' u'^2 is positive here (no
  // mirror, light along +z), so Delta z > 0 means behind the paraxial image plane.
  // Real fields 0.5 and 1 deg (H = 0.1, 0.2): fifth-order terms are of relative size H^2 tan^2
  // (5 deg) <= 3.1e-4, the neighbour-ray error O(delta^2) is far below; tolerance relative 1e-2.
  System s = load("m1/singlet_const.rtt.json");
  const MaterialLibrary lib;
  const CompiledSystem first = compile(s, lib);
  element(s, 2).pose = Pose::along_z(*rtt::paraxial::first_order(first, PathId{0}, 1).rear_focal_z);
  const CompiledSystem cs = compile(s, lib);
  const auto seidel = rtt::paraxial::seidel(cs, PathId{0}, 1);
  const auto marginal = rtt::paraxial::trace_ray(cs, PathId{0}, 1, seidel.marginal.z,
                                                 seidel.marginal.y, seidel.marginal.u);
  const double nu2 = marginal.back().n * marginal.back().u * marginal.back().u;
  const double w220 = (seidel.sum.s4 + seidel.sum.s3) / 4.0;
  const double w222 = seidel.sum.s3 / 2.0;
  std::vector<double> angle;
  std::vector<double> tangential;
  for (const double theta : {0.25, 0.5, 1.0, 2.0}) {
    const auto p = rtt::analysis::field_curvature_at(cs, PathId{0}, Field{0.0, theta, 1.0}, 1);
    const double h = std::tan(theta * kDeg) / std::tan(5.0 * kDeg);
    const double expected_s = -2.0 * w220 * h * h / nu2;
    const double expected_t = -2.0 * (w220 + w222) * h * h / nu2;
    INFO("theta " << theta << ": T " << p.tangential << " (" << expected_t << "), S " << p.sagittal
                  << " (" << expected_s << ")");
    if (theta == 0.5 || theta == 1.0) {
      REQUIRE(std::abs(p.sagittal - expected_s) <= 1e-2 * std::abs(expected_s));
      REQUIRE(std::abs(p.tangential - expected_t) <= 1e-2 * std::abs(expected_t));
    }
    REQUIRE(std::abs(p.astigmatism - (p.tangential - p.sagittal)) <= 1e-15);
    angle.push_back(theta);
    tangential.push_back(p.tangential);
  }
  const double exponent = fit_exponent(angle, tangential);
  INFO("fit exponent " << exponent);
  REQUIRE(std::abs(exponent - 2.0) <= 0.05);
}

TEST_CASE("field sweeps run along +y of the largest field", "[distortion][field-curvature]") {
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(load("m1/singlet_const.rtt.json"), lib);
  rtt::analysis::FieldCurvatureOptions options;
  options.samples = 3;
  const auto fc = rtt::analysis::field_curvature(cs, PathId{0}, 1, options);
  REQUIRE(fc.size() == 3);
  REQUIRE(fc[0].field.y == 0.0);
  REQUIRE(fc[1].field.y == 2.5);
  REQUIRE(fc[2].field.y == 5.0);
  REQUIRE(fc[2].field.x == 0.0);
}

TEST_CASE("invalid field-analysis input", "[distortion][field-curvature]") {
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(load("m1/singlet_const.rtt.json"), lib);
  rtt::analysis::FieldSweepOptions sweep;
  sweep.samples = 1;
  REQUIRE_THROWS_AS(rtt::analysis::distortion(cs, PathId{0}, 1, sweep), std::invalid_argument);
  REQUIRE_THROWS_AS(rtt::analysis::distortion(cs, PathId{0}, 7), std::invalid_argument);
  rtt::analysis::FieldCurvatureOptions bad;
  bad.delta = 0.0;
  REQUIRE_THROWS_AS(rtt::analysis::field_curvature(cs, PathId{0}, 1, bad), std::invalid_argument);
  // Only an on-axis field point: no sweep direction.
  System s = load("m1/singlet_const.rtt.json");
  s.fields.points = {{0.0, 0.0, 1.0}};
  const CompiledSystem axis = compile(s, lib);
  REQUIRE_THROWS_AS(rtt::analysis::distortion(axis, PathId{0}, 1), std::invalid_argument);
}
