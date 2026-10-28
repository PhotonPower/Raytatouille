// object.medium (#282, ADR 0037): rays start in the medium of the object space. Reference cases
// with closed formulas: Snell's law n1 sin(theta1) = n2 sin(theta2) at the first surface (the
// same law and source as in test_sequential.cpp, "plane-parallel plate shifts the ray ...": de
// Greve, Reflections and Refractions in Ray Tracing, 2006, Eq. (23)), the optical path
// n * d / cos(theta) of the first leg (definition of the optical path length), and the paraxial
// reading of the object-space NA, tan U = NA / n (docs/quellen.md, row "Lagrange-Invariante, NA,
// Arbeits-F/#", Greivenkamp OPTI-502, Sec. 9; the same reading as the test in test_sources.cpp
// "object_na: paraxial reading"). The environment is VACUUM (AIR is Ciddor air since #25).

#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/sequential.hpp"
#include "rtt/trace/sources.hpp"

using rtt::compile::compile;
using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::material::MaterialLibrary;
using rtt::model::Element;
using rtt::model::ElementKind;
using rtt::model::Param;
using rtt::model::Pose;
using rtt::model::Surface;
using rtt::model::SurfaceId;
using rtt::model::System;
using rtt::trace::RayBatch;
using rtt::trace::RayStatus;
using rtt::trace::SequentialTracer;

namespace {

constexpr double kPi = 3.14159265358979323846;

Surface surface(const std::string& id, double z_mm = 0.0) {
  Surface s;
  s.id = SurfaceId(id);
  s.pose = Pose::along_z(z_mm);
  return s;
}

/// Finite object at `distance` mm in front of z = 0, in a medium of index 1.5, VACUUM around.
System glass_object(double distance) {
  System s;
  s.name = "object medium";
  s.environment.medium = "VACUUM";
  s.object.at_infinity = false;
  s.object.distance = Param(distance);
  s.object.medium = "CONST:1.5";
  s.wavelengths = {{0.5876, 1.0, true}};
  s.aperture = {rtt::model::SystemApertureType::EntrancePupilDiameter, Param(10.0)};
  s.fields = {rtt::model::FieldType::ObjectHeight, {{0.0, 0.0, 1.0}}};
  s.root.name = "root";
  s.paths = {{"main", true, {}}};
  return s;
}

/// Axial rays from (0, 0, z0) in the y-z plane at the given angles (rad) against +z.
RayBatch axial_rays(const std::vector<double>& angles, double z0) {
  RayBatch rays(angles.size());
  for (std::size_t i = 0; i < angles.size(); ++i) {
    rays.pos_x()[i] = 0.0;
    rays.pos_y()[i] = 0.0;
    rays.pos_z()[i] = z0;
    rays.dir_x()[i] = 0.0;
    rays.dir_y()[i] = std::sin(angles[i]);
    rays.dir_z()[i] = std::cos(angles[i]);
  }
  return rays;
}

}  // namespace

TEST_CASE("object medium: Snell's law at the first surface, glass to vacuum",
          "[sequential][object-medium]") {
  // A plate of VACUUM behind the object space: the first surface is the interface 1.5 -> 1, the
  // second one 1 -> 1 (no change), so the rays leave the plate at the refracted direction
  // sin(theta') = 1.5 sin(theta).
  System s = glass_object(2.0);
  s.root.children.push_back({Element{"P",
                                     ElementKind::Plate,
                                     Pose::along_z(0.0),
                                     "VACUUM",
                                     {surface("P.S1"), surface("P.S2", 5.0)}}});
  s.root.children.push_back(
      {Element{"D", ElementKind::Detector, Pose::along_z(40.0), std::nullopt, {surface("D")}}});
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(s, lib);

  const std::vector<double> degrees{10.0, 30.0, 45.0};  // 45 degrees: beyond the critical angle
  std::vector<double> angles;
  for (const double d : degrees) angles.push_back(d * kPi / 180.0);
  RayBatch rays = axial_rays(angles, -2.0);
  [[maybe_unused]] const auto stats = SequentialTracer().trace(cs, PathId{0}, rays);

  for (std::size_t i = 0; i < 2; ++i) {
    INFO("theta = " << degrees[i] << " deg");
    REQUIRE(rays.status()[i] == RayStatus::Alive);
    const double expected = 1.5 * std::sin(angles[i]);
    REQUIRE(std::abs(rays.dir_y()[i] - expected) <= 1e-12);
    // Content guard: the direction really changed (without the object medium it would not).
    REQUIRE(std::abs(rays.dir_y()[i] - std::sin(angles[i])) > 0.05);
  }
  // Reference values: asin(1.5 sin 10 deg) = 15.0980866 deg, asin(1.5 sin 30 deg) = 48.5903779 deg.
  REQUIRE(std::abs(std::asin(rays.dir_y()[0]) * 180.0 / kPi - 15.0980866) <= 1e-6);
  REQUIRE(std::abs(std::asin(rays.dir_y()[1]) * 180.0 / kPi - 48.5903779) <= 1e-6);
  // 1.5 sin 45 deg = 1.0607 > 1: total internal reflection at the first surface (critical angle
  // asin(1 / 1.5) = 41.8103 deg).
  REQUIRE(rays.status()[2] == RayStatus::Tir);
}

TEST_CASE("object medium: the first leg adds n * d / cos(theta) to the optical path",
          "[sequential][object-medium]") {
  // Only a detector at z = 0, 2 mm behind the start point: the ray travels 2 / cos(theta) in the
  // object medium of index 1.5. Without the object medium the path would be 2 / cos(theta).
  System s = glass_object(2.0);
  s.root.children.push_back(
      {Element{"D", ElementKind::Detector, Pose::along_z(0.0), std::nullopt, {surface("D")}}});
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(s, lib);

  const std::vector<double> angles{0.0, 30.0 * kPi / 180.0};
  RayBatch rays = axial_rays(angles, -2.0);
  [[maybe_unused]] const auto stats = SequentialTracer().trace(cs, PathId{0}, rays);
  REQUIRE(rays.status()[0] == RayStatus::Alive);
  REQUIRE(rays.status()[1] == RayStatus::Alive);
  REQUIRE(std::abs(rays.opl()[0] - 3.0) <= 1e-12);
  REQUIRE(std::abs(rays.opl()[1] - 3.4641016151377544) <= 1e-12);  // 1.5 * 2 / cos 30 deg
}

TEST_CASE("object medium: the object-space NA is read as tan U = NA / n",
          "[sources][aiming][object-medium]") {
  // Stop at z = 0 as the first element, object point at z = -20 in a medium of index 1.5,
  // NA = 0.9: the paraxial marginal slope is u = NA / n = 0.6, so the real marginal ray meets
  // the stop at r_s = 0.6 * 20 = 12 mm and has tan U = 0.6 (a straight line in the homogeneous
  // object space). In vacuum the same NA would give tan U = 0.9.
  System s = glass_object(20.0);
  s.aperture = {rtt::model::SystemApertureType::ObjectSpaceNA, Param(0.9)};
  Surface stop = surface("STO");
  stop.aperture = rtt::model::CircularAperture{30.0, 0.0};
  s.root.children.push_back(
      {Element{"stop", ElementKind::Stop, Pose::along_z(0.0), std::nullopt, {stop}}});
  s.root.children.push_back(
      {Element{"D", ElementKind::Detector, Pose::along_z(100.0), std::nullopt, {surface("IMG")}}});
  const MaterialLibrary lib;
  const CompiledSystem cs = compile(s, lib);

  const std::vector<std::uint16_t> axis{0};
  const RayBatch marginal = rtt::trace::make_rays(
      cs, PathId{0}, axis, 0, rtt::trace::SinglePupilPoint{0.0, 1.0}, rtt::trace::Aiming::Real);
  REQUIRE(std::abs(marginal.dir_y()[0] - 0.6 / std::sqrt(1.0 + 0.6 * 0.6)) <= 1e-11);
  REQUIRE(std::abs(marginal.dir_y()[0] - 0.9 / std::sqrt(1.0 + 0.9 * 0.9)) > 0.05);
}
