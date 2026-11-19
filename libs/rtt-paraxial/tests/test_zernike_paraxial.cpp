// The paraxial rule for Zernike sag terms (#316, ADR 0051, point 6; ADR 0041, points 1 and 2):
// a rotationally symmetric term enters first order and Seidel like the even asphere with the
// same vertex curvature and cap, and an astigmatic term is a ParaxialError.

#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <string>
#include <variant>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/paraxial/paraxial.hpp"
#include "rtt/paraxial/seidel.hpp"

using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::model::Element;
using rtt::model::Param;
using rtt::model::System;

namespace {

/// m1/singlet_const with the shape of L1.S1 (/root/children/1/surfaces/0) replaced.
System singlet(const rtt::model::ShapeStack& shape) {
  System s = rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m1/singlet_const.rtt.json");
  std::get<Element>(s.root.children[1].value).surfaces[0].shape = shape;
  return s;
}

CompiledSystem compiled(const System& s) {
  const rtt::material::MaterialLibrary lib;
  return rtt::compile::compile(s, lib);
}

}  // namespace

TEST_CASE("a spherical Zernike term on a plane equals the even asphere of its expansion",
          "[paraxial][zernike]") {
  // Noll Z11 = sqrt(5) (6 rho^4 - 6 rho^2 + 1) with a, and Z1 = 1 with -sqrt(5) a, so that the
  // vertex stays at the origin (value = sqrt(5) a - sqrt(5) a = 0, the same expression as the
  // code). Then c = 2 sqrt(5) a (-6) / R^2 and b4 = sqrt(5) a 6 / R^4: the even asphere with
  // radius 1/c, k = 0 and A4 = b4 - c^3 / 8 has the same first order and the same cap (ADR 0041,
  // point 1: K c^3 + 8 A4 = 8 b4 - c^3).
  const double rad = 12.7;
  const double a = 2e-3;
  const double coef = std::sqrt(5.0) * a;
  std::vector<Param> z(11, Param(0.0));
  z[0] = Param(-coef);
  z[10] = Param(a);
  rtt::model::ShapeStack zernike;
  zernike.terms = {rtt::model::ZernikeSag{Param(rad), z}};
  const double c = 2.0 * coef * -6.0 / (rad * rad);
  const double b4 = coef * 6.0 / (rad * rad * rad * rad);
  rtt::model::ShapeStack even;
  even.base = rtt::model::EvenAsphere{Param(1.0 / c), Param(0.0), {Param(b4 - c * c * c / 8.0)}};
  const CompiledSystem cz = compiled(singlet(zernike));
  const CompiledSystem ce = compiled(singlet(even));
  const auto fz = rtt::paraxial::first_order(cz, PathId{0}, 1);
  const auto fe = rtt::paraxial::first_order(ce, PathId{0}, 1);
  REQUIRE(fz.efl.has_value());
  REQUIRE(fe.efl.has_value());
  const double efl_z = fz.efl.value_or(0.0);
  const double efl_e = fe.efl.value_or(0.0);
  // c and 1 / (1 / c) differ by rounding only.
  CHECK(std::abs(efl_z - efl_e) <= 1e-12 * std::abs(efl_e));
  const auto sz = rtt::paraxial::seidel(cz, PathId{0}, 1);
  const auto se = rtt::paraxial::seidel(ce, PathId{0}, 1);
  CHECK(std::abs(sz.sum.s1 - se.sum.s1) <= 1e-10 * std::abs(se.sum.s1));
  // Content guard: the term changes the cap, so S_I differs from the sphere of the same c.
  rtt::model::ShapeStack sphere;
  sphere.base = rtt::model::Conic{Param(1.0 / c), Param(0.0)};
  const auto ss = rtt::paraxial::seidel(compiled(singlet(sphere)), PathId{0}, 1);
  CHECK(std::abs(ss.sum.s1 - se.sum.s1) > 1e-6 * std::abs(se.sum.s1));
}

TEST_CASE("an astigmatic Zernike term is a ParaxialError (ADR 0041, point 2)",
          "[paraxial][zernike]") {
  rtt::model::ShapeStack astig;
  astig.terms = {rtt::model::ZernikeSag{
      Param(12.7), {Param(0.0), Param(0.0), Param(0.0), Param(0.0), Param(0.0), Param(1e-3)}}};
  const CompiledSystem cs = compiled(singlet(astig));
  CHECK_THROWS_AS((void)rtt::paraxial::first_order(cs, PathId{0}, 1), rtt::paraxial::ParaxialError);
}
