// The paraxial rule for the XY polynomial (#316, ADR 0051, points 4 and 9; ADR 0041, points 1
// and 2): with c_x2 = c_y2 first order is that of the sphere with the same curvature, bit for
// bit when 1 / radius == 2 c / (R R); the term alone is the paraboloid z = c (x^2 + y^2) / R^2,
// so Seidel agrees with the conic k = -1 (cap 8 b4 - c^3 = -c^3 = K c^3), and with the fourth-order
// part b (x^2 + y^2)^2 / R^4 (c_x4 = c_y4 = b, c_x2y2 = 2 b) with the even asphere k = -1,
// A4 = b / R^4, which is the same surface. c_x2 != c_y2 is a ParaxialError.
//
// Tolerance of S_I: both sides form c and the cap from the same values in a different order
// (XY: 2 c / (R R) and b / R^4 - c^3 / 8; conic or asphere: 1 / radius and K c^3 + 8 A4), so they
// differ by a few ulp of the summands of S_I; 1e-10 relative covers that with a wide margin. The
// content guards differ by more than 1e-6 relative.

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

rtt::model::ShapeStack xy(double radius, double cx2, double cy2) {
  rtt::model::ShapeStack s;
  s.terms = {rtt::model::XyPolynomial{
      Param(radius), {Param(0.0), Param(0.0), Param(cx2), Param(0.0), Param(cy2)}}};
  return s;
}

}  // namespace

TEST_CASE("an XY polynomial with c_x2 = c_y2 has the first order of the sphere", "[paraxial][xy]") {
  // R = 10, c = 0.0625: 2 c / (R R) = 1/800 exactly, the sphere radius 800 mm.
  constexpr double kR = 10.0;
  constexpr double kC = 0.0625;
  constexpr double kRadius = 800.0;
  REQUIRE(1.0 / kRadius == 2.0 * kC / (kR * kR));
  rtt::model::ShapeStack sphere;
  sphere.base = rtt::model::Conic{Param(kRadius), Param(0.0)};
  const CompiledSystem cx = compiled(singlet(xy(kR, kC, kC)));
  const CompiledSystem cs = compiled(singlet(sphere));
  const auto fx = rtt::paraxial::first_order(cx, PathId{0}, 1);
  const auto fs = rtt::paraxial::first_order(cs, PathId{0}, 1);
  REQUIRE(fx.efl.has_value());
  REQUIRE(fs.efl.has_value());
  CHECK(fx.efl.value_or(0.0) == fs.efl.value_or(1.0));
  CHECK(fx.bfl == fs.bfl);
  // Seidel: the term is the paraboloid, the conic k = -1 with the same radius.
  rtt::model::ShapeStack paraboloid;
  paraboloid.base = rtt::model::Conic{Param(kRadius), Param(-1.0)};
  const auto sx = rtt::paraxial::seidel(cx, PathId{0}, 1);
  const auto sp = rtt::paraxial::seidel(compiled(singlet(paraboloid)), PathId{0}, 1);
  CHECK(std::abs(sx.sum.s1 - sp.sum.s1) <= 1e-10 * std::abs(sp.sum.s1));
  // Content guard: the sphere's S_I differs (its cap is 0, not -c^3).
  const auto ss = rtt::paraxial::seidel(cs, PathId{0}, 1);
  CHECK(std::abs(ss.sum.s1 - sp.sum.s1) > 1e-6 * std::abs(sp.sum.s1));
  // Content guard: L1.S2 is a plane, so the term alone gives the singlet its power; without it
  // the singlet is afocal.
  const auto fp =
      rtt::paraxial::first_order(compiled(singlet(rtt::model::ShapeStack{})), PathId{0}, 1);
  CHECK_FALSE(fp.efl.has_value());
}

TEST_CASE("the fourth-order XY part enters Seidel like A4 of the even asphere", "[paraxial][xy]") {
  // c (x^2 + y^2) / R^2 + b (x^2 + y^2)^2 / R^4 = r^2 / (2 R_p) + (b / R^4) r^4 with
  // R_p = R^2 / (2 c): the even asphere with k = -1 (paraboloid) and A4 = b / R^4.
  constexpr double kR = 10.0;
  constexpr double kC = 0.0625;
  constexpr double kRadius = 800.0;
  constexpr double kB = 0.5;
  const auto stack = [&](double b) {
    rtt::model::ShapeStack s;
    std::vector<Param> c(14, Param(0.0));
    c[2] = Param(kC);        // x^2
    c[4] = Param(kC);        // y^2
    c[9] = Param(b);         // x^4
    c[11] = Param(2.0 * b);  // x^2 y^2
    c[13] = Param(b);        // y^4
    s.terms = {rtt::model::XyPolynomial{Param(kR), c}};
    return s;
  };
  rtt::model::ShapeStack asphere;
  asphere.base =
      rtt::model::EvenAsphere{Param(kRadius), Param(-1.0), {Param(kB / (kR * kR * kR * kR))}};
  const auto sx = rtt::paraxial::seidel(compiled(singlet(stack(kB))), PathId{0}, 1);
  const auto sa = rtt::paraxial::seidel(compiled(singlet(asphere)), PathId{0}, 1);
  CHECK(std::abs(sx.sum.s1 - sa.sum.s1) <= 1e-10 * std::abs(sa.sum.s1));
  // Content guard: without the fourth-order part S_I is that of the paraboloid, clearly different.
  const auto s0 = rtt::paraxial::seidel(compiled(singlet(stack(0.0))), PathId{0}, 1);
  CHECK(std::abs(s0.sum.s1 - sa.sum.s1) > 1e-6 * std::abs(sa.sum.s1));
}

TEST_CASE("an XY polynomial with c_x2 != c_y2 is a ParaxialError (ADR 0041, point 2)",
          "[paraxial][xy]") {
  const CompiledSystem cs = compiled(singlet(xy(10.0, 0.0625, 0.03125)));
  CHECK_THROWS_AS((void)rtt::paraxial::first_order(cs, PathId{0}, 1), rtt::paraxial::ParaxialError);
}
