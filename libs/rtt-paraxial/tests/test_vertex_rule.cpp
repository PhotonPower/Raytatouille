#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <optional>
#include <string>
#include <utility>

#include "../src/vertex_rule.hpp"
#include "rtt/geom/asphere.hpp"
#include "rtt/geom/conic.hpp"
#include "rtt/geom/plane.hpp"
#include "rtt/geom/shape.hpp"

// ADR 0041, point 2: the rule for the rotationally symmetric paraxial model of one surface, a
// pure function over the vertex expansion. The test shapes exist only here, not in the model.

using rtt::geom::VertexExpansion;
using rtt::paraxial::detail::vertex_rule;
using rtt::paraxial::detail::vertex_z;
using rtt::paraxial::detail::VertexRule;

namespace {

/// Smooth, rotationally symmetric expansion with curvature c and cap (K, A4).
VertexExpansion<double> smooth(double c, double conic = 0.0, double a4 = 0.0) {
  VertexExpansion<double> e;
  e.hessian = {c, 0.0, c};
  e.order = 4;
  e.rotational4 = true;
  e.conic = conic;
  e.a4 = a4;
  return e;
}

/// A shape whose starting value of the intersection (base_conic) is not its vertex curvature,
/// as for Forbes Q-bfs or a Zernike defocus term (ADR 0041, Kontext).
struct OffsetStart {
  [[nodiscard]] std::pair<double, double> base_conic() const { return {0.5, -1.0}; }
  [[nodiscard]] VertexExpansion<double> vertex_expansion() const noexcept { return smooth(0.01); }
};

void require_first_order_error(const VertexRule& r, const std::string& part) {
  REQUIRE(r.first_order_error.has_value());
  CHECK(r.first_order_error->find(part) != std::string::npos);
  // Without a first-order model there is no fourth-order model either.
  CHECK(r.fourth_order_error.has_value());
}

}  // namespace

TEST_CASE("the existing shapes keep their curvature and cap", "[paraxial][vertex]") {
  const VertexRule plane = vertex_rule(rtt::geom::Plane<double>{}.vertex_expansion());
  CHECK_FALSE(plane.first_order_error);
  CHECK_FALSE(plane.fourth_order_error);
  CHECK(plane.c == 0.0);
  CHECK(plane.value == 0.0);

  const VertexRule conic = vertex_rule(rtt::geom::Conic<double>(-0.04, 2.0).vertex_expansion());
  CHECK_FALSE(conic.first_order_error);
  CHECK_FALSE(conic.fourth_order_error);
  CHECK(conic.c == -0.04);
  CHECK(conic.conic == 2.0);
  CHECK(conic.a4 == 0.0);

  const VertexRule asphere =
      vertex_rule(rtt::geom::EvenAsphere<double>(0.03, -0.7, {4e-6, 1e-9}).vertex_expansion());
  CHECK_FALSE(asphere.first_order_error);
  CHECK_FALSE(asphere.fourth_order_error);
  CHECK(asphere.c == 0.03);
  CHECK(asphere.conic == -0.7);
  CHECK(asphere.a4 == 4e-6);
}

TEST_CASE("the paraxial curvature comes from the Hessian, not from base_conic",
          "[paraxial][vertex]") {
  const OffsetStart shape;
  REQUIRE(shape.base_conic().first != 0.01);
  const VertexRule r = vertex_rule(shape.vertex_expansion());
  CHECK_FALSE(r.first_order_error);
  CHECK(r.c == 0.01);
}

TEST_CASE("unequal diagonals of the Hessian have no rotationally symmetric model",
          "[paraxial][vertex]") {
  VertexExpansion<double> toric = smooth(0.02);
  toric.hessian = {0.02, 0.0, 0.025};
  require_first_order_error(vertex_rule(toric), "not rotationally symmetric in second order");

  VertexExpansion<double> twisted = smooth(0.02);
  twisted.hessian = {0.02, 1e-6, 0.02};
  require_first_order_error(vertex_rule(twisted), "not rotationally symmetric in second order");
}

TEST_CASE("a slope at the vertex has no rotationally symmetric model", "[paraxial][vertex]") {
  for (const auto& slope : {std::pair{1e-4, 0.0}, std::pair{0.0, -1e-4}}) {
    VertexExpansion<double> tilted = smooth(0.02);
    tilted.slope = slope;
    require_first_order_error(vertex_rule(tilted), "slope at the vertex");
  }
}

TEST_CASE("an expansion below second order has no paraxial model", "[paraxial][vertex]") {
  for (const int order : {0, 1}) {
    VertexExpansion<double> cone = smooth(0.0);
    cone.order = order;
    require_first_order_error(vertex_rule(cone), "no paraxial model at the vertex");
  }
}

TEST_CASE("the order is checked before the slope, the slope before the Hessian",
          "[paraxial][vertex]") {
  // Slope and Hessian have no meaning for order 0 (ADR 0041, point 1), e.g. an odd asphere with
  // A1 != 0: the message names the missing paraxial model, not the slope or the asymmetry.
  VertexExpansion<double> cone = smooth(0.02);
  cone.order = 0;
  cone.slope = {1e-3, 0.0};
  cone.hessian = {0.02, 0.0, 0.03};
  require_first_order_error(vertex_rule(cone), "no paraxial model at the vertex");

  VertexExpansion<double> tilted_toric = smooth(0.02);
  tilted_toric.slope = {0.0, 1e-3};
  tilted_toric.hessian = {0.02, 0.0, 0.03};
  require_first_order_error(vertex_rule(tilted_toric), "slope at the vertex");
}

TEST_CASE("the Seidel cap needs order 4 and rotational symmetry to fourth order",
          "[paraxial][vertex]") {
  // |r|^3 (A3 != 0): order 2, first order works, Seidel does not.
  VertexExpansion<double> cubic = smooth(0.02);
  cubic.order = 2;
  const VertexRule r_cubic = vertex_rule(cubic);
  CHECK_FALSE(r_cubic.first_order_error);
  CHECK(r_cubic.c == 0.02);
  REQUIRE(r_cubic.fourth_order_error.has_value());
  CHECK(r_cubic.fourth_order_error->find("fourth order") != std::string::npos);

  // A fourth-order term that is not rotationally symmetric (e.g. Zernike trefoil or x^4 alone).
  VertexExpansion<double> quartic = smooth(0.02);
  quartic.rotational4 = false;
  const VertexRule r_quartic = vertex_rule(quartic);
  CHECK_FALSE(r_quartic.first_order_error);
  REQUIRE(r_quartic.fourth_order_error.has_value());
  CHECK(r_quartic.fourth_order_error->find("rotationally symmetric") != std::string::npos);
}

TEST_CASE("the constant part shifts the vertex along the local z axis", "[paraxial][vertex]") {
  VertexExpansion<double> shifted = smooth(0.02);
  shifted.value = 0.25;
  const VertexRule r = vertex_rule(shifted);
  CHECK_FALSE(r.first_order_error);
  CHECK(r.value == 0.25);
  // Local z along global +z and along -z (a surface after a mirror or a flipped element).
  CHECK(vertex_z(10.0, 1.0, r.value) == 10.25);
  CHECK(vertex_z(10.0, -1.0, r.value) == 9.75);
  CHECK(vertex_z(10.0, 1.0, -0.5) == 9.5);
  CHECK(vertex_z(10.0, -1.0, -0.5) == 10.5);
}

TEST_CASE("without a constant part the vertex stays at the origin, also for -0",
          "[paraxial][vertex]") {
  CHECK(vertex_z(12.5, 1.0, 0.0) == 12.5);
  CHECK(vertex_z(12.5, -1.0, 0.0) == 12.5);
  for (const double axis_z : {1.0, -1.0}) {
    for (const double value : {0.0, -0.0}) {
      const double z = vertex_z(-0.0, axis_z, value);
      CHECK(z == 0.0);
      CHECK(std::signbit(z));
    }
  }
}
