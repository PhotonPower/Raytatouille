#include <array>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <utility>
#include <vector>

#include "rtt/geom/asphere.hpp"
#include "rtt/geom/conic.hpp"
#include "rtt/geom/plane.hpp"
#include "rtt/geom/shape.hpp"

// ADR 0041, point 1: every shape gives the Taylor expansion of its sag at the vertex, and a stack
// of a base and additive terms sums the expansions.

using Catch::Matchers::WithinRel;
using rtt::geom::Conic;
using rtt::geom::EvenAsphere;
using rtt::geom::ParaxialShape;
using rtt::geom::Plane;
using rtt::geom::stack_expansion;
using rtt::geom::VertexExpansion;

static_assert(ParaxialShape<Plane<double>, double>);
static_assert(ParaxialShape<Conic<double>, double>);
static_assert(ParaxialShape<EvenAsphere<double>, double>);

namespace {

/// Every value of an expansion that the paraxial rule reads (ADR 0041, point 2).
void require_rotational(const VertexExpansion<double>& e, double c, double conic, double a4) {
  CHECK(e.value == 0.0);
  CHECK(e.slope == std::pair{0.0, 0.0});
  CHECK(e.hessian == std::array{c, 0.0, c});
  CHECK(e.order == 4);
  CHECK(e.rotational4);
  CHECK(e.conic == conic);
  CHECK(e.a4 == a4);
}

/// Second derivative of the meridional sag at the vertex by central differences, an independent
/// check of the closed curvature: z(h) - 2 z(0) + z(-h) = z'' h^2 + O(h^4).
template <class S>
double numeric_curvature(const S& shape, double h) {
  return (shape.sag(0.0, h) - 2.0 * shape.sag(0.0, 0.0) + shape.sag(0.0, -h)) / (h * h);
}

}  // namespace

TEST_CASE("the plane has a flat vertex expansion", "[geom][vertex]") {
  require_rotational(Plane<double>{}.vertex_expansion(), 0.0, 0.0, 0.0);
}

TEST_CASE("a conic gives its curvature and conic constant at the vertex", "[geom][vertex]") {
  const Conic<double> conic(1.0 / 25.0, -0.5);
  const VertexExpansion<double> e = conic.vertex_expansion();
  require_rotational(e, 1.0 / 25.0, -0.5, 0.0);
  // Independent check against the sag (rounding of the difference quotient ~1e-9 relative).
  CHECK_THAT(numeric_curvature(conic, 1e-3), WithinRel(e.hessian[0], 1e-6));
  CHECK(e.value == conic.sag(0.0, 0.0));
}

TEST_CASE("an even asphere adds A4 to the cap but not to the curvature", "[geom][vertex]") {
  const double c = -1.0 / 40.0;
  const EvenAsphere<double> plain(c, 0.3, {});
  const EvenAsphere<double> aspheric(c, 0.3, {2e-6, -3e-9});
  require_rotational(plain.vertex_expansion(), c, 0.3, 0.0);
  require_rotational(aspheric.vertex_expansion(), c, 0.3, 2e-6);
  // Content guard: A4 changes a4, not c or the vertex position.
  CHECK(aspheric.vertex_expansion().a4 != plain.vertex_expansion().a4);
  CHECK(aspheric.vertex_expansion().hessian == plain.vertex_expansion().hessian);
  CHECK_THAT(numeric_curvature(aspheric, 1e-3), WithinRel(c, 1e-6));
}

TEST_CASE("a stack without terms keeps the expansion of its base", "[geom][vertex]") {
  const VertexExpansion<double> base = EvenAsphere<double>(0.02, -1.2, {5e-7}).vertex_expansion();
  const VertexExpansion<double> s = stack_expansion<double>(base, {});
  CHECK(s.value == base.value);
  CHECK(s.slope == base.slope);
  CHECK(s.hessian == base.hessian);
  CHECK(s.order == base.order);
  CHECK(s.rotational4 == base.rotational4);
  // Bitwise the (K, A4) of the base, not the equivalent K = 0 form.
  CHECK(s.conic == -1.2);
  CHECK(s.a4 == 5e-7);
}

TEST_CASE("a stack sums the expansions; its cap is K = 0 with A4 = b4 - c^3 / 8",
          "[geom][vertex]") {
  const double c0 = 0.02;
  const double k0 = -0.5;
  const double a0 = 1e-6;
  const VertexExpansion<double> base = EvenAsphere<double>(c0, k0, {a0}).vertex_expansion();
  // A test term (only here, not in the model): z = 0.1 + 0.005 r^2 + 3e-7 r^4, so c = 0.01.
  VertexExpansion<double> term;
  term.value = 0.1;
  term.hessian = {0.01, 0.0, 0.01};
  term.order = 4;
  term.rotational4 = true;
  term.conic = 0.0;
  term.a4 = 3e-7 - 0.01 * 0.01 * 0.01 / 8.0;  // b4 = 3e-7 with K = 0
  const std::vector<VertexExpansion<double>> terms{term};
  const VertexExpansion<double> s = stack_expansion<double>(base, terms);

  CHECK(s.value == 0.1);
  CHECK(s.slope == std::pair{0.0, 0.0});
  CHECK(s.hessian[0] == c0 + 0.01);
  CHECK(s.hessian[1] == 0.0);
  CHECK(s.hessian[2] == s.hessian[0]);
  CHECK(s.order == 4);
  CHECK(s.rotational4);
  CHECK(s.conic == 0.0);
  // b4 of the sum: (1 + K) c0^3 / 8 + A4 of the base plus 3e-7 of the term (ADR 0041, point 1).
  const double c = c0 + 0.01;
  const double b4 = (1.0 + k0) * c0 * c0 * c0 / 8.0 + a0 + 3e-7;
  CHECK_THAT(s.a4, WithinRel(b4 - c * c * c / 8.0, 1e-12));
  // Content guard: the stack differs from its base.
  CHECK(s.a4 != base.a4);
}

TEST_CASE("a stack takes the lowest order and needs every part rotational", "[geom][vertex]") {
  const VertexExpansion<double> base = Conic<double>(0.01, 0.0).vertex_expansion();
  VertexExpansion<double> astigmatic;
  astigmatic.hessian = {0.002, 0.0, -0.002};
  astigmatic.order = 4;
  astigmatic.rotational4 = false;
  VertexExpansion<double> cubic;
  cubic.hessian = {0.0, 0.0, 0.0};
  cubic.order = 2;
  cubic.rotational4 = true;
  const std::vector<VertexExpansion<double>> terms{astigmatic, cubic};
  const VertexExpansion<double> s = stack_expansion<double>(base, terms);
  CHECK(s.order == 2);
  CHECK_FALSE(s.rotational4);
  CHECK(s.hessian[0] == 0.01 + 0.002);
  CHECK(s.hessian[2] == 0.01 + -0.002);
}
