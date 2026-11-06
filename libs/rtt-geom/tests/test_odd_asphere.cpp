#include <algorithm>
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <utility>
#include <vector>

#include "rtt/geom/asphere.hpp"
#include "rtt/geom/conic.hpp"
#include "rtt/geom/intersect.hpp"
#include "rtt/geom/odd_asphere.hpp"
#include "rtt/geom/shape.hpp"

// Odd asphere of ADR 0043, point 1 (OpticStudio User Guide 2025 R1, "Odd Asphere"):
// z = c r^2 / (1 + phi) + sum_{n >= 1} A_n r^n with phi = sqrt(1 - (1 + k) c^2 r^2) (Forbes,
// Opt. Express 19(10), 9923-9942 (2011), Eqs. (2.1), (2.2)). Acceptance 6 of ADR 0043.

using rtt::geom::Conic;
using rtt::geom::DomainShape;
using rtt::geom::EvenAsphere;
using rtt::geom::HitStatus;
using rtt::geom::IntersectableShape;
using rtt::geom::OddAsphere;
using rtt::geom::ParaxialShape;
using rtt::geom::VertexExpansion;
using Vec3 = rtt::math::Vec3T<double>;

static_assert(IntersectableShape<OddAsphere<double>, double>);
static_assert(ParaxialShape<OddAsphere<double>, double>);
// The domain is the disk of the base conic, as for the even asphere (ADR 0043, point 1).
static_assert(!DomainShape<OddAsphere<double>, double>);

namespace {

constexpr double kEps = std::numeric_limits<double>::epsilon();

/// Every field of two expansions compared with == (ADR 0043, acceptance 6: bit for bit).
void require_same(const VertexExpansion<double>& a, const VertexExpansion<double>& b) {
  CHECK(a.value == b.value);
  CHECK(a.slope == b.slope);
  CHECK(a.hessian == b.hessian);
  CHECK(a.order == b.order);
  CHECK(a.rotational4 == b.rotational4);
  CHECK(a.conic == b.conic);
  CHECK(a.a4 == b.a4);
}

/// The odd asphere written out independently of the library: the conic sag plus sum A_n r^n,
/// term by term with std::pow, and its radial derivative c r / phi + sum n A_n r^(n - 1).
double sag_written_out(double c, double k, const std::vector<double>& a, double r) {
  double z = c * r * r / (1.0 + std::sqrt(1.0 - (1.0 + k) * c * c * r * r));
  for (std::size_t n = 1; n <= a.size(); ++n) z += a[n - 1] * std::pow(r, static_cast<double>(n));
  return z;
}

}  // namespace

TEST_CASE("odd asphere with only even terms equals the even asphere to 4 ulp at dyadic r",
          "[geom][odd_asphere]") {
  // ADR 0043, acceptance 6: Horner in r (odd asphere) against Horner in r^2 (even asphere);
  // c = 1/32 is exact. A2 = 0, so both have the base conic as vertex curvature.
  const double c = 1.0 / 32.0;
  const double k = -1.3;
  const EvenAsphere<double> even(c, k, {2e-6, -3e-9, 4e-12});
  const OddAsphere<double> odd(c, k, {0.0, 0.0, 0.0, 2e-6, 0.0, -3e-9, 0.0, 4e-12});
  double worst = 0.0;
  for (int i = 0; i <= 960; ++i) {
    const double r = i / 64.0;  // dyadic, up to 15 mm
    for (const auto& [x, y] : {std::pair{r, 0.0}, std::pair{0.0, -r}}) {
      const double ze = even.sag(x, y);
      const double zo = odd.sag(x, y);
      REQUIRE(std::isfinite(ze));
      const double scale = std::max(std::abs(ze), std::numeric_limits<double>::min());
      worst = std::max(worst, std::abs(ze - zo) / (kEps * scale));
    }
  }
  CHECK(worst <= 4.0);
  // Content guard: the polynomial part is not negligible at the rim.
  CHECK(std::abs(odd.sag(15.0, 0.0) - Conic<double>(c, k).sag(15.0, 0.0)) > 0.05);
}

TEST_CASE("odd asphere with only even terms has the vertex expansion of the even asphere",
          "[geom][odd_asphere][vertex]") {
  // ADR 0043, point 3: with A2 = 0 the cap is (k, A4) directly, bit for bit as EvenAsphere.
  for (const auto& [c, k] : {std::pair{1.0 / 32.0, -1.3}, std::pair{-0.02, 0.0},
                             std::pair{0.0, 0.0}, std::pair{-0.0, 0.5}}) {
    const EvenAsphere<double> even(c, k, {2e-6, -3e-9});
    const OddAsphere<double> odd(c, k, {0.0, 0.0, 0.0, 2e-6, 0.0, -3e-9});
    require_same(odd.vertex_expansion(), even.vertex_expansion());
    CHECK(odd.vertex_expansion().order == 4);
    CHECK(odd.base_conic() == even.base_conic());
    CHECK(odd.max_radius() == even.max_radius());
  }
  // Without A4 the cap has A4 = 0 as the conic.
  const OddAsphere<double> odd(0.03, -0.5, {0.0, 0.0, 0.0, 0.0, 1e-9});
  require_same(odd.vertex_expansion(), Conic<double>(0.03, -0.5).vertex_expansion());
}

TEST_CASE("odd asphere with A2: vertex curvature c + 2 A2 and the cap with native K",
          "[geom][odd_asphere][vertex]") {
  // ADR 0043, point 3: c_v = c + 2 A2, K = k, A4 = A4 + (1 + k)(c^3 - c_v^3) / 8. Dyadic values:
  // c = 1/32, A2 = 1/128, so c_v = 3/64 exactly.
  const double c = 1.0 / 32.0;
  const double k = -0.75;
  const double a2 = 1.0 / 128.0;
  const double a4 = 2e-6;
  const OddAsphere<double> odd(c, k, {0.0, a2, 0.0, a4, 0.0, 1e-9});
  const VertexExpansion<double> e = odd.vertex_expansion();
  const double cv = 3.0 / 64.0;
  CHECK(e.hessian == std::array{cv, 0.0, cv});
  CHECK(e.hessian == Conic<double>(cv, k).vertex_expansion().hessian);
  CHECK(e.value == 0.0);
  CHECK(e.slope == std::pair{0.0, 0.0});
  CHECK(e.order == 4);
  CHECK(e.rotational4);
  CHECK(e.conic == k);
  // The start value of the intersection stays the base conic (ADR 0043, point 3).
  CHECK(odd.base_conic() == std::pair{c, k});
  // Independent check of the cap: a fit z(r) = a r^2 + b r^4 + e r^6 at r = 1/4, 1/2, 1 of the
  // written-out sag gives the r^4 coefficient b4; the cap obeys K c_v^3 + 8 A4 = 8 b4 - c_v^3
  // (ADR 0041, point 1). Neglected r^8 terms are below 1e-12 of b4 at these r.
  const std::vector<double> coefficients = {0.0, a2, 0.0, a4, 0.0, 1e-9};
  const std::array<double, 3> h = {0.25, 0.5, 1.0};
  std::array<double, 3> z{};
  for (std::size_t i = 0; i < 3; ++i) z[i] = sag_written_out(c, k, coefficients, h[i]);
  // Solve [h^2 h^4 h^6] (a, b, e) = z by Cramer's rule.
  const auto det3 = [](const std::array<std::array<double, 3>, 3>& m) {
    return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
           m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
           m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
  };
  std::array<std::array<double, 3>, 3> m{};
  for (std::size_t i = 0; i < 3; ++i) {
    m[i] = {std::pow(h[i], 2.0), std::pow(h[i], 4.0), std::pow(h[i], 6.0)};
  }
  std::array<std::array<double, 3>, 3> mb = m;
  for (std::size_t i = 0; i < 3; ++i) mb[i][1] = z[i];
  const double b4 = det3(mb) / det3(m);
  const double lhs = e.conic * cv * cv * cv + 8.0 * e.a4;
  const double rhs = 8.0 * b4 - cv * cv * cv;
  CHECK(std::abs(lhs - rhs) <= 1e-6 * std::abs(rhs));
  // Content guard: A2 changes the cap away from A4 alone.
  CHECK(std::abs(e.a4 - a4) > 1e-6);
}

TEST_CASE("odd asphere order: 0 with A1, 2 with A3, else 4", "[geom][odd_asphere][vertex]") {
  // ADR 0041, point 1 and ADR 0043, point 3, decided from the parameters.
  const auto order = [](std::vector<double> a) {
    return OddAsphere<double>(0.02, -1.0, std::move(a)).vertex_expansion().order;
  };
  CHECK(order({}) == 4);
  CHECK(order({0.0, 0.0, 0.0}) == 4);
  CHECK(order({1e-3}) == 0);
  CHECK(order({1e-3, 0.0, 1e-5}) == 0);
  CHECK(order({0.0, 0.0, 1e-5}) == 2);
  CHECK(order({0.0, 1e-3, 1e-5, 1e-6}) == 2);
  CHECK(order({0.0, 0.0, 0.0, 1e-6, 1e-8}) == 4);  // |r|^5 is C^4
  CHECK(order({0.0, 1e-3}) == 4);
  for (const std::vector<double>& a :
       {std::vector<double>{1e-3}, std::vector<double>{0.0, 0.0, 1e-5}}) {
    CHECK(OddAsphere<double>(0.02, -1.0, a).vertex_expansion().rotational4);
  }
}

TEST_CASE("odd asphere without non-zero coefficients is the conic bit for bit",
          "[geom][odd_asphere]") {
  // ADR 0043, point 5.1, applied to the odd asphere: the polynomial part is only added when a
  // coefficient is non-zero, so that the sign of a zero sag survives (-0 + 0 would be +0).
  for (const double c : {-1.0 / 32.0, 0.025}) {
    const Conic<double> conic(c, -0.4);
    for (const std::vector<double>& a :
         {std::vector<double>{}, std::vector<double>{0.0, 0.0, 0.0, -0.0}}) {
      const OddAsphere<double> odd(c, -0.4, a);
      for (const auto& [x, y] : {std::pair{0.0, 0.0}, std::pair{-0.0, 0.0}, std::pair{3.0, -4.0},
                                 std::pair{1e-200, 0.0}, std::pair{-9.5, 2.25}}) {
        const double zc = conic.sag(x, y);
        const double zo = odd.sag(x, y);
        CHECK(zo == zc);
        CHECK(std::signbit(zo) == std::signbit(zc));
        const auto [gx, gy] = odd.grad(x, y);
        const auto [hx, hy] = conic.grad(x, y);
        CHECK(gx == hx);
        CHECK(gy == hy);
        CHECK(std::signbit(gx) == std::signbit(hx));
        CHECK(std::signbit(gy) == std::signbit(hy));
      }
      require_same(odd.vertex_expansion(), conic.vertex_expansion());
    }
  }
  // Content guard: the case c < 0 at the vertex has the sag -0.
  REQUIRE(std::signbit(Conic<double>(-1.0 / 32.0, -0.4).sag(0.0, 0.0)));
}

TEST_CASE("odd asphere sag and gradient against the written-out form", "[geom][odd_asphere]") {
  const double c = 1.0 / 40.0;
  const double k = -1.3;
  const std::vector<double> a = {1e-3, 2e-4, -3e-5, 1e-6, 2e-8};
  const OddAsphere<double> odd(c, k, a);
  int checked = 0;
  for (const auto& [x, y] : {std::pair{0.05, 0.0}, std::pair{3.0, -4.0}, std::pair{-0.3, 0.2},
                             std::pair{-9.0, 11.0}, std::pair{0.0, -14.5}}) {
    const double r = std::hypot(x, y);
    const double z = sag_written_out(c, k, a, r);
    CHECK(std::abs(odd.sag(x, y) - z) <= 8.0 * kEps * std::max(1.0, std::abs(z)));
    // Central differences, step 1e-6 mm: truncation ~1e-12, rounding ~eps |z| / h < 1e-8.
    const double hh = 1e-6;
    const double nx = (odd.sag(x + hh, y) - odd.sag(x - hh, y)) / (2.0 * hh);
    const double ny = (odd.sag(x, y + hh) - odd.sag(x, y - hh)) / (2.0 * hh);
    const auto [gx, gy] = odd.grad(x, y);
    CHECK(std::abs(gx - nx) <= 1e-8);
    CHECK(std::abs(gy - ny) <= 1e-8);
    ++checked;
  }
  REQUIRE(checked == 5);
}

TEST_CASE("odd asphere gradient at the vertex and at subnormal r", "[geom][odd_asphere]") {
  // ADR 0043, point 1 (F8): the A1 term is A1 x / hypot(x, y), 0 at r = 0, and does not
  // overflow for subnormal r.
  const double a1 = std::sqrt(3.0);
  const OddAsphere<double> odd(0.01, 0.0, {a1, 0.0, 1e-4});
  CHECK(odd.grad(0.0, 0.0) == std::pair{0.0, 0.0});
  CHECK(odd.sag(0.0, 0.0) == 0.0);
  const double tiny = std::numeric_limits<double>::denorm_min();
  const auto [gx, gy] = odd.grad(tiny, 0.0);
  CHECK(std::isfinite(gx));
  CHECK(gx == a1);  // the conic and A3 parts vanish next to A1 x / r = A1
  CHECK(gy == 0.0);
  const auto [hx, hy] = odd.grad(0.0, -tiny);
  CHECK(hx == 0.0);
  CHECK(hy == -a1);
  // Off the vertex the slope of the cone is A1 plus the smooth part.
  const auto [px, py] = odd.grad(0.6, 0.8);
  // g(r) = c / phi + 3 A3 r + A1 / r at r = 1 (up to rounding of hypot(0.6, 0.8)).
  CHECK(std::abs(px - 0.6 * (0.01 / std::sqrt(1.0 - 1e-4) + a1 + 3e-4)) <= 1e-14);
  CHECK(std::abs(py - 0.8 * (0.01 / std::sqrt(1.0 - 1e-4) + a1 + 3e-4)) <= 1e-14);
}

TEST_CASE("Newton on an odd asphere with A1 = tan 60 deg near the apex",
          "[geom][odd_asphere][newton]") {
  // ADR 0043, acceptance 6: rays that cross the axis 1e-12 to 1e-3 mm from the apex, tilted up
  // to 25 deg from the axis (the generatrix of the cone lies 30 deg from it), converge on the
  // surface. Fixed seed (AGENTS.md, rule 7).
  const double a1 = std::tan(60.0 * std::acos(-1.0) / 180.0);
  const OddAsphere<double> odd(0.01, 0.0, {a1, 0.0, 1e-4});
  std::mt19937_64 rng(60);
  std::uniform_real_distribution<double> unit(0.0, 1.0);
  int hits = 0;
  int max_iterations = 0;
  double max_residual = 0.0;
  for (int i = 0; i < 4000; ++i) {
    const double phi = 2.0 * std::acos(-1.0) * unit(rng);
    const double theta = unit(rng) * 25.0 * std::acos(-1.0) / 180.0;
    const Vec3 d(std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta));
    const double scale = std::pow(10.0, -12.0 + 9.0 * unit(rng));
    const Vec3 cross((2.0 * unit(rng) - 1.0) * scale, (2.0 * unit(rng) - 1.0) * scale, 0.0);
    const Vec3 o = cross - 7.0 * d;
    const auto hit = rtt::geom::intersect(odd, o, d);
    if (hit.status != HitStatus::Hit) continue;
    ++hits;
    max_iterations = std::max(max_iterations, hit.iterations);
    max_residual =
        std::max(max_residual, std::abs(hit.point.z() - odd.sag(hit.point.x(), hit.point.y())));
    // The hit lies next to the axis crossing: the cone is crossed within ~1e-3 mm of it.
    CHECK(std::abs(hit.t - 7.0) < 1e-2);
  }
  CHECK(hits == 4000);
  CHECK(max_residual <= 1e-12);
  CHECK(max_iterations <= rtt::geom::kMaxNewtonIterations);
  // An on-axis ray converges at the apex with the normal +z (gradient 0 at r = 0).
  const auto axial = rtt::geom::intersect(odd, Vec3(0.0, 0.0, -5.0), Vec3(0.0, 0.0, 1.0));
  REQUIRE(axial.status == HitStatus::Hit);
  CHECK(std::abs(axial.t - 5.0) <= 1e-12);
  CHECK(axial.normal == Vec3(0.0, 0.0, 1.0));
}
