// Tests of the sources of the non-sequential tracer (ADR 0057, point 6, Nachtrag Review PR 1, d,
// Nachtrag PR 4, f).

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <cstdint>
#include <numbers>
#include <stdexcept>

#include "rtt/math/isometry.hpp"
#include "rtt/math/philox.hpp"
#include "rtt/nonseq/source.hpp"

using rtt::math::Vec3;
using namespace rtt::nonseq;
using Catch::Matchers::WithinAbs;

namespace {

constexpr double kPi = std::numbers::pi;

Source point(double max_angle, std::uint64_t n, Sampling s = Sampling::Random) {
  Source src;
  src.emitter = PointEmitter{max_angle};
  src.power = 1.0;
  src.rays = n;
  src.sampling = s;
  return src;
}

/// The two uniform numbers of draw j of start ray `index` (ADR 0057, Nachtrag PR 4, f), computed
/// here independently of source.cpp.
std::pair<double, double> xi(std::uint32_t source,
                             std::uint64_t index,
                             std::uint64_t seed,
                             std::uint32_t j) {
  const auto w = rtt::math::philox4x32_10(
      {source << 16, static_cast<std::uint32_t>(index), static_cast<std::uint32_t>(index >> 32), j},
      {static_cast<std::uint32_t>(seed), static_cast<std::uint32_t>(seed >> 32)});
  return {rtt::math::uniform53(w[0], w[1]), rtt::math::uniform53(w[2], w[3])};
}

}  // namespace

TEST_CASE("Ray counts of random and grid sampling", "[nonseq][source]") {
  CHECK(source_ray_count(point(1.0, 10)) == 10);
  // Two dimensions: n = ceil(sqrt(N)), n^2 rays (Nachtrag Review PR 1, d).
  CHECK(source_ray_count(point(1.0, 10, Sampling::Grid)) == 16);
  CHECK(source_ray_count(point(1.0, 16, Sampling::Grid)) == 16);
  CHECK(source_ray_count(point(1.0, 17, Sampling::Grid)) == 25);
  CHECK(source_ray_count(point(1.0, 1, Sampling::Grid)) == 1);
  Source col = point(1.0, 1000001, Sampling::Grid);
  col.emitter = CollimatedEmitter{Disc{1.0}};
  CHECK(source_ray_count(col) == 1002001);  // n = 1001
  // Lambert with a shape: four dimensions, n = ceil(N^(1/4)), n^4 rays.
  Source lam = point(1.0, 16, Sampling::Grid);
  lam.emitter = LambertEmitter{Disc{1.0}, kPi / 2};
  CHECK(source_ray_count(lam) == 16);
  lam.rays = 17;
  CHECK(source_ray_count(lam) == 81);
  lam.rays = 1;
  CHECK(source_ray_count(lam) == 1);
  lam.sampling = Sampling::Random;
  lam.rays = 17;
  CHECK(source_ray_count(lam) == 17);
}

TEST_CASE("Grid start rays at the cell centres", "[nonseq][source]") {
  // n = 2: index = i1 + 2 i2, xi = ((i1 + 1/2)/2, (i2 + 1/2)/2).
  SECTION("point, cone of 90 degrees: 1 - cos theta = 2 xi_1 sin^2(45 deg) = xi_1") {
    const Source s = point(kPi / 2, 4, Sampling::Grid);
    // index 0: xi = (0.25, 0.25): cos theta = 0.75, phi = pi/2.
    SourceRay r = source_ray(s, 0, 0, 0);
    const double sin0 = std::sqrt(1.0 - 0.75 * 0.75);
    CHECK_THAT(r.direction.x(), WithinAbs(0.0, 1e-15));
    CHECK_THAT(r.direction.y(), WithinAbs(sin0, 1e-15));
    CHECK_THAT(r.direction.z(), WithinAbs(0.75, 1e-15));
    CHECK(r.position.isZero(0.0));
    // index 3: xi = (0.75, 0.75): cos theta = 0.25, phi = 3 pi / 2.
    r = source_ray(s, 0, 3, 0);
    CHECK_THAT(r.direction.y(), WithinAbs(-std::sqrt(1.0 - 0.0625), 1e-15));
    CHECK_THAT(r.direction.z(), WithinAbs(0.25, 1e-15));
    // index 1: xi = (0.75, 0.25).
    r = source_ray(s, 0, 1, 0);
    CHECK_THAT(r.direction.z(), WithinAbs(0.25, 1e-15));
    CHECK(r.direction.y() > 0.0);
  }
  SECTION("full sphere: 1 - cos theta = 2 xi_1") {
    const SourceRay r = source_ray(point(kPi, 4, Sampling::Grid), 0, 1, 0);
    CHECK_THAT(r.direction.z(), WithinAbs(-0.5, 1e-15));
  }
  SECTION("collimated disc: r = R sqrt(xi_1), phi = 2 pi xi_2") {
    Source s = point(1.0, 4, Sampling::Grid);
    s.emitter = CollimatedEmitter{Disc{2.0}};
    SourceRay r = source_ray(s, 0, 0, 0);
    CHECK_THAT(r.position.x(), WithinAbs(0.0, 1e-15));
    CHECK_THAT(r.position.y(), WithinAbs(1.0, 1e-15));  // 2 sqrt(0.25)
    CHECK(r.position.z() == 0.0);
    CHECK(r.direction == Vec3::UnitZ());
    r = source_ray(s, 0, 2, 0);  // xi = (0.25, 0.75): phi = 3 pi/2
    CHECK_THAT(r.position.y(), WithinAbs(-1.0, 1e-15));
  }
  SECTION("collimated rectangle") {
    Source s = point(1.0, 4, Sampling::Grid);
    s.emitter = CollimatedEmitter{Rectangle{2.0, 3.0}};
    const SourceRay r = source_ray(s, 0, 1, 0);  // xi = (0.75, 0.25)
    CHECK_THAT(r.position.x(), WithinAbs(1.0, 1e-15));
    CHECK_THAT(r.position.y(), WithinAbs(-1.5, 1e-15));
  }
  SECTION("Lambert: index = p n^2 + q, directions on the disc of radius sin(max_angle)") {
    Source s = point(1.0, 16, Sampling::Grid);
    s.emitter = LambertEmitter{Rectangle{2.0, 2.0}, kPi / 6};
    // p = 1 (xi = (0.75, 0.25)), q = 2 (xi = (0.25, 0.75)).
    const SourceRay r = source_ray(s, 0, 6, 0);
    CHECK_THAT(r.position.x(), WithinAbs(1.0, 1e-15));
    CHECK_THAT(r.position.y(), WithinAbs(-1.0, 1e-15));
    const double rho = 0.5 * std::sqrt(0.25);  // sin 30 deg sqrt(xi_3)
    CHECK_THAT(r.direction.x(), WithinAbs(0.0, 1e-15));
    CHECK_THAT(r.direction.y(), WithinAbs(-rho, 1e-15));
    CHECK_THAT(r.direction.z(), WithinAbs(std::sqrt(1.0 - rho * rho), 1e-15));
  }
}

TEST_CASE("Random start rays follow the Philox counter of ADR 0057", "[nonseq][source]") {
  const std::uint64_t seed = 0x123456789abcdefULL;
  const std::uint64_t index = (std::uint64_t{5} << 32) + 7;
  SECTION("point") {
    const double tmax = 0.7;
    const auto [x1, x2] = xi(3, index, seed, 0);
    const double one_minus = 2.0 * x1 * std::pow(std::sin(tmax / 2), 2);
    const double c = 1.0 - one_minus;
    const double sn = std::sqrt(one_minus * (2.0 - one_minus));
    const SourceRay r = source_ray(point(tmax, 10), 3, index, seed);
    CHECK_THAT(r.direction.x(), WithinAbs(sn * std::cos(2 * kPi * x2), 1e-15));
    CHECK_THAT(r.direction.y(), WithinAbs(sn * std::sin(2 * kPi * x2), 1e-15));
    CHECK_THAT(r.direction.z(), WithinAbs(c, 1e-15));
  }
  SECTION("Lambert takes xi_3, xi_4 from draw 1") {
    Source s = point(1.0, 10);
    s.emitter = LambertEmitter{Disc{3.0}, kPi / 2};
    const auto [x1, x2] = xi(0, index, seed, 0);
    const auto [x3, x4] = xi(0, index, seed, 1);
    const SourceRay r = source_ray(s, 0, index, seed);
    CHECK_THAT(r.position.x(), WithinAbs(3.0 * std::sqrt(x1) * std::cos(2 * kPi * x2), 1e-14));
    CHECK_THAT(r.position.y(), WithinAbs(3.0 * std::sqrt(x1) * std::sin(2 * kPi * x2), 1e-14));
    const double rho = std::sqrt(x3);
    CHECK_THAT(r.direction.x(), WithinAbs(rho * std::cos(2 * kPi * x4), 1e-15));
    CHECK_THAT(r.direction.y(), WithinAbs(rho * std::sin(2 * kPi * x4), 1e-15));
    CHECK_THAT(r.direction.z(), WithinAbs(std::sqrt(1.0 - rho * rho), 1e-15));
  }
  SECTION("another seed or source gives another ray") {
    const SourceRay a = source_ray(point(kPi, 10), 0, 1, 1);
    CHECK(a.direction != source_ray(point(kPi, 10), 0, 1, 2).direction);
    CHECK(a.direction != source_ray(point(kPi, 10), 1, 1, 1).direction);
    CHECK(a.direction == source_ray(point(kPi, 10), 0, 1, 1).direction);
  }
}

TEST_CASE("Random start rays have the moments of their distributions", "[nonseq][source]") {
  // Means of 20000 samples within 4 sigma, sigma = sqrt(Var / N).
  const std::uint64_t n = 20000;
  const auto mean = [&](const Source& s, auto f) {
    double m = 0.0;
    double m2 = 0.0;
    for (std::uint64_t i = 0; i < n; ++i) {
      const double v = f(source_ray(s, 0, i, 11));
      m += v;
      m2 += v * v;
    }
    m /= static_cast<double>(n);
    const double var = m2 / static_cast<double>(n) - m * m;
    return std::pair{m, 4.0 * std::sqrt(var / static_cast<double>(n))};
  };
  SECTION("cone: cos theta uniform in [cos theta_max, 1] (PBRT 4ed, §A.5.4)") {
    const double tmax = 1.2;
    const auto [m, tol] = mean(point(tmax, n), [](const SourceRay& r) { return r.direction.z(); });
    CHECK_THAT(m, WithinAbs((1.0 + std::cos(tmax)) / 2.0, tol));
    const auto [mx, tolx] =
        mean(point(tmax, n), [](const SourceRay& r) { return r.direction.x(); });
    CHECK_THAT(mx, WithinAbs(0.0, tolx));
  }
  SECTION("disc: r^2 uniform in [0, R^2] (§A.5.1)") {
    Source s = point(1.0, n);
    s.emitter = CollimatedEmitter{Disc{2.0}};
    const auto [m, tol] = mean(s, [](const SourceRay& r) { return r.position.squaredNorm(); });
    CHECK_THAT(m, WithinAbs(2.0, tol));
  }
  SECTION("Lambert: sin^2 theta uniform in [0, sin^2 theta_max] (Malley, §A.5.3)") {
    Source s = point(1.0, n);
    const double tmax = 0.9;
    s.emitter = LambertEmitter{Rectangle{1.0, 1.0}, tmax};
    const auto [m, tol] =
        mean(s, [](const SourceRay& r) { return 1.0 - r.direction.z() * r.direction.z(); });
    CHECK_THAT(m, WithinAbs(std::pow(std::sin(tmax), 2) / 2.0, tol));
    const auto [mx, tolx] = mean(s, [](const SourceRay& r) { return r.position.x(); });
    CHECK_THAT(mx, WithinAbs(0.0, tolx));
  }
}

TEST_CASE("The pose maps start rays to global coordinates", "[nonseq][source]") {
  Source s = point(kPi / 2, 4, Sampling::Grid);
  s.emitter = CollimatedEmitter{Disc{2.0}};
  // Rotation by 90 degrees about x: local +z -> global -y, local +y -> global +z.
  const rtt::math::Isometry3 pose((rtt::math::Mat3() << 1, 0, 0, 0, 0, -1, 0, 1, 0).finished(),
                                  Vec3(1.0, 2.0, 3.0));
  s.pose = pose;
  const SourceRay r = source_ray(s, 0, 0, 0);  // local position (0, 1, 0)
  CHECK_THAT(r.position.x(), WithinAbs(1.0, 1e-15));
  CHECK_THAT(r.position.y(), WithinAbs(2.0, 1e-15));
  CHECK_THAT(r.position.z(), WithinAbs(4.0, 1e-15));
  CHECK_THAT(r.direction.y(), WithinAbs(-1.0, 1e-15));
}

TEST_CASE("check_source rejects invalid sources", "[nonseq][source]") {
  CHECK_NOTHROW(check_source(point(kPi, 1), 2));
  CHECK_THROWS_AS(check_source(point(0.0, 1), 2), std::invalid_argument);
  CHECK_THROWS_AS(check_source(point(kPi + 1e-9, 1), 2), std::invalid_argument);
  CHECK_THROWS_AS(check_source(point(1.0, 0), 2), std::invalid_argument);
  CHECK_THROWS_AS(check_source(point(1.0, (std::uint64_t{1} << 40) + 1), 2), std::invalid_argument);
  Source s = point(1.0, 1);
  s.power = 0.0;
  CHECK_THROWS_AS(check_source(s, 2), std::invalid_argument);
  s.power = std::nan("");
  CHECK_THROWS_AS(check_source(s, 2), std::invalid_argument);
  s.power = 1.0;
  s.wavelengths = {0, 2};
  CHECK_THROWS_AS(check_source(s, 2), std::invalid_argument);
  s.wavelengths = {1, 1};
  CHECK_THROWS_AS(check_source(s, 2), std::invalid_argument);
  s.wavelengths = {1, 0};
  CHECK_NOTHROW(check_source(s, 2));
  s.emitter = LambertEmitter{Disc{1.0}, kPi / 2 + 1e-9};
  CHECK_THROWS_AS(check_source(s, 2), std::invalid_argument);
  s.emitter = LambertEmitter{Disc{0.0}, kPi / 2};
  CHECK_THROWS_AS(check_source(s, 2), std::invalid_argument);
  s.emitter = CollimatedEmitter{Rectangle{1.0, -1.0}};
  CHECK_THROWS_AS(check_source(s, 2), std::invalid_argument);
  s.emitter = CollimatedEmitter{Rectangle{1.0, 1.0}};
  s.pose = rtt::math::Isometry3(rtt::math::Mat3::Identity(), Vec3(std::nan(""), 0.0, 0.0));
  CHECK_THROWS_AS(check_source(s, 2), std::invalid_argument);
}
