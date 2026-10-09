// Gaussian quadrature in the pupil (#168, ADR 0030 point 4): rings at the Gauss-Legendre nodes
// in u = rho^2 and evenly spaced arms. Source: NIST DLMF, Sec. 3.5(v), Eqs. (3.5.15),
// (3.5.20_1), (3.5.21), Table 3.5.1 (docs/quellen.md). The mean of f over the unit disk is
// (1/pi) int int f dA = (1/(2 pi)) int_0^{2 pi} int_0^1 f du dtheta, so R Gauss nodes in u are
// exact for polynomials in u of degree <= 2R - 1 (3.5.20_1), and A arms for e^{i m theta} with
// 0 < |m| < A (geometric sum).

#include <oneapi/tbb/task_arena.h>

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/sources.hpp"

using rtt::trace::gauss_pupil_weights;
using rtt::trace::GaussPupil;
using rtt::trace::pupil_points;
using rtt::trace::PupilPoint;

namespace {

constexpr double kEps = std::numeric_limits<double>::epsilon();

/// Mean of f(px, py) over the unit disk with the Gauss rule (weights sum to 1).
double mean(const GaussPupil& g, const std::function<double(double, double)>& f) {
  const std::vector<PupilPoint> p = pupil_points(g);
  const std::vector<double> q = gauss_pupil_weights(g);
  REQUIRE(p.size() == q.size());
  double sum = 0.0;
  for (std::size_t i = 0; i < p.size(); ++i) sum += q[i] * f(p[i].px, p[i].py);
  return sum;
}

/// RMS of f about its weighted mean.
double rms(const GaussPupil& g, const std::function<double(double, double)>& f) {
  const double m = mean(g, f);
  return std::sqrt(mean(g, [&](double x, double y) { return (f(x, y) - m) * (f(x, y) - m); }));
}

/// Angle of a pupil point from +y towards +x (the orientation of the samplings).
double theta(double px, double py) {
  return std::atan2(px, py);
}

}  // namespace

TEST_CASE("gauss pupil: nodes and weights of DLMF Table 3.5.1 (5 points)", "[sources][gauss]") {
  // One arm (theta = 0, so px = 0, py = rho): the Gauss-Legendre node is x_k = 2 rho^2 - 1 and
  // the weight w_k = 2 A q_k = 2 q_k. The table has 15 decimals (rounding <= 5e-16); recovering
  // x from rho adds at most a few ulp of 1: 2e-15.
  const GaussPupil g{5, 1};
  const std::vector<PupilPoint> p = pupil_points(g);
  const std::vector<double> q = gauss_pupil_weights(g);
  REQUIRE(p.size() == 5);
  // Table 3.5.1: +-x_k and w_k; the rings come with ascending u, so x ascending.
  const std::array<double, 5> x{-0.906179845938664, -0.538469310105683, 0.0, 0.538469310105683,
                                0.906179845938664};
  const std::array<double, 5> w{0.236926885056189, 0.478628670499366, 0.568888888888889,
                                0.478628670499366, 0.236926885056189};
  for (std::size_t k = 0; k < 5; ++k) {
    INFO("node " << k);
    CHECK(p[k].px == 0.0);
    CHECK(std::abs((2.0 * p[k].py * p[k].py - 1.0) - x[k]) <= 2e-15);
    CHECK(std::abs(2.0 * q[k] - w[k]) <= 2e-15);
  }
}

TEST_CASE("gauss pupil: polynomials in u = rho^2 up to degree 2R - 1 are exact",
          "[sources][gauss]") {
  // Mean of u^p over the disk = int_0^1 u^p du = 1 / (p + 1). At most 6 x 6 = 36 terms <= 1:
  // rounding below 64 eps. Counter-check: p = 2R is not exact (DLMF (3.5.19); the error is
  // 9e-8 or more for R <= 6, far above the rounding).
  for (int rings = 1; rings <= 6; ++rings) {
    for (const int arms : {1, 6}) {
      const GaussPupil g{rings, arms};
      INFO("rings " << rings << ", arms " << arms);
      const std::vector<double> q = gauss_pupil_weights(g);
      double total = 0.0;
      for (const double v : q) total += v;
      CHECK(std::abs(total - 1.0) <= 64.0 * kEps);
      for (int p = 0; p <= 2 * rings; ++p) {
        INFO("p " << p);
        const double m = mean(g, [p](double x, double y) { return std::pow(x * x + y * y, p); });
        const double expected = 1.0 / (p + 1);
        if (p <= 2 * rings - 1) {
          CHECK(std::abs(m - expected) <= 64.0 * kEps);
        } else {
          CHECK(std::abs(m - expected) > 1e-9);
        }
      }
    }
  }
}

TEST_CASE("gauss pupil: A arms integrate cos(m theta) exactly for 0 < m < A", "[sources][gauss]") {
  // sum_j e^{i m 2 pi j / A} = 0 for 0 < m < A and = A for m = A (geometric sum), with any radial
  // factor: mean of u cos(m theta) is 0 below m = A and mean(u) = 1/2 at m = A.
  for (const int arms : {3, 6, 8}) {
    const GaussPupil g{3, arms};
    INFO("arms " << arms);
    for (int m = 1; m <= arms; ++m) {
      INFO("m " << m);
      const double v =
          mean(g, [m](double x, double y) { return (x * x + y * y) * std::cos(m * theta(x, y)); });
      CHECK(std::abs(v - (m < arms ? 0.0 : 0.5)) <= 64.0 * kEps);
    }
  }
}

TEST_CASE("gauss pupil: RMS of the Seidel wavefronts about their mean", "[sources][gauss]") {
  // Disk means (u = rho^2, theta from +y, py = rho cos theta):
  // - defocus rho^2: <u> = 1/2, <u^2> = 1/3, sigma^2 = 1/12;
  // - spherical rho^4: <u^2> = 1/3, <u^4> = 1/5, sigma^2 = 1/5 - 1/9 = 4/45;
  // - coma rho^3 cos theta = u py: mean 0, <u^3 cos^2> = 1/4 * 1/2 = 1/8;
  // - astigmatism (rho cos theta)^2 = py^2 = u (1 + cos 2 theta) / 2: mean 1/4, <py^4> =
  //   <u^2 (1 + cos 2t)^2> / 4 = (1/3)(3/2)/4 = 1/8, sigma^2 = 1/8 - 1/16 = 1/16.
  // Exact with the rule when the radial degree of f^2 is <= 2R - 1 and its angular orders are
  // < A (astigmatism: order 4, so A >= 5). Rounding: 64 eps relative.
  const auto near = [](double a, double b) { return std::abs(a - b) <= 64.0 * kEps * b; };
  const auto u = [](double x, double y) { return x * x + y * y; };
  CHECK(near(rms({2, 6}, u), 1.0 / std::sqrt(12.0)));
  CHECK(near(rms({3, 6}, [&](double x, double y) { return u(x, y) * u(x, y); }),
             2.0 / (3.0 * std::sqrt(5.0))));
  CHECK(near(rms({2, 6}, [&](double x, double y) { return u(x, y) * y; }), 1.0 / std::sqrt(8.0)));
  CHECK(near(rms({2, 6}, [](double, double y) { return y * y; }), 0.25));
  // Counter-check: with A = 4 the order-4 term of py^4 is not cancelled (0.3227, not 0.25).
  CHECK(std::abs(rms({2, 4}, [](double, double y) { return y * y; }) - 0.25) > 0.05);
}

TEST_CASE("gauss pupil: order, defaults and input errors", "[sources][gauss]") {
  // Ring-major, arms from +y towards +x (as HexapolarPupil); rings with ascending rho.
  const GaussPupil g{2, 4};
  const std::vector<PupilPoint> p = pupil_points(g);
  REQUIRE(p.size() == 8);
  CHECK(p[0].px == 0.0);
  CHECK(p[0].py > 0.0);
  CHECK(p[1].px > 0.0);  // second arm: theta = pi/2 towards +x
  CHECK(std::abs(p[1].py) <= 1e-16);
  CHECK(std::hypot(p[4].px, p[4].py) > std::hypot(p[0].px, p[0].py));
  const std::vector<double> q = gauss_pupil_weights(g);
  CHECK(q[0] == q[1]);  // same ring, same weight
  CHECK(GaussPupil{}.rings == 3);
  CHECK(GaussPupil{}.arms == 6);
  CHECK_THROWS_AS(pupil_points(GaussPupil{0, 6}), std::invalid_argument);
  CHECK_THROWS_AS(pupil_points(GaussPupil{3, 0}), std::invalid_argument);
  CHECK_THROWS_AS(gauss_pupil_weights(GaussPupil{0, 6}), std::invalid_argument);
  CHECK_THROWS_AS(gauss_pupil_weights(GaussPupil{3, -1}), std::invalid_argument);
}

TEST_CASE("gauss pupil: make_rays aims the points, weights stay power, any thread count",
          "[sources][gauss]") {
  // Each ray equals the ray aimed at the same pupil point with SinglePupilPoint; the ray weight
  // is the power (ADR 0021), the quadrature weights are gauss_pupil_weights(); bitwise the same
  // for 1 and 4 threads (#119).
  const rtt::model::System s =
      rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m1/singlet_const.rtt.json");
  const rtt::material::MaterialLibrary lib;
  const rtt::compile::CompiledSystem cs = rtt::compile::compile(s, lib);
  const std::array<std::uint16_t, 1> fields{1};
  const GaussPupil g{3, 6};
  const auto run = [&](int threads) {
    oneapi::tbb::task_arena arena(threads);
    return arena.execute(
        [&] { return rtt::trace::make_rays(cs, rtt::compile::PathId{0}, fields, 0, g); });
  };
  const rtt::trace::RayBatch one = run(1);
  const rtt::trace::RayBatch four = run(4);
  const std::vector<PupilPoint> p = pupil_points(g);
  REQUIRE(one.size() == p.size());
  for (std::size_t i = 0; i < p.size(); ++i) {
    INFO("ray " << i);
    const rtt::trace::RayBatch single = rtt::trace::make_rays(
        cs, rtt::compile::PathId{0}, fields, 0, rtt::trace::SinglePupilPoint{p[i].px, p[i].py});
    CHECK(one.pos_x()[i] == single.pos_x()[0]);
    CHECK(one.pos_y()[i] == single.pos_y()[0]);
    CHECK(one.dir_x()[i] == single.dir_x()[0]);
    CHECK(one.dir_y()[i] == single.dir_y()[0]);
    CHECK(one.dir_z()[i] == single.dir_z()[0]);
    CHECK(one.weight()[i] == single.weight()[0]);
    CHECK(one.pos_x()[i] == four.pos_x()[i]);
    CHECK(one.dir_y()[i] == four.dir_y()[i]);
    CHECK(one.status()[i] == rtt::trace::RayStatus::Alive);
  }
}
