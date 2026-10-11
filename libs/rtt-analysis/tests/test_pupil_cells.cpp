// Cells of the pupil function (ADR 0032, point 1; #202): N x N cells of edge a = 2/N on
// [-1, 1]^2, area fraction f inside the unit circle in closed form, sample point = centre of a
// full cell or centroid of the part inside.
//
// Tolerances:
// - sum f a^2 = pi: each f is a sum of a few closed-form terms, each rounded to a few u
//   (u = 2^-53 = 1.1e-16); N^2 cells add up to N^2 such errors of size about u a^2, so the sum is
//   off by about N^2 u a^2 = 4 u. 1e-13 relative is more than 100 times that (ADR 0032).
// - Fine midpoint rule over l(t) = max(0, min(y1, h) - max(y0, -h)), h = sqrt(1 - t^2), with M
//   steps of width s = a / M: l is continuous and piecewise smooth, kinks give O(s^2) per kink
//   and the square root at t = +-1 gives O(s^1.5) (the integrand behaves like sqrt(1 - |t|)). With
//   M = 200000 and a <= 1: s^1.5 <= 1.2e-8, relative to a^2 below 1e-6 for N <= 64
//   (a = 1/32: s^1.5 = 4.4e-11, a^2 = 9.8e-4, ratio 4.5e-8). Tolerance on f: 1e-6. Centroids
//   divide by the area: an error e of a moment and of the area moves the centroid by at most
//   about 2 e / (f a^2); for f >= 1e-2 and e <= 1e-6 f a^2-scale this is below 1e-4, so the
//   centroid tolerance is 1e-4 for cells with f >= 1e-2.

#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <stdexcept>
#include <vector>

#include "rtt/analysis/pupil.hpp"

using rtt::analysis::kMaxPupilGrid;
using rtt::analysis::pupil_cells;
using rtt::analysis::PupilCells;

namespace {

/// Index of cell (i, k): i along x, k along y; rows from py = -1, px fastest.
std::size_t cell(int n, int i, int k) {
  return static_cast<std::size_t>(k) * static_cast<std::size_t>(n) + static_cast<std::size_t>(i);
}

struct Moments {
  double area = 0.0;
  double mx = 0.0;
  double my = 0.0;
};

/// Area and first moments of [x0, x1] x [y0, y1] inside the unit circle by the midpoint rule
/// over x with `steps` steps (independent of the closed form of the code).
Moments midpoint(double x0, double x1, double y0, double y1, int steps) {
  Moments m;
  const double s = (x1 - x0) / steps;
  for (int j = 0; j < steps; ++j) {
    const double t = x0 + (j + 0.5) * s;
    if (std::abs(t) >= 1.0) continue;
    const double h = std::sqrt(1.0 - t * t);
    const double top = std::min(y1, h);
    const double bottom = std::max(y0, -h);
    if (top <= bottom) continue;
    m.area += (top - bottom) * s;
    m.mx += t * (top - bottom) * s;
    m.my += 0.5 * (top * top - bottom * bottom) * s;
  }
  return m;
}

}  // namespace

TEST_CASE("pupil cells: the area fractions sum to pi / a^2", "[pupil][cells]") {
  for (const int n : {2, 3, 7, 64, 65}) {
    INFO("N = " << n);
    const PupilCells c = pupil_cells(n);
    REQUIRE(c.grid == n);
    const auto size = static_cast<std::size_t>(n) * static_cast<std::size_t>(n);
    REQUIRE(c.fraction.size() == size);
    REQUIRE(c.px.size() == size);
    REQUIRE(c.py.size() == size);
    const double a = 2.0 / n;
    double sum = 0.0;
    for (const double f : c.fraction) sum += f;
    CHECK(std::abs(sum * a * a - std::numbers::pi) <= 1e-13 * std::numbers::pi);
  }
}

TEST_CASE("pupil cells: order, full cells, empty cells and sample points", "[pupil][cells]") {
  for (const int n : {2, 3, 7, 64, 65}) {
    INFO("N = " << n);
    const PupilCells c = pupil_cells(n);
    int partial = 0;
    for (int k = 0; k < n; ++k) {
      for (int i = 0; i < n; ++i) {
        INFO("cell (" << i << ", " << k << ")");
        const std::size_t j = cell(n, i, k);
        const double f = c.fraction[j];
        REQUIRE(f >= 0.0);
        REQUIRE(f <= 1.0);
        const double cx = static_cast<double>(2 * i + 1 - n) / n;
        const double cy = static_cast<double>(2 * k + 1 - n) / n;
        // A cell is full if its corner farthest from the origin lies in the circle (the disk
        // is convex), empty if its point nearest to the origin lies outside.
        const double a = 2.0 / n;
        const double far_x = std::abs(cx) + a / 2.0;
        const double far_y = std::abs(cy) + a / 2.0;
        const double near_x = std::max(0.0, std::abs(cx) - a / 2.0);
        const double near_y = std::max(0.0, std::abs(cy) - a / 2.0);
        if (far_x * far_x + far_y * far_y <= 1.0) {
          CHECK(f == 1.0);
          CHECK(c.px[j] == cx);  // the centre, exactly
          CHECK(c.py[j] == cy);
        } else if (near_x * near_x + near_y * near_y >= 1.0) {
          CHECK(f == 0.0);
          CHECK(c.px[j] == cx);  // documented: the centre, not traced
          CHECK(c.py[j] == cy);
        } else {
          ++partial;
          CHECK(f > 0.0);
          CHECK(f < 1.0);
          // The centroid of the convex part inside lies strictly inside the circle and inside
          // the cell.
          CHECK(c.px[j] * c.px[j] + c.py[j] * c.py[j] < 1.0);
          CHECK(std::abs(c.px[j] - cx) <= a / 2.0);
          CHECK(std::abs(c.py[j] - cy) <= a / 2.0);
        }
      }
    }
    CHECK(partial > 0);  // content guard: there are boundary cells
  }
}

TEST_CASE("pupil cells: symmetric under x -> -x, y -> -y and x <-> y", "[pupil][cells]") {
  // The code computes mirrored cells with mirrored inputs and integrates along the coordinate
  // nearer to the axis, so the mirror and swap images agree bit for bit; 1e-15 leaves room for
  // the diagonal cells, which are their own swap images (px against py of the same cell).
  for (const int n : {2, 3, 7, 64, 65}) {
    INFO("N = " << n);
    const PupilCells c = pupil_cells(n);
    for (int k = 0; k < n; ++k) {
      for (int i = 0; i < n; ++i) {
        INFO("cell (" << i << ", " << k << ")");
        const std::size_t j = cell(n, i, k);
        const std::size_t mx = cell(n, n - 1 - i, k);
        const std::size_t my = cell(n, i, n - 1 - k);
        const std::size_t sw = cell(n, k, i);  // x <-> y
        CHECK(std::abs(c.fraction[j] - c.fraction[mx]) <= 1e-15);
        CHECK(std::abs(c.fraction[j] - c.fraction[my]) <= 1e-15);
        CHECK(std::abs(c.fraction[j] - c.fraction[sw]) <= 1e-15);
        CHECK(std::abs(c.px[j] + c.px[mx]) <= 1e-15);
        CHECK(std::abs(c.py[j] - c.py[mx]) <= 1e-15);
        CHECK(std::abs(c.py[j] + c.py[my]) <= 1e-15);
        CHECK(std::abs(c.px[j] - c.px[my]) <= 1e-15);
        CHECK(std::abs(c.px[j] - c.py[sw]) <= 1e-15);
        CHECK(std::abs(c.py[j] - c.px[sw]) <= 1e-15);
      }
    }
  }
}

TEST_CASE("pupil cells: area fractions and centroids agree with a fine midpoint rule",
          "[pupil][cells]") {
  constexpr int kSteps = 200000;
  for (const int n : {3, 7, 64}) {
    INFO("N = " << n);
    const PupilCells c = pupil_cells(n);
    const double a = 2.0 / n;
    int checked = 0;
    for (int k = 0; k < n; ++k) {
      for (int i = 0; i < n; ++i) {
        const std::size_t j = cell(n, i, k);
        const double f = c.fraction[j];
        if (f == 0.0 || f == 1.0) continue;
        // N = 64 has about 250 boundary cells; every 7th keeps the test fast.
        if (n == 64 && j % 7 != 0) continue;
        INFO("cell (" << i << ", " << k << "), f = " << f);
        const double x0 = -1.0 + i * a;
        const double y0 = -1.0 + k * a;
        const Moments m = midpoint(x0, x0 + a, y0, y0 + a, kSteps);
        CHECK(std::abs(f - m.area / (a * a)) <= 1e-6);
        if (f >= 1e-2) {
          CHECK(std::abs(c.px[j] - m.mx / m.area) <= 1e-4);
          CHECK(std::abs(c.py[j] - m.my / m.area) <= 1e-4);
        }
        ++checked;
      }
    }
    CHECK(checked > 3);
  }
}

TEST_CASE("pupil cells: N = 3, the edge y = 1 touches the circle", "[pupil][cells]") {
  // The middle cell of the top row, [-1/3, 1/3] x [1/3, 1]: the circle lies below y = 1 except
  // at t = 0, where it touches (ADR 0032, point 1: with "<" instead of "<=" the form choice
  // would be wrong there, sum f a^2 - pi = 1.3e-2). Its area is
  // int_{-1/3}^{1/3} (h - 1/3) dt = 2 P(1/3) - 2/9 with P(t) = (t sqrt(1 - t^2) + asin t) / 2.
  const PupilCells c = pupil_cells(3);
  const double t = 1.0 / 3.0;
  const double p = 0.5 * (t * std::sqrt(1.0 - t * t) + std::asin(t));
  const double area = 2.0 * p - 2.0 / 9.0;
  const double a = 2.0 / 3.0;
  CHECK(std::abs(c.fraction[cell(3, 1, 2)] - area / (a * a)) <= 1e-15);
  CHECK(std::abs(c.fraction[cell(3, 1, 0)] - area / (a * a)) <= 1e-15);
  CHECK(c.fraction[cell(3, 1, 1)] == 1.0);  // the centre cell lies inside
}

TEST_CASE("pupil cells: grid out of range is invalid_argument", "[pupil][cells]") {
  CHECK_THROWS_AS(pupil_cells(1), std::invalid_argument);
  CHECK_THROWS_AS(pupil_cells(0), std::invalid_argument);
  CHECK_THROWS_AS(pupil_cells(-64), std::invalid_argument);
  CHECK_THROWS_AS(pupil_cells(kMaxPupilGrid + 1), std::invalid_argument);
  CHECK(kMaxPupilGrid >= 256);  // the acceptance of M6 needs N up to 256
}
