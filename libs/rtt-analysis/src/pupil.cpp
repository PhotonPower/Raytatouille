#include "rtt/analysis/pupil.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace rtt::analysis {

namespace {

/// Antiderivative of h(t) = sqrt(1 - t^2) on [-1, 1]: P(t) = (t sqrt(1 - t^2) + asin t) / 2.
double area_antiderivative(double t) {
  return 0.5 * (t * std::sqrt(1.0 - t * t) + std::asin(t));
}

/// Antiderivative of t h(t): -(1 - t^2)^(3/2) / 3.
double moment_antiderivative(double t) {
  const double r = 1.0 - t * t;
  return -r * std::sqrt(r) / 3.0;
}

/// Antiderivative of h(t)^2 = 1 - t^2: t - t^3 / 3.
double square_antiderivative(double t) {
  return t - t * t * t / 3.0;
}

/// Area and first moments of the part of [x0, x1] x [y0, y1] inside the unit circle.
struct Part {
  double area = 0.0;
  double mx = 0.0;  ///< integral of x over the part
  double my = 0.0;  ///< integral of y over the part
};

/// Closed form of ADR 0032, point 1 (derived there): for t in [x0, x1] with |t| < 1 the cut in y
/// has length l(t) = max(0, min(y1, h) - max(y0, -h)), h = sqrt(1 - t^2). The kinks
/// t = +-sqrt(1 - y0^2), +-sqrt(1 - y1^2) split [x0, x1] into pieces on which the upper end is
/// either y1 or h and the lower end either y0 or -h; the form is decided at the midpoint m of the
/// piece: upper h if h(m) <= y1, lower -h if -h(m) >= y0. The "<=" matters where an edge y = +-1
/// touches the circle in t = 0 (odd N): there h <= y1 on the whole piece with equality in one
/// point only. A piece with upper <= lower at m lies outside (its sign cannot change inside the
/// piece, the crossings are kinks). Area = int (upper - lower) dt, x moment = int t (upper -
/// lower) dt, y moment = int (upper^2 - lower^2) / 2 dt, each with the antiderivatives above.
Part part_inside(double x0, double x1, double y0, double y1) {
  Part part;
  const double lo = std::max(x0, -1.0);
  const double hi = std::min(x1, 1.0);
  if (!(lo < hi)) return part;
  // Ends and at most four kinks (a vector: std::sort on a partly filled std::array trips
  // -Warray-bounds of GCC 16 at -O2).
  std::vector<double> t;
  t.reserve(6);
  t.push_back(lo);
  for (const double y : {y0, y1}) {
    if (std::abs(y) >= 1.0) continue;
    const double kink = std::sqrt(1.0 - y * y);
    for (const double k : {-kink, kink}) {
      if (k > lo && k < hi) t.push_back(k);
    }
  }
  t.push_back(hi);
  std::sort(t.begin(), t.end());
  for (std::size_t j = 0; j + 1 < t.size(); ++j) {
    const double ta = t[j];
    const double tb = t[j + 1];
    if (!(ta < tb)) continue;
    const double m = 0.5 * (ta + tb);
    const double hm = std::sqrt(1.0 - m * m);
    const bool upper_h = hm <= y1;
    const bool lower_h = -hm >= y0;
    const double upper_m = upper_h ? hm : y1;
    const double lower_m = lower_h ? -hm : y0;
    if (!(upper_m > lower_m)) continue;
    const double dt = tb - ta;
    const double dh = area_antiderivative(tb) - area_antiderivative(ta);
    const double dth = moment_antiderivative(tb) - moment_antiderivative(ta);
    const double dt2 = 0.5 * (tb * tb - ta * ta);
    const double dh2 = square_antiderivative(tb) - square_antiderivative(ta);
    const double upper = upper_h ? dh : y1 * dt;
    const double lower = lower_h ? -dh : y0 * dt;
    const double upper_t = upper_h ? dth : y1 * dt2;
    const double lower_t = lower_h ? -dth : y0 * dt2;
    const double upper_sq = upper_h ? dh2 : y1 * y1 * dt;
    const double lower_sq = lower_h ? dh2 : y0 * y0 * dt;
    part.area += upper - lower;
    part.mx += upper_t - lower_t;
    part.my += 0.5 * (upper_sq - lower_sq);
  }
  return part;
}

}  // namespace

PupilCells pupil_cells(int grid) {
  if (grid < 2 || grid > kMaxPupilGrid) {
    throw std::invalid_argument("analysis: pupil grid must lie in [2, " +
                                std::to_string(kMaxPupilGrid) + "], got " + std::to_string(grid));
  }
  const int n = grid;
  const auto size = static_cast<std::size_t>(n) * static_cast<std::size_t>(n);
  PupilCells cells;
  cells.grid = n;
  cells.fraction.resize(size);
  cells.px.resize(size);
  cells.py.resize(size);
  const double nd = static_cast<double>(n);
  const double a = 2.0 / nd;
  // Edges and centres as (2 i - N) / N and (2 i + 1 - N) / N: the mirror cell has the negated
  // numerator, so the values are mirror images bit for bit.
  const auto edge = [nd, n](int i) { return static_cast<double>(2 * i - n) / nd; };
  for (int k = 0; k < n; ++k) {
    for (int i = 0; i < n; ++i) {
      const std::size_t j =
          static_cast<std::size_t>(k) * static_cast<std::size_t>(n) + static_cast<std::size_t>(i);
      const double cx = static_cast<double>(2 * i + 1 - n) / nd;
      const double cy = static_cast<double>(2 * k + 1 - n) / nd;
      cells.px[j] = cx;
      cells.py[j] = cy;
      double x0 = edge(i);
      double x1 = edge(i + 1);
      double y0 = edge(k);
      double y1 = edge(k + 1);
      // The disk is convex: full if the farthest corner lies inside, empty if the nearest point
      // of the cell lies outside.
      const double far_x = std::max(std::abs(x0), std::abs(x1));
      const double far_y = std::max(std::abs(y0), std::abs(y1));
      const double near_x = x0 <= 0.0 && x1 >= 0.0 ? 0.0 : std::min(std::abs(x0), std::abs(x1));
      const double near_y = y0 <= 0.0 && y1 >= 0.0 ? 0.0 : std::min(std::abs(y0), std::abs(y1));
      if (far_x * far_x + far_y * far_y <= 1.0) {
        cells.fraction[j] = 1.0;
        continue;
      }
      if (near_x * near_x + near_y * near_y >= 1.0) {
        cells.fraction[j] = 0.0;
        continue;
      }
      // Boundary cell. Cells on the negative side are computed mirrored, so that mirror cells
      // give the same numbers (x -> -x, y -> -y).
      const double sx = x1 <= 0.0 ? -1.0 : 1.0;
      const double sy = y1 <= 0.0 ? -1.0 : 1.0;
      if (sx < 0.0) {
        x0 = -edge(i + 1);
        x1 = -edge(i);
      }
      if (sy < 0.0) {
        y0 = -edge(k + 1);
        y1 = -edge(k);
      }
      // Integrated along the coordinate nearer to the axis: the antiderivatives are then evaluated
      // away from t = +-1, where asin amplifies the rounding of t by 1 / sqrt(1 - t^2) (up to
      // 5e-13 in f at N = 64 the other way round). The swapped cell (x <-> y) takes the same
      // path with the roles exchanged, so f and the centroid are symmetric bit for bit.
      // A cell on a diagonal (the same x and y range after mirroring) is its own swap image, so
      // its two moments are equal; their mean puts the centroid on the diagonal bit for bit (the
      // two formulas alone differ by up to 3e-14 in the centroid at N = 64).
      Part part;
      if (x0 == y0 && x1 == y1) {
        part = part_inside(x0, x1, y0, y1);
        part.mx = part.my = 0.5 * (part.mx + part.my);
      } else if (std::max(std::abs(x0), std::abs(x1)) > std::max(std::abs(y0), std::abs(y1))) {
        const Part swapped = part_inside(y0, y1, x0, x1);
        part = {swapped.area, swapped.my, swapped.mx};
      } else {
        part = part_inside(x0, x1, y0, y1);
      }
      if (!(part.area > 0.0)) {
        cells.fraction[j] = 0.0;  // touches the circle in a point only
        continue;
      }
      cells.fraction[j] = std::min(1.0, part.area / (a * a));
      cells.px[j] = sx * (part.mx / part.area);
      cells.py[j] = sy * (part.my / part.area);
    }
  }
  return cells;
}

}  // namespace rtt::analysis
