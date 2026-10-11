#pragma once

/// @file pupil.hpp
/// Pupil function and its cells (ADR 0032; #202).
///
/// Conventions:
/// - Normalised pupil coordinates p = (px, py) as in rtt-trace (sources.hpp): the unit circle is
///   the rim of the paraxial entrance pupil, +y meridional.
/// - The cells divide the square [-1, 1]^2 into N x N squares of edge a = 2/N, in the order of
///   trace::GridPupil: rows from py = -1 upwards, px fastest; cell (i, k), i along x and k along
///   y, has index k N + i and the centre (-1 + (i + 1/2) a, -1 + (k + 1/2) a).

#include <vector>

namespace rtt::analysis {

/// Largest grid N of the pupil function. Memory bound: each traced cell is a ray of a
/// trace::RayBatch (13 doubles, 9 complex PRT elements, ids and status: about 260 bytes) plus up
/// to about 390 bytes in the pupil function (cells, end states, PRT copy), so N = 2048 needs
/// about 0.79 N^2 = 3.3 million traced cells (pi/4 of N^2) or 2.7 GB at most. A typo such as
/// 20480 would ask for 100 times that; it is an input error instead. The acceptance of M6 needs
/// N <= 256.
inline constexpr int kMaxPupilGrid = 2048;

/// Cells of the pupil function (ADR 0032, point 1), N x N entries each in the order above.
struct PupilCells {
  int grid = 0;  ///< N
  /// Area fraction f in [0, 1] of the cell inside the unit circle, in closed form: exactly 1 for
  /// a cell inside the circle, exactly 0 for a cell outside.
  std::vector<double> fraction;
  /// Sample point: the centre of a full cell (f = 1) and of an empty one (f = 0, not traced);
  /// the centroid of the part inside the circle for a boundary cell (0 < f < 1), which lies
  /// strictly inside the circle.
  std::vector<double> px;
  std::vector<double> py;
};

/// Cells of an N x N grid, N = `grid`.
/// @throws std::invalid_argument for grid < 2 or grid > kMaxPupilGrid
[[nodiscard]] PupilCells pupil_cells(int grid);

}  // namespace rtt::analysis
