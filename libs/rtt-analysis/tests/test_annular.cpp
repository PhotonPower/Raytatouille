// Annular pupil (#229, PR C; ADR 0032 and ADR 0033, addenda #229): cells with the exact area of
// the annulus epsilon <= |p| <= 1, and the PSF and MTF of the annular pupil against the analytic
// expected values.
//
// System: tests/reference/m2/paraboloid_stop (f = 100 mm, object at infinity, on axis) with EPD
// 10 mm and the stop aperture {5 mm, inner 2 mm}: obscuration epsilon = 2 / 5 = 0.4. The
// paraboloid images the axial point without aberration, so W = 0 and the PSF and MTF are those of
// the annular pupil in the linear model.
//
// Expected values (docs/quellen.md):
// - PSF in Strehl units, I(v) = [(2 J1(v) / v - eps^2 2 J1(eps v) / (eps v)) / (1 - eps^2)]^2,
//   v = pi r / (lambda F#_w) (Hakobyan 2018, Eq. (1); derived like the Airy pattern of ADR 0033,
//   point 6, from DLMF (10.22.1)); the first zero is the first root of J1(v) - eps J1(eps v),
//   found here by Newton on std::cyl_bessel_j with J1' = J0 - J1 / v (DLMF (10.6.2), (10.6.3)).
// - MTF(x) = [A(1, 1, d) - 2 A(1, eps, d) + A(eps, eps, d)] / (pi (1 - eps^2)), d = 2 x, x the
//   frequency over the cut-off, with the lens area A of two circles (Weisstein, MathWorld
//   "Circle-Circle Intersection", Eqs. (11), (12)): the overlap of two annuli, the same as
//   Hakobyan, Eq. (4). On the cell grid the shift by l cells is x = l / N.
// Tolerances (ADR 0033, point 10, with the addendum #229): X at N and 2N, X_2N compared. A single
// quantity (the zero) with the bound e = 2 |X_N - X_2N|; a curve (profile, MTF) with the bound
// over the curve e_max = max 2 |X_N - X_2N|, absolute here because both go down to small values;
// e and e_max at most 1e-3.

#include <algorithm>
#include <array>
#include <bit>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#include "rtt/analysis/psf.hpp"
#include "rtt/analysis/pupil.hpp"
#include "rtt/analysis/zernike.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"

using rtt::analysis::Point2;
using rtt::analysis::PupilCells;
using rtt::analysis::PupilFunction;
using rtt::analysis::PupilFunctionOptions;
using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::model::Element;
using rtt::model::System;
using rtt::trace::RayStatus;

namespace {

constexpr double kEpsilon = 0.4;
constexpr double kLambdaMm = 0.5876e-3;

CompiledSystem paraboloid(double inner_radius) {
  System s = rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m2/paraboloid_stop.rtt.json");
  s.aperture.value = rtt::model::Param(10.0);
  std::get<Element>(s.root.children[0].value).surfaces[0].aperture =
      rtt::model::CircularAperture{5.0, inner_radius};
  return rtt::compile::compile(s, rtt::material::MaterialLibrary{});
}

PupilFunction pupil(const CompiledSystem& cs, int n) {
  PupilFunctionOptions o;
  o.grid = n;
  return rtt::analysis::pupil_function(cs, PathId{0}, 0, 0, o);
}

bool same_bits(double a, double b) {
  return std::bit_cast<std::uint64_t>(a) == std::bit_cast<std::uint64_t>(b);
}

/// F#_w = 1 / (2 |M_11|) of the mapping.
double f_number(const PupilFunction& p) {
  if (!p.mapping.has_value()) {
    FAIL("no pupil mapping");
    return std::numeric_limits<double>::quiet_NaN();
  }
  return 1.0 / (2.0 * std::abs(p.mapping->m[0]));
}

/// The annular PSF in Strehl units (Hakobyan, Eq. (1)).
double annular_psf(double v, double eps) {
  if (v == 0.0) return 1.0;
  const double a = 2.0 * std::cyl_bessel_j(1.0, v) / v -
                   eps * eps * 2.0 * std::cyl_bessel_j(1.0, eps * v) / (eps * v);
  const double s = a / (1.0 - eps * eps);
  return s * s;
}

/// First positive root of J1(v) - eps J1(eps v) by Newton from the Airy zero (DLMF (10.6.2):
/// J1' = J0 - J1 / v).
double first_zero(double eps) {
  double v = 3.83;
  for (int i = 0; i < 50; ++i) {
    const double g = std::cyl_bessel_j(1.0, v) - eps * std::cyl_bessel_j(1.0, eps * v);
    const double d1 = std::cyl_bessel_j(0.0, v) - std::cyl_bessel_j(1.0, v) / v;
    const double d2 = std::cyl_bessel_j(0.0, eps * v) - std::cyl_bessel_j(1.0, eps * v) / (eps * v);
    v -= g / (d1 - eps * eps * d2);
  }
  return v;
}

/// Lens area of two circles of radii R and r with centres d apart (MathWorld, Eqs. (11), (12)).
double lens(double big, double small, double d) {
  if (d >= big + small) return 0.0;
  if (d <= big - small) return std::numbers::pi * small * small;
  const double a =
      small * small * std::acos((d * d + small * small - big * big) / (2.0 * d * small)) +
      big * big * std::acos((d * d + big * big - small * small) / (2.0 * d * big));
  return a - 0.5 * std::sqrt((-d + small + big) * (d + small - big) * (d - small + big) *
                             (d + small + big));
}

/// MTF of the annular pupil at x = nu / nu_c.
double annular_mtf(double x, double eps) {
  const double d = 2.0 * x;
  return (lens(1.0, 1.0, d) - 2.0 * lens(1.0, eps, d) + lens(eps, eps, d)) /
         (std::numbers::pi * (1.0 - eps * eps));
}

/// Intensity in Strehl units at distance r (mm) along e_x.
double along_x(const PupilFunction& p, double r) {
  const std::array<Point2, 1> x{Point2{r, 0.0}};
  return rtt::analysis::psf_at(p, x)[0];
}

/// First minimum of the PSF along e_x in [lo, hi] lambda F#_w by golden section, mm.
double first_null(const PupilFunction& p, double lo, double hi) {
  const double s = kLambdaMm * f_number(p);
  double a = lo * s;
  double b = hi * s;
  const double g = (std::sqrt(5.0) - 1.0) / 2.0;
  double c = b - g * (b - a);
  double d = a + g * (b - a);
  while (b - a > 1e-12 * s) {
    if (along_x(p, c) < along_x(p, d)) {
      b = d;
    } else {
      a = c;
    }
    c = b - g * (b - a);
    d = a + g * (b - a);
  }
  return 0.5 * (a + b);
}

}  // namespace

TEST_CASE("pupil_cells with an obscuration: the exact area of the annulus", "[annular]") {
  for (const int n : {16, 33, 64}) {
    for (const double eps : {0.25, kEpsilon}) {
      INFO("N = " << n << ", epsilon = " << eps);
      const PupilCells full = rtt::analysis::pupil_cells(n);
      const PupilCells ring = rtt::analysis::pupil_cells(n, eps);
      const double a2 = (2.0 / n) * (2.0 / n);
      double area = 0.0;
      std::size_t rim = 0;
      for (const double f : ring.fraction) area += f * a2;
      CHECK(ring.obscuration == eps);
      CHECK(std::abs(area / (std::numbers::pi * (1.0 - eps * eps)) - 1.0) <= 1e-12);
      for (int k = 0; k < n; ++k) {
        for (int i = 0; i < n; ++i) {
          const auto at = [n](int ii, int kk) {
            return static_cast<std::size_t>(kk) * static_cast<std::size_t>(n) +
                   static_cast<std::size_t>(ii);
          };
          const std::size_t j = at(i, k);
          // Mirror and swap images give the same numbers (ADR 0032, point 1); "==" is the bit
          // pattern up to the sign of zero, as in test_pupil_cells.cpp (a centre on an axis is
          // +0.0, its negation -0.0).
          const std::size_t mx = at(n - 1 - i, k);
          const std::size_t sw = at(k, i);
          CHECK(ring.fraction[j] == ring.fraction[mx]);
          CHECK(ring.px[j] == -ring.px[mx]);
          CHECK(ring.fraction[j] == ring.fraction[sw]);
          CHECK(ring.px[j] == ring.py[sw]);
          // Corners of the cell: inside the inner circle, outside it, or on its rim.
          const double x0 = -1.0 + 2.0 * i / n;
          const double y0 = -1.0 + 2.0 * k / n;
          const double x1 = x0 + 2.0 / n;
          const double y1 = y0 + 2.0 / n;
          const double far2 = std::max(x0 * x0, x1 * x1) + std::max(y0 * y0, y1 * y1);
          const double nx = x0 <= 0.0 && x1 >= 0.0 ? 0.0 : std::min(std::abs(x0), std::abs(x1));
          const double ny = y0 <= 0.0 && y1 >= 0.0 ? 0.0 : std::min(std::abs(y0), std::abs(y1));
          const double near2 = nx * nx + ny * ny;
          if (far2 < eps * eps - 1e-12) {
            CHECK(ring.fraction[j] == 0.0);  // obscured: not traced
            CHECK(ring.inner_rim[j] == 0);
          } else if (near2 > eps * eps + 1e-12) {
            // Untouched by the obscuration: bitwise the cell of the full pupil.
            CHECK(same_bits(ring.fraction[j], full.fraction[j]));
            CHECK(same_bits(ring.px[j], full.px[j]));
            CHECK(same_bits(ring.py[j], full.py[j]));
            CHECK(ring.inner_rim[j] == 0);
          } else if (ring.inner_rim[j] != 0) {
            ++rim;
            CHECK(ring.fraction[j] > 0.0);
            CHECK(ring.fraction[j] < 1.0);
          }
        }
      }
      CHECK(rim > 0);  // content guard
    }
  }
  // Obscuration 0: bitwise the cells of the full pupil, no rim.
  const PupilCells full = rtt::analysis::pupil_cells(33);
  const PupilCells zero = rtt::analysis::pupil_cells(33, 0.0);
  CHECK(zero.obscuration == 0.0);
  for (std::size_t j = 0; j < full.fraction.size(); ++j) {
    CHECK(same_bits(zero.fraction[j], full.fraction[j]));
    CHECK(same_bits(zero.px[j], full.px[j]));
    CHECK(same_bits(zero.py[j], full.py[j]));
    CHECK(zero.inner_rim[j] == 0);
  }
  for (const double bad : {-0.1, 1.0, std::numeric_limits<double>::quiet_NaN()}) {
    CHECK_THROWS_AS(rtt::analysis::pupil_cells(16, bad), std::invalid_argument);
  }
}

TEST_CASE("pupil function of the annular paraboloid: every cell of the annulus arrives",
          "[annular]") {
  const CompiledSystem cs = paraboloid(2.0);
  const PupilFunction p = pupil(cs, 32);
  CHECK(p.obscuration == kEpsilon);  // 2 mm over the stop radius of |p| = 1, 5 mm
  const PupilCells ring = rtt::analysis::pupil_cells(32, kEpsilon);
  std::size_t rim = 0;
  for (std::size_t j = 0; j < p.fraction.size(); ++j) {
    CHECK(same_bits(p.fraction[j], ring.fraction[j]));
    if (!(p.fraction[j] > 0.0)) continue;
    INFO("cell " << j);
    CHECK(p.status[j] == RayStatus::Alive);  // the inner rim too, as virtual rays
    if (ring.inner_rim[j] != 0) ++rim;
  }
  CHECK(rim > 0);  // content guard
  CHECK(p.arrived_fraction == 1.0);
  CHECK(p.warnings.empty());  // no rays.lost, no stop.clips_beam
  // The open system has no obscuration and the full cells.
  CHECK(pupil(paraboloid(0.0), 32).obscuration == 0.0);
  // A rim cell whose sample point lies inside the inner circle arrives only as a virtual ray.
  // N = 33, epsilon = 0.39375 (inner radius 1.96875 mm = 63/32 mm): the circle passes 1.9e-4
  // below the upper edge 13/33 of the centre cell of column 16, so its annular part is two
  // corner slivers whose centroid lies on the axis below the circle.
  const double eps_rim = 1.96875 / 5.0;
  const PupilCells sliver = rtt::analysis::pupil_cells(33, eps_rim);
  const PupilFunction q = pupil(paraboloid(1.96875), 33);
  CHECK(q.obscuration == eps_rim);
  std::size_t inside = 0;
  for (std::size_t j = 0; j < q.fraction.size(); ++j) {
    if (!(q.fraction[j] > 0.0)) continue;
    INFO("cell " << j << " at (" << q.px[j] << ", " << q.py[j] << ")");
    CHECK(q.status[j] == RayStatus::Alive);
    if (sliver.inner_rim[j] != 0 && std::hypot(q.px[j], q.py[j]) < eps_rim) ++inside;
  }
  CHECK(inside > 0);  // content guard: the case the virtual rim rays are for
  // Zernike: the Noll polynomials are not orthogonal on the annulus.
  const auto fit = rtt::analysis::zernike_fit(p);
  REQUIRE(fit.warnings.size() == 1);
  CHECK(fit.warnings[0].code == "zernike.pupil_incomplete");
  CHECK(fit.warnings[0].message.find("annular") != std::string::npos);
}

TEST_CASE("annular PSF: first zero, profile, energy and Strehl ratio", "[annular]") {
  const CompiledSystem cs = paraboloid(2.0);
  const PupilFunction p64 = pupil(cs, 64);
  const PupilFunction p128 = pupil(cs, 128);
  const double v0 = first_zero(kEpsilon);
  // Content guard: the obscuration moves the zero well inside the Airy zero (3.8317).
  CHECK(v0 < 3.4);
  const double scale = kLambdaMm * f_number(p128) / std::numbers::pi;  // r per unit of v
  const double expected = v0 * scale;
  const double r64 = first_null(p64, 0.9, 1.25);
  const double r128 = first_null(p128, 0.9, 1.25);
  const double e = 2.0 * std::abs(r64 - r128) / expected;
  INFO("r0: N=64 " << r64 << ", N=128 " << r128 << ", expected " << expected << ", e " << e);
  CHECK(e <= 1e-3);
  CHECK(std::abs(r128 / expected - 1.0) <= e);
  CHECK(along_x(p128, r128) < 1e-4);  // a zero, not any minimum
  // Profile, a curve: the bound over the curve (absolute, Strehl units).
  const std::array<double, 4> vs{1.0, 2.0, 3.0, 5.0};
  std::array<double, 4> i64{};
  std::array<double, 4> i128{};
  double e_max = 0.0;
  for (std::size_t k = 0; k < vs.size(); ++k) {
    i64[k] = along_x(p64, vs[k] * scale);
    i128[k] = along_x(p128, vs[k] * scale);
    e_max = std::max(e_max, 2.0 * std::abs(i64[k] - i128[k]));
  }
  e_max += 1e-12;
  CHECK(e_max <= 1e-3);
  for (std::size_t k = 0; k < vs.size(); ++k) {
    INFO("v = " << vs[k] << ": N=128 " << i128[k] << ", expected " << annular_psf(vs[k], kEpsilon)
                << ", e_max " << e_max);
    CHECK(std::abs(i128[k] - annular_psf(vs[k], kEpsilon)) <= e_max);
  }
  // Content guard: the annular profile differs from the Airy pattern at v = 3.
  const double airy3 = std::pow(2.0 * std::cyl_bessel_j(1.0, 3.0) / 3.0, 2.0);
  CHECK(std::abs(annular_psf(3.0, kEpsilon) - airy3) > 0.01);
  // Without aberration the Strehl ratio is 1 and the lattice holds 1 - edge_energy_deficit.
  rtt::analysis::PsfOptions o;
  o.pupil.grid = 64;
  const rtt::analysis::Psf psf = rtt::analysis::psf(cs, PathId{0}, 0, o);
  CHECK(std::abs(psf.strehl - 1.0) <= 1e-12);
  double sum = 0.0;
  for (const double v : psf.intensity) sum += v;
  CHECK(std::abs(sum - (1.0 - psf.edge_energy_deficit)) <= 1e-12);
  CHECK(psf.edge_energy_deficit > 0.0);
}

TEST_CASE("annular MTF against the overlap of two annuli", "[annular]") {
  // N = 128 and 256, samples x = l / 128 in [0.1, 0.9] (l = 13 ... 115) and the same x at 256.
  const CompiledSystem cs = paraboloid(2.0);
  const PupilFunction p128 = pupil(cs, 128);
  const rtt::analysis::Mtf m128 = rtt::analysis::mtf(p128);
  const rtt::analysis::Mtf m256 = rtt::analysis::mtf(pupil(cs, 256));
  REQUIRE(m128.tangential.modulation.size() == 129);
  REQUIRE(m256.tangential.modulation.size() == 257);
  // The frequency axis: sample l lies at x = l / N of the cut-off nu_c = 1 / (lambda F#_w), F#_w
  // from the mapping M that the guard of #205 ties to the geometry (review of #246, f2): a wrong
  // scale of the frequencies (n', a factor 2) would otherwise go unnoticed.
  const double nu_c = 1.0 / (kLambdaMm * f_number(p128));
  for (std::size_t l = 1; l <= 128; ++l) {
    const double x = static_cast<double>(l) / 128.0;
    CHECK(std::abs(m128.tangential.frequency[l] / nu_c - x) <= 1e-12);
    CHECK(std::abs(m128.sagittal.frequency[l] / nu_c - x) <= 1e-12);
  }
  double e_max = 0.0;
  for (std::size_t l = 13; l <= 115; ++l) {
    for (const auto* curve : {&m128.tangential, &m128.sagittal}) {
      const auto& fine = curve == &m128.tangential ? m256.tangential : m256.sagittal;
      e_max = std::max(e_max, 2.0 * std::abs(curve->modulation[l] - fine.modulation[2 * l]));
    }
  }
  e_max += 1e-12;
  CHECK(e_max <= 1e-3);
  for (std::size_t l = 13; l <= 115; ++l) {
    const double x = static_cast<double>(l) / 128.0;
    INFO("x = " << x << ": N=256 " << m256.tangential.modulation[2 * l] << ", expected "
                << annular_mtf(x, kEpsilon) << ", e_max " << e_max);
    CHECK(std::abs(m256.tangential.modulation[2 * l] - annular_mtf(x, kEpsilon)) <= e_max);
    CHECK(std::abs(m256.sagittal.modulation[2 * l] - annular_mtf(x, kEpsilon)) <= e_max);
  }
  // Content guard: at x = 0.25 the annulus transfers clearly less than the full pupil. It compares
  // the expected values only, but it is effective: a computation without the annulus would lie
  // more than 0.05 from the annular expected value, far beyond e_max.
  CHECK(annular_mtf(0.25, 0.0) - annular_mtf(0.25, kEpsilon) > 0.05);
}
