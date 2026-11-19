// Zernike sag terms in the analyses (#316, ADR 0051, point 9, PR 2): a defocus term changes the
// OPD as the sign rule predicts (the T1 acceptance of P3), and a fit of the OPD of a mirror with
// Zernike terms returns the surface coefficients (normalisation, angle and index agree with
// rtt::analysis::zernike()).
//
// Sign rule (ADR 0051, point 9): a surface moved by dz > 0 (local +z) changes the OPL of a ray
// at normal incidence by (n - n') dz, with n' = -n for a reflection; with W = OPL_ref - OPL
// (ADR 0032) the OPD changes by -(n - n') (dz(rho) - dz(0)). VACUUM keeps the environment index
// exactly 1; W is in waves, 1e3 mm / um.
//
// Bounds:
// - Second order of the term: it tilts the surface by at most s = max |grad dz|, the ray then by
//   alpha <= |n - n'| s. Along a ray W does not change (the wavefront is a surface of constant
//   OPL, theorem of Malus and Dupin); what changes is where the ray meets the exit pupil, by at
//   most L alpha with L the distance from the surface to the exit pupil. That moves W by at most
//   |grad W| L alpha <= L alpha^2 (in mm; |grad W| <= alpha). The exit pupil lies between the
//   stop and the ideal lens, so L <= (surface to lens) + (stop to lens).
// - The ideal lens is free of aberration only for collimated light (ADR 0031): the term makes the
//   beam slightly convergent, and the lens then adds spherical aberration of about NA^2 times W.
//   This is an empirical safety margin, not a derived bound: measured 3.9 % at NA = 0.2 with this
//   code. The focal length 5000 mm keeps NA = 0.002, so 2 NA^2 |c| is about 1e-7 waves against a
//   signal of 1e-2 waves; the result does not depend on the margin.
// The fit is a projection: with orthonormal Z_j (mean Z_j^2 = 1) Cauchy-Schwarz gives
// |dc_j| = |mean(dW Z_j)| <= max |dW| on a sampling that keeps the orthonormality, and the factor 2
// covers the discrete grid. The tests use eps = 2 * 1e3 L alpha^2 / lambda + 2 NA^2 |c| (waves);
// max_slope samples the slope on a fine polar grid (an estimate, with the factor 2 as margin).

#include <algorithm>
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

#include "rtt/analysis/nominal_reference.hpp"
#include "rtt/analysis/opd.hpp"
#include "rtt/analysis/pupil.hpp"
#include "rtt/analysis/zernike.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
#include "rtt/paraxial/paraxial.hpp"
#include "rtt/trace/sources.hpp"

using rtt::compile::CompiledSystem;
using rtt::compile::PathId;

namespace {

constexpr double kLambda = 0.5876;  // um
constexpr double kR = 10.0;         // normalisation radius = stop radius at the surface, mm

CompiledSystem compiled(const std::string& json) {
  return rtt::compile::compile(rtt::io::parse_system(json), rtt::material::MaterialLibrary{});
}

/// Collimated light through a stop (radius 10) onto a lens of two planes (CONST:1.5, 5 mm)
/// whose exit surface carries `terms`, then an ideal lens f = 5000 (NA = 0.002) and the focus.
std::string plate_system(const std::string& terms) {
  return R"({
  "schema_version": "0.11.0",
  "units": {"length": "mm", "wavelength": "um"},
  "environment": {"medium": "VACUUM"},
  "wavelengths": [{"um": 0.5876, "reference": true}],
  "aperture": {"type": "epd", "value": 20.0},
  "fields": {"points": [{}]},
  "root": {"type": "assembly", "name": "system", "children": [
    {"type": "stop", "name": "stop", "pose": {"position": [0.0, 0.0, -1.0]},
     "surfaces": [{"id": "STO", "aperture": {"type": "circular", "radius": 10.0}}]},
    {"type": "lens", "name": "L", "material": "CONST:1.5", "surfaces": [
      {"id": "A", "aperture": {"type": "circular", "radius": 12.0}},
      {"id": "B", "pose": {"position": [0.0, 0.0, 5.0], "reference": "relative_to_preceding"},
       "shape": {"base": {"type": "plane"}, "terms": [)" +
         terms + R"(]},
       "aperture": {"type": "circular", "radius": 12.0}}]},
    {"type": "thin_element", "name": "IL", "pose": {"position": [0.0, 0.0, 6.0]},
     "surfaces": [{"id": "IL", "aperture": {"type": "circular", "radius": 15.0},
                   "interaction": {"type": "ideal_lens", "focal_length": 5000.0}}]},
    {"type": "detector", "name": "F", "pose": {"position": [0.0, 0.0, 5006.0]},
     "surfaces": [{"id": "F"}]}
  ]},
  "paths": [{"name": "main", "events": "auto"}]
})";
}

/// Collimated light through a stop (radius 10) onto a plane mirror with `terms`, back through an
/// ideal lens f = 5000 at z = -20 to the focus at z = -5020.
std::string mirror_system(const std::string& terms) {
  return R"({
  "schema_version": "0.11.0",
  "units": {"length": "mm", "wavelength": "um"},
  "environment": {"medium": "VACUUM"},
  "wavelengths": [{"um": 0.5876, "reference": true}],
  "aperture": {"type": "epd", "value": 20.0},
  "fields": {"points": [{}]},
  "root": {"type": "assembly", "name": "system", "children": [
    {"type": "stop", "name": "stop", "pose": {"position": [0.0, 0.0, -1.0]},
     "surfaces": [{"id": "STO", "aperture": {"type": "circular", "radius": 10.0}}]},
    {"type": "mirror", "name": "M", "surfaces": [
      {"id": "M", "shape": {"base": {"type": "plane"}, "terms": [)" +
         terms + R"(]},
       "aperture": {"type": "circular", "radius": 12.0}}]},
    {"type": "thin_element", "name": "IL", "pose": {"position": [0.0, 0.0, -20.0]},
     "surfaces": [{"id": "IL", "aperture": {"type": "circular", "radius": 15.0},
                   "interaction": {"type": "ideal_lens", "focal_length": 5000.0}}]},
    {"type": "detector", "name": "F", "pose": {"position": [0.0, 0.0, -5020.0]},
     "surfaces": [{"id": "F"}]}
  ]},
  "paths": [{"name": "main", "events": "auto"}]
})";
}

std::string noll_term(const std::vector<double>& a) {
  std::ostringstream s;
  s.precision(17);
  s << R"({"type": "zernike_sag", "normalization_radius": 10.0, "coefficients": [)";
  for (std::size_t i = 0; i < a.size(); ++i) s << (i > 0 ? ", " : "") << a[i];
  s << "]}";
  return s.str();
}

/// Largest |grad dz| of the Noll term over the disk of radius kR, from the polynomials.
double max_slope(const std::vector<double>& a) {
  double s = 0.0;
  for (int i = 0; i <= 40; ++i) {
    for (int k = 0; k < 72; ++k) {
      const double rho = i / 40.0;
      const double th = k * 2.0 * 3.141592653589793 / 72.0;
      const double h = 1e-6;
      const auto z = [&](double r, double t) {
        double v = 0.0;
        for (std::size_t j = 0; j < a.size(); ++j) {
          v += a[j] * rtt::analysis::zernike(static_cast<int>(j) + 1, r, t);
        }
        return v;
      };
      const double x = rho * std::cos(th) * kR;
      const double y = rho * std::sin(th) * kR;
      const auto at = [&](double px, double py) {
        return z(std::hypot(px, py) / kR, std::atan2(py, px));
      };
      const double gx = (at(x + h, y) - at(x - h, y)) / (2.0 * h);
      const double gy = (at(x, y + h) - at(x, y - h)) / (2.0 * h);
      s = std::max(s, std::hypot(gx, gy));
    }
  }
  return s;
}

}  // namespace

TEST_CASE("a Zernike defocus term changes the OPD by -(n - 1) dz (ADR 0051, T1 acceptance)",
          "[zernike][opd]") {
  // a4 Z4 (Noll) on the exit surface (glass n = 1.5 -> vacuum): dW = -(n - 1)(dz - dz(0)) =
  // -(n - 1) 2 sqrt(3) a4 rho^2, so the fit gives c4 = -1e3 (n - 1) a4 / lambda.
  const double a4 = 2e-5;  // mm
  const CompiledSystem with = compiled(plate_system(noll_term({0.0, 0.0, 0.0, a4})));
  const CompiledSystem without = compiled(plate_system(noll_term({0.0})));
  rtt::analysis::PupilFunctionOptions o;
  o.grid = 32;
  rtt::analysis::ZernikeOptions zo;
  zo.terms = 15;
  const auto fit_with =
      rtt::analysis::zernike_fit(rtt::analysis::pupil_function(with, PathId{0}, 0, 0, o), zo);
  const auto fit_without =
      rtt::analysis::zernike_fit(rtt::analysis::pupil_function(without, PathId{0}, 0, 0, o), zo);
  const double expected = -1e3 * 0.5 * a4 / kLambda;
  // L = 1 mm (surface B to the ideal lens) + 7 mm (stop to the lens); alpha = (n - 1) s;
  // NA = 10 / 5000.
  const double alpha = 0.5 * max_slope({0.0, 0.0, 0.0, a4});
  const double na = 10.0 / 5000.0;
  const double eps =
      2.0 * 1e3 * 8.0 * alpha * alpha / kLambda + 2.0 * na * na * std::abs(expected) + 1e-9;
  INFO("c4 " << fit_with.coefficients[3] << ", expected " << expected << ", eps " << eps);
  CHECK(std::abs(fit_with.coefficients[3] - expected) <= eps);
  CHECK(fit_with.coefficients[3] < 0.0);  // content guard: the sign of the rule
  CHECK(std::abs(fit_without.coefficients[3]) <= eps);
  CHECK(std::abs(expected) > 50.0 * eps);  // the test separates the effect from the bound
}

TEST_CASE("the OPD fit of a mirror returns its Zernike coefficients (ADR 0051, point 9)",
          "[zernike][opd]") {
  // a_j for j = 4 ... 15 with alternating signs on a plane mirror (n = 1, n' = -1):
  // dW = -2 (dz - dz(0)), so c_j = -2e3 a_j / lambda. Astigmatism, coma and the higher terms are
  // not rotationally symmetric: the rays are aimed with the reference of the nominal mirror
  // without terms (ADR 0041, point 3); without it the trial throws ParaxialError.
  std::vector<double> a(15, 0.0);
  for (std::size_t j = 3; j < 15; ++j) {
    a[j] = (j % 2 == 0 ? 1.0 : -1.0) * (1.0 + 0.1 * static_cast<double>(j)) * 1e-6;
  }
  const CompiledSystem trial = compiled(mirror_system(noll_term(a)));
  const CompiledSystem nominal = compiled(mirror_system(noll_term({0.0})));
  CHECK_THROWS_AS(rtt::analysis::opd_points(trial, PathId{0}, 0, 0, rtt::trace::GridPupil{24}),
                  rtt::paraxial::ParaxialError);
  rtt::analysis::OpdOptions o;
  o.reference = rtt::analysis::nominal_reference(nominal, PathId{0});
  const auto opd = rtt::analysis::opd_points(trial, PathId{0}, 0, 0, rtt::trace::GridPupil{24}, o);
  // Least squares over the arrived points with Z_1 ... Z_15 (normal equations in long double).
  constexpr std::size_t kJ = 15;
  std::vector<long double> ata(kJ * kJ, 0.0L);
  std::vector<long double> atw(kJ, 0.0L);
  std::size_t used = 0;
  for (const auto& p : opd.points) {
    const double rho = std::hypot(p.px, p.py);
    if (rho > 1.0) continue;
    REQUIRE(p.status == rtt::trace::RayStatus::Alive);
    ++used;
    std::array<double, kJ> z{};
    for (std::size_t j = 0; j < kJ; ++j) {
      z[j] = rtt::analysis::zernike(static_cast<int>(j) + 1, rho, std::atan2(p.py, p.px));
    }
    for (std::size_t i = 0; i < kJ; ++i) {
      atw[i] += static_cast<long double>(z[i]) * p.w;
      for (std::size_t j = 0; j < kJ; ++j) ata[i * kJ + j] += static_cast<long double>(z[i]) * z[j];
    }
  }
  REQUIRE(used > 300);
  // Gaussian elimination (small, well conditioned on the uniform grid).
  std::vector<long double> c = atw;
  for (std::size_t k = 0; k < kJ; ++k) {
    for (std::size_t i = k + 1; i < kJ; ++i) {
      const long double f = ata[i * kJ + k] / ata[k * kJ + k];
      for (std::size_t j = k; j < kJ; ++j) ata[i * kJ + j] -= f * ata[k * kJ + j];
      c[i] -= f * c[k];
    }
  }
  for (std::size_t k = kJ; k-- > 0;) {
    for (std::size_t j = k + 1; j < kJ; ++j) c[k] -= ata[k * kJ + j] * c[j];
    c[k] /= ata[k * kJ + k];
  }
  // L = 20 mm (mirror to the ideal lens) + 21 mm (stop to the lens along the light, 1 + 20);
  // alpha = 2 s;
  // NA = 10 / 5000.
  const double alpha = 2.0 * max_slope(a);
  const double na = 10.0 / 5000.0;
  double largest = 0.0;
  for (std::size_t j = 3; j < kJ; ++j) largest = std::max(largest, 2e3 * std::abs(a[j]) / kLambda);
  const double eps = 2.0 * 1e3 * 41.0 * alpha * alpha / kLambda + 2.0 * na * na * largest + 1e-9;
  for (std::size_t j = 3; j < kJ; ++j) {
    const double expected = -2e3 * a[j] / kLambda;
    INFO("j " << j + 1 << ": c " << static_cast<double>(c[j]) << ", expected " << expected
              << ", eps " << eps);
    CHECK(std::abs(static_cast<double>(c[j]) - expected) <= eps);
    // content guard: the sign alternates with a_j and the effect exceeds the bound
    CHECK((c[j] < 0.0L) == (a[j] > 0.0));
    CHECK(std::abs(expected) > 10.0 * eps);
  }
}
