// FFT diffraction PSF (#205; ADR 0033, points 2-4, 7, 8, 10, 13). Reference system: the
// paraboloid with stop (tests/reference/m2/paraboloid_stop, f = 100 mm, object at infinity, on
// axis) with the entrance pupil diameter set in the test: F/10 (D = 10 mm) and F/50 (D = 2 mm).
// The paraboloid images the axial point without aberration, so W = 0 and the PSF is the Airy
// pattern of the linear model.
//
// Expected values (derived, ADR 0033 and docs/quellen.md):
// - Airy: intensity [2 J1(v)/v]^2 with v = pi r / (lambda F#_w), first zero r0 = j11 lambda F#_w
//   / pi; j11 found here by Newton on std::cyl_bessel_j (DLMF (10.6.2) for J1').
// - Defocus: S = [sin(pi w)/(pi w)]^2 for W = w rho^2 (Wyant & Creath (64)).
// - F#_w = 1 / (2 |M_11|) from the mapping, which the guard of ADR 0033, point 10 ties to the
//   geometry (1 / (2 F#) with F# = f / D) and to first_order (n' |u'|).
// Tolerances of ADR 0033, point 10: X at N and 2N, compare X_2N, bound e = 2 |X_N - X_2N|, and e
// itself must stay below the architecture's 1e-3.

#include <oneapi/tbb/task_arena.h>

#include <algorithm>
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstring>
#include <numbers>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include "rtt/analysis/psf.hpp"
#include "rtt/analysis/pupil.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
#include "rtt/paraxial/paraxial.hpp"

using rtt::analysis::Point2;
using rtt::analysis::Psf;
using rtt::analysis::PsfOptions;
using rtt::analysis::PupilFunction;
using rtt::analysis::PupilFunctionOptions;
using rtt::compile::CompiledSystem;
using rtt::compile::PathId;

namespace {

constexpr double kLambdaMm = 0.5876e-3;
constexpr double kF = 100.0;  // focal length of the paraboloid, R / 2

CompiledSystem paraboloid(double epd) {
  rtt::model::System s =
      rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m2/paraboloid_stop.rtt.json");
  s.aperture.value = rtt::model::Param(epd);
  return rtt::compile::compile(s, rtt::material::MaterialLibrary{});
}

PupilFunction pupil(const CompiledSystem& cs, int n, double focus_shift = 0.0) {
  PupilFunctionOptions o;
  o.grid = n;
  o.focus_shift = focus_shift;
  return rtt::analysis::pupil_function(cs, PathId{0}, 0, 0, o);
}

/// F#_w = 1 / (2 |M_11|) of the mapping.
double f_number(const PupilFunction& p) {
  REQUIRE(p.mapping.has_value());
  return 1.0 / (2.0 * std::abs(p.mapping->m[0]));
}

/// First positive zero of J1 by Newton, J1' = J0 - J1 / x (DLMF (10.6.2), (10.6.3)).
double j11() {
  double x = 3.83;
  for (int i = 0; i < 50; ++i) {
    const double j1 = std::cyl_bessel_j(1.0, x);
    const double d = std::cyl_bessel_j(0.0, x) - j1 / x;
    x -= j1 / d;
  }
  return x;
}

double airy(double v) {
  if (v == 0.0) return 1.0;
  const double a = 2.0 * std::cyl_bessel_j(1.0, v) / v;
  return a * a;
}

/// Intensity in Strehl units at distance r (mm) along e_x.
double on_axis_x(const PupilFunction& p, double r) {
  const std::array<Point2, 1> x{Point2{r, 0.0}};
  return rtt::analysis::psf_at(p, x)[0];
}

/// First minimum of the PSF along e_x in [1.0, 1.45] lambda F#_w by golden section, mm.
double first_null(const PupilFunction& p) {
  const double s = kLambdaMm * f_number(p);
  double a = 1.0 * s;
  double b = 1.45 * s;
  const double g = (std::sqrt(5.0) - 1.0) / 2.0;
  double c = b - g * (b - a);
  double d = a + g * (b - a);
  while (b - a > 1e-12 * s) {
    if (on_axis_x(p, c) < on_axis_x(p, d)) {
      b = d;
    } else {
      a = c;
    }
    c = b - g * (b - a);
    d = a + g * (b - a);
  }
  return 0.5 * (a + b);
}

/// Index of lattice cell (k1, k2) in Psf::intensity.
std::size_t cell(const Psf& psf, int k1, int k2) {
  const int h = psf.length / 2;
  return static_cast<std::size_t>(k2 + h) * static_cast<std::size_t>(psf.length) +
         static_cast<std::size_t>(k1 + h);
}

double sinc2(double w) {
  if (w == 0.0) return 1.0;
  const double t = std::sin(std::numbers::pi * w) / (std::numbers::pi * w);
  return t * t;
}

bool same_bits(double a, double b) {
  return std::memcmp(&a, &b, sizeof(double)) == 0;
}

}  // namespace

TEST_CASE("psf: the mapping M against the geometry and first_order (guard, ADR 0033, 10)",
          "[psf]") {
  // M is the scale of every expected value below; a wrong scale (factor 2, n') in M would be
  // shared by expected and computed values, so it is tied here to quantities that do not use M.
  // Tolerance: the central difference at delta = 1e-3 of s(p) = m p (1 - p^2 / (16 F#^2) + ...)
  // has the relative error (delta / (4 F#))^2, plus rounding.
  const double epd = 10.0;
  const CompiledSystem cs = paraboloid(epd);
  const PupilFunction p = pupil(cs, 16);
  REQUIRE(p.mapping.has_value());
  const auto& m = p.mapping->m;
  const double f_number_geometry = kF / epd;
  const double expected = 1.0 / (2.0 * f_number_geometry);
  const double tol = 2.0 * std::pow(1e-3 / (4.0 * f_number_geometry), 2) + 1e-12;
  INFO("M = [" << m[0] << ", " << m[1] << "; " << m[2] << ", " << m[3] << "]");
  CHECK(std::abs(std::abs(m[0]) / expected - 1.0) <= tol);
  CHECK(std::abs(std::abs(m[3]) / expected - 1.0) <= tol);
  CHECK(std::abs(m[1]) <= 1e-12 * expected);
  CHECK(std::abs(m[2]) <= 1e-12 * expected);
  // Concave mirror with an unrotated image surface: M = n'|u'| diag(-1, 1) (ADR 0032, point 4).
  CHECK(m[0] < 0.0);
  CHECK(m[3] > 0.0);
  // first_order: n' u' = -Phi y1 for a ray from infinity, y1 = D / 2.
  const auto fo = rtt::paraxial::first_order(cs, PathId{0}, 0);
  REQUIRE(fo.entrance_pupil.has_value());
  REQUIRE(fo.entrance_pupil->diameter.has_value());
  const double nu = std::abs(fo.power) * *fo.entrance_pupil->diameter / 2.0;
  CHECK(std::abs(std::abs(m[0]) / nu - 1.0) <= tol);
  CHECK(std::abs(p.mapping->s0[0]) <= 1e-15);
  CHECK(std::abs(p.mapping->s0[1]) <= 1e-15);
}

TEST_CASE("psf: the FFT lattice equals the direct sum psf_at", "[psf]") {
  // I_k = S(x_k) |sum f A|^2 / (L^2 E) (psf.hpp); compared up to 1e-12 of the peak.
  const CompiledSystem cs = paraboloid(10.0);
  PsfOptions o;
  o.pupil.grid = 32;
  o.padding = 4;
  const Psf psf = rtt::analysis::psf(cs, PathId{0}, 0, o);
  const PupilFunction p = pupil(cs, 32);
  double sum_fa = 0.0;
  double energy = 0.0;
  for (std::size_t j = 0; j < p.fraction.size(); ++j) {
    sum_fa += p.fraction[j] * p.amplitude[j];
    energy += p.fraction[j] * p.amplitude[j] * p.amplitude[j];
  }
  const double l = static_cast<double>(psf.length);
  const double scale = sum_fa * sum_fa / (l * l * energy);
  const double peak = psf.intensity[cell(psf, 0, 0)];
  REQUIRE(peak > 0.0);  // content guard
  for (int k2 = -6; k2 <= 6; k2 += 2) {
    for (int k1 = -6; k1 <= 6; ++k1) {
      const Point2 x{psf.lattice[0] * k1 + psf.lattice[1] * k2,
                     psf.lattice[2] * k1 + psf.lattice[3] * k2};
      const std::array<Point2, 1> xs{x};
      const double direct = rtt::analysis::psf_at(p, xs)[0] * scale;
      INFO("k = (" << k1 << ", " << k2 << ")");
      CHECK(std::abs(psf.intensity[cell(psf, k1, k2)] - direct) <= 1e-12 * peak);
    }
  }
}

TEST_CASE("psf: first zero of the Airy pattern at j11 lambda F#_w / pi", "[psf]") {
  const CompiledSystem cs = paraboloid(10.0);
  const PupilFunction p64 = pupil(cs, 64);
  const PupilFunction p128 = pupil(cs, 128);
  const double expected = j11() / std::numbers::pi * kLambdaMm * f_number(p128);
  const double r64 = first_null(p64);
  const double r128 = first_null(p128);
  const double e = 2.0 * std::abs(r64 - r128) / expected;
  INFO("r0: N=64 " << r64 << ", N=128 " << r128 << ", expected " << expected << ", e " << e);
  CHECK(e <= 1e-3);
  CHECK(std::abs(r128 / expected - 1.0) <= e);
  // Content guard: it is a zero of the pattern, not any minimum.
  CHECK(on_axis_x(p128, r128) < 1e-4);
}

TEST_CASE("psf: Airy profile [2 J1(v) / v]^2", "[psf]") {
  const CompiledSystem cs = paraboloid(10.0);
  const PupilFunction p64 = pupil(cs, 64);
  const PupilFunction p128 = pupil(cs, 128);
  const double s = kLambdaMm * f_number(p128) / std::numbers::pi;  // r per unit of v
  for (const double v : {1.0, 2.0, 3.0, 5.0}) {
    const double i64 = on_axis_x(p64, v * s);
    const double i128 = on_axis_x(p128, v * s);
    const double e = 2.0 * std::abs(i64 - i128) + 1e-12;
    INFO("v = " << v << ": N=64 " << i64 << ", N=128 " << i128 << ", Airy " << airy(v));
    CHECK(e <= 1e-3);
    CHECK(std::abs(i128 - airy(v)) <= e);
  }
}

TEST_CASE("psf: defocus Strehl sinc^2(w), through_focus and psf agree", "[psf]") {
  // F/50 (D = 2 mm): the wavefront of a focus shift is w rho^2 up to a relative rho^4 term below
  // alpha_m^2 / 12 + 2 (h / 2f)^2 ~ 6e-5, which changes S by less than |dS/dw| w 6e-5.
  // w from the exact sphere formula at the rim (ADR 0032, point 5), with n' = 1:
  // W = Delta (1 - cos alpha_m) + sqrt(R'^2 - Delta^2 sin^2 alpha_m) - R', both terms ~ rho^2;
  // sin alpha_m = h / (f + h^2 / (4 f)) of the paraboloid, h = D / 2. R' = |C' - E| with the focus
  // C at z = -f, C' = C + Delta (0, 0, -1) (the chief ray runs in -z after the mirror) and the exit
  // pupil E from first_order (on the side of the incoming light, so the formula's branch holds).
  // At Delta ~ 10 mm the second term is ~10 % of W: it must not be dropped.
  const double epd = 2.0;
  const CompiledSystem cs = paraboloid(epd);
  const double h = epd / 2.0;
  const double sin_a = h / (kF + h * h / (4.0 * kF));
  const double one_minus_cos = sin_a * sin_a / (1.0 + std::sqrt(1.0 - sin_a * sin_a));
  const auto fo = rtt::paraxial::first_order(cs, PathId{0}, 0);
  REQUIRE(fo.exit_pupil.has_value());
  REQUIRE(fo.exit_pupil->z.has_value());
  const double z_exit = *fo.exit_pupil->z;
  REQUIRE(z_exit > -kF);  // upstream of the focus along the reflected chief ray
  std::vector<double> shifts;
  std::vector<double> ws;
  for (const double w_target : {0.1, 0.3, 0.5, 0.8}) {
    const double delta = w_target * kLambdaMm / one_minus_cos;
    const double r = std::abs(-kF - delta - z_exit);
    const double t = delta * delta * sin_a * sin_a;
    const double w_mm = delta * one_minus_cos - t / (std::sqrt(r * r - t) + r);
    shifts.push_back(delta);
    ws.push_back(w_mm / kLambdaMm);
  }
  PupilFunctionOptions o64;
  o64.grid = 64;
  PupilFunctionOptions o128;
  o128.grid = 128;
  const std::vector<double> s64 = rtt::analysis::through_focus(cs, PathId{0}, 0, 0, shifts, o64);
  const std::vector<double> s128 = rtt::analysis::through_focus(cs, PathId{0}, 0, 0, shifts, o128);
  REQUIRE(s128.size() == shifts.size());
  for (std::size_t i = 0; i < shifts.size(); ++i) {
    const double w = ws[i];
    const double ds_dw = 2.0 * std::numbers::pi;  // bound of |d sinc^2 / dw| on [0, 1]
    const double e = 2.0 * std::abs(s64[i] - s128[i]) + ds_dw * w * 6e-5 + 1e-12;
    INFO("w = " << w << ": N=64 " << s64[i] << ", N=128 " << s128[i] << ", sinc^2 " << sinc2(w));
    CHECK(e <= 1e-3);
    CHECK(std::abs(s128[i] - sinc2(w)) <= e);
    CHECK(s128[i] < 1.0);  // content guard
    // psf().strehl with the same focus shift: the same formula on the same cells.
    PsfOptions po;
    po.pupil = o128;
    po.pupil.focus_shift = shifts[i];
    CHECK(same_bits(rtt::analysis::psf(cs, PathId{0}, 0, po).strehl, s128[i]));
  }
}

TEST_CASE("psf: energy, deficit and the ideal peak", "[psf]") {
  const CompiledSystem cs = paraboloid(10.0);
  PsfOptions o;
  o.pupil.grid = 64;
  const Psf psf = rtt::analysis::psf(cs, PathId{0}, 0, o);
  const PupilFunction p = pupil(cs, 64);
  double f1 = 0.0;
  double f2 = 0.0;
  for (const double f : p.fraction) {
    f1 += f;
    f2 += f * f;
  }
  const double deficit = 1.0 - f2 / f1;  // A = 1: ideal mirror, weight 1
  CHECK(std::abs(psf.edge_energy_deficit - deficit) <= 1e-14);
  CHECK(psf.edge_energy_deficit > 0.0);
  double sum = 0.0;
  for (const double v : psf.intensity) sum += v;
  CHECK(std::abs(sum - (1.0 - deficit)) <= 1e-12);
  // Ideal pupil: Strehl 1, peak in the centre cell.
  CHECK(std::abs(psf.strehl - 1.0) <= 1e-12);
  CHECK(psf.peak.k1 == 0);
  CHECK(psf.peak.k2 == 0);
}

TEST_CASE("psf: two wavelengths sum with weights w_i E_i on the reference lattice", "[psf]") {
  // lambda_2 = 1.25 lambda_1 with Q = 4: L_2 = 5 N gives exactly the lattice of the reference
  // (lambda / (L a m) equal), so the polychromatic PSF is the weighted sum of the single PSFs, the
  // second computed alone with padding 5.
  rtt::model::System s =
      rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m2/paraboloid_stop.rtt.json");
  s.aperture.value = rtt::model::Param(10.0);
  s.wavelengths = {{0.5876, 1.0, true}, {0.5876 * 1.25, 3.0, false}};
  const CompiledSystem cs = rtt::compile::compile(s, rtt::material::MaterialLibrary{});
  PsfOptions both;
  both.pupil.grid = 32;
  const Psf poly = rtt::analysis::psf(cs, PathId{0}, 0, both);
  PsfOptions first = both;
  first.wavelengths = std::vector<std::uint16_t>{0};
  PsfOptions second = both;
  second.wavelengths = std::vector<std::uint16_t>{1};
  second.padding = 5;
  const Psf a = rtt::analysis::psf(cs, PathId{0}, 0, first);
  const Psf b = rtt::analysis::psf(cs, PathId{0}, 0, second);
  REQUIRE(poly.wavelengths.size() == 2);
  CHECK(poly.wavelengths[1].length == 5 * 32);
  CHECK(poly.polychromatic_scale_error <= 1e-15);
  const double wa = 0.25 * a.wavelengths[0].energy;
  const double wb = 0.75 * b.wavelengths[0].energy;
  for (int k2 = -8; k2 <= 8; k2 += 4) {
    for (int k1 = -8; k1 <= 8; ++k1) {
      const double expected =
          (wa * a.intensity[cell(a, k1, k2)] + wb * b.intensity[cell(b, k1, k2)]) / (wa + wb);
      INFO("k = (" << k1 << ", " << k2 << ")");
      CHECK(std::abs(poly.intensity[cell(poly, k1, k2)] - expected) <=
            1e-14 * poly.intensity[cell(poly, 0, 0)]);
    }
  }
  // A ratio without a whole L: the scale error stays within the rounding of L_i.
  s.wavelengths = {{0.5876, 1.0, true}, {0.5876 * 1.1, 1.0, false}};
  const Psf rounded = rtt::analysis::psf(rtt::compile::compile(s, rtt::material::MaterialLibrary{}),
                                         PathId{0}, 0, both);
  const double l2 = static_cast<double>(rounded.wavelengths[1].length);
  CHECK(rounded.polychromatic_scale_error > 0.0);
  CHECK(rounded.polychromatic_scale_error <= 1.0 / (2.0 * l2) + 1e-15);
}

TEST_CASE("psf: warnings psf.aliasing and psf.pupil_undersampled", "[psf]") {
  const CompiledSystem cs = paraboloid(10.0);
  const auto has = [](const Psf& p, const char* code) {
    return std::any_of(p.warnings.begin(), p.warnings.end(),
                       [code](const auto& d) { return d.code == code; });
  };
  for (const int n : {32, 64, 128}) {
    PsfOptions o;
    o.pupil.grid = n;
    const Psf ideal = rtt::analysis::psf(cs, PathId{0}, 0, o);
    INFO("ideal pupil, N = " << n);
    CHECK_FALSE(has(ideal, "psf.aliasing"));
    CHECK_FALSE(has(ideal, "psf.pupil_undersampled"));
  }
  // Strong defocus at N = 32: w = 20 waves, the blur exceeds the field and the phase steps by
  // about 2 w rho (2 / N) = 2.5 waves between rim cells.
  const double h = 5.0;
  const double sin_a = h / (kF + h * h / (4.0 * kF));
  const double one_minus_cos = sin_a * sin_a / (1.0 + std::sqrt(1.0 - sin_a * sin_a));
  PsfOptions o;
  o.pupil.grid = 32;
  o.pupil.focus_shift = 20.0 * kLambdaMm / one_minus_cos;
  const Psf blurred = rtt::analysis::psf(cs, PathId{0}, 0, o);
  CHECK(has(blurred, "psf.aliasing"));
  CHECK(has(blurred, "psf.pupil_undersampled"));
}

TEST_CASE("psf: bitwise the same for 1, 4 and all threads", "[psf]") {
  const CompiledSystem cs = paraboloid(10.0);
  PsfOptions o;
  o.pupil.grid = 64;
  o.pupil.focus_shift = 0.5;
  const auto run = [&](int threads) {
    oneapi::tbb::task_arena arena(threads);
    return arena.execute([&] { return rtt::analysis::psf(cs, PathId{0}, 0, o); });
  };
  const Psf all = rtt::analysis::psf(cs, PathId{0}, 0, o);
  for (const int t : {1, 4}) {
    const Psf p = run(t);
    INFO("threads " << t);
    REQUIRE(p.intensity.size() == all.intensity.size());
    CHECK(std::memcmp(p.intensity.data(), all.intensity.data(),
                      all.intensity.size() * sizeof(double)) == 0);
    CHECK(same_bits(p.strehl, all.strehl));
  }
  CHECK(all.strehl < 1.0);  // content guard: the focus shift lowers the peak
}

TEST_CASE("psf: input errors", "[psf]") {
  const CompiledSystem cs = paraboloid(10.0);
  PsfOptions o;
  o.padding = 1;
  CHECK_THROWS_AS(rtt::analysis::psf(cs, PathId{0}, 0, o), std::invalid_argument);
  o.padding = 16;
  o.pupil.grid = 1024;  // L = 16384 > 8192
  CHECK_THROWS_AS(rtt::analysis::psf(cs, PathId{0}, 0, o), std::invalid_argument);
  PsfOptions w;
  w.wavelengths = std::vector<std::uint16_t>{};
  CHECK_THROWS_AS(rtt::analysis::psf(cs, PathId{0}, 0, w), std::invalid_argument);
  w.wavelengths = std::vector<std::uint16_t>{7};
  CHECK_THROWS_AS(rtt::analysis::psf(cs, PathId{0}, 0, w), std::invalid_argument);
  PsfOptions f;
  f.aliasing_factor = 0.5;
  CHECK_THROWS_AS(rtt::analysis::psf(cs, PathId{0}, 0, f), std::invalid_argument);
  const std::array<double, 1> bad{std::nan("")};
  CHECK_THROWS_AS(rtt::analysis::through_focus(cs, PathId{0}, 0, 0, bad), std::invalid_argument);
}

TEST_CASE("psf: an afocal image space has no PSF lattice", "[psf]") {
  // singlet_const with the object in its front focal point: the rays leave collimated, the
  // mapping is singular (ADR 0032, point 4), psf() gives AnalysisError.
  rtt::model::System s =
      rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m1/singlet_const.rtt.json");
  const CompiledSystem infinite = rtt::compile::compile(s, rtt::material::MaterialLibrary{});
  const auto fo = rtt::paraxial::first_order(infinite, PathId{0}, 0);
  REQUIRE(fo.front_focal_z.has_value());
  s.object.at_infinity = false;
  s.object.distance = rtt::model::Param(-*fo.front_focal_z);
  s.fields.points.resize(1);
  const CompiledSystem cs = rtt::compile::compile(s, rtt::material::MaterialLibrary{});
  CHECK_THROWS_AS(rtt::analysis::psf(cs, PathId{0}, 0), rtt::analysis::AnalysisError);
}
