// Zernike polynomials (Noll), the wavefront fit and the Strehl approximation (ADR 0032, points 6,
// 7 and 12; #203).
//
// Sources (read in full text, docs/quellen.md): K. Niu, C. Tian, "Zernike polynomials and their
// applications", J. Opt. 24 (2022) 123001, Eqs. (3), (4), (7), (8), (9), (54) and Table 10; J. C.
// Wyant, K. Creath, "Basic Wavefront Aberration Theory for Optical Metrology" (1992), Sec. IX,
// Eqs. (64), (67).

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numbers>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "rtt/analysis/pupil.hpp"
#include "rtt/analysis/zernike.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
#include "rtt/paraxial/paraxial.hpp"
#include "rtt/paraxial/seidel.hpp"
#include "rtt/trace/sources.hpp"

using rtt::analysis::AnalysisError;
using rtt::analysis::PupilFunction;
using rtt::analysis::PupilFunctionOptions;
using rtt::analysis::zernike;
using rtt::analysis::ZernikeOptions;
using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::model::Element;
using rtt::model::Param;
using rtt::model::System;

namespace {

constexpr PathId kPath{0};
constexpr double kU = std::numeric_limits<double>::epsilon() / 2.0;  // unit roundoff 2^-53

System load(const std::string& name) {
  return rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/" + name);
}

CompiledSystem compiled(const System& s) {
  rtt::material::MaterialLibrary lib;
  lib.add_catalog(std::string(RTT_CATALOG_DIR) + "/m2/schott.agf");
  return rtt::compile::compile(s, lib);
}

/// Exact coefficient c_s of R_n^m, Niu/Tian Eq. (4), as a double (exact below 2^53 for n <= 20).
double radial_coefficient(int n, int m, int s) {
  // (n - s)! / (s! ((n + m)/2 - s)! ((n - m)/2 - s)!) in uint64: 20! < 2^64, the quotient is an
  // integer below 2^53 for n <= 20 (ADR 0032, point 6), so exact as a double.
  std::uint64_t num = 1;
  for (int i = 1; i <= n - s; ++i) num *= static_cast<std::uint64_t>(i);
  std::uint64_t den = 1;
  for (int i = 1; i <= s; ++i) den *= static_cast<std::uint64_t>(i);
  for (int i = 1; i <= (n + m) / 2 - s; ++i) den *= static_cast<std::uint64_t>(i);
  for (int i = 1; i <= (n - m) / 2 - s; ++i) den *= static_cast<std::uint64_t>(i);
  const std::uint64_t quotient = num / den;
  const auto c = static_cast<double>(quotient);
  return s % 2 == 0 ? c : -c;
}

// Double-double arithmetic (error-free transformations with fma): a value is hi + lo with
// |lo| <= u |hi|. Each operation below has a relative error of at most about 4 u^2 = 5e-32.
struct Dd {
  double hi = 0.0;
  double lo = 0.0;
};

Dd two_sum(double a, double b) {
  const double s = a + b;
  const double bb = s - a;
  return {s, (a - (s - bb)) + (b - bb)};
}

Dd dd_mul(const Dd& a, const Dd& b) {
  const double p = a.hi * b.hi;
  const double e = std::fma(a.hi, b.hi, -p);
  const Dd r = two_sum(p, e + (a.hi * b.lo + a.lo * b.hi));
  return r;
}

Dd dd_add(const Dd& a, double b) {
  const Dd s = two_sum(a.hi, b);
  return two_sum(s.hi, s.lo + a.lo);
}

/// R_n^m(rho) in double-double: x = rho^2 and every Horner step without the double rounding.
Dd radial_dd(int n, int m, double rho) {
  const Dd r{rho, 0.0};
  const Dd x = dd_mul(r, r);
  const int k = (n - m) / 2;
  Dd acc{radial_coefficient(n, m, 0), 0.0};
  for (int s = 1; s <= k; ++s) acc = dd_add(dd_mul(acc, x), radial_coefficient(n, m, s));
  for (int i = 0; i < m; ++i) acc = dd_mul(acc, r);
  return acc;
}

/// A pupil function of the paraboloid (W = 0 up to rounding, all cells arrive) as carrier of
/// synthetic wavefronts; the fit reads only fraction, px, py, status and w.
PupilFunction carrier(const std::string& file,
                      std::uint16_t field,
                      std::uint16_t wavelength,
                      int grid) {
  PupilFunctionOptions o;
  o.grid = grid;
  return rtt::analysis::pupil_function(compiled(load(file)), kPath, field, wavelength, o);
}

/// W = sum c_j Z_j on the arrived cells of `pf`.
void set_wavefront(PupilFunction& pf, const std::vector<double>& c) {
  for (std::size_t k = 0; k < pf.fraction.size(); ++k) {
    if (pf.fraction[k] == 0.0 || pf.status[k] != rtt::trace::RayStatus::Alive) continue;
    const double rho = std::hypot(pf.px[k], pf.py[k]);
    const double theta = std::atan2(pf.py[k], pf.px[k]);
    double w = 0.0;
    for (std::size_t j = 0; j < c.size(); ++j)
      w += c[j] * zernike(static_cast<int>(j + 1), rho, theta);
    pf.w[k] = w;
  }
}

/// sigma_BF by hand: weighted least squares with Z1 ... Z4 (4 x 4 normal equations by Gaussian
/// elimination with partial pivoting, independent of the Cholesky of the code), then the
/// f-weighted RMS of the remainder.
double sigma_bf_by_hand(const PupilFunction& pf) {
  double a[4][5] = {};
  for (std::size_t k = 0; k < pf.fraction.size(); ++k) {
    if (pf.fraction[k] == 0.0 || pf.status[k] != rtt::trace::RayStatus::Alive) continue;
    const double rho = std::hypot(pf.px[k], pf.py[k]);
    const double theta = std::atan2(pf.py[k], pf.px[k]);
    double z[4];
    for (int j = 0; j < 4; ++j) z[j] = zernike(j + 1, rho, theta);
    for (int r = 0; r < 4; ++r) {
      for (int col = 0; col < 4; ++col) a[r][col] += pf.fraction[k] * z[r] * z[col];
      a[r][4] += pf.fraction[k] * z[r] * pf.w[k];
    }
  }
  for (int i = 0; i < 4; ++i) {
    int p = i;
    for (int r = i + 1; r < 4; ++r) {
      if (std::abs(a[r][i]) > std::abs(a[p][i])) p = r;
    }
    for (int col = 0; col < 5; ++col) std::swap(a[i][col], a[p][col]);
    for (int r = i + 1; r < 4; ++r) {
      const double q = a[r][i] / a[i][i];
      for (int col = i; col < 5; ++col) a[r][col] -= q * a[i][col];
    }
  }
  double c[4];
  for (int i = 3; i >= 0; --i) {
    double s = a[i][4];
    for (int col = i + 1; col < 4; ++col) s -= a[i][col] * c[col];
    c[i] = s / a[i][i];
  }
  double sum_f = 0.0;
  double sum = 0.0;
  for (std::size_t k = 0; k < pf.fraction.size(); ++k) {
    if (pf.fraction[k] == 0.0 || pf.status[k] != rtt::trace::RayStatus::Alive) continue;
    const double rho = std::hypot(pf.px[k], pf.py[k]);
    const double theta = std::atan2(pf.py[k], pf.px[k]);
    double r = pf.w[k];
    for (int j = 0; j < 4; ++j) r -= c[j] * zernike(j + 1, rho, theta);
    sum_f += pf.fraction[k];
    sum += pf.fraction[k] * r * r;
  }
  return std::sqrt(sum / sum_f);
}

}  // namespace

TEST_CASE("Noll indices: j -> (n, m) after Niu/Tian Eq. (8) and back after Eq. (9)", "[zernike]") {
  for (int j = 1; j <= rtt::analysis::kMaxZernikeTerms; ++j) {
    INFO("j = " << j);
    const auto nm = rtt::analysis::noll_to_nm(j);
    // Eq. (8) in double: n = floor(sqrt(2j - 1) + 0.5) - 1, m from the floors below.
    const int n = static_cast<int>(std::floor(std::sqrt(2.0 * j - 1.0) + 0.5)) - 1;
    const int m = n % 2 == 0 ? 2 * ((2 * j + 1 - n * (n + 1)) / 4)
                             : 2 * ((2 * (j + 1) - n * (n + 1)) / 4) - 1;
    CHECK(nm.n == n);
    CHECK(nm.m == m);
    CHECK((n - m) % 2 == 0);
    CHECK(m <= n);
    // Eq. (9) with the parity: even j cosine, odd j sine (Eq. (3)).
    CHECK(rtt::analysis::noll_from_nm(nm.n, nm.m, m != 0 && j % 2 == 1) == j);
  }
  CHECK(rtt::analysis::noll_to_nm(37).n == 8);
  CHECK(rtt::analysis::noll_to_nm(231).n == 20);
  CHECK_THROWS_AS(rtt::analysis::noll_to_nm(0), std::invalid_argument);
  CHECK_THROWS_AS(rtt::analysis::noll_to_nm(232), std::invalid_argument);
}

TEST_CASE("Zernike polynomials Z1 ... Z11 against the explicit Noll formulas", "[zernike]") {
  // Niu/Tian Eq. (3), (4); theta from +x counter-clockwise (p. 5-6). Tolerance 1e-15: a few
  // roundings of values below 7.
  for (const double rho : {0.0, 0.3, 0.71, 1.0}) {
    for (const double t : {0.0, 0.4, 2.0, -1.3}) {
      INFO("rho " << rho << ", theta " << t);
      const double r2 = rho * rho;
      const double expected[11] = {1.0,
                                   2.0 * rho * std::cos(t),
                                   2.0 * rho * std::sin(t),
                                   std::sqrt(3.0) * (2.0 * r2 - 1.0),
                                   std::sqrt(6.0) * r2 * std::sin(2.0 * t),
                                   std::sqrt(6.0) * r2 * std::cos(2.0 * t),
                                   std::sqrt(8.0) * (3.0 * r2 * rho - 2.0 * rho) * std::sin(t),
                                   std::sqrt(8.0) * (3.0 * r2 * rho - 2.0 * rho) * std::cos(t),
                                   std::sqrt(8.0) * r2 * rho * std::sin(3.0 * t),
                                   std::sqrt(8.0) * r2 * rho * std::cos(3.0 * t),
                                   std::sqrt(5.0) * (6.0 * r2 * r2 - 6.0 * r2 + 1.0)};
      for (int j = 1; j <= 11; ++j) {
        CHECK(std::abs(zernike(j, rho, t) - expected[j - 1]) <= 1e-14);
      }
    }
  }
}

TEST_CASE("Zernike polynomials are orthonormal over the unit disk (j, k <= 37)", "[zernike]") {
  // Niu/Tian Eq. (7): the mean of Z_j Z_k over the disk is delta_jk. The Gaussian pupil quadrature
  // of #168 is exact for u^p e^{i M theta} with p <= 2 rings - 1 and |M| < arms (sources.hpp). For
  // n <= 8 the product Z_j Z_k has frequencies |M| <= 16 and, at M = 0, a polynomial in u = rho^2
  // of degree <= 8: rings = 6, arms = 20 integrate it exactly. Tolerance 1e-13: 120 points, values
  // up to 3, products up to 9, rounding about 120 u 9.
  const rtt::trace::GaussPupil gauss{6, 20};
  const auto points = rtt::trace::pupil_points(gauss);
  const auto weights = rtt::trace::gauss_pupil_weights(gauss);
  for (int j = 1; j <= 37; ++j) {
    for (int k = j; k <= 37; ++k) {
      double mean = 0.0;
      for (std::size_t i = 0; i < points.size(); ++i) {
        const double rho = std::hypot(points[i].px, points[i].py);
        const double theta = std::atan2(points[i].py, points[i].px);
        mean += weights[i] * zernike(j, rho, theta) * zernike(k, rho, theta);
      }
      INFO("j " << j << ", k " << k);
      CHECK(std::abs(mean - (j == k ? 1.0 : 0.0)) <= 1e-13);
    }
  }
}

TEST_CASE("R_n^m within its rounding bound for n <= 20 (ADR 0032, point 6)", "[zernike]") {
  // Bound of the ADR: |R^ - R| <= gamma_k sum_s |c_s| rho^(n - 2s), k = 3K + m + 1, K = (n - m)/2,
  // gamma_k = k u / (1 - k u). Reference: the same rho (a double, taken as exact) in double-double
  // (radial_dd above). Its own error: 2K + m + 1 double-double operations, each with a relative
  // error of at most 4 u^2 against magnitudes bounded by sum_s |c_s| rho^(n - 2s), i.e. at most
  // (2K + m + 1) 4 u^2 S <= 100 u^2 S = 1.2e-30 S, against gamma_k S >= u S: below 1e-14 of the
  // bound. The check adds it to the bound anyway.
  std::size_t checked = 0;
  double worst = 0.0;  // largest error relative to the bound, reported
  for (int n = 0; n <= 20; ++n) {
    for (int m = n % 2; m <= n; m += 2) {
      const int big_k = (n - m) / 2;
      const int k = 3 * big_k + m + 1;
      const double gamma = k * kU / (1.0 - k * kU);
      for (int i = 0; i <= 1000; ++i) {
        const double rho = i * 1e-3;
        double s = 0.0;  // sum |c_s| rho^(n - 2s)
        for (int t = 0; t <= big_k; ++t) {
          s += std::abs(radial_coefficient(n, m, t)) * std::pow(rho, n - 2 * t);
        }
        const Dd exact = radial_dd(n, m, rho);
        const double got = rtt::analysis::zernike_radial(n, m, rho);
        const double err = std::abs((got - exact.hi) - exact.lo);
        const double bound = gamma * s * (1.0 + 1e-12) + 100.0 * kU * kU * s;
        INFO("n " << n << ", m " << m << ", rho " << rho);
        CHECK(err <= bound);
        if (bound > 0.0) worst = std::max(worst, err / bound);
        ++checked;
      }
    }
  }
  INFO("largest error / bound " << worst);
  CHECK(checked == 121u * 1001u);
  CHECK(worst > 0.0);  // content guard: the cancellation at large n is visible
}

TEST_CASE("zernike fit recovers synthetic coefficients (j <= 37)", "[zernike][fit]") {
  // W = sum c_j Z_j is in the span of the 37 fitted terms, so the weighted least squares
  // reproduce c up to rounding: about cond(A^T F A) u max|c|, cond near 1 on the full disk.
  // Tolerance 1e-12 max|c| (ADR 0032, point 12).
  PupilFunction pf = carrier("m2/paraboloid_stop.rtt.json", 0, 0, 32);
  std::vector<double> c(37);
  for (std::size_t j = 0; j < c.size(); ++j) {
    c[j] = 0.1 * ((j % 3 == 0) ? -1.0 : 1.0) / static_cast<double>(j + 1);
  }
  set_wavefront(pf, c);
  const auto fit = rtt::analysis::zernike_fit(pf);
  REQUIRE(fit.coefficients.size() == 37);
  for (std::size_t j = 0; j < c.size(); ++j) {
    INFO("j = " << j + 1);
    CHECK(std::abs(fit.coefficients[j] - c[j]) <= 1e-12 * 0.1);
  }
  CHECK(fit.residual_rms <= 1e-12);
  CHECK(fit.cells > 700);
  CHECK(fit.arrived_fraction == 1.0);
  CHECK(fit.path == kPath);
  CHECK(fit.warnings.empty());
}

TEST_CASE("zernike fit of the plano-convex lens: c11 and c4 from W040 = S_I / 8",
          "[zernike][fit]") {
  // m1/singlet_const in vacuum at EPD 10 mm (NA 0.05; the M2 Seidel case of #34 used EPD 2, where
  // W040 is only 5e-4 waves), image surface in the paraxial focus, axis, both orientations. On the
  // axis W = W040 rho^4 + W060 rho^6 + W080 rho^8 + ..., W040 = S_I / 8 exactly (Taylor
  // coefficient). Exact decompositions (Niu/Tian Eqs. (3), (4)): rho^4 = Z11/(6 sqrt5) + Z4/(2
  // sqrt3) + 1/3 and rho^6 = Z22/(20 sqrt7) + Z11/(4 sqrt5) + 9 Z4/(20 sqrt3) + 1/4, so with W060 =
  // 20 sqrt7 c22:
  //   c11 = W040/(6 sqrt5) + sqrt35 c22,  c4 = W040/(2 sqrt3) + 3 sqrt21 c22
  // up to the terms from rho^8 on. Tolerance (ADR 0032, point 12): the rho^8 part, estimated from
  // c37 = W080 <rho^8, Z37> / ... through the exact projections <rho^8, Z_j> (the mean of rho^8
  // Z_j over the disk, computed below from Eq. (4)), times 2, plus 1e-6 |W040| for the rounding of
  // S_I and of the OPD (about 1e-9 waves against W040 of order 0.1 wave).
  const auto projection = [](int p, int j) {
    // <rho^(2p), Z_j> for m = 0: sqrt(n + 1) sum_s c_s 2 / (2p + n - 2s + 2).
    const auto nm = rtt::analysis::noll_to_nm(j);
    double sum = 0.0;
    for (int s = 0; s <= nm.n / 2; ++s) {
      sum += radial_coefficient(nm.n, 0, s) * 2.0 / (2.0 * p + nm.n - 2.0 * s + 2.0);
    }
    return std::sqrt(nm.n + 1.0) * sum;
  };
  // Check of the projections against the decompositions of the ADR.
  CHECK(std::abs(projection(2, 11) - 1.0 / (6.0 * std::sqrt(5.0))) <= 1e-15);
  CHECK(std::abs(projection(3, 22) - 1.0 / (20.0 * std::sqrt(7.0))) <= 1e-15);
  CHECK(std::abs(projection(3, 4) - 9.0 / (20.0 * std::sqrt(3.0))) <= 1e-15);
  const auto check = [&](System s) {
    s.environment.medium = "VACUUM";
    s.aperture = {rtt::model::SystemApertureType::EntrancePupilDiameter, Param(10.0)};
    const CompiledSystem initial = compiled(s);
    const std::uint16_t ref = initial.reference_wavelength();
    std::get<Element>(s.root.children[2].value).pose =
        rtt::model::Pose::along_z(*rtt::paraxial::first_order(initial, kPath, ref).rear_focal_z);
    const CompiledSystem cs = compiled(s);
    const double lambda_mm = cs.wavelengths_um()[ref] * 1e-3;
    const double w040 = rtt::paraxial::seidel(cs, kPath, ref).sum.s1 / 8.0 / lambda_mm;  // waves
    PupilFunctionOptions o;
    o.grid = 64;
    const auto fit =
        rtt::analysis::zernike_fit(rtt::analysis::pupil_function(cs, kPath, 0, ref, o));
    const auto& c = fit.coefficients;
    const double c22 = c[21];
    const double w080 = c[36] / projection(4, 37);
    const double e11 = w080 * (projection(4, 11) - std::sqrt(35.0) * projection(4, 22));
    const double e4 = w080 * (projection(4, 4) - 3.0 * std::sqrt(21.0) * projection(4, 22));
    const double expected11 = w040 / (6.0 * std::sqrt(5.0)) + std::sqrt(35.0) * c22;
    const double expected4 = w040 / (2.0 * std::sqrt(3.0)) + 3.0 * std::sqrt(21.0) * c22;
    INFO("W040 " << w040 << " waves, c11 " << c[10] << ", c4 " << c[3] << ", c22 " << c22
                 << ", c37 " << c[36]);
    CHECK(std::abs(c[10] - expected11) <= 2.0 * std::abs(e11) + 1e-6 * std::abs(w040));
    CHECK(std::abs(c[3] - expected4) <= 2.0 * std::abs(e4) + 1e-6 * std::abs(w040));
    CHECK(std::abs(w040) > 0.1);  // content guard: a measurable spherical aberration
    return w040;
  };
  const double curved_first = check(load("m1/singlet_const.rtt.json"));
  System flipped = load("m1/singlet_const.rtt.json");
  Element& lens = std::get<Element>(flipped.root.children[1].value);
  lens.surfaces[0].shape.base = rtt::model::Plane{};
  lens.surfaces[1].shape.base = rtt::model::Conic{Param(-51.68), Param(0.0)};
  const double flat_first = check(flipped);
  CHECK(std::abs(curved_first - flat_first) > 1e-3);  // the orientation changes W040
}

TEST_CASE("Strehl approximation after Marechal from sigma about the diffraction focus",
          "[zernike][strehl]") {
  // Wyant/Creath Eq. (67): S = exp(-(2 pi sigma)^2), sigma about the diffraction focus (piston,
  // tilt and defocus removed, Sec. IX), in waves at the evaluated wavelength lambda, i.e.
  // sigma_BF lambda_ref / lambda for W in waves at lambda_ref. Cooke triplet, wavelength 0 (the
  // reference is 1), so the factor lambda_ref / lambda differs from 1.
  PupilFunction pf = carrier("m2/cooke_triplet.rtt.json", 0, 0, 24);
  std::vector<double> spherical(11, 0.0);
  spherical[10] = 0.05;  // a Z11
  set_wavefront(pf, spherical);
  const auto a = rtt::analysis::zernike_fit(pf);
  const double sigma_bf = sigma_bf_by_hand(pf);
  CHECK(std::abs(a.sigma_bf - sigma_bf) <= 1e-13 * sigma_bf);
  const double ratio = pf.reference_wavelength_um / pf.wavelength_um;
  REQUIRE(std::abs(ratio - 1.0) > 0.1);
  const double expected = std::exp(-std::pow(2.0 * std::numbers::pi * sigma_bf * ratio, 2.0));
  CHECK(std::abs(a.strehl - expected) <= 1e-14);
  // Tilt and defocus fall out exactly: the same sigma_BF with them added.
  std::vector<double> with_tilt = spherical;
  with_tilt[1] = 0.3;   // Z2
  with_tilt[2] = -0.2;  // Z3
  with_tilt[3] = 0.4;   // Z4
  set_wavefront(pf, with_tilt);
  const auto b = rtt::analysis::zernike_fit(pf);
  CHECK(std::abs(b.sigma_bf - a.sigma_bf) <= 1e-13 * a.sigma_bf);
  CHECK(b.sigma > 2.0 * a.sigma);  // content guard: sigma itself grows with tilt and defocus
  // sigma: the f-weighted standard deviation of W (piston removed), checked by hand.
  double sum_f = 0.0;
  double sum_fw = 0.0;
  for (std::size_t k = 0; k < pf.fraction.size(); ++k) {
    if (pf.fraction[k] == 0.0 || pf.status[k] != rtt::trace::RayStatus::Alive) continue;
    sum_f += pf.fraction[k];
    sum_fw += pf.fraction[k] * pf.w[k];
  }
  const double mean = sum_fw / sum_f;
  double var = 0.0;
  for (std::size_t k = 0; k < pf.fraction.size(); ++k) {
    if (pf.fraction[k] == 0.0 || pf.status[k] != rtt::trace::RayStatus::Alive) continue;
    var += pf.fraction[k] * (pf.w[k] - mean) * (pf.w[k] - mean);
  }
  CHECK(std::abs(b.sigma - std::sqrt(var / sum_f)) <= 1e-14);
}

TEST_CASE("zernike fit: errors and the incomplete pupil", "[zernike][fit]") {
  SECTION("fewer arrived cells than terms") {
    const PupilFunction pf = carrier("m2/paraboloid_stop.rtt.json", 0, 0, 4);
    try {
      static_cast<void>(rtt::analysis::zernike_fit(pf));
      FAIL("no AnalysisError");
    } catch (const AnalysisError& e) {
      CHECK(e.location() == std::optional<std::string>("/paths/0"));
    }
  }
  SECTION("cells on one row: the normal equations are singular") {
    // Only the row py = 1/N arrives: there Z3 = 2 py is a multiple of Z1 (pivot below 1e-8 of the
    // largest diagonal element, ADR 0032, point 6).
    PupilFunction pf = carrier("m2/paraboloid_stop.rtt.json", 0, 0, 64);
    std::size_t kept = 0;
    for (std::size_t k = 0; k < pf.fraction.size(); ++k) {
      if (k / 64 != 32) {
        pf.status[k] = rtt::trace::RayStatus::Vignetted;
      } else if (pf.fraction[k] > 0.0) {
        ++kept;
      }
    }
    REQUIRE(kept >= 37);
    try {
      static_cast<void>(rtt::analysis::zernike_fit(pf));
      FAIL("no AnalysisError");
    } catch (const AnalysisError& e) {
      CHECK(e.location() == std::optional<std::string>("/paths/0"));
    }
  }
  SECTION("a nearly dependent column: a pivot below the floor") {
    // Three terms on the row py = 1/N, where Z3 = 2 py is a multiple of Z1, plus one cell of
    // another row with the tiny weight f = 1e-12: the normal equations are regular, but their
    // pivot for Z3 is about 1e-12 of the largest diagonal element, below the floor 1e-8 of
    // ADR 0032, point 6 (without the floor the fit would return coefficients of about 1e12).
    PupilFunction pf = carrier("m2/paraboloid_stop.rtt.json", 0, 0, 64);
    for (std::size_t k = 0; k < pf.fraction.size(); ++k) {
      if (k / 64 != 32 && k != 40 * 64 + 32) pf.status[k] = rtt::trace::RayStatus::Vignetted;
    }
    pf.fraction[40 * 64 + 32] = 1e-12;
    ZernikeOptions o;
    o.terms = 3;
    try {
      static_cast<void>(rtt::analysis::zernike_fit(pf, o));
      FAIL("no AnalysisError");
    } catch (const AnalysisError& e) {
      CHECK(e.location() == std::optional<std::string>("/paths/0"));
    }
  }
  SECTION("a clipped pupil gives the warning, the full one none") {
    // Stop 28 mm instead of 30: about 13 % of the area is lost, above 1 %. (A strongly clipped
    // pupil, e.g. 10 mm, leaves a disk of radius 1/3, on which 37 terms are too badly
    // conditioned: an AnalysisError, as it should.)
    System s = load("m2/paraboloid_stop.rtt.json");
    std::get<Element>(s.root.children[0].value).surfaces[0].aperture =
        rtt::model::CircularAperture{28.0, 0.0};
    PupilFunctionOptions o;
    o.grid = 64;
    const auto fit =
        rtt::analysis::zernike_fit(rtt::analysis::pupil_function(compiled(s), kPath, 0, 0, o));
    REQUIRE(fit.warnings.size() == 1);
    CHECK(fit.warnings[0].code == "zernike.pupil_incomplete");
    CHECK(fit.warnings[0].location == "/paths/0");
    CHECK(rtt::analysis::zernike_fit(carrier("m2/paraboloid_stop.rtt.json", 0, 0, 32))
              .warnings.empty());
  }
  SECTION("options out of range") {
    const PupilFunction pf = carrier("m2/paraboloid_stop.rtt.json", 0, 0, 16);
    for (const int terms : {0, -1, rtt::analysis::kMaxZernikeTerms + 1}) {
      ZernikeOptions o;
      o.terms = terms;
      CHECK_THROWS_AS(rtt::analysis::zernike_fit(pf, o), std::invalid_argument);
    }
    for (const double f : {-0.1, 1.5, std::numeric_limits<double>::quiet_NaN()}) {
      ZernikeOptions o;
      o.incomplete_warning_fraction = f;
      CHECK_THROWS_AS(rtt::analysis::zernike_fit(pf, o), std::invalid_argument);
    }
    CHECK_THROWS_AS(rtt::analysis::zernike(0, 0.5, 0.0), std::invalid_argument);
    CHECK_THROWS_AS(rtt::analysis::zernike_radial(21, 1, 0.5), std::invalid_argument);
    CHECK_THROWS_AS(rtt::analysis::zernike_radial(4, 1, 0.5), std::invalid_argument);
  }
}
