// M6 acceptance (#212): every criterion of the roadmap row through the tracer, on the reference
// systems under tests/reference/m6/ and r2/, with expected values that do not come from the code
// under test and tolerances fixed in advance (plan A1, docs/architecture.md "Abnahme M6").
//
// Tolerance rule of ADR 0033, point 10 (with the erratum #206): X on the grids N and 2N, compare
// X_2N, bound e = 2 |X_N - X_2N|, and e itself must stay below 1e-3; for a curve the largest e over
// its points (relative where the acceptance is relative). Every expected value that uses F#_w is
// preceded by a guard that ties the mapping M to quantities that do not use it.
//
// Cases: 1 Airy zero, 2 encircled energy and 3 MTF (F/10), 4 defocus Strehl (F/50), 5 Zernike
// from W040 (plano-convex), 6 and 7 Jones pupil of a 45 deg metal fold (collimated, converging),
// 8 ideal lens (Airy zero, Strehl 1), and the reproducibility at 1, 4 and all threads. The
// annular pupil (9) follows with #229.

#include <tbb/global_control.h>

#include <algorithm>
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <numbers>
#include <optional>
#include <string>
#include <variant>
#include <vector>

#include "rtt/analysis/jones.hpp"
#include "rtt/analysis/psf.hpp"
#include "rtt/analysis/pupil.hpp"
#include "rtt/analysis/zernike.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/math/types.hpp"
#include "rtt/model/model.hpp"
#include "rtt/paraxial/paraxial.hpp"
#include "rtt/paraxial/seidel.hpp"
#include "rtt/polar/fresnel.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/ray_paths.hpp"
#include "rtt/trace/sequential.hpp"

using Cx = std::complex<double>;
using rtt::analysis::JonesPupil;
using rtt::analysis::Point2;
using rtt::analysis::PsfOptions;
using rtt::analysis::PupilFunction;
using rtt::analysis::PupilFunctionOptions;
using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::math::Vec3;
using rtt::model::Element;
using rtt::model::System;
using rtt::trace::RayBatch;
using rtt::trace::RayStatus;

namespace {

constexpr double kPi = std::numbers::pi;
constexpr double kF = 100.0;             // focal length of the paraboloid, R / 2, mm
constexpr double kLambdaMm = 0.5876e-3;  // reference wavelength of the m6 paraboloids, mm
const Cx kMetal(1.2, 7.26);              // test metal of the m6 folds (as M3)

System load(const std::string& relative) {
  return rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/" + relative);
}

CompiledSystem compiled(const System& s) {
  return rtt::compile::compile(s, rtt::material::MaterialLibrary{});
}

PupilFunction pupil(const CompiledSystem& cs, PathId path, int n) {
  PupilFunctionOptions o;
  o.grid = n;
  return rtt::analysis::pupil_function(cs, path, 0, cs.reference_wavelength(), o);
}

/// n' = |Re n| of the image medium, read from the compiled media (not through the analysis).
double image_index(const CompiledSystem& cs, PathId path) {
  const auto& last = cs.path(path).events.back();
  return std::abs(cs.media()[last.medium_after].index[cs.reference_wavelength()].real());
}

/// F#_w = 1 / (2 |M_11|) of the mapping.
double f_number(const PupilFunction& p) {
  REQUIRE(p.mapping.has_value());
  return 1.0 / (2.0 * std::abs(p.mapping->m[0]));
}

/// Guard of ADR 0033, point 10: |M_11| = |M_22| = n' |u'| with the paraxial image-space slope u'
/// of the marginal ray, given from the geometry of the system; M diagonal. The central difference
/// at delta = 1e-3 of s(p) = m p (1 + O(p^2 / F#^2)) has the relative error (delta / (4 F#))^2.
void guard_mapping(const PupilFunction& p, double n_u) {
  REQUIRE(p.mapping.has_value());
  const auto& m = p.mapping->m;
  const double f_number_geometry = 1.0 / (2.0 * n_u);
  const double tol = 2.0 * std::pow(1e-3 / (4.0 * f_number_geometry), 2) + 1e-12;
  INFO("M = [" << m[0] << ", " << m[1] << "; " << m[2] << ", " << m[3] << "], n'u' = " << n_u);
  REQUIRE(std::abs(std::abs(m[0]) / n_u - 1.0) <= tol);
  REQUIRE(std::abs(std::abs(m[3]) / n_u - 1.0) <= tol);
  REQUIRE(std::abs(m[1]) <= 1e-12 * n_u);
  REQUIRE(std::abs(m[2]) <= 1e-12 * n_u);
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

/// Intensity in Strehl units at distance r (mm) along e_x.
double along_x(const PupilFunction& p, double r) {
  const std::array<Point2, 1> x{Point2{r, 0.0}};
  return rtt::analysis::psf_at(p, x)[0];
}

/// First minimum of the PSF along e_x in [1.0, 1.45] lambda F#_w by golden section, mm.
double first_null(const PupilFunction& p, double lambda_mm) {
  const double s = lambda_mm * f_number(p);
  double a = 1.0 * s;
  double b = 1.45 * s;
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

/// The Airy zero on the grids 64 and 128 against j11 lambda F#_w / pi (case 1 and 8).
void check_airy_zero(const CompiledSystem& cs, PathId path, double n_u) {
  const double lambda_mm = cs.wavelengths_um()[cs.reference_wavelength()] * 1e-3;
  const PupilFunction p64 = pupil(cs, path, 64);
  const PupilFunction p128 = pupil(cs, path, 128);
  guard_mapping(p128, n_u);
  const double expected = j11() / kPi * lambda_mm * f_number(p128);
  const double r64 = first_null(p64, lambda_mm);
  const double r128 = first_null(p128, lambda_mm);
  const double e = 2.0 * std::abs(r64 - r128) / expected;
  INFO("r0: N=64 " << r64 << ", N=128 " << r128 << ", expected " << expected << ", e " << e);
  CHECK(e <= 1e-3);
  CHECK(std::abs(r128 / expected - 1.0) <= e);
  CHECK(along_x(p128, r128) < 1e-4);  // content guard: a zero of the pattern, not any minimum
}

/// Diffraction-limited MTF of the circular pupil, (2 / pi) [acos x - x sqrt(1 - x^2)].
double mtf_circle(double x) {
  return 2.0 / kPi * (std::acos(x) - x * std::sqrt(1.0 - x * x));
}

/// Encircled energy of the Airy pattern, 1 - J0(v)^2 - J1(v)^2, v = pi r / (lambda F#_w).
double ee_airy(double v) {
  const double j0 = std::cyl_bessel_j(0.0, v);
  const double j1 = std::cyl_bessel_j(1.0, v);
  return 1.0 - j0 * j0 - j1 * j1;
}

double sinc2(double w) {
  if (w == 0.0) return 1.0;
  const double t = std::sin(kPi * w) / (kPi * w);
  return t * t;
}

/// (n - s)! / (s! ((n + m)/2 - s)! ((n - m)/2 - s)!) with the sign (-1)^s (Noll; ADR 0032,
/// point 6), exact in uint64 for n <= 20.
double radial_coefficient(int n, int m, int s) {
  std::uint64_t num = 1;
  for (int i = 1; i <= n - s; ++i) num *= static_cast<std::uint64_t>(i);
  std::uint64_t den = 1;
  for (int i = 1; i <= s; ++i) den *= static_cast<std::uint64_t>(i);
  for (int i = 1; i <= (n + m) / 2 - s; ++i) den *= static_cast<std::uint64_t>(i);
  for (int i = 1; i <= (n - m) / 2 - s; ++i) den *= static_cast<std::uint64_t>(i);
  const auto c = static_cast<double>(num / den);
  return s % 2 == 0 ? c : -c;
}

/// <rho^(2p), Z_j> over the unit disk for m = 0: sqrt(n + 1) sum_s c_s 2 / (2p + n - 2s + 2).
double projection(int p, int j) {
  const auto nm = rtt::analysis::noll_to_nm(j);
  double sum = 0.0;
  for (int s = 0; s <= nm.n / 2; ++s) {
    sum += radial_coefficient(nm.n, 0, s) * 2.0 / (2.0 * p + nm.n - 2.0 * s + 2.0);
  }
  return std::sqrt(nm.n + 1.0) * sum;
}

Element& child(System& s, std::size_t i) {
  return std::get<Element>(s.root.children[i].value);
}

bool same_bits(const std::vector<double>& a, const std::vector<double>& b) {
  return a.size() == b.size() &&
         (a.empty() || std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0);
}

bool same_bits(const std::vector<Cx>& a, const std::vector<Cx>& b) {
  return a.size() == b.size() &&
         (a.empty() || std::memcmp(a.data(), b.data(), a.size() * sizeof(Cx)) == 0);
}

/// 5 x 5 parallel rays along +z from z = 0 at x, y in {-3, -1.5, 0, 1.5, 3} (inside D = 6),
/// labelled with pupil coordinates (x, y) / 3; the chief ray (0, 0) has index 12.
RayBatch parallel_rays() {
  RayBatch rays(25);
  for (std::size_t i = 0; i < 25; ++i) {
    const double x = -3.0 + 1.5 * static_cast<double>(i % 5);
    const double y = -3.0 + 1.5 * static_cast<double>(i / 5);
    rays.pos_x()[i] = x;
    rays.pos_y()[i] = y;
    rays.pupil_x()[i] = x / 3.0;
    rays.pupil_y()[i] = y / 3.0;
  }
  return rays;
}

constexpr std::size_t kChief = 12;

/// Collimated start rays on a grid of pitch 1.25 mm inside r = 5 (F/5 after f = 50), from
/// z = -1 in front of the ideal lens at z = 0 (AGENTS.md); `chief` is set to the ray at (0, 0).
RayBatch converging_rays(std::size_t& chief) {
  RayBatch start(0);
  for (int iy = -4; iy <= 4; ++iy) {
    for (int ix = -4; ix <= 4; ++ix) {
      const double x = 1.25 * ix;
      const double y = 1.25 * iy;
      if (x * x + y * y > 25.0) continue;
      const std::size_t i = start.size();
      start.resize(i + 1);
      start.pos_x()[i] = x;
      start.pos_y()[i] = y;
      start.pos_z()[i] = -1.0;
      start.pupil_x()[i] = x / 5.0;
      start.pupil_y()[i] = y / 5.0;
      if (ix == 0 && iy == 0) chief = i;
    }
  }
  return start;
}

}  // namespace

TEST_CASE("M6 case 1: Airy zero at j11 lambda F#_w / pi, paraboloid F/10", "[m6]") {
  // tests/reference/m6/paraboloid_f10.rtt.json: f = 100 mm, D = 10 mm, AIR, on axis; the
  // paraboloid images the axial point without aberration (W = 0), so the PSF is the Airy pattern
  // of the linear model. Guard: n'|u'| = n' D / (2 f) for the object at infinity, and first_order
  // gives the same with n' u' = -Phi y1. Grids 64/128 (ADR 0033, point 10).
  const CompiledSystem cs = compiled(load("m6/paraboloid_f10.rtt.json"));
  const double n_image = image_index(cs, PathId{0});
  REQUIRE(n_image > 1.0);  // AIR: n' carries into M and into F#_w
  const auto fo = rtt::paraxial::first_order(cs, PathId{0}, cs.reference_wavelength());
  REQUIRE(fo.entrance_pupil.has_value());
  REQUIRE(fo.entrance_pupil->diameter.has_value());
  const double n_u = n_image * 10.0 / (2.0 * kF);
  CHECK(std::abs(std::abs(fo.power) * *fo.entrance_pupil->diameter / 2.0 / n_u - 1.0) <= 1e-12);
  check_airy_zero(cs, PathId{0}, n_u);
}

TEST_CASE("M6 case 2: encircled energy 1 - J0^2 - J1^2, paraboloid F/10", "[m6]") {
  // tests/reference/m6/paraboloid_f10.rtt.json, W = 0: EE(r) of the Airy pattern with
  // v = pi r / (lambda F#_w) at r = 0.5, 1, r0 = j11 / pi, 2, 3 and 5 lambda F#_w. Bound per
  // curve, relative as the acceptance (ADR 0033, point 10 with the addendum of #229): e_max =
  // max 2 |X_N - X_2N| / X <= 1e-3 on the grids 128/256 (erratum #206), and every radius against
  // e_max; EE(r0) on its own lies in (0.83, 0.84) (content guard).
  const CompiledSystem cs = compiled(load("m6/paraboloid_f10.rtt.json"));
  const double n_u = image_index(cs, PathId{0}) * 10.0 / (2.0 * kF);
  const PupilFunction p = pupil(cs, PathId{0}, 16);
  guard_mapping(p, n_u);
  const double unit = kLambdaMm * f_number(p);  // lambda F#_w, mm
  const std::vector<double> units{0.5, 1.0, j11() / kPi, 2.0, 3.0, 5.0};
  std::vector<double> radii;
  for (const double u : units) radii.push_back(u * unit);
  rtt::analysis::EncircledEnergyOptions o128;
  o128.pupil.grid = 128;
  rtt::analysis::EncircledEnergyOptions o256;
  o256.pupil.grid = 256;
  const auto a = rtt::analysis::encircled_energy(cs, PathId{0}, 0, radii, o128);
  const auto b = rtt::analysis::encircled_energy(cs, PathId{0}, 0, radii, o256);
  REQUIRE(a.energy.size() == radii.size());
  REQUIRE(b.energy.size() == radii.size());
  double e_max = 0.0;
  for (std::size_t i = 0; i < radii.size(); ++i) {
    e_max = std::max(e_max, 2.0 * std::abs(a.energy[i] - b.energy[i]) / ee_airy(kPi * units[i]));
  }
  INFO("e_max (relative) " << e_max);
  CHECK(e_max <= 1e-3);
  for (std::size_t i = 0; i < radii.size(); ++i) {
    const double expected = ee_airy(kPi * units[i]);
    INFO("r = " << units[i] << " lambda F#_w: N=256 " << b.energy[i] << ", expected " << expected);
    CHECK(std::abs(b.energy[i] - expected) / expected <= e_max + 1e-12);
  }
  CHECK(b.energy[2] > 0.83);  // 83.8 % at the first zero
  CHECK(b.energy[2] < 0.84);
}

TEST_CASE("M6 case 3: MTF (2 / pi) [acos x - x sqrt(1 - x^2)], paraboloid F/10", "[m6]") {
  // tests/reference/m6/paraboloid_f10.rtt.json, W = 0: tangential and sagittal against the
  // circular MTF at every lattice point with x = nu / nu_c in [0.1, 0.9], nu_c = 1 / (lambda
  // F#_w) from the guarded mapping. Grids 256/512 (erratum #206), samples l of N and 2 l of 2 N
  // share x; bound per curve, relative (as case 2); x taken from the result's frequencies, so a
  // wrong frequency scale fails too.
  const CompiledSystem cs = compiled(load("m6/paraboloid_f10.rtt.json"));
  const double n_u = image_index(cs, PathId{0}) * 10.0 / (2.0 * kF);
  const PupilFunction p = pupil(cs, PathId{0}, 16);
  guard_mapping(p, n_u);
  const double nu_c = 1.0 / (kLambdaMm * f_number(p));
  rtt::analysis::MtfOptions o256;
  o256.pupil.grid = 256;
  rtt::analysis::MtfOptions o512;
  o512.pupil.grid = 512;
  const auto a = rtt::analysis::mtf(cs, PathId{0}, 0, o256);
  const auto b = rtt::analysis::mtf(cs, PathId{0}, 0, o512);
  int checked = 0;
  for (int d = 0; d < 2; ++d) {
    INFO((d == 0 ? "tangential" : "sagittal"));
    const auto& ca = d == 0 ? a.tangential : a.sagittal;
    const auto& cb = d == 0 ? b.tangential : b.sagittal;
    REQUIRE(ca.modulation.size() == 257);
    REQUIRE(ca.frequency.size() == 257);
    REQUIRE(cb.modulation.size() == 513);
    CHECK(std::abs(ca.frequency.back() / nu_c - 1.0) <= 1e-12);
    double e_max = 0.0;
    for (std::size_t l = 0; l <= 256; ++l) {
      const double x = ca.frequency[l] / nu_c;
      if (x < 0.1 || x > 0.9) continue;
      e_max =
          std::max(e_max, 2.0 * std::abs(ca.modulation[l] - cb.modulation[2 * l]) / mtf_circle(x));
    }
    INFO("e_max (relative) " << e_max);
    CHECK(e_max <= 1e-3);
    for (std::size_t l = 0; l <= 256; ++l) {
      const double x = ca.frequency[l] / nu_c;
      if (x < 0.1 || x > 0.9) continue;
      INFO("x = " << x << ": N=512 " << cb.modulation[2 * l] << ", expected " << mtf_circle(x));
      CHECK(std::abs(cb.modulation[2 * l] - mtf_circle(x)) / mtf_circle(x) <= e_max + 1e-12);
      ++checked;
    }
    CHECK(ca.modulation[26] > 0.8);    // content guard: x = 0.10
    CHECK(ca.modulation[230] < 0.05);  // x = 0.90
  }
  CHECK(checked > 400);
}

TEST_CASE("M6 case 4: defocus Strehl sinc^2(w), paraboloid F/50", "[m6]") {
  // tests/reference/m6/paraboloid_f50.rtt.json (D = 2 mm): the wavefront of a focus shift Delta
  // is w rho^2 up to a relative rho^4 term below 6e-5, whose effect on S is at most
  // 4 pi 6e-5 w (derivation in test_psf.cpp, "defocus Strehl"). w from the exact sphere formula
  // at the rim (ADR 0032, point 5): W = n' [Delta (1 - cos a) + sqrt(R'^2 - Delta^2 sin^2 a) - R'],
  // sin a = h / (f + h^2 / (4 f)), R' = |C' - E| with the exit pupil E from first_order and
  // n' of AIR. Expected S = sinc^2(w) (Wyant & Creath, Eq. (64)). Grids 128/256 (as F2).
  const CompiledSystem cs = compiled(load("m6/paraboloid_f50.rtt.json"));
  const double n_image = image_index(cs, PathId{0});
  const double h = 1.0;
  const double sin_a = h / (kF + h * h / (4.0 * kF));
  const double one_minus_cos = sin_a * sin_a / (1.0 + std::sqrt(1.0 - sin_a * sin_a));
  const auto fo = rtt::paraxial::first_order(cs, PathId{0}, 0);
  REQUIRE(fo.exit_pupil.has_value());
  REQUIRE(fo.exit_pupil->z.has_value());
  const double z_exit = *fo.exit_pupil->z;
  REQUIRE(z_exit > -kF);
  std::vector<double> shifts;
  std::vector<double> ws;
  for (const double w_target : {0.1, 0.3, 0.5, 0.8}) {
    const double delta = w_target * kLambdaMm / one_minus_cos;
    const double r = std::abs(-kF - delta - z_exit);
    const double t = delta * delta * sin_a * sin_a;
    const double w_mm = n_image * (delta * one_minus_cos - t / (std::sqrt(r * r - t) + r));
    shifts.push_back(delta);
    ws.push_back(w_mm / kLambdaMm);
  }
  PupilFunctionOptions o128;
  o128.grid = 128;
  PupilFunctionOptions o256;
  o256.grid = 256;
  const auto t128 = rtt::analysis::through_focus(cs, PathId{0}, 0, 0, shifts, o128);
  const auto t256 = rtt::analysis::through_focus(cs, PathId{0}, 0, 0, shifts, o256);
  REQUIRE(t128.strehl.size() == shifts.size());
  REQUIRE(t256.strehl.size() == shifts.size());
  double e_max = 0.0;
  for (std::size_t i = 0; i < shifts.size(); ++i) {
    const double w = ws[i];
    const double e = 2.0 * std::abs(t128.strehl[i] - t256.strehl[i]) + 4.0 * kPi * 6e-5 * w + 1e-12;
    INFO("w = " << w << ": N=128 " << t128.strehl[i] << ", N=256 " << t256.strehl[i] << ", sinc^2 "
                << sinc2(w));
    CHECK(std::abs(t256.strehl[i] - sinc2(w)) <= e);
    CHECK(t256.strehl[i] < 1.0);  // content guard
    e_max = std::max(e_max, e);
  }
  CHECK(e_max <= 1e-3);  // the curve bound itself (ADR 0033, point 10)
}

TEST_CASE("M6 case 5: Zernike c11 and c4 from W040 = S_I / 8, plano-convex lens", "[m6]") {
  // tests/reference/m6/planoconvex_f10.rtt.json: the lens of m1/singlet_const in VACUUM at
  // D = 10 mm (NA 0.05), the image surface in the paraxial focus (checked against first_order),
  // on axis. W = W040 rho^4 + W060 rho^6 + ..., W040 = S_I / 8 (Seidel, rtt::paraxial::seidel);
  // the projection of rho^4, rho^6 and rho^8 on Z_4, Z_11, Z_22 gives (as V2, #230)
  //   c11 = W040 / (6 sqrt5) + sqrt35 c22,  c4 = W040 / (2 sqrt3) + 3 sqrt21 c22,
  // up to the rho^8 share e = W080 (<rho^8, Z_j> - ...), W080 from c37; tolerance 2 |e| plus
  // 1e-6 |W040| for the rounding of S_I and the OPD. Grid 64 (as V2).
  const CompiledSystem cs = compiled(load("m6/planoconvex_f10.rtt.json"));
  const std::uint16_t ref = cs.reference_wavelength();
  const auto fo = rtt::paraxial::first_order(cs, PathId{0}, ref);
  REQUIRE(fo.rear_focal_z.has_value());
  const auto image = cs.find_surface(rtt::model::SurfaceId("IMG"));
  REQUIRE(image.has_value());
  const Vec3 vertex = cs.surfaces()[*image].to_global.apply_point(Vec3::Zero());
  CHECK(std::abs(vertex.z() - *fo.rear_focal_z) <= 1e-9);  // the file holds the paraxial focus
  const double lambda_mm = cs.wavelengths_um()[ref] * 1e-3;
  const double w040 = rtt::paraxial::seidel(cs, PathId{0}, ref).sum.s1 / 8.0 / lambda_mm;
  const auto fit = rtt::analysis::zernike_fit(pupil(cs, PathId{0}, 64));
  const auto& c = fit.coefficients;
  REQUIRE(c.size() >= 37);
  const double c22 = c[21];
  const double w080 = c[36] / projection(4, 37);
  const double e11 = w080 * (projection(4, 11) - std::sqrt(35.0) * projection(4, 22));
  const double e4 = w080 * (projection(4, 4) - 3.0 * std::sqrt(21.0) * projection(4, 22));
  const double expected11 = w040 / (6.0 * std::sqrt(5.0)) + std::sqrt(35.0) * c22;
  const double expected4 = w040 / (2.0 * std::sqrt(3.0)) + 3.0 * std::sqrt(21.0) * c22;
  INFO("W040 " << w040 << " waves, c11 " << c[10] << ", c4 " << c[3] << ", c22 " << c22);
  CHECK(std::abs(c[10] - expected11) <= 2.0 * std::abs(e11) + 1e-6 * std::abs(w040));
  CHECK(std::abs(c[3] - expected4) <= 2.0 * std::abs(e4) + 1e-6 * std::abs(w040));
  CHECK(std::abs(w040) > 0.1);  // content guard: a measurable spherical aberration
}

TEST_CASE("M6 case 6: Jones pupil of a collimated 45 deg fold, diag(r_s, r_p)", "[m6]") {
  // tests/reference/m6/fold45_metal.rtt.json: a collimated beam along +z on a mirror tilted
  // 45 deg about x (normal (0, -sin 45, cos 45)) into +y, the detector facing +y. In the s/p basis
  // J = diag(r_s, r_p) at every point (ADR 0034, points 5 and 9; Fresnel after Byrnes, Eq. (6),
  // at 45 deg, rtt::polar::fresnel_at_angle), D = (|r_s|^2 - |r_p|^2) / (|r_s|^2 + |r_p|^2) and
  // the physical retardance |arg(-r_p / r_s)|. The ideal conductor (no material) gives r_s = -1,
  // r_p = 1, D = delta = 0. Tolerances as the component test: J 1e-9, D and delta 1e-10.
  const auto a = rtt::polar::fresnel_at_angle<double>(Cx(1.0), kMetal, std::cos(kPi / 4.0));
  struct Case {
    const char* name;
    bool metal;
    Cx rs, rp;
  };
  for (const Case& k :
       {Case{"metal", true, a.rs, a.rp}, Case{"ideal conductor", false, -1.0, 1.0}}) {
    INFO(k.name);
    System s = load("m6/fold45_metal.rtt.json");
    if (!k.metal) child(s, 0).material = std::nullopt;
    const CompiledSystem cs = compiled(s);
    const JonesPupil jp = rtt::analysis::jones_pupil(cs, PathId{0}, parallel_rays(), kChief);
    REQUIRE(jp.status.size() == 25);
    REQUIRE(jp.j.size() == 100);
    const double d = (std::norm(k.rs) - std::norm(k.rp)) / (std::norm(k.rs) + std::norm(k.rp));
    const double delta = std::abs(std::arg(-k.rp / k.rs));
    for (std::size_t r = 0; r < 25; ++r) {
      INFO("point " << r);
      REQUIRE(jp.status[r] == RayStatus::Alive);
      CHECK(std::abs(jp.j[4 * r + 0] - k.rs) <= 1e-9);
      CHECK(std::abs(jp.j[4 * r + 1]) <= 1e-9);
      CHECK(std::abs(jp.j[4 * r + 2]) <= 1e-9);
      CHECK(std::abs(jp.j[4 * r + 3] - k.rp) <= 1e-9);
      CHECK(std::abs(jp.diattenuation[r] - d) <= 1e-10);
      CHECK(std::abs(jp.retardance[r] - delta) <= 1e-10);
    }
    if (k.metal) {
      CHECK(d > 1e-2);      // content guard: a real diattenuation
      CHECK(delta > 1e-2);  // and a real physical retardance
    }
  }
}

TEST_CASE("M6 case 7: Jones pupil of a converging 45 deg fold, per point against Fresnel", "[m6]") {
  // tests/reference/m6/fold45_focus.rtt.json: an ideal lens f = 50 (P = Q = a rotation without
  // retardance or diattenuation, ADR 0031), the metal fold at z = 20, the detector 2 mm behind
  // the focus. Per point the basis-free values: D and the physical retardance are those of one
  // Fresnel reflection at the local angle of incidence cos(theta) = |k1 . N| (k1 the direction
  // after the lens, N the mirror normal), and W = (OPL_c - OPL) / lambda from the recorded OPL
  // (ADR 0034, point 6). The full J per point is the component test (test_jones_pupil.cpp,
  // "converging beam over a 45 deg fold"). Tolerances: D and delta 1e-10, W 1e-9 waves.
  const CompiledSystem cs = compiled(load("m6/fold45_focus.rtt.json"));
  std::size_t chief = 0;
  const RayBatch start = converging_rays(chief);
  const JonesPupil jp = rtt::analysis::jones_pupil(cs, PathId{0}, start, chief);
  RayBatch traced = start;
  rtt::trace::RayPaths paths;
  [[maybe_unused]] const auto stats =
      rtt::trace::SequentialTracer().trace(cs, PathId{0}, traced, paths);
  const Vec3 n(0.0, -std::sin(kPi / 4.0), std::cos(kPi / 4.0));
  const double lambda_mm = 0.55e-3;
  const double opl_c = paths.opl[chief * paths.slots + 3];  // slots: start, IL, M, IMG
  REQUIRE(jp.status.size() == start.size());
  double w_max = 0.0;
  double d_min = 1.0;
  double d_max = 0.0;
  for (std::size_t r = 0; r < start.size(); ++r) {
    INFO("point " << r);
    REQUIRE(jp.status[r] == RayStatus::Alive);
    const Vec3 k1 = paths.direction_at(r, 1);
    const auto f = rtt::polar::fresnel_at_angle<double>(Cx(1.0), kMetal, std::abs(k1.dot(n)));
    const double d = (std::norm(f.rs) - std::norm(f.rp)) / (std::norm(f.rs) + std::norm(f.rp));
    const double delta = std::abs(std::arg(-f.rp / f.rs));
    CHECK(std::abs(jp.diattenuation[r] - d) <= 1e-10);
    CHECK(std::abs(jp.retardance[r] - delta) <= 1e-10);
    const double w = (opl_c - paths.opl[r * paths.slots + 3]) / lambda_mm;
    CHECK(std::abs(jp.w[r] - w) <= 1e-9);
    w_max = std::max(w_max, std::abs(w));
    d_min = std::min(d_min, d);
    d_max = std::max(d_max, d);
  }
  CHECK(w_max > 1.0);           // content guard: the plane OPD of a converging wave
  CHECK(d_max - d_min > 1e-4);  // and the angle of incidence varies over the pupil
}

TEST_CASE("M6 case 8: ideal lens, Airy zero and Strehl 1", "[m6]") {
  // tests/reference/r2/ideal_lens.rtt.json, path "ideal lens": f = 50 mm, object at 75 mm, image
  // at 150 mm, stop of radius 10 at z = -1 (the entrance pupil); for the design conjugate W = 0
  // exactly (ADR 0031), so the PSF is the Airy pattern of the linear model and S = 1. Guard:
  // the paraxial marginal ray from the axial object point at z = -75 through the stop rim meets
  // the lens at h = 10 * 75 / 74 and the image at z = 150, n'|u'| = n' h / 150 (AIR).
  const CompiledSystem cs = compiled(load("r2/ideal_lens.rtt.json"));
  const auto id = cs.find_path("ideal lens");
  REQUIRE(id.has_value());
  const double h = 10.0 * 75.0 / 74.0;
  check_airy_zero(cs, *id, image_index(cs, *id) * h / 150.0);
  PsfOptions o;
  o.pupil.grid = 64;
  const auto psf = rtt::analysis::psf(cs, *id, 0, o);
  CHECK(std::abs(psf.strehl - 1.0) <= 1e-12);
}

TEST_CASE("M6: bitwise the same at 1, 4 and all threads", "[m6]") {
  // Determinism (rule 7, ADR 0033, point 12): the PSF, the through-focus Strehl and the Jones
  // pupil of the acceptance systems, at the default thread count and with 1 and 4 threads.
  const CompiledSystem f10 = compiled(load("m6/paraboloid_f10.rtt.json"));
  const CompiledSystem fold = compiled(load("m6/fold45_focus.rtt.json"));
  std::size_t chief = 0;
  const RayBatch start = converging_rays(chief);
  PsfOptions o;
  o.pupil.grid = 32;
  const std::vector<double> shifts{0.0, 0.05};
  const auto psf_all = rtt::analysis::psf(f10, PathId{0}, 0, o);
  const auto tf_all = rtt::analysis::through_focus(f10, PathId{0}, 0, 0, shifts, o.pupil);
  const JonesPupil jp_all = rtt::analysis::jones_pupil(fold, PathId{0}, start, chief);
  REQUIRE(psf_all.intensity.size() > 1);
  REQUIRE(jp_all.j.size() == 4 * start.size());
  for (const int threads : {1, 4}) {
    INFO(threads << " threads");
    const tbb::global_control limit(tbb::global_control::max_allowed_parallelism,
                                    static_cast<std::size_t>(threads));
    CHECK(same_bits(rtt::analysis::psf(f10, PathId{0}, 0, o).intensity, psf_all.intensity));
    CHECK(same_bits(rtt::analysis::through_focus(f10, PathId{0}, 0, 0, shifts, o.pupil).strehl,
                    tf_all.strehl));
    const JonesPupil jp = rtt::analysis::jones_pupil(fold, PathId{0}, start, chief);
    CHECK(same_bits(jp.j, jp_all.j));
    CHECK(same_bits(jp.diattenuation, jp_all.diattenuation));
    CHECK(same_bits(jp.retardance, jp_all.retardance));
  }
}
