// FFT MTF and encircled energy (#206; ADR 0033, points 5-7, 10, 13, with the erratum on the
// signs). Reference system: the paraboloid with stop (tests/reference/m2/paraboloid_stop, f = 100
// mm, object at infinity) at F/10 (EPD 10 mm), as in test_psf.cpp; on axis it images without
// aberration, so W = 0 and MTF and EE are those of the circular pupil.
//
// Expected values (derived in ADR 0033, docs/quellen.md):
// - MTF(x) = (2/pi)[acos x - x sqrt(1 - x^2)], x = nu / nu_c, nu_c = 1 / (lambda F#_w)
//   (Shaw, EELE 582 "OTF/MTF", slides 10 and 11, with 2/pi instead of the misprinted 2/rho).
// - EE(v) = 1 - J0(v)^2 - J1(v)^2, v = pi r / (lambda F#_w) (from DLMF (10.6.2), (10.6.3),
//   (10.22.1)); first zero r0 = j11 lambda F#_w / pi, j11 by Newton on std::cyl_bessel_j.
// - F#_w = 1 / (2 |M_11|) from the mapping; the guard that ties M to the geometry and to
//   first_order is in test_psf.cpp (ADR 0033, point 10).
// Tolerances of ADR 0033, point 10: X at N and 2N, compare X_2N, bound e = 2 |X_N - X_2N|, and e
// itself below the architecture's 1e-3.
//
// The exact cases inject a known wavefront into a pupil function and use the overloads on pupil
// level: a cylindrical defocus w py^2 (tangential and sagittal not swapped) and a tilt (the
// centre EeCentre::LinearPhase and the signs of the erratum; psf_at() finds the moved peak
// independently).

#include <oneapi/tbb/task_arena.h>

#include <algorithm>
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
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

using rtt::analysis::EeCentre;
using rtt::analysis::EncircledEnergy;
using rtt::analysis::EncircledEnergyOptions;
using rtt::analysis::Mtf;
using rtt::analysis::MtfCurve;
using rtt::analysis::MtfOptions;
using rtt::analysis::MtfOutputs;
using rtt::analysis::Point2;
using rtt::analysis::PupilFunction;
using rtt::analysis::PupilFunctionOptions;
using rtt::compile::CompiledSystem;
using rtt::compile::PathId;

namespace {

constexpr double kLambdaMm = 0.5876e-3;
constexpr double kF = 100.0;  // focal length of the paraboloid, R / 2

rtt::model::System paraboloid_model(double epd) {
  rtt::model::System s =
      rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m2/paraboloid_stop.rtt.json");
  s.aperture.value = rtt::model::Param(epd);
  return s;
}

CompiledSystem compiled(const rtt::model::System& s) {
  return rtt::compile::compile(s, rtt::material::MaterialLibrary{});
}

CompiledSystem paraboloid(double epd) {
  return compiled(paraboloid_model(epd));
}

PupilFunction pupil(const CompiledSystem& cs, int n, std::uint16_t field = 0) {
  PupilFunctionOptions o;
  o.grid = n;
  return rtt::analysis::pupil_function(cs, PathId{0}, field, 0, o);
}

/// n' = |Re n| of the image medium at `wavelength`, from the compiled media (not through the code
/// under test). The paraboloid lies in AIR: n' = 1.00027 at 0.5876 um, carried by M.
double image_index(const CompiledSystem& cs, std::uint16_t wavelength) {
  const auto& last = cs.path(PathId{0}).events.back();
  return std::abs(cs.media()[last.medium_after].index[wavelength].real());
}

/// lambda F#_w in mm, F#_w = 1 / (2 |M_11|) of the mapping (M includes n').
double airy_unit(const PupilFunction& p) {
  REQUIRE(p.mapping.has_value());
  return kLambdaMm / (2.0 * std::abs(p.mapping->m[0]));
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

double ee_airy(double v) {
  const double j0 = std::cyl_bessel_j(0.0, v);
  const double j1 = std::cyl_bessel_j(1.0, v);
  return 1.0 - j0 * j0 - j1 * j1;
}

double mtf_circle(double x) {
  return 2.0 / std::numbers::pi * (std::acos(x) - x * std::sqrt(1.0 - x * x));
}

Mtf mtf_n(const CompiledSystem& cs, int n, std::uint16_t field = 0) {
  MtfOptions o;
  o.pupil.grid = n;
  return rtt::analysis::mtf(cs, PathId{0}, field, o);
}

EncircledEnergy ee_n(const CompiledSystem& cs, int n, std::span<const double> radii) {
  EncircledEnergyOptions o;
  o.pupil.grid = n;
  return rtt::analysis::encircled_energy(cs, PathId{0}, 0, radii, o);
}

bool same_bits(double a, double b) {
  return std::memcmp(&a, &b, sizeof(double)) == 0;
}

/// The curves exist with the samples l = 0 ... n before they are indexed, so that the tests fail
/// cleanly against a stub instead of reading out of range.
void require_curves(const Mtf& m, int n) {
  const auto size = static_cast<std::size_t>(n) + 1;
  REQUIRE(m.tangential.frequency.size() == size);
  REQUIRE(m.tangential.modulation.size() == size);
  REQUIRE(m.sagittal.frequency.size() == size);
  REQUIRE(m.sagittal.modulation.size() == size);
}

/// The energies exist for every radius before they are indexed.
void require_energy(const EncircledEnergy& e, std::size_t radii) {
  REQUIRE(e.energy.size() == radii);
}

/// Cell value u_j = f A exp(-2 pi i W lambda_ref / lambda) (ADR 0032), W in waves of lambda_ref.
std::complex<double> cell_value(const PupilFunction& p, std::size_t j) {
  const double phase =
      -2.0 * std::numbers::pi * p.w[j] * p.reference_wavelength_um / p.wavelength_um;
  return p.fraction[j] * p.amplitude[j] * std::polar(1.0, phase);
}

/// Adds a tilt to W so that the PSF moves by t (mm in e_x, e_y): u_j gets the factor
/// exp(+2 pi i t . M p_j / lambda), since U(x - t) = sum_j u_j exp(2 pi i t . M p_j / lambda)
/// exp(-2 pi i x . M p_j / lambda) with the kernel of ADR 0033, point 2.
PupilFunction tilted(const PupilFunction& p, const Point2& t) {
  REQUIRE(p.mapping.has_value());
  const auto& m = p.mapping->m;
  PupilFunction q = p;
  for (std::size_t j = 0; j < q.w.size(); ++j) {
    const double sx = m[0] * q.px[j] + m[1] * q.py[j];
    const double sy = m[2] * q.px[j] + m[3] * q.py[j];
    // -2 pi i W' lambda_ref / lambda = -2 pi i W lambda_ref / lambda + 2 pi i t . s / lambda
    q.w[j] -= (t.x * sx + t.y * sy) / (p.reference_wavelength_um * 1e-3);
  }
  return q;
}

bool has_code(const std::vector<rtt::model::Diagnostic>& warnings, const char* code) {
  return std::any_of(warnings.begin(), warnings.end(),
                     [code](const auto& d) { return d.code == code; });
}

}  // namespace

TEST_CASE("mtf: diffraction-limited MTF of the paraboloid at N = 128 and 256", "[mtf]") {
  // Samples l of N and 2 l of 2 N have the same x = l / N; x is taken from the result's
  // frequency, so a wrong frequency scale fails too.
  const CompiledSystem cs = paraboloid(10.0);
  const Mtf a = mtf_n(cs, 128);
  const Mtf b = mtf_n(cs, 256);
  const double nu_c = 1.0 / airy_unit(pupil(cs, 16));
  int checked = 0;
  for (const MtfCurve* ca : {&a.tangential, &a.sagittal}) {
    const MtfCurve& cb = ca == &a.tangential ? b.tangential : b.sagittal;
    REQUIRE(ca->frequency.size() == 129);
    REQUIRE(ca->modulation.size() == 129);
    REQUIRE(cb.frequency.size() == 257);
    REQUIRE(cb.modulation.size() == 257);
    CHECK(same_bits(ca->modulation.front(), 1.0));  // C(0) := E
    CHECK(ca->modulation.back() == 0.0);            // no overlap at l = N
    CHECK(std::abs(ca->frequency.back() / nu_c - 1.0) <= 1e-12);
    CHECK(std::is_sorted(ca->frequency.begin(), ca->frequency.end()));
    for (std::size_t l = 0; l <= 128; ++l) {
      const double x = ca->frequency[l] / nu_c;
      CHECK(std::abs(x - static_cast<double>(l) / 128.0) <= 1e-12);
      CHECK(std::abs(cb.frequency[2 * l] - ca->frequency[l]) <= 1e-12 * ca->frequency.back());
      if (x < 0.1 || x > 0.9) continue;
      const double expected = mtf_circle(x);
      const double xa = ca->modulation[l];
      const double xb = cb.modulation[2 * l];
      const double e = 2.0 * std::abs(xa - xb) + 1e-12;
      INFO("x = " << x << ": N=128 " << xa << ", N=256 " << xb << ", expected " << expected);
      CHECK(e <= 1e-3 * expected);
      CHECK(std::abs(xb - expected) <= e);
      ++checked;
    }
    // Content guard: the curve falls from near 1 to near 0 over the checked range.
    CHECK(ca->modulation[13] > 0.8);    // x = 0.10
    CHECK(ca->modulation[115] < 0.05);  // x = 0.90
  }
  CHECK(checked > 200);
}

TEST_CASE("mtf: tangential equals sagittal on the axis, along e_y and e_x", "[mtf]") {
  // Not bit for bit: the rays of (px, py) and (py, px) are traced with the coordinates in
  // different roles, so W and M agree up to rounding only.
  const Mtf m = mtf_n(paraboloid(10.0), 64);
  REQUIRE(m.tangential.modulation.size() == 65);
  REQUIRE(m.sagittal.modulation.size() == 65);
  for (std::size_t l = 0; l <= 64; ++l) {
    INFO("l = " << l);
    CHECK(std::abs(m.tangential.modulation[l] - m.sagittal.modulation[l]) <= 1e-12);
    CHECK(std::abs(m.tangential.frequency[l] - m.sagittal.frequency[l]) <=
          1e-12 * m.sagittal.frequency.back());
  }
  CHECK(std::abs(std::abs(m.tangential.direction[1]) - 1.0) <= 1e-12);
  CHECK(std::abs(m.tangential.direction[0]) <= 1e-12);
  CHECK(std::abs(std::abs(m.sagittal.direction[0]) - 1.0) <= 1e-12);
  CHECK(std::abs(m.sagittal.direction[1]) <= 1e-12);
  CHECK(m.tangential.otf.empty());  // only on request
  CHECK(m.map.empty());
}

TEST_CASE("mtf: a cylindrical defocus w py^2 lowers only the tangential curve", "[mtf]") {
  // For shifts along px the phases of u_(j+l) and u_j cancel (same py), so the sagittal curve is
  // that of the undisturbed pupil; the tangential one is checked against the autocorrelation
  // summed here directly, without FFT. Swapped axes in the code fail both checks.
  const int n = 64;
  const PupilFunction p = pupil(paraboloid(10.0), n);
  PupilFunction q = p;
  for (std::size_t j = 0; j < q.w.size(); ++j) q.w[j] += 0.5 * q.py[j] * q.py[j];
  const Mtf m0 = rtt::analysis::mtf(p);
  const Mtf m1 = rtt::analysis::mtf(q);
  require_curves(m0, n);
  require_curves(m1, n);
  for (std::size_t l = 0; l <= 64; ++l) {
    CHECK(std::abs(m1.sagittal.modulation[l] - m0.sagittal.modulation[l]) <= 1e-13);
  }
  double energy = 0.0;
  for (std::size_t j = 0; j < q.fraction.size(); ++j) {
    energy += q.fraction[j] * q.amplitude[j] * q.amplitude[j];
  }
  const auto size = static_cast<std::size_t>(n);
  for (const std::size_t l : {1U, 5U, 20U, 40U}) {
    std::complex<double> c = 0.0;
    for (std::size_t k = 0; k + l < size; ++k) {  // rows (py)
      for (std::size_t i = 0; i < size; ++i) {    // columns (px)
        c += cell_value(q, (k + l) * size + i) * std::conj(cell_value(q, k * size + i));
      }
    }
    INFO("l = " << l);
    CHECK(std::abs(m1.tangential.modulation[l] - std::abs(c) / energy) <= 1e-12);
  }
  // Content guard: x = 19 / 64 = 0.30.
  CHECK(m1.tangential.modulation[19] < m1.sagittal.modulation[19] - 0.05);
}

TEST_CASE("mtf: fields +y and -y give mirrored curves", "[mtf]") {
  // 2 degrees: the tangential coma of the paraboloid is about 3 theta / (16 F#^2) f = 6.5 um,
  // close to the Airy radius 1.22 lambda F# = 7.2 um, so tangential and sagittal differ visibly
  // (estimate; the content guard below checks it).
  rtt::model::System s = paraboloid_model(10.0);
  s.fields.points = {{0.0, 2.0, 1.0}, {0.0, -2.0, 1.0}};
  const CompiledSystem cs = compiled(s);
  const Mtf up = mtf_n(cs, 64, 0);
  const Mtf down = mtf_n(cs, 64, 1);
  double largest = 0.0;
  for (const auto curve : {&Mtf::tangential, &Mtf::sagittal}) {
    const MtfCurve& a = up.*curve;
    const MtfCurve& b = down.*curve;
    REQUIRE(a.modulation.size() == 65);
    REQUIRE(b.modulation.size() == 65);
    for (std::size_t l = 0; l <= 64; ++l) {
      INFO("l = " << l);
      CHECK(std::abs(a.modulation[l] - b.modulation[l]) <= 1e-12);
      CHECK(std::abs(a.frequency[l] - b.frequency[l]) <= 1e-12 * a.frequency.back());
    }
  }
  for (std::size_t l = 0; l <= 64; ++l) {
    largest = std::max(largest, std::abs(up.tangential.modulation[l] - up.sagittal.modulation[l]));
  }
  // Tangential in the meridional plane (e_x is the global x axis for a field along y).
  for (const Mtf* m : {&up, &down}) {
    CHECK(std::abs(std::abs(m->tangential.direction[1]) - 1.0) <= 1e-12);
    CHECK(std::abs(std::abs(m->sagittal.direction[0]) - 1.0) <= 1e-12);
    CHECK(std::abs(std::abs(m->e_x.x()) - 1.0) <= 1e-12);
  }
  CHECK(largest > 0.02);  // content guard: off the axis the two curves differ
}

TEST_CASE("encircled energy: 1 - J0^2 - J1^2 of the circular pupil at N = 64 and 128", "[ee]") {
  const CompiledSystem cs = paraboloid(10.0);
  const double s = airy_unit(pupil(cs, 16));
  const double r0 = j11() / std::numbers::pi;
  const std::vector<double> units{0.5, 1.0, r0, 2.0, 3.0, 5.0};
  std::vector<double> radii;
  for (const double u : units) radii.push_back(u * s);
  const EncircledEnergy a = ee_n(cs, 64, radii);
  const EncircledEnergy b = ee_n(cs, 128, radii);
  REQUIRE(a.energy.size() == radii.size());
  REQUIRE(b.energy.size() == radii.size());
  CHECK(a.radii == radii);
  for (std::size_t i = 0; i < radii.size(); ++i) {
    const double expected = ee_airy(std::numbers::pi * units[i]);
    const double e = 2.0 * std::abs(a.energy[i] - b.energy[i]) + 1e-12;
    INFO("r = " << units[i] << " lambda F#: N=64 " << a.energy[i] << ", N=128 " << b.energy[i]
                << ", expected " << expected);
    CHECK(e <= 1e-3 * expected);
    CHECK(std::abs(b.energy[i] - expected) <= e);
  }
  // Content guard: 83.8 % at the first zero.
  CHECK(b.energy[2] > 0.83);
  CHECK(b.energy[2] < 0.84);
  CHECK(b.offset.x == 0.0);
  CHECK(b.offset.y == 0.0);
}

TEST_CASE("encircled energy: up to radius_max, a quarter of the period", "[ee]") {
  // The periodic PSF of the discrete model has the period N lambda F#_w on the axis; within
  // r_max = N lambda F#_w / 4 the neighbouring periods add only their Airy tails, about
  // 4 r^2 / (pi^2 N^3) (ADR 0033, point 6): at r = 16 lambda F#_w that is 4e-4 for N = 64 and
  // 5e-5 for N = 128, falling like N^-3, so the rule e = 2 |X_N - X_2N| covers it; e stays below
  // 2e-3 (twice the tail plus the discretisation of point 9).
  const CompiledSystem cs = paraboloid(10.0);
  const double s = airy_unit(pupil(cs, 16));
  const std::array<double, 0> none{};
  const EncircledEnergy probe = ee_n(cs, 64, none);
  CHECK(probe.energy.empty());  // empty radii: an empty result, no error
  CHECK(std::abs(probe.radius_max / (16.0 * s) - 1.0) <= 1e-12);
  REQUIRE(probe.edge_energy_deficit > 0.0);
  const std::array<double, 1> at_max{probe.radius_max};
  const EncircledEnergy full = ee_n(cs, 64, at_max);
  require_energy(full, 1);
  CHECK(full.energy[0] <= 1.0 - full.edge_energy_deficit + 1e-12);
  const std::array<double, 1> r16{16.0 * s};
  const EncircledEnergy at64 = ee_n(cs, 64, r16);
  const EncircledEnergy at128 = ee_n(cs, 128, r16);
  require_energy(at64, 1);
  require_energy(at128, 1);
  const double a = at64.energy[0];
  const double b = at128.energy[0];
  const double expected = ee_airy(16.0 * std::numbers::pi);
  const double e = 2.0 * std::abs(a - b) + 1e-12;
  INFO("N=64 " << a << ", N=128 " << b << ", expected " << expected);
  CHECK(e <= 2e-3);
  CHECK(std::abs(b - expected) <= e);
  const std::array<double, 1> beyond{
      std::nextafter(probe.radius_max, std::numeric_limits<double>::infinity())};
  CHECK_THROWS_AS(ee_n(cs, 64, beyond), std::invalid_argument);
}

TEST_CASE("encircled energy: LinearPhase centre of a tilted wavefront", "[ee]") {
  // A tilt moves the PSF by t exactly (ADR 0033, point 2). EeCentre::LinearPhase must find t, the
  // circle about it must hold the energy of the untilted PSF about C', and psf_at() - which
  // shares no code with the lag sum - must find the peak there. This fixes both signs of the
  // erratum: a wrong sign in the centre or in exp(-2 pi i nu . x0) fails one of the checks.
  const PupilFunction p = pupil(paraboloid(10.0), 64);
  const double s = airy_unit(p);
  const Point2 t{1.5 * s, -0.7 * s};
  const PupilFunction q = tilted(p, t);
  const double r0 = j11() / std::numbers::pi * s;
  const std::vector<double> radii{0.5 * s, r0, 3.0 * s};
  const EncircledEnergy plain = rtt::analysis::encircled_energy(p, radii);
  const EncircledEnergy moved = rtt::analysis::encircled_energy(q, radii, EeCentre::LinearPhase);
  const EncircledEnergy fixed = rtt::analysis::encircled_energy(q, radii);
  for (const EncircledEnergy* e : {&plain, &moved, &fixed}) require_energy(*e, radii.size());
  const double norm_t = std::hypot(t.x, t.y);
  CHECK(std::abs(moved.offset.x - t.x) <= 1e-12 * norm_t);
  CHECK(std::abs(moved.offset.y - t.y) <= 1e-12 * norm_t);
  REQUIRE(q.mapping.has_value());
  const rtt::math::Vec3 expected_centre =
      q.mapping->centre + moved.offset.x * moved.e_x + moved.offset.y * moved.e_y;
  CHECK((moved.centre - expected_centre).norm() <= 1e-12 * norm_t);
  const std::array<Point2, 2> at{moved.offset, Point2{0.0, 0.0}};
  const std::vector<double> strehl = rtt::analysis::psf_at(q, at);
  CHECK(std::abs(strehl[0] - 1.0) <= 1e-12);  // the peak of the moved Airy pattern
  CHECK(strehl[1] < 0.5);
  for (std::size_t i = 0; i < radii.size(); ++i) {
    INFO("r = " << radii[i] / s << " lambda F#");
    CHECK(std::abs(moved.energy[i] - plain.energy[i]) <= 1e-12);
  }
  // Content guard: about C' the moved PSF loses most of the energy of the first ring.
  CHECK(fixed.energy[1] < 0.5);
  CHECK(plain.energy[1] > 0.83);
}

TEST_CASE("mtf: complex OTF, phase of a moved PSF and the 2D map", "[mtf]") {
  // OTF(nu) = integral I exp(-2 pi i nu . x) / integral I: a PSF moved by t gets the factor
  // exp(-2 pi i nu . t) (erratum to ADR 0033, point 5).
  const PupilFunction p = pupil(paraboloid(10.0), 32);
  const double s = airy_unit(p);
  const Point2 t{1.5 * s, -0.7 * s};
  const PupilFunction q = tilted(p, t);
  MtfOutputs out;
  out.complex_otf = true;
  const Mtf m0 = rtt::analysis::mtf(p, out);
  const Mtf m1 = rtt::analysis::mtf(q, out);
  const auto& f = m1.frequency_matrix;
  for (const Mtf* m : {&m0, &m1}) {
    require_curves(*m, 32);
    REQUIRE(m->tangential.otf.size() == 33);
    REQUIRE(m->sagittal.otf.size() == 33);
  }
  double largest_phase = 0.0;
  for (std::size_t l = 0; l <= 32; ++l) {
    const auto ld = static_cast<double>(l);
    // sagittal l = (l, 0), tangential l = (0, l); nu = F l
    const std::array<std::array<double, 2>, 2> nu{{{f[0] * ld, f[2] * ld}, {f[1] * ld, f[3] * ld}}};
    const std::array<const MtfCurve*, 2> c0{&m0.sagittal, &m0.tangential};
    const std::array<const MtfCurve*, 2> c1{&m1.sagittal, &m1.tangential};
    for (std::size_t d = 0; d < 2; ++d) {
      const double arg = -2.0 * std::numbers::pi * (nu[d][0] * t.x + nu[d][1] * t.y);
      const std::complex<double> expected = c0[d]->otf[l] * std::polar(1.0, arg);
      INFO("l = " << l << ", " << (d == 0 ? "sagittal" : "tangential"));
      CHECK(std::abs(c1[d]->otf[l] - expected) <= 1e-12);
      CHECK(same_bits(std::abs(c1[d]->otf[l]), c1[d]->modulation[l]));
      largest_phase = std::max(largest_phase, std::abs(std::arg(c1[d]->otf[l])));
    }
  }
  CHECK(largest_phase > 0.5);  // content guard: the phase transfer is not trivial
  CHECK(m1.sagittal.otf[0] == std::complex<double>(1.0, 0.0));
  // 2D map: the rows and columns through l = 0 are the curves.
  out.map = true;
  const Mtf mm = rtt::analysis::mtf(q, out);
  require_curves(mm, 32);
  REQUIRE(mm.sagittal.otf.size() == 33);
  const std::size_t side = 2 * 32 - 1;
  REQUIRE(mm.map.size() == side * side);
  REQUIRE(mm.otf_map.size() == side * side);
  const auto at = [side](int l1, int l2) {
    return static_cast<std::size_t>(l2 + 31) * side + static_cast<std::size_t>(l1 + 31);
  };
  CHECK(mm.map[at(0, 0)] == 1.0);
  for (int l = 0; l < 32; ++l) {
    const auto lu = static_cast<std::size_t>(l);
    CHECK(same_bits(mm.map[at(l, 0)], mm.sagittal.modulation[lu]));
    CHECK(same_bits(mm.map[at(0, l)], mm.tangential.modulation[lu]));
    CHECK(mm.otf_map[at(l, 0)] == mm.sagittal.otf[lu]);
    CHECK(std::abs(mm.map[at(-l, -l)] - mm.map[at(l, l)]) <= 1e-14);  // |C(-l)| = |C(l)|
  }
}

TEST_CASE("mtf, encircled energy: two wavelengths combine with weights w_i E_i", "[mtf][ee]") {
  // lambda_2 = 1.25 lambda_1, weights 1 : 3. EE: the weighted sum at the same radius. MTF: the
  // complex OTF of lambda_2, interpolated linearly in |nu| onto the reference frequencies (0
  // above its cut-off), summed with the OTF of lambda_1 (ADR 0033, point 7); the interpolation is
  // repeated here.
  rtt::model::System sys = paraboloid_model(10.0);
  sys.wavelengths = {{0.5876, 1.0, true}, {0.5876 * 1.25, 3.0, false}};
  const CompiledSystem cs = compiled(sys);
  const double s = airy_unit(pupil(cs, 16));
  MtfOptions both;
  both.pupil.grid = 32;
  both.outputs.complex_otf = true;
  MtfOptions first = both;
  first.wavelengths = std::vector<std::uint16_t>{0};
  MtfOptions second = both;
  second.wavelengths = std::vector<std::uint16_t>{1};
  const Mtf poly = rtt::analysis::mtf(cs, PathId{0}, 0, both);
  const Mtf a = rtt::analysis::mtf(cs, PathId{0}, 0, first);
  const Mtf b = rtt::analysis::mtf(cs, PathId{0}, 0, second);
  for (const Mtf* m : {&poly, &a, &b}) require_curves(*m, 32);
  REQUIRE(a.tangential.otf.size() == 33);
  REQUIRE(b.tangential.otf.size() == 33);
  REQUIRE(a.wavelengths.size() == 1);
  REQUIRE(b.wavelengths.size() == 1);
  REQUIRE(poly.wavelengths.size() == 2);
  const double wa = 0.25 * a.wavelengths[0].energy;
  const double wb = 0.75 * b.wavelengths[0].energy;
  const auto interp = [](const MtfCurve& c, double nu) -> std::complex<double> {
    if (nu >= c.frequency.back()) return 0.0;
    const auto hi = std::upper_bound(c.frequency.begin(), c.frequency.end(), nu);
    const auto k = static_cast<std::size_t>(hi - c.frequency.begin());
    const double g = (nu - c.frequency[k - 1]) / (c.frequency[k] - c.frequency[k - 1]);
    return (1.0 - g) * c.otf[k - 1] + g * c.otf[k];
  };
  // The second curve ends at the cut-off of lambda_2. nu_c,i = 2 |M_i| / lambda_i with
  // M_i = (n'_i / n'_1) M_1 (the mirror reflects the same rays, the direction sines carry the
  // index of AIR; as in test_psf.cpp), so nu_c,2 / nu_c,1 = (n'_2 / n'_1) / 1.25.
  const double n_ratio = image_index(cs, 1) / image_index(cs, 0);
  CHECK(std::abs(b.tangential.frequency.back() * 1.25 / (n_ratio * a.tangential.frequency.back()) -
                 1.0) <= 1e-12);
  double largest = 0.0;
  for (std::size_t l = 0; l <= 32; ++l) {
    const double nu = a.tangential.frequency[l];
    CHECK(same_bits(poly.tangential.frequency[l], nu));
    const std::complex<double> otf =
        (wa * a.tangential.otf[l] + wb * interp(b.tangential, nu)) / (wa + wb);
    INFO("l = " << l);
    CHECK(std::abs(poly.tangential.modulation[l] - std::abs(otf)) <= 1e-14);
    largest =
        std::max(largest, std::abs(poly.tangential.modulation[l] - a.tangential.modulation[l]));
  }
  CHECK(largest > 0.01);  // content guard: lambda_2 changes the curve
  EncircledEnergyOptions eo;
  eo.pupil.grid = 32;
  const std::vector<double> radii{1.0 * s, 2.0 * s, 3.0 * s};
  const EncircledEnergy ep = rtt::analysis::encircled_energy(cs, PathId{0}, 0, radii, eo);
  eo.wavelengths = std::vector<std::uint16_t>{0};
  const EncircledEnergy ea = rtt::analysis::encircled_energy(cs, PathId{0}, 0, radii, eo);
  eo.wavelengths = std::vector<std::uint16_t>{1};
  const EncircledEnergy eb = rtt::analysis::encircled_energy(cs, PathId{0}, 0, radii, eo);
  for (const EncircledEnergy* e : {&ep, &ea, &eb}) require_energy(*e, radii.size());
  REQUIRE(ea.wavelengths.size() == 1);
  REQUIRE(eb.wavelengths.size() == 1);
  const double va = 0.25 * ea.wavelengths[0].energy;
  const double vb = 0.75 * eb.wavelengths[0].energy;
  for (std::size_t i = 0; i < radii.size(); ++i) {
    CHECK(std::abs(ep.energy[i] - (va * ea.energy[i] + vb * eb.energy[i]) / (va + vb)) <= 1e-15);
    CHECK(std::abs(ea.energy[i] - eb.energy[i]) > 0.01);  // content guard
  }
  CHECK(same_bits(ep.radius_max, std::min(ea.radius_max, eb.radius_max)));
}

TEST_CASE("mtf, encircled energy: bitwise the same for 1, 4 and all threads", "[mtf][ee]") {
  const CompiledSystem cs = paraboloid(10.0);
  const double s = airy_unit(pupil(cs, 16));
  MtfOptions mo;
  mo.pupil.grid = 64;
  mo.pupil.focus_shift = 0.5;
  EncircledEnergyOptions eo;
  eo.pupil = mo.pupil;
  eo.centre = EeCentre::LinearPhase;
  const std::vector<double> radii{1.0 * s, 2.0 * s, 5.0 * s};
  const Mtf m_all = rtt::analysis::mtf(cs, PathId{0}, 0, mo);
  const EncircledEnergy e_all = rtt::analysis::encircled_energy(cs, PathId{0}, 0, radii, eo);
  require_curves(m_all, 64);  // content guard: there are curves to compare
  require_energy(e_all, radii.size());
  for (const int threads : {1, 4}) {
    oneapi::tbb::task_arena arena(threads);
    const Mtf m = arena.execute([&] { return rtt::analysis::mtf(cs, PathId{0}, 0, mo); });
    const EncircledEnergy e =
        arena.execute([&] { return rtt::analysis::encircled_energy(cs, PathId{0}, 0, radii, eo); });
    INFO("threads " << threads);
    REQUIRE(m.tangential.modulation.size() == m_all.tangential.modulation.size());
    CHECK(std::memcmp(m.tangential.modulation.data(), m_all.tangential.modulation.data(),
                      m_all.tangential.modulation.size() * sizeof(double)) == 0);
    CHECK(std::memcmp(m.sagittal.modulation.data(), m_all.sagittal.modulation.data(),
                      m_all.sagittal.modulation.size() * sizeof(double)) == 0);
    REQUIRE(e.energy.size() == e_all.energy.size());
    CHECK(std::memcmp(e.energy.data(), e_all.energy.data(), e_all.energy.size() * sizeof(double)) ==
          0);
    CHECK(same_bits(e.offset.x, e_all.offset.x));
    CHECK(same_bits(e.offset.y, e_all.offset.y));
  }
  // Content guard: the focus shift lowers the MTF below the diffraction limit.
  CHECK(m_all.tangential.modulation[16] < mtf_circle(0.25) - 0.05);
}

TEST_CASE("mtf, encircled energy: warning psf.pupil_undersampled, never psf.aliasing",
          "[mtf][ee]") {
  const CompiledSystem cs = paraboloid(10.0);
  MtfOptions ideal;
  ideal.pupil.grid = 32;
  const Mtf m0 = rtt::analysis::mtf(cs, PathId{0}, 0, ideal);
  CHECK_FALSE(has_code(m0.warnings, "psf.pupil_undersampled"));
  // Strong defocus at N = 32 as in test_psf.cpp: w = 20 waves, steps of about 2.5 waves.
  const double h = 5.0;
  const double sin_a = h / (kF + h * h / (4.0 * kF));
  const double one_minus_cos = sin_a * sin_a / (1.0 + std::sqrt(1.0 - sin_a * sin_a));
  MtfOptions blurred = ideal;
  blurred.pupil.focus_shift = 20.0 * kLambdaMm / one_minus_cos;
  const Mtf m1 = rtt::analysis::mtf(cs, PathId{0}, 0, blurred);
  CHECK(has_code(m1.warnings, "psf.pupil_undersampled"));
  CHECK_FALSE(has_code(m1.warnings, "psf.aliasing"));  // no image lattice
  EncircledEnergyOptions eo;
  eo.pupil = blurred.pupil;
  const std::array<double, 1> r{0.01};
  const EncircledEnergy e1 = rtt::analysis::encircled_energy(cs, PathId{0}, 0, r, eo);
  CHECK(has_code(e1.warnings, "psf.pupil_undersampled"));
  CHECK_FALSE(has_code(e1.warnings, "psf.aliasing"));
}

TEST_CASE("mtf, encircled energy: input errors", "[mtf][ee]") {
  const CompiledSystem cs = paraboloid(10.0);
  for (const double bad : {-1e-3, std::nan(""), std::numeric_limits<double>::infinity()}) {
    const std::array<double, 2> radii{1e-3, bad};
    INFO("radius " << bad);
    CHECK_THROWS_AS(ee_n(cs, 32, radii), std::invalid_argument);
  }
  EncircledEnergyOptions eo;
  eo.wavelengths = std::vector<std::uint16_t>{};
  const std::array<double, 1> r{1e-3};
  CHECK_THROWS_AS(rtt::analysis::encircled_energy(cs, PathId{0}, 0, r, eo), std::invalid_argument);
  eo.wavelengths = std::vector<std::uint16_t>{7};
  CHECK_THROWS_AS(rtt::analysis::encircled_energy(cs, PathId{0}, 0, r, eo), std::invalid_argument);
  MtfOptions mo;
  mo.undersampled_waves = 0.0;
  CHECK_THROWS_AS(rtt::analysis::mtf(cs, PathId{0}, 0, mo), std::invalid_argument);
  // The 2D map is monochromatic only.
  rtt::model::System two = paraboloid_model(10.0);
  two.wavelengths = {{0.5876, 1.0, true}, {0.6563, 1.0, false}};
  MtfOptions map;
  map.outputs.map = true;
  CHECK_THROWS_AS(rtt::analysis::mtf(compiled(two), PathId{0}, 0, map), std::invalid_argument);
  map.wavelengths = std::vector<std::uint16_t>{1};
  CHECK(rtt::analysis::mtf(compiled(two), PathId{0}, 0, map).map.size() == 127U * 127U);
  // An afocal image space (singlet_const with the object in its front focal point) has no MTF.
  rtt::model::System s =
      rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m1/singlet_const.rtt.json");
  const auto fo = rtt::paraxial::first_order(compiled(s), PathId{0}, 0);
  REQUIRE(fo.front_focal_z.has_value());
  s.object.at_infinity = false;
  s.object.distance = rtt::model::Param(-*fo.front_focal_z);
  s.fields.points.resize(1);
  const CompiledSystem afocal = compiled(s);
  CHECK_THROWS_AS(rtt::analysis::mtf(afocal, PathId{0}, 0), rtt::analysis::AnalysisError);
  CHECK_THROWS_AS(rtt::analysis::encircled_energy(afocal, PathId{0}, 0, r),
                  rtt::analysis::AnalysisError);
}
