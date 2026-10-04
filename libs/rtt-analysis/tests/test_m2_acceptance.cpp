// Acceptance protocol of milestone M2 (#34): the reference cases of docs/architecture.md
// ("Referenzfälle" and the M2 row of the roadmap) as tests. Tag [m2].
//
// 1. Spherical aberration of a plano-convex lens: Seidel S_I = fit to real rays at small NA,
//    relative 1e-3.
// 2. Longitudinal colour of a cemented achromat (catalogue glass of the test catalogue) against
//    an independent computation.
// 3. Cooke triplet as golden test: design data from J. Sasian, OPTI 517 lecture L20 "Cooke
//    triplet" (p. 7 design, p. 8 aberration coefficients; docs/quellen.md), reference values
//    from Sasian and from an independent computation with ray-optics (BSD-3).
//
// Every tolerance is justified before the comparison (rule 2), in the comment of its test.

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <numbers>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "rtt/analysis/chromatic.hpp"
#include "rtt/analysis/opd.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/air.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
#include "rtt/paraxial/paraxial.hpp"
#include "rtt/paraxial/seidel.hpp"
#include "rtt/trace/sequential.hpp"
#include "rtt/trace/sources.hpp"

using rtt::compile::compile;
using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::material::MaterialLibrary;
using rtt::math::Vec3;
using rtt::model::Element;
using rtt::model::Param;
using rtt::model::Pose;
using rtt::model::System;
using rtt::trace::RayStatus;

namespace {

constexpr double kDeg = std::numbers::pi / 180.0;

System load(const std::string& relative) {
  return rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/" + relative);
}

Element& element(System& s, std::size_t i) {
  return std::get<Element>(s.root.children[i].value);
}

/// n'_K u'_K of the paraxial marginal ray of seidel() in image space.
double image_nu(const CompiledSystem& cs,
                const rtt::paraxial::Seidel& seidel,
                std::uint16_t wavelength) {
  const auto ray = rtt::paraxial::trace_ray(cs, PathId{0}, wavelength, seidel.marginal.z,
                                            seidel.marginal.y, seidel.marginal.u);
  return ray.back().n * ray.back().u;
}

/// Traces one ray from a given start (global) along path 0; the ray after the last event.
rtt::trace::RayBatch trace_from(const CompiledSystem& cs,
                                const Vec3& pos,
                                const Vec3& dir,
                                std::uint16_t wavelength) {
  rtt::trace::RayBatch rays(1);
  rays.pos_x()[0] = pos.x();
  rays.pos_y()[0] = pos.y();
  rays.pos_z()[0] = pos.z();
  rays.dir_x()[0] = dir.x();
  rays.dir_y()[0] = dir.y();
  rays.dir_z()[0] = dir.z();
  rays.wl()[0] = wavelength;
  [[maybe_unused]] const auto stats = rtt::trace::SequentialTracer().trace(cs, PathId{0}, rays);
  return rays;
}

/// Least squares v = a x^p + b x^q; returns (a, b).
std::pair<double, double> fit_two(const std::vector<double>& x,
                                  const std::vector<double>& v,
                                  int p,
                                  int q) {
  double spp = 0.0, spq = 0.0, sqq = 0.0, spv = 0.0, sqv = 0.0;
  for (std::size_t i = 0; i < x.size(); ++i) {
    const double xp = std::pow(x[i], p);
    const double xq = std::pow(x[i], q);
    spp += xp * xp;
    spq += xp * xq;
    sqq += xq * xq;
    spv += xp * v[i];
    sqv += xq * v[i];
  }
  const double det = spp * sqq - spq * spq;
  return {(spv * sqq - sqv * spq) / det, (sqv * spp - spv * spq) / det};
}

}  // namespace

// ---------------------------------------------------------------------------------------------
// 1. Plano-convex lens
// ---------------------------------------------------------------------------------------------

TEST_CASE("M2: plano-convex lens, Seidel S_I equals a fit to real rays at small NA", "[m2]") {
  // m1/singlet_const (n = 1.5168, f = 100 mm) in vacuum, EPD 2 mm (NA about 0.01), image surface
  // in the paraxial focus, axial field; both orientations of the lens. The stop lies in air in
  // front of the lens and is the entrance pupil, so the normalised pupil coordinate rho of the
  // aimed rays is the paraxial one of seidel().
  //
  // Transverse fit: in the paraxial image plane a wavefront W(rho) gives the transverse
  // aberration eps = (dW/drho) / (n'_K u'_K) (seidel.hpp, checked against real rays in #30).
  // With W = S_I rho^4 / 8 + W060 rho^6 + ... (Sasian, OPTI 517 L4 p. 23) this is
  // eps = S_I rho^3 / (2 n'u') + c5 rho^5 + O(rho^7). Least squares on rho^3 and rho^5 for
  // rho = 0.2 ... 1: S_I,fit = 2 n'u' a.
  // OPD fit: opd_fan() (#29) refers W to a sphere centred on the chief ray in the paraxial focus,
  // so there is no defocus term: W = S_I rho^4 / 8 + W060 rho^6 + ...; least squares on rho^4
  // and rho^6: S_I,fit = 8 a lambda. This also settles the S_I reference value for the OPD that
  // #35 asked for.
  // Tolerance relative 1e-3 (architecture, acceptance table). After the two-term fits the
  // neglected terms are of relative order NA^4 ~ 1e-8.
  const auto check = [](System s) {
    const MaterialLibrary lib;
    s.environment.medium = "VACUUM";
    s.aperture = {rtt::model::SystemApertureType::EntrancePupilDiameter, Param(2.0)};
    const CompiledSystem initial = compile(s, lib);
    const std::uint16_t ref = initial.reference_wavelength();
    element(s, 2).pose =
        Pose::along_z(*rtt::paraxial::first_order(initial, PathId{0}, ref).rear_focal_z);
    const CompiledSystem cs = compile(s, lib);
    const auto seidel = rtt::paraxial::seidel(cs, PathId{0}, ref);
    const double nu = image_nu(cs, seidel, ref);

    std::vector<double> rho;
    std::vector<double> eps;
    for (int i = 1; i <= 5; ++i) {
      rho.push_back(0.2 * i);
      const auto aimed = rtt::trace::aim_ray(cs, PathId{0}, 0, ref, 0.0, rho.back());
      const auto rays = trace_from(cs, aimed.ray.pos, aimed.ray.dir, ref);
      REQUIRE(rays.status()[0] == RayStatus::Alive);
      eps.push_back(rays.pos_y()[0]);
    }
    const double s1_rays = 2.0 * nu * fit_two(rho, eps, 3, 5).first;

    rtt::analysis::OpdOptions options;
    options.fan_points = 21;
    const auto fan = rtt::analysis::opd_fan(cs, PathId{0}, 0, ref, options);
    const double lambda_mm = cs.wavelengths_um()[ref] * 1e-3;
    std::vector<double> r;
    std::vector<double> w;
    for (const auto& p : fan.tangential) {
      if (p.py <= 0.0) continue;
      REQUIRE(p.status == RayStatus::Alive);
      r.push_back(p.py);
      w.push_back(p.w * lambda_mm);
    }
    const double s1_opd = 8.0 * fit_two(r, w, 4, 6).first;

    const double s1 = seidel.sum.s1;
    INFO("S_I " << s1 << ", ray fit " << s1_rays << ", OPD fit " << s1_opd);
    REQUIRE(std::abs(s1_rays - s1) <= 1e-3 * std::abs(s1));
    REQUIRE(std::abs(s1_opd - s1) <= 1e-3 * std::abs(s1));
    return s1;
  };
  // Curved side towards the object (as in the file).
  const double curved_first = check(load("m1/singlet_const.rtt.json"));
  // Flat side towards the object: S1 plane, S2 with R = -51.68 mm.
  System flipped = load("m1/singlet_const.rtt.json");
  Element& lens = element(flipped, 1);
  lens.surfaces[0].shape.base = rtt::model::Plane{};
  lens.surfaces[1].shape.base = rtt::model::Conic{Param(-51.68), Param(0.0)};
  const double flat_first = check(flipped);
  // Both under-corrected (S_I > 0, seidel.hpp); the orientation changes S_I.
  REQUIRE(curved_first > 0.0);
  REQUIRE(flat_first > curved_first);
}

// ---------------------------------------------------------------------------------------------
// 2. Achromat
// ---------------------------------------------------------------------------------------------

namespace {

/// Independent real trace in the meridional plane, written here without rtt-trace, of the
/// cemented doublet in m2/achromat.rtt.json (vertices z = 5, 11, 14 mm; R = 61.5, -35.2,
/// -136.9 mm). Start parallel to the axis at height y0 in the stop plane z = 0; returns the
/// global z where the ray crosses the axis. Intersection of the line with the circle (centre
/// z_v + R) on the cap through the vertex; refraction with the vector form of B. de Greve,
/// "Reflections and Refractions in Ray Tracing" (2006), Eqs. (22), (23), (28) (docs/quellen.md):
/// t = mu d + (mu cos_i - cos_t) N with mu = n1/n2, cos_i = -d.N, N oriented into the incident
/// medium, sin_t^2 = mu^2 (1 - cos_i^2).
double achromat_axis_crossing(double y0, const std::array<double, 4>& n) {
  const std::array<double, 3> zv{5.0, 11.0, 14.0};
  const std::array<double, 3> radius{61.5, -35.2, -136.9};
  double py = y0;
  double pz = 0.0;
  double dy = 0.0;
  double dz = 1.0;
  for (std::size_t i = 0; i < 3; ++i) {
    const double zc = zv[i] + radius[i];
    const double oz = pz - zc;
    const double b = py * dy + oz * dz;
    const double c = py * py + oz * oz - radius[i] * radius[i];
    const double root = std::sqrt(b * b - c);
    const double t1 = -b - root;
    const double t2 = -b + root;
    // The cap through the vertex: the crossing closer to the vertex plane.
    const double t = std::abs(pz + t1 * dz - zv[i]) < std::abs(pz + t2 * dz - zv[i]) ? t1 : t2;
    py += t * dy;
    pz += t * dz;
    double ny = py / radius[i];
    double nz = (pz - zc) / radius[i];
    if (ny * dy + nz * dz > 0.0) {
      ny = -ny;
      nz = -nz;
    }
    const double mu = n[i] / n[i + 1];
    const double cos_i = -(ny * dy + nz * dz);
    const double sin_t2 = mu * mu * (1.0 - cos_i * cos_i);
    const double cos_t = std::sqrt(1.0 - sin_t2);
    const double ty = mu * dy + (mu * cos_i - cos_t) * ny;
    const double tz = mu * dz + (mu * cos_i - cos_t) * nz;
    dy = ty;
    dz = tz;
  }
  return pz - py * dz / dy;
}

double medium_index(const CompiledSystem& cs, const std::string& reference, std::uint16_t wl) {
  for (const auto& m : cs.media()) {
    if (m.reference == reference) return m.index[wl].real();
  }
  FAIL("medium " << reference << " not found");
  return 0.0;
}

void add_schott(MaterialLibrary& lib) {
  lib.add_catalog(std::string(RTT_CATALOG_DIR) + "/schott.agf");
}

}  // namespace

TEST_CASE("M2: achromat, real longitudinal colour against an independent meridional trace",
          "[m2]") {
  // m2/achromat.rtt.json (N-BK7/F2 from the test catalogue). The zone ray of
  // longitudinal_colour() is the axial ray at normalised pupil (0, zone); the stop lies in air in
  // front of the lens and is the entrance pupil (radius 10 mm), so the ray runs parallel to the
  // axis at height 10 zone. The paraxial part was checked against an independent y-nu trace in
  // #31; here the real foci are checked against achromat_axis_crossing() for F, d, C.
  // Tolerance 1e-9 mm, a priori: real aiming meets the stop to < kAimTolerance = 1e-9 mm, which
  // moves the crossing by |dz/dy| * 1e-9 mm with |dz/dy| ~ 2 LSA / y ~ 0.02 here; rounding in
  // three refractions and the final division by the exit slope (~0.1) stays below 1e-12 mm.
  MaterialLibrary lib;
  add_schott(lib);
  const CompiledSystem cs = compile(load("m2/achromat.rtt.json"), lib);
  for (const double zone : {1.0, 0.5}) {
    rtt::analysis::ChromaticOptions options;
    options.zone = zone;
    const auto lc = rtt::analysis::longitudinal_colour(cs, PathId{0}, options);
    std::array<double, 3> crossing{};
    for (std::uint16_t wl = 0; wl < 3; ++wl) {
      const double n0 = cs.media()[cs.environment_medium()].index[wl].real();
      crossing[wl] = achromat_axis_crossing(10.0 * zone, {n0, medium_index(cs, "SCHOTT:N-BK7", wl),
                                                          medium_index(cs, "SCHOTT:F2", wl), n0});
      INFO("zone " << zone << ", wavelength " << wl << ": " << lc.foci[wl].real_z << " / "
                   << crossing[wl]);
      REQUIRE(std::abs(lc.foci[wl].real_z - crossing[wl]) <= 1e-9);
    }
    REQUIRE(std::abs(lc.real - (crossing[lc.pair.first] - crossing[lc.pair.second])) <= 1e-9);
  }
}

TEST_CASE("M2: achromat, Seidel C_L against the paraxial colour of a close wavelength pair",
          "[m2]") {
  // W contains 1/2 C_L rho^2 for the pair (first, second), dn = n(first) - n(second)
  // (seidel.hpp; Sasian, OPTI 517 L4 p. 23: d_lambda W020 = C_L / 2; C_L per surface: OPTI 518
  // L6). A defocus term W020 moves the focus by eps_z = -2 W020 / (n'u'^2) along the light
  // (Wyant & Creath, Eq. (18), as in the field-curvature test of #31), so to first order in dn
  // the paraxial focus of `first` lies dz = -C_L / (n'_K u'_K^2) from that of `second`.
  // F and C are unsuitable for this first-order relation in an achromat: the contributions of
  // the two glasses (each about the colour of a singlet, ~1.6 mm here) cancel, and the
  // second-order term per glass, about 1/2 (B''/B') (dF + dC) dn with B the focus as a function
  // of the index and dF, dC the index offsets of F and C from d, is comparable to the remaining
  // first-order sum. With the symmetric pair lambda_d -+ 1e-4 um, dF + dC = n'' dlambda^2 and the
  // second-order remainder is <~ 1e-6 relative. Tolerance relative 1e-4, a priori.
  MaterialLibrary lib;
  add_schott(lib);
  System s = load("m2/achromat.rtt.json");
  s.wavelengths = {{0.5875, 1.0, false}, {0.5876, 1.0, true}, {0.5877, 1.0, false}};
  const CompiledSystem cs = compile(s, lib);
  const rtt::paraxial::ChromaticPair pair{0, 2};
  const auto seidel = rtt::paraxial::seidel(cs, PathId{0}, 1, pair);
  const auto marginal = rtt::paraxial::trace_ray(cs, PathId{0}, 1, seidel.marginal.z,
                                                 seidel.marginal.y, seidel.marginal.u);
  const double nu2 = marginal.back().n * marginal.back().u * marginal.back().u;
  const double expected = -seidel.sum.c_l / nu2;
  rtt::analysis::ChromaticOptions options;
  options.pair = pair;
  const auto lc = rtt::analysis::longitudinal_colour(cs, PathId{0}, options);
  INFO("paraxial " << lc.paraxial << ", from C_L " << expected);
  REQUIRE(expected != 0.0);
  REQUIRE(std::abs(lc.paraxial - expected) <= 1e-4 * std::abs(expected));
}

// ---------------------------------------------------------------------------------------------
// 3. Cooke triplet (golden test)
// ---------------------------------------------------------------------------------------------
//
// tests/reference/m2/cooke_triplet.rtt.json: Sasian L20 p. 7 (design attributed to Geiser):
// f = 50 mm, f/4 (EPD 12.5 mm), field 20 deg, R/d = 23.713/4.831 (LAK9), 7331.288/5.86,
// -24.456/0.975 (SF5, stop), 21.896/4.822, 86.759/3.127 (LAK9), -20.4942/41.10346.
// Glasses N-LAK9 and N-SF5 from the free SCHOTT AGF (tests/catalogs/m2/schott.agf, coefficients
// identical to opticalglass 2.0.2 used by ray-optics). LAK9 of Sasian has the same coefficients
// as N-LAK9 in that AGF; SF5 differs from N-SF5 (nd 1.67270 / 1.67271).
//
// Air (decision for #34, variant A): the sources compute relative to air (air = 1) with
// wavelengths in air. The file uses AIR (Ciddor) at 20 degC / 1 atm, the reference temperature
// of the catalogue, and the vacuum wavelengths lambda_vac = lambda_air n_air(lambda_vac) of the
// F, d, C air wavelengths 486.1327, 587.5618, 656.2725 nm. A catalogue glass is evaluated at
// lambda_vac / n_air = lambda_air and delivered as n_rel n_air (#25), so every index ratio equals
// that of the sources: geometric quantities (BFL, ray intercepts) agree directly, the rear focal
// length f'_R = n'/Phi equals the relative focal length (EFL_abs = 1/Phi = f_rel / n_air), the
// Seidel sums are S_abs = n_air S_rel, and in waves W = S_abs / (k lambda_vac) = S_rel /
// (k lambda_air).
//
// Stop: Sasian's stop is the surface L2.S1. A stop is an element of its own here, and a stop
// surface coincident with L2.S1 would be hit at t = 0 < t_min = 1e-9 mm. The plane STO
// therefore lies delta = 1e-8 mm behind the vertex of L2.S1, inside the glass, on an explicit
// path (Transmit keeps the medium, #5). This shifts the paraxial stop; with the stop shifting
// parameter S = (y_bar,new - y_bar,old) / y (Sasian, OPTI 518 L14 p. 22) the sums change by
// dS_II = S S_I, dS_III = 2 S S_II + S^2 S_I, dS_V = S (S_IV + 3 S_III) + ... (L14 p. 23).
// At the stop, y_bar,new - y_bar,old = u_bar delta with the chief slope u_bar in the glass, so
// |S| = |u_bar| delta / y ~ 0.3 * 1e-8 / 5 < 1e-9 (checked below from the paraxial data). The
// largest change of a single value is then |dS_V| <= |S| (|S_IV| + 3 |S_III|) + O(S^2) <
// 0.65e-9 mm (|S_IV| + 3 |S_III| < 0.65 mm for every surface of this design); dS_II <= |S| |S_I|
// and dS_III <= 2 |S| |S_II| + S^2 |S_I| are smaller, S_I and S_IV do not change.

namespace {

constexpr std::array<double, 3> kAirWavelengthsUm{0.4861327, 0.5875618, 0.6562725};
constexpr double kStopOffset = 1e-8;  // mm, STO behind the vertex of L2.S1

/// Vacuum wavelength whose wavelength in air (Ciddor, 20 degC, 1 atm) is `air_um`.
double vacuum_wavelength(double air_um) {
  double v = air_um;
  for (int i = 0; i < 10; ++i) v = air_um * rtt::material::ciddor_air_index(v, 20.0, 1.0);
  return v;
}

void add_triplet_glasses(MaterialLibrary& lib) {
  lib.add_catalog(std::string(RTT_CATALOG_DIR) + "/m2/schott.agf");
}

/// Path events of the six refracting surfaces (event 3 is STO, event 7 the image).
constexpr std::array<std::size_t, 6> kSurfaceEvents{0, 1, 2, 4, 5, 6};

}  // namespace

TEST_CASE("M2: Cooke triplet, file wavelengths are the vacuum wavelengths of F, d, C in air",
          "[m2]") {
  // Rounding the file values by 1e-12 um moves lambda_air by the same amount and the glass
  // indices by < 1e-12 * |dn/dlambda| < 1e-12; negligible against every tolerance below.
  const System s = load("m2/cooke_triplet.rtt.json");
  REQUIRE(s.wavelengths.size() == 3);
  REQUIRE(s.environment.medium == "AIR");
  REQUIRE(s.environment.temperature_c == 20.0);
  REQUIRE(s.environment.pressure_atm == 1.0);
  for (std::size_t i = 0; i < 3; ++i) {
    const double v = vacuum_wavelength(kAirWavelengthsUm[i]);
    INFO(std::setprecision(17) << "wavelength " << i << ": file " << s.wavelengths[i].um
                               << ", vacuum " << v);
    REQUIRE(std::abs(s.wavelengths[i].um - v) <= 1e-12);
  }
}

TEST_CASE("M2: Cooke triplet, first-order data", "[m2]") {
  // References:
  // - ray-optics 0.9.10 run of the same prescription (script in the PR of #34): efl, bfl, ffl
  //   and enp_dist at full double precision. Tolerance 1e-9 mm for f, BFL and FFL, a priori:
  //   both are double y-nu traces with the same index ratios, rounding ~1e-13 mm, and they do
  //   not depend on the stop. The entrance pupil depends on the paraxial chief ray, which
  //   ray-optics computes from an object at 1e10 mm (its model of an object at infinity): the
  //   chief ray starts at |y_bar| = 1e10 tan(20 deg) = 3.64e9 mm and loses digits up to the first
  //   surface, about 4 ulp * 3.64e9 mm = 3.2e-6 mm in y_bar, i.e. 3.2e-6 / tan(20 deg) < 1e-5 mm
  //   in z_EP (found in #34, see the Seidel test); tolerance 1e-5 mm against ray-optics.
  // - Independent y-nu trace of the stop centre through L2.S1 and L1 (Greivenkamp, OPTI-201/202
  //   Sec. 9, p. 9-2; Python script in the PR of #34, without rtt): z_EP = 11.679209268268448 mm
  //   with the stop 1e-8 mm behind the vertex of L2.S1 (11.679209257540240 mm without the
  //   offset). Tolerance 1e-9 mm (double rounding ~1e-14 mm).
  // - ray-optics documentation, "SasianTriplet" example (printed with 4 significant digits):
  //   efl 50, bfl 41.24, ffl -37.1, enp_dist 11.68; tolerance half the last digit.
  // - Sasian L20 p. 7: f = 50 mm; tolerance 0.5 mm (integer).
  // ray-optics measures ffl from the first vertex along +z (negative = F in front),
  // rtt::paraxial positive if F lies upstream.
  MaterialLibrary lib;
  add_triplet_glasses(lib);
  const CompiledSystem cs = compile(load("m2/cooke_triplet.rtt.json"), lib);
  const std::uint16_t ref = cs.reference_wavelength();
  const double n_air = cs.media()[cs.environment_medium()].index[ref].real();
  const auto fo = rtt::paraxial::first_order(cs, PathId{0}, ref);
  REQUIRE(fo.rear_focal_length);
  REQUIRE(fo.efl);
  REQUIRE(fo.bfl);
  REQUIRE(fo.ffl);
  REQUIRE(fo.entrance_pupil);
  const double f_rel = *fo.rear_focal_length;
  const double z_ep = *fo.entrance_pupil->z;
  INFO(std::setprecision(17) << "f'_R " << f_rel << ", EFL " << *fo.efl << ", BFL " << *fo.bfl
                             << ", FFL " << *fo.ffl << ", z_EP " << z_ep);
  REQUIRE(std::abs(f_rel - 50.00098928993641) <= 1e-9);
  REQUIRE(std::abs(*fo.efl * n_air - 50.00098928993641) <= 1e-9);
  REQUIRE(std::abs(*fo.bfl - 41.23764895362289) <= 1e-9);
  REQUIRE(std::abs(*fo.ffl - 37.102124967857094) <= 1e-9);
  REQUIRE(std::abs(z_ep - 11.679207956830547) <= 1e-5);
  REQUIRE(std::abs(z_ep - 11.679209268268448) <= 1e-9);
  REQUIRE(std::abs(*fo.entrance_pupil->diameter - 12.5) <= 1e-12);
  REQUIRE(std::abs(f_rel - 50.0) <= 0.005);
  REQUIRE(std::abs(*fo.bfl - 41.24) <= 0.005);
  REQUIRE(std::abs(*fo.ffl - 37.1) <= 0.05);
  REQUIRE(std::abs(z_ep - 11.68) <= 0.005);
  REQUIRE(std::abs(f_rel - 50.0) <= 0.5);
}

TEST_CASE("M2: Cooke triplet, Seidel sums against ray-optics and Sasian L20", "[m2]") {
  MaterialLibrary lib;
  add_triplet_glasses(lib);
  const CompiledSystem cs = compile(load("m2/cooke_triplet.rtt.json"), lib);
  const std::uint16_t ref = cs.reference_wavelength();
  const double n_air = cs.media()[cs.environment_medium()].index[ref].real();
  const auto seidel = rtt::paraxial::seidel(cs, PathId{0}, ref);
  REQUIRE(seidel.surfaces.size() == 8);

  // Stop shifting parameter of the offset STO (see the comment above the triplet tests).
  const auto chief =
      rtt::paraxial::trace_ray(cs, PathId{0}, ref, seidel.chief.z, seidel.chief.y, seidel.chief.u);
  const auto marginal = rtt::paraxial::trace_ray(cs, PathId{0}, ref, seidel.marginal.z,
                                                 seidel.marginal.y, seidel.marginal.u);
  const double stop_shift = chief[3].u * kStopOffset / marginal[3].y;
  INFO("stop shifting parameter " << stop_shift);
  REQUIRE(std::abs(stop_shift) < 1e-9);

  // S_I ... S_V relative to air, mm, per refracting surface and the sum.
  const auto terms = [&](std::size_t k) {
    const auto& t = k < 6 ? seidel.surfaces[kSurfaceEvents[k]].terms : seidel.sum;
    return std::array<double, 5>{t.s1 / n_air, t.s2 / n_air, t.s3 / n_air, t.s4 / n_air,
                                 t.s5 / n_air};
  };

  // (a) ray-optics 0.9.10 run (script in the PR of #34), printed to 16 digits.
  // S_I and S_IV (marginal ray and Lagrange invariant only): tolerance 1e-9 mm, a priori: the
  // stop offset does not change them (L14 p. 23), the two double implementations differ by
  // rounding (~1e-13 mm).
  // S_II, S_III and S_V depend on the paraxial chief ray. ray-optics models the object at
  // infinity at 1e10 mm, so its chief ray starts at |y_bar| = 1e10 tan(20 deg) = 3.64e9 mm and
  // loses digits up to the first surface: the error in y_bar is about 4 ulp * 3.64e9 mm =
  // 4 * 2.2e-16 * 3.64e9 mm = 3.2e-6 mm. Such an error acts like a stop shift with
  // S = dy_bar / y <= 3.2e-6 / 6.25 = 5e-7 (y = r_EP; Sasian, OPTI 518 L14 p. 22), which changes
  // a value by at most |dS| <= S (|S_IV| + 3 |S_III|) + O(S^2) <= 5e-7 * 0.65 mm = 3.3e-7 mm
  // (L14 p. 23, as above). Tolerance 5e-7 mm. Found in #34: before this bound the comparison at
  // 1e-9 mm failed by up to 1e-8 mm; an independent y-nu trace confirmed rtt (entrance pupil, see
  // the first-order test), and the stop-shift formulas with the measured pupil offset of
  // ray-optics (1.31e-6 mm) predict the differences of surface 1 (2.1e-9, 3.0e-9, 9.9e-9 mm).
  // (b) ray-optics documentation, same table printed with 6 decimals: tolerance 5e-7 mm (half
  // the last digit) + the bound of (a) (5e-7 mm for S_II, S_III, S_V; 1e-9 mm otherwise).
  constexpr std::array<std::array<double, 5>, 7> kRayOptics{{
      {2.765359752522226e-02, 1.937938497952812e-02, 1.358089347479001e-02, 8.917445783138278e-02,
       7.200999833214197e-02},
      {2.208227484777396e-02, -5.950122499892305e-02, 1.603274934660710e-01, -2.884341630768809e-04,
       -4.312291257059224e-01},
      {-1.051559239169821e-01, 1.376919806815802e-01, -1.802949452375537e-01,
       -8.509664370932608e-02, 3.475058006924731e-01},
      {-4.535808963116495e-02, -7.679623722947593e-02, -1.300244807610623e-01,
       -9.504583113606498e-02, -3.810688061670026e-01},
      {7.941672506770307e-03, 2.838184712564030e-02, 1.014306804487971e-01, 2.437319377304464e-02,
       4.495962686617573e-01},
      {1.038095714711020e-01, -5.006777394154603e-02, 2.414788879231219e-02, 1.031801152792292e-01,
       -6.141080859828842e-02},
      {1.097310280272150e-02, -9.120233831964100e-04, -1.083246981664562e-02, 3.629685787518870e-02,
       -4.596672784840904e-03},
  }};
  constexpr std::array<std::array<double, 5>, 7> kRayOpticsDoc{{
      {0.027654, 0.019379, 0.013581, 0.089174, 0.072010},
      {0.022082, -0.059501, 0.160327, -0.000288, -0.431229},
      {-0.105156, 0.137692, -0.180295, -0.085097, 0.347506},
      {-0.045358, -0.076796, -0.130024, -0.095046, -0.381069},
      {0.007942, 0.028382, 0.101431, 0.024373, 0.449596},
      {0.103810, -0.050068, 0.024148, 0.103180, -0.061411},
      {0.010973, -0.000912, -0.010832, 0.036297, -0.004597},
  }};
  for (std::size_t k = 0; k < 7; ++k) {
    const auto s = terms(k);
    for (std::size_t j = 0; j < 5; ++j) {
      INFO(std::setprecision(17) << "row " << k << " (6 = sum), S_" << j + 1 << ": " << s[j]);
      const bool chief_dependent = j == 1 || j == 2 || j == 4;  // S_II, S_III, S_V
      REQUIRE(std::abs(s[j] - kRayOptics[k][j]) <= (chief_dependent ? 5e-7 : 1e-9));
      REQUIRE(std::abs(s[j] - kRayOpticsDoc[k][j]) <= 5e-7 + (chief_dependent ? 5e-7 : 1e-9));
    }
  }

  // Wavefront coefficients in waves at 587.5618 nm (in air): W040 = S_I/8, W131 = S_II/2,
  // W222 = S_III/2, W220 = (S_IV + S_III)/4, W311 = S_V/2 (Sasian, OPTI 517 L4 p. 23).
  constexpr double kLambdaAirMm = 0.5875618e-3;
  const auto sum = terms(6);
  const std::array<double, 5> w_sum{
      sum[0] / 8.0 / kLambdaAirMm, sum[1] / 2.0 / kLambdaAirMm, sum[2] / 2.0 / kLambdaAirMm,
      (sum[3] + sum[2]) / 4.0 / kLambdaAirMm, sum[4] / 2.0 / kLambdaAirMm};
  // (c) ray-optics documentation, wavefront coefficients printed with 6 decimals: tolerance
  // 5e-7 waves (half the last digit) + the chief-ray error of ray-optics in the sums, with
  // S <= 5e-7 from (a): dS_II = S S_I <= 5.5e-9 mm, dS_III <= 2 S |S_II| < 1e-9 mm,
  // dS_V <= S |S_IV + 3 S_III| < 2e-9 mm (sums), i.e. <= 5.5e-9 mm / (2 lambda) < 5e-6 waves;
  // the offset STO adds < 1e-9 mm. Tolerance 5e-7 + 5e-6 waves.
  constexpr std::array<double, 5> kWaveDoc{2.334457, -0.776108, -9.218154, 10.834770, -3.911650};
  for (std::size_t j = 0; j < 5; ++j) {
    INFO("W sum " << j << ": " << w_sum[j]);
    REQUIRE(std::abs(w_sum[j] - kWaveDoc[j]) <= 5e-7 + 5e-6);
  }

  // (d) Sasian L20 p. 8: W040, W131, W222, W220 and W311 per surface and total, in waves. The
  // column W220 there is the Petzval part S_IV / 4 (surface 1: 0.0891745 / 4 / lambda =
  // 37.9424); the W220 with S_III is the 10.8356 listed above the table. Sasian uses LAK9 and
  // SF5: LAK9 has the coefficients of N-LAK9, SF5 differs from N-SF5. Tolerance relative 1e-3
  // for this glass change (decision for #34; ray-optics with N-SF5 and Sasian differ by at most
  // 3.2e-4 in the printed totals).
  constexpr std::array<std::array<double, 5>, 7> kSasian{{
      {5.883061, 16.491222, 11.556926, 37.942379, 61.278483},
      {4.697811, -50.633600, 136.433840, -0.122724, -366.963936},
      {-22.370883, 117.170758, -153.424726, -36.207024, 295.715864},
      {-9.649013, -65.348394, -110.643768, -40.440216, -324.276663},
      {1.689360, 24.150389, 86.310980, 10.370424, 382.592103},
      {22.084875, -42.606408, 20.549199, 43.901573, -52.258664},
      {2.335211, -0.776033, -9.217549, 15.444412, -3.912814},
  }};
  for (std::size_t k = 0; k < 7; ++k) {
    const auto s = terms(k);
    const std::array<double, 5> w{s[0] / 8.0 / kLambdaAirMm, s[1] / 2.0 / kLambdaAirMm,
                                  s[2] / 2.0 / kLambdaAirMm, s[3] / 4.0 / kLambdaAirMm,
                                  s[4] / 2.0 / kLambdaAirMm};
    for (std::size_t j = 0; j < 5; ++j) {
      INFO("Sasian row " << k << " (6 = total), column " << j << ": " << w[j] << " / "
                         << kSasian[k][j]);
      REQUIRE(std::abs(w[j] - kSasian[k][j]) <= 1e-3 * std::abs(kSasian[k][j]));
    }
  }
}

TEST_CASE("M2: Cooke triplet, real rays against ray-optics", "[m2]") {
  // Fixed start conditions, no aiming (decision for #34): start plane z = -10 mm (10 mm before
  // the first vertex), direction (0, sin t, cos t) for the field angle t, start point
  // (px * 6.25, -(10 + 11.68) tan t + py * 6.25) for chief (0, 0), upper (0, 1), lower (0, -1)
  // and sagittal (1, 0) rays; fields 0, 14, 20 deg; F, d, C. Reference: intercept (x, y) on the
  // image plane and direction cosines in image space from ray-optics 0.9.10 (script in the PR of
  // #34, trace.trace with the same start), printed to 17 digits.
  // Tolerance 1e-8 mm and 1e-8 for the direction cosines (decision for #34), a priori: same
  // prescription and start, index ratios equal up to rounding (lambda_rel = lambda_air to
  // 1e-12 um); two correct double ray tracers agree to about 1e-9 mm after six refractions and
  // ~60 mm of path. STO inside L2 only transmits.
  // The RMS spot is not compared: no ray-optics pupil sampling is identical to the real aiming
  // on the stop used by rtt-analysis.
  struct Reference {
    double field_deg;
    std::uint16_t wavelength;
    int px, py;
    double x, y, l, m, n;
  };
  // clang-format off
  constexpr std::array<Reference, 36> kRays{{
    {0.0, 0, 0, 0, 0, 0, 0, 0, 1},
    {0.0, 0, 0, 1, 0, 0.0045451578786283205, 0, -0.12515525690852666, 0.9921371687766567},
    {0.0, 0, 0, -1, 0, -0.0045451578786283205, 0, 0.12515525690852666, 0.9921371687766567},
    {0.0, 0, 1, 0, 0.0045451578786283205, 0, -0.12515525690852666, 0, 0.9921371687766567},
    {0.0, 1, 0, 0, 0, 0, 0, 0, 1},
    {0.0, 1, 0, 1, 0, -0.0011614290405511025, 0, -0.12527951809357171, 0.99212148567916969},
    {0.0, 1, 0, -1, 0, 0.0011614290405511025, 0, 0.12527951809357171, 0.99212148567916969},
    {0.0, 1, 1, 0, -0.0011614290405511025, 0, -0.12527951809357171, 0, 0.99212148567916969},
    {0.0, 2, 0, 0, 0, 0, 0, 0, 1},
    {0.0, 2, 0, 1, 0, 0.0025692994720563837, 0, -0.12521173458297458, 0.99213004264699245},
    {0.0, 2, 0, -1, 0, -0.0025692994720563837, 0, 0.12521173458297458, 0.99213004264699245},
    {0.0, 2, 1, 0, 0.0025692994720563837, 0, -0.12521173458297458, 0, 0.99213004264699245},
    {14.0, 0, 0, 0, 0, 12.439531353421028, 0, 0.23271929326724505, 0.9725439478713509},
    {14.0, 0, 0, 1, 0, 12.500929045740621, 0, 0.11937765295025005, 0.99284891900837047},
    {14.0, 0, 0, -1, 0, 12.408013239430229, 0, 0.34662876074334936, 0.93800239990392786},
    {14.0, 0, 1, 0, -0.0098206300195425311, 12.441178247893465, -0.12171647685041223, 0.2326262876320212, 0.96491974255150792},
    {14.0, 1, 0, 0, 0, 12.438188066064237, 0, 0.23298384282732224, 0.97248060596672747},
    {14.0, 1, 0, 1, 0, 12.499956888873884, 0, 0.11971349348090388, 0.99280848076484451},
    {14.0, 1, 0, -1, 0, 12.405670949683309, 0, 0.34696640632310594, 0.93787755751122992},
    {14.0, 1, 1, 0, -0.015277218115981095, 12.439837956916877, -0.12183195460340025, 0.23289094361405044, 0.96484132541059431},
    {14.0, 2, 0, 0, 0, 12.437786953229574, 0, 0.2331026873195419, 0.97245212589844621},
    {14.0, 2, 0, 1, 0, 12.505718163181095, 0, 0.11997273590666047, 0.99277718680430527},
    {14.0, 2, 0, -1, 0, 12.398506726964508, 0, 0.34700988951549594, 0.93786146982293905},
    {14.0, 2, 1, 0, -0.011356057636187467, 12.439390433360883, -0.12176303562006803, 0.23301157990775007, 0.96482089881256217},
    {20.0, 0, 0, 0, 0, 18.151529745461772, 0, 0.32299012339169347, 0.94640233526308415},
    {20.0, 0, 0, 1, 0, 18.187378972063488, 0, 0.21781492284778475, 0.97599009184766494},
    {20.0, 0, 0, -1, 0, 18.143562887863414, 0, 0.42604354598217181, 0.9047026566374925},
    {20.0, 0, 1, 0, 0.0062686607782368453, 18.147329724740402, -0.11768122526381655, 0.32317883885714083, 0.93899231484360768},
    {20.0, 1, 0, 0, 0, 18.150978511797796, 0, 0.32355514267829499, 0.94620931597951907},
    {20.0, 1, 0, 1, 0, 18.204100634369489, 0, 0.21895196032191744, 0.9757356399513093},
    {20.0, 1, 0, -1, 0, 18.122113982370166, 0, 0.42651900093676282, 0.90447860220123832},
    {20.0, 1, 1, 0, 0.0015439528904211963, 18.146560398014994, -0.11777031017379426, 0.32374635703758164, 0.93878562534077137},
    {20.0, 2, 0, 0, 0, 18.151132717237097, 0, 0.32380660879056189, 0.9461232901179194},
    {20.0, 2, 0, 1, 0, 18.217221704358629, 0, 0.21953207832204946, 0.97560528216466713},
    {20.0, 2, 0, -1, 0, 18.107338508203689, 0, 0.42664197794551062, 0.90442060052540996},
    {20.0, 2, 1, 0, 0.0058177118405455808, 18.146547957231785, -0.11769406277596504, 0.32400046471775468, 0.93870751911869033}
  }};
  // clang-format on
  MaterialLibrary lib;
  add_triplet_glasses(lib);
  const CompiledSystem cs = compile(load("m2/cooke_triplet.rtt.json"), lib);
  for (const Reference& r : kRays) {
    const double t = r.field_deg * kDeg;
    const Vec3 start(r.px * 6.25, -(10.0 + 11.68) * std::tan(t) + r.py * 6.25, -10.0);
    const Vec3 dir(0.0, std::sin(t), std::cos(t));
    const auto rays = trace_from(cs, start, dir, r.wavelength);
    INFO(std::setprecision(17) << "field " << r.field_deg << ", wavelength " << r.wavelength
                               << ", pupil (" << r.px << ", " << r.py << "): x " << rays.pos_x()[0]
                               << ", y " << rays.pos_y()[0] << ", z " << rays.pos_z()[0]);
    REQUIRE(rays.status()[0] == RayStatus::Alive);
    REQUIRE(std::abs(rays.pos_z()[0] - 60.71846) <= 1e-12);
    REQUIRE(std::abs(rays.pos_x()[0] - r.x) <= 1e-8);
    REQUIRE(std::abs(rays.pos_y()[0] - r.y) <= 1e-8);
    REQUIRE(std::abs(rays.dir_x()[0] - r.l) <= 1e-8);
    REQUIRE(std::abs(rays.dir_y()[0] - r.m) <= 1e-8);
    REQUIRE(std::abs(rays.dir_z()[0] - r.n) <= 1e-8);
  }
}
