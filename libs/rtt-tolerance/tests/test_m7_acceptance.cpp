// M7 acceptance (#262): tolerancing end to end on the reference systems tests/reference/m7, read
// through rtt-io. Every perturbed trial runs through prepare/evaluate with the nominal reference
// (ADR 0036); no case computes the perturbed system by hand past the tolerance layer. The cases
// carry the numbers of the plan in #262 (plan v3 of #252, sections 4 and 5a). The target values
// and bands were derived before the code and agreed with the coordinator; a red case is a
// finding, not a reason to change the test (AGENTS.md, rule 2).
//
// Bounds:
// - traced coordinates: kEps = 1e-12 mm, as in test_sensitivity.cpp and test_monte_carlo.cpp;
// - paraxial values (efl): in a trial, also at the nominal values, they come from parabasal rays
//   (ADR 0030, addendum #256), so the target is the real parabasal ray and the bound the
//   rounding term of the Richardson bound, as in rtt-optim/tests/test_parabasal.cpp;
// - statistics: the bands of test_monte_carlo.cpp (mean 4 sigma / sqrt(N), Var(s^2) =
//   (mu_4 - (N - 3) / (N - 1) sigma^4) / N with MathWorld "Sample Variance Distribution", Wilson
//   interval with z = 3.8906), each derived in its case.

#include <oneapi/tbb/task_arena.h>

#include <algorithm>
#include <bit>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <initializer_list>
#include <limits>
#include <memory>
#include <numbers>
#include <optional>
#include <string>
#include <variant>
#include <vector>

#include "rtt/coating/catalog.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/compile/perturbation.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/math/units.hpp"
#include "rtt/model/model.hpp"
#include "rtt/model/tolerance.hpp"
#include "rtt/model/tolerance_classes.hpp"
#include "rtt/optim/levenberg_marquardt.hpp"
#include "rtt/optim/optimize.hpp"
#include "rtt/tolerance/classes.hpp"
#include "rtt/tolerance/expand.hpp"
#include "rtt/tolerance/monte_carlo.hpp"
#include "rtt/tolerance/sensitivity.hpp"
#include "rtt/tolerance/trial.hpp"

namespace fs = std::filesystem;

using rtt::model::System;
using rtt::model::ToleranceDistribution;
using rtt::model::ToleranceItem;
using rtt::tolerance::MonteCarloOptions;
using rtt::tolerance::MonteCarloResult;
using rtt::tolerance::TrialContext;
using rtt::tolerance::TrialResult;
using rtt::tolerance::TrialStatus;

namespace {

constexpr double kEps = 1e-12;  // mm, see the head of the file

bool same(double a, double b) {
  return std::bit_cast<std::uint64_t>(a) == std::bit_cast<std::uint64_t>(b);
}

System load(const std::string& name) {
  return rtt::io::load_system(fs::path(RTT_REFERENCE_DIR) / "m7" / name);
}

/// Glass-free systems (CONST and ideal elements).
const rtt::material::MaterialLibrary& no_glass() {
  static const rtt::material::MaterialLibrary lib;
  return lib;
}

const rtt::material::MaterialLibrary& schott() {
  struct Schott {  // MaterialLibrary is not movable
    rtt::material::MaterialLibrary lib;
    Schott() { lib.add_catalog(std::string(RTT_CATALOG_DIR) + "/schott.agf"); }
  };
  static const Schott s;
  return s.lib;
}

/// The items of the file in the order of expand(), restricted to `ids` in that order.
std::vector<ToleranceItem> items(const System& s,
                                 std::initializer_list<const char*> ids,
                                 const rtt::material::MaterialLibrary& lib = no_glass()) {
  const rtt::tolerance::ToleranceItems all = rtt::tolerance::expand(s, lib);
  std::vector<ToleranceItem> out;
  for (const char* id : ids) {
    const auto it = std::find_if(all.items.begin(), all.items.end(),
                                 [&](const ToleranceItem& i) { return i.id == id; });
    REQUIRE(it != all.items.end());
    out.push_back(*it);
  }
  return out;
}

TrialContext context(const System& s,
                     const std::vector<ToleranceItem>& subset,
                     const rtt::material::MaterialLibrary& lib = no_glass()) {
  return rtt::tolerance::prepare(s, lib, nullptr, subset);
}

/// Index of the criterion `id` in the order of the context.
std::size_t criterion(const TrialContext& ctx, const std::string& id) {
  const auto c = ctx.criteria();
  const auto it = std::find_if(c.begin(), c.end(), [&](const auto& x) { return x.id == id; });
  REQUIRE(it != c.end());
  return static_cast<std::size_t>(it - c.begin());
}

/// The value of criterion `id` in a valid trial with `values`.
double value(const TrialContext& ctx, const std::vector<double>& values, const std::string& id) {
  const TrialResult r = rtt::tolerance::evaluate(ctx, values);
  INFO("reason: " << r.reason);
  REQUIRE(r.status == TrialStatus::Valid);
  const std::optional<double> v = r.criteria[criterion(ctx, id)].value;
  REQUIRE(v.has_value());
  return *v;
}

double nominal_value(const TrialContext& ctx, const std::string& id) {
  const std::optional<double> v = ctx.nominal().criteria[criterion(ctx, id)].value;
  REQUIRE(v.has_value());
  return *v;
}

/// The global z of the origin of surface `id` in `system` with `perturbation`.
double global_z(const System& system,
                const rtt::compile::Perturbation& perturbation,
                const std::string& id) {
  const rtt::coating::CoatingLibrary coatings;
  const rtt::compile::CompiledSystem cs =
      rtt::compile::compile(system, no_glass(), coatings, 0, perturbation);
  const std::optional<std::uint32_t> k = cs.find_surface(rtt::model::SurfaceId(id));
  REQUIRE(k.has_value());
  return cs.surfaces()[*k].to_global.translation().z();
}

/// z with erfc(z / sqrt(2)) = alpha (bisection; as in test_monte_carlo.cpp, NIST DLMF, Eqs.
/// (7.2.1), (7.2.2)).
double z_of(double alpha) {
  double lo = 0.0;
  double hi = 10.0;
  for (int i = 0; i < 200; ++i) {
    const double mid = (lo + hi) / 2.0;
    (std::erfc(mid / std::numbers::sqrt2) > alpha ? lo : hi) = mid;
  }
  return lo;
}

MonteCarloResult run(const TrialContext& ctx, int threads, MonteCarloOptions options = {}) {
  oneapi::tbb::task_arena arena(threads);
  return arena.execute([&] { return rtt::tolerance::monte_carlo(ctx, options); });
}

/// Both empty or both equal bitwise.
void check_same(const std::optional<double>& a, const std::optional<double>& b) {
  REQUIRE(a.has_value() == b.has_value());
  if (a) CHECK(same(*a, *b));
}

/// Bitwise equality of two Monte Carlo runs: per trial the status, LM status, draws, criterion
/// values and compensators; per criterion mean, standard deviation and the yield with its
/// bounds; the total yield.
void check_bitwise(const MonteCarloResult& a, const MonteCarloResult& b) {
  REQUIRE(a.trials.size() == b.trials.size());
  for (std::size_t n = 0; n < a.trials.size(); ++n) {
    const auto& x = a.trials[n];
    const auto& y = b.trials[n];
    REQUIRE(x.trial == y.trial);
    CHECK(x.result.status == y.result.status);
    CHECK(x.result.lm_status == y.result.lm_status);
    CHECK(x.good == y.good);
    REQUIRE(x.values.size() == y.values.size());
    for (std::size_t i = 0; i < x.values.size(); ++i) CHECK(same(x.values[i], y.values[i]));
    REQUIRE(x.result.criteria.size() == y.result.criteria.size());
    for (std::size_t k = 0; k < x.result.criteria.size(); ++k) {
      const auto& u = x.result.criteria[k].value;
      const auto& v = y.result.criteria[k].value;
      REQUIRE(u.has_value() == v.has_value());
      if (u) CHECK(same(*u, *v));
    }
    REQUIRE(x.result.compensators.size() == y.result.compensators.size());
    for (std::size_t c = 0; c < x.result.compensators.size(); ++c) {
      CHECK(same(x.result.compensators[c].end, y.result.compensators[c].end));
    }
  }
  REQUIRE(a.criteria.size() == b.criteria.size());
  for (std::size_t k = 0; k < a.criteria.size(); ++k) {
    check_same(a.criteria[k].mean, b.criteria[k].mean);
    check_same(a.criteria[k].standard_deviation, b.criteria[k].standard_deviation);
    CHECK(a.criteria[k].yield.good == b.criteria[k].yield.good);
    check_same(a.criteria[k].yield.lower, b.criteria[k].yield.lower);
    check_same(a.criteria[k].yield.upper, b.criteria[k].yield.upper);
  }
  CHECK(a.yield.good == b.yield.good);
  check_same(a.yield.lower, b.yield.lower);
  check_same(a.yield.upper, b.yield.upper);
}

}  // namespace

TEST_CASE("m7 case 1: an assembly tilt moves all its children", "[m7]") {
  // assembly_motion: the assembly "module" holds the ideal lens f = 50 (ADR 0031) and, relative
  // to it, its detector in the focal plane; tilt theta about x at the pivot in the lens. The
  // axial bundle keeps its global direction z; in the frame of the assembly it runs along
  // R_x(-theta) z = (0, sin theta, cos theta), because R_x(+theta) turns +y toward +z
  // (docs/architecture.md, conventions; checked in the review of #254). ray_y is read in the
  // local frame of the detector, a child of the assembly, so in that frame: the ideal lens maps
  // the direction exactly onto f tan(theta) in its focal plane (case 11 for the sign).
  const System s = load("assembly_motion.rtt.json");
  const TrialContext module = context(s, items(s, {"module-tilt"}));
  volatile double theta_in = 1.0;  // volatile: no folding of tan (AGENTS.md, #218)
  const double theta_deg = theta_in;
  const double expected = 50.0 * std::tan(rtt::math::deg_to_rad(theta_deg));
  REQUIRE(nominal_value(module, "chief-y") == 0.0);
  const double y = value(module, {theta_deg}, "chief-y");
  INFO("ray_y " << y << " against f tan(theta) " << expected);
  CHECK(std::abs(y - expected) <= kEps);

  // The tilt of the assembly is the same tilt of each child about the same global pivot: the
  // pivots of lens and detector are the lens centre in their own coordinates (file overrides).
  // The detector is not in the subtree of the lens, so tilting the lens alone does not carry it.
  const TrialContext parts = context(s, items(s, {"lens-tilt", "image-tilt"}));
  CHECK(std::abs(value(parts, {theta_deg, theta_deg}, "chief-y") - y) <= kEps);
  // Content guard: the lens alone moves the image less than both together.
  CHECK(std::abs(value(parts, {theta_deg, 0.0}, "chief-y") - y) > 1e3 * kEps);

  // A node outside the assembly that is placed relative to its last surface stays at its
  // nominal place (ADR 0036, point 1b; compile::Perturbation).
  const rtt::tolerance::Applied a = rtt::tolerance::apply(module, std::vector<double>{theta_deg});
  const rtt::compile::Perturbation none;
  CHECK(std::abs(global_z(a.system, a.perturbation, "W.S1") - global_z(s, none, "W.S1")) <= kEps);
  CHECK(global_z(a.system, a.perturbation, "IMG") != global_z(s, none, "IMG"));  // content
}

TEST_CASE("m7 case 2: decenter of a thin lens moves the image by delta (1 - m)", "[m7]") {
  // Ideal lens f = 50 (ADR 0031), decenter delta along y. Object at infinity: the axial bundle
  // focuses on the lens axis, y' = +delta. Object at s = -150 (the lens's design conjugate):
  // s' = 75, m = s'/s = -0.5 (docs/quellen.md, row "Pupillen, Winkel- und Lateralvergrößerung");
  // the object point lies at -delta from the lens axis, its image at m (-delta) from it, so y' =
  // delta (1 - m) = 1.5 delta. The detector does not move.
  volatile double delta_in = 0.1;
  const double delta = delta_in;
  const System inf = load("decentred_lens.rtt.json");
  const TrialContext at_inf = context(inf, items(inf, {"lens-dy"}));
  CHECK(std::abs(value(at_inf, {delta}, "image-y") - delta) <= kEps);

  const System fin = load("decentred_lens_finite.rtt.json");
  const TrialContext at_fin = context(fin, items(fin, {"lens-dy"}));
  const double m = -0.5;
  CHECK(std::abs(value(at_fin, {delta}, "image-y") - delta * (1.0 - m)) <= kEps);

  // efl of the decentered trial from the parabasal rays around the virtual chief ray (ADR 0030,
  // addendum #256), so the nominal basis takes part. Hand value as in test_parabasal.cpp ("a
  // decentered ideal lens"), there for x, here for y: the chief ray crosses the lens plane at
  // the origin and runs to F = (0, delta, f), L = sqrt(f^2 + delta^2); the y pair (decenter
  // direction) gives k_y = f L^2 / (L^4 - delta^2 h^2) with h = 1e-3 * 5 mm (parabasal_step
  // times the pupil radius), the x pair k_x = 1/L; efl = 2 / (k_x + k_y), about
  // f + 3 delta^2 / (4 f) (1/L = (1 - delta^2/(2 f^2))/f, f/L^2 = (1 - delta^2/f^2)/f). Bound: the
  // Richardson rounding term of that test.
  const double f = 50.0;
  const double h = 1e-3 * 5.0;
  const double l = std::hypot(f, delta);
  const double ky = f * l * l / (l * l * l * l - delta * delta * h * h);
  const double kx = 1.0 / l;
  const double efl_expected = 2.0 / (kx + ky);
  const double r = 20.0 * std::numeric_limits<double>::epsilon() / (1e-3 * 5.0 / f) * f;
  const double efl = value(at_inf, {delta}, "efl");
  INFO("efl " << efl << " against " << efl_expected << ", bound " << r);
  CHECK(std::abs(efl - efl_expected) <= r);
  CHECK(std::abs(efl - f) > 1e3 * r);  // content guard: the decenter is in the value
}

TEST_CASE("m7 case 3: a tilted plane-parallel plate offsets the axial ray", "[m7]") {
  // Plate t = 10 mm, n = 1.5 (CONST in VACUUM), tilted by theta = 10 degree about x through its
  // centre. Snell's law at both faces (docs/quellen.md, row "Versatz der planparallelen Platte")
  // gives the lateral offset of the emerging parallel ray s = t sin(theta) [1 - cos(theta) /
  // sqrt(n^2 - sin^2(theta))]. Sign: the normal z of the plate turns to R_x(theta) z = (0, -sin
  // theta, cos theta); inside, the ray bends toward the normal, i.e. toward -y, so ray_y = -s on
  // the detector.
  const System s = load("tilted_plate.rtt.json");
  const TrialContext ctx = context(s, items(s, {"plate-tilt"}));
  volatile double theta_in = 10.0;  // degree
  const double theta = rtt::math::deg_to_rad(theta_in);
  const double t = 10.0;
  const double n = 1.5;
  const double offset =
      t * std::sin(theta) *
      (1.0 - std::cos(theta) / std::sqrt(n * n - std::sin(theta) * std::sin(theta)));
  REQUIRE(nominal_value(ctx, "chief-y") == 0.0);
  const double y = value(ctx, {static_cast<double>(theta_in)}, "chief-y");
  INFO("ray_y " << y << ", offset " << offset);
  CHECK(std::abs(y + offset) <= kEps);
}

TEST_CASE("m7 case 4: index and Abbe number end to end through first_order", "[m7]") {
  // Plano-convex N-BK7, curved face first, in VACUUM. In a trial efl comes from the parabasal
  // rays around the chief ray (ADR 0030, addendum #256; rtt-optim parabasal.cpp): the parallel
  // pair at the pupil offsets +-delta (delta = 1e-3 of the pupil radius 5 mm, h = 5e-3 mm) gives
  // phi = n' u' / y_1 with u' the tangent of the image-side angle, so efl = h / tan(beta(h)) of
  // the real ray at height h. That differs from the paraxial R / (n - 1) by the spherical
  // aberration at h (about 3e-7 mm here), so the target is the real ray itself: refraction at the
  // sphere and at the plane by Snell's law (docs/quellen.md, row "Brechung in Vektorform"): sin I =
  // h / R, sin I' = sin I / n, the ray runs at alpha = I - I' to the axis, and leaves the plane at
  // sin beta = n sin alpha. The thickness does not enter. Bound per value: the rounding term of
  // the Richardson bound, r = 20 eps / (delta u_m) |efl| with u_m = r_ep / |efl| (addendum #256,
  // as in rtt-optim/tests/test_parabasal.cpp), plus the aim tolerance of 1e-9 mm on h through the
  // h^2 dependence of the aberration term (2 |efl - R/(n-1)| 1e-9 / h, about 1.7e-13 mm).
  const System s = load("planoconvex_bk7.rtt.json");
  const TrialContext ctx = context(s, items(s, {"index", "abbe"}, schott()), schott());
  const double radius = 51.68;
  const double h = 1e-3 * 5.0;
  const auto efl_of = [&](double n) {
    const double sin_i = h / radius;
    const double alpha = std::asin(sin_i) - std::asin(sin_i / n);
    return h / std::tan(std::asin(n * std::sin(alpha)));
  };
  const auto bound = [&](double n) {
    const double efl = efl_of(n);
    const double r = 20.0 * std::numeric_limits<double>::epsilon() / (1e-3 * 5.0 / efl) * efl;
    return r + 2.0 * std::abs(efl - radius / (n - 1.0)) * 1e-9 / h;
  };
  const std::shared_ptr<const rtt::material::Material> glass = schott().resolve("SCHOTT:N-BK7");
  const double tc = s.environment.temperature_c;
  const double p = s.environment.pressure_atm;
  // The index at the wavelengths of the system, evaluated as compile does (compile.cpp). The
  // system's 0.5876 um is not the d line (air 0.5875618 um, in vacuum via
  // line_vacuum_wavelength_um), so x_d = n(0.5876) - n_d != 0.
  std::vector<double> n;
  for (const auto& w : s.wavelengths) n.push_back(glass->index(w.um, tc, p).real());
  REQUIRE(n.size() == 3);
  // n_d, n_F, n_C at the Fraunhofer lines, from the same function as the overlay (ADR 0036,
  // point 1c; AGENTS.md: the same function and input as the code).
  const rtt::compile::FraunhoferIndices fi = rtt::compile::fraunhofer_indices(*glass, tc, p);
  REQUIRE(fi.status == rtt::compile::FraunhoferIndices::Status::Ok);
  // ADR 0036, point 1c: D = n_F - n_C, nu_d = (n_d - 1)/D, D' = (n_d + Delta n_d - 1) /
  // (nu_d + Delta nu_d), n'(lambda) = n(lambda) + Delta n_d + (D'/D - 1)(n(lambda) - n_d).
  const double d = fi.nf - fi.nc;
  const double vd = (fi.nd - 1.0) / d;
  const auto perturbed = [&](double nl, double dn, double dv) {
    const double dp = (fi.nd + dn - 1.0) / (vd + dv);
    return nl + dn + (dp / d - 1.0) * (nl - fi.nd);
  };
  const char* ids[3] = {"efl-F", "efl-d", "efl-C"};
  for (std::size_t k = 0; k < 3; ++k) {
    const double v = nominal_value(ctx, ids[k]);
    INFO(ids[k] << " " << v << " against " << efl_of(n[k]));
    CHECK(std::abs(v - efl_of(n[k])) <= bound(n[k]));
    // Content guard: the aberration term (about 4.3e-7 mm, 49 times the bound) is in the value.
    CHECK(std::abs(v - radius / (n[k] - 1.0)) > 10.0 * bound(n[k]));
  }

  SECTION("Delta n_d") {
    // With Delta nu_d = 0, D'/D = (n_d + Delta n_d - 1)/(n_d - 1) != 1, so n' = n + Delta n_d (1 +
    // x / (n_d - 1)), x = n - n_d (not n + Delta n_d: x_d != 0 at the system's 0.5876 um).
    volatile double dn_in = 0.001;
    const double dn = dn_in;
    for (std::size_t k = 0; k < 3; ++k) {
      CHECK(std::abs(value(ctx, {dn, 0.0}, ids[k]) - efl_of(perturbed(n[k], dn, 0.0))) <=
            bound(perturbed(n[k], dn, 0.0)));
    }
    // First order in the effective change x = Delta n' / (n - 1) of the paraxial efl ratio
    // (n - 1)/(n' - 1) = 1/(1 + x): the series is alternating, the remainder at most x^2 (about
    // 3.7e-6). The code's ratio differs from the paraxial one by the change of the aberration
    // term (of the order 1e-11, from efl_of) and the rounding of both values.
    const double n1 = perturbed(n[1], dn, 0.0);
    const double ratio = value(ctx, {dn, 0.0}, "efl-d") / nominal_value(ctx, "efl-d");
    const double x = (n1 - n[1]) / (n[1] - 1.0);
    const double aberration = std::abs(efl_of(n1) / efl_of(n[1]) - (n[1] - 1.0) / (n1 - 1.0));
    const double rounding = (bound(n1) + bound(n[1])) / efl_of(n[1]);
    CHECK(std::abs((ratio - 1.0) + x) <= x * x + aberration + rounding);
    CHECK(std::abs((ratio - 1.0) + x) > 0.5 * x * x);  // content: the second order is there
  }

  SECTION("Delta nu_d") {
    volatile double dv_in = 1.0;
    const double dv = dv_in;
    std::vector<double> v(3);
    std::vector<double> v0(3);
    for (std::size_t k = 0; k < 3; ++k) {
      const double nk = perturbed(n[k], 0.0, dv);
      v[k] = value(ctx, {0.0, dv}, ids[k]);
      v0[k] = nominal_value(ctx, ids[k]);
      CHECK(std::abs(v[k] - efl_of(nk)) <= bound(nk));
    }
    // Linear in eps = D'/D - 1 (n' is exactly linear in eps): for the paraxial f = R / (n - 1),
    // d f = -f^2 / R * eps x with x = n - n_d, so the change of f_F - f_C is lin = -eps (f_F^2 x_F
    // - f_C^2 x_C) / R. The remainder per wavelength is exactly f t^2 / (1 + t), t = eps x / (n -
    // 1) (f / (1 + t) = f (1 - t + t^2 / (1 + t))). The code's values add the change of the
    // aberration term (from efl_of) and their rounding. Sign: a larger nu_d is less dispersion,
    // f_F - f_C (< 0) grows.
    const double eps = ((fi.nd - 1.0) / (vd + dv)) / d - 1.0;
    const double change = (v[0] - v[2]) - (v0[0] - v0[2]);
    double lin = 0.0;
    double remainder = 0.0;
    for (const std::size_t k : {std::size_t{0}, std::size_t{2}}) {
      const double f0 = radius / (n[k] - 1.0);
      const double x = n[k] - fi.nd;
      const double t = eps * x / (n[k] - 1.0);
      const double nk = perturbed(n[k], 0.0, dv);
      lin += (k == 0 ? -1.0 : 1.0) * eps * f0 * f0 * x / radius;
      remainder += f0 * t * t / (1.0 - std::abs(t));
      remainder += std::abs((efl_of(nk) - radius / (nk - 1.0)) - (efl_of(n[k]) - f0));
      remainder += bound(nk) + bound(n[k]);
    }
    INFO("change " << change << ", linear " << lin << ", bound " << remainder);
    CHECK(change > 0.0);
    CHECK(std::abs(change - lin) <= remainder);
  }
}

TEST_CASE("m7 case 5: the compensator takes back a shift of the mirror", "[m7]") {
  // Paraboloid R = -200 (f = 100), shifted by dz = +0.2 along z: its focus moves rigidly from
  // z = -100 to -99.8. Compensator: the detector position z (absolute), travel +-1 around the
  // perturbed value c0 = -100; merit: ray_y of the two rim rays, target 0. The paraboloid
  // focuses the axial bundle exactly, so r_i = J_i (c - c*) is exactly linear in c, c* = -99.8.
  //
  // Bound (derivation agreed with the coordinator): the travel is a bound on both sides, so LM
  // works in the MINUIT variable theta, c = a + (b - a)/2 (sin theta + 1) (rtt-optim
  // bounds.cpp), c' = (b - a)/2 cos theta. The Jacobian is the central difference in theta with
  // h = eps_f^(1/3) max(|theta|, 1) (levenberg_marquardt.cpp, Run::jacobian), exact for sin up to
  // s = sin(h)/h <= 1: J = J_c c' s. It is recomputed after every accepted step
  // (levenberg_marquardt.cpp, step 6, "a new Jacobian after an accepted step"); a rejected step
  // keeps theta and J. Scaling d = |J(theta_0)| (column norm of the start Jacobian, fixed). The
  // solve gives h = -J r / (J^2 + mu d^2), e_c = r / J_c, so |e_c| = |h| s (c'^2 + mu c'_0^2
  // s_0^2 / s^2) / |c'|. For |theta| <= 1 and theta_0 = 0 the step h_fd is the same, s_0 = s, and
  // since s <= 1 the factor is at most (c'^2 + mu c'_0^2) / |c'|. The run stops with
  // ConvergedStep, ||D h|| <= xtol (||D theta|| + xtol), the step h NOT applied
  // (levenberg_marquardt.cpp, step 3), so |e_c| <= xtol (|theta| + xtol/d) (c'^2 + mu c'_0^2) /
  // |c'|; or ConvergedGradient with gtol = 0, g = 0 and e_c = 0. About 2.9e-9 mm here.
  const System s = load("paraboloid_compensator.rtt.json");
  const TrialContext ctx = context(s, items(s, {"mirror-dz"}));
  volatile double dz_in = 0.2;
  const double dz = dz_in;
  const TrialResult r = rtt::tolerance::evaluate(ctx, std::vector<double>{dz});
  INFO("reason: " << r.reason);
  REQUIRE(r.status == TrialStatus::Valid);
  REQUIRE(r.compensators.size() == 1);
  REQUIRE(r.lm_status.has_value());
  REQUIRE(r.mu.has_value());
  const rtt::optim::LmStatus status = *r.lm_status;
  REQUIRE((status == rtt::optim::LmStatus::ConvergedStep ||
           status == rtt::optim::LmStatus::ConvergedGradient));
  const auto& comp = r.compensators[0];
  CHECK(!comp.at_bound);
  CHECK(!comp.infeasible);
  const double c0 = comp.start;
  REQUIRE(c0 == -100.0);
  const double target = -100.0 + dz;
  const double lo = c0 - 1.0;
  const double hi = c0 + 1.0;
  const double theta = std::asin(std::clamp(2.0 * (comp.end - lo) / (hi - lo) - 1.0, -1.0, 1.0));
  REQUIRE(std::abs(theta) <= 1.0);  // premise of s_0 = s above
  const double c0_prime = (hi - lo) / 2.0;
  const double c_prime = c0_prime * std::cos(theta);
  // |J_c|: each rim ray (height 5 at the mirror, sag 5^2 / (2 R) for the paraboloid) runs to the
  // focus with slope 5 / (distance); d >= 0.99 |J_c| c'_0, since s_0 >= 1 - h^2/6 with
  // h = cbrt(function_precision) = cbrt(1e-9) = 1e-3 (OptimizeOptions default kMeritPrecision,
  // which optimize.cpp passes to LmOptions; the compensator run uses TrialOptions::compensators).
  const double sag = 25.0 / (2.0 * -200.0);
  const double slope = 5.0 / ((sag + dz) - target);
  const double d = 0.99 * std::sqrt(2.0) * std::abs(slope) * c0_prime;
  const rtt::optim::OptimizeOptions o;
  const double bound = o.xtol * (std::abs(theta) + o.xtol / d) *
                       (c_prime * c_prime + *r.mu * c0_prime * c0_prime) / std::abs(c_prime);
  INFO("end " << comp.end << ", target " << target << ", mu " << *r.mu << ", bound " << bound);
  CHECK(std::abs(comp.end - target) <= bound);
  // The spot is back at its nominal value: a defocus e gives an RMS radius of at most |e| times
  // the largest slope (Burge, slide 8: the compensator restores the criterion).
  const double spot = r.criteria[criterion(ctx, "spot")].value.value_or(1.0);
  CHECK(spot <= nominal_value(ctx, "spot") + bound * std::abs(slope));
  // Content guard: without the compensator the shift defocuses.
  CHECK(std::abs(comp.end - comp.start) > 0.1);
}

TEST_CASE("m7 case 6: sensitivity and RSS against the Monte Carlo", "[m7]") {
  // decentred_lens_finite: m = -0.5. Linear: image-y = 1.5 delta_lens - delta_image (the
  // detector decenter moves the local origin), exactly, so the sensitivity deltas (per limit T,
  // ADR 0036, point 6) are 1.5 T_lens and -T_image, and predicted_sigma = sqrt(sum (r_i delta_i)^2)
  // with r_i = sigma_ratio (Burge, slides 5 and 9) is the true standard deviation.
  const System s = load("decentred_lens_finite.rtt.json");
  const TrialContext lin = context(s, items(s, {"lens-dy", "image-dy"}));
  rtt::tolerance::SensitivityOptions so;
  so.confirm = false;
  const rtt::tolerance::ToleranceSensitivity sens = rtt::tolerance::sensitivity(lin, so);
  const std::size_t ky = criterion(lin, "image-y");
  const auto& cy = sens.criteria[ky];
  REQUIRE(cy.items.size() == 2);
  REQUIRE(cy.items[0].delta.has_value());
  REQUIRE(cy.items[1].delta.has_value());
  CHECK(std::abs(*cy.items[0].delta - 1.5 * 0.1) <= kEps);  // 1 - m = 1.5
  CHECK(std::abs(*cy.items[1].delta + 0.05) <= kEps);       // the detector: -1
  const double ru = rtt::tolerance::sigma_ratio(ToleranceDistribution::Uniform);
  const double rn = rtt::tolerance::sigma_ratio(ToleranceDistribution::Normal);
  const double sigma = std::hypot(ru * 0.15, rn * 0.05);  // 0.08935
  REQUIRE(cy.predicted_sigma.has_value());
  CHECK(std::abs(*cy.predicted_sigma - sigma) <= 1e-12);

  // Monte Carlo, N = 1000 (file): Var(s^2) = (mu_4 - (N - 3)/(N - 1) sigma^4)/N with
  // mu_4 <= 3 sigma^4, a sum of independent uniform and truncated normal terms (each with a
  // kurtosis <= 3) has kurtosis <= 3; 4 standard deviations: s^2/sigma^2 in [0.821, 1.179].
  const MonteCarloResult mc = rtt::tolerance::monte_carlo(lin, {});
  REQUIRE(mc.trials.size() == 1000);
  const auto& st = mc.criteria[ky];
  REQUIRE(st.standard_deviation.has_value());
  const double n = 1000.0;
  const double band = 4.0 * std::sqrt((3.0 - (n - 3.0) / (n - 1.0)) / n);
  const double s2 = *st.standard_deviation * *st.standard_deviation / (sigma * sigma);
  INFO("s^2 / sigma^2 = " << s2 << ", band " << band);
  CHECK(std::abs(s2 - 1.0) <= band);

  // Quadratic: the detector shift dz alone. Behind the ideal lens every ray is straight, so at
  // the detector it lies at dz tan U_i from the image point: the RMS about the centroid is
  // Phi = k |dz| exactly, and the centroid stays at the chief ray (guard below). With dz uniform
  // on [-T, T]: q = Phi(T)^2 - Phi_0^2, predicted_mean_square = Phi_0^2 + q / 3, E Phi^2 = q/3
  // (+ Phi_0^2), Var(Phi^2) = k^4 T^4 (1/5 - 1/9) = (4/45) q^2; band 4 sqrt(4/45) q / sqrt(N),
  // 11.3 % of q/3 at N = 1000.
  const TrialContext quad = context(s, items(s, {"image-dz"}));
  const double up = value(quad, {0.1}, "rim-up");
  const double down = value(quad, {0.1}, "rim-down");
  CHECK(std::abs(up + down) <= kEps);  // guard: symmetric about the chief ray
  CHECK(std::abs(up) > 1e-3);          // content: the shift defocuses
  const rtt::tolerance::ToleranceSensitivity qs = rtt::tolerance::sensitivity(quad, so);
  const std::size_t kq = criterion(quad, "spot");
  const auto& cq = qs.criteria[kq];
  REQUIRE(cq.items.size() == 1);
  REQUIRE(cq.items[0].q.has_value());
  REQUIRE(cq.predicted_mean_square.has_value());
  const double q = *cq.items[0].q;
  const MonteCarloResult mq = rtt::tolerance::monte_carlo(quad, {});
  REQUIRE(mq.trials.size() == 1000);
  double sum = 0.0;
  for (const auto& t : mq.trials) {
    REQUIRE(t.result.criteria[kq].value.has_value());
    sum += *t.result.criteria[kq].value * *t.result.criteria[kq].value;
  }
  const double mean_sq = sum / n;
  const double qband = 4.0 * std::sqrt(4.0 / 45.0) * q / std::sqrt(n);
  INFO("mean Phi^2 " << mean_sq << ", predicted " << *cq.predicted_mean_square << ", band "
                     << qband);
  CHECK(std::abs(mean_sq - *cq.predicted_mean_square) <= qband);
}

TEST_CASE("m7 case 7: the Monte Carlo is reproducible, also with a compensator", "[m7]") {
  // paraboloid_compensator: two items, LM in every trial (P3-3 of the review of #301). Bitwise
  // the same for 1, 4 and all threads and on repetition (ADR 0004; ADR 0036, point 4).
  const System s = load("paraboloid_compensator.rtt.json");
  const TrialContext ctx = context(s, items(s, {"mirror-dz", "mirror-tilt"}));
  const MonteCarloResult one = run(ctx, 1);
  REQUIRE(one.trials.size() == 50);
  REQUIRE(std::all_of(one.trials.begin(), one.trials.end(), [](const auto& t) {
    return t.result.lm_status.has_value();
  }));  // content: the compensator ran in every trial
  check_bitwise(one, run(ctx, 4));
  check_bitwise(one, run(ctx, oneapi::tbb::task_arena::automatic));
  check_bitwise(one, run(ctx, oneapi::tbb::task_arena::automatic));

  // An item added at the end leaves the draws of the others unchanged (counter-based Philox,
  // ADR 0036, point 4).
  const TrialContext first = context(s, items(s, {"mirror-dz"}));
  const MonteCarloResult fewer = run(first, 1);
  REQUIRE(fewer.trials.size() == one.trials.size());
  for (std::size_t n = 0; n < one.trials.size(); ++n) {
    CHECK(same(fewer.trials[n].values[0], one.trials[n].values[0]));
  }
  // Content guard: another seed draws other values.
  MonteCarloOptions other;
  other.seed = 26206;
  const MonteCarloResult m = run(ctx, 1, other);
  CHECK(!same(m.trials[0].values[0], one.trials[0].values[0]));
}

TEST_CASE("m7 case 8: statistics against the analytic distribution", "[m7]") {
  // decentred_lens: image-y = delta exactly (case 2). N = 200 (file).
  const System s = load("decentred_lens.rtt.json");
  const double t = 0.1;
  const double n = 200.0;
  MonteCarloOptions o;
  o.z = z_of(1e-4);
  REQUIRE(std::abs(o.z - 3.8906) < 1e-4);

  SECTION("uniform") {
    // sigma = T/sqrt(3); mean within 4 sigma/sqrt(N); mu_4 = 9/5 sigma^4, so Var(s^2)/sigma^4 =
    // (9/5 - (N-3)/(N-1))/N = 4.05e-3, s^2/sigma^2 within +-0.2546; yield of |y| <= L = 0.05 is
    // L/T = 0.5 and lies in the Wilson interval (z = 3.8906).
    const TrialContext ctx = context(s, items(s, {"lens-dy"}));
    const MonteCarloResult mc = rtt::tolerance::monte_carlo(ctx, o);
    REQUIRE(mc.trials.size() == 200);
    const auto& st = mc.criteria[criterion(ctx, "image-y")];
    REQUIRE(st.mean.has_value());
    REQUIRE(st.standard_deviation.has_value());
    const double sigma = t / std::sqrt(3.0);
    CHECK(std::abs(*st.mean) <= 4.0 * sigma / std::sqrt(n));
    const double band = 4.0 * std::sqrt((9.0 / 5.0 - (n - 3.0) / (n - 1.0)) / n);
    const double s2 = *st.standard_deviation * *st.standard_deviation / (sigma * sigma);
    INFO("s^2/sigma^2 " << s2 << ", band " << band);
    CHECK(std::abs(s2 - 1.0) <= band);
    REQUIRE(st.yield.lower.has_value());
    REQUIRE(st.yield.upper.has_value());
    INFO("yield [" << *st.yield.lower << ", " << *st.yield.upper << "]");
    CHECK(*st.yield.lower <= 0.5);
    CHECK(0.5 <= *st.yield.upper);
  }

  SECTION("truncated normal") {
    // sigma = T/2, truncated at +-T = +-2 sigma: Var = sigma^2 (1 - 2 k phi(k) / (2 Phi(k) - 1))
    // with k = 2 (0.773741 sigma^2; sigma_T / T = 0.439813 = sigma_ratio(Normal)); its kurtosis
    // is below 3, so mu_4 <= 3 sigma_T^4 in the band.
    std::vector<ToleranceItem> it = items(s, {"lens-dy"});
    it[0].distribution = ToleranceDistribution::Normal;
    const TrialContext ctx = context(s, it);
    const MonteCarloResult mc = rtt::tolerance::monte_carlo(ctx, o);
    const auto& st = mc.criteria[criterion(ctx, "image-y")];
    REQUIRE(st.standard_deviation.has_value());
    const double k = 2.0;
    const double phi = std::exp(-k * k / 2.0) / std::sqrt(2.0 * std::numbers::pi);
    const double big_phi = std::erfc(-k / std::numbers::sqrt2) / 2.0;
    const double var = (t / 2.0) * (t / 2.0) * (1.0 - 2.0 * k * phi / (2.0 * big_phi - 1.0));
    REQUIRE(std::abs(std::sqrt(var) / t -
                     rtt::tolerance::sigma_ratio(ToleranceDistribution::Normal)) <= 1e-12);
    const double band = 4.0 * std::sqrt((3.0 - (n - 3.0) / (n - 1.0)) / n);
    const double s2 = *st.standard_deviation * *st.standard_deviation / var;
    INFO("s^2/sigma^2 " << s2 << ", band " << band);
    CHECK(std::abs(s2 - 1.0) <= band);
  }
}

TEST_CASE("m7 case 9: a tolerance on a row of the parameter table", "[m7]") {
  // configuration_row: the detector lies T behind the ideal lens f = 50, T = 50 (near) and 55
  // (far). The rim ray (5 mm) runs behind the lens with tan U' = -5/50 in both configurations, so
  // an offset Delta T of the row moves it by Delta T tan U' = -0.01 mm in each (straight rays).
  const System s = load("configuration_row.rtt.json");
  const TrialContext ctx = context(s, items(s, {"T"}));
  volatile double dt_in = 0.1;
  const double dt = dt_in;
  const double tan_u = -5.0 / 50.0;
  for (const char* id : {"rim-near", "rim-far"}) {
    INFO(id);
    CHECK(std::abs(value(ctx, {dt}, id) - nominal_value(ctx, id) - dt * tan_u) <= kEps);
  }
  // Content guard: the two configurations differ (0 and -0.5 mm).
  CHECK(std::abs(nominal_value(ctx, "rim-near")) <= kEps);
  CHECK(std::abs(nominal_value(ctx, "rim-far") + 0.5) <= kEps);
}

TEST_CASE("m7 case 10: double pass, a mirror tilt deviates by 2 theta", "[m7]") {
  // double_pass: the axial ray passes the plate, meets the flat mirror 100 mm behind P.S2 and
  // returns. Mirror tilt theta about x at its vertex: the normal -z turns to (0, sin theta,
  // -cos theta), d' = d - 2 (d.n) n (de Greve, Eq. (13); docs/quellen.md) = (0, sin 2theta, -cos
  // 2theta); after 100 mm back along -z the ray is at y = 100 tan(2 theta) on P.S2 (second pass).
  const System s = load("double_pass.rtt.json");
  const TrialContext tilt = context(s, items(s, {"mirror-tilt"}));
  volatile double theta_in = 0.5;
  const double theta = rtt::math::deg_to_rad(theta_in);
  const double y = value(tilt, {static_cast<double>(theta_in)}, "return-y");
  INFO("y " << y << " against " << 100.0 * std::tan(2.0 * theta));
  CHECK(std::abs(y - 100.0 * std::tan(2.0 * theta)) <= kEps);

  // The element met twice sees the same perturbation in both passes: decenter delta of the
  // plate moves the local hit of the axial ray to -delta on P.S1, forward and back.
  const TrialContext plate = context(s, items(s, {"plate-dy"}));
  volatile double delta_in = 0.1;
  const double delta = delta_in;
  CHECK(std::abs(value(plate, {delta}, "plate-in") + delta) <= kEps);
  CHECK(std::abs(value(plate, {delta}, "plate-out") + delta) <= kEps);
}

TEST_CASE("m7 case 11: the sign of the assembly tilt", "[m7]") {
  // As case 1: ray_y = +f tan(theta) in the frame of the assembly; -theta gives the mirror image.
  const System s = load("assembly_motion.rtt.json");
  const TrialContext module = context(s, items(s, {"module-tilt"}));
  volatile double theta_in = 1.0;
  const double theta_deg = theta_in;
  const double plus = value(module, {theta_deg}, "chief-y");
  const double minus = value(module, {-theta_deg}, "chief-y");
  CHECK(plus > 0.0);
  CHECK(std::abs(plus - 50.0 * std::tan(rtt::math::deg_to_rad(theta_deg))) <= kEps);
  CHECK(std::abs(minus + plus) <= kEps);
}

TEST_CASE("m7 case 12: relative placement follows a thickness, not a rigid motion", "[m7]") {
  // L3 (#252): a thickness offset changes the model (ADR 0036, point 1a), so the detector placed
  // relative to the lens and the witness placed relative to the detector move with it. A rigid
  // motion of the assembly is an overlay (point 1b): its subtree moves, the witness outside it
  // stays, although it is placed relative to the assembly's last surface.
  const System s = load("assembly_motion.rtt.json");
  const rtt::compile::Perturbation none;
  const double il = global_z(s, none, "IL");
  const double img = global_z(s, none, "IMG");
  const double w = global_z(s, none, "W.S1");
  volatile double dt_in = 0.05;
  const double dt = dt_in;
  const TrialContext thick = context(s, items(s, {"lens-thickness"}));
  const rtt::tolerance::Applied a = rtt::tolerance::apply(thick, std::vector<double>{dt});
  CHECK(std::abs(global_z(a.system, a.perturbation, "IL") - il - dt) <= kEps);
  CHECK(std::abs(global_z(a.system, a.perturbation, "IMG") - img - dt) <= kEps);
  CHECK(std::abs(global_z(a.system, a.perturbation, "W.S1") - w - dt) <= kEps);

  volatile double dz_in = 0.1;
  const double dz = dz_in;
  const TrialContext shift = context(s, items(s, {"module-shift"}));
  const rtt::tolerance::Applied b = rtt::tolerance::apply(shift, std::vector<double>{dz});
  CHECK(std::abs(global_z(b.system, b.perturbation, "IL") - il - dz) <= kEps);
  CHECK(std::abs(global_z(b.system, b.perturbation, "IMG") - img - dz) <= kEps);
  CHECK(std::abs(global_z(b.system, b.perturbation, "W.S1") - w) <= kEps);
}

TEST_CASE("m7 case 13: tolerance classes with overrides and explicit items", "[m7]") {
  // achromat_tolerances: class normal, the assembly "objective" precise, a single value for
  // L1.S2 curvature_relative, the image without class, three explicit items. Precedence explicit
  // > override > class (ADR 0035, point 5), deterministic order, stable ids; the Monte Carlo over
  // the generated items is bitwise the same on repetition.
  const System s = load("achromat_tolerances.rtt.json");
  const rtt::tolerance::ToleranceItems a = rtt::tolerance::expand(s, schott());
  const rtt::tolerance::ToleranceItems b = rtt::tolerance::expand(s, schott());
  REQUIRE(a.items.size() == b.items.size());
  for (std::size_t i = 0; i < a.items.size(); ++i) CHECK(a.items[i] == b.items[i]);

  std::size_t precise = 0;
  std::size_t single_values = 0;
  for (std::size_t i = 0; i < a.items.size(); ++i) {
    const ToleranceItem& item = a.items[i];
    const rtt::tolerance::ItemOrigin& o = a.origins[i];
    if (o.source != rtt::tolerance::ItemOrigin::Source::Class) continue;
    const auto& c = a.classes.items[o.index];
    if (const std::optional<std::size_t> v = c.value_override) {
      // The single value of the file replaces the class limit (two here: L1.S2
      // curvature_relative 0.005 and L1 segment 1 index 0.001).
      CHECK(item.limit ==
            std::get<rtt::model::ToleranceValueOverride>(s.tolerances.overrides[*v]).limit);
      ++single_values;
      continue;
    }
    if (c.kind == rtt::model::ToleranceKind::Abbe) continue;  // relative to nu_d (#260)
    CHECK(item.limit == rtt::tolerance::class_limit(c.tolerance_class, c.kind, c.level));
    if (c.tolerance_class == rtt::model::ToleranceClass::Precise) ++precise;
  }
  CHECK(single_values == 2);
  CHECK(precise > 0);            // the override of the objective took effect
  REQUIRE(a.items.size() >= 3);  // before the iterator arithmetic below
  // Explicit before class: the class entry L1.S1 conic is replaced by the explicit item
  // S1-conic (index 1 of the file's items), and no generated item perturbs the same quantity.
  const auto replaced =
      std::find_if(a.classes.items.begin(), a.classes.items.end(),
                   [](const rtt::model::ClassTolerance& c) { return c.replaced_by.has_value(); });
  REQUIRE(replaced != a.classes.items.end());
  REQUIRE(*replaced->replaced_by < s.tolerances.items.size());
  CHECK(s.tolerances.items[*replaced->replaced_by].id == "S1-conic");
  CHECK(std::none_of(a.items.begin(), a.items.end() - 3, [&](const ToleranceItem& i) {
    return i.target == replaced->target && i.kind == replaced->kind;
  }));
  // The explicit items follow the generated ones, in file order.
  CHECK(a.items[a.items.size() - 3].id == "S3-wedge");
  CHECK(a.items[a.items.size() - 2].id == "S1-conic");
  CHECK(a.items[a.items.size() - 1].id == "R1");

  const TrialContext ctx = rtt::tolerance::prepare(s, schott(), nullptr, a.items);
  MonteCarloOptions o;
  o.trials = 10;
  check_bitwise(run(ctx, oneapi::tbb::task_arena::automatic, o),
                run(ctx, oneapi::tbb::task_arena::automatic, o));
}
