// Acceptance of M5 (#170): the optimizer against independent reference values, tolerances derived
// beforehand. Case 2 here: tests/reference/m5/two_lens_gap.rtt.json, a variable air gap D (row of
// the parameter table, bound to the z position of L2) with L2 and the image placed
// relative_to_preceding (ADR 0028, 0029); merit: EFL = 55 mm and the marginal ray (py = 1) at
// the image height 0 (two equations, two unknowns D and the image distance).

#include <tbb/global_control.h>
#include <tbb/info.h>

#include <algorithm>
#include <array>
#include <bit>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <iomanip>
#include <limits>
#include <optional>
#include <span>
#include <sstream>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
#include "rtt/model/parameters.hpp"
#include "rtt/optim/bounds.hpp"
#include "rtt/optim/merit.hpp"
#include "rtt/optim/optimize.hpp"
#include "rtt/optim/variables.hpp"
#include "rtt/paraxial/paraxial.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/run_control.hpp"
#include "rtt/trace/sequential.hpp"
#include "rtt/trace/sources.hpp"

using rtt::compile::CompiledSystem;
using rtt::material::MaterialLibrary;
using rtt::model::System;
using rtt::optim::LmStatus;
using rtt::optim::OptimizeOptions;
using rtt::optim::OptimResult;

namespace {

constexpr double kEps = std::numeric_limits<double>::epsilon();
constexpr double kTargetEfl = 55.0;  // the efl operand of two_lens_gap.rtt.json

System load(const std::string& relative) {
  return rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/" + relative);
}

double z_of(const CompiledSystem& cs, const char* id) {
  return cs.surfaces()[*cs.find_surface(rtt::model::SurfaceId{id})].to_global.translation().z();
}

/// The system with only the children `keep` of the root (by index), without merit function.
System only(const System& s, std::initializer_list<std::size_t> keep) {
  System out = s;
  out.optimization = {};
  out.root.children.clear();
  for (const std::size_t i : keep) out.root.children.push_back(s.root.children[i]);
  return out;
}

/// Run options of the acceptance: the step test ends the run (ftol = 0, gtol = 0 by default), and
/// xtol = 1e-12 makes its bound small against 1e-10 of the EFL (derivation below). The default
/// difference step (kMeritPrecision = 1e-9, h = 1e-3 max(|theta|, 1)) is kept: this is a
/// zero-residual problem, so the Jacobian error does not move the solution (J~^T f = 0 at f = 0),
/// it only slows the convergence.
OptimizeOptions acceptance_options() {
  OptimizeOptions o;
  o.ftol = 0.0;
  o.xtol = 1e-12;
  return o;
}

/// Jacobian of the residuals with respect to the internal variables theta (bounds.hpp) at the
/// external values p, by central differences with a small step (only for the tolerance bound:
/// magnitudes, not reference values). Rows: residuals, columns: variables.
std::vector<std::array<double, 2>> jacobian_theta(const rtt::optim::MeritFunction& merit,
                                                  std::span<const double> p) {
  const auto& vars = merit.variables();
  REQUIRE(vars.size() == 2);
  std::vector<std::array<double, 2>> j(merit.size());
  for (std::size_t c = 0; c < 2; ++c) {
    const rtt::optim::Bounds& b = vars[c].bounds;
    const double theta = rtt::optim::to_internal(p[c], b);
    const double h = 1e-6 * std::max(std::abs(theta), 1.0);
    std::vector<double> plus(p.begin(), p.end());
    std::vector<double> minus(p.begin(), p.end());
    plus[c] = rtt::optim::to_external(theta + h, b);
    minus[c] = rtt::optim::to_external(theta - h, b);
    const rtt::optim::MeritEvaluation fp = merit.evaluate(plus);
    const rtt::optim::MeritEvaluation fm = merit.evaluate(minus);
    REQUIRE(fp.valid());
    REQUIRE(fm.valid());
    for (std::size_t i = 0; i < merit.size(); ++i) {
      j[i][c] = (fp.residuals[i] - fm.residuals[i]) / (2.0 * h);
    }
  }
  return j;
}

// --- Case 1: the independent 1-D reference for m5/singlet_optim ---------------------------------

/// The singlet of m5/singlet_optim with R1 given and R2 from the EFL condition, and the best
/// image plane for the RMS spot about the centroid of the Gauss rays (rms_spot of the file:
/// GaussPupil{3, 6}, field 0, reference wavelength, quadrature weights).
struct Bending {
  double r1 = 0.0;
  double r2 = 0.0;
  double z_image = 0.0;  ///< global z of the best image plane, mm (= the detector position)
  double s = 0.0;        ///< RMS^2 about the centroid there, mm^2 (the generator's sum of squares)
};

/// R2 of a lens of index n and centre thickness d in a medium n_a for the power phi (EFL 1/phi):
/// the lens as two components (Greivenkamp, OPTI-201/202, Sec. 7 p. 7-4, docs/quellen.md):
/// phi1 = (n - n_a) c1, phi2 = (n_a - n) c2 at the vertices (their principal planes), tau = d / n,
/// phi = phi1 + phi2 - phi1 phi2 tau, so phi2 = (phi - phi1) / (1 - tau phi1).
double thick_r2(double r1, double n, double n_a, double d, double phi) {
  const double phi1 = (n - n_a) / r1;
  const double phi2 = (phi - phi1) / (1.0 - d / n * phi1);
  return (n_a - n) / phi2;
}

Bending bending(const System& base, const MaterialLibrary& lib, double r1) {
  System s = base;
  s.optimization = {};
  auto& lens = std::get<rtt::model::Element>(s.root.children[1].value);
  const CompiledSystem probe = rtt::compile::compile(s, lib);
  const std::uint16_t ref = probe.reference_wavelength();
  const double n_a = probe.media()[probe.environment_medium()].index[ref].real();
  const double d = z_of(probe, "L1.S2") - z_of(probe, "L1.S1");
  constexpr double kN = 1.5168;  // CONST:1.5168 of the file
  Bending b;
  b.r1 = r1;
  b.r2 = thick_r2(r1, kN, n_a, d, 1.0 / 100.0);
  std::get<rtt::model::Conic>(lens.surfaces[0].shape.base).radius.value = b.r1;
  std::get<rtt::model::Conic>(lens.surfaces[1].shape.base).radius.value = b.r2;
  const CompiledSystem cs = rtt::compile::compile(s, lib);
  const rtt::compile::PathId main = *cs.find_path("main");
  const rtt::trace::GaussPupil gauss{3, 6};
  const std::uint16_t field = 0;
  rtt::trace::RayBatch rays = rtt::trace::make_rays(cs, main, std::span(&field, 1), ref, gauss);
  static_cast<void>(rtt::trace::SequentialTracer().trace(cs, main, rays));
  const std::vector<double> q = rtt::trace::gauss_pupil_weights(gauss);
  REQUIRE(q.size() == rays.size());
  // After the last surface every ray is straight: at the plane z0 + s the point is a + s b with
  // a = (x, y) on the detector (z0) and b = (dx/dz, dy/dz). With weights q (sum 1) the centroid
  // is a_c + s b_c, and RMS^2(s) = sum q |(a - a_c) + s (b - b_c)|^2, a parabola in s with its
  // minimum at s* = -sum q (a - a_c).(b - b_c) / sum q |b - b_c|^2.
  const double z0 = z_of(cs, "IMG");
  std::array<double, 2> ac{};
  std::array<double, 2> bc{};
  for (std::size_t k = 0; k < rays.size(); ++k) {
    REQUIRE(rays.status()[k] == rtt::trace::RayStatus::Alive);
    ac[0] += q[k] * rays.pos_x()[k];
    ac[1] += q[k] * rays.pos_y()[k];
    bc[0] += q[k] * rays.dir_x()[k] / rays.dir_z()[k];
    bc[1] += q[k] * rays.dir_y()[k] / rays.dir_z()[k];
  }
  double ab = 0.0;
  double bb = 0.0;
  for (std::size_t k = 0; k < rays.size(); ++k) {
    const double ax = rays.pos_x()[k] - ac[0];
    const double ay = rays.pos_y()[k] - ac[1];
    const double bx = rays.dir_x()[k] / rays.dir_z()[k] - bc[0];
    const double by = rays.dir_y()[k] / rays.dir_z()[k] - bc[1];
    ab += q[k] * (ax * bx + ay * by);
    bb += q[k] * (bx * bx + by * by);
  }
  REQUIRE(bb > 0.0);
  const double s_star = -ab / bb;
  b.z_image = z0 + s_star;
  for (std::size_t k = 0; k < rays.size(); ++k) {
    const double x = rays.pos_x()[k] - ac[0] + s_star * (rays.dir_x()[k] / rays.dir_z()[k] - bc[0]);
    const double y = rays.pos_y()[k] - ac[1] + s_star * (rays.dir_y()[k] / rays.dir_z()[k] - bc[1]);
    b.s += q[k] * (x * x + y * y);
  }
  return b;
}

/// Golden-section search for the minimum of S(R1) on [lo, hi] (bracket from a coarse grid), until
/// the bracket is narrower than `width`.
Bending golden_section(
    const System& base, const MaterialLibrary& lib, double lo, double hi, double width) {
  const double g = 0.5 * (std::sqrt(5.0) - 1.0);
  double a = lo;
  double b = hi;
  Bending c = bending(base, lib, b - g * (b - a));
  Bending d = bending(base, lib, a + g * (b - a));
  while (b - a > width) {
    if (c.s < d.s) {
      b = d.r1;
      d = c;
      c = bending(base, lib, b - g * (b - a));
    } else {
      a = c.r1;
      c = d;
      d = bending(base, lib, a + g * (b - a));
    }
  }
  return c.s < d.s ? c : d;
}

}  // namespace

TEST_CASE("M5 acceptance, case 1: the 1-D bending reference and Coddington's shape factor",
          "[optim][m5]") {
  // The reference of case 1, independent of rtt-optim: R2(R1) keeps the EFL at 100 mm, the image
  // plane is the best one for each R1, and a golden-section search over R1 finds the minimum of
  // the RMS spot. Checked here: the EFL of the reference system (thick-lens relation against the
  // paraxial engine) and, as plausibility only, the shape factor of minimal spherical aberration
  // of a thin lens (Sasian, OPTI 518 (2019), Lecture 14: p. 9 X = (c1 + c2) / (c1 - c2), p. 15
  // object at infinity Y = 1, p. 30 minimum of sigma_I at X = 2 (n + 1)(n - 1) / (n + 2) Y with the
  // relative index n; docs/quellen.md). The stop 5 mm before the lens does not change S_I (stop
  // shift, S*_I = S_I). Band +-0.1, set beforehand: thick lens (d = 4 mm), RMS spot with real rays
  // instead of the third-order S_I.
  const System s = load("m5/singlet_optim.rtt.json");
  const MaterialLibrary lib;
  // Coarse grid for the bracket: R1 from 40 to 100 mm in steps of 5 (the thin-lens minimum lies
  // near R1 = 59 mm).
  double best_r1 = 40.0;
  double best_s = std::numeric_limits<double>::infinity();
  for (double r1 = 40.0; r1 <= 100.0; r1 += 5.0) {
    const Bending b = bending(s, lib, r1);
    if (b.s < best_s) {
      best_s = b.s;
      best_r1 = r1;
    }
  }
  REQUIRE(best_r1 > 40.0);
  REQUIRE(best_r1 < 100.0);  // an interior minimum
  const Bending ref = golden_section(s, lib, best_r1 - 5.0, best_r1 + 5.0, 1e-6);
  INFO("R1* = " << ref.r1 << ", R2* = " << ref.r2 << ", z* = " << ref.z_image
                << ", RMS* = " << std::sqrt(ref.s) << " mm");

  // The EFL of the reference system from the paraxial engine.
  System at_ref = s;
  auto& lens = std::get<rtt::model::Element>(at_ref.root.children[1].value);
  std::get<rtt::model::Conic>(lens.surfaces[0].shape.base).radius.value = ref.r1;
  std::get<rtt::model::Conic>(lens.surfaces[1].shape.base).radius.value = ref.r2;
  const CompiledSystem cs = rtt::compile::compile(at_ref, lib);
  const std::uint16_t wl = cs.reference_wavelength();
  CHECK(std::abs(*rtt::paraxial::first_order(cs, *cs.find_path("main"), wl).efl - 100.0) <=
        1e-12 * 100.0);

  // Plausibility: Coddington's shape factor.
  const double n_rel = 1.5168 / cs.media()[cs.environment_medium()].index[wl].real();
  const double x_thin = 2.0 * (n_rel + 1.0) * (n_rel - 1.0) / (n_rel + 2.0);
  const double x_ref = (1.0 / ref.r1 + 1.0 / ref.r2) / (1.0 / ref.r1 - 1.0 / ref.r2);
  INFO("X* = " << x_ref << ", thin lens X = " << x_thin);
  CHECK(std::abs(x_ref - x_thin) <= 0.1);
  CHECK(ref.s > 0.0);  // content guard: a real spot, not a perfect image
}

namespace {

// --- Small dense helpers for case 1 (3 variables) ------------------------------------------------

using Vec = std::vector<double>;
using Mat = std::vector<std::vector<double>>;  // rows

Vec residuals_at(const rtt::optim::MeritFunction& merit, const Vec& p) {
  const rtt::optim::MeritEvaluation e = merit.evaluate(p);
  REQUIRE(e.valid());
  return e.residuals;
}

/// Jacobian (m x n) of the residuals at p by central differences with the steps `h` (no bounds:
/// theta = p).
Mat jacobian_at(const rtt::optim::MeritFunction& merit, const Vec& p, const Vec& h) {
  const std::size_t n = p.size();
  Mat j(merit.size(), Vec(n, 0.0));
  for (std::size_t c = 0; c < n; ++c) {
    Vec plus = p;
    Vec minus = p;
    plus[c] += h[c];
    minus[c] -= h[c];
    const Vec fp = residuals_at(merit, plus);
    const Vec fm = residuals_at(merit, minus);
    const double denominator = plus[c] - minus[c];
    for (std::size_t i = 0; i < j.size(); ++i) j[i][c] = (fp[i] - fm[i]) / denominator;
  }
  return j;
}

Vec steps(const Vec& p, double relative) {
  Vec h;
  for (const double v : p) h.push_back(relative * std::max(std::abs(v), 1.0));
  return h;
}

/// J^T f.
Vec gradient(const Mat& j, const Vec& f) {
  Vec g(j[0].size(), 0.0);
  for (std::size_t i = 0; i < j.size(); ++i) {
    for (std::size_t c = 0; c < g.size(); ++c) g[c] += j[i][c] * f[i];
  }
  return g;
}

/// Solves the small system a x = b (Gaussian elimination with partial pivoting; the test's own
/// solver for n <= 3).
Vec solve_n(Mat a, Vec b) {
  const std::size_t n = b.size();
  for (std::size_t c = 0; c < n; ++c) {
    std::size_t pivot = c;
    for (std::size_t r = c + 1; r < n; ++r) {
      if (std::abs(a[r][c]) > std::abs(a[pivot][c])) pivot = r;
    }
    REQUIRE(a[pivot][c] != 0.0);
    std::swap(a[c], a[pivot]);
    std::swap(b[c], b[pivot]);
    for (std::size_t r = c + 1; r < n; ++r) {
      const double factor = a[r][c] / a[c][c];
      for (std::size_t k = c; k < n; ++k) a[r][k] -= factor * a[c][k];
      b[r] -= factor * b[c];
    }
  }
  Vec x(n, 0.0);
  for (std::size_t c = n; c-- > 0;) {
    double s = b[c];
    for (std::size_t k = c + 1; k < n; ++k) s -= a[c][k] * x[k];
    x[c] = s / a[c][c];
  }
  return x;
}

double frobenius(const Mat& a) {
  double s = 0.0;
  for (const Vec& row : a) {
    for (const double v : row) s += v * v;
  }
  return std::sqrt(s);
}

/// Hessian of F = 1/2 f^T f at p: J^T J + sum_i f_i d2f_i, the second part from central
/// differences of J contracted with f (so the large EFL row, whose residual is about 0 at the
/// reference, does not enter through second differences of F).
Mat hessian_at(const rtt::optim::MeritFunction& merit,
               const Vec& p,
               const Mat& j,
               const Vec& f,
               const Vec& h_j,
               const Vec& h_h) {
  const std::size_t n = p.size();
  Mat a(n, Vec(n, 0.0));
  for (std::size_t r = 0; r < n; ++r) {
    for (std::size_t c = 0; c < n; ++c) {
      for (std::size_t i = 0; i < j.size(); ++i) a[r][c] += j[i][r] * j[i][c];
    }
  }
  for (std::size_t k = 0; k < n; ++k) {
    Vec plus = p;
    Vec minus = p;
    plus[k] += h_h[k];
    minus[k] -= h_h[k];
    const Mat jp = jacobian_at(merit, plus, h_j);
    const Mat jm = jacobian_at(merit, minus, h_j);
    for (std::size_t r = 0; r < n; ++r) {
      double s = 0.0;
      for (std::size_t i = 0; i < j.size(); ++i) {
        s += f[i] * (jp[i][r] - jm[i][r]) / (plus[k] - minus[k]);
      }
      a[r][k] += s;
    }
  }
  // Symmetrize the difference quotient.
  for (std::size_t r = 0; r < n; ++r) {
    for (std::size_t c = r + 1; c < n; ++c) {
      const double m = 0.5 * (a[r][c] + a[c][r]);
      a[r][c] = m;
      a[c][r] = m;
    }
  }
  return a;
}

/// The 1-D bending reference of m5/singlet_optim (radii free in the model): a coarse grid over
/// R1 for the bracket, then the golden section.
Bending bending_reference(const System& base, const MaterialLibrary& lib) {
  double best_r1 = 40.0;
  double best_s = std::numeric_limits<double>::infinity();
  for (double r1 = 40.0; r1 <= 100.0; r1 += 5.0) {
    const Bending b = bending(base, lib, r1);
    if (b.s < best_s) {
      best_s = b.s;
      best_r1 = r1;
    }
  }
  return golden_section(base, lib, best_r1 - 5.0, best_r1 + 5.0, 1e-6);
}

/// Runs optimize() on `s` (ftol = 0, default xtol) and compares the result with the least-squares
/// optimum predicted from the reference point p_ref, with the tolerance derived before the run
/// (see the test of case 1). `w_efl`: weight of the efl operand in row 0 if it enters the merit
/// (then its predicted residual must lie below 1e-9 mm); none if it is only an observer.
/// `golden_width0`: for a reference from the golden section over variable 0 with the other
/// variable (n = 2) in closed form for it, the width of the final bracket in the units of
/// variable 0. Then the reference is checked to lie at the optimum of the merit function
/// before the run (review of #197, P2): |delta_n| against a bound from that width and the
/// rounding of a flat minimum, and |delta_2| << |delta_n| for the Newton prediction.
/// Returns the result for further checks.
OptimResult check_against_prediction(const System& s,
                                     const MaterialLibrary& lib,
                                     const Vec& p_ref,
                                     std::optional<double> w_efl,
                                     std::optional<double> golden_width0 = std::nullopt) {
  const std::size_t n = p_ref.size();
  const rtt::optim::MeritFunction merit(s, lib, nullptr);
  const auto& vars = merit.variables();
  REQUIRE(vars.size() == n);
  for (const auto& v : vars) REQUIRE(!v.bounds.min.has_value());
  // Newton correction from the reference.
  const Vec h_acc_ref = steps(p_ref, 1e-6);
  const Vec f_ref = residuals_at(merit, p_ref);
  const Mat j_ref = jacobian_at(merit, p_ref, h_acc_ref);
  const Mat hess = hessian_at(merit, p_ref, j_ref, f_ref, h_acc_ref, steps(p_ref, 1e-4));
  const Vec delta_n = solve_n(hess, gradient(j_ref, f_ref));
  Vec p_pred = p_ref;
  for (std::size_t c = 0; c < n; ++c) p_pred[c] -= delta_n[c];
  const Vec f_pred = residuals_at(merit, p_pred);
  const Mat j_pred = jacobian_at(merit, p_pred, steps(p_pred, 1e-6));
  const Vec delta_2 = solve_n(hess, gradient(j_pred, f_pred));
  // Solver's difference step and its bias.
  const double h_rel = std::cbrt(rtt::optim::kMeritPrecision);
  const Mat j_solver = jacobian_at(merit, p_pred, steps(p_pred, h_rel));
  Vec bias_g(n, 0.0);
  for (std::size_t i = 0; i < f_pred.size(); ++i) {
    for (std::size_t c = 0; c < n; ++c) bias_g[c] += (j_solver[i][c] - j_pred[i][c]) * f_pred[i];
  }
  const Vec delta_fd = solve_n(hess, bias_g);
  std::ostringstream terms;
  terms << std::setprecision(6);
  for (std::size_t c = 0; c < n; ++c) {
    terms << "variable " << c << ": reference " << p_ref[c] << ", Newton correction " << delta_n[c]
          << ", remainder " << delta_2[c] << ", FD bias " << delta_fd[c] << "; ";
  }
  INFO(terms.str());
  if (golden_width0) {
    // The independent reference lies at the optimum of the merit function (review of #197, P2):
    // otherwise a fault that moves the optimum (e.g. a wrong weighting in the generator) would
    // move p_pred and the run alike. Bounds, all from the reference before the run:
    // - variable 0 (golden section): the bracket width, plus the uncertainty of the minimum of a
    //   flat function sqrt(2 sigma_F / H_00), with the noise of F = 1/2 f^T f at most
    //   sigma_F = ||f|| sigma_f and sigma_f = kMeritPrecision ||f(p0)||inf (ADR 0030, point 9;
    //   ten times the measured noise);
    // - variable 1 (closed form for the reference's variable 0): an error dc of variable 0 moves
    //   its optimum along the valley dF/dtheta_1 = 0 by (H_10 / H_11) dc, plus its own flat term.
    REQUIRE(n == 2);
    double f_inf0 = 0.0;
    for (const double v : residuals_at(merit, merit.start()))
      f_inf0 = std::max(f_inf0, std::abs(v));
    double f_norm = 0.0;
    for (const double v : f_ref) f_norm += v * v;
    const double sigma_f_merit = std::sqrt(f_norm) * rtt::optim::kMeritPrecision * f_inf0;
    const double flat0 = std::sqrt(2.0 * sigma_f_merit / hess[0][0]);
    const double flat1 = std::sqrt(2.0 * sigma_f_merit / hess[1][1]);
    const double bound0 = *golden_width0 + flat0;
    const double bound1 = std::abs(hess[1][0] / hess[1][1]) * bound0 + flat1;
    INFO("reference at the optimum: |delta_n| <= " << bound0 << ", " << bound1 << " (golden width "
                                                   << *golden_width0 << ", flat " << flat0 << ", "
                                                   << flat1 << ")");
    CHECK(std::abs(delta_n[0]) <= bound0);
    CHECK(std::abs(delta_n[1]) <= bound1);
    // The Newton prediction is valid: its second step is small against the first.
    for (std::size_t c = 0; c < n; ++c) {
      CHECK(std::abs(delta_2[c]) <= 0.1 * std::abs(delta_n[c]) + 64.0 * kEps * std::abs(p_ref[c]));
    }
  }
  if (w_efl) {
    // The EFL residual of the least-squares optimum: e = f_0 / sqrt(w) at p_pred. The issue
    // requires EFL relative 1e-10 (1e-8 mm): the predicted residual must lie well below.
    const double e_pred = f_pred[0] / std::sqrt(*w_efl);
    INFO("predicted EFL residual e = " << e_pred << " mm");
    REQUIRE(std::abs(e_pred) <= 1e-9);
  }

  // --- The run ---
  OptimizeOptions options;
  options.ftol = 0.0;
  OptimResult r = optimize(s, lib, nullptr, options);
  INFO("status " << static_cast<int>(r.status) << ", iterations " << r.iterations
                 << ", evaluations " << r.evaluations);
  REQUIRE(r.status == LmStatus::ConvergedStep);
  REQUIRE(r.variables.size() == n);

  // Stopping bound from the last solve: D from the start Jacobian with the solver's step.
  const Vec p0 = merit.start();
  const Mat j0 = jacobian_at(merit, p0, steps(p0, h_rel));
  Vec d(n, 0.0);
  for (const Vec& row : j0) {
    for (std::size_t c = 0; c < n; ++c) d[c] += row[c] * row[c];
  }
  for (double& v : d) v = std::sqrt(v);
  // In the scaled variables z = D theta of the solver: h_z = -(A_z + mu I)^-1 g~_z at the stop and
  // e_z = H_z^-1 g~_z, so e_z = -H_z^-1 (A_z + mu I) h_z and ||e_z|| <= ||M||_F ||D h|| with
  // M = H_z^-1 (A_z + mu I), H_z = D^-1 H D^-1, A_z = D^-1 J^T J D^-1; per variable
  // |e_j| <= ||e_z|| / d_j.
  Mat h_z(n, Vec(n, 0.0));
  Mat am_z(n, Vec(n, 0.0));
  for (std::size_t rr = 0; rr < n; ++rr) {
    for (std::size_t c = 0; c < n; ++c) {
      h_z[rr][c] = hess[rr][c] / (d[rr] * d[c]);
      double a = 0.0;
      for (const Vec& row : j_pred) a += row[rr] * row[c];
      am_z[rr][c] = a / (d[rr] * d[c]);
    }
    am_z[rr][rr] += r.history.back().mu;
  }
  Mat m_z(n, Vec(n, 0.0));
  for (std::size_t c = 0; c < n; ++c) {
    Vec column(n, 0.0);
    for (std::size_t rr = 0; rr < n; ++rr) column[rr] = am_z[rr][c];
    const Vec x = solve_n(h_z, column);
    for (std::size_t rr = 0; rr < n; ++rr) m_z[rr][c] = x[rr];
  }
  const double e_z = frobenius(m_z) * r.history.back().step_norm;
  INFO("stop: ||M||_F = " << frobenius(m_z) << ", ||D h|| = " << r.history.back().step_norm
                          << ", ||e_z|| <= " << e_z);
  for (std::size_t c = 0; c < n; ++c) {
    const double tol = e_z / d[c] + std::abs(delta_fd[c]) + 2.0 * std::abs(delta_2[c]) +
                       64.0 * kEps * std::abs(p_pred[c]);
    INFO("variable " << c << ": result " << r.variables[c].end << ", predicted " << p_pred[c]
                     << ", tolerance " << tol);
    CHECK(std::abs(r.variables[c].end - p_pred[c]) <= tol);
  }
  // The generator table: every ray arrived, a real spot (content guard).
  REQUIRE(r.generators.size() == 1);
  CHECK(r.generators[0].rays_lost == 0);
  CHECK(r.generators[0].rms > 0.0);
  return r;
}

double efl_of(const System& s, const MaterialLibrary& lib) {
  const CompiledSystem cs = rtt::compile::compile(s, lib);
  return *rtt::paraxial::first_order(cs, *cs.find_path("main"), cs.reference_wavelength()).efl;
}

}  // namespace

TEST_CASE("M5 acceptance, case 1: singlet bent for minimal RMS spot at EFL 100 mm", "[optim][m5]") {
  // m5/singlet_solve: the EFL is held exactly by the parameter table (a "solve", decision of the
  // maintainer for #170; ADR 0030, addendum #170): rows D = n - n_a, K = 4 D / n, P = 1/100, C1
  // (variable), C2 = (P - D C1) / (D (K C1 - 1)) from the thick-lens power (Greivenkamp,
  // OPTI-201/202, p. 7-4: phi = phi1 + phi2 - phi1 phi2 t / n with phi1 = D c1, phi2 = -D c2),
  // R1 = 1 / C1 and R2 = 1 / C2 bound to the radii. Variables: C1 and the detector position;
  // merit: rms_spot (efl only as an observer, weight 0). The reference is the 1-D bending of
  // m5/singlet_optim (the same lens with free radii, independent of rtt-optim and of the table).
  // Tolerances are derived from the reference before the run (check_against_prediction):
  // - one Newton step on F = 1/2 f^T f from the reference predicts the least-squares optimum
  //   (golden-section width), a second step bounds the remainder;
  // - the solver's difference step h = eps_f^(1/3) max(|theta|, 1), eps_f = kMeritPrecision =
  //   1e-9 (ADR 0030, point 9), moves the fixed point by delta_FD = -H^-1 (J~ - J)^T f;
  // - the step test (ftol = 0, default xtol) leaves at most ||H_z^-1 (A_z + mu I)||_F ||D h|| / d_j
  //   in the scaled variables z = D theta: an a posteriori bound with a formula derived before
  //   the run (mu and ||D h|| come from the last solve, as in the step-test bounds of #166);
  // - the reference itself lies at the optimum of the merit function: |delta_n| within the
  //   golden-section width and the rounding of a flat minimum, |delta_2| << |delta_n|.
  // The default xtol: with a residual left (the spot) the gradient noise gives steps in the
  // flat bending direction that a step test of 1e-12 would never accept.
  const System s = load("m5/singlet_solve.rtt.json");
  const MaterialLibrary lib;
  // D in the file is n - n_a of the compiled environment at the reference wavelength, in the
  // shortest round-trip form: then the binding holds the EFL to rounding.
  {
    const CompiledSystem cs = rtt::compile::compile(s, lib);
    const double n_a = cs.media()[cs.environment_medium()].index[cs.reference_wavelength()].real();
    REQUIRE(std::get<double>(s.parameters[0].form) == 1.5168 - n_a);
  }
  const Bending ref = bending_reference(load("m5/singlet_optim.rtt.json"), lib);
  INFO("reference R1 " << ref.r1 << " R2 " << ref.r2 << " z " << ref.z_image);
  // The golden section over R1 ends with a bracket of 1e-6 mm, in C1 = 1/R1 a width of
  // 1e-6 / R1*^2.
  const OptimResult r = check_against_prediction(s, lib, Vec{1.0 / ref.r1, ref.z_image},
                                                 std::nullopt, 1e-6 / (ref.r1 * ref.r1));

  // The EFL of the result: relative 1e-10 (the issue), from the analysis and in the operand
  // table (the observer).
  const double efl = efl_of(r.system, lib);
  CHECK(std::abs(efl - 100.0) <= 1e-10 * 100.0);
  REQUIRE(r.operands.size() == 1);
  CHECK(r.operands[0].value == efl);
  CHECK(r.operands[0].contribution == 0.0);  // weight 0
  // The radii of the result follow the table (bound Params): R1 = 1 / C1.
  const auto table = rtt::model::evaluate_parameters(r.system);
  CHECK(table.at(5, 0) == 1.0 / r.variables[0].end);
  INFO("R1 result " << table.at(5, 0) << ", R2 result " << table.at(6, 0));
}

TEST_CASE("M5 acceptance, case 1 with the EFL as a penalty weight (known limit, #196)",
          "[.][optim][m5]") {
  // Hidden: m5/singlet_optim states the EFL as an operand with weight 1e4 and R1, R2 and the
  // detector free. Levenberg-Marquardt (and MINPACK lmdif, trf, BFGS and the MNT hybrid of
  // Sec. 3.4, measured for #170) crawls along the narrow, curved valley that the penalty makes
  // and stops far from the bending optimum: ADR 0030, addendum #170 "Bekannte Grenze:
  // Gleichheitsbedingungen als Strafgewicht". This test states the real goal (the reference
  // optimum R1*) and is enabled by the solution of #196 (exact equality constraints).
  const System s = load("m5/singlet_optim.rtt.json");
  const MaterialLibrary lib;
  const Bending ref = bending_reference(s, lib);
  const double w_efl =
      std::get<rtt::model::FirstOrderOperand>(s.optimization.operands[0]).common.weight;
  REQUIRE(w_efl == 1e4);
  const OptimResult r = check_against_prediction(s, lib, Vec{ref.r1, ref.r2, ref.z_image}, w_efl);
  CHECK(std::abs(efl_of(r.system, lib) - 100.0) <= 1e-10 * 100.0);
}

TEST_CASE("M5 acceptance, case 2: air gap D from Gullstrand's equation", "[optim][m5]") {
  const System s = load("m5/two_lens_gap.rtt.json");
  const MaterialLibrary lib;
  REQUIRE(s.optimization.operands.size() == 2);

  // --- Reference D* (independent of rtt-optim) -------------------------------------------------
  // Greivenkamp, OPTI-201/202, Sec. 7 "Gaussian Reduction", p. 7-4 (docs/quellen.md): two
  // components of powers phi1, phi2 (omega = n u, phi = -omega'/y) separated by t from the rear
  // principal plane P'1 of the first to the front principal plane P2 of the second in a medium of
  // index n2 have phi = phi1 + phi2 - phi1 phi2 tau, tau = t / n2. Here n2 = n_a (air) and phi is
  // first_order().power, EFL = 1 / phi. Solved for t: t = (phi1 + phi2 - phi) n_a / (phi1 phi2).
  // The powers and principal planes of L1 and L2 come from first_order of systems with only that
  // lens (children of the root: 0 stop, 1 L1, 2 L2, 3 image). The vertex distance D follows from
  // t = z_H2 - z_H'1 with the offsets of the principal planes from their vertices:
  // D = z_L2.S1 - z_L1.S2 = t + (z_H'1 - z_L1.S2) - (z_H2 - z_L2.S1).
  const CompiledSystem l1 = rtt::compile::compile(only(s, {0, 1, 3}), lib);
  const CompiledSystem l2 = rtt::compile::compile(only(s, {0, 2, 3}), lib);
  const std::uint16_t ref = l1.reference_wavelength();
  const rtt::paraxial::FirstOrder f1 = rtt::paraxial::first_order(l1, *l1.find_path("main"), ref);
  const rtt::paraxial::FirstOrder f2 = rtt::paraxial::first_order(l2, *l2.find_path("main"), ref);
  const double n_a = l1.media()[l1.environment_medium()].index[ref].real();
  const double phi = 1.0 / kTargetEfl;
  const double t = (f1.power + f2.power - phi) * n_a / (f1.power * f2.power);
  const double d_star =
      t + (*f1.rear_principal_z - z_of(l1, "L1.S2")) - (*f2.front_principal_z - z_of(l2, "L2.S1"));
  INFO("phi1 = " << f1.power << ", phi2 = " << f2.power << ", t = " << t << ", D* = " << d_star);
  // Content guards: the lenses converge, D* lies inside the bound (D >= 1) and away from the
  // start (D = 10), and the full system at D* has the target EFL (paraxial engine against the
  // two-component formula; both are y-nu traces, so they agree to rounding).
  REQUIRE(f1.power > 0.0);
  REQUIRE(d_star > 15.0);  // start D = 10, bound D >= 1
  System at_star = s;
  at_star.parameters[0].form = d_star;
  const CompiledSystem full = rtt::compile::compile(at_star, lib);
  const double efl_star = *rtt::paraxial::first_order(full, *full.find_path("main"), ref).efl;
  CHECK(std::abs(efl_star - kTargetEfl) <= 1e-12 * kTargetEfl);

  // --- Reference image distance (independent of rtt-optim) -------------------------------------
  // The marginal ray (field 0, py = 1, real aiming) of the system at D*: after the last surface it
  // is straight, so it crosses the axis at z_f = z - y dz/dy from any point (y, z) on it with
  // direction (dy, dz). Its distance behind L2.S2 is the image distance that zeroes the ray
  // operand.
  const std::uint16_t field = 0;
  rtt::trace::RayBatch rays =
      rtt::trace::make_rays(full, *full.find_path("main"), std::span(&field, 1), ref,
                            rtt::trace::SinglePupilPoint{0.0, 1.0});
  static_cast<void>(rtt::trace::SequentialTracer().trace(full, *full.find_path("main"), rays));
  REQUIRE(rays.status()[0] == rtt::trace::RayStatus::Alive);
  REQUIRE(rays.dir_y()[0] < 0.0);  // converging towards the axis
  const double z_focus = rays.pos_z()[0] - rays.pos_y()[0] * rays.dir_z()[0] / rays.dir_y()[0];
  const double image_star = z_focus - z_of(full, "L2.S2");
  INFO("image distance* = " << image_star);

  // --- Tolerances, derived before the run -----------------------------------------------------
  // ADR 0030, point 8, step test (MNT eq. (3.15b)), as in test_levenberg_marquardt.cpp: the run
  // ends when ||D h|| <= xtol (||D theta|| + xtol) and h is not applied. For a zero-residual
  // problem near the solution the damped step in z = D theta is h_z = -(A_z + mu I)^-1 A_z e_z,
  // A_z = J_z^T J_z, so ||e_z|| <= ||h_z|| (1 + mu / lambda_min); the test requires
  // mu <= lambda_min / 2 at the last solve and uses the factor 2. D: column norms of the Jacobian
  // in theta at the start (ADR 0030, point 6). In theta: |e_j| <= ||e_z|| / d_j; in p:
  // |dp/dtheta| times that (bounds.hpp, MINUIT transform for D >= 1; the image distance has no
  // bound, p = theta).
  const OptimizeOptions options = acceptance_options();
  const rtt::optim::MeritFunction merit(s, lib, nullptr);
  const auto& vars = merit.variables();
  REQUIRE(vars.size() == 2);
  REQUIRE(vars[0].pointer == "/parameters/0/value");  // D (table first, ADR 0030 point 5)
  const std::vector<double> p0 = merit.start();
  const std::vector<double> p_star = {d_star, image_star};
  const auto j0 = jacobian_theta(merit, p0);
  const auto js = jacobian_theta(merit, p_star);
  std::array<double, 2> d{};
  for (std::size_t c = 0; c < 2; ++c) {
    for (const auto& row : j0) d[c] += row[c] * row[c];
    d[c] = std::sqrt(d[c]);
    REQUIRE(d[c] > 0.0);
  }
  // A_z at the solution and its smallest eigenvalue (2 x 2, closed form).
  std::array<std::array<double, 2>, 2> a{};
  for (const auto& row : js) {
    for (std::size_t r = 0; r < 2; ++r) {
      for (std::size_t c = 0; c < 2; ++c) a[r][c] += row[r] / d[r] * (row[c] / d[c]);
    }
  }
  const double lambda_min =
      0.5 * (a[0][0] + a[1][1]) - std::hypot(0.5 * (a[0][0] - a[1][1]), a[0][1]);
  REQUIRE(lambda_min > 0.0);
  double scaled = 0.0;
  std::array<double, 2> theta_star{};
  for (std::size_t c = 0; c < 2; ++c) {
    theta_star[c] = rtt::optim::to_internal(p_star[c], vars[c].bounds);
    scaled += d[c] * theta_star[c] * d[c] * theta_star[c];
  }
  const double e_z = 2.0 * options.xtol * (std::sqrt(scaled) + options.xtol);
  std::array<double, 2> tol{};
  for (std::size_t c = 0; c < 2; ++c) {
    const double ht = 1e-6 * std::max(std::abs(theta_star[c]), 1.0);
    const double dp_dtheta = (rtt::optim::to_external(theta_star[c] + ht, vars[c].bounds) -
                              rtt::optim::to_external(theta_star[c] - ht, vars[c].bounds)) /
                             (2.0 * ht);
    tol[c] = std::abs(dp_dtheta) * e_z / d[c];
  }
  // The EFL residual depends on D only: its error is J_theta[0][0] times the error of theta_D,
  // at most e_z / d_0, plus the rounding of the EFL itself.
  const double tol_efl = std::abs(js[0][0]) * e_z / d[0] + 64.0 * kEps * kTargetEfl;
  INFO("lambda_min = " << lambda_min << ", tol D = " << tol[0] << " mm, tol image = " << tol[1]
                       << " mm, tol EFL = " << tol_efl << " mm");
  // The issue requires EFL relative 1e-10: the derived bound must lie below it.
  REQUIRE(tol_efl <= 1e-10 * kTargetEfl);

  // --- The run ----------------------------------------------------------------------------------
  const OptimResult r = optimize(s, lib, nullptr, options);
  INFO("status " << static_cast<int>(r.status) << ", iterations " << r.iterations);
  REQUIRE((r.status == LmStatus::ConvergedStep || r.status == LmStatus::ConvergedGradient));
  REQUIRE(!r.history.empty());
  REQUIRE(r.history.back().mu <= lambda_min / 2.0);
  REQUIRE(r.variables.size() == 2);
  CHECK(std::abs(r.variables[0].end - d_star) <= tol[0]);
  CHECK(std::abs(r.variables[1].end - image_star) <= tol[1]);
  CHECK(r.diagnostics.empty());

  // The EFL of the result, from the analysis itself.
  const CompiledSystem result = rtt::compile::compile(r.system, lib);
  const double efl = *rtt::paraxial::first_order(result, *result.find_path("main"), ref).efl;
  CHECK(std::abs(efl - kTargetEfl) <= tol_efl);
  CHECK(std::abs(efl - kTargetEfl) <= 1e-10 * kTargetEfl);

  // Guards of the issue: L2.S1 lies D behind L1.S2, and the detector moved with L2 (its distance
  // behind L2.S2 is the image distance variable). Sums of a few doubles: within rounding of the
  // global z (a few ulp).
  const double d_end = std::get<double>(r.system.parameters[0].form);
  CHECK(d_end == r.variables[0].end);
  const double z_s1 = z_of(result, "L2.S1");
  CHECK(std::abs((z_s1 - z_of(result, "L1.S2")) - d_end) <= 4.0 * kEps * z_s1);
  const double z_img = z_of(result, "IMG");
  CHECK(std::abs((z_img - z_of(result, "L2.S2")) - r.variables[1].end) <= 4.0 * kEps * z_img);
  CHECK(z_s1 != z_of(rtt::compile::compile(s, lib), "L2.S1"));  // content guard: L2 moved
}

TEST_CASE("M5 acceptance, case 2: bitwise for 1, 4 and all threads, with RunControl, repeated",
          "[optim][m5]") {
  const System s = load("m5/two_lens_gap.rtt.json");
  const MaterialLibrary lib;
  const OptimizeOptions options = acceptance_options();
  const auto run = [&](std::size_t threads, bool control) {
    const tbb::global_control limit(tbb::global_control::max_allowed_parallelism, threads);
    rtt::trace::RunControl c;
    if (control) c.cancel = rtt::trace::CancelToken{};  // active, never requested
    return optimize(s, lib, nullptr, options, c);
  };
  const OptimResult one = run(1, false);
  REQUIRE(one.patch != "[]");  // content guard
  for (const std::size_t threads : {std::size_t{1}, std::size_t{4},
                                    static_cast<std::size_t>(tbb::info::default_concurrency())}) {
    for (const bool control : {false, true}) {
      INFO("threads " << threads << ", control " << control);
      const OptimResult other = run(threads, control);
      CHECK(other.patch == one.patch);
      CHECK(other.evaluations == one.evaluations);
      REQUIRE(other.history.size() == one.history.size());
      for (std::size_t k = 0; k < one.history.size(); ++k) {
        CHECK(std::bit_cast<std::uint64_t>(other.history[k].phi) ==
              std::bit_cast<std::uint64_t>(one.history[k].phi));
        CHECK(std::bit_cast<std::uint64_t>(other.history[k].mu) ==
              std::bit_cast<std::uint64_t>(one.history[k].mu));
      }
    }
  }
}
