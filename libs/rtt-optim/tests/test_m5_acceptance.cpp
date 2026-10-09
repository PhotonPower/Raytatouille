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
#include <limits>
#include <optional>
#include <span>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
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

}  // namespace

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
