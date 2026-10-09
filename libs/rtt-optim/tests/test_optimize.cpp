// Operands, merit evaluation and optimize() (ADR 0030, points 1-3, 5, 10-12; #167, tests M1-M8).
// Tolerances are derived beforehand (solution_tolerance below), not measured.

#include <tbb/global_control.h>
#include <tbb/info.h>

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

#include "rtt/analysis/opd.hpp"
#include "rtt/analysis/spot.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/compile/errors.hpp"
#include "rtt/io/edit.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/model/optimization.hpp"
#include "rtt/model/parameters.hpp"
#include "rtt/optim/merit.hpp"
#include "rtt/optim/optimize.hpp"
#include "rtt/paraxial/paraxial.hpp"
#include "rtt/paraxial/prescription.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/run_control.hpp"
#include "rtt/trace/sequential.hpp"
#include "rtt/trace/sources.hpp"

using rtt::material::MaterialLibrary;
using rtt::model::FirstOrderOperand;
using rtt::model::FirstOrderQuantity;
using rtt::model::OperandCommon;
using rtt::model::Param;
using rtt::model::System;
using rtt::optim::LmStatus;
using rtt::optim::MeritEvaluation;
using rtt::optim::MeritFunction;
using rtt::optim::OptimError;
using rtt::optim::optimize;
using rtt::optim::OptimizeOptions;
using rtt::optim::OptimResult;

namespace {

constexpr double kN = 1.5168;  // CONST index of m1/singlet_const

System load(const std::string& relative) {
  return rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/" + relative);
}

rtt::model::Element& element(System& s, std::size_t child) {
  return std::get<rtt::model::Element>(s.root.children[child].value);
}

Param& radius(System& s, std::size_t child, std::size_t surface) {
  return std::get<rtt::model::Conic>(element(s, child).surfaces[surface].shape.base).radius;
}

OperandCommon common(double target,
                     double weight = 1.0,
                     std::optional<std::string> configuration = std::nullopt) {
  return {target, weight, std::move(configuration)};
}

FirstOrderOperand efl(double target, std::optional<std::string> configuration = std::nullopt) {
  FirstOrderOperand op;
  op.common = common(target, 1.0, std::move(configuration));
  op.quantity = FirstOrderQuantity::Efl;
  op.path = "main";
  return op;
}

/// Makes the first unbound Param whose pointer ends with `suffix` variable; returns its pointer.
std::string make_variable(System& s, std::string_view suffix) {
  std::string found;
  rtt::model::for_each_param(s, [&](std::string_view pointer, Param& p) {
    if (found.empty() && pointer.ends_with(suffix) && !p.is_bound()) {
      p.variable = true;
      found = pointer;
    }
  });
  REQUIRE(!found.empty());
  return found;
}

/// m1/singlet_const: stop, plano-convex L1 (R1 = 51.68 on L1.S1 at z = 5, L1.S2 plane at z = 9),
/// detector IMG at z = 106.363; object at infinity, EPD 20, fields y = 0, 3.5, 5 mm, three
/// wavelengths (reference 0.5876 um).
System singlet() {
  return load("m1/singlet_const.rtt.json");
}

// A priori bound on the error of the solution when the run stops with the step test (ADR 0030,
// point 8; MNT eq. (3.15b)), as in test_levenberg_marquardt.cpp: ||D h|| <= xtol (||D theta|| +
// xtol), and with mu <= lambda_min / 2 at the last solve ||e_z|| <= 1.5 ||h_z||; the bound uses
// 2. Here every Jacobian is diagonal with entries j_i, d_i = |j_i| (no bounds, theta = p), so
// A_z = I and lambda_min = 1. In p: |e_i| <= ||e_z|| / d_i.
double solution_tolerance(std::span<const double> d,
                          std::span<const double> p_star,
                          double mu_last,
                          double xtol) {
  INFO("mu at the last solve = " << mu_last);
  REQUIRE(mu_last <= 0.5);
  double scaled = 0.0;
  double d_min = d[0];
  for (std::size_t i = 0; i < d.size(); ++i) {
    scaled += d[i] * p_star[i] * d[i] * p_star[i];
    d_min = std::min(d_min, d[i]);
  }
  return 2.0 * xtol * (std::sqrt(scaled) + xtol) / d_min;
}

}  // namespace

TEST_CASE("merit M9: kMeritPrecision is a power of ten in [eps_M, 1e-6]", "[optim][merit]") {
  // The value comes from the measurement in test_merit_noise.cpp ([.noise]); above 1e-6 the
  // finding goes to the coordinator first (plan #167).
  constexpr double p = rtt::optim::kMeritPrecision;
  CHECK(p >= std::numeric_limits<double>::epsilon());
  CHECK(p <= 1e-6);
  const double exponent = std::log10(p);
  CHECK(std::pow(10.0, std::round(exponent)) == p);
}

TEST_CASE("merit M1: every operand is exactly its analysis", "[optim][merit]") {
  System s = singlet();
  s.parameters = {{"K", 2.5}};
  auto& ops = s.optimization.operands;
  ops.emplace_back(efl(0.0));  // 0
  FirstOrderOperand bfl = efl(0.0);
  bfl.quantity = FirstOrderQuantity::Bfl;
  bfl.wavelength = 0;
  ops.emplace_back(bfl);  // 1
  FirstOrderOperand fno = efl(0.0);
  fno.quantity = FirstOrderQuantity::ImageFNumber;
  ops.emplace_back(fno);  // 2
  rtt::model::RayOperand ray_img;
  ray_img.common = common(0.0);
  ray_img.path = "main";
  ray_img.surface = rtt::model::SurfaceId{"IMG"};
  ray_img.field = 2;
  ray_img.px = 0.25;
  ray_img.py = 0.75;
  ops.emplace_back(ray_img);  // 3: ray_y at the image
  rtt::model::RayOperand ray_stop = ray_img;
  ray_stop.surface = rtt::model::SurfaceId{"STO"};
  ray_stop.coordinate = rtt::model::RayCoordinate::X;
  ops.emplace_back(ray_stop);  // 4: ray_x at the stop
  rtt::model::SpotRmsOperand spot;
  spot.common = common(0.0);
  spot.path = "main";
  spot.field = 1;
  spot.rings = 4;
  ops.emplace_back(spot);  // 5: centroid, reference wavelength
  rtt::model::SpotRmsOperand poly = spot;
  poly.field = 2;
  poly.polychromatic = true;
  poly.reference = rtt::model::SpotReference::Chief;
  ops.emplace_back(poly);  // 6: chief, all wavelengths
  rtt::model::OpdRmsOperand opd;
  opd.common = common(0.0);
  opd.path = "main";
  opd.field = 1;
  opd.wavelength = 2;
  opd.grid = 17;
  ops.emplace_back(opd);  // 7
  rtt::model::ParamValueOperand param;
  param.common = common(0.0);
  param.parameter = "K";
  ops.emplace_back(param);  // 8

  const MaterialLibrary lib;
  const MeritFunction merit(s, lib, nullptr);
  REQUIRE(merit.variables().empty());
  REQUIRE(merit.size() == 9);
  const MeritEvaluation e = merit.evaluate({});
  REQUIRE(e.valid());

  const rtt::compile::CompiledSystem cs = rtt::compile::compile(s, lib);
  const rtt::compile::PathId main = *cs.find_path("main");
  const std::uint16_t ref = cs.reference_wavelength();
  const rtt::paraxial::FirstOrder fo = rtt::paraxial::first_order(cs, main, ref);
  CHECK(e.values[0] == *fo.efl);
  CHECK(e.values[1] == *rtt::paraxial::first_order(cs, main, 0).bfl);
  CHECK(e.values[2] == *rtt::paraxial::prescription(cs, main, ref).paraxial_working_f_number);
  // The ray at the last event: the state the plain trace leaves in the batch, in local
  // coordinates of IMG (the detector is not tilted, so local = global - vertex).
  const std::uint16_t field = 2;
  rtt::trace::RayBatch rays = rtt::trace::make_rays(cs, main, std::span(&field, 1), ref,
                                                    rtt::trace::SinglePupilPoint{0.25, 0.75});
  static_cast<void>(rtt::trace::SequentialTracer().trace(cs, main, rays));
  REQUIRE(rays.status()[0] == rtt::trace::RayStatus::Alive);
  const std::uint32_t img = *cs.find_surface(rtt::model::SurfaceId{"IMG"});
  const rtt::math::Vec3 local =
      cs.surfaces()[img].to_local.apply_point({rays.pos_x()[0], rays.pos_y()[0], rays.pos_z()[0]});
  CHECK(e.values[3] == local.y());
  // At the stop real aiming hits the target px R_s = 0.25 * 10 mm within kAimTolerance (the
  // stop is the first surface, so the entrance pupil is the stop, R_s = EPD / 2).
  CHECK(std::abs(e.values[4] - 2.5) <= rtt::trace::kAimTolerance);
  rtt::analysis::SpotOptions so;
  so.sampling = rtt::trace::HexapolarPupil{4};
  CHECK(e.values[5] == rtt::analysis::spot(cs, main, 1, ref, so).stats.rms_centroid);
  CHECK(e.values[6] == rtt::analysis::spot(cs, main, 2, std::nullopt, so).stats.rms_chief);
  rtt::analysis::OpdOptions oo;
  oo.grid = 17;
  CHECK(e.values[7] == rtt::analysis::opd_map(cs, main, 1, 2, oo).rms);
  CHECK(e.values[8] == 2.5);
  // Target 0, weight 1: the residual is the value.
  for (std::size_t i = 0; i < e.values.size(); ++i) CHECK(e.residuals[i] == e.values[i]);
}

TEST_CASE("merit M1: residual sqrt(w) (v - t)", "[optim][merit]") {
  System s = singlet();
  FirstOrderOperand op = efl(80.0);
  op.common.weight = 4.0;
  s.optimization.operands.emplace_back(op);
  const MaterialLibrary lib;
  const MeritEvaluation e = MeritFunction(s, lib, nullptr).evaluate({});
  CHECK(e.residuals[0] == 2.0 * (e.values[0] - 80.0));
}

TEST_CASE("merit M2: an evaluation equals the system changed by hand", "[optim][merit]") {
  System s = singlet();
  s.configurations = {{"near"}, {"far"}};
  s.parameters = {{"Z", 5.0, true}, {"R", std::vector<double>{51.68, 60.0}, true}};
  element(s, 1).pose.position[2] = Param::bound("Z");
  radius(s, 1, 0) = Param::bound("R");
  element(s, 2).pose.position[2].variable = true;  // detector z
  s.aperture.value.variable = true;                // EPD
  rtt::model::RayOperand ray;
  ray.common = common(0.0, 1.0, "far");
  ray.path = "main";
  ray.surface = rtt::model::SurfaceId{"IMG"};
  ray.field = 2;
  ray.py = 1.0;
  s.optimization.operands = {efl(0.0, "near"), efl(0.0, "far"), ray,
                             rtt::model::ParamValueOperand{common(0.0, 1.0, "far"), "R"},
                             rtt::model::ParamValueOperand{common(0.0), "Z"}};
  const MaterialLibrary lib;
  const MeritFunction merit(s, lib, nullptr);
  // ADR 0030, point 5: table rows first (R per configuration), then the model Params.
  const auto& vars = merit.variables();
  REQUIRE(vars.size() == 5);
  CHECK(vars[0].pointer == "/parameters/0/value");
  CHECK(vars[1].pointer == "/parameters/1/values/0");
  CHECK(vars[2].pointer == "/parameters/1/values/1");
  CHECK(vars[3].row.empty());
  CHECK(vars[4].row.empty());
  const std::vector<double> p = {5.5, 50.0, 62.5, 21.0, 104.0};
  // vars[3] and vars[4] in edit-form order: the aperture comes before the tree.
  REQUIRE(vars[3].pointer == "/aperture/value/value");
  const MeritEvaluation e = merit.evaluate(p);
  REQUIRE(e.valid());

  System hand = s;
  hand.parameters[0].form = 5.5;
  hand.parameters[1].form = std::vector<double>{50.0, 62.5};
  hand.aperture.value.value = 21.0;
  element(hand, 2).pose.position[2].value = 104.0;
  const rtt::compile::CompiledSystem near = rtt::compile::compile(hand, lib, 0);
  const rtt::compile::CompiledSystem far = rtt::compile::compile(hand, lib, 1);
  CHECK(e.values[0] == *rtt::paraxial::first_order(near, *near.find_path("main"), 1).efl);
  CHECK(e.values[1] == *rtt::paraxial::first_order(far, *far.find_path("main"), 1).efl);
  CHECK(e.values[0] != e.values[1]);  // content guard: the configurations differ
  const std::uint16_t field = 2;
  const rtt::compile::PathId main = *far.find_path("main");
  rtt::trace::RayBatch rays =
      rtt::trace::make_rays(far, main, std::span(&field, 1), 1, rtt::trace::SinglePupilPoint{0, 1});
  static_cast<void>(rtt::trace::SequentialTracer().trace(far, main, rays));
  const std::uint32_t img = *far.find_surface(rtt::model::SurfaceId{"IMG"});
  CHECK(e.values[2] ==
        far.surfaces()[img]
            .to_local.apply_point({rays.pos_x()[0], rays.pos_y()[0], rays.pos_z()[0]})
            .y());
  CHECK(e.values[3] == 62.5);
  CHECK(e.values[4] == 5.5);
}

TEST_CASE("optimize M3: EFL of a plano-convex lens, analytic R1", "[optim][optimize]") {
  // Paraxial y-nu trace (Greivenkamp, OPTI-201/202, p. 9-2): a ray parallel to the axis at
  // height y gets n u' = -y (n - 1) / R1 at L1.S1; L1.S2 is plane (phi = 0), so after it
  // u'' = -y (n - 1) / R1 for every thickness, and Phi = -u'' / y = (n - 1) / R1 (the focal
  // length of a single refracting surface, docs/quellen.md). EFL = R1 / (n - 1), linear in R1:
  // target 80 mm gives R1* = 80 (n - 1) = 41.344 mm, and d = |dEFL/dR1| = 1 / (n - 1).
  System s = singlet();
  radius(s, 1, 0) = Param(60.0);
  radius(s, 1, 0).variable = true;
  s.optimization.operands = {efl(80.0)};
  OptimizeOptions options;
  options.ftol = 0.0;  // only the step test (or exact zero) ends the run: the bound below holds
  const MaterialLibrary lib;
  OptimResult r = optimize(s, lib, nullptr, options);
  REQUIRE((r.status == LmStatus::ConvergedStep || r.status == LmStatus::ConvergedGradient));
  const double r_star = 80.0 * (kN - 1.0);
  const std::vector<double> d = {1.0 / (kN - 1.0)};
  const double tol = solution_tolerance(
      d, std::vector<double>{r_star}, r.history.empty() ? 0.0 : r.history.back().mu, options.xtol);
  INFO("tolerance = " << tol);
  REQUIRE(r.variables.size() == 1);
  CHECK(std::abs(r.variables[0].end - r_star) <= tol);
  CHECK(radius(r.system, 1, 0).value == r.variables[0].end);
  CHECK(r.variables[0].changed);
  CHECK(r.diagnostics.empty());
  // The final state is the solver's state: residual bitwise from the same evaluation.
  REQUIRE(r.operands.size() == 1);
  CHECK(r.operands[0].contribution == (r.operands[0].value == 80.0 ? 0.0 : 100.0));
}

TEST_CASE("optimize M4: two configurations, two independent unknowns", "[optim][optimize]") {
  System s = singlet();
  s.configurations = {{"long"}, {"short"}};
  s.parameters = {{"R", std::vector<double>{60.0, 60.0}, true}};
  radius(s, 1, 0) = Param::bound("R");
  s.optimization.operands = {efl(100.0, "long"), efl(80.0, "short")};
  OptimizeOptions options;
  options.ftol = 0.0;
  const MaterialLibrary lib;
  OptimResult r = optimize(s, lib, nullptr, options);
  REQUIRE((r.status == LmStatus::ConvergedStep || r.status == LmStatus::ConvergedGradient));
  const std::vector<double> p_star = {100.0 * (kN - 1.0), 80.0 * (kN - 1.0)};
  const std::vector<double> d = {1.0 / (kN - 1.0), 1.0 / (kN - 1.0)};
  REQUIRE(!r.history.empty());
  const double tol = solution_tolerance(d, p_star, r.history.back().mu, options.xtol);
  INFO("tolerance = " << tol);
  REQUIRE(r.variables.size() == 2);
  CHECK(r.variables[0].pointer == "/parameters/0/values/0");
  CHECK(r.variables[1].pointer == "/parameters/0/values/1");
  CHECK(std::abs(r.variables[0].end - p_star[0]) <= tol);
  CHECK(std::abs(r.variables[1].end - p_star[1]) <= tol);
  // The binding stays; the table holds the result.
  CHECK(radius(r.system, 1, 0).param == std::optional<std::string>("R"));
  CHECK(std::get<std::vector<double>>(r.system.parameters[0].form) ==
        std::vector<double>{r.variables[0].end, r.variables[1].end});
}

TEST_CASE("optimize M5: the result patch", "[optim][optimize]") {
  System s = singlet();
  s.configurations = {{"long"}, {"short"}};
  s.parameters = {{"R", std::vector<double>{60.0, 60.0}, true}};
  radius(s, 1, 0) = Param::bound("R");
  element(s, 2).pose.position[2].variable = true;
  rtt::model::RayOperand ray;
  ray.common = common(0.0);
  ray.path = "main";
  ray.surface = rtt::model::SurfaceId{"IMG"};
  ray.py = 1.0;
  s.optimization.operands = {efl(100.0, "long"), efl(80.0, "short"), ray};
  const MaterialLibrary lib;
  const OptimResult r = optimize(s, lib);
  REQUIRE(r.patch != "[]");
  const rtt::io::PatchResult applied = rtt::io::apply_patch(s, r.patch);
  CHECK(rtt::io::to_json(applied.system) == rtt::io::to_json(r.system));
  const rtt::io::PatchResult back =
      rtt::io::apply_patch(applied.system, applied.inverse, rtt::io::EditCheck::StructureOnly);
  CHECK(rtt::io::to_json(back.system) == rtt::io::to_json(s));

  SECTION("no accepted step: empty patch, bitwise the input") {
    OptimizeOptions none;
    none.max_iterations = 0;
    const OptimResult z = optimize(s, lib, nullptr, none);
    CHECK(z.status == LmStatus::MaxIterations);
    CHECK(z.patch == "[]");
    CHECK(z.system == s);
    CHECK(rtt::io::to_json(z.system) == rtt::io::to_json(s));
  }
}

TEST_CASE("optimize M6: bitwise for 1, 4 and all threads, with and without RunControl",
          "[optim][optimize]") {
  System s = singlet();
  s.configurations = {{"long"}, {"short"}};
  s.parameters = {{"R", std::vector<double>{60.0, 60.0}, true}};
  radius(s, 1, 0) = Param::bound("R");
  element(s, 2).pose.position[2].variable = true;
  rtt::model::SpotRmsOperand spot;
  spot.common = common(0.0, 1.0, "long");
  spot.path = "main";
  s.optimization.operands = {efl(100.0, "long"), efl(80.0, "short"), spot};
  const MaterialLibrary lib;
  const auto run = [&](std::size_t threads, bool control) {
    const tbb::global_control limit(tbb::global_control::max_allowed_parallelism, threads);
    rtt::trace::RunControl c;
    if (control) c.cancel = rtt::trace::CancelToken{};  // active, never requested
    return optimize(s, lib, nullptr, {}, c);
  };
  const OptimResult one = run(1, false);
  for (const std::size_t threads : {std::size_t{1}, std::size_t{4},
                                    static_cast<std::size_t>(tbb::info::default_concurrency())}) {
    for (const bool control : {false, true}) {
      INFO("threads " << threads << ", control " << control);
      const OptimResult other = run(threads, control);
      CHECK(other.patch == one.patch);
      CHECK(other.evaluations == one.evaluations);
      REQUIRE(other.history.size() == one.history.size());
      for (std::size_t k = 0; k < one.history.size(); ++k) {
        CHECK(other.history[k].phi == one.history[k].phi);
      }
    }
  }
  CHECK(one.patch != "[]");  // content guard
}

TEST_CASE("optimize M7: input errors at the start", "[optim][optimize]") {
  const MaterialLibrary lib;
  SECTION("no variable") {
    System s = singlet();
    s.optimization.operands = {efl(80.0)};
    try {
      static_cast<void>(optimize(s, lib));
      FAIL("no OptimError");
    } catch (const OptimError& e) {
      REQUIRE(e.diagnostics().size() == 1);
      CHECK(e.diagnostics()[0].code == "optim.no_variables");
    }
  }
  SECTION("no operand and no generator") {
    System s = singlet();
    radius(s, 1, 0).variable = true;
    try {
      static_cast<void>(optimize(s, lib));
      FAIL("no OptimError");
    } catch (const OptimError& e) {
      REQUIRE(e.diagnostics().size() == 1);
      CHECK(e.diagnostics()[0].code == "optim.no_operands");
    }
  }
  SECTION("neither variables nor operands: both codes") {
    try {
      static_cast<void>(optimize(singlet(), lib));
      FAIL("no OptimError");
    } catch (const OptimError& e) {
      REQUIRE(e.diagnostics().size() == 2);
      CHECK(e.diagnostics()[0].code == "optim.no_variables");
      CHECK(e.diagnostics()[1].code == "optim.no_operands");
    }
  }
  SECTION("magnification with the object at infinity") {
    System s = singlet();
    radius(s, 1, 0).variable = true;
    FirstOrderOperand m = efl(-1.0);
    m.quantity = FirstOrderQuantity::Magnification;
    s.optimization.operands = {efl(80.0), m};
    try {
      static_cast<void>(optimize(s, lib));
      FAIL("no OptimError");
    } catch (const OptimError& e) {
      REQUIRE(e.diagnostics().size() == 1);
      CHECK(e.diagnostics()[0].code == "merit.operand_unsupported");
      CHECK(e.diagnostics()[0].location == "/optimization/operands/1");
    }
  }
  SECTION("a generator before #168") {
    System s = singlet();
    radius(s, 1, 0).variable = true;
    rtt::model::SpotGenerator g;
    g.path = "main";
    s.optimization.generators = {g};
    try {
      static_cast<void>(optimize(s, lib));
      FAIL("no OptimError");
    } catch (const OptimError& e) {
      REQUIRE(e.diagnostics().size() == 1);
      CHECK(e.diagnostics()[0].code == "merit.operand_unsupported");
      CHECK(e.diagnostics()[0].location == "/optimization/generators/0");
    }
  }
  SECTION("an operand on a crystal path (ordinary and extraordinary events)") {
    System s = load("m4/calcite_walkoff.rtt.json");
    make_variable(s, "/position/2");
    FirstOrderOperand op = efl(80.0);
    op.path = "o";
    s.optimization.operands = {op};
    try {
      static_cast<void>(optimize(s, lib));
      FAIL("no OptimError");
    } catch (const OptimError& e) {
      REQUIRE(e.diagnostics().size() == 1);
      CHECK(e.diagnostics()[0].code == "merit.operand_unsupported");
      CHECK(e.diagnostics()[0].location == "/optimization/operands/0");
    }
  }
  SECTION("an undefined value at the start names the operand") {
    // m2/telecentric_4f is afocal: first_order accepts the path, but there is no EFL.
    System s = load("m2/telecentric_4f.rtt.json");
    make_variable(s, "/position/2");
    s.optimization.operands = {efl(80.0)};
    try {
      static_cast<void>(optimize(s, lib));
      FAIL("no exception");
    } catch (const std::invalid_argument& e) {
      CHECK(std::string(e.what()).find("/optimization/operands/0") != std::string::npos);
    }
  }
}

TEST_CASE("optimize M7: an invalid trial is a rejected step with a warning", "[optim][optimize]") {
  // EFL = R1 / (n - 1) with target 15 mm wants R1* = 15 (n - 1) = 7.752 mm, but the rim ray
  // (py = 1, height 10 mm at the stop = L1.S1 aperture side) cannot meet a sphere of radius
  // below 10 mm. The first trial, close to the undamped step R1 = 60 -> 7.75 (tau = 1e-3, the
  // problem is linear in R1), therefore has an undefined ray operand: an invalid evaluation,
  // rejected (ADR 0030, point 10). Content guard: the ray operand is valid at the start.
  System s = singlet();
  radius(s, 1, 0) = Param(60.0);
  radius(s, 1, 0).variable = true;
  rtt::model::RayOperand ray;
  ray.common = common(0.0, 1e-6);
  ray.path = "main";
  ray.surface = rtt::model::SurfaceId{"IMG"};
  ray.py = 1.0;
  s.optimization.operands = {efl(15.0), ray};
  const MaterialLibrary lib;
  REQUIRE(MeritFunction(s, lib, nullptr).evaluate(std::vector<double>{60.0}).valid());
  const OptimResult r = optimize(s, lib);
  CHECK(r.failed_evaluations > 0);
  std::size_t warnings = 0;
  for (const auto& d : r.diagnostics) warnings += d.code == "optim.evaluation_failed" ? 1 : 0;
  CHECK(warnings == 1);
}

TEST_CASE("optimize M7: a Jacobian column without a valid evaluation fails the run",
          "[optim][optimize]") {
  // Row E = X * 1e308 overflows for X > 1.7976931348623157 (DBL_MAX = 1.7976931348623157e308),
  // and validate reports a non-finite row (parameters.not_finite), so compile throws. Start
  // X = 1.79769 is valid; the difference step h = eps_f^(1/3) max(|X|, 1) >= 1e-5 (eps_f <=
  // 1e-6, ADR 0030 point 9) puts X + h above the limit: the column of X has no valid pair, the
  // run ends with status Failed and the input (ADR 0030, point 10).
  System s = singlet();
  s.parameters = {{"X", 1.79769, true}, {"E", rtt::model::ParameterExpression{"X * 1e308"}}};
  s.optimization.operands = {efl(80.0)};
  const MaterialLibrary lib;
  const OptimResult r = optimize(s, lib);
  CHECK(r.status == LmStatus::Failed);
  CHECK(r.patch == "[]");
  CHECK(r.system == s);
  CHECK(r.failed_evaluations >= 1);
  bool reported = false;
  for (const auto& d : r.diagnostics) {
    reported =
        reported || (d.code == "optim.jacobian_failed" && d.location == "/parameters/0/value");
  }
  CHECK(reported);
}

TEST_CASE("optimize M7: cancellation keeps the last accepted state", "[optim][optimize]") {
  System s = singlet();
  radius(s, 1, 0) = Param(60.0);
  radius(s, 1, 0).variable = true;
  s.optimization.operands = {efl(80.0)};
  const MaterialLibrary lib;
  rtt::trace::CancelToken token;
  token.request_cancel();  // before the start: the run ends with the input
  rtt::trace::RunControl c;
  c.cancel = token;
  const OptimResult r = optimize(s, lib, nullptr, {}, c);
  CHECK(r.status == LmStatus::Cancelled);
  CHECK(r.patch == "[]");
  CHECK(r.system == s);
}

TEST_CASE("optimize M8: a bound that excludes the optimum", "[optim][optimize]") {
  System s = singlet();
  radius(s, 1, 0) = Param(60.0);
  radius(s, 1, 0).variable = true;
  radius(s, 1, 0).min = 45.0;  // R1* = 41.344 lies below
  s.optimization.operands = {efl(80.0)};
  const MaterialLibrary lib;
  const OptimResult r = optimize(s, lib);
  REQUIRE(r.variables.size() == 1);
  CHECK(r.variables[0].at_bound);
  CHECK(std::abs(r.variables[0].end - 45.0) <= 1e-6);
  bool warned = false;
  for (const auto& d : r.diagnostics) {
    warned =
        warned || (d.code == "optim.parameter_at_bound" && d.location == r.variables[0].pointer);
  }
  CHECK(warned);
}
