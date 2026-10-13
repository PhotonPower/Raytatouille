// The merit function and optimize with the overlay of a tolerance trial (ADR 0030, addendum #256;
// ADR 0036, points 2 and 3): the perturbation is compiled into every evaluation, aiming and
// reference spheres come from the nominal basis. With an empty perturbation and the basis of the
// same system everything but the paraxial operands is bitwise the run without overlay.

#include <oneapi/tbb/info.h>
#include <oneapi/tbb/task_arena.h>

#include <algorithm>
#include <bit>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <optional>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <variant>
#include <vector>

#include "rtt/coating/catalog.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/compile/perturbation.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
#include "rtt/optim/merit.hpp"
#include "rtt/optim/optimize.hpp"
#include "rtt/optim/overlay.hpp"

namespace fs = std::filesystem;
using rtt::model::Element;
using rtt::model::ElementKind;
using rtt::model::FieldType;
using rtt::model::Param;
using rtt::model::Pose;
using rtt::model::Surface;
using rtt::model::SurfaceId;
using rtt::model::System;
using rtt::optim::MeritEvaluation;
using rtt::optim::MeritFunction;
using rtt::optim::TrialOverlay;

namespace {

bool same(double a, double b) {
  return std::bit_cast<std::uint64_t>(a) == std::bit_cast<std::uint64_t>(b);
}

const rtt::material::MaterialLibrary& library() {
  static const rtt::material::MaterialLibrary materials;
  return materials;
}

Surface plane(const std::string& id, double z_mm = 0.0) {
  Surface s;
  s.id = SurfaceId(id);
  s.pose = Pose::along_z(z_mm);
  return s;
}

/// Ideal lens f = 100 mm in vacuum: stop r = 10 at z = -5, lens at z = 0, detector at z = 100 with
/// a variable position x; object at infinity, EPD 10, field 0.
System ideal_lens() {
  System s;
  s.name = "ideal lens";
  s.environment.medium = "VACUUM";
  s.wavelengths = {{0.5876, 1.0, true}};
  s.aperture = {rtt::model::SystemApertureType::EntrancePupilDiameter, Param(10.0)};
  s.fields = {FieldType::AngleDeg, {{0.0, 0.0, 1.0}}};
  s.root.name = "root";
  Surface stop = plane("STO");
  stop.aperture = rtt::model::CircularAperture{10.0, 0.0};
  Surface lens = plane("IL");
  lens.aperture = rtt::model::CircularAperture{30.0, 0.0};
  lens.interaction = rtt::model::IdealLens{Param(100.0), std::nullopt};
  Pose detector = Pose::along_z(100.0);
  detector.position[0].variable = true;
  s.root.children = {
      {Element{"stop", ElementKind::Stop, Pose::along_z(-5.0), std::nullopt, {stop}}},
      {Element{"lens", ElementKind::ThinElement, Pose::along_z(0.0), std::nullopt, {lens}}},
      {Element{"D", ElementKind::Detector, detector, std::nullopt, {plane("IMG")}}}};
  s.paths = {{"main", true, {}}};
  rtt::model::RayOperand chief;
  chief.coordinate = rtt::model::RayCoordinate::X;
  chief.path = "main";
  chief.surface = SurfaceId("IMG");
  s.optimization.operands = {chief};
  return s;
}

TrialOverlay overlay_of(const System& nominal,
                        rtt::compile::Perturbation perturbation = {},
                        const rtt::material::MaterialLibrary& materials = library(),
                        const rtt::coating::CoatingLibrary* coatings = nullptr) {
  return TrialOverlay{std::move(perturbation),
                      rtt::optim::nominal_basis(nominal, materials, coatings), 1e-3};
}

/// Every operand kind but the paraxial ones, and both generators, on the system's first path in
/// configuration 0: spot about the centroid and the chief ray (polychromatic), OPD, a ray at the
/// image surface.
rtt::model::Optimization battery(const System& s, const rtt::compile::CompiledSystem& cs) {
  rtt::model::Optimization o;
  const std::string& path = s.paths.at(0).name;
  const auto last = static_cast<std::uint16_t>(s.fields.points.size() - 1);
  for (const std::uint16_t field : {std::uint16_t{0}, last}) {
    rtt::model::SpotRmsOperand centroid;
    centroid.path = path;
    centroid.field = field;
    centroid.rings = 3;
    rtt::model::SpotRmsOperand chief = centroid;
    chief.reference = rtt::model::SpotReference::Chief;
    chief.polychromatic = true;
    rtt::model::OpdRmsOperand opd;
    opd.path = path;
    opd.field = field;
    opd.grid = 9;
    rtt::model::RayOperand ray;
    ray.path = path;
    ray.field = field;
    ray.surface = cs.surfaces()[cs.paths().at(0).events.back().surface].id;
    ray.py = 0.7;
    o.operands.insert(o.operands.end(), {centroid, chief, opd, ray});
  }
  rtt::model::SpotGenerator spot;
  spot.path = path;
  spot.reference = rtt::model::SpotReference::Chief;
  rtt::model::WavefrontGenerator wavefront;
  wavefront.path = path;
  o.generators = {spot, wavefront};
  return o;
}

/// The evaluation at the start values, or the name of the exception's type.
std::variant<MeritEvaluation, std::string> outcome(const System& s,
                                                   const rtt::material::MaterialLibrary& materials,
                                                   const rtt::coating::CoatingLibrary& coatings,
                                                   const TrialOverlay* overlay) {
  try {
    const MeritFunction merit(s, materials, &coatings, overlay);
    return merit.evaluate(merit.start());
  } catch (const std::exception& e) {
    return std::string(typeid(e).name());
  }
}

std::vector<fs::path> reference_files() {
  std::vector<fs::path> files;
  for (const auto& entry : fs::recursive_directory_iterator(RTT_REFERENCE_DIR)) {
    const std::string name = entry.path().filename().string();
    if (entry.is_regular_file() && name.ends_with(".rtt.json")) files.push_back(entry.path());
  }
  std::sort(files.begin(), files.end());
  return files;
}

}  // namespace

TEST_CASE("trial overlay: empty perturbation and the basis of the same system are bitwise",
          "[optim][overlay]") {
  // Addendum #256: spot, OPD, rays and both generators aim with the nominal aiming data and take
  // the nominal exit pupil; from the same system these are the same numbers from the same
  // computation (#255, #265), and an empty perturbation compiles bitwise (#254). Every reference
  // system that compiles with the catalogues of the reference files (as test_nominal_reference.cpp
  // in rtt-analysis).
  rtt::material::MaterialLibrary m1;
  m1.add_catalog(std::string(RTT_CATALOG_DIR) + "/schott.agf");
  rtt::material::MaterialLibrary m2;
  m2.add_catalog(std::string(RTT_CATALOG_DIR) + "/m2/schott.agf");
  rtt::coating::CoatingLibrary coatings;
  coatings.add_catalog(std::string(RTT_CATALOG_DIR) + "/coatings/demo.json");
  coatings.add_catalog(std::string(RTT_CATALOG_DIR) + "/coatings/m3.json");
  std::size_t compared = 0;
  for (const fs::path& file : reference_files()) {
    INFO(file.string());
    System s = rtt::io::load_system(file);
    std::optional<rtt::compile::CompiledSystem> cs;
    const rtt::material::MaterialLibrary* materials = nullptr;
    for (const rtt::material::MaterialLibrary* lib : {&m1, &m2}) {
      try {
        cs.emplace(rtt::compile::compile(s, *lib, coatings));
        materials = lib;
        break;
      } catch (const rtt::compile::CompileError& e) {
        static_cast<void>(e);  // another catalogue
      }
    }
    if (materials == nullptr) continue;
    s.optimization = battery(s, *cs);
    const TrialOverlay overlay = overlay_of(s, {}, *materials, &coatings);
    const auto a = outcome(s, *materials, coatings, nullptr);
    const auto b = outcome(s, *materials, coatings, &overlay);
    REQUIRE(a.index() == b.index());
    if (a.index() == 1) {
      CHECK(std::get<1>(a) == std::get<1>(b));
      continue;
    }
    const MeritEvaluation& x = std::get<0>(a);
    const MeritEvaluation& y = std::get<0>(b);
    REQUIRE(x.values.size() == y.values.size());
    REQUIRE(x.residuals.size() == y.residuals.size());
    for (std::size_t i = 0; i < x.values.size(); ++i) CHECK(same(x.values[i], y.values[i]));
    for (std::size_t i = 0; i < x.residuals.size(); ++i) {
      CHECK(same(x.residuals[i], y.residuals[i]));
    }
    CHECK(x.undefined == y.undefined);
    CHECK(x.warnings.size() == y.warnings.size());
    ++compared;
  }
  CHECK(compared >= 20U);
}

TEST_CASE("trial overlay: the perturbation is in every evaluation of optimize",
          "[optim][overlay]") {
  // Acceptance #256: the ideal lens decentered by d = 0.1 mm in x through the overlay. The chief
  // ray (stop centre, axis-parallel) crosses the lens plane at the origin and runs to the focus on
  // the lens axis, (d, 0, 100) (collinear mapping, ADR 0031): on the detector at x_det its local x
  // is d - x_det. Variable x_det, operand ray_x = 0: the run must end at x_det = d, which it can
  // only if the decenter is compiled into every evaluation, the Jacobian's included.
  // Bound (ADR 0036, point 3): the residual r = d - x is linear in x with J^T J = 1 and scaling
  // D = 1; the run stops with ||D h|| <= xtol (||D theta|| + xtol) (ConvergedStep, MNT eq.
  // (3.15b)), and for a linear residual ||D e|| <= ||D h|| (1 + mu / lambda_min), mu the damping at
  // the end; or with r = 0 exactly (ConvergedGradient). So |r_end| <= xtol (d + xtol) (1 + mu).
  const double d = 0.1;
  const System s = ideal_lens();
  rtt::compile::Perturbation p;
  p.nodes["lens"].translation_mm = rtt::math::Vec3(d, 0.0, 0.0);
  const TrialOverlay overlay = overlay_of(s, p);
  const rtt::optim::OptimizeOptions options;
  const rtt::optim::OptimResult r =
      rtt::optim::optimize(s, library(), nullptr, options, {}, overlay);
  REQUIRE(r.variables.size() == 1);
  REQUIRE(!r.history.empty());
  REQUIRE(r.operands.size() == 1);
  const double mu = r.history.back().mu;
  const double bound = options.xtol * (d + options.xtol) * (1.0 + mu);
  INFO("x_det " << r.variables[0].end << ", ray_x " << r.operands[0].value << ", mu " << mu);
  CHECK(std::abs(r.operands[0].value) <= bound);
  CHECK(std::abs(r.variables[0].end - d) <= bound);
  // The result model is the input model with the result value, without the overlay.
  CHECK(r.system.root.children.size() == s.root.children.size());

  // Without the overlay the chief ray meets x = 0, and the detector stays.
  const rtt::optim::OptimResult plain = rtt::optim::optimize(s, library(), nullptr, options);
  CHECK(plain.variables[0].end == 0.0);
}

TEST_CASE("trial overlay: optimize is bitwise the same for 1, 4 and all threads",
          "[optim][overlay]") {
  System s = ideal_lens();
  rtt::model::FirstOrderOperand efl;
  efl.path = "main";
  efl.common.weight = 0.0;  // watched only, through the parabasal rays
  s.optimization.operands.push_back(efl);
  rtt::compile::Perturbation p;
  p.nodes["lens"].translation_mm = rtt::math::Vec3(0.1, -0.05, 0.0);
  p.nodes["lens"].rotation_rad = rtt::math::Vec3(0.01, 0.0, 0.0);
  const TrialOverlay overlay = overlay_of(s, p);
  const auto run = [&](int threads) {
    oneapi::tbb::task_arena arena(threads);
    return arena.execute(
        [&] { return rtt::optim::optimize(s, library(), nullptr, {}, {}, overlay); });
  };
  const rtt::optim::OptimResult one = run(1);
  REQUIRE(one.operands.size() == 2);
  CHECK(std::isfinite(one.operands[1].value));
  for (const int threads : {4, oneapi::tbb::info::default_concurrency()}) {
    INFO(threads << " threads");
    const rtt::optim::OptimResult other = run(threads);
    REQUIRE(other.variables.size() == one.variables.size());
    CHECK(same(other.variables[0].end, one.variables[0].end));
    REQUIRE(other.operands.size() == one.operands.size());
    for (std::size_t i = 0; i < one.operands.size(); ++i) {
      CHECK(same(other.operands[i].value, one.operands[i].value));
    }
    CHECK(other.history.size() == one.history.size());
  }
}

TEST_CASE("trial overlay: a nominal basis of another system is an invalid argument",
          "[optim][overlay]") {
  const System s = ideal_lens();
  // Another path name.
  System renamed = s;
  renamed.paths[0].name = "other";
  renamed.optimization.operands.clear();
  const TrialOverlay other_path = overlay_of(renamed);
  CHECK_THROWS_AS(MeritFunction(s, library(), nullptr, &other_path), std::invalid_argument);
  // Another number of configurations.
  System configured = s;
  configured.configurations = {{"a"}, {"b"}};
  const TrialOverlay other_count = overlay_of(configured);
  CHECK_THROWS_AS(MeritFunction(s, library(), nullptr, &other_count), std::invalid_argument);
  // Another number of wavelengths: checked by the constructor.
  System more = s;
  more.wavelengths.push_back({0.6563, 1.0, false});
  const TrialOverlay fewer = overlay_of(s);
  CHECK_THROWS_AS(MeritFunction(more, library(), nullptr, &fewer), std::invalid_argument);
  // Another wavelength: the aiming data do not fit (trace::AimReference, #255), at the first
  // evaluation.
  System shifted = s;
  shifted.wavelengths[0].um = 0.6563;
  const TrialOverlay other_wavelength = overlay_of(shifted);
  const MeritFunction merit(s, library(), nullptr, &other_wavelength);
  CHECK_THROWS_AS(merit.evaluate(merit.start()), std::invalid_argument);
}
