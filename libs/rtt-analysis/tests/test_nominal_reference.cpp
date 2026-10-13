// Spot and OPD with the nominal reference (ADR 0036, point 2; ADR 0032, addendum #265): the rays
// are aimed with the paraxial data of the nominal system, the reference sphere has its centre on
// the chief ray of the trial and its radius to the exit pupil of the nominal system. With the
// reference of the same system the results are bitwise those without reference.
//
// The perturbed systems are built directly with poses (a decentered ideal lens, a decentered and
// tilted real lens); rtt::paraxial::first_order rejects them.

#include <oneapi/tbb/info.h>
#include <oneapi/tbb/task_arena.h>

#include <algorithm>
#include <bit>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <numbers>
#include <optional>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <utility>
#include <variant>
#include <vector>

#include "rtt/analysis/nominal_reference.hpp"
#include "rtt/analysis/opd.hpp"
#include "rtt/analysis/spot.hpp"
#include "rtt/coating/catalog.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
#include "rtt/paraxial/paraxial.hpp"

namespace fs = std::filesystem;
using rtt::analysis::NominalReference;
using rtt::analysis::OpdMap;
using rtt::analysis::OpdOptions;
using rtt::analysis::OpdPupilPoints;
using rtt::analysis::SpotDiagram;
using rtt::analysis::SpotOptions;
using rtt::compile::compile;
using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::material::MaterialLibrary;
using rtt::model::Element;
using rtt::model::ElementKind;
using rtt::model::FieldType;
using rtt::model::Param;
using rtt::model::Pose;
using rtt::model::Surface;
using rtt::model::SurfaceId;
using rtt::model::System;

namespace {

bool same(double a, double b) {
  return std::bit_cast<std::uint64_t>(a) == std::bit_cast<std::uint64_t>(b);
}

void require_same(const SpotDiagram& a, const SpotDiagram& b) {
  REQUIRE(a.points.size() == b.points.size());
  for (std::size_t i = 0; i < a.points.size(); ++i) {
    CHECK(same(a.points[i].x, b.points[i].x));
    CHECK(same(a.points[i].y, b.points[i].y));
    CHECK(same(a.points[i].weight, b.points[i].weight));
    CHECK(a.points[i].wavelength == b.points[i].wavelength);
  }
  CHECK(same(a.chief.x, b.chief.x));
  CHECK(same(a.chief.y, b.chief.y));
  CHECK(same(a.stats.rms_centroid, b.stats.rms_centroid));
  CHECK(same(a.stats.rms_chief, b.stats.rms_chief));
  CHECK(a.rays_launched == b.rays_launched);
  CHECK(a.rays_arrived == b.rays_arrived);
  CHECK(a.warnings.size() == b.warnings.size());
}

void require_same(const rtt::analysis::ReferenceSphere& a,
                  const rtt::analysis::ReferenceSphere& b) {
  CHECK(same(a.centre.x(), b.centre.x()));
  CHECK(same(a.centre.y(), b.centre.y()));
  CHECK(same(a.centre.z(), b.centre.z()));
  CHECK(same(a.radius, b.radius));
}

void require_same(const std::vector<rtt::analysis::OpdPoint>& a,
                  const std::vector<rtt::analysis::OpdPoint>& b) {
  REQUIRE(a.size() == b.size());
  for (std::size_t i = 0; i < a.size(); ++i) {
    CHECK(same(a[i].w, b[i].w));
    CHECK(same(a[i].px, b[i].px));
    CHECK(same(a[i].py, b[i].py));
    CHECK(a[i].status == b[i].status);
  }
}

void require_same(const OpdMap& a, const OpdMap& b) {
  require_same(a.sphere, b.sphere);
  require_same(a.points, b.points);
  CHECK(same(a.rms, b.rms));
  CHECK(same(a.pv, b.pv));
  CHECK(a.arrived == b.arrived);
  CHECK(a.warnings.size() == b.warnings.size());
}

void require_same(const OpdPupilPoints& a, const OpdPupilPoints& b) {
  require_same(a.sphere, b.sphere);
  require_same(a.points, b.points);
  CHECK(a.warnings.size() == b.warnings.size());
}

/// Runs `f`; the result, or the name of the type of the exception it threw.
template <class F>
auto outcome(F f) -> std::variant<decltype(f()), std::string> {
  try {
    return f();
  } catch (const std::exception& e) {
    return std::string(typeid(e).name());
  }
}

/// Asserts that the analysis with and without reference ends the same: both bitwise the same
/// result, or both an exception of the same type. True if both gave a result.
template <class F, class G>
bool same_outcome(F without, G with) {
  const auto a = outcome(without);
  const auto b = outcome(with);
  REQUIRE(a.index() == b.index());
  if (a.index() == 1) {
    CHECK(std::get<1>(a) == std::get<1>(b));
    return false;
  }
  require_same(std::get<0>(a), std::get<0>(b));
  return true;
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

Surface plane(const std::string& id, double z_mm = 0.0, std::optional<double> radius = {}) {
  Surface s;
  s.id = SurfaceId(id);
  s.pose = Pose::along_z(z_mm);
  if (radius) s.shape.base = rtt::model::Conic{Param(*radius), Param(0.0)};
  return s;
}

/// In vacuum: the stop (r = 10 mm) at z = -5, an ideal lens f = 100 mm (design conjugate at
/// infinity, ADR 0031) at z = 0, the detector at z = 100 (rear focal plane); object at infinity,
/// EPD 10 mm, fields 0 and 2 degree. With `decenter` the lens is shifted by 0.1 mm along x.
System ideal_lens(bool decenter = false) {
  System s;
  s.name = "decentered ideal lens";
  s.environment.medium = "VACUUM";
  s.wavelengths = {{0.5876, 1.0, true}};
  s.aperture = {rtt::model::SystemApertureType::EntrancePupilDiameter, Param(10.0)};
  s.fields = {FieldType::AngleDeg, {{0.0, 0.0, 1.0}, {0.0, 2.0, 1.0}}};
  s.root.name = "root";
  Surface stop = plane("STO");
  stop.aperture = rtt::model::CircularAperture{10.0, 0.0};
  Surface lens = plane("IL");
  lens.aperture = rtt::model::CircularAperture{30.0, 0.0};
  lens.interaction = rtt::model::IdealLens{Param(100.0), std::nullopt};
  Pose lens_pose = Pose::along_z(0.0);
  if (decenter) lens_pose.position[0] = Param(0.1);
  s.root.children = {
      {Element{"stop", ElementKind::Stop, Pose::along_z(-5.0), std::nullopt, {stop}}},
      {Element{"lens", ElementKind::ThinElement, lens_pose, std::nullopt, {lens}}},
      {Element{"D", ElementKind::Detector, Pose::along_z(100.0), std::nullopt, {plane("IMG")}}}};
  s.paths = {{"main", true, {}}};
  return s;
}

/// Two lenses in CONST:1.5 with the stop between them (as in rtt-trace test_aim_reference.cpp);
/// with `perturb` L1 is decentered by 0.2 mm in x and tilted by 0.5 degree about x.
System two_lenses(bool perturb = false) {
  System s;
  s.name = "two lenses";
  s.environment.medium = "VACUUM";
  s.wavelengths = {{0.5876, 1.0, true}, {0.4861, 1.0, false}};
  s.aperture = {rtt::model::SystemApertureType::EntrancePupilDiameter, Param(4.0)};
  s.fields = {FieldType::AngleDeg, {{0.0, 0.0, 1.0}, {0.0, 3.0, 1.0}}};
  s.root.name = "root";
  Surface stop = plane("STO");
  stop.aperture = rtt::model::CircularAperture{50.0, 0.0};
  Pose l1_pose = Pose::along_z(0.0);
  if (perturb) {
    l1_pose.position[0] = Param(0.2);
    l1_pose.rotation_deg[0] = Param(0.5);
  }
  s.root.children = {
      {Element{"L1",
               ElementKind::Lens,
               l1_pose,
               "CONST:1.5",
               {plane("L1.S1", 0.0, 40.0), plane("L1.S2", 3.0)}}},
      {Element{"stop", ElementKind::Stop, Pose::along_z(10.0), std::nullopt, {stop}}},
      {Element{"L2",
               ElementKind::Lens,
               Pose::along_z(20.0),
               "CONST:1.5",
               {plane("L2.S1", 0.0, 30.0), plane("L2.S2", 3.0)}}},
      {Element{"D", ElementKind::Detector, Pose::along_z(60.0), std::nullopt, {plane("IMG")}}}};
  s.paths = {{"main", true, {}}};
  return s;
}

CompiledSystem compiled(const System& s) {
  return compile(s, MaterialLibrary{});
}

SpotOptions spot_with(const NominalReference& reference) {
  SpotOptions o;
  o.reference = reference;
  return o;
}

OpdOptions opd_with(const NominalReference& reference) {
  OpdOptions o;
  o.reference = reference;
  return o;
}

}  // namespace

TEST_CASE("nominal reference: the reference of the same system gives bitwise the same results",
          "[analysis][nominal_reference]") {
  // Acceptance: every reference system that compiles, every path, field and wavelength; spot
  // per wavelength and polychromatic, OPD map and OPD points (Gauss sampling).
  MaterialLibrary m1;
  m1.add_catalog(std::string(RTT_CATALOG_DIR) + "/schott.agf");
  MaterialLibrary m2;
  m2.add_catalog(std::string(RTT_CATALOG_DIR) + "/m2/schott.agf");
  rtt::coating::CoatingLibrary coatings;
  coatings.add_catalog(std::string(RTT_CATALOG_DIR) + "/coatings/demo.json");
  coatings.add_catalog(std::string(RTT_CATALOG_DIR) + "/coatings/m3.json");
  std::size_t compared = 0;
  for (const fs::path& file : reference_files()) {
    INFO(file.string());
    const System s = rtt::io::load_system(file);
    std::optional<CompiledSystem> cs;
    for (const MaterialLibrary* lib : {&m1, &m2}) {
      try {
        cs = compile(s, *lib, coatings);
        break;
      } catch (const rtt::compile::CompileError&) {
      }
    }
    if (!cs) continue;
    for (std::uint32_t p = 0; p < cs->paths().size(); ++p) {
      const PathId path{p};
      const NominalReference reference = rtt::analysis::nominal_reference(*cs, path);
      OpdOptions opd;
      opd.grid = 9;
      OpdOptions opd_ref = opd;
      opd_ref.reference = reference;
      for (std::uint16_t field = 0; field < cs->fields().points.size(); ++field) {
        INFO("path " << p << ", field " << field);
        std::vector<std::optional<std::uint16_t>> spots{std::nullopt};
        for (std::uint16_t wl = 0; wl < cs->wavelengths_um().size(); ++wl) spots.emplace_back(wl);
        for (const auto& wl : spots) {
          if (same_outcome([&] { return rtt::analysis::spot(*cs, path, field, wl); },
                           [&] {
                             return rtt::analysis::spot(*cs, path, field, wl, spot_with(reference));
                           })) {
            ++compared;
          }
        }
        for (std::uint16_t wl = 0; wl < cs->wavelengths_um().size(); ++wl) {
          INFO("wavelength " << wl);
          if (same_outcome([&] { return rtt::analysis::opd_map(*cs, path, field, wl, opd); },
                           [&] { return rtt::analysis::opd_map(*cs, path, field, wl, opd_ref); })) {
            ++compared;
          }
          const rtt::trace::GaussPupil gauss{3, 6};
          same_outcome(
              [&] { return rtt::analysis::opd_points(*cs, path, field, wl, gauss, opd); },
              [&] { return rtt::analysis::opd_points(*cs, path, field, wl, gauss, opd_ref); });
        }
      }
    }
  }
  CHECK(compared >= 100U);
}

TEST_CASE("nominal reference: a decentered ideal lens, hand values",
          "[analysis][nominal_reference]") {
  // The ideal lens maps collinearly (ADR 0031): a parallel bundle of direction d meets the rear
  // focal plane at the lens axis plus f (d_x/d_z, d_y/d_z). Decentered by 0.1 mm, that point moves
  // rigidly: the spots of 0 and 2 degree are points at (0.1, 0) and (0.1, 100 tan 2 deg) on the
  // detector. Tolerance 1e-9 mm as in test_ideal_lens_analysis.cpp (more than 100 times the
  // rounding).
  const CompiledSystem nominal = compiled(ideal_lens());
  const CompiledSystem perturbed = compiled(ideal_lens(true));
  const NominalReference reference = rtt::analysis::nominal_reference(nominal, PathId{0});
  // Without reference the perturbed system has no paraxial data.
  CHECK_THROWS_AS(rtt::analysis::spot(perturbed, PathId{0}, 1, std::uint16_t{0}),
                  rtt::paraxial::ParaxialError);
  CHECK_THROWS_AS(rtt::analysis::opd_map(perturbed, PathId{0}, 1, 0), rtt::paraxial::ParaxialError);

  const double t = std::tan(2.0 * std::numbers::pi / 180.0);
  for (const std::uint16_t field : {std::uint16_t{0}, std::uint16_t{1}}) {
    INFO("field " << field);
    const double y = field == 0 ? 0.0 : 100.0 * t;
    const SpotDiagram spot =
        rtt::analysis::spot(perturbed, PathId{0}, field, std::uint16_t{0}, spot_with(reference));
    CHECK(spot.rays_arrived > 100);  // content guard: hexapolar 6 rings
    CHECK(std::abs(spot.chief.x - 0.1) < 1e-9);
    CHECK(std::abs(spot.chief.y - y) < 1e-9);
    CHECK(spot.stats.rms_chief < 1e-9);
    CHECK(spot.stats.geo_chief < 1e-9);

    // OPD: every ray meets C, so the wavefront is a sphere about C and every sphere about C is
    // parallel to it, whatever its radius (here to the nominal exit pupil): W = 0 exactly up to
    // rounding. OPL about 200 mm, a few roundings per event: |dOPL| <~ 10 * 2.2e-16 * 200 mm =
    // 4.4e-13 mm, over lambda = 5.876e-4 mm about 7.5e-10 waves; bound 1e-8 waves. The tilt is
    // zero too, since C is the point of convergence.
    const OpdMap map = rtt::analysis::opd_map(perturbed, PathId{0}, field, 0, opd_with(reference));
    CHECK(map.arrived > 100);
    double worst = 0.0;
    for (const auto& q : map.points) {
      if (q.status == rtt::trace::RayStatus::Alive) worst = std::max(worst, std::abs(q.w));
    }
    CHECK(worst <= 1e-8);
    // The sphere: centre on the chief ray of the trial (the spot's chief on the detector at
    // z = 100), radius to the exit pupil of the nominal system (rtt-paraxial).
    CHECK(std::abs(map.sphere.centre.x() - 0.1) < 1e-9);
    CHECK(std::abs(map.sphere.centre.y() - y) < 1e-9);
    CHECK(std::abs(map.sphere.centre.z() - 100.0) < 1e-9);
    const auto fo = rtt::paraxial::first_order(nominal, PathId{0}, 0);
    REQUIRE(fo.exit_pupil);
    REQUIRE(fo.exit_pupil->z);
    const rtt::math::Vec3 e(0.0, 0.0, *fo.exit_pupil->z);
    CHECK(reference.exit_pupil() == e);
    CHECK(same(map.sphere.radius, (map.sphere.centre - e).norm()));
  }
}

TEST_CASE("nominal reference: a decentered and tilted real lens runs",
          "[analysis][nominal_reference]") {
  const CompiledSystem nominal = compiled(two_lenses());
  const CompiledSystem perturbed = compiled(two_lenses(true));
  const NominalReference reference = rtt::analysis::nominal_reference(nominal, PathId{0});
  CHECK_THROWS_AS(rtt::analysis::opd_map(perturbed, PathId{0}, 1, 0), rtt::paraxial::ParaxialError);
  for (const std::uint16_t field : {std::uint16_t{0}, std::uint16_t{1}}) {
    INFO("field " << field);
    const SpotDiagram spot =
        rtt::analysis::spot(perturbed, PathId{0}, field, std::nullopt, spot_with(reference));
    CHECK(spot.rays_arrived > 100);
    for (std::uint16_t wl = 0; wl < 2; ++wl) {
      const OpdMap map =
          rtt::analysis::opd_map(perturbed, PathId{0}, field, wl, opd_with(reference));
      CHECK(map.arrived > 50);
      // The chief ray of the same wavelength is the reference of W (W(0, 0) = 0, #29).
      const OpdPupilPoints centre =
          rtt::analysis::opd_points(perturbed, PathId{0}, field, wl,
                                    rtt::trace::SinglePupilPoint{0.0, 0.0}, opd_with(reference));
      REQUIRE(centre.points.size() == 1);
      CHECK(centre.points[0].status == rtt::trace::RayStatus::Alive);
      CHECK(centre.points[0].w == 0.0);
      require_same(centre.sphere, map.sphere);
    }
  }
}

TEST_CASE("nominal reference: bitwise the same for 1, 4 and all threads",
          "[analysis][nominal_reference]") {
  const CompiledSystem nominal = compiled(two_lenses());
  const CompiledSystem perturbed = compiled(two_lenses(true));
  const NominalReference reference = rtt::analysis::nominal_reference(nominal, PathId{0});
  const auto run = [&](int threads) {
    oneapi::tbb::task_arena arena(threads);
    return arena.execute([&] {
      return std::pair{
          rtt::analysis::spot(perturbed, PathId{0}, 1, std::nullopt, spot_with(reference)),
          rtt::analysis::opd_points(perturbed, PathId{0}, 1, 1, rtt::trace::GaussPupil{4, 8},
                                    opd_with(reference))};
    });
  };
  const auto one = run(1);
  for (const int threads : {4, oneapi::tbb::info::default_concurrency()}) {
    INFO(threads << " threads");
    const auto other = run(threads);
    require_same(one.first, other.first);
    require_same(one.second, other.second);
  }
}

TEST_CASE("nominal reference: a reference that does not fit is an invalid argument",
          "[analysis][nominal_reference]") {
  const CompiledSystem nominal = compiled(two_lenses());
  const NominalReference reference = rtt::analysis::nominal_reference(nominal, PathId{0});
  // Another path than the reference's.
  System two_paths = two_lenses();
  two_paths.paths.push_back({"second", true, {}});
  const CompiledSystem other_path = compiled(two_paths);
  CHECK_THROWS_AS(
      rtt::analysis::spot(other_path, PathId{1}, 0, std::uint16_t{0}, spot_with(reference)),
      std::invalid_argument);
  CHECK_THROWS_AS(rtt::analysis::opd_map(other_path, PathId{1}, 0, 0, opd_with(reference)),
                  std::invalid_argument);
  // Other wavelengths: the aiming data do not fit (trace::AimReference, #255).
  System other = two_lenses();
  other.wavelengths[1].um = 0.6563;
  const CompiledSystem other_wl = compiled(other);
  CHECK_THROWS_AS(
      rtt::analysis::spot(other_wl, PathId{0}, 0, std::uint16_t{0}, spot_with(reference)),
      std::invalid_argument);
  CHECK_THROWS_AS(rtt::analysis::opd_map(other_wl, PathId{0}, 0, 0, opd_with(reference)),
                  std::invalid_argument);
  CHECK_THROWS_AS(rtt::analysis::nominal_reference(nominal, PathId{5}), std::invalid_argument);
}
