// Aiming with the paraxial data of another system (ADR 0007, addendum #252; ADR 0036, point 2;
// #255): the entrance pupil, stop scale, sigma and the field conversion come from an
// AimReference of the nominal system, the rays run through the aimed (perturbed) system and
// converge onto its stop. With the reference of the same system the rays are bitwise those
// without reference.
//
// The perturbed systems are built directly with poses (a decentered stop, a decentered and tilted
// lens), which makes them not rotationally symmetric: rtt::paraxial::first_order rejects them.
// Expected values: the aim target sigma (px R_s, py R_s) in the local coordinates of the stop
// surface (aim_ray(), ADR 0007) with R_s of the nominal system, met within kAimTolerance.

#include <oneapi/tbb/task_arena.h>

#include <algorithm>
#include <array>
#include <bit>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <numbers>
#include <numeric>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "rtt/coating/catalog.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
#include "rtt/paraxial/paraxial.hpp"
#include "rtt/trace/apply_event.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/run_control.hpp"
#include "rtt/trace/sequential.hpp"
#include "rtt/trace/sources.hpp"

namespace fs = std::filesystem;
using rtt::compile::compile;
using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::material::MaterialLibrary;
using rtt::math::Vec3;
using rtt::model::Element;
using rtt::model::ElementKind;
using rtt::model::FieldType;
using rtt::model::Param;
using rtt::model::Pose;
using rtt::model::Surface;
using rtt::model::SurfaceId;
using rtt::model::System;
using rtt::trace::Aiming;
using rtt::trace::AimReference;
using rtt::trace::kAimTolerance;
using rtt::trace::RayBatch;
using rtt::trace::RayState;
using rtt::trace::RayStatus;

namespace {

bool same(double a, double b) {
  return std::bit_cast<std::uint64_t>(a) == std::bit_cast<std::uint64_t>(b);
}

/// Asserts that two ray batches are bitwise the same.
void require_same(const RayBatch& a, const RayBatch& b) {
  REQUIRE(a.size() == b.size());
  for (std::size_t i = 0; i < a.size(); ++i) {
    INFO("ray " << i);
    CHECK(same(a.pos_x()[i], b.pos_x()[i]));
    CHECK(same(a.pos_y()[i], b.pos_y()[i]));
    CHECK(same(a.pos_z()[i], b.pos_z()[i]));
    CHECK(same(a.dir_x()[i], b.dir_x()[i]));
    CHECK(same(a.dir_y()[i], b.dir_y()[i]));
    CHECK(same(a.dir_z()[i], b.dir_z()[i]));
    CHECK(a.wl()[i] == b.wl()[i]);
    CHECK(a.field()[i] == b.field()[i]);
    CHECK(same(a.pupil_x()[i], b.pupil_x()[i]));
    CHECK(same(a.pupil_y()[i], b.pupil_y()[i]));
    CHECK(a.status()[i] == b.status()[i]);
  }
}

bool same(const RayState& a, const RayState& b) {
  return same(a.pos.x(), b.pos.x()) && same(a.pos.y(), b.pos.y()) && same(a.pos.z(), b.pos.z()) &&
         same(a.dir.x(), b.dir.x()) && same(a.dir.y(), b.dir.y()) && same(a.dir.z(), b.dir.z()) &&
         a.status == b.status;
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

Surface surface(const std::string& id, double z_mm = 0.0, std::optional<double> radius = {}) {
  Surface s;
  s.id = SurfaceId(id);
  s.pose = Pose::along_z(z_mm);
  if (radius) s.shape.base = rtt::model::Conic{Param(*radius), Param(0.0)};
  return s;
}

/// In vacuum: L1 (CONST:1.5, R 40 / plane, 3 mm) at z = 0, the stop (r = 50 mm, never
/// vignetting) at z = 10, L2 (CONST:1.5, R 30 / plane, 3 mm) at z = 20, detector at z = 60;
/// EPD 4 mm, fields 0 and 3 degree (object at infinity). With `stop_first` the stop is at z = -5
/// in front of L1 (the entrance pupil is the stop: R_s = 2 mm).
System two_lenses(bool stop_first = false) {
  System s;
  s.name = "aim reference";
  s.environment.medium = "VACUUM";
  s.wavelengths = {{0.5876, 1.0, true}, {0.4861, 1.0, false}};
  s.aperture = {rtt::model::SystemApertureType::EntrancePupilDiameter, Param(4.0)};
  s.fields = {FieldType::AngleDeg, {{0.0, 0.0, 1.0}, {0.0, 3.0, 1.0}}};
  s.root.name = "root";
  Surface stop = surface("STO");
  stop.aperture = rtt::model::CircularAperture{50.0, 0.0};
  const Element stop_element{
      "stop", ElementKind::Stop, Pose::along_z(stop_first ? -5.0 : 10.0), std::nullopt, {stop}};
  const Element l1{"L1",
                   ElementKind::Lens,
                   Pose::along_z(0.0),
                   "CONST:1.5",
                   {surface("L1.S1", 0.0, 40.0), surface("L1.S2", 3.0)}};
  const Element l2{"L2",
                   ElementKind::Lens,
                   Pose::along_z(20.0),
                   "CONST:1.5",
                   {surface("L2.S1", 0.0, 30.0), surface("L2.S2", 3.0)}};
  if (stop_first) {
    s.root.children = {{stop_element}, {l1}, {l2}};
  } else {
    s.root.children = {{l1}, {stop_element}, {l2}};
  }
  s.root.children.push_back(
      {Element{"D", ElementKind::Detector, Pose::along_z(60.0), std::nullopt, {surface("IMG")}}});
  s.paths = {{"main", true, {}}};
  return s;
}

/// The element `name` among the children of the root.
Element& element(System& s, const std::string& name) {
  for (auto& child : s.root.children) {
    if (auto* e = std::get_if<Element>(&child.value); e != nullptr && e->name == name) return *e;
  }
  FAIL("no element " << name);
  return std::get<Element>(s.root.children[0].value);
}

/// The stop of `s` decentered by (0.1, -0.05) mm.
System decentered_stop(System s) {
  Pose& pose = element(s, "stop").pose;
  pose.position[0] = Param(0.1);
  pose.position[1] = Param(-0.05);
  return s;
}

const MaterialLibrary& materials() {
  static const MaterialLibrary library;
  return library;
}

CompiledSystem compiled(const System& s) {
  return compile(s, materials());
}

/// Index of the first Stop event of path 0.
std::size_t stop_event(const CompiledSystem& cs) {
  const auto& events = cs.path(PathId{0}).events;
  for (std::size_t i = 0; i < events.size(); ++i) {
    if (cs.surfaces()[events[i].surface].element_kind == ElementKind::Stop) return i;
  }
  FAIL("no stop");
  return 0;
}

/// The ray traced to the stop surface of `cs` (apertures before the stop included), hit in local
/// stop coordinates (as in test_sources.cpp).
Vec3 hit_on_stop(const CompiledSystem& cs, RayState ray, std::uint16_t wl) {
  const auto& events = cs.path(PathId{0}).events;
  const std::size_t stop = stop_event(cs);
  for (std::size_t i = 0; i < stop; ++i) {
    const auto& e = events[i];
    ray = rtt::trace::sequential_step(ray, cs.surfaces()[e.surface], e.surface, e.kind,
                                      cs.media()[e.medium_before].index[wl].real(),
                                      cs.media()[e.medium_after].index[wl].real());
  }
  REQUIRE(ray.status == RayStatus::Alive);
  const auto hit = rtt::trace::intersect_surface(ray, cs.surfaces()[events[stop].surface]);
  REQUIRE(hit.status == rtt::geom::HitStatus::Hit);
  return hit.point;
}

/// Signed R_s of `cs` at `wl`: local y of the hit of the real-aimed rim ray (0, 1) on the stop.
double signed_stop_radius(const CompiledSystem& cs, std::uint16_t wl) {
  const auto rim = rtt::trace::aim_ray(cs, PathId{0}, std::uint16_t{0}, wl, 0.0, 1.0);
  REQUIRE(rim.ray.status == RayStatus::Alive);
  return hit_on_stop(cs, rim.ray, wl).y();
}

constexpr std::array<std::array<double, 2>, 4> kPupil{
    {{0.0, 0.0}, {0.0, 1.0}, {1.0, 0.0}, {-0.7, 0.7}}};

/// Real aiming of `perturbed` with the reference of `nominal` meets sigma p R_s of the nominal
/// R_s on the stop of `perturbed`, for every field, wavelength and the points of kPupil.
void require_on_target(const CompiledSystem& nominal, const CompiledSystem& perturbed) {
  for (std::uint16_t wl = 0; wl < nominal.wavelengths_um().size(); ++wl) {
    const AimReference reference = rtt::trace::aim_reference(nominal, PathId{0}, wl);
    const double r_s = signed_stop_radius(nominal, wl);
    CHECK(std::abs(std::abs(r_s) - rtt::trace::stop_radius(reference)) <= kAimTolerance);
    for (std::uint16_t field = 0; field < nominal.fields().points.size(); ++field) {
      for (const auto& [px, py] : kPupil) {
        INFO("wavelength " << wl << ", field " << field << ", p = (" << px << ", " << py << ")");
        const auto aimed =
            rtt::trace::aim_ray(perturbed, PathId{0}, field, wl, px, py, Aiming::Real, reference);
        REQUIRE(aimed.ray.status == RayStatus::Alive);
        CHECK(aimed.residual < kAimTolerance);
        const Vec3 hit = hit_on_stop(perturbed, aimed.ray, wl);
        CHECK(std::hypot(hit.x() - px * r_s, hit.y() - py * r_s) < kAimTolerance);
      }
    }
  }
}

}  // namespace

TEST_CASE("aim reference: the reference of the same system gives bitwise the same rays",
          "[sources][aim_reference]") {
  // Acceptance "ohne Kontext bitgleich": every reference system that aims, every path with a
  // stop, every wavelength and field, real and paraxial aiming.
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
    std::vector<std::uint16_t> fields(cs->fields().points.size());
    std::iota(fields.begin(), fields.end(), std::uint16_t{0});
    for (std::uint32_t p = 0; p < cs->paths().size(); ++p) {
      for (std::uint16_t wl = 0; wl < cs->wavelengths_um().size(); ++wl) {
        INFO("path " << p << ", wavelength " << wl);
        const PathId path{p};
        std::optional<AimReference> reference;
        try {
          reference = rtt::trace::aim_reference(*cs, path, wl);
        } catch (const std::exception&) {
          // The same paraxial computation fails without reference.
          CHECK_THROWS(rtt::trace::stop_radius(*cs, path, wl));
          continue;
        }
        CHECK(same(rtt::trace::stop_radius(*reference), rtt::trace::stop_radius(*cs, path, wl)));
        for (const Aiming aiming : {Aiming::Real, Aiming::Paraxial}) {
          const rtt::trace::HexapolarPupil sampling{3};
          std::optional<RayBatch> plain;
          try {
            plain = rtt::trace::make_rays(*cs, path, fields, wl, sampling, aiming);
          } catch (const std::invalid_argument&) {
            CHECK_THROWS_AS(
                rtt::trace::make_rays(*cs, path, fields, wl, sampling, aiming, *reference),
                std::invalid_argument);
            continue;
          }
          require_same(*plain,
                       rtt::trace::make_rays(*cs, path, fields, wl, sampling, aiming, *reference));
          for (const auto& field : cs->fields().points) {
            const auto a = rtt::trace::aim_ray(*cs, path, field, wl, 0.3, -0.6, aiming);
            const auto b = rtt::trace::aim_ray(*cs, path, field, wl, 0.3, -0.6, aiming, *reference);
            CHECK(same(a.ray, b.ray));
            CHECK(a.iterations == b.iterations);
            CHECK(same(a.residual, b.residual));
          }
          ++compared;
        }
      }
    }
  }
  CHECK(compared >= 40U);
}

TEST_CASE("aim reference: real aiming converges onto a decentered stop",
          "[sources][aim_reference]") {
  // Acceptance: the stop decentered by (0.1, -0.05) mm; with the nominal reference the chief ray
  // passes through the centre of the perturbed stop and pupil point p meets it at p R_s.
  const CompiledSystem nominal = compiled(two_lenses());
  const CompiledSystem perturbed = compiled(decentered_stop(two_lenses()));
  // Without reference the perturbed system has no paraxial data.
  CHECK_THROWS_AS(rtt::trace::aim_ray(perturbed, PathId{0}, std::uint16_t{1}, 0, 0.0, 0.0),
                  rtt::paraxial::ParaxialError);
  require_on_target(nominal, perturbed);
}

TEST_CASE("aim reference: hand values with the stop in front", "[sources][aim_reference]") {
  // Stop first in vacuum: the entrance pupil is the stop, R_s = EPD / 2 = 2 mm, and the rays
  // run straight to it. Decentered by (0.1, -0.05) mm, the chief ray of 3 degree crosses the stop
  // plane z = -5 at (0.1, -0.05) and the rim ray p = (0, 1) at (0.1, 1.95).
  const CompiledSystem nominal = compiled(two_lenses(true));
  const CompiledSystem perturbed = compiled(decentered_stop(two_lenses(true)));
  const AimReference reference = rtt::trace::aim_reference(nominal, PathId{0}, 0);
  CHECK(std::abs(rtt::trace::stop_radius(reference) - 2.0) <= 1e-15);
  const auto on_stop_plane = [](const RayState& ray) {
    return ray.pos + ((-5.0 - ray.pos.z()) / ray.dir.z()) * ray.dir;
  };
  const auto chief = rtt::trace::aim_ray(perturbed, PathId{0}, std::uint16_t{1}, 0, 0.0, 0.0,
                                         Aiming::Real, reference);
  REQUIRE(chief.ray.status == RayStatus::Alive);
  const Vec3 c = on_stop_plane(chief.ray);
  CHECK(std::abs(c.x() - 0.1) < kAimTolerance);
  CHECK(std::abs(c.y() - -0.05) < kAimTolerance);
  const auto rim = rtt::trace::aim_ray(perturbed, PathId{0}, std::uint16_t{1}, 0, 0.0, 1.0,
                                       Aiming::Real, reference);
  REQUIRE(rim.ray.status == RayStatus::Alive);
  const Vec3 r = on_stop_plane(rim.ray);
  CHECK(std::abs(r.x() - 0.1) < kAimTolerance);
  CHECK(std::abs(r.y() - 1.95) < kAimTolerance);
  // The field direction is that of the nominal field point (3 degree).
  const double t = std::tan(3.0 * std::numbers::pi / 180.0);
  CHECK(std::abs(chief.ray.dir.y() / chief.ray.dir.z() - t) <= 1e-15);
}

TEST_CASE("aim reference: a decentered and tilted lens before the stop",
          "[sources][aim_reference]") {
  System s = two_lenses();
  Pose& pose = element(s, "L1").pose;
  pose.position[0] = Param(0.2);
  pose.rotation_deg[0] = Param(0.5);
  const CompiledSystem nominal = compiled(two_lenses());
  const CompiledSystem perturbed = compiled(s);
  require_on_target(nominal, perturbed);
}

TEST_CASE("aim reference: fields mean the nominal directions and object points",
          "[sources][aim_reference]") {
  SECTION("paraxial image height, object at infinity: the direction of the nominal field") {
    System s = two_lenses();
    s.fields = {FieldType::ParaxialImageHeight, {{0.0, 0.0, 1.0}, {0.0, 1.5, 1.0}}};
    const CompiledSystem nominal = compiled(s);
    const CompiledSystem perturbed = compiled(decentered_stop(s));
    const AimReference reference = rtt::trace::aim_reference(nominal, PathId{0}, 0);
    const auto a = rtt::trace::aim_ray(nominal, PathId{0}, std::uint16_t{1}, 0, 0.0, 0.0);
    const auto b = rtt::trace::aim_ray(perturbed, PathId{0}, std::uint16_t{1}, 0, 0.0, 0.0,
                                       Aiming::Real, reference);
    CHECK(same(a.ray.dir.x(), b.ray.dir.x()));
    CHECK(same(a.ray.dir.y(), b.ray.dir.y()));
    CHECK(same(a.ray.dir.z(), b.ray.dir.z()));
    require_on_target(nominal, perturbed);
  }
  SECTION("object height, finite object: the ray starts in the nominal object point") {
    System s = two_lenses();
    s.object.at_infinity = false;
    s.object.distance = Param(100.0);
    s.fields = {FieldType::ObjectHeight, {{0.0, 0.0, 1.0}, {0.0, 2.0, 1.0}}};
    const CompiledSystem nominal = compiled(s);
    const CompiledSystem perturbed = compiled(decentered_stop(s));
    const AimReference reference = rtt::trace::aim_reference(nominal, PathId{0}, 0);
    const auto b = rtt::trace::aim_ray(perturbed, PathId{0}, std::uint16_t{1}, 0, 0.0, 0.0,
                                       Aiming::Real, reference);
    CHECK(b.ray.pos == Vec3(0.0, 2.0, -100.0));
    require_on_target(nominal, perturbed);
  }
}

TEST_CASE("aim reference: make_rays is bitwise the same for every number of threads",
          "[sources][aim_reference]") {
  const CompiledSystem nominal = compiled(two_lenses());
  const CompiledSystem perturbed = compiled(decentered_stop(two_lenses()));
  const AimReference reference = rtt::trace::aim_reference(nominal, PathId{0}, 1);
  const std::vector<std::uint16_t> fields{0, 1};
  const auto run = [&](int threads, const rtt::trace::RunControl* control) {
    oneapi::tbb::task_arena arena(threads);
    return arena.execute([&] {
      const rtt::trace::HexapolarPupil sampling{4};
      return control != nullptr ? rtt::trace::make_rays(perturbed, PathId{0}, fields, 1, sampling,
                                                        Aiming::Real, reference, *control)
                                : rtt::trace::make_rays(perturbed, PathId{0}, fields, 1, sampling,
                                                        Aiming::Real, reference);
    });
  };
  const RayBatch one = run(1, nullptr);
  for (std::size_t i = 0; i < one.size(); ++i) CHECK(one.status()[i] == RayStatus::Alive);
  require_same(one, run(2, nullptr));
  require_same(one, run(8, nullptr));
  const rtt::trace::RunControl empty;
  require_same(one, run(2, &empty));
}

TEST_CASE("aim reference: a reference that does not fit is an invalid argument",
          "[sources][aim_reference]") {
  const System s = two_lenses();
  const CompiledSystem nominal = compiled(s);
  const AimReference reference = rtt::trace::aim_reference(nominal, PathId{0}, 0);
  const auto throws = [&](const CompiledSystem& system, std::uint16_t wl) {
    CHECK_THROWS_AS(rtt::trace::aim_ray(system, PathId{0}, std::uint16_t{0}, wl, 0.0, 0.0,
                                        Aiming::Real, reference),
                    std::invalid_argument);
  };
  throws(nominal, 1);  // another wavelength than the reference's
  CHECK_THROWS_AS(rtt::trace::aim_ray(nominal, PathId{1}, std::uint16_t{0}, 0, 0.0, 0.0,
                                      Aiming::Real, reference),
                  std::invalid_argument);
  System other = s;
  other.wavelengths[1].um = 0.6563;
  throws(compiled(other), 0);
  other = s;
  other.fields.points[1].y = 4.0;
  throws(compiled(other), 0);
  other = s;
  other.aperture.value = Param(5.0);
  throws(compiled(other), 0);
  other = s;
  other.object.at_infinity = false;
  other.object.distance = Param(1000.0);
  throws(compiled(other), 0);
  other = s;  // another path: an extra plate before the detector (rebuilt in order, no insert)
  const rtt::model::Node detector = other.root.children.back();
  other.root.children.pop_back();
  other.root.children.push_back({Element{"W",
                                         ElementKind::Plate,
                                         Pose::along_z(40.0),
                                         "CONST:1.5",
                                         {surface("W.S1"), surface("W.S2", 1.0)}}});
  other.root.children.push_back(detector);
  throws(compiled(other), 0);
  // The field weight does not aim and does not matter.
  other = s;
  other.fields.points[1].weight = 3.0;
  CHECK_NOTHROW(rtt::trace::aim_ray(compiled(other), PathId{0}, std::uint16_t{1}, 0, 0.0, 0.0,
                                    Aiming::Real, reference));
}
