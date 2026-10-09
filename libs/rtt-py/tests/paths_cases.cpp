// Path evaluation (#122), ghost generator (#123) and ghost ranking (#124) cases of
// rtt_py_reference (#133): the same calls as CASES in test_bitwise_paths.py, every result
// written as named .npy arrays (the Python side flattens its results with the same names).

#include "paths_cases.hpp"

#include <oneapi/tbb/task_arena.h>

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <optional>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "npy_writer.hpp"
#include "rtt/analysis/ghosts.hpp"
#include "rtt/analysis/paths.hpp"
#include "rtt/coating/catalog.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/compile/ghosts.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/sources.hpp"

namespace rtt::py::reference {
namespace {

namespace fs = std::filesystem;
using namespace rtt::analysis;
using compile::CompiledSystem;
using compile::PathId;
using trace::RayBatch;

constexpr std::uint64_t kNone = std::numeric_limits<std::uint64_t>::max();
constexpr double kNan = std::numeric_limits<double>::quiet_NaN();  // None in float arrays

/// Writes the arrays of one case as <out>/<case>.<name>.npy.
class Writer {
 public:
  Writer(fs::path out, std::string name) : out_(std::move(out)), name_(std::move(name)) {}

  void f8(const std::string& array, const std::vector<double>& v) const {
    write_column<double>(file(array), "<f8", v);
  }
  void u8(const std::string& array, const std::vector<std::uint64_t>& v) const {
    write_column<std::uint64_t>(file(array), "<u8", v);
  }
  template <typename Points, typename Get>
  void f8_of(const std::string& array, const Points& points, Get get) const {
    std::vector<double> v;
    for (const auto& p : points) v.push_back(get(p));
    f8(array, v);
  }
  template <typename Points, typename Get>
  void u8_of(const std::string& array, const Points& points, Get get) const {
    std::vector<std::uint64_t> v;
    for (const auto& p : points) v.push_back(static_cast<std::uint64_t>(get(p)));
    u8(array, v);
  }

 private:
  [[nodiscard]] fs::path file(const std::string& array) const {
    return out_ / (name_ + "." + array + ".npy");
  }
  fs::path out_;
  std::string name_;
};

// Flattening of the results; the names match the flatten_* functions in
// test_bitwise_paths.py. Integers (indices, counts, status) as uint64, None as kNone or NaN.

/// RayLosses as one array: launched, by_status..., worst_surface (kNone), worst_count.
std::vector<std::uint64_t> losses_values(const RayLosses& l) {
  std::vector<std::uint64_t> v{l.launched};
  for (const std::size_t n : l.by_status) v.push_back(n);
  v.push_back(l.worst_surface ? *l.worst_surface : kNone);
  v.push_back(l.worst_surface_count);
  return v;
}

void write_transmission(const Writer& w, const PathTransmission& t) {
  w.f8_of("px", t.rays, [](const PathRay& r) { return r.px; });
  w.f8_of("py", t.rays, [](const PathRay& r) { return r.py; });
  w.f8_of("weight", t.rays, [](const PathRay& r) { return r.weight; });
  w.u8_of("status", t.rays, [](const PathRay& r) { return r.status; });
  w.f8("scalars", {t.mean, t.min, t.max});
  w.u8("ints", {t.path.index, t.image_surface, t.rays_launched, t.rays_arrived, t.warnings.size()});
  w.u8("losses", losses_values(t.losses));
}

void write_difference(const Writer& w, const PathOplDifference& d) {
  w.f8_of("px", d.points, [](const OplDifferencePoint& p) { return p.px; });
  w.f8_of("py", d.points, [](const OplDifferencePoint& p) { return p.py; });
  w.f8_of("delta", d.points, [](const OplDifferencePoint& p) { return p.delta; });
  w.u8_of("status", d.points, [](const OplDifferencePoint& p) { return p.status; });
  w.f8("scalars", {d.chief.value_or(kNan)});
  w.u8("ints", {d.path_a.index, d.path_b.index, d.image_surface, d.warnings.size()});
  w.u8("losses_a", losses_values(d.losses_a));
  w.u8("losses_b", losses_values(d.losses_b));
}

void write_ghosts(const Writer& w, const compile::GhostSystem& g) {
  using compile::GhostPath;
  w.u8_of("path", g.ghosts, [](const GhostPath& p) { return p.path.index; });
  w.u8_of("base", g.ghosts, [](const GhostPath& p) { return p.base.index; });
  w.u8_of("surface_j", g.ghosts, [](const GhostPath& p) { return p.surface_j; });
  w.u8_of("surface_i", g.ghosts, [](const GhostPath& p) { return p.surface_i; });
  w.u8_of("event_j", g.ghosts, [](const GhostPath& p) { return p.event_j; });
  w.u8_of("event_i", g.ghosts, [](const GhostPath& p) { return p.event_i; });
  w.u8("ints", {g.system.paths().size()});
}

void write_ranking(const Writer& w, const GhostRanking& r) {
  w.u8_of("path", r.entries, [](const GhostEntry& e) { return e.path.index; });
  w.u8_of("surface_j", r.entries, [](const GhostEntry& e) { return e.surface_j; });
  w.u8_of("surface_i", r.entries, [](const GhostEntry& e) { return e.surface_i; });
  w.f8_of("power", r.entries, [](const GhostEntry& e) { return e.power; });
  w.f8_of("relative_power", r.entries, [](const GhostEntry& e) { return e.relative_power; });
  w.f8_of("rms_radius", r.entries, [](const GhostEntry& e) { return e.rms_radius; });
  w.u8_of("rays_arrived", r.entries, [](const GhostEntry& e) { return e.rays_arrived; });
  w.f8_of("relative_irradiance", r.entries,
          [](const GhostEntry& e) { return e.relative_irradiance; });
  w.f8_of("focus_offset", r.entries,
          [](const GhostEntry& e) { return e.focus_offset.value_or(kNan); });
  w.f8_of("paraxial_blur_radius", r.entries,
          [](const GhostEntry& e) { return e.paraxial_blur_radius.value_or(kNan); });
  std::vector<std::uint64_t> losses;
  for (const GhostEntry& e : r.entries) {
    for (const std::uint64_t n : losses_values(e.losses)) losses.push_back(n);
  }
  w.u8("losses", losses);
  w.f8("scalars", {r.base_power, r.base_rms_radius, r.resolution_radius});
  w.u8("ints", {r.base.index, r.field, r.wavelength, r.warnings.size()});
}

/// Collimated start bundle along +z at z = 0, as collimated_bundle() in
/// test_bitwise_paths.py: the axis ray first, then a square grid with a pitch of 0.5 mm inside
/// r <= 2.25 mm; pupil coordinates r / 2.5.
RayBatch collimated_bundle() {
  std::vector<std::pair<double, double>> points{{0.0, 0.0}};
  for (int i = -4; i <= 4; ++i) {
    for (int j = -4; j <= 4; ++j) {
      const double x = 0.5 * i;
      const double y = 0.5 * j;
      if ((i != 0 || j != 0) && x * x + y * y <= 2.25 * 2.25) points.emplace_back(x, y);
    }
  }
  RayBatch rays(points.size());
  for (std::size_t k = 0; k < points.size(); ++k) {
    rays.pos_x()[k] = points[k].first;
    rays.pos_y()[k] = points[k].second;
    rays.pupil_x()[k] = points[k].first / 2.5;
    rays.pupil_y()[k] = points[k].second / 2.5;
  }
  return rays;
}

/// Systems of the cases, as in systems() of test_bitwise_paths.py.
struct Systems {
  CompiledSystem michelson;
  CompiledSystem michelson_small;  // test mirror radius 1 mm: most rays end vignetted at M2
  CompiledSystem singlet;
  compile::GhostSystem plate;
  compile::GhostSystem cooke;
  compile::GhostSystem ghost_plates;   // M4 acceptance (#135)
  compile::GhostSystem ghost_singlet;  // M4 acceptance (#135)
  CompiledSystem tour;                 // feature tour, refract with order at G.S1 (#135)
};

model::Element& element(model::System& s, std::size_t child) {
  return std::get<model::Element>(s.root.children[child].value);
}

Systems load_systems(const fs::path& reference_dir, const fs::path& catalog_dir) {
  const material::MaterialLibrary plain;
  const coating::CoatingLibrary no_coatings;
  material::MaterialLibrary schott;
  schott.add_catalog(catalog_dir / "m2" / "schott.agf");
  material::MaterialLibrary schott_bk7;  // N-BK7 of the feature tour
  schott_bk7.add_catalog(catalog_dir / "schott.agf");
  // The feature tour shows every format element and does not compile as it is: without the
  // Zernike term of A.S1 (M8) and with a bare Fresnel A.S1 instead of the coating "AR_VIS"
  // (no catalog), as in test_m4_acceptance.cpp.
  model::System tour = io::load_system(reference_dir / "m0/feature_tour.rtt.json");
  element(tour, 2).surfaces[0].shape.terms.clear();
  element(tour, 2).surfaces[0].interaction = model::Fresnel{};
  const model::System michelson = io::load_system(reference_dir / "m4/michelson_offset.rtt.json");
  model::System small = michelson;
  element(small, 3).surfaces[0].aperture = model::CircularAperture{1.0, 0.0};  // test mirror
  return {compile::compile(michelson, plain),
          compile::compile(small, plain),
          compile::compile(io::load_system(reference_dir / "m1/singlet_const.rtt.json"), plain),
          compile::compile_with_ghosts(io::load_system(reference_dir / "m3/fresnel_bk7.rtt.json"),
                                       "main", plain, no_coatings),
          compile::compile_with_ghosts(io::load_system(reference_dir / "m2/cooke_triplet.rtt.json"),
                                       "main", schott, no_coatings),
          compile::compile_with_ghosts(io::load_system(reference_dir / "m4/ghost_plates.rtt.json"),
                                       "main", plain, no_coatings),
          compile::compile_with_ghosts(io::load_system(reference_dir / "m4/ghost_singlet.rtt.json"),
                                       "main", plain, no_coatings),
          compile::compile(tour, schott_bk7)};
}

PathId path(const CompiledSystem& s, const char* name) {
  return *s.find_path(name);
}

}  // namespace

void run_paths_cases(const fs::path& reference_dir,
                     const fs::path& catalog_dir,
                     const fs::path& out,
                     int threads) {
  const Systems s = load_systems(reference_dir, catalog_dir);
  const auto w = [&](const std::string& name) { return Writer(out, name); };
  const RayBatch start = collimated_bundle();
  const char* const arms[] = {"reference arm", "test arm", "reference arm, return",
                              "test arm, return"};

  oneapi::tbb::task_arena arena(threads);
  arena.execute([&] {
    for (std::size_t k = 0; k < 4; ++k) {
      write_transmission(w("michelson_t" + std::to_string(k)),
                         path_transmission(s.michelson, path(s.michelson, arms[k]), start));
    }
    write_difference(w("michelson_opl"),
                     opl_difference(s.michelson, path(s.michelson, "reference arm"),
                                    path(s.michelson, "test arm"), start));
    write_transmission(
        w("michelson_small_t"),
        path_transmission(s.michelson_small, path(s.michelson_small, "test arm"), start));
    write_difference(w("michelson_small_opl"),
                     opl_difference(s.michelson_small, path(s.michelson_small, "test arm"),
                                    path(s.michelson_small, "reference arm"), start));
    PathOptions hex6;
    hex6.sampling = trace::HexapolarPupil{6};
    write_transmission(w("singlet_t"), path_transmission(s.singlet, PathId{0}, 1, 0, hex6));
    PathOptions grid9;
    grid9.sampling = trace::GridPupil{9};
    grid9.aiming = trace::Aiming::Paraxial;
    write_difference(w("singlet_opl"),
                     opl_difference(s.singlet, PathId{0}, PathId{0}, 2, 0, grid9));
    write_ghosts(w("plate_ghosts"), s.plate);
    GhostRankingOptions plate_options;
    plate_options.sampling = trace::HexapolarPupil{4};
    write_ranking(w("plate_ranking"),
                  ghost_ranking(s.plate, 0, s.plate.system.reference_wavelength(), plate_options));
    write_ghosts(w("cooke_ghosts"), s.cooke);
    GhostRankingOptions cooke_options;
    cooke_options.resolution_radius = 0.01;
    write_ranking(w("cooke_ranking"),
                  ghost_ranking(s.cooke, 1, s.cooke.system.reference_wavelength(), cooke_options));
    // M4 acceptance (#135): the ghost reference files with the default options, and the path
    // "first order" of the feature tour (refract with order 1 at G.S1) for the start rays. Its
    // hand-made ghost path loses every ray of this bundle at A.S2, so it is not compared.
    write_ghosts(w("ghost_plates_ghosts"), s.ghost_plates);
    write_ranking(w("ghost_plates_ranking"),
                  ghost_ranking(s.ghost_plates, 0, s.ghost_plates.system.reference_wavelength()));
    write_ghosts(w("ghost_singlet_ghosts"), s.ghost_singlet);
    write_ranking(w("ghost_singlet_ranking"),
                  ghost_ranking(s.ghost_singlet, 0, s.ghost_singlet.system.reference_wavelength()));
    write_transmission(w("tour_first_order"),
                       path_transmission(s.tour, path(s.tour, "first order"), start));
  });
}

}  // namespace rtt::py::reference
