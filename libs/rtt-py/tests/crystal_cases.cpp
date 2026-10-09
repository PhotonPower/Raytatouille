// Crystal cases of rtt_py_reference (#134, ADR 0026): the same calls as CASES in
// test_bitwise_crystal.py. Every case traces start rays set by hand through the calcite plate of
// tests/reference/m4/calcite_walkoff.rtt.json, on its paths "o" and "e" to the detector and on
// two paths added here that end right after the entry ("o entry", "e entry"), where the rays
// are inside the crystal (wave != dir, mode_index > 0). The traced columns including wave and
// mode_index, the PRT, the trace statistics and the path transmission are written as named .npy
// arrays (the Python side flattens its results with the same names).

#include "crystal_cases.hpp"

#include <oneapi/tbb/task_arena.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "npy_writer.hpp"
#include "rtt/analysis/paths.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/sequential.hpp"

namespace rtt::py::reference {
namespace {

namespace fs = std::filesystem;
using compile::CompiledSystem;
using compile::PathId;
using trace::RayBatch;

constexpr std::uint64_t kNone = std::numeric_limits<std::uint64_t>::max();

/// Writes the arrays of one case as <out>/<case>.<name>.npy.
class Writer {
 public:
  Writer(fs::path out, std::string name) : out_(std::move(out)), name_(std::move(name)) {}

  template <typename T>
  void column(const std::string& array, const char* descr, std::span<const T> v) const {
    write_column<T>(out_ / (name_ + "." + array + ".npy"), descr, v);
  }
  void u8(const std::string& array, const std::vector<std::uint64_t>& v) const {
    column<std::uint64_t>(array, "<u8", v);
  }
  void f8(const std::string& array, const std::vector<double>& v) const {
    column<double>(array, "<f8", v);
  }

 private:
  fs::path out_;
  std::string name_;
};

/// The columns of a traced batch, as flatten() in test_bitwise_crystal.py.
void write_rays(const Writer& w, const RayBatch& r, const trace::TraceStats& stats) {
  w.column<double>("pos_x", "<f8", r.pos_x());
  w.column<double>("pos_y", "<f8", r.pos_y());
  w.column<double>("pos_z", "<f8", r.pos_z());
  w.column<double>("dir_x", "<f8", r.dir_x());
  w.column<double>("dir_y", "<f8", r.dir_y());
  w.column<double>("dir_z", "<f8", r.dir_z());
  w.column<double>("wave_x", "<f8", r.wave_x());
  w.column<double>("wave_y", "<f8", r.wave_y());
  w.column<double>("wave_z", "<f8", r.wave_z());
  w.column<double>("mode_index", "<f8", r.mode_index());
  w.column<double>("opl", "<f8", r.opl());
  w.column<double>("weight", "<f8", r.weight());
  w.column<std::uint32_t>("last_surface", "<u4", r.last_surface());
  w.column<trace::RayStatus>("status", "|u1", r.status());
  for (int row = 0; row < 3; ++row) {
    for (int col = 0; col < 3; ++col) {
      w.column<math::Complex>("prt" + std::to_string(row) + std::to_string(col), "<c16",
                              r.prt(row, col));
    }
  }
  w.u8("stats", std::vector<std::uint64_t>(stats.rays.begin(), stats.rays.end()));
}

/// path_transmission of the same start rays, as flatten() in test_bitwise_crystal.py.
void write_transmission(const Writer& w, const analysis::PathTransmission& t) {
  std::vector<double> weight;
  std::vector<std::uint64_t> status;
  for (const analysis::PathRay& r : t.rays) {
    weight.push_back(r.weight);
    status.push_back(static_cast<std::uint64_t>(r.status));
  }
  w.f8("t_weight", weight);
  w.u8("t_status", status);
  w.f8("t_scalars", {t.mean, t.min, t.max});
  std::vector<std::uint64_t> losses{t.losses.launched};
  for (const std::size_t n : t.losses.by_status) losses.push_back(n);
  losses.push_back(t.losses.worst_surface ? *t.losses.worst_surface : kNone);
  losses.push_back(t.losses.worst_surface_count);
  w.u8("t_losses", losses);
  w.u8("t_ints", {t.rays_launched, t.rays_arrived, t.warnings.size()});
}

/// Start rays on a square grid with a pitch of 0.5 mm inside r <= 1 mm at z = 0 (the plate
/// starts at z = 10), direction (dx, dy, dz), as start_rays() in test_bitwise_crystal.py.
RayBatch start_rays(double dx, double dy, double dz) {
  std::vector<std::pair<double, double>> points;
  for (int i = -2; i <= 2; ++i) {
    for (int j = -2; j <= 2; ++j) {
      const double x = 0.5 * i;
      const double y = 0.5 * j;
      if (x * x + y * y <= 1.0) points.emplace_back(x, y);
    }
  }
  RayBatch rays(points.size());
  for (std::size_t k = 0; k < points.size(); ++k) {
    rays.pos_x()[k] = points[k].first;
    rays.pos_y()[k] = points[k].second;
    rays.dir_x()[k] = dx;
    rays.dir_y()[k] = dy;
    rays.dir_z()[k] = dz;
  }
  return rays;
}

/// The calcite plate with the paths "o entry" and "e entry", which end at P.S1 inside the
/// crystal.
model::System with_entry_paths(model::System s) {
  s.paths.push_back(
      {"o entry", false, {{model::SurfaceId("P.S1"), model::EventKind::Ordinary, 0}}});
  s.paths.push_back(
      {"e entry", false, {{model::SurfaceId("P.S1"), model::EventKind::Extraordinary, 0}}});
  return s;
}

}  // namespace

void run_crystal_cases(const fs::path& reference_dir,
                       const fs::path& /*catalog_dir*/,
                       const fs::path& out,
                       int threads) {
  const material::MaterialLibrary plain;
  const CompiledSystem plate = compile::compile(
      with_entry_paths(io::load_system(reference_dir / "m4/calcite_walkoff.rtt.json")), plain);
  // Normal incidence, and oblique with the direction (1, 1, 5) / sqrt(27) (no trigonometry, so
  // that both sides build the same bits): out of the principal plane of the optic axis.
  const double norm = std::sqrt(27.0);
  const RayBatch normal = start_rays(0.0, 0.0, 1.0);
  const RayBatch oblique = start_rays(1.0 / norm, 1.0 / norm, 5.0 / norm);

  struct Case {
    const char* name;
    const char* path;
    const RayBatch* start;
  };
  const Case cases[] = {
      {"calcite_o_normal", "o", &normal},
      {"calcite_e_normal", "e", &normal},
      {"calcite_o_oblique", "o", &oblique},
      {"calcite_e_oblique", "e", &oblique},
      {"calcite_o_entry_oblique", "o entry", &oblique},
      {"calcite_e_entry_normal", "e entry", &normal},
      {"calcite_e_entry_oblique", "e entry", &oblique},
  };
  oneapi::tbb::task_arena arena(threads);
  arena.execute([&] {
    for (const Case& c : cases) {
      const PathId path = *plate.find_path(c.path);
      RayBatch rays = *c.start;
      const trace::TraceStats stats = trace::SequentialTracer{}.trace(plate, path, rays);
      const Writer w(out, c.name);
      write_rays(w, rays, stats);
      write_transmission(w, analysis::path_transmission(plate, path, *c.start));
    }
  });
}

}  // namespace rtt::py::reference
