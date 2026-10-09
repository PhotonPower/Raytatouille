// Diffraction cases of rtt_py_reference (#134, ADR 0025): the same calls as CASES in
// test_bitwise_diffraction.py. Every case traces start rays set by hand through the grating of
// tests/reference/m4/grating_transmission.rtt.json (and two variants: a reflection grating and
// diffraction efficiencies) and writes the traced columns and the path transmission as named
// .npy arrays (the Python side flattens its results with the same names).

#include "diffraction_cases.hpp"

#include <oneapi/tbb/task_arena.h>

#include <cmath>
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

/// The columns of a traced batch, as flatten_rays() in test_bitwise_diffraction.py.
void write_rays(const Writer& w, const RayBatch& r, const trace::TraceStats& stats) {
  w.column<double>("pos_x", "<f8", r.pos_x());
  w.column<double>("pos_y", "<f8", r.pos_y());
  w.column<double>("pos_z", "<f8", r.pos_z());
  w.column<double>("dir_x", "<f8", r.dir_x());
  w.column<double>("dir_y", "<f8", r.dir_y());
  w.column<double>("dir_z", "<f8", r.dir_z());
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

/// path_transmission of the same start rays, as flatten_transmission() in
/// test_bitwise_diffraction.py.
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

/// Start rays on a square grid with a pitch of 0.5 mm inside r <= 2 mm at z = -5 mm (in front of
/// the stop), direction (dx, 0, dz), as start_rays() in test_bitwise_diffraction.py.
RayBatch start_rays(double dx, double dz) {
  std::vector<std::pair<double, double>> points;
  for (int i = -4; i <= 4; ++i) {
    for (int j = -4; j <= 4; ++j) {
      const double x = 0.5 * i;
      const double y = 0.5 * j;
      if (x * x + y * y <= 4.0) points.emplace_back(x, y);
    }
  }
  RayBatch rays(points.size());
  for (std::size_t k = 0; k < points.size(); ++k) {
    rays.pos_x()[k] = points[k].first;
    rays.pos_y()[k] = points[k].second;
    rays.pos_z()[k] = -5.0;
    rays.dir_x()[k] = dx;
    rays.dir_z()[k] = dz;
  }
  return rays;
}

model::Element& grating(model::System& s) {
  return std::get<model::Element>(s.root.children[1].value);
}

/// The grating bench with an ideal mirror as the grating and the path "reflect +1": the order
/// +1 reflected back through the stop (a reflection grating, ADR 0025, point 4).
model::System reflection_grating(model::System s) {
  grating(s).surfaces[0].interaction = model::IdealMirror{};
  s.paths.push_back({"reflect +1",
                     false,
                     {{model::SurfaceId("STO"), model::EventKind::Transmit, 0},
                      {model::SurfaceId("G"), model::EventKind::Reflect, 1},
                      {model::SurfaceId("STO"), model::EventKind::Transmit, 0}}});
  return s;
}

/// The grating bench with efficiencies: order 0 0.3, order +1 0.4, all others 0 (ADR 0025,
/// point 5).
model::System with_efficiencies(model::System s) {
  grating(s).surfaces[0].diffraction_efficiency =
      std::vector<model::DiffractionEfficiency>{{0, 0.3}, {1, 0.4}};
  return s;
}

}  // namespace

void run_diffraction_cases(const fs::path& reference_dir,
                           const fs::path& /*catalog_dir*/,
                           const fs::path& out,
                           int threads) {
  const material::MaterialLibrary plain;
  const model::System bench = io::load_system(reference_dir / "m4/grating_transmission.rtt.json");
  const CompiledSystem transmission = compile::compile(bench, plain);
  const CompiledSystem reflection = compile::compile(reflection_grating(bench), plain);
  const CompiledSystem efficiency = compile::compile(with_efficiencies(bench), plain);
  // Normal incidence, and oblique with the direction (-1, 0, 5) / sqrt(26) (no trigonometry,
  // so that both sides build the same bits): order +6 is evanescent at normal incidence and
  // propagates obliquely, off the screen.
  const double norm = std::sqrt(26.0);
  const RayBatch normal = start_rays(0.0, 1.0);
  const RayBatch oblique = start_rays(-1.0 / norm, 5.0 / norm);

  struct Case {
    const char* name;
    const CompiledSystem* system;
    const char* path;
    const RayBatch* start;
  };
  const Case cases[] = {
      {"grating_m1_normal", &transmission, "order -1", &normal},
      {"grating_0_normal", &transmission, "order 0", &normal},
      {"grating_p1_normal", &transmission, "order +1", &normal},
      {"grating_p6_normal", &transmission, "order +6", &normal},
      {"grating_p1_oblique", &transmission, "order +1", &oblique},
      {"grating_p6_oblique", &transmission, "order +6", &oblique},
      {"reflection_p1", &reflection, "reflect +1", &normal},
      {"efficiency_m1", &efficiency, "order -1", &normal},
      {"efficiency_0", &efficiency, "order 0", &normal},
      {"efficiency_p1", &efficiency, "order +1", &normal},
  };
  oneapi::tbb::task_arena arena(threads);
  arena.execute([&] {
    for (const Case& c : cases) {
      const PathId path = *c.system->find_path(c.path);
      RayBatch rays = *c.start;
      const trace::TraceStats stats = trace::SequentialTracer{}.trace(*c.system, path, rays);
      const Writer w(out, c.name);
      write_rays(w, rays, stats);
      write_transmission(w, analysis::path_transmission(*c.system, path, *c.start));
    }
  });
}

}  // namespace rtt::py::reference
