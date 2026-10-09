// Configuration cases of rtt_py_reference (#169, ADR 0028 and 0029): the zoom of
// tests/reference/m5/zoom.rtt.json compiled for each configuration (by index and by name, as
// Python calls it), traced with start rays set by hand, and the frames of every node. Written as
// named .npy arrays (the Python side flattens its results with the same names).

#include "configuration_cases.hpp"

#include <oneapi/tbb/task_arena.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>
#include <utility>
#include <vector>

#include "npy_writer.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/sequential.hpp"

namespace rtt::py::reference {
namespace {

namespace fs = std::filesystem;
using compile::CompiledSystem;
using trace::RayBatch;

/// Three rays from z = -5 at y = -1, 0, 1 mm along +z and one tilted ray, as start_rays() in
/// test_bitwise_configurations.py.
RayBatch start_rays() {
  RayBatch rays(4);
  for (std::size_t k = 0; k < 3; ++k) {
    rays.pos_y()[k] = static_cast<double>(k) - 1.0;
    rays.pos_z()[k] = -5.0;
  }
  rays.pos_z()[3] = -5.0;
  rays.dir_y()[3] = 0.0625;
  rays.dir_z()[3] = std::sqrt(1.0 - 0.0625 * 0.0625);  // correctly rounded on both sides
  return rays;
}

/// The global transforms of all node frames: 16 values per node (reference, then to_global),
/// row-major 3 x 4 each with the translation in the last column (12 values), in node order.
std::vector<double> frames(const CompiledSystem& cs) {
  std::vector<double> out;
  for (const compile::NodeFrame& f : cs.node_frames()) {
    for (const math::Isometry3* t : {&f.reference, &f.to_global}) {
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) out.push_back(t->rotation()(i, j));
        out.push_back(t->translation()(i));
      }
    }
  }
  return out;
}

void write_case(const fs::path& out,
                const std::string& name,
                const CompiledSystem& cs,
                const trace::TraceStats& stats,
                const RayBatch& r) {
  const auto file = [&](const std::string& array) { return out / (name + "." + array + ".npy"); };
  write_column<double>(file("pos_x"), "<f8", r.pos_x());
  write_column<double>(file("pos_y"), "<f8", r.pos_y());
  write_column<double>(file("pos_z"), "<f8", r.pos_z());
  write_column<double>(file("dir_y"), "<f8", r.dir_y());
  write_column<double>(file("dir_z"), "<f8", r.dir_z());
  write_column<double>(file("opl"), "<f8", r.opl());
  write_column<trace::RayStatus>(file("status"), "|u1", r.status());
  const std::vector<std::uint64_t> counts(stats.rays.begin(), stats.rays.end());
  write_column<std::uint64_t>(file("stats"), "<u8", counts);
  const std::vector<std::uint64_t> config{cs.configuration()};
  write_column<std::uint64_t>(file("configuration"), "<u8", config);
  write_column<double>(file("frames"), "<f8", frames(cs));
}

}  // namespace

void run_configuration_cases(const fs::path& reference_dir,
                             const fs::path& /*catalog_dir*/,
                             const fs::path& out,
                             int threads) {
  const model::System zoom = io::load_system(reference_dir / "m5/zoom.rtt.json");
  const material::MaterialLibrary plain;
  struct Case {
    const char* name;
    CompiledSystem system;
  };
  const Case cases[] = {
      {"zoom_wide", compile::compile(zoom, plain, std::size_t{0})},
      {"zoom_tele", compile::compile(zoom, plain, "tele")},
  };
  oneapi::tbb::task_arena arena(threads);
  arena.execute([&] {
    for (const Case& c : cases) {
      RayBatch rays = start_rays();
      const trace::TraceStats stats =
          trace::SequentialTracer{}.trace(c.system, compile::PathId{0}, rays);
      write_case(out, c.name, c.system, stats, rays);
    }
  });
}

}  // namespace rtt::py::reference
