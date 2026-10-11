// Q and P per slot cases of rtt_py_reference (#207, ADR 0034, points 2, 3 and 8):
// tests/reference/m3/metal_periscope.rtt.json along both paths with three rays set by hand with
// the same operations as start_rays() in test_bitwise_q.py, Q requested, P and Q recorded per
// slot. Written as named .npy arrays (the per-slot arrays flattened in C order).

#include "q_cases.hpp"

#include <oneapi/tbb/task_arena.h>

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

#include "npy_writer.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/ray_paths.hpp"
#include "rtt/trace/sequential.hpp"

namespace rtt::py::reference {
namespace {

namespace fs = std::filesystem;
using trace::RayBatch;

/// Three rays along +z from (0, 0, 0), (0.5, 0, 0) and one along +x from the origin, which
/// runs parallel to the plane of the first mirror (Missed there); Q requested.
RayBatch start_rays() {
  RayBatch rays(3);
  rays.pos_x()[1] = 0.5;
  rays.dir_x()[2] = 1.0;
  rays.dir_z()[2] = 0.0;
  rays.enable_q();
  return rays;
}

void write_case(const fs::path& out,
                const std::string& name,
                const trace::TraceStats& stats,
                const RayBatch& r,
                const trace::RayPaths& paths) {
  const auto file = [&](const std::string& array) { return out / (name + "." + array + ".npy"); };
  write_column<double>(file("dir_x"), "<f8", r.dir_x());
  write_column<double>(file("dir_y"), "<f8", r.dir_y());
  write_column<double>(file("dir_z"), "<f8", r.dir_z());
  write_column<double>(file("weight"), "<f8", r.weight());
  write_column<trace::RayStatus>(file("status"), "|u1", r.status());
  for (int row = 0; row < 3; ++row) {
    for (int col = 0; col < 3; ++col) {
      const std::string rc = std::to_string(row) + std::to_string(col);
      write_column<math::Complex>(file("prt" + rc), "<c16", r.prt(row, col));
      write_column<double>(file("q" + rc), "<f8", r.q(row, col));
    }
  }
  write_column<math::Complex>(file("paths_prt"), "<c16", paths.prt);
  write_column<double>(file("paths_q"), "<f8", paths.q);
  write_column<std::uint32_t>(file("paths_count"), "<u4", paths.count);
  const std::vector<std::uint64_t> counts(stats.rays.begin(), stats.rays.end());
  write_column<std::uint64_t>(file("stats"), "<u8", counts);
}

}  // namespace

void run_q_cases(const fs::path& reference_dir, const fs::path& out, int threads) {
  const compile::CompiledSystem cs = compile::compile(
      io::load_system(reference_dir / "m3/metal_periscope.rtt.json"), material::MaterialLibrary{});
  oneapi::tbb::task_arena arena(threads);
  arena.execute([&] {
    for (const char* path : {"single", "periscope"}) {
      RayBatch rays = start_rays();
      trace::RayPaths paths;
      paths.polarization = true;
      const trace::TraceStats stats =
          trace::SequentialTracer{}.trace(cs, *cs.find_path(path), rays, paths);
      write_case(out, std::string("q_") + path, stats, rays, paths);
    }
  });
}

}  // namespace rtt::py::reference
