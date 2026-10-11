// Ideal lens cases of rtt_py_reference (#178, ADR 0031): tests/reference/r2/ideal_lens.rtt.json
// traced along both paths with start rays from an object point, set by hand with the same
// operations as start_rays() in test_bitwise_ideal_lens.py. Written as named .npy arrays.

#include "ideal_lens_cases.hpp"

#include <oneapi/tbb/task_arena.h>

#include <cmath>
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
#include "rtt/trace/sequential.hpp"

namespace rtt::py::reference {
namespace {

namespace fs = std::filesystem;
using trace::RayBatch;

/// Rays from the object point (x0, y0, -75) towards the points (x, y, -1) of the stop with
/// x, y in {-6, -3, 0, 3, 6} (inside its radius 10); directions normalised with
/// sqrt(dx^2 + dy^2 + dz^2) in this order.
RayBatch start_rays(double x0, double y0) {
  RayBatch rays(25);
  std::size_t i = 0;
  for (int j = -2; j <= 2; ++j) {
    for (int k = -2; k <= 2; ++k, ++i) {
      const double dx = 3.0 * k - x0;
      const double dy = 3.0 * j - y0;
      const double dz = 74.0;
      const double n = std::sqrt(dx * dx + dy * dy + dz * dz);
      rays.pos_x()[i] = x0;
      rays.pos_y()[i] = y0;
      rays.pos_z()[i] = -75.0;
      rays.dir_x()[i] = dx / n;
      rays.dir_y()[i] = dy / n;
      rays.dir_z()[i] = dz / n;
    }
  }
  return rays;
}

void write_case(const fs::path& out,
                const std::string& name,
                const trace::TraceStats& stats,
                const RayBatch& r) {
  const auto file = [&](const std::string& array) { return out / (name + "." + array + ".npy"); };
  write_column<double>(file("pos_x"), "<f8", r.pos_x());
  write_column<double>(file("pos_y"), "<f8", r.pos_y());
  write_column<double>(file("pos_z"), "<f8", r.pos_z());
  write_column<double>(file("dir_x"), "<f8", r.dir_x());
  write_column<double>(file("dir_y"), "<f8", r.dir_y());
  write_column<double>(file("dir_z"), "<f8", r.dir_z());
  write_column<double>(file("opl"), "<f8", r.opl());
  write_column<double>(file("weight"), "<f8", r.weight());
  write_column<trace::RayStatus>(file("status"), "|u1", r.status());
  for (int row = 0; row < 3; ++row) {
    for (int col = 0; col < 3; ++col) {
      write_column<math::Complex>(file("prt" + std::to_string(row) + std::to_string(col)), "<c16",
                                  r.prt(row, col));
    }
  }
  const std::vector<std::uint64_t> counts(stats.rays.begin(), stats.rays.end());
  write_column<std::uint64_t>(file("stats"), "<u8", counts);
}

}  // namespace

void run_ideal_lens_cases(const fs::path& reference_dir, const fs::path& out, int threads) {
  const compile::CompiledSystem cs = compile::compile(
      io::load_system(reference_dir / "r2/ideal_lens.rtt.json"), material::MaterialLibrary{});
  struct Case {
    const char* name;
    const char* path;
    double x0;
    double y0;
  };
  const Case cases[] = {
      {"ideal_lens_image", "ideal lens", 3.0, -4.0},
      {"cylinder_lens_image", "cylinder lens", 3.0, 2.0},
  };
  oneapi::tbb::task_arena arena(threads);
  arena.execute([&] {
    for (const Case& c : cases) {
      RayBatch rays = start_rays(c.x0, c.y0);
      const trace::TraceStats stats =
          trace::SequentialTracer{}.trace(cs, *cs.find_path(c.path), rays);
      write_case(out, c.name, stats, rays);
    }
  });
}

}  // namespace rtt::py::reference
