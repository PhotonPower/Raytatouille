// Odd asphere and axicon cases of rtt_py_reference (#314, ADR 0043, acceptance 12): the
// examples under examples/rotational. The odd asphere singlet (A2, A3: first order, no Seidel
// sums) is aimed with a hexapolar pupil, real and paraxial, at three fields and traced; the
// axicon has no paraxial model at its apex (ADR 0043, F6), so its own start rays are traced
// (axicon_rays(), as in test_bitwise_rotational.py). Written as named .npy arrays.

#include "rotational_cases.hpp"

#include <oneapi/tbb/task_arena.h>

#include <cmath>
#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

#include "npy_writer.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/sequential.hpp"
#include "rtt/trace/sources.hpp"

namespace rtt::py::reference {
namespace {

namespace fs = std::filesystem;
using trace::RayBatch;

/// As axicon_rays() in test_bitwise_rotational.py: 41 rays from z = -1 at heights y = 0.1 j mm,
/// x = 0 for j <= 20 and x = 0.03 (j - 20) mm above; every fourth ray (j % 4 == 3) tilted along
/// (0.01, -0.02, sqrt(1 - 0.0005)). Ray 0 runs along the axis through the apex. sqrt is
/// correctly rounded in IEEE 754 on both sides.
RayBatch axicon_rays() {
  RayBatch rays(41);
  for (std::size_t j = 0; j < rays.size(); ++j) {
    const double i = static_cast<double>(j);
    rays.pos_x()[j] = j > 20 ? 0.03 * (i - 20.0) : 0.0;
    rays.pos_y()[j] = 0.1 * i;
    rays.pos_z()[j] = -1.0;
    if (j % 4 == 3) {
      rays.dir_x()[j] = 0.01;
      rays.dir_y()[j] = -0.02;
      rays.dir_z()[j] = std::sqrt(1.0 - 0.0005);
    } else {
      rays.dir_x()[j] = 0.0;
      rays.dir_y()[j] = 0.0;
      rays.dir_z()[j] = 1.0;
    }
  }
  return rays;
}

void write_case(const fs::path& out,
                const std::string& name,
                const trace::TraceStats& stats,
                const RayBatch& start,
                const RayBatch& r) {
  const auto file = [&](const std::string& array) { return out / (name + "." + array + ".npy"); };
  write_column<double>(file("start_pos_x"), "<f8", start.pos_x());
  write_column<double>(file("start_pos_y"), "<f8", start.pos_y());
  write_column<double>(file("start_pos_z"), "<f8", start.pos_z());
  write_column<double>(file("start_dir_x"), "<f8", start.dir_x());
  write_column<double>(file("start_dir_y"), "<f8", start.dir_y());
  write_column<double>(file("start_dir_z"), "<f8", start.dir_z());
  write_column<double>(file("pos_x"), "<f8", r.pos_x());
  write_column<double>(file("pos_y"), "<f8", r.pos_y());
  write_column<double>(file("pos_z"), "<f8", r.pos_z());
  write_column<double>(file("dir_x"), "<f8", r.dir_x());
  write_column<double>(file("dir_y"), "<f8", r.dir_y());
  write_column<double>(file("dir_z"), "<f8", r.dir_z());
  write_column<double>(file("opl"), "<f8", r.opl());
  write_column<std::uint16_t>(file("field"), "<u2", r.field());
  write_column<std::uint32_t>(file("last_surface"), "<u4", r.last_surface());
  write_column<trace::RayStatus>(file("status"), "|u1", r.status());
  const std::vector<std::uint64_t> counts(stats.rays.begin(), stats.rays.end());
  write_column<std::uint64_t>(file("stats"), "<u8", counts);
}

}  // namespace

void run_rotational_cases(const fs::path& reference_dir, const fs::path& out, int threads) {
  const fs::path examples = reference_dir / ".." / ".." / "examples" / "rotational";
  oneapi::tbb::task_arena arena(threads);
  arena.execute([&] {
    {
      model::System s = io::load_system(examples / "odd_asphere_singlet.rtt.json");
      // Field angles in x and y (degree), as FIELDS in test_bitwise_rotational.py.
      s.fields.points = {{0.0, 0.0, 1.0}, {0.0, 2.0, 1.0}, {1.5, -1.0, 1.0}};
      const compile::CompiledSystem cs = compile::compile(s, material::MaterialLibrary{});
      const std::vector<std::uint16_t> fields{0, 1, 2};
      for (const trace::Aiming aiming : {trace::Aiming::Real, trace::Aiming::Paraxial}) {
        const RayBatch start =
            trace::make_rays(cs, compile::PathId{0}, fields, cs.reference_wavelength(),
                             trace::HexapolarPupil{3}, aiming);
        RayBatch rays = start;
        const trace::TraceStats stats =
            trace::SequentialTracer{}.trace(cs, compile::PathId{0}, rays);
        write_case(out,
                   std::string("rotational_odd_asphere_") +
                       (aiming == trace::Aiming::Real ? "real" : "paraxial"),
                   stats, start, rays);
      }
    }
    {
      const model::System s = io::load_system(examples / "axicon_planoconvex.rtt.json");
      const compile::CompiledSystem cs = compile::compile(s, material::MaterialLibrary{});
      const RayBatch start = axicon_rays();
      RayBatch rays = start;
      const trace::TraceStats stats = trace::SequentialTracer{}.trace(cs, compile::PathId{0}, rays);
      write_case(out, "rotational_axicon", stats, start, rays);
    }
  });
}

}  // namespace rtt::py::reference
