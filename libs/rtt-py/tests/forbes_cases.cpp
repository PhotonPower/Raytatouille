// Q-con and Q-bfs cases of rtt_py_reference (#314, ADR 0043, acceptance 12): qcon_singlet and
// qbfs_singlet under examples/rotational, aimed with a hexapolar pupil, real and paraxial, at three
// fields and traced. Written as named .npy arrays, as in test_bitwise_forbes.py.

#include "forbes_cases.hpp"

#include <oneapi/tbb/task_arena.h>

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

void run_forbes_cases(const fs::path& reference_dir, const fs::path& out, int threads) {
  const fs::path examples = reference_dir / ".." / ".." / "examples" / "rotational";
  oneapi::tbb::task_arena arena(threads);
  arena.execute([&] {
    for (const std::string stem : {"qcon_singlet", "qbfs_singlet"}) {
      model::System s = io::load_system(examples / (stem + ".rtt.json"));
      // Field angles in x and y (degree), as FIELDS in test_bitwise_forbes.py.
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
        write_case(out, "forbes_" + stem + (aiming == trace::Aiming::Real ? "_real" : "_paraxial"),
                   stats, start, rays);
      }
    }
  });
}

}  // namespace rtt::py::reference
