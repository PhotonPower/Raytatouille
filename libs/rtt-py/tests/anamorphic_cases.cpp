// Anamorphic cases of rtt_py_reference (#315, ADR 0042, acceptance 9): the examples under
// examples/anamorphic (biconic, cylinder, crossed cylinders) with four field angles in x and y,
// aimed per section (ADR 0042, point 3) with a hexapolar pupil, real and paraxial, and traced.
// Written as named .npy arrays, as in test_bitwise_anamorphic.py.

#include "anamorphic_cases.hpp"

#include <oneapi/tbb/task_arena.h>

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <string>
#include <variant>
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

/// The case "toroid_before_stop", as toroid_before_stop() in test_bitwise_anamorphic.py: the
/// cylinder example with a toroid (profile R = 51.68, radius of rotation 77.52) as L1.S1 and the
/// stop (r = 4) moved behind the lens to z = 20, system aperture stop_size. The entrance pupil
/// then differs per section in position and size (ADR 0042, acceptance 6 and 8).
void make_toroid_before_stop(model::System& s) {
  auto& children = s.root.children;
  std::rotate(children.begin(), children.begin() + 1, children.begin() + 2);  // lens, stop, ...
  auto& lens = std::get<model::Element>(children[0].value);
  lens.surfaces[0].shape.base =
      model::Toroid{model::Param(51.68), model::Param(0.0), {}, model::Param(77.52)};
  auto& stop = std::get<model::Element>(children[1].value);
  stop.pose = model::Pose::along_z(20.0);
  stop.surfaces[0].aperture = model::CircularAperture{4.0, 0.0};
  s.aperture = {model::SystemApertureType::StopSize, model::Param(0.0)};
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

void run_anamorphic_cases(const fs::path& reference_dir, const fs::path& out, int threads) {
  const fs::path examples = reference_dir / ".." / ".." / "examples" / "anamorphic";
  // The three examples, and the cylinder example turned into a toroid in front of the stop.
  struct File {
    const char* stem;
    const char* file;
  };
  const File files[] = {{"biconic_singlet", "biconic_singlet"},
                        {"crossed_cylinders", "crossed_cylinders"},
                        {"cylinder_planoconvex", "cylinder_planoconvex"},
                        {"toroid_before_stop", "cylinder_planoconvex"}};
  oneapi::tbb::task_arena arena(threads);
  arena.execute([&] {
    for (const File& f : files) {
      const std::string stem = f.stem;
      model::System s = io::load_system(examples / (std::string(f.file) + ".rtt.json"));
      if (stem == "toroid_before_stop") make_toroid_before_stop(s);
      // Field angles in x and y (degree), as FIELDS in test_bitwise_anamorphic.py.
      s.fields.points = {{0.0, 0.0, 1.0}, {1.5, 0.0, 1.0}, {0.0, 1.5, 1.0}, {-1.0, 0.5, 1.0}};
      const compile::CompiledSystem cs = compile::compile(s, material::MaterialLibrary{});
      const std::vector<std::uint16_t> fields{0, 1, 2, 3};
      for (const trace::Aiming aiming : {trace::Aiming::Real, trace::Aiming::Paraxial}) {
        const RayBatch start =
            trace::make_rays(cs, compile::PathId{0}, fields, cs.reference_wavelength(),
                             trace::HexapolarPupil{3}, aiming);
        RayBatch rays = start;
        const trace::TraceStats stats =
            trace::SequentialTracer{}.trace(cs, compile::PathId{0}, rays);
        const std::string name =
            "anamorphic_" + stem + (aiming == trace::Aiming::Real ? "_real" : "_paraxial");
        write_case(out, name, stats, start, rays);
      }
    }
  });
}

}  // namespace rtt::py::reference
