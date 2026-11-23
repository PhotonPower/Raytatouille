// M8 cases of rtt_py_reference (#324, plan v1, section 3): the freeform mirrors of S3
// (zernike_mirror, xy_mirror, grid_sag_mirror) and the aperture systems of F1 (hexagon_mirror,
// spider_telescope) under tests/reference/m8, with three fields (FIELDS in test_bitwise_m8.py),
// aimed with a 21 x 21 grid pupil and traced; for the aimable systems also the spot of field 1,
// first_order and seidel. The grid-sag mirror has no paraxial model (ADR 0051, point 5) and
// Python has no aiming reference: its start rays are those of its nominal design (the conic
// without the term), traced through the grid mirror. Written as named .npy arrays, as in
// test_bitwise_m8.py.

#include "m8_cases.hpp"

#include <oneapi/tbb/task_arena.h>

#include <cstdint>
#include <filesystem>
#include <limits>
#include <optional>
#include <span>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "npy_writer.hpp"
#include "rtt/analysis/spot.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
#include "rtt/paraxial/paraxial.hpp"
#include "rtt/paraxial/seidel.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/sequential.hpp"
#include "rtt/trace/sources.hpp"

namespace rtt::py::reference {
namespace {

namespace fs = std::filesystem;
using compile::CompiledSystem;
using compile::PathId;
using trace::RayBatch;

constexpr std::uint64_t kNone = std::numeric_limits<std::uint64_t>::max();
constexpr double kNan = std::numeric_limits<double>::quiet_NaN();

class Writer {
 public:
  Writer(fs::path out, std::string name) : out_(std::move(out)), name_(std::move(name)) {}
  void f8(const std::string& array, const std::vector<double>& v) const {
    write_column<double>(file(array), "<f8", v);
  }
  void f8(const std::string& array, std::span<const double> v) const {
    write_column<double>(file(array), "<f8", v);
  }
  void u8(const std::string& array, const std::vector<std::uint64_t>& v) const {
    write_column<std::uint64_t>(file(array), "<u8", v);
  }
  [[nodiscard]] fs::path file(const std::string& array) const {
    return out_ / (name_ + "." + array + ".npy");
  }

 private:
  fs::path out_;
  std::string name_;
};

/// Start rays and traced columns, as write_case() in anamorphic_cases.cpp.
void write_rays(const Writer& w,
                const trace::TraceStats& stats,
                const RayBatch& start,
                const RayBatch& r) {
  w.f8("start_pos_x", start.pos_x());
  w.f8("start_pos_y", start.pos_y());
  w.f8("start_pos_z", start.pos_z());
  w.f8("start_dir_x", start.dir_x());
  w.f8("start_dir_y", start.dir_y());
  w.f8("start_dir_z", start.dir_z());
  w.f8("pos_x", r.pos_x());
  w.f8("pos_y", r.pos_y());
  w.f8("pos_z", r.pos_z());
  w.f8("dir_x", r.dir_x());
  w.f8("dir_y", r.dir_y());
  w.f8("dir_z", r.dir_z());
  w.f8("opl", r.opl());
  write_column<std::uint16_t>(w.file("field"), "<u2", r.field());
  write_column<std::uint32_t>(w.file("last_surface"), "<u4", r.last_surface());
  write_column<trace::RayStatus>(w.file("status"), "|u1", r.status());
  w.u8("stats", std::vector<std::uint64_t>(stats.rays.begin(), stats.rays.end()));
}

/// As write_spot() in analysis_cases.cpp and flatten_spot() in test_bitwise_analysis.py.
void write_spot(const Writer& w, const analysis::SpotDiagram& s) {
  std::vector<double> x;
  std::vector<double> y;
  std::vector<std::uint64_t> wl;
  std::vector<double> weight;
  for (const analysis::SpotPoint& p : s.points) {
    x.push_back(p.x);
    y.push_back(p.y);
    wl.push_back(p.wavelength);
    weight.push_back(p.weight);
  }
  w.f8("spot_x", x);
  w.f8("spot_y", y);
  w.u8("spot_wavelengths", wl);
  w.f8("spot_weight", weight);
  w.f8("spot_scalars",
       {s.chief.x, s.chief.y, s.stats.centroid.x, s.stats.centroid.y, s.stats.rms_centroid,
        s.stats.rms_chief, s.stats.geo_centroid, s.stats.geo_chief, s.vignetted_fraction});
  w.u8("spot_ints", {s.field, s.wavelength ? *s.wavelength : kNone, s.image_surface,
                     s.rays_launched, s.rays_arrived});
}

/// As first_order_values() in rtt_py_reference.cpp and test_bitwise.py.
void write_first_order(const Writer& w, const paraxial::FirstOrder& fo) {
  const auto value = [](const std::optional<double>& v) { return v.value_or(kNan); };
  const auto pupil = [&](const std::optional<paraxial::Pupil>& p, bool diameter) {
    if (!p) return kNan;
    return value(diameter ? p->diameter : p->z);
  };
  w.f8("first_order", {fo.object_index,
                       fo.image_index,
                       static_cast<double>(fo.image_direction),
                       fo.power,
                       value(fo.efl),
                       value(fo.front_focal_length),
                       value(fo.rear_focal_length),
                       value(fo.ffl),
                       value(fo.bfl),
                       value(fo.front_focal_z),
                       value(fo.rear_focal_z),
                       value(fo.front_principal_z),
                       value(fo.rear_principal_z),
                       value(fo.image_z),
                       value(fo.lateral_magnification),
                       value(fo.angular_magnification),
                       pupil(fo.entrance_pupil, false),
                       pupil(fo.entrance_pupil, true),
                       pupil(fo.exit_pupil, false),
                       pupil(fo.exit_pupil, true)});
}

/// Sums of seidel(), as the "scalars" of write_seidel() in analysis_cases.cpp.
void write_seidel(const Writer& w, const paraxial::Seidel& s) {
  w.f8("seidel",
       {s.sum.s1, s.sum.s2, s.sum.s3, s.sum.s4, s.sum.s5, s.sum.c_l, s.sum.c_t, s.lagrange,
        s.marginal.z, s.marginal.y, s.marginal.u, s.chief.z, s.chief.y, s.chief.u});
}

model::System load(const fs::path& reference_dir, const std::string& stem) {
  model::System s = io::load_system(reference_dir / "m8" / (stem + ".rtt.json"));
  // Field angles in x and y (degree), as FIELDS in test_bitwise_m8.py.
  s.fields.points = {{0.0, 0.0, 1.0}, {0.0, 0.5, 1.0}, {0.3, 0.0, 1.0}};
  return s;
}

}  // namespace

void run_m8_cases(const fs::path& reference_dir, const fs::path& out, int threads) {
  const std::vector<std::uint16_t> fields{0, 1, 2};
  const trace::GridPupil grid{21};
  oneapi::tbb::task_arena arena(threads);
  arena.execute([&] {
    for (const char* stem : {"zernike_mirror", "xy_mirror", "hexagon_mirror", "spider_telescope"}) {
      const CompiledSystem cs =
          compile::compile(load(reference_dir, stem), material::MaterialLibrary{});
      const std::uint16_t wl = cs.reference_wavelength();
      const Writer w(out, std::string("m8_") + stem);
      const RayBatch start = trace::make_rays(cs, PathId{0}, fields, wl, grid);
      RayBatch rays = start;
      const trace::TraceStats stats = trace::SequentialTracer{}.trace(cs, PathId{0}, rays);
      write_rays(w, stats, start, rays);
      write_spot(w, analysis::spot(cs, PathId{0}, 1, std::nullopt, {grid}));
      write_first_order(w, paraxial::first_order(cs, PathId{0}, wl));
      write_seidel(w, paraxial::seidel(cs, PathId{0}, wl));
    }
    // Grid sag: start rays of the nominal design, traced through the grid mirror.
    const model::System grid_model = load(reference_dir, "grid_sag_mirror");
    model::System nominal_model = grid_model;
    std::get<model::Element>(nominal_model.root.children[1].value).surfaces[0].shape.terms.clear();
    const CompiledSystem cs = compile::compile(grid_model, material::MaterialLibrary{});
    const CompiledSystem nominal = compile::compile(nominal_model, material::MaterialLibrary{});
    const RayBatch start =
        trace::make_rays(nominal, PathId{0}, fields, nominal.reference_wavelength(), grid);
    RayBatch rays = start;
    const trace::TraceStats stats = trace::SequentialTracer{}.trace(cs, PathId{0}, rays);
    write_rays(Writer(out, "m8_grid_sag_mirror"), stats, start, rays);
  });
}

}  // namespace rtt::py::reference
