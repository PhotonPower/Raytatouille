// Report cases of rtt_py_reference (#177): the same calls as CASES in test_bitwise_reports.py,
// every report written as named .npy arrays (the Python side flattens its reports with the same
// names).

#include "reports_cases.hpp"

#include <oneapi/tbb/task_arena.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "npy_writer.hpp"
#include "rtt/analysis/reports.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/paraxial/prescription.hpp"
#include "rtt/trace/ray_batch.hpp"

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

  template <typename T, typename Items, typename Get>
  void column(const std::string& array, const char* descr, const Items& items, Get get) const {
    std::vector<T> v;
    for (const auto& item : items) v.push_back(static_cast<T>(get(item)));
    write_column<T>(out_ / (name_ + "." + array + ".npy"), descr, v);
  }
  void f8(const std::string& array, const std::vector<double>& v) const {
    column<double>(array, "<f8", v, [](double x) { return x; });
  }
  void u8(const std::string& array, const std::vector<std::uint64_t>& v) const {
    column<std::uint64_t>(array, "<u8", v, [](std::uint64_t x) { return x; });
  }

 private:
  fs::path out_;
  std::string name_;
};

// Flattening of the reports; the names match the flatten_* functions in
// test_bitwise_reports.py. Integers (indices, counts) as uint64, None as kNone or NaN.

void write_raytrace(const Writer& w, const RaytraceReport& r) {
  const auto& rows = r.rows;
  w.column<std::uint64_t>("ray", "<u8", rows, [](const RaytraceRow& x) { return x.ray; });
  w.column<std::uint64_t>("slot", "<u8", rows, [](const RaytraceRow& x) { return x.slot; });
  w.column<std::uint32_t>("surface", "<u4", rows, [](const RaytraceRow& x) { return x.surface; });
  w.column<double>("x", "<f8", rows, [](const RaytraceRow& x) { return x.x; });
  w.column<double>("y", "<f8", rows, [](const RaytraceRow& x) { return x.y; });
  w.column<double>("z", "<f8", rows, [](const RaytraceRow& x) { return x.z; });
  w.column<double>("dx", "<f8", rows, [](const RaytraceRow& x) { return x.dx; });
  w.column<double>("dy", "<f8", rows, [](const RaytraceRow& x) { return x.dy; });
  w.column<double>("dz", "<f8", rows, [](const RaytraceRow& x) { return x.dz; });
  w.column<double>("local_x", "<f8", rows, [](const RaytraceRow& x) { return x.local_x; });
  w.column<double>("local_y", "<f8", rows, [](const RaytraceRow& x) { return x.local_y; });
  w.column<double>("local_z", "<f8", rows, [](const RaytraceRow& x) { return x.local_z; });
  w.column<double>("local_dx", "<f8", rows, [](const RaytraceRow& x) { return x.local_dx; });
  w.column<double>("local_dy", "<f8", rows, [](const RaytraceRow& x) { return x.local_dy; });
  w.column<double>("local_dz", "<f8", rows, [](const RaytraceRow& x) { return x.local_dz; });
  w.column<double>("opl", "<f8", rows, [](const RaytraceRow& x) { return x.opl; });
  w.column<double>("weight", "<f8", rows, [](const RaytraceRow& x) { return x.weight; });
  w.column<std::uint8_t>("status", "|u1", rows,
                         [](const RaytraceRow& x) { return static_cast<std::uint8_t>(x.status); });
  w.u8("ints", {r.path.index, r.rays, r.slots});
}

double value(const std::optional<double>& v) {
  return v.value_or(kNan);
}

void write_system(const Writer& w, const SystemReport& r) {
  w.u8("ints",
       {r.path.index, r.wavelength, r.reference_wavelength, r.field_count, r.surface_count,
        r.event_count, r.stop ? *r.stop : kNone, r.prescription ? 1U : 0U, r.warnings.size()});
  w.f8("wavelengths_um", r.wavelengths_um);
  if (!r.prescription) return;
  const paraxial::Prescription& p = *r.prescription;
  using paraxial::PrescriptionSurface;
  const auto ray = [](const std::optional<paraxial::PrescriptionRay>& x, int k) {
    if (!x) return kNan;
    return k == 0 ? x->y : (k == 1 ? x->u : x->i);
  };
  w.column<std::uint64_t>("p_surface", "<u8", p.surfaces,
                          [](const PrescriptionSurface& s) { return s.surface; });
  w.column<double>("p_z", "<f8", p.surfaces, [](const PrescriptionSurface& s) { return s.z; });
  w.column<double>("p_n", "<f8", p.surfaces, [](const PrescriptionSurface& s) { return s.n; });
  w.column<double>("p_y", "<f8", p.surfaces,
                   [&](const PrescriptionSurface& s) { return ray(s.marginal, 0); });
  w.column<double>("p_u", "<f8", p.surfaces,
                   [&](const PrescriptionSurface& s) { return ray(s.marginal, 1); });
  w.column<double>("p_i", "<f8", p.surfaces,
                   [&](const PrescriptionSurface& s) { return ray(s.marginal, 2); });
  w.column<double>("p_y_bar", "<f8", p.surfaces,
                   [&](const PrescriptionSurface& s) { return ray(s.chief, 0); });
  w.column<double>("p_u_bar", "<f8", p.surfaces,
                   [&](const PrescriptionSurface& s) { return ray(s.chief, 1); });
  w.column<double>("p_i_bar", "<f8", p.surfaces,
                   [&](const PrescriptionSurface& s) { return ray(s.chief, 2); });
  w.column<double>("p_lagrange", "<f8", p.surfaces,
                   [](const PrescriptionSurface& s) { return value(s.lagrange); });
  w.f8("p_scalars", {p.total_track, value(p.object_distance), value(p.paraxial_working_f_number),
                     value(p.paraxial_image_na), value(p.lagrange_invariant),
                     value(p.first_order.efl), value(p.first_order.bfl)});
}

void write_dimensions(const Writer& w, const DimensionReport& d) {
  const auto& s = d.segments;
  w.column<std::uint32_t>("element", "<u4", s,
                          [](const SegmentDimensions& x) { return x.element; });
  w.column<std::uint32_t>("first_surface", "<u4", s,
                          [](const SegmentDimensions& x) { return x.first_surface; });
  w.column<std::uint8_t>("coaxial", "|b1", s,
                         [](const SegmentDimensions& x) { return x.coaxial ? 1 : 0; });
  w.column<double>("centre_thickness", "<f8", s,
                   [](const SegmentDimensions& x) { return x.centre_thickness; });
  w.column<double>("semi_diameter_first", "<f8", s,
                   [](const SegmentDimensions& x) { return x.semi_diameter_first; });
  w.column<double>("semi_diameter_second", "<f8", s,
                   [](const SegmentDimensions& x) { return x.semi_diameter_second; });
  w.column<std::uint8_t>("aperture_first", "|u1", s, [](const SegmentDimensions& x) {
    return static_cast<std::uint8_t>(x.aperture_first);
  });
  w.column<std::uint8_t>("aperture_second", "|u1", s, [](const SegmentDimensions& x) {
    return static_cast<std::uint8_t>(x.aperture_second);
  });
  w.column<double>("edge_thickness", "<f8", s,
                   [](const SegmentDimensions& x) { return x.edge_thickness; });
  w.column<double>("diameter", "<f8", s, [](const SegmentDimensions& x) { return x.diameter; });
}

/// As start_rays() in test_bitwise_reports.py: an axial ray, a ray at y = 5 mm and a ray
/// steeply off axis (lost), all from z = -1 (in front of the stop plane z = 0).
RayBatch start_rays() {
  RayBatch rays(3);
  for (std::size_t k = 0; k < rays.size(); ++k) rays.pos_z()[k] = -1.0;
  rays.pos_y()[1] = 5.0;
  rays.dir_y()[2] = 0.8;
  rays.dir_z()[2] = 0.6;
  return rays;
}

}  // namespace

void run_reports_cases(const fs::path& reference_dir,
                       const fs::path& catalog_dir,
                       const fs::path& out,
                       int threads) {
  const auto w = [&](const std::string& name) { return Writer(out, name); };
  material::MaterialLibrary plain;
  material::MaterialLibrary schott;
  schott.add_catalog(catalog_dir / "schott.agf");
  const CompiledSystem singlet =
      compile::compile(io::load_system(reference_dir / "m1/singlet_const.rtt.json"), plain);
  const CompiledSystem cooke =
      compile::compile(io::load_system(reference_dir / "m2/cooke_triplet.rtt.json"), schott);
  const CompiledSystem grating =
      compile::compile(io::load_system(reference_dir / "m4/grating_transmission.rtt.json"), plain);

  oneapi::tbb::task_arena arena(threads);
  arena.execute([&] {
    write_raytrace(w("singlet_raytrace"), raytrace_report(singlet, PathId{0}, start_rays()));
    write_raytrace(w("cooke_raytrace"), raytrace_report(cooke, PathId{0}, start_rays()));
  });
  write_system(w("singlet_system"), system_report(singlet, PathId{0}, 1));
  write_system(w("cooke_system"), system_report(cooke, PathId{0}, cooke.reference_wavelength()));
  write_system(w("grating_system"), system_report(grating, *grating.find_path("order +1"), 0));
  write_dimensions(w("singlet_dimensions"), dimension_report(singlet));
  write_dimensions(w("cooke_dimensions"), dimension_report(cooke));
}

}  // namespace rtt::py::reference
