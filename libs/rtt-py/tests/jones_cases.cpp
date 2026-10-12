// Jones pupil cases of rtt_py_reference (#208 B2): the same calls as CASES in
// test_bitwise_jones.py, every result written as named .npy arrays (the Python side flattens
// with the same names).

#include "jones_cases.hpp"

#include <oneapi/tbb/task_arena.h>

#include <complex>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "npy_writer.hpp"
#include "rtt/analysis/jones.hpp"
#include "rtt/analysis/stokes.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/sources.hpp"

namespace rtt::py::reference {
namespace {

namespace fs = std::filesystem;
using namespace rtt::analysis;
using compile::CompiledSystem;
using compile::PathId;

/// Writes the arrays of one case as <out>/<case>.<name>.npy.
class Writer {
 public:
  Writer(fs::path out, std::string name) : out_(std::move(out)), name_(std::move(name)) {}
  void f8(const std::string& array, const std::vector<double>& v) const {
    write_column<double>(file(array), "<f8", std::span<const double>(v));
  }
  void u1(const std::string& array, const std::vector<std::uint8_t>& v) const {
    write_column<std::uint8_t>(file(array), "|u1", std::span<const std::uint8_t>(v));
  }
  void u8(const std::string& array, const std::vector<std::uint64_t>& v) const {
    write_column<std::uint64_t>(file(array), "<u8", std::span<const std::uint64_t>(v));
  }
  void c16(const std::string& array, const std::vector<std::complex<double>>& v) const {
    write_column<std::complex<double>>(file(array), "<c16",
                                       std::span<const std::complex<double>>(v));
  }

 private:
  [[nodiscard]] fs::path file(const std::string& array) const {
    return out_ / (name_ + "." + array + ".npy");
  }
  fs::path out_;
  std::string name_;
};

std::vector<std::uint64_t> losses(const RayLosses& l, std::size_t warnings) {
  std::vector<std::uint64_t> v{l.launched};
  for (const std::size_t n : l.by_status) v.push_back(n);
  v.push_back(l.worst_surface ? *l.worst_surface : std::numeric_limits<std::uint64_t>::max());
  v.push_back(l.worst_surface_count);
  v.push_back(warnings);
  return v;
}

// Flattening; names and order as flatten_jones() in test_bitwise_jones.py.
void write_jones(const Writer& w, const JonesPupil& p) {
  w.c16("j", p.j);
  w.f8("px", p.px);
  w.f8("py", p.py);
  w.f8("w", p.w);
  w.f8("transmission", p.transmission);
  w.f8("diattenuation", p.diattenuation);
  w.f8("diattenuation_orientation", p.diattenuation_orientation);
  w.f8("diattenuation_ellipticity", p.diattenuation_ellipticity);
  w.f8("retardance", p.retardance);
  w.f8("retardance_orientation", p.retardance_orientation);
  w.f8("retardance_ellipticity", p.retardance_ellipticity);
  std::vector<std::uint8_t> status;
  status.reserve(p.status.size());
  for (const trace::RayStatus s : p.status) status.push_back(static_cast<std::uint8_t>(s));
  w.u1("status", status);
  std::vector<double> frames;
  for (const math::Vec3* v : {&p.x_in, &p.y_in, &p.k_in, &p.x_out, &p.y_out, &p.k_out}) {
    frames.insert(frames.end(), {v->x(), v->y(), v->z()});
  }
  w.f8("frames", frames);
  w.f8("sphere", {p.sphere.centre.x(), p.sphere.centre.y(), p.sphere.centre.z(), p.sphere.radius});
  w.f8("scalars", {p.wavelength_um, p.reference_wavelength_um, p.global_phase});
  w.u8("ints", {p.path.index, p.field, p.wavelength, static_cast<std::uint64_t>(p.basis),
                p.remove_global_phase ? 1U : 0U, p.field_amplitude ? 1U : 0U});
  w.u8("losses", losses(p.losses, p.warnings.size()));
}

// Names and order as flatten_stokes() in test_bitwise_jones.py.
void write_stokes(const Writer& w, const StokesPupil& p) {
  w.c16("e_in", {p.e1, p.e2});
  for (const auto& [name, column] :
       {std::pair{"s0", &p.s0}, std::pair{"s1", &p.s1}, std::pair{"s2", &p.s2},
        std::pair{"s3", &p.s3}, std::pair{"orientation", &p.orientation},
        std::pair{"ellipticity", &p.ellipticity}}) {
    w.f8(name, *column);
  }
  w.u8("ints", {p.path.index, p.field, p.wavelength, p.warnings.size()});
}

/// 5 x 5 rays along +z from z = 0 at x, y in {-3, -1.5, 0, 1.5, 3}, labelled (x, y) / 3; the
/// chief ray (0, 0) has index 12 (as parallel_rays() in test_bitwise_jones.py).
trace::RayBatch parallel_rays() {
  trace::RayBatch rays(25);
  for (std::size_t i = 0; i < 25; ++i) {
    const double x = -3.0 + 1.5 * static_cast<double>(i % 5);
    const double y = -3.0 + 1.5 * static_cast<double>(i / 5);
    rays.pos_x()[i] = x;
    rays.pos_y()[i] = y;
    rays.pupil_x()[i] = x / 3.0;
    rays.pupil_y()[i] = y / 3.0;
  }
  return rays;
}

}  // namespace

void run_jones_cases(const fs::path& reference_dir, const fs::path& out, int threads) {
  const auto w = [&](const std::string& name) { return Writer(out, name); };
  const material::MaterialLibrary plain;
  const CompiledSystem singlet =
      compile::compile(io::load_system(reference_dir / "m1/singlet_const.rtt.json"), plain);
  const CompiledSystem periscope =
      compile::compile(io::load_system(reference_dir / "m3/metal_periscope.rtt.json"), plain);
  oneapi::tbb::task_arena arena(threads);
  arena.execute([&] {
    const std::uint16_t ref = singlet.reference_wavelength();
    JonesPupilOptions grid;
    grid.sampling = trace::GridPupil{9};
    write_jones(w("singlet_jones"), jones_pupil(singlet, PathId{0}, 1, ref, grid));
    JonesPupilOptions dipole = grid;
    dipole.basis = JonesBasis::Dipole;
    dipole.remove_global_phase = true;
    write_jones(w("singlet_jones_dipole"), jones_pupil(singlet, PathId{0}, 1, 0, dipole));
    JonesPupilOptions field = grid;
    field.field_amplitude = true;
    write_jones(w("singlet_jones_field"), jones_pupil(singlet, PathId{0}, 1, ref, field));
    JonesPupilOptions given;
    given.remove_global_phase = true;
    write_jones(w("periscope_jones"),
                jones_pupil(periscope, PathId{1}, parallel_rays(), 12, given));
    write_stokes(w("periscope_stokes"),
                 stokes_pupil(jones_pupil(periscope, PathId{1}, parallel_rays(), 12, {}),
                              std::complex<double>(1.0, 0.0), std::complex<double>(0.0, 1.0)));
  });
}

}  // namespace rtt::py::reference
