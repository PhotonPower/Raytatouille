// M6 acceptance cases of rtt_py_reference (#212): the same calls as CASES in test_bitwise_m6.py
// on tests/reference/m6/ and r2/ideal_lens (the systems of test_m6_acceptance.cpp), every result
// written as named .npy arrays (the Python side flattens with the same names).

#include "m6_cases.hpp"

#include <oneapi/tbb/task_arena.h>

#include <array>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "npy_writer.hpp"
#include "rtt/analysis/jones.hpp"
#include "rtt/analysis/psf.hpp"
#include "rtt/analysis/pupil.hpp"
#include "rtt/analysis/zernike.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/trace/ray_batch.hpp"

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

// Flattening; names as the flatten_* functions in test_bitwise_m6.py.

void write_psf(const Writer& w, const Psf& p) {
  w.f8("intensity", p.intensity);
  w.f8("scalars",
       {p.strehl, p.edge_energy_deficit, p.lattice[0], p.lattice[1], p.lattice[2], p.lattice[3]});
}

void write_mtf(const Writer& w, const Mtf& m) {
  w.f8("t_frequency", m.tangential.frequency);
  w.f8("t_modulation", m.tangential.modulation);
  w.f8("s_modulation", m.sagittal.modulation);
}

void write_jones(const Writer& w, const JonesPupil& j) {
  w.c16("j", j.j);
  w.f8("w", j.w);
  w.f8("diattenuation", j.diattenuation);
  w.f8("retardance", j.retardance);
}

CompiledSystem load(const fs::path& reference_dir, const char* file) {
  return compile::compile(io::load_system(reference_dir / file), material::MaterialLibrary{});
}

/// 5 x 5 parallel rays along +z from z = 0 at x, y in {-3, -1.5, 0, 1.5, 3}, pupil labels
/// (x, y) / 3, chief index 12 (as parallel_rays() in test_m6_acceptance.cpp).
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

/// Collimated rays of pitch 1.25 mm inside r = 5 from z = -1 (as converging_rays() in
/// test_m6_acceptance.cpp); `chief` is set to the ray at (0, 0).
trace::RayBatch converging_rays(std::size_t& chief) {
  trace::RayBatch start(0);
  for (int iy = -4; iy <= 4; ++iy) {
    for (int ix = -4; ix <= 4; ++ix) {
      const double x = 1.25 * ix;
      const double y = 1.25 * iy;
      if (x * x + y * y > 25.0) continue;
      const std::size_t i = start.size();
      start.resize(i + 1);
      start.pos_x()[i] = x;
      start.pos_y()[i] = y;
      start.pos_z()[i] = -1.0;
      start.pupil_x()[i] = x / 5.0;
      start.pupil_y()[i] = y / 5.0;
      if (ix == 0 && iy == 0) chief = i;
    }
  }
  return start;
}

}  // namespace

void run_m6_cases(const fs::path& reference_dir, const fs::path& out, int threads) {
  const auto w = [&](const std::string& name) { return Writer(out, name); };
  const CompiledSystem f10 = load(reference_dir, "m6/paraboloid_f10.rtt.json");
  const CompiledSystem f50 = load(reference_dir, "m6/paraboloid_f50.rtt.json");
  const CompiledSystem plano = load(reference_dir, "m6/planoconvex_f10.rtt.json");
  const CompiledSystem metal = load(reference_dir, "m6/fold45_metal.rtt.json");
  const CompiledSystem focus = load(reference_dir, "m6/fold45_focus.rtt.json");
  const CompiledSystem ideal = load(reference_dir, "r2/ideal_lens.rtt.json");
  oneapi::tbb::task_arena arena(threads);
  arena.execute([&] {
    PsfOptions po;
    po.pupil.grid = 32;
    write_psf(w("m6_f10_psf"), psf(f10, PathId{0}, 0, po));
    MtfOptions mo;
    mo.pupil.grid = 32;
    write_mtf(w("m6_f10_mtf"), mtf(f10, PathId{0}, 0, mo));
    EncircledEnergyOptions eo;
    eo.pupil.grid = 32;
    const std::array<double, 3> radii{0.002, 0.004, 0.008};
    w("m6_f10_ee").f8("energy", encircled_energy(f10, PathId{0}, 0, radii, eo).energy);
    PupilFunctionOptions g32;
    g32.grid = 32;
    const std::array<double, 3> shifts{0.0, 0.5, 1.0};
    w("m6_f50_through_focus").f8("strehl", through_focus(f50, PathId{0}, 0, 0, shifts, g32).strehl);
    const ZernikeFit fit = zernike_fit(pupil_function(plano, PathId{0}, 0, 0, g32));
    w("m6_planoconvex_zernike").f8("coefficients", fit.coefficients);
    write_jones(w("m6_fold_metal_jones"), jones_pupil(metal, PathId{0}, parallel_rays(), 12));
    std::size_t chief = 0;
    const trace::RayBatch start = converging_rays(chief);
    write_jones(w("m6_fold_focus_jones"), jones_pupil(focus, PathId{0}, start, chief));
    write_psf(w("m6_ideal_lens_psf"), psf(ideal, PathId{0}, 0, po));
  });
}

}  // namespace rtt::py::reference
