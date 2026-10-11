// Diffraction cases of rtt_py_reference (#209, P1b): the same calls as CASES in
// test_bitwise_psf.py, every result written as named .npy arrays (the Python side flattens with
// the same names).

#include "psf_cases.hpp"

#include <oneapi/tbb/task_arena.h>

#include <array>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <span>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "npy_writer.hpp"
#include "rtt/analysis/psf.hpp"
#include "rtt/analysis/pupil.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"

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

/// centre, e_x, e_y (9 values).
std::vector<double> frame(const math::Vec3& c, const math::Vec3& x, const math::Vec3& y) {
  return {c.x(), c.y(), c.z(), x.x(), x.y(), x.z(), y.x(), y.y(), y.z()};
}

/// The wavelengths as (weight, energy, outside_field) per entry and (wavelength, length).
void write_wavelengths(const Writer& w, const std::vector<PsfWavelength>& wls) {
  std::vector<double> f;
  std::vector<std::uint64_t> i;
  for (const PsfWavelength& x : wls) {
    f.insert(f.end(), {x.weight, x.energy, x.outside_field});
    i.insert(i.end(), {x.wavelength, static_cast<std::uint64_t>(x.length)});
  }
  w.f8("wavelengths", f);
  w.u8("wavelength_ints", i);
}

// Flattening; names as flatten_psf() and friends in test_bitwise_psf.py.

void write_psf(const Writer& w, const Psf& p) {
  w.f8("intensity", p.intensity);
  w.f8("frame", frame(p.centre, p.e_x, p.e_y));
  w.f8("lattice", {p.lattice.begin(), p.lattice.end()});
  w.f8("scalars",
       {p.strehl, p.edge_energy_deficit, p.polychromatic_scale_error, p.peak.value,
        static_cast<double>(p.peak.k1), static_cast<double>(p.peak.k2), p.peak.x, p.peak.y});
  w.u8("ints", {p.field, static_cast<std::uint64_t>(p.grid), static_cast<std::uint64_t>(p.length)});
  write_wavelengths(w, p.wavelengths);
  w.u8("losses", losses(p.losses, p.warnings.size()));
}

void write_through_focus(const Writer& w, const ThroughFocus& t) {
  w.f8("shifts", t.shifts);
  w.f8("strehl", t.strehl);
  w.f8("arrived_fraction", t.arrived_fraction);
  w.u8("ints", {t.field, t.wavelength});
  w.u8("losses", losses(t.losses, t.warnings.size()));
}

void write_curve(const Writer& w, const std::string& prefix, const MtfCurve& c) {
  w.f8(prefix + "_direction", {c.direction[0], c.direction[1]});
  w.f8(prefix + "_frequency", c.frequency);
  w.f8(prefix + "_modulation", c.modulation);
  w.c16(prefix + "_otf", c.otf);
}

void write_mtf(const Writer& w, const Mtf& m) {
  write_curve(w, "t", m.tangential);
  write_curve(w, "s", m.sagittal);
  w.f8("frequency_matrix", {m.frequency_matrix.begin(), m.frequency_matrix.end()});
  w.f8("map", m.map);
  w.c16("otf_map", m.otf_map);
  w.f8("frame", frame(m.centre, m.e_x, m.e_y));
  w.f8("scalars", {m.edge_energy_deficit});
  w.u8("ints", {m.field, static_cast<std::uint64_t>(m.grid)});
  write_wavelengths(w, m.wavelengths);
  w.u8("losses", losses(m.losses, m.warnings.size()));
}

void write_ee(const Writer& w, const EncircledEnergy& e) {
  w.f8("radii", e.radii);
  w.f8("energy", e.energy);
  w.f8("frame", frame(e.centre, e.e_x, e.e_y));
  w.f8("scalars", {e.offset.x, e.offset.y, e.radius_max, e.edge_energy_deficit});
  w.u8("ints", {e.field, static_cast<std::uint64_t>(e.grid)});
  write_wavelengths(w, e.wavelengths);
  w.u8("losses", losses(e.losses, e.warnings.size()));
}

/// m2/paraboloid_stop (f = 100 mm, AIR) with the entrance pupil diameter `epd` and the fields
/// 0 and +2 deg, as systems() in test_bitwise_psf.py; `two_wavelengths` replaces the wavelengths
/// by 0.5876 um (reference, weight 1) and 0.7345 um (weight 3).
CompiledSystem paraboloid(const fs::path& reference_dir,
                          double epd,
                          bool two_wavelengths,
                          bool obscured = false) {
  model::System s = io::load_system(reference_dir / "m2/paraboloid_stop.rtt.json");
  s.aperture.value = model::Param(epd);
  s.fields.points.push_back({0.0, 2.0, 1.0});
  if (two_wavelengths) s.wavelengths = {{0.5876, 1.0, true}, {0.7345, 3.0, false}};
  if (obscured) {  // stop aperture 5 mm with inner radius 2 mm: epsilon = 0.4 (#229 C)
    std::get<model::Element>(s.root.children[0].value).surfaces[0].aperture =
        model::CircularAperture{5.0, 2.0};
  }
  return compile::compile(s, material::MaterialLibrary{});
}

}  // namespace

void run_psf_cases(const fs::path& reference_dir, const fs::path& out, int threads) {
  const auto w = [&](const std::string& name) { return Writer(out, name); };
  const CompiledSystem f10 = paraboloid(reference_dir, 10.0, false);
  const CompiledSystem f50 = paraboloid(reference_dir, 2.0, false);
  const CompiledSystem poly = paraboloid(reference_dir, 10.0, true);
  const CompiledSystem ring = paraboloid(reference_dir, 10.0, false, true);
  oneapi::tbb::task_arena arena(threads);
  arena.execute([&] {
    PsfOptions po;
    po.pupil.grid = 32;
    write_psf(w("f10_psf"), psf(f10, PathId{0}, 1, po));
    PsfOptions pp;
    pp.pupil.grid = 16;
    write_psf(w("poly_psf"), psf(poly, PathId{0}, 0, pp));
    PupilFunctionOptions g32;
    g32.grid = 32;
    const std::array<double, 3> shifts{0.0, 0.5, 1.0};
    write_through_focus(w("f50_through_focus"), through_focus(f50, PathId{0}, 0, 0, shifts, g32));
    MtfOptions mo;
    mo.pupil.grid = 16;
    mo.outputs.map = true;
    mo.outputs.complex_otf = true;
    write_mtf(w("f10_mtf"), mtf(f10, PathId{0}, 1, mo));
    MtfOptions mp;
    mp.pupil.grid = 16;
    write_mtf(w("poly_mtf"), mtf(poly, PathId{0}, 0, mp));
    const std::array<double, 3> radii{0.002, 0.004, 0.008};
    EncircledEnergyOptions eo;
    eo.pupil.grid = 32;
    write_ee(w("f10_ee_reference"), encircled_energy(f10, PathId{0}, 1, radii, eo));
    eo.centre = EeCentre::LinearPhase;
    write_ee(w("f10_ee_linear_phase"), encircled_energy(f10, PathId{0}, 1, radii, eo));
    const PupilFunction pupil = pupil_function(f10, PathId{0}, 1, 0, g32);
    const PupilFunction moved = refocus(f10, pupil, 0.05);
    MtfOutputs both;
    both.complex_otf = true;
    write_mtf(w("refocus_mtf"), mtf(moved, both));
    write_ee(w("refocus_ee"), encircled_energy(moved, radii, EeCentre::LinearPhase));
    const std::vector<Point2> points{{0.0, 0.0}, {0.001, 0.0}, {0.0, 0.002}, {-0.003, 0.001}};
    w("psf_at").f8("values", psf_at(pupil, points));
    MtfOptions mr;
    mr.pupil.grid = 32;
    write_mtf(w("annular_mtf"), mtf(ring, PathId{0}, 0, mr));
  });
}

}  // namespace rtt::py::reference
