// Pupil cases of rtt_py_reference (#209, P1a): the same calls as CASES in test_bitwise_pupil.py,
// every result written as named .npy arrays (the Python side flattens with the same names).

#include "pupil_cases.hpp"

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
#include "rtt/analysis/opd.hpp"
#include "rtt/analysis/pupil.hpp"
#include "rtt/analysis/zernike.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/trace/sources.hpp"

namespace rtt::py::reference {
namespace {

namespace fs = std::filesystem;
using namespace rtt::analysis;
using compile::CompiledSystem;
using compile::PathId;

constexpr double kNan = std::numeric_limits<double>::quiet_NaN();

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

std::vector<double> sphere(const ReferenceSphere& s) {
  return {s.centre.x(), s.centre.y(), s.centre.z(), s.radius};
}

std::vector<std::uint64_t> losses(const RayLosses& l, std::size_t warnings) {
  std::vector<std::uint64_t> v{l.launched};
  for (const std::size_t n : l.by_status) v.push_back(n);
  v.push_back(l.worst_surface ? *l.worst_surface : std::numeric_limits<std::uint64_t>::max());
  v.push_back(l.worst_surface_count);
  v.push_back(warnings);
  return v;
}

// Flattening; names as flatten_pupil() and friends in test_bitwise_pupil.py.

void write_pupil(const Writer& w, const PupilFunction& p) {
  w.f8("fraction", p.fraction);
  w.f8("px", p.px);
  w.f8("py", p.py);
  w.f8("amplitude", p.amplitude);
  w.f8("w", p.w);
  std::vector<std::uint8_t> status;
  status.reserve(p.status.size());
  for (const trace::RayStatus s : p.status) status.push_back(static_cast<std::uint8_t>(s));
  w.u1("status", status);
  w.c16("prt", p.prt);
  w.f8("sphere", sphere(p.sphere));
  std::vector<double> mapping(18, kNan);  // m (4), s0 (2), centre, e_x, e_y, d_c (12), deltas
  if (p.mapping) {
    const PupilMapping& m = *p.mapping;
    mapping = {m.m[0], m.m[1], m.m[2], m.m[3], m.s0[0], m.s0[1]};
    for (const math::Vec3* v : {&m.centre, &m.e_x, &m.e_y, &m.d_c}) {
      mapping.insert(mapping.end(), {v->x(), v->y(), v->z()});
    }
    mapping.insert(mapping.end(), {m.delta_max, m.delta_rms});
  }
  w.f8("mapping", mapping);
  w.f8("scalars", {p.wavelength_um, p.reference_wavelength_um, p.focus_shift, p.apodization,
                   p.arrived_fraction});
  w.u8("losses", losses(p.losses, p.warnings.size()));
}

void write_zernike(const Writer& w, const ZernikeFit& f) {
  w.f8("coefficients", f.coefficients);
  w.f8("scalars", {f.residual_rms, f.sigma, f.sigma_bf, f.strehl, f.arrived_fraction});
  w.u8("ints", {f.cells, f.warnings.size()});
}

void write_points(const Writer& w, const OpdPupilPoints& p) {
  std::vector<double> px;
  std::vector<double> py;
  std::vector<double> wv;
  std::vector<std::uint8_t> status;
  px.reserve(p.points.size());
  py.reserve(p.points.size());
  wv.reserve(p.points.size());
  status.reserve(p.points.size());
  for (const OpdPoint& q : p.points) {
    px.push_back(q.px);
    py.push_back(q.py);
    wv.push_back(q.w);
    status.push_back(static_cast<std::uint8_t>(q.status));
  }
  w.f8("px", px);
  w.f8("py", py);
  w.f8("w", wv);
  w.u1("status", status);
  w.f8("sphere", sphere(p.sphere));
  w.u8("losses", losses(p.losses, p.warnings.size()));
}

}  // namespace

void run_pupil_cases(const fs::path& reference_dir,
                     const fs::path& catalog_dir,
                     const fs::path& out,
                     int threads) {
  const auto w = [&](const std::string& name) { return Writer(out, name); };
  const material::MaterialLibrary plain;
  material::MaterialLibrary schott;  // the full catalogue with N-LAK9 and N-SF5
  schott.add_catalog(catalog_dir / "m2" / "schott.agf");
  const CompiledSystem paraboloid =
      compile::compile(io::load_system(reference_dir / "m2/paraboloid_stop.rtt.json"), plain);
  const CompiledSystem cooke =
      compile::compile(io::load_system(reference_dir / "m2/cooke_triplet.rtt.json"), schott);
  const CompiledSystem singlet =
      compile::compile(io::load_system(reference_dir / "m1/singlet_const.rtt.json"), plain);
  oneapi::tbb::task_arena arena(threads);
  arena.execute([&] {
    PupilFunctionOptions para;
    para.grid = 24;
    para.with_prt = true;
    para.apodization = 0.5;
    write_pupil(w("paraboloid_pupil"), pupil_function(paraboloid, PathId{0}, 0, 0, para));
    PupilFunctionOptions o;
    o.grid = 24;
    o.with_prt = true;
    o.focus_shift = 0.3;
    const PupilFunction cooke_pf = pupil_function(cooke, PathId{0}, 2, 0, o);
    write_pupil(w("cooke_pupil"), cooke_pf);
    write_pupil(w("cooke_refocus"), refocus(cooke, cooke_pf, -0.2));
    write_zernike(w("cooke_zernike"), zernike_fit(cooke_pf));
    OpdOptions oo;
    oo.focus_shift = 0.2;
    write_points(w("singlet_points"),
                 opd_points(singlet, PathId{0}, 1, 0, trace::HexapolarPupil{4}, oo));
  });
}

}  // namespace rtt::py::reference
