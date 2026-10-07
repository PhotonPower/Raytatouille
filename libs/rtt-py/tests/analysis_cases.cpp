// Analysis cases of rtt_py_reference (#33): the same calls as CASES in
// test_bitwise_analysis.py, every result written as named .npy arrays (see the write_*
// functions; the Python side flattens its results with the same names).

#include <oneapi/tbb/task_arena.h>

#include <cstdint>
#include <filesystem>
#include <functional>
#include <limits>
#include <optional>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "npy_writer.hpp"
#include "rtt/analysis/chromatic.hpp"
#include "rtt/analysis/field.hpp"
#include "rtt/analysis/opd.hpp"
#include "rtt/analysis/spot.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/paraxial/prescription.hpp"
#include "rtt/paraxial/seidel.hpp"

namespace rtt::py::reference {
namespace {

namespace fs = std::filesystem;
using namespace rtt::analysis;
using compile::CompiledSystem;
using compile::PathId;

constexpr std::uint64_t kNone = std::numeric_limits<std::uint64_t>::max();
constexpr double kNan = std::numeric_limits<double>::quiet_NaN();  // None in float arrays

/// Writes the arrays of one case as <out>/<case>.<name>.npy.
class Writer {
 public:
  Writer(fs::path out, std::string name) : out_(std::move(out)), name_(std::move(name)) {}

  void f8(const std::string& array, const std::vector<double>& v) const {
    write_column<double>(file(array), "<f8", v);
  }
  void u8(const std::string& array, const std::vector<std::uint64_t>& v) const {
    write_column<std::uint64_t>(file(array), "<u8", v);
  }
  template <typename Points, typename Get>
  void f8_of(const std::string& array, const Points& points, Get get) const {
    std::vector<double> v;
    for (const auto& p : points) v.push_back(get(p));
    f8(array, v);
  }
  template <typename Points, typename Get>
  void u8_of(const std::string& array, const Points& points, Get get) const {
    std::vector<std::uint64_t> v;
    for (const auto& p : points) v.push_back(static_cast<std::uint64_t>(get(p)));
    u8(array, v);
  }

 private:
  [[nodiscard]] fs::path file(const std::string& array) const {
    return out_ / (name_ + "." + array + ".npy");
  }
  fs::path out_;
  std::string name_;
};

// Flattening of the results; the names match flatten() in test_bitwise_analysis.py. Integers
// (indices, counts, status) as uint64, None as kNone.

void write_spot(const Writer& w, const SpotDiagram& s) {
  w.f8_of("x", s.points, [](const SpotPoint& p) { return p.x; });
  w.f8_of("y", s.points, [](const SpotPoint& p) { return p.y; });
  w.u8_of("wavelengths", s.points, [](const SpotPoint& p) { return p.wavelength; });
  w.f8_of("weight", s.points, [](const SpotPoint& p) { return p.weight; });
  w.f8("scalars",
       {s.chief.x, s.chief.y, s.stats.centroid.x, s.stats.centroid.y, s.stats.rms_centroid,
        s.stats.rms_chief, s.stats.geo_centroid, s.stats.geo_chief, s.vignetted_fraction});
  w.u8("ints", {s.field, s.wavelength ? *s.wavelength : kNone, s.image_surface, s.rays_launched,
                s.rays_arrived});
}

void write_fan_points(const Writer& w,
                      const std::string& prefix,
                      const std::vector<FanPoint>& points) {
  w.f8_of(prefix + "p", points, [](const FanPoint& p) { return p.p; });
  w.f8_of(prefix + "ex", points, [](const FanPoint& p) { return p.ex; });
  w.f8_of(prefix + "ey", points, [](const FanPoint& p) { return p.ey; });
  w.u8_of(prefix + "status", points, [](const FanPoint& p) { return p.status; });
}

void write_fan(const Writer& w, const RayFan& f) {
  write_fan_points(w, "t_", f.tangential);
  write_fan_points(w, "s_", f.sagittal);
  w.f8("scalars", {f.chief.x, f.chief.y});
  w.u8("ints", {f.field, f.wavelength, f.image_surface});
}

void write_opd_points(const Writer& w,
                      const std::string& prefix,
                      const std::vector<OpdPoint>& points) {
  w.f8_of(prefix + "px", points, [](const OpdPoint& p) { return p.px; });
  w.f8_of(prefix + "py", points, [](const OpdPoint& p) { return p.py; });
  w.f8_of(prefix + "w", points, [](const OpdPoint& p) { return p.w; });
  w.u8_of(prefix + "status", points, [](const OpdPoint& p) { return p.status; });
}

std::vector<double> sphere_values(const ReferenceSphere& s) {
  return {s.centre.x(), s.centre.y(), s.centre.z(), s.radius};
}

void write_opd_map(const Writer& w, const OpdMap& m) {
  write_opd_points(w, "", m.points);
  std::vector<double> scalars = sphere_values(m.sphere);
  scalars.push_back(m.rms);
  scalars.push_back(m.pv);
  w.f8("scalars", scalars);
  w.u8("ints", {m.field, m.wavelength, m.arrived, m.vignetted});
}

void write_opd_fan(const Writer& w, const OpdFan& f) {
  write_opd_points(w, "t_", f.tangential);
  write_opd_points(w, "s_", f.sagittal);
  w.f8("scalars", sphere_values(f.sphere));
  w.u8("ints", {f.field, f.wavelength});
}

void write_longitudinal(const Writer& w, const LongitudinalColour& l) {
  w.u8_of("f_wavelength", l.foci, [](const FocusPosition& f) { return f.wavelength; });
  w.f8_of("f_paraxial_z", l.foci, [](const FocusPosition& f) { return f.paraxial_z; });
  w.f8_of("f_real_z", l.foci, [](const FocusPosition& f) { return f.real_z; });
  w.f8("scalars", {l.paraxial, l.real});
  w.u8("ints", {l.pair.first, l.pair.second});
}

void write_lateral(const Writer& w, const LateralColour& l) {
  w.f8_of("chief_x", l.chief, [](const Point2& p) { return p.x; });
  w.f8_of("chief_y", l.chief, [](const Point2& p) { return p.y; });
  w.f8_of("offset_x", l.offset, [](const Point2& p) { return p.x; });
  w.f8_of("offset_y", l.offset, [](const Point2& p) { return p.y; });
  w.u8("ints", {l.field});
}

std::vector<double> distortion_values(const DistortionPoint& p) {
  return {p.fraction, p.field.x, p.field.y, p.real_height, p.paraxial_height, p.percent};
}

std::vector<double> field_curvature_values(const FieldCurvaturePoint& p) {
  return {p.fraction, p.field.x, p.field.y, p.tangential, p.sagittal, p.astigmatism};
}

/// One array per quantity over the sweep: "q0" ... "q5" in the order of the *_values().
template <typename Point, typename Values>
void write_sweep(const Writer& w, const std::vector<Point>& sweep, Values values) {
  for (std::size_t q = 0; q < 6; ++q) {
    w.f8_of("q" + std::to_string(q), sweep, [&](const Point& p) { return values(p)[q]; });
  }
}

void write_seidel(const Writer& w, const paraxial::Seidel& s) {
  using paraxial::SeidelSurface;
  w.u8_of("surface", s.surfaces, [](const SeidelSurface& q) { return q.surface; });
  w.f8_of("y", s.surfaces, [](const SeidelSurface& q) { return q.y; });
  w.f8_of("y_bar", s.surfaces, [](const SeidelSurface& q) { return q.y_bar; });
  w.f8_of("a", s.surfaces, [](const SeidelSurface& q) { return q.a; });
  w.f8_of("a_bar", s.surfaces, [](const SeidelSurface& q) { return q.a_bar; });
  w.f8_of("lagrange", s.surfaces, [](const SeidelSurface& q) { return q.lagrange; });
  w.f8_of("s1", s.surfaces, [](const SeidelSurface& q) { return q.terms.s1; });
  w.f8_of("s2", s.surfaces, [](const SeidelSurface& q) { return q.terms.s2; });
  w.f8_of("s3", s.surfaces, [](const SeidelSurface& q) { return q.terms.s3; });
  w.f8_of("s4", s.surfaces, [](const SeidelSurface& q) { return q.terms.s4; });
  w.f8_of("s5", s.surfaces, [](const SeidelSurface& q) { return q.terms.s5; });
  w.f8_of("c_l", s.surfaces, [](const SeidelSurface& q) { return q.terms.c_l; });
  w.f8_of("c_t", s.surfaces, [](const SeidelSurface& q) { return q.terms.c_t; });
  w.f8("scalars",
       {s.sum.s1, s.sum.s2, s.sum.s3, s.sum.s4, s.sum.s5, s.sum.c_l, s.sum.c_t, s.lagrange,
        s.marginal.z, s.marginal.y, s.marginal.u, s.chief.z, s.chief.y, s.chief.u});
  w.u8("ints",
       {s.chromatic ? s.chromatic->first : kNone, s.chromatic ? s.chromatic->second : kNone});
}

void write_prescription(const Writer& w, const paraxial::Prescription& p) {
  using paraxial::PrescriptionSurface;
  w.u8_of("surface", p.surfaces, [](const PrescriptionSurface& q) { return q.surface; });
  w.f8_of("z", p.surfaces, [](const PrescriptionSurface& q) { return q.z; });
  w.f8_of("n", p.surfaces, [](const PrescriptionSurface& q) { return q.n; });
  w.f8_of("y", p.surfaces,
          [](const PrescriptionSurface& q) { return q.marginal ? q.marginal->y : kNan; });
  w.f8_of("u", p.surfaces,
          [](const PrescriptionSurface& q) { return q.marginal ? q.marginal->u : kNan; });
  w.f8_of("i", p.surfaces,
          [](const PrescriptionSurface& q) { return q.marginal ? q.marginal->i : kNan; });
  w.f8_of("y_bar", p.surfaces,
          [](const PrescriptionSurface& q) { return q.chief ? q.chief->y : kNan; });
  w.f8_of("u_bar", p.surfaces,
          [](const PrescriptionSurface& q) { return q.chief ? q.chief->u : kNan; });
  w.f8_of("i_bar", p.surfaces,
          [](const PrescriptionSurface& q) { return q.chief ? q.chief->i : kNan; });
  w.f8_of("lagrange", p.surfaces,
          [](const PrescriptionSurface& q) { return q.lagrange.value_or(kNan); });
  w.f8("scalars",
       {p.total_track, p.object_distance.value_or(kNan), p.paraxial_working_f_number.value_or(kNan),
        p.paraxial_image_na.value_or(kNan), p.lagrange_invariant.value_or(kNan)});
}

/// Systems of the cases, as in systems() of test_bitwise_analysis.py.
struct Systems {
  CompiledSystem singlet;
  CompiledSystem achromat;
  CompiledSystem paraboloid21;  // mirror aperture radius 21 mm < beam radius 30 mm: vignetting
};

Systems load_systems(const fs::path& reference_dir, const fs::path& catalog_dir) {
  const material::MaterialLibrary plain;
  material::MaterialLibrary schott;
  schott.add_catalog(catalog_dir / "schott.agf");
  model::System paraboloid = io::load_system(reference_dir / "m2/paraboloid_stop.rtt.json");
  std::get<model::Element>(paraboloid.root.children[1].value).surfaces[0].aperture =
      model::CircularAperture{21.0, 0.0};
  return {compile::compile(io::load_system(reference_dir / "m1/singlet_const.rtt.json"), plain),
          compile::compile(io::load_system(reference_dir / "m2/achromat.rtt.json"), schott),
          compile::compile(paraboloid, plain)};
}

}  // namespace

void run_analysis_cases(const fs::path& reference_dir,
                        const fs::path& catalog_dir,
                        const fs::path& out,
                        int threads) {
  const Systems s = load_systems(reference_dir, catalog_dir);
  const PathId path{0};
  const auto ref = [](const CompiledSystem& c) { return c.reference_wavelength(); };
  const auto w = [&](const char* name) { return Writer(out, name); };
  const trace::Aiming real = trace::Aiming::Real;

  oneapi::tbb::task_arena arena(threads);
  arena.execute([&] {
    write_spot(w("spot_mono"),
               spot(s.singlet, path, 2, std::uint16_t{0}, {trace::HexapolarPupil{6}, real}));
    write_spot(w("spot_poly"),
               spot(s.achromat, path, 2, std::nullopt, {trace::RandomPupil{300, 7}, real}));
    write_spot(w("spot_vignetted"),
               spot(s.paraboloid21, path, 0, std::nullopt, {trace::GridPupil{21}, real}));
    write_fan(w("fan"), ray_fan(s.singlet, path, 1, 0, {21, real}));
    write_fan(w("fan_vignetted"),
              ray_fan(s.paraboloid21, path, 0, ref(s.paraboloid21), {41, trace::Aiming::Paraxial}));
    OpdOptions grid17;
    grid17.grid = 17;
    write_opd_map(w("opd_map"), opd_map(s.singlet, path, 2, ref(s.singlet), grid17));
    OpdOptions grid21;
    grid21.grid = 21;
    write_opd_map(w("opd_map_vignetted"),
                  opd_map(s.paraboloid21, path, 0, ref(s.paraboloid21), grid21));
    write_opd_fan(w("opd_fan"), opd_fan(s.achromat, path, 1, 2, {}));
    write_longitudinal(w("longitudinal"),
                       longitudinal_colour(s.achromat, path, {std::nullopt, 0.7, real}));
    write_lateral(w("lateral"), lateral_colour(s.achromat, path, 2, real));
    write_sweep(w("distortion"), distortion(s.singlet, path, ref(s.singlet), {11, real}),
                distortion_values);
    w("distortion_at")
        .f8("scalars", distortion_values(distortion_at(s.singlet, path, {0.0, 4.0, 1.0}, 0, real)));
    write_sweep(w("field_curvature"),
                field_curvature(s.achromat, path, ref(s.achromat), {7, 1e-3, real}),
                field_curvature_values);
    w("field_curvature_at")
        .f8("scalars", field_curvature_values(field_curvature_at(
                           s.achromat, path, {1.0, 3.0, 1.0}, ref(s.achromat), {11, 2e-3, real})));
    write_seidel(w("seidel"), paraxial::seidel(s.achromat, path, ref(s.achromat),
                                               paraxial::ChromaticPair{0, 2}));
    write_prescription(w("prescription"),
                       paraxial::prescription(s.achromat, path, ref(s.achromat)));
  });
}

}  // namespace rtt::py::reference
