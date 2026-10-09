// Writes the C++ results of the cases in test_bitwise.py (rays, #32) and
// test_bitwise_analysis.py (analyses, #33) as .npy files, so that pytest can compare the Python
// results bit for bit. Built in the same CMake tree with the same
// flags as the extension module raytatouille._core.
//
// Usage: rtt_py_reference <reference dir> <catalog dir> <output dir> <threads>
// Output: <output dir>/<case>.<column>.npy for every RayBatch column after the trace, the
// stats and the first-order values (see first_order_values()); the analysis cases are in
// analysis_cases.cpp.

#include <oneapi/tbb/task_arena.h>

#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "configuration_cases.hpp"
#include "crystal_cases.hpp"
#include "diffraction_cases.hpp"
#include "npy_writer.hpp"
#include "paths_cases.hpp"
#include "polar_batch.hpp"
#include "reports_cases.hpp"
#include "rtt/coating/catalog.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/paraxial/paraxial.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/ray_paths.hpp"
#include "rtt/trace/sequential.hpp"
#include "rtt/trace/sources.hpp"

namespace fs = std::filesystem;
using rtt::py::reference::write_column;
using rtt::trace::RayBatch;

namespace {

/// First-order values in the order of FIRST_ORDER_FIELDS in test_bitwise.py; NaN for None.
std::vector<double> first_order_values(const rtt::paraxial::FirstOrder& fo) {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const auto value = [&](const std::optional<double>& v) { return v.value_or(nan); };
  const auto pupil = [&](const std::optional<rtt::paraxial::Pupil>& p, bool diameter) {
    if (!p) return nan;
    return value(diameter ? p->diameter : p->z);
  };
  return {fo.object_index,
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
          pupil(fo.exit_pupil, true)};
}

/// Rays set by hand when a case has no pupil sampling.
enum class Hand { Singlet, Plate, Michelson };

struct Case {
  std::string name;
  std::string file;                                   // relative to the reference dir
  bool catalog;                                       // load <catalog dir>/schott.agf
  std::optional<rtt::trace::PupilSampling> sampling;  // none: rays set by hand, see `hand`
  std::optional<std::uint16_t> wavelength;            // none: reference wavelength
  rtt::trace::Aiming aiming;
  std::uint32_t path = 0;
  bool coatings = false;  // load <catalog dir>/coatings (ADR 0019)
  Hand hand = Hand::Singlet;
  bool paraxial = true;  // write the first-order values
};

/// Same cases as CASES in test_bitwise.py.
std::vector<Case> cases() {
  using rtt::trace::Aiming;
  return {
      {"singlet_hex", "m1/singlet_const.rtt.json", false, rtt::trace::HexapolarPupil{6},
       std::nullopt, Aiming::Real},
      {"paraboloid_random", "m2/paraboloid_stop.rtt.json", false, rtt::trace::RandomPupil{500, 42},
       std::nullopt, Aiming::Paraxial},
      {"achromat_grid", "m2/achromat.rtt.json", true, rtt::trace::GridPupil{15}, std::uint16_t{0},
       Aiming::Real},
      {"singlet_hand_filled", "m1/singlet_const.rtt.json", false, std::nullopt, std::nullopt,
       Aiming::Real},
      {"ar_singlet_hex", "m3/ar_singlet.rtt.json", false, rtt::trace::HexapolarPupil{6},
       std::uint16_t{0}, Aiming::Real, 0, true},
      {"absorbing_ar_plate_angles", "m3/absorbing_ar_plate.rtt.json", false, std::nullopt,
       std::nullopt, Aiming::Real, 0, true, Hand::Plate, false},
      {"michelson_reference_arm", "m0/michelson.rtt.json", false, std::nullopt, std::nullopt,
       Aiming::Real, 0, false, Hand::Michelson, false},
      {"michelson_test_arm", "m0/michelson.rtt.json", false, std::nullopt, std::nullopt,
       Aiming::Real, 1, false, Hand::Michelson, false},
      {"polarizer_qwp_hex", "m3/polarizer_qwp.rtt.json", false, rtt::trace::HexapolarPupil{4},
       std::nullopt, Aiming::Real},
  };
}

/// Rays set column by column, as hand_filled_rays() in test_bitwise.py: status and last_surface
/// vary (ALIVE at IMG, VIGNETTED at the stop for heights above its radius of 10 mm, MISSED
/// without any surface for rays travelling towards -z). All inputs are exact binary values or
/// the same decimal literals on both sides.
RayBatch hand_filled_rays() {
  RayBatch rays(48);
  for (std::size_t i = 0; i < rays.size(); ++i) {
    const double k = static_cast<double>(i);
    rays.pos_z()[i] = -10.0;
    if (i < 40) {  // parallel to the axis at y = 0 ... 19.5 mm, all three wavelengths
      rays.pos_y()[i] = 0.5 * k;
      rays.wl()[i] = static_cast<std::uint16_t>(i % 3);
    } else if (i < 44) {  // oblique in the x-z plane
      rays.pos_y()[i] = 2.0 * (k - 40.0);
      rays.dir_x()[i] = 0.6;
      rays.dir_z()[i] = 0.8;
      rays.wl()[i] = 1;
    } else {  // away from the system
      rays.pos_y()[i] = k - 44.0;
      rays.dir_z()[i] = -1.0;
      rays.wl()[i] = 1;
    }
  }
  return rays;
}

/// As plate_rays() in test_bitwise.py: 9 rays in the y-z plane at direction sines -0.4 ... 0.4.
RayBatch plate_rays() {
  RayBatch rays(9);
  for (std::size_t j = 0; j < rays.size(); ++j) {
    const double s = (static_cast<double>(j) - 4.0) / 10.0;
    rays.pos_y()[j] = 0.25 * (static_cast<double>(j) - 4.0);
    rays.dir_y()[j] = s;
    rays.dir_z()[j] = std::sqrt(1.0 - s * s);
  }
  return rays;
}

/// As michelson_rays() in test_bitwise.py: 8 rays tilted in the x-z plane.
RayBatch michelson_rays() {
  RayBatch rays(8);
  for (std::size_t j = 0; j < rays.size(); ++j) {
    const double k = static_cast<double>(j);
    const double s = (k - 4.0) / 200.0;
    rays.pos_x()[j] = 0.5 * k - 1.75;
    rays.pos_y()[j] = 0.125 * k;
    rays.dir_x()[j] = s;
    rays.dir_z()[j] = std::sqrt(1.0 - s * s);
  }
  return rays;
}

RayBatch hand_rays(Hand hand) {
  switch (hand) {
    case Hand::Singlet:
      return hand_filled_rays();
    case Hand::Plate:
      return plate_rays();
    case Hand::Michelson:
      return michelson_rays();
  }
  return hand_filled_rays();
}

/// Writes rows of `width` values as an (N, width) array.
template <typename T, typename Rows>
void write_rows(const fs::path& file, std::string_view descr, const Rows& rows, std::size_t width) {
  std::vector<T> data;
  data.reserve(width * rows.size());
  for (const auto& row : rows) {
    for (std::size_t c = 0; c < width; ++c) data.push_back(row[static_cast<int>(c)]);
  }
  const std::span<const std::byte> bytes = std::as_bytes(std::span<const T>(data));
  rtt::py::reference::write_npy(file, descr, {rows.size(), width},
                                {reinterpret_cast<const char*>(bytes.data()), bytes.size()});
}

/// raytatouille.polar on the traced rays, as polar_results() in test_bitwise.py.
template <typename File>
void write_polar(const RayBatch& rays, const File& file) {
  namespace polar = rtt::py::polar_batch;
  const std::vector<polar::CVec3> x =
      polar::transverse_polarization(rays, polar::CVec3(1.0, 0.0, 0.0));
  write_rows<double>(file("k0"), "<f8", polar::initial_directions(rays), 3);
  write_rows<rtt::math::Complex>(file("e_x"), "<c16", x, 3);
  const std::vector<double> unpolarized = polar::transmission(rays);
  write_column<double>(file("transmission"), "<f8", unpolarized);
  const std::vector<double> polarized = polar::transmission(rays, x);
  write_column<double>(file("transmission_x"), "<f8", polarized);
  const polar::Diattenuations d = polar::diattenuation(rays);
  write_column<double>(file("d_value"), "<f8", d.value);
  write_column<double>(file("d_maximum"), "<f8", d.maximum);
  write_column<double>(file("d_minimum"), "<f8", d.minimum);
  write_rows<rtt::math::Complex>(file("d_axis"), "<c16", d.axis, 3);
  const polar::Retardances r = polar::retardance(rays);
  write_column<double>(file("r_value"), "<f8", r.value);
  write_rows<rtt::math::Complex>(file("r_fast_axis"), "<c16", r.fast_axis, 3);
  write_rows<double>(file("stokes_x"), "<f8",
                     polar::stokes(rays, x, rtt::math::Vec3(1.0, 0.0, 0.0)), 4);
}

/// RayPaths of the recorded trace, as path_results() in test_bitwise.py.
template <typename File>
void write_paths(const rtt::trace::RayPaths& p, const File& file) {
  using rtt::py::reference::write_npy;
  const auto bytes = [](const auto& v) {
    const std::span<const std::byte> b = std::as_bytes(std::span(v));
    return std::span<const char>(reinterpret_cast<const char*>(b.data()), b.size());
  };
  const std::size_t n = p.ray_count();
  const std::size_t s = p.slots;
  static_assert(sizeof(std::size_t) == 8);
  write_npy(file("path_ray_indices"), "<u8", {n}, bytes(p.ray_indices));
  write_npy(file("path_event_surfaces"), "<u4", {p.event_surfaces.size()}, bytes(p.event_surfaces));
  write_npy(file("path_position"), "<f8", {n, s, 3}, bytes(p.position));
  write_npy(file("path_direction"), "<f8", {n, s, 3}, bytes(p.direction));
  write_npy(file("path_opl"), "<f8", {n, s}, bytes(p.opl));
  write_npy(file("path_weight"), "<f8", {n, s}, bytes(p.weight));
  write_npy(file("path_status"), "|u1", {n, s}, bytes(p.status));
  write_npy(file("path_count"), "<u4", {n}, bytes(p.count));
  write_npy(file("path_lost_at"), "<i4", {n}, bytes(p.lost_at));
}

void run(const Case& c,
         const fs::path& reference_dir,
         const fs::path& catalog_dir,
         const fs::path& out,
         int threads) {
  rtt::material::MaterialLibrary materials;
  if (c.catalog) materials.add_catalog(catalog_dir / "schott.agf");
  rtt::coating::CoatingLibrary coatings;
  if (c.coatings) coatings.add_catalog(catalog_dir / "coatings");
  const auto system =
      rtt::compile::compile(rtt::io::load_system(reference_dir / c.file), materials, coatings);
  const rtt::compile::PathId path{c.path};
  const std::uint16_t wl = c.wavelength.value_or(system.reference_wavelength());
  std::vector<std::uint16_t> fields(system.fields().points.size());
  std::iota(fields.begin(), fields.end(), std::uint16_t{0});

  oneapi::tbb::task_arena arena(threads);
  RayBatch rays;
  rtt::trace::TraceStats stats;
  rtt::trace::RayPaths paths;
  arena.execute([&] {
    rays = c.sampling ? rtt::trace::make_rays(system, path, fields, wl, *c.sampling, c.aiming)
                      : hand_rays(c.hand);
    // The paths come from a second, recorded trace of the same rays (#80); the columns below
    // from the plain trace, so the Python side (which records) also checks that recording
    // leaves the rays bitwise unchanged.
    RayBatch recorded = rays;
    static_cast<void>(rtt::trace::SequentialTracer{}.trace(system, path, recorded, paths));
    stats = rtt::trace::SequentialTracer{}.trace(system, path, rays);
  });

  const auto file = [&](std::string_view column) {
    return out / (c.name + "." + std::string(column) + ".npy");
  };
  const RayBatch& r = rays;
  write_column<double>(file("pos_x"), "<f8", r.pos_x());
  write_column<double>(file("pos_y"), "<f8", r.pos_y());
  write_column<double>(file("pos_z"), "<f8", r.pos_z());
  write_column<double>(file("dir_x"), "<f8", r.dir_x());
  write_column<double>(file("dir_y"), "<f8", r.dir_y());
  write_column<double>(file("dir_z"), "<f8", r.dir_z());
  write_column<std::uint16_t>(file("wl"), "<u2", r.wl());
  write_column<double>(file("opl"), "<f8", r.opl());
  write_column<double>(file("weight"), "<f8", r.weight());
  write_column<std::uint16_t>(file("field"), "<u2", r.field());
  write_column<double>(file("pupil_x"), "<f8", r.pupil_x());
  write_column<double>(file("pupil_y"), "<f8", r.pupil_y());
  write_column<std::uint32_t>(file("last_surface"), "<u4", r.last_surface());
  write_column<rtt::trace::RayStatus>(file("status"), "|u1", r.status());
  for (int row = 0; row < 3; ++row) {
    for (int col = 0; col < 3; ++col) {
      write_column<rtt::math::Complex>(file("prt" + std::to_string(row) + std::to_string(col)),
                                       "<c16", r.prt(row, col));
    }
  }
  std::vector<std::uint64_t> counts(stats.rays.begin(), stats.rays.end());
  write_column<std::uint64_t>(file("stats"), "<u8", counts);
  if (c.paraxial) {
    const std::vector<double> fo = first_order_values(rtt::paraxial::first_order(system, path, wl));
    write_column<double>(file("first_order"), "<f8", fo);
  }
  write_polar(r, file);
  write_paths(paths, file);
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const std::vector<std::string> args(argv + 1, argv + argc);
    if (args.size() != 4) {
      std::cerr << "usage: rtt_py_reference <reference dir> <catalog dir> <output dir> "
                   "<threads>\n";
      return 2;
    }
    const int threads = std::stoi(args[3]);
    if (threads < 1) throw std::invalid_argument("threads must be at least 1");
    fs::create_directories(args[2]);
    for (const Case& c : cases()) run(c, args[0], args[1], args[2], threads);
    rtt::py::reference::run_analysis_cases(args[0], args[1], args[2], threads);
    rtt::py::reference::run_paths_cases(args[0], args[1], args[2], threads);
    rtt::py::reference::run_crystal_cases(args[0], args[1], args[2], threads);
    rtt::py::reference::run_diffraction_cases(args[0], args[1], args[2], threads);
    rtt::py::reference::run_configuration_cases(args[0], args[1], args[2], threads);
    rtt::py::reference::run_reports_cases(args[0], args[1], args[2], threads);
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "rtt_py_reference: " << e.what() << '\n';
    return 1;
  }
}
