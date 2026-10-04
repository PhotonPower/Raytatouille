// Writes the C++ results of the cases in test_bitwise.py as .npy files, so that pytest can
// compare the Python results bit for bit (issue #32). Built in the same CMake tree with the same
// flags as the extension module raytatouille._core.
//
// Usage: rtt_py_reference <reference dir> <catalog dir> <output dir> <threads>
// Output: <output dir>/<case>.<column>.npy for every RayBatch column after the trace, the
// stats and the first-order values (see first_order_values()).

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

#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/paraxial/paraxial.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/sequential.hpp"
#include "rtt/trace/sources.hpp"

namespace fs = std::filesystem;
using rtt::trace::RayBatch;

namespace {

static_assert(std::endian::native == std::endian::little, ".npy descriptors below are '<'");

/// Writes a C-contiguous array as NumPy .npy file, format version 1.0.
void write_npy(const fs::path& file,
               std::string_view descr,
               const std::vector<std::size_t>& shape,
               std::span<const char> bytes) {
  std::string dims;
  for (const std::size_t d : shape) dims += std::to_string(d) + ", ";
  if (shape.size() > 1) dims.resize(dims.size() - 2);  // "(N,)" for 1-D, "(N, 3, 3)" else
  std::string header = "{'descr': '" + std::string(descr) + "', 'fortran_order': False, " +
                       "'shape': (" + dims + "), }";
  // Magic (6) + version (2) + header length (2) + header + '\n' is a multiple of 64.
  const std::size_t total = 10 + header.size() + 1;
  header.append((64 - total % 64) % 64, ' ');
  header += '\n';
  std::ofstream out(file, std::ios::binary);
  out.write("\x93NUMPY\x01\x00", 8);
  const auto length = static_cast<std::uint16_t>(header.size());
  const char length_bytes[2] = {static_cast<char>(length & 0xffU), static_cast<char>(length >> 8U)};
  out.write(length_bytes, 2);
  out.write(header.data(), static_cast<std::streamsize>(header.size()));
  out.write(bytes.data(), static_cast<std::streamsize>(bytes.size()));
  if (!out) throw std::runtime_error("cannot write " + file.string());
}

template <typename T>
void write_column(const fs::path& file, std::string_view descr, std::span<const T> values) {
  const std::span<const std::byte> bytes = std::as_bytes(values);
  write_npy(file, descr, {values.size()},
            {reinterpret_cast<const char*>(bytes.data()), bytes.size()});
}

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

struct Case {
  std::string name;
  std::string file;  // relative to the reference dir
  bool catalog;      // load <catalog dir>/schott.agf
  rtt::trace::PupilSampling sampling;
  std::optional<std::uint16_t> wavelength;  // none: reference wavelength
  rtt::trace::Aiming aiming;
};

/// Same cases as CASES in test_bitwise.py.
std::vector<Case> cases() {
  return {
      {"singlet_hex", "m1/singlet_const.rtt.json", false, rtt::trace::HexapolarPupil{6},
       std::nullopt, rtt::trace::Aiming::Real},
      {"paraboloid_random", "m2/paraboloid_stop.rtt.json", false, rtt::trace::RandomPupil{500, 42},
       std::nullopt, rtt::trace::Aiming::Paraxial},
      {"achromat_grid", "m2/achromat.rtt.json", true, rtt::trace::GridPupil{15}, std::uint16_t{0},
       rtt::trace::Aiming::Real},
  };
}

void run(const Case& c,
         const fs::path& reference_dir,
         const fs::path& catalog_dir,
         const fs::path& out,
         int threads) {
  rtt::material::MaterialLibrary materials;
  if (c.catalog) materials.add_catalog(catalog_dir / "schott.agf");
  const auto system =
      rtt::compile::compile(rtt::io::load_system(reference_dir / c.file), materials);
  const rtt::compile::PathId path{0};
  const std::uint16_t wl = c.wavelength.value_or(system.reference_wavelength());
  std::vector<std::uint16_t> fields(system.fields().points.size());
  std::iota(fields.begin(), fields.end(), std::uint16_t{0});

  oneapi::tbb::task_arena arena(threads);
  RayBatch rays;
  rtt::trace::TraceStats stats;
  arena.execute([&] {
    rays = rtt::trace::make_rays(system, path, fields, wl, c.sampling, c.aiming);
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
  const std::vector<double> fo = first_order_values(rtt::paraxial::first_order(system, path, wl));
  write_column<double>(file("first_order"), "<f8", fo);
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
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "rtt_py_reference: " << e.what() << '\n';
    return 1;
  }
}
