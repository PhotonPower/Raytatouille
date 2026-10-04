#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/variant.h>
#include <nanobind/stl/vector.h>
#include <oneapi/tbb/task_arena.h>

#include <cstddef>
#include <cstdint>
#include <memory>
#include <numeric>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include "bindings.hpp"
#include "rtt/math/types.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/sequential.hpp"
#include "rtt/trace/sources.hpp"

namespace nb = nanobind;
using namespace nb::literals;

namespace rtt::py {
namespace {

using trace::RayBatch;
using trace::RayStatus;

template <typename T>
using Column = nb::ndarray<nb::numpy, T, nb::ndim<1>, nb::c_contig>;

/// NumPy view on one column of a RayBatch, without a copy. The array holds a reference to the
/// Python RayBatch object `owner`, so the batch lives at least as long as the array. A RayBatch
/// never changes its size from Python, so the view cannot dangle.
template <typename T>
Column<T> column(nb::handle owner, std::span<T> data) {
  const std::size_t shape[1] = {data.size()};
  return Column<T>(data.data(), 1, shape, owner);
}

/// Binds a column getter: name, member function of RayBatch, docstring (string literal).
template <typename T, typename Get>
void def_column(nb::class_<RayBatch>& cls, const char* name, Get get, const char* doc) {
  cls.def_prop_ro(
      name,
      [get](nb::pointer_and_handle<RayBatch> self) -> Column<T> {
        return column<T>(self.h, (self.p->*get)());
      },
      doc);
}

trace::TraceStats run_trace(const compile::CompiledSystem& system,
                            RayBatch& rays,
                            const PathArg& path,
                            std::optional<int> threads) {
  const compile::PathId id = path_id(system, path);
  const trace::SequentialTracer tracer;
  if (!threads) return tracer.trace(system, id, rays);
  if (*threads < 1) throw std::invalid_argument("threads must be at least 1");
  // Limits the threads like the C++ tests (test_sequential.cpp); the result does not depend on
  // it (ADR 0004, rule 7).
  oneapi::tbb::task_arena arena(*threads);
  return arena.execute([&] { return tracer.trace(system, id, rays); });
}

}  // namespace

void bind_trace(nb::module_& m) {
  nb::enum_<RayStatus>(m, "RayStatus", nb::is_arithmetic(),
                       "State of a ray; RayBatch.status holds these values as uint8.")
      .value("ALIVE", RayStatus::Alive, "still being traced")
      .value("MISSED", RayStatus::Missed, "did not intersect the next surface")
      .value("VIGNETTED", RayStatus::Vignetted, "hit a surface outside its aperture")
      .value("TIR", RayStatus::Tir, "total internal reflection where refraction was requested")
      .value("NO_CONVERGENCE", RayStatus::NoConvergence,
             "iterative intersection or aiming did not converge")
      .value("ABSORBED", RayStatus::Absorbed, "stopped by an absorbing surface")
      .value("EVENT_IMPOSSIBLE", RayStatus::EventImpossible, "the requested event cannot happen");
  m.attr("NO_SURFACE") = trace::kNoSurface;

  nb::enum_<trace::Aiming>(m, "Aiming", "How a ray is aimed into the stop.")
      .value("REAL", trace::Aiming::Real,
             "Newton iteration on the real ray onto the target point on the stop surface")
      .value("PARAXIAL", trace::Aiming::Paraxial,
             "straight line through the target point on the paraxial entrance pupil");

  // Pupil samplings; (px, py) are normalised pupil coordinates, the unit circle is the rim of
  // the paraxial entrance pupil and +y the meridional direction.
  nb::class_<trace::SinglePupilPoint>(m, "SinglePupilPoint",
                                      "One ray at normalised pupil coordinates (px, py).")
      .def(
          "__init__",
          [](trace::SinglePupilPoint* t, double px, double py) {
            new (t) trace::SinglePupilPoint{px, py};
          },
          "px"_a = 0.0, "py"_a = 0.0)
      .def_rw("px", &trace::SinglePupilPoint::px)
      .def_rw("py", &trace::SinglePupilPoint::py);
  nb::class_<trace::HexapolarPupil>(
      m, "HexapolarPupil",
      "Centre ray plus rings k = 1..rings with 6k rays on radius k/rings, starting at +y and "
      "turning towards +x.")
      .def(
          "__init__",
          [](trace::HexapolarPupil* t, int rings) { new (t) trace::HexapolarPupil{rings}; },
          "rings"_a = 6)
      .def_rw("rings", &trace::HexapolarPupil::rings);
  nb::class_<trace::GridPupil>(
      m, "GridPupil",
      "n x n points on [-1, 1]^2 (rows from py = -1, px fastest) inside the unit circle.")
      .def(
          "__init__", [](trace::GridPupil* t, int n) { new (t) trace::GridPupil{n}; }, "n"_a = 11)
      .def_rw("n", &trace::GridPupil::n);
  nb::class_<trace::FanXPupil>(m, "FanXPupil", "n points on px in [-1, 1] at py = 0.")
      .def(
          "__init__", [](trace::FanXPupil* t, int n) { new (t) trace::FanXPupil{n}; }, "n"_a = 11)
      .def_rw("n", &trace::FanXPupil::n);
  nb::class_<trace::FanYPupil>(m, "FanYPupil", "n points on py in [-1, 1] at px = 0.")
      .def(
          "__init__", [](trace::FanYPupil* t, int n) { new (t) trace::FanYPupil{n}; }, "n"_a = 11)
      .def_rw("n", &trace::FanYPupil::n);
  nb::class_<trace::RandomPupil>(
      m, "RandomPupil",
      "`count` points uniformly in the unit disk from std::mt19937_64 seeded with `seed` "
      "(reproducible; see rtt/trace/sources.hpp).")
      .def(
          "__init__",
          [](trace::RandomPupil* t, std::size_t count, std::uint64_t seed) {
            new (t) trace::RandomPupil{count, seed};
          },
          "count"_a = 100, "seed"_a = 0)
      .def_rw("count", &trace::RandomPupil::count)
      .def_rw("seed", &trace::RandomPupil::seed);

  nb::class_<RayBatch> batch(
      m, "RayBatch",
      "Rays as columns (structure of arrays). Positions in mm and unit directions in global "
      "coordinates (right-handed, optical axis +z), OPL in mm, weight as power (source = 1). "
      "Every column (pos_x, ..., last_surface, status, prt(row, col)) is a writable NumPy "
      "view on the batch without a copy and keeps the batch alive; the size is fixed from "
      "Python. Do not change a column from another Python thread while trace() runs on the "
      "batch (the GIL is released during the trace).");
  batch
      .def(nb::init<std::size_t>(), "size"_a,
           "Batch of `size` rays at the origin along +z, wavelength 0, OPL 0, weight 1, "
           "P = identity, status ALIVE, field 0, pupil (0, 0), last_surface NO_SURFACE.")
      .def("__len__", &RayBatch::size)
      .def_prop_ro("size", &RayBatch::size, "Number of rays.");
  using Double = std::span<double> (RayBatch::*)() noexcept;
  using U16 = std::span<std::uint16_t> (RayBatch::*)() noexcept;
  def_column<double>(batch, "pos_x", static_cast<Double>(&RayBatch::pos_x), "x in mm, global.");
  def_column<double>(batch, "pos_y", static_cast<Double>(&RayBatch::pos_y), "y in mm, global.");
  def_column<double>(batch, "pos_z", static_cast<Double>(&RayBatch::pos_z), "z in mm, global.");
  def_column<double>(batch, "dir_x", static_cast<Double>(&RayBatch::dir_x),
                     "x of the unit direction, global.");
  def_column<double>(batch, "dir_y", static_cast<Double>(&RayBatch::dir_y),
                     "y of the unit direction, global.");
  def_column<double>(batch, "dir_z", static_cast<Double>(&RayBatch::dir_z),
                     "z of the unit direction, global.");
  def_column<std::uint16_t>(batch, "wl", static_cast<U16>(&RayBatch::wl),
                            "Wavelength index into CompiledSystem.wavelengths_um (uint16).");
  def_column<double>(batch, "opl", static_cast<Double>(&RayBatch::opl),
                     "Accumulated optical path length in mm.");
  def_column<double>(batch, "weight", static_cast<Double>(&RayBatch::weight),
                     "Power weight, dimensionless (source normalised to 1).");
  def_column<std::uint16_t>(batch, "field", static_cast<U16>(&RayBatch::field),
                            "Index of the field point of the ray (uint16).");
  def_column<double>(batch, "pupil_x", static_cast<Double>(&RayBatch::pupil_x),
                     "Normalised pupil x of the ray's origin.");
  def_column<double>(batch, "pupil_y", static_cast<Double>(&RayBatch::pupil_y),
                     "Normalised pupil y of the ray's origin.");
  def_column<std::uint32_t>(
      batch, "last_surface",
      static_cast<std::span<std::uint32_t> (RayBatch::*)() noexcept>(&RayBatch::last_surface),
      "Index into CompiledSystem.surface_ids of the last surface hit, or NO_SURFACE (uint32).");
  batch.def_prop_ro(
      "status",
      [](nb::pointer_and_handle<RayBatch> self) -> Column<std::uint8_t> {
        // RayStatus has the underlying type std::uint8_t; access through an unsigned char type
        // is allowed for any object.
        const std::span<RayStatus> status = self.p->status();
        return column<std::uint8_t>(
            self.h, {reinterpret_cast<std::uint8_t*>(status.data()), status.size()});
      },
      "Ray status as uint8, values of RayStatus.");
  batch.def(
      "prt",
      [](nb::pointer_and_handle<RayBatch> self, int row, int col) -> Column<math::Complex> {
        if (row < 0 || row > 2 || col < 0 || col > 2) {
          throw std::out_of_range("prt: row and col must be 0, 1 or 2");
        }
        return column<math::Complex>(self.h, self.p->prt(row, col));
      },
      "row"_a, "col"_a,
      "Element (row, col) of the 3x3 polarization ray-tracing matrix P (Yun, McClain, "
      "Chipman 2011), global coordinates, complex128; view without a copy.");
  batch.def(
      "prt_matrices",
      [](const RayBatch& rays) {
        const std::size_t n = rays.size();
        auto data = std::make_unique<std::vector<math::Complex>>(n * 9);
        for (int row = 0; row < 3; ++row) {
          for (int col = 0; col < 3; ++col) {
            const auto element = rays.prt(row, col);
            for (std::size_t i = 0; i < n; ++i) (*data)[9 * i + 3 * row + col] = element[i];
          }
        }
        const std::size_t shape[3] = {n, 3, 3};
        math::Complex* values = data->data();
        // The capsule owns the copy from here on and frees it with the array.
        nb::capsule owner(data.get(), [](void* p) noexcept {
          const std::unique_ptr<std::vector<math::Complex>> owned(
              static_cast<std::vector<math::Complex>*>(p));
        });
        static_cast<void>(data.release());
        return nb::ndarray<nb::numpy, math::Complex, nb::shape<-1, 3, 3>, nb::c_contig>(
            values, 3, shape, owner);
      },
      "P of all rays as an (N, 3, 3) complex128 array. This is a copy: P is stored as nine "
      "separate columns.");

  nb::class_<trace::TraceStats>(m, "TraceStats", "Number of rays per status after a trace.")
      .def_ro("rays", &trace::TraceStats::rays, "Counts indexed by RayStatus.")
      .def("count", &trace::TraceStats::count, "status"_a, "Number of rays with `status`.");

  m.def(
      "make_rays",
      [](const compile::CompiledSystem& system, const trace::PupilSampling& sampling,
         const PathArg& path, std::optional<std::vector<std::uint16_t>> fields,
         std::optional<std::uint16_t> wavelength, trace::Aiming aiming) {
        std::vector<std::uint16_t> all;
        if (!fields) {
          all.resize(system.fields().points.size());
          std::iota(all.begin(), all.end(), std::uint16_t{0});
        }
        const std::vector<std::uint16_t>& selected = fields ? *fields : all;
        return trace::make_rays(system, path_id(system, path), selected,
                                wavelength_index(system, wavelength), sampling, aiming);
      },
      "system"_a, "sampling"_a, nb::kw_only(), "path"_a = 0, "fields"_a.none() = nb::none(),
      "wavelength"_a.none() = nb::none(), "aiming"_a = trace::Aiming::Real,
      nb::call_guard<nb::gil_scoped_release>(),
      "Rays for the field indices `fields` (None: all fields) at wavelength index "
      "`wavelength` (None: reference): for every field every pupil point of `sampling`, aimed "
      "with `aiming`.\n\nRaises ValueError for unknown paths, fields or wavelengths and "
      "ParaxialError for paths without paraxial data.");
  m.def("trace", &run_trace, "system"_a, "rays"_a, nb::kw_only(), "path"_a = 0,
        "threads"_a.none() = nb::none(), nb::call_guard<nb::gil_scoped_release>(),
        "Traces `rays` in place along `path` (index or name) with the sequential tracer and "
        "returns the counts per status. `threads` limits the worker threads (None: all); the "
        "result is bitwise the same for every number of threads. The GIL is released; do not "
        "change the columns of `rays` from another thread meanwhile.\n\nRaises ValueError for "
        "an unknown path name, a wavelength index that is not a system wavelength or an invalid "
        "status, IndexError for an unknown path index.");
}

}  // namespace rtt::py
