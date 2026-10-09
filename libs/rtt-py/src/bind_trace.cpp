#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
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
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "bindings.hpp"
#include "rtt/math/types.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/ray_paths.hpp"
#include "rtt/trace/run_control.hpp"
#include "rtt/trace/sequential.hpp"
#include "rtt/trace/sources.hpp"

namespace nb = nanobind;
using namespace nb::literals;

namespace rtt::py {
namespace {

using trace::RayBatch;
using trace::RayStatus;

// RayBatch.status is exposed as uint8 through an unsigned char pointer (bind_trace below).
static_assert(std::is_same_v<std::underlying_type_t<RayStatus>, std::uint8_t> &&
              std::is_same_v<std::uint8_t, unsigned char>);

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

/// Read-only NumPy view of a RayPaths array with the given shape. The array holds a reference
/// to the Python RayPaths object `owner`; RayPaths is immutable from Python, so the view cannot
/// dangle.
template <typename T>
nb::ndarray<nb::numpy, const T, nb::c_contig> paths_view(nb::handle owner,
                                                         const T* data,
                                                         std::initializer_list<std::size_t> shape) {
  return nb::ndarray<nb::numpy, const T, nb::c_contig>(data, shape.size(), shape.begin(), owner);
}

using IndexArray = nb::ndarray<const std::int64_t, nb::ndim<1>, nb::c_contig, nb::device::cpu>;

std::tuple<trace::TraceStats, trace::RayPaths> run_trace_recorded(
    const compile::CompiledSystem& system,
    RayBatch& rays,
    const PathArg& path,
    std::optional<int> threads,
    const std::optional<IndexArray>& record_rays,
    std::size_t max_recorded_rays) {
  std::vector<std::size_t> selection;
  std::optional<std::span<const std::size_t>> selected;
  if (record_rays) {
    selection.reserve(record_rays->shape(0));
    for (std::size_t r = 0; r < record_rays->shape(0); ++r) {
      const std::int64_t i = record_rays->data()[r];
      if (i < 0) {
        throw std::invalid_argument("record_rays[" + std::to_string(r) +
                                    "] = " + std::to_string(i) + " is negative");
      }
      selection.push_back(static_cast<std::size_t>(i));
    }
    selected = selection;  // an empty selection is rejected by the tracer
  }
  const compile::PathId id = path_id(system, path);
  const trace::SequentialTracer tracer;
  trace::RayPaths paths;
  const trace::TraceStats stats = with_threads(
      threads, [&] { return tracer.trace(system, id, rays, paths, selected, max_recorded_rays); });
  return {stats, std::move(paths)};
}

trace::TraceStats run_trace(const compile::CompiledSystem& system,
                            RayBatch& rays,
                            const PathArg& path,
                            std::optional<int> threads,
                            const std::optional<trace::CancelToken>& cancel,
                            const std::optional<nb::callable>& progress) {
  const compile::PathId id = path_id(system, path);
  const trace::SequentialTracer tracer;
  const trace::RunControl control = run_control(cancel, progress);
  return released(threads, [&] { return tracer.trace(system, id, rays, control); });
}

}  // namespace

trace::RunControl run_control(const std::optional<trace::CancelToken>& cancel,
                              const std::optional<nb::callable>& progress) {
  trace::RunControl control;
  control.cancel = cancel;
  if (progress) {
    // By pointer: copying the callable in a worker would change its reference count without
    // the GIL. The temporaries of the call are destroyed before the GIL is released again.
    control.progress = [callable = &*progress](const trace::Progress& p) {
      const nb::gil_scoped_acquire acquire;
      (*callable)(p.done, p.total, nb::str(p.stage.data(), p.stage.size()));
    };
  }
  return control;
}

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
      .value("EVENT_IMPOSSIBLE", RayStatus::EventImpossible,
             "the requested event cannot happen (a modelling problem; TIR and EVANESCENT are "
             "physical limits)")
      .value("EVANESCENT", RayStatus::Evanescent,
             "the diffraction order has no real direction (ADR 0025)");
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
  nb::class_<trace::GaussPupil>(
      m, "GaussPupil",
      "Gaussian quadrature in the pupil (#168): `rings` rings at the Gauss-Legendre nodes in "
      "u = rho^2 (ascending) and `arms` arms per ring at 2 pi j / arms from +y towards +x. Use "
      "the weights of gauss_pupil_weights() for means over the pupil; spot, the path analyses "
      "and the ghost ranking reject this sampling.")
      .def(
          "__init__",
          [](trace::GaussPupil* t, int rings, int arms) { new (t) trace::GaussPupil{rings, arms}; },
          "rings"_a = 3, "arms"_a = 6)
      .def_rw("rings", &trace::GaussPupil::rings)
      .def_rw("arms", &trace::GaussPupil::arms);
  m.def(
      "gauss_pupil_weights", &trace::gauss_pupil_weights, "gauss"_a,
      "Quadrature weights of the points of `gauss` in the order of its rays (they sum to 1): "
      "sum q f approximates the mean of f over the unit pupil, exactly for polynomials in "
      "rho^2 up to degree 2 rings - 1 and angular orders below arms. Not the ray weights (power). "
      "Raises ValueError for rings or arms < 1.");

  nb::class_<RayBatch> batch(
      m, "RayBatch",
      "Rays as columns (structure of arrays). Positions in mm and unit directions in global "
      "coordinates (right-handed, optical axis +z), OPL in mm. weight is the power for an "
      "unpolarized source (source = 1); P (prt) is power-normalised, see raytatouille.polar and "
      "ADR 0021. "
      "Every column (pos_x, ..., wave_x, ..., mode_index, last_surface, status, "
      "prt(row, col)) is a writable NumPy "
      "view on the batch without a copy and keeps the batch alive; the size is fixed from "
      "Python. While trace() or make_rays() runs on the batch (the GIL is released), do not "
      "read or change its columns from another Python thread and do not start a second "
      "trace() on it.");
  batch
      .def(nb::init<std::size_t>(), "size"_a,
           "Batch of `size` rays at the origin along +z, wavelength 0, OPL 0, weight 1, "
           "P = identity, status ALIVE, field 0, pupil (0, 0), last_surface NO_SURFACE, "
           "wave +z and mode_index 0.")
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
  def_column<double>(
      batch, "weight", static_cast<Double>(&RayBatch::weight),
      "Power for an unpolarized source, dimensionless (source = 1): weight = s ||P_T||^2 / 2 "
      "with P_T = (I - k k^T) P and s the polarization-independent factors "
      "(volume absorption, absorber); ADR 0021. For a polarized state use "
      "raytatouille.polar.transmission.");
  def_column<std::uint16_t>(batch, "field", static_cast<U16>(&RayBatch::field),
                            "Index of the field point of the ray (uint16).");
  def_column<double>(batch, "pupil_x", static_cast<Double>(&RayBatch::pupil_x),
                     "Normalised pupil x of the ray's origin.");
  def_column<double>(batch, "pupil_y", static_cast<Double>(&RayBatch::pupil_y),
                     "Normalised pupil y of the ray's origin.");
  def_column<double>(batch, "wave_x", static_cast<Double>(&RayBatch::wave_x),
                     "x of the unit wave normal k, global (ADR 0026, point 3). dir is the "
                     "energy direction S; in an isotropic medium wave = dir. Counts only where "
                     "mode_index > 0: trace() sets wave := dir for every ray with mode_index 0, "
                     "so a batch that sets only dir is valid (reading rule).");
  def_column<double>(batch, "wave_y", static_cast<Double>(&RayBatch::wave_y),
                     "y of the unit wave normal k, global; see wave_x.");
  def_column<double>(batch, "wave_z", static_cast<Double>(&RayBatch::wave_z),
                     "z of the unit wave normal k, global; see wave_x.");
  def_column<double>(batch, "mode_index", static_cast<Double>(&RayBatch::mode_index),
                     "Index n of the crystal mode along the wave normal, real, dimensionless: "
                     "> 0 inside a uniaxial crystal, 0 in an isotropic medium (no mode; "
                     "ADR 0026, point 3). OPL in the crystal grows by n l (k . S).");
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
      "Chipman 2011), global coordinates, complex128; view without a copy. P is "
      "power-normalised (ADR 0021): |P E|^2 is the power fraction for an incident state E "
      "(|E| = 1, transverse to the initial direction k0) without the polarization-independent "
      "factors in weight; the phases are those of the field PRT matrix, the field amplitude "
      "differs by the factors sqrt(c) of the interfaces. P k0 = k.");
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
      "P of all rays as an (N, 3, 3) complex128 array, global coordinates, power-normalised as "
      "described at prt() (ADR 0021). This is a copy: P is stored as nine separate columns.");

  nb::class_<trace::RayPaths>(
      m, "RayPaths",
      "Recorded paths of traced rays (#80): the state before the first event (slot 0) and after "
      "each event of the path (slot k after event k - 1), for drawing rays. Positions in mm and "
      "unit directions in global coordinates, OPL in mm, weight as power for an unpolarized "
      "source (ADR 0021). A ray lost at event j has count = j + 2 valid slots and lost_at = j; "
      "slot j + 1 holds the state at the loss with the loss status, later slots are NaN with the "
      "loss status. A ray that did not start has count 1, a ray that passed all events count = "
      "slots; lost_at is -1 for both. Slot count - 1 is bitwise the final state in the "
      "RayBatch. All arrays are read-only views without a copy and keep the object alive.")
      .def_prop_ro(
          "slots", [](const trace::RayPaths& p) { return p.slots; },
          "Number of slots S = number of path events + 1.")
      .def_prop_ro(
          "n_events", [](const trace::RayPaths& p) { return p.event_surfaces.size(); },
          "Number of path events (S - 1).")
      .def_prop_ro("ray_count", &trace::RayPaths::ray_count, "Number of recorded rays N.")
      .def("__len__", &trace::RayPaths::ray_count)
      .def_prop_ro(
          "ray_indices",
          [](nb::pointer_and_handle<trace::RayPaths> self) {
            return paths_view<std::size_t>(self.h, self.p->ray_indices.data(),
                                           {self.p->ray_count()});
          },
          "RayBatch index of every recorded ray, (N,) uint64.")
      .def_prop_ro(
          "event_surfaces",
          [](nb::pointer_and_handle<trace::RayPaths> self) {
            return paths_view<std::uint32_t>(self.h, self.p->event_surfaces.data(),
                                             {self.p->event_surfaces.size()});
          },
          "Surface index (into CompiledSystem.surface_ids) of every path event, (S - 1,) uint32.")
      .def_prop_ro(
          "position",
          [](nb::pointer_and_handle<trace::RayPaths> self) {
            return paths_view<double>(self.h, self.p->position.data(),
                                      {self.p->ray_count(), self.p->slots, 3});
          },
          "Position per ray and slot, (N, S, 3) float64, mm, global; NaN after a loss.")
      .def_prop_ro(
          "direction",
          [](nb::pointer_and_handle<trace::RayPaths> self) {
            return paths_view<double>(self.h, self.p->direction.data(),
                                      {self.p->ray_count(), self.p->slots, 3});
          },
          "Unit direction per ray and slot, (N, S, 3) float64, global; NaN after a loss.")
      .def_prop_ro(
          "opl",
          [](nb::pointer_and_handle<trace::RayPaths> self) {
            return paths_view<double>(self.h, self.p->opl.data(),
                                      {self.p->ray_count(), self.p->slots});
          },
          "Accumulated optical path length per ray and slot, (N, S) float64, mm; NaN after a "
          "loss.")
      .def_prop_ro(
          "weight",
          [](nb::pointer_and_handle<trace::RayPaths> self) {
            return paths_view<double>(self.h, self.p->weight.data(),
                                      {self.p->ray_count(), self.p->slots});
          },
          "Power for an unpolarized source per ray and slot, (N, S) float64 (ADR 0021); NaN "
          "after a loss.")
      .def_prop_ro(
          "status",
          [](nb::pointer_and_handle<trace::RayPaths> self) {
            // RayStatus has the underlying type std::uint8_t; access through an unsigned char
            // type is allowed for any object.
            const auto* data = reinterpret_cast<const std::uint8_t*>(self.p->status.data());
            return paths_view<std::uint8_t>(self.h, data, {self.p->ray_count(), self.p->slots});
          },
          "Status per ray and slot, (N, S) uint8, values of RayStatus; after a loss the loss "
          "status.")
      .def_prop_ro(
          "count",
          [](nb::pointer_and_handle<trace::RayPaths> self) {
            return paths_view<std::uint32_t>(self.h, self.p->count.data(), {self.p->ray_count()});
          },
          "Number of valid slots per ray, (N,) uint32.")
      .def_prop_ro(
          "lost_at",
          [](nb::pointer_and_handle<trace::RayPaths> self) {
            return paths_view<std::int32_t>(self.h, self.p->lost_at.data(), {self.p->ray_count()});
          },
          "Index of the event at which each ray stopped, (N,) int32; -1 if it passed all events "
          "or did not start. The surface is event_surfaces[lost_at].");

  nb::class_<trace::CancelToken>(
      m, "CancelToken",
      "Cancellation request for a running trace, make_rays or analysis (#83). Pass it as "
      "`cancel=`; cancel() from any thread (e.g. a GUI) makes the run stop after at most one "
      "block of rays per worker and raise raytatouille.errors.Cancelled.")
      .def(nb::init<>())
      .def("cancel", &trace::CancelToken::request_cancel,
           "Requests cancellation; idempotent, never blocks.")
      .def_prop_ro("cancelled", &trace::CancelToken::cancelled,
                   "Whether cancellation was requested.");

  m.attr("DEFAULT_MAX_RECORDED_RAYS") = trace::kDefaultMaxRecordedRays;
  m.def("trace_recorded", &run_trace_recorded, "system"_a, "rays"_a, nb::kw_only(), "path"_a = 0,
        "threads"_a.none() = nb::none(), "record_rays"_a.none() = nb::none(),
        "max_recorded_rays"_a = trace::kDefaultMaxRecordedRays,
        nb::call_guard<nb::gil_scoped_release>(),
        "trace() that also records the paths of the rays `record_rays` (int64 indices into "
        "`rays`, no duplicates; None: all rays) and returns (TraceStats, RayPaths). The rays end "
        "bitwise as without recording. Use raytatouille.trace.trace(..., record_path=True).\n\n"
        "Raises ValueError for an index outside the batch, a duplicate, an empty selection, or "
        "more than `max_recorded_rays` rays without a selection (the message gives the memory "
        "need), and as trace().");

  nb::class_<trace::TraceStats>(m, "TraceStats", "Number of rays per status after a trace.")
      .def_ro("rays", &trace::TraceStats::rays, "Counts indexed by RayStatus.")
      .def("count", &trace::TraceStats::count, "status"_a, "Number of rays with `status`.");

  m.def(
      "make_rays",
      [](const compile::CompiledSystem& system, const trace::PupilSampling& sampling,
         const PathArg& path, std::optional<std::vector<std::uint16_t>> fields,
         std::optional<std::uint16_t> wavelength, trace::Aiming aiming, std::optional<int> threads,
         const std::optional<trace::CancelToken>& cancel,
         const std::optional<nb::callable>& progress) {
        std::vector<std::uint16_t> all;
        if (!fields) {
          all.resize(system.fields().points.size());
          std::iota(all.begin(), all.end(), std::uint16_t{0});
        }
        const std::vector<std::uint16_t>& selected = fields ? *fields : all;
        const trace::RunControl control = run_control(cancel, progress);
        return released(threads, [&] {
          return trace::make_rays(system, path_id(system, path), selected,
                                  wavelength_index(system, wavelength), sampling, aiming, control);
        });
      },
      "system"_a, "sampling"_a, nb::kw_only(), "path"_a = 0, "fields"_a.none() = nb::none(),
      "wavelength"_a.none() = nb::none(), "aiming"_a = trace::Aiming::Real,
      "threads"_a.none() = nb::none(), "cancel"_a.none() = nb::none(),
      "progress"_a.none() = nb::none(),
      "Rays for the field indices `fields` (None: all fields) at wavelength index "
      "`wavelength` (None: reference): for every field every pupil point of `sampling`, aimed "
      "with `aiming`, in parallel (#119). `threads` limits the worker threads (None: all); the "
      "rays are bitwise the same for every number of threads. `cancel` (CancelToken) and "
      "`progress(done, total, stage)` (stage "
      "'aim', called from any thread with the GIL) control the run; neither changes the "
      "rays.\n\nRaises ValueError for unknown paths, fields or wavelengths, ParaxialError for "
      "paths without paraxial data, ValueError for threads < 1, raytatouille.errors.Cancelled "
      "after a cancellation and the exception of `progress`.");
  m.def("trace", &run_trace, "system"_a, "rays"_a, nb::kw_only(), "path"_a = 0,
        "threads"_a.none() = nb::none(), "cancel"_a.none() = nb::none(),
        "progress"_a.none() = nb::none(),
        "Traces `rays` in place along `path` (index or name) with the sequential tracer and "
        "returns the counts per status. `threads` limits the worker threads (None: all); the "
        "result is bitwise the same for every number of threads. The GIL is released; do not "
        "read or change the columns of `rays` from another thread meanwhile. `cancel` and "
        "`progress` as in make_rays (stage 'trace'); neither changes the result.\n\nRaises "
        "ValueError for an unknown path name, a wavelength index that is not a system "
        "wavelength or an invalid status, IndexError for an unknown path index, "
        "raytatouille.errors.Cancelled after a cancellation and the exception of `progress`.");
}

}  // namespace rtt::py
