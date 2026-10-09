#pragma once

/// @file bindings.hpp
/// Internal declarations of the extension module raytatouille._core (ADR 0002). Each bind_*()
/// adds one area of the C++ API to the module; the Python modules in python/raytatouille
/// re-export the names.

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <oneapi/tbb/task_arena.h>

#include <cstddef>
#include <cstdint>
#include <iterator>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/model/system.hpp"
#include "rtt/model/validate.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/run_control.hpp"

namespace rtt::py {

/// Issues every warning of `diagnostics` as a raytatouille.errors.RaytatouilleWarning (ADR 0022).
/// Called from C++, stacklevel 1 is already the Python line that called the bound function. An
/// error raised by a warnings filter ("error") propagates as that Python exception.
void warn(const std::vector<model::Diagnostic>& diagnostics);

/// A path given from Python: index into CompiledSystem paths or path name.
using PathArg = std::variant<std::uint32_t, std::string>;

/// Resolves a path argument.
/// @throws std::invalid_argument if a name is not a path of `system` (an unknown index is
///         reported by the called C++ function)
[[nodiscard]] compile::PathId path_id(const compile::CompiledSystem& system, const PathArg& path);

/// Wavelength index from Python: None means the reference wavelength.
[[nodiscard]] std::uint16_t wavelength_index(const compile::CompiledSystem& system,
                                             std::optional<std::uint16_t> wavelength);

/// Runs `f` with at most `threads` worker threads (None: all). Limits the threads like the C++
/// tests (test_sequential.cpp); results do not depend on it (ADR 0004, rule 7).
/// @throws std::invalid_argument if threads < 1
template <typename F>
auto with_threads(std::optional<int> threads, F&& f) {
  if (!threads) return std::forward<F>(f)();
  if (*threads < 1) throw std::invalid_argument("threads must be at least 1");
  oneapi::tbb::task_arena arena(*threads);
  return arena.execute(std::forward<F>(f));
}

/// Run-time control from Python (#83): the token `cancel` and the callable
/// `progress(done, total, stage)`, both optional. The control refers to `progress` by pointer,
/// so the workers never touch its reference count without the GIL: `progress` must outlive the
/// run (an argument of the bound function does). The callback acquires the GIL; a Python
/// exception from it ends the run and propagates (nanobind::python_error).
[[nodiscard]] trace::RunControl run_control(const std::optional<trace::CancelToken>& cancel,
                                            const std::optional<nanobind::callable>& progress);

/// Releases the GIL and runs `f` with at most `threads` worker threads (with_threads). For the
/// functions with a RunControl instead of a nanobind::call_guard: the control is built with the
/// GIL held (run_control) and only the computation runs without it, which the explicit release
/// makes visible. The control holds no Python reference.
template <typename F>
auto released(std::optional<int> threads, F&& f) {
  const nanobind::gil_scoped_release release;
  return with_threads(threads, std::forward<F>(f));
}

/// Read-only one-dimensional NumPy array.
template <typename T>
using ReadOnlyArray =
    nanobind::ndarray<nanobind::numpy, const T, nanobind::ndim<1>, nanobind::c_contig>;

/// Read-only NumPy copy of `get(item)` for every item of `items` (analysis results, ADR 0002:
/// small, so a copy; RayBatch columns are views instead).
template <typename T, typename Items, typename Get>
ReadOnlyArray<T> read_only_array(const Items& items, Get get) {
  auto data = std::make_unique<std::vector<T>>();
  data->reserve(std::size(items));
  for (const auto& item : items) data->push_back(static_cast<T>(get(item)));
  const std::size_t shape[1] = {data->size()};
  const T* values = data->data();
  // The capsule owns the copy from here on and frees it with the array.
  nanobind::capsule owner(data.get(), [](void* p) noexcept {
    const std::unique_ptr<std::vector<T>> owned(static_cast<std::vector<T>*>(p));
  });
  static_cast<void>(data.release());
  return ReadOnlyArray<T>(values, 1, shape, owner);
}

/// NaN for "none" in float arrays of analysis results (ADR 0023, as Prescription, #108).
inline constexpr double kNan = std::numeric_limits<double>::quiet_NaN();

/// Binds a list class (a struct with a member `points`); its attributes are read-only NumPy
/// copies of the C++ result (ADR 0002: analysis results are small), one array per member of
/// the point type.
template <typename C>
nanobind::class_<C> columns(nanobind::module_& m, const char* name, const char* doc) {
  nanobind::class_<C> cls(m, name, doc);
  cls.def("__len__", [](const C& c) { return c.points.size(); });
  return cls;
}

/// Adds a read-only array attribute `name` with values get(point) of type T.
template <typename T, typename C, typename Get>
void column(nanobind::class_<C>& cls, const char* name, Get get, const char* doc) {
  // reference: the NumPy array views the copy and keeps its capsule alive (move would copy
  // again into a writable array; reference_internal needs an array without owner).
  cls.def_prop_ro(
      name, [get](const C& c) { return read_only_array<T>(c.points, get); },
      nanobind::rv_policy::reference, doc);
}

/// RayStatus as the uint8 of the status columns.
[[nodiscard]] inline std::uint8_t status_value(trace::RayStatus s) {
  return static_cast<std::uint8_t>(s);
}

/// __eq__ with Python semantics: comparing with another type gives False, not TypeError.
template <typename T>
bool equal(const T& self, nanobind::handle other) {
  return nanobind::isinstance<T>(other) && self == nanobind::cast<const T&>(other);
}

/// Registers the translation of C++ exceptions into the classes of raytatouille.errors.
void register_errors(nanobind::module_& m);

void bind_model(nanobind::module_& m);
void bind_model_tree(nanobind::module_& m);
void bind_edit(nanobind::module_& m);

/// The edit form of `s` (ADR 0024) as Python objects (json.loads of rtt::io::to_edit_json).
[[nodiscard]] nanobind::object edit_dict(const model::System& s);
/// The value at a JSON pointer (RFC 6901) in the edit form of `s`.
/// Raises KeyError if it does not exist and ValueError for an invalid pointer or array index.
[[nodiscard]] nanobind::object json_at(const model::System& s, std::string_view pointer);
/// Current JSON pointer of a surface id, a node (assembly or element) name or a path name;
/// std::nullopt (None in Python) if there is none (ADR 0024 point 2: separate name spaces).
[[nodiscard]] std::optional<std::string> locate_surface(const model::System& s,
                                                        std::string_view id);
[[nodiscard]] std::optional<std::string> locate_node(const model::System& s, std::string_view name);
[[nodiscard]] std::optional<std::string> locate_path(const model::System& s, std::string_view name);
void bind_material(nanobind::module_& m);
void bind_compile(nanobind::module_& m);
void bind_paraxial(nanobind::module_& m);
void bind_trace(nanobind::module_& m);
void bind_analysis(nanobind::module_& m);
/// Path evaluation (#122) and ghost ranking (#124); after bind_analysis (RayLosses).
void bind_paths(nanobind::module_& m);
/// Reports as data (#177); after bind_analysis (Prescription).
void bind_reports(nanobind::module_& m);
/// Optimization (#169, ADR 0030); after bind_model (Diagnostic) and bind_trace (CancelToken).
void bind_optim(nanobind::module_& m);
void bind_polar(nanobind::module_& m);
void bind_layout(nanobind::module_& m);
/// to_dict()/to_json() on the result classes (ADR 0023); call after all bind_* functions.
void bind_result_methods(nanobind::module_& m);

}  // namespace rtt::py
