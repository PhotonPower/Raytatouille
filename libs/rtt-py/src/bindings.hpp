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
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "rtt/compile/compiled_system.hpp"

namespace rtt::py {

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

/// Registers the translation of C++ exceptions into the classes of raytatouille.errors.
void register_errors(nanobind::module_& m);

void bind_model(nanobind::module_& m);
void bind_material(nanobind::module_& m);
void bind_compile(nanobind::module_& m);
void bind_paraxial(nanobind::module_& m);
void bind_trace(nanobind::module_& m);
void bind_analysis(nanobind::module_& m);
void bind_polar(nanobind::module_& m);
void bind_layout(nanobind::module_& m);

}  // namespace rtt::py
