#pragma once

/// @file bindings.hpp
/// Internal declarations of the extension module raytatouille._core (ADR 0002). Each bind_*()
/// adds one area of the C++ API to the module; the Python modules in python/raytatouille
/// re-export the names.

#include <nanobind/nanobind.h>

#include <cstdint>
#include <optional>
#include <string>
#include <variant>

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

/// Registers the translation of C++ exceptions into the classes of raytatouille.errors.
void register_errors(nanobind::module_& m);

void bind_model(nanobind::module_& m);
void bind_compile(nanobind::module_& m);
void bind_paraxial(nanobind::module_& m);
void bind_trace(nanobind::module_& m);

}  // namespace rtt::py
