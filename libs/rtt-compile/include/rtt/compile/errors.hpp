#pragma once

/// @file errors.hpp
/// Errors of the tracing layer that rtt-paraxial, rtt-trace and rtt-analysis share (ADR 0022).

#include <stdexcept>
#include <string>

#include "rtt/compile/compiled_system.hpp"

namespace rtt::compile {

/// The path has no stop: no surface of a Stop element on it. Aiming at the stop, the pupils,
/// Seidel sums and OPD need one (first_order() does not: its pupils are then empty).
/// A std::invalid_argument, because the input system does not allow the requested analysis.
class NoStopError : public std::invalid_argument {
 public:
  /// `path_name` as in the model, `location` the JSON pointer of the path, e.g. "/paths/0".
  NoStopError(std::string path_name, std::string location);

  /// Name of the path without stop.
  [[nodiscard]] const std::string& path_name() const noexcept { return path_name_; }
  /// JSON pointer of the path in the system file, e.g. "/paths/0".
  [[nodiscard]] const std::string& location() const noexcept { return location_; }

 private:
  std::string path_name_;
  std::string location_;
};

/// The one check for "the path has a stop" (ADR 0022): a path has a stop if one of its events
/// is at a surface of a Stop element. Called by everything that needs the stop (aimed ray
/// sources, Seidel sums, OPD) before any ray is traced, so that all of them throw the same
/// error. Their own argument checks (path, field, wavelength, rotational symmetry, an off-axis
/// field for distortion) may report first.
/// @throws NoStopError if the path has no stop
/// @throws std::out_of_range if `path` does not belong to `system`
void require_stop(const CompiledSystem& system, PathId path);

}  // namespace rtt::compile
