#pragma once

/// @file result_names.hpp
/// Names of the enums in the result format (ADR 0023; ADR 0052, point 1): the C++ result writer
/// (result_json.hpp) and the Python binding take them from here, so that both spell every member
/// the same way. Python writes an enum as its member name (raytatouille.results, to_dict).

#include <array>
#include <cstddef>
#include <string_view>

#include "rtt/diagnostics/codes.hpp"
#include "rtt/optim/levenberg_marquardt.hpp"
#include "rtt/trace/ray_batch.hpp"

namespace rtt::io {

/// Name of a status of the Levenberg-Marquardt run (Python: raytatouille.OptimStatus).
[[nodiscard]] constexpr std::string_view result_name(optim::LmStatus s) noexcept {
  switch (s) {
    case optim::LmStatus::ConvergedGradient:
      return "CONVERGED_GRADIENT";
    case optim::LmStatus::ConvergedStep:
      return "CONVERGED_STEP";
    case optim::LmStatus::ConvergedMerit:
      return "CONVERGED_MERIT";
    case optim::LmStatus::MaxIterations:
      return "MAX_ITERATIONS";
    case optim::LmStatus::Cancelled:
      return "CANCELLED";
    case optim::LmStatus::Failed:
      return "FAILED";
  }
  return "?";
}

/// Name of a severity (Python: raytatouille.Severity).
[[nodiscard]] constexpr std::string_view result_name(diagnostics::Severity s) noexcept {
  switch (s) {
    case diagnostics::Severity::Error:
      return "ERROR";
    case diagnostics::Severity::Warning:
      return "WARNING";
  }
  return "?";
}

/// Name of a ray status (Python: raytatouille.RayStatus); in the results the status is a uint8
/// value, the names serve the CSV of the command line (ADR 0052, point 5).
[[nodiscard]] constexpr std::string_view result_name(trace::RayStatus s) noexcept {
  switch (s) {
    case trace::RayStatus::Alive:
      return "ALIVE";
    case trace::RayStatus::Missed:
      return "MISSED";
    case trace::RayStatus::Vignetted:
      return "VIGNETTED";
    case trace::RayStatus::Tir:
      return "TIR";
    case trace::RayStatus::NoConvergence:
      return "NO_CONVERGENCE";
    case trace::RayStatus::Absorbed:
      return "ABSORBED";
    case trace::RayStatus::EventImpossible:
      return "EVENT_IMPOSSIBLE";
    case trace::RayStatus::Evanescent:
      return "EVANESCENT";
  }
  return "?";
}

}  // namespace rtt::io
