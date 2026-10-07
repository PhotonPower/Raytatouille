#pragma once

/// @file validate.hpp
/// Semantic checks of a model that structural parsing cannot express.

#include <cstdint>
#include <string>
#include <vector>

#include "rtt/diagnostics/codes.hpp"
#include "rtt/model/system.hpp"

namespace rtt::model {

/// Severity of a Diagnostic; defined with the code registry in rtt-diagnostics (ADR 0022).
using Severity = diagnostics::Severity;

/// Result of a semantic check.
struct Diagnostic {
  Severity severity = Severity::Error;
  std::string location;  ///< JSON pointer into the system file, e.g. "/root/children/0"
  std::string message;   ///< human-readable text; may change between versions, the code not
  /// Stable code from the registry rtt/diagnostics/codes.hpp, e.g. "material.unknown"
  /// (docs/diagnostics.md, ADR 0022); the severity is the one of the registry entry. Last
  /// member, so that {severity, location, message} keeps its meaning (code empty).
  std::string code;
};

/// Runs all semantic checks. An empty result means the model is consistent.
[[nodiscard]] std::vector<Diagnostic> validate(const System& system);

/// True if any diagnostic is an error.
[[nodiscard]] bool has_errors(const std::vector<Diagnostic>& diagnostics);

/// Single-line description "error [code] /root/...: message" (without brackets if the code is
/// empty).
[[nodiscard]] std::string to_string(const Diagnostic& diagnostic);

}  // namespace rtt::model
