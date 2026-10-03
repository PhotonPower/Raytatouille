#pragma once

/// @file validate.hpp
/// Semantic checks of a model that structural parsing cannot express.

#include <cstdint>
#include <string>
#include <vector>

#include "rtt/model/system.hpp"

namespace rtt::model {

enum class Severity : std::uint8_t { Error, Warning };

struct Diagnostic {
  Severity severity = Severity::Error;
  std::string location;  ///< JSON pointer into the system file, e.g. "/root/children/0"
  std::string message;
};

/// Runs all semantic checks. An empty result means the model is consistent.
[[nodiscard]] std::vector<Diagnostic> validate(const System& system);

/// True if any diagnostic is an error.
[[nodiscard]] bool has_errors(const std::vector<Diagnostic>& diagnostics);

/// Single-line description "error /root/...: message".
[[nodiscard]] std::string to_string(const Diagnostic& diagnostic);

}  // namespace rtt::model
