#pragma once

/// @file load_warning.hpp
/// Warnings of the catalogue readers (ADR 0022, #71): data, not exceptions. Errors stay
/// exceptions (AgfError); a warning reports a documented, narrow exception to the strict
/// reading (ADR 0008 addendum) that the reader applied.

#include <cstddef>
#include <string>

namespace rtt::material {

/// A warning while reading a catalogue: a stable code of the registry rtt/diagnostics/codes.hpp
/// (group agf.*, producer "agf", docs/diagnostics.md) and the place in the file.
struct LoadWarning {
  std::string code;      ///< e.g. "agf.duplicate_glass"
  std::string file;      ///< file name, or the source name of an in-memory catalogue
  std::size_t line = 0;  ///< 1-based line of the file, 0 for the whole file
  std::string message;   ///< human-readable text; may change between versions, the code not
  bool operator==(const LoadWarning&) const = default;
};

}  // namespace rtt::material
