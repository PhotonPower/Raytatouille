#pragma once

/// @file json_io.hpp
/// Reading and writing `.rtt.json` system files (ADR 08, schema/raytatouille.schema.json).
///
/// Parsing is strict: unknown keys, wrong types, wrong units and incompatible schema
/// versions are errors. Writing is canonical: the same model always yields the same bytes,
/// so `to_json(parse(text)) == text` holds for every canonical file. Files of an older supported
/// schema version (0.1) are migrated on reading; writing always uses `model::kSchemaVersion`.

#include <filesystem>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

#include "rtt/model/system.hpp"

namespace rtt::io {

/// Structural error in a system file. `pointer()` is a JSON pointer to the offending value.
class ParseError : public std::runtime_error {
 public:
  ParseError(std::string pointer, const std::string& message)
      : std::runtime_error((pointer.empty() ? std::string("/") : pointer) + ": " + message),
        pointer_(std::move(pointer)) {}

  [[nodiscard]] const std::string& pointer() const noexcept { return pointer_; }

 private:
  std::string pointer_;
};

/// Parses a system from JSON text. Throws ParseError.
[[nodiscard]] model::System parse_system(std::string_view json_text);

/// Reads a system file. Throws ParseError or std::runtime_error on I/O failure.
[[nodiscard]] model::System load_system(const std::filesystem::path& file);

/// Canonical JSON text of a system (2-space indent, LF line ends, trailing newline,
/// default values omitted).
[[nodiscard]] std::string to_json(const model::System& system);

/// Writes the canonical JSON text. Throws std::runtime_error on I/O failure.
void save_system(const model::System& system, const std::filesystem::path& file);

}  // namespace rtt::io
