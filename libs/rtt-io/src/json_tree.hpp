#pragma once

/// @file json_tree.hpp
/// The JSON trees behind parse_system and to_json, for the edit functions (edit.cpp). Internal to
/// rtt-io.

#include <cstdint>
#include <nlohmann/json.hpp>

#include "rtt/model/system.hpp"

namespace rtt::io::detail {

/// Which values the writer writes.
enum class Form : std::uint8_t {
  Canonical,  ///< defaults omitted, a plain Param as a number (to_json)
  Edit,       ///< every value, every Param as an object (to_edit_json, ADR 0024)
};

/// Reads a system from a parsed JSON value (migration, strict checks).
/// @throws ParseError with the JSON pointer of the offending value
[[nodiscard]] model::System read_system(const nlohmann::json& j);

/// Writes a system as a JSON tree in the given form, keys in file order.
/// @throws std::invalid_argument for a non-finite number, or material and segment_materials both
///         set
[[nodiscard]] nlohmann::ordered_json write_system(const model::System& s, Form form);

/// The same JSON value with sorted object keys (for patching and reading).
[[nodiscard]] nlohmann::json plain(const nlohmann::ordered_json& j);

}  // namespace rtt::io::detail
