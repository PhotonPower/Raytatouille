#pragma once

#include <nlohmann/json.hpp>
#include <string>

namespace rtt::io::detail {

/// Pretty printer with 2-space indent that keeps arrays of scalars on one line.
[[nodiscard]] std::string format_canonical(const nlohmann::ordered_json& j);

}  // namespace rtt::io::detail
