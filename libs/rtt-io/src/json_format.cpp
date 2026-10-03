#include "json_format.hpp"

#include <algorithm>
#include <optional>

namespace rtt::io::detail {

namespace {

using Json = nlohmann::ordered_json;

/// Maximum length of an object or array written on a single line.
constexpr std::size_t kMaxInlineLength = 72;

bool is_scalar(const Json& j) {
  return !j.is_object() && !j.is_array();
}

bool is_scalar_array(const Json& j) {
  return j.is_array() && std::all_of(j.begin(), j.end(), is_scalar);
}

/// Single-line rendering for objects whose values are scalars or scalar arrays,
/// e.g. {"um": 0.5876, "reference": true}. Empty optional if not applicable.
std::optional<std::string> inline_object(const Json& j) {
  if (!j.is_object() || j.empty()) return std::nullopt;
  std::string out = "{";
  bool first = true;
  for (const auto& [key, value] : j.items()) {
    if (!is_scalar(value) && !is_scalar_array(value)) return std::nullopt;
    if (!first) out += ", ";
    first = false;
    out += Json(key).dump() + ": ";
    if (value.is_array()) {
      out += "[";
      for (std::size_t i = 0; i < value.size(); ++i) out += (i ? ", " : "") + value[i].dump();
      out += "]";
    } else {
      out += value.dump();
    }
  }
  out += "}";
  if (out.size() > kMaxInlineLength) return std::nullopt;
  return out;
}

void write(const Json& j, int indent, std::string& out) {
  const std::string pad(static_cast<std::size_t>(indent) * 2, ' ');
  const std::string inner(static_cast<std::size_t>(indent + 1) * 2, ' ');
  if (j.is_object()) {
    if (j.empty()) {
      out += "{}";
      return;
    }
    if (indent > 0) {
      if (const auto line = inline_object(j)) {
        out += *line;
        return;
      }
    }
    out += "{\n";
    bool first = true;
    for (const auto& [key, value] : j.items()) {
      if (!first) out += ",\n";
      first = false;
      out += inner;
      out += Json(key).dump();
      out += ": ";
      write(value, indent + 1, out);
    }
    out += "\n" + pad + "}";
  } else if (j.is_array()) {
    if (j.empty()) {
      out += "[]";
    } else if (is_scalar_array(j)) {
      out += "[";
      bool first = true;
      for (const auto& e : j) {
        if (!first) out += ", ";
        first = false;
        out += e.dump();
      }
      out += "]";
    } else {
      out += "[\n";
      bool first = true;
      for (const auto& e : j) {
        if (!first) out += ",\n";
        first = false;
        out += inner;
        write(e, indent + 1, out);
      }
      out += "\n" + pad + "]";
    }
  } else {
    out += j.dump();
  }
}

}  // namespace

std::string format_canonical(const Json& j) {
  std::string out;
  write(j, 0, out);
  out += "\n";
  return out;
}

}  // namespace rtt::io::detail
