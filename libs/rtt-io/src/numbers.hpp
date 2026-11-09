#pragma once

/// @file numbers.hpp
/// Numbers in the text formats of rtt-io (ZMX import and export, surface table; ADR 0045,
/// point 4): written in the shortest form that reads back to the same double, read without
/// locale.

#include <array>
#include <charconv>
#include <cmath>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>

namespace rtt::io::detail {

/// The shortest text that reads back to the same double ([charconv.to.chars], paragraphs 2, 7
/// and 8; ADR 0045, F3).
inline std::string format_number(double v) {
  std::array<char, 32> buf{};
  const auto [end, ec] = std::to_chars(buf.data(), buf.data() + buf.size(), v);
  if (ec != std::errc()) throw std::logic_error("to_chars failed");
  return {buf.data(), end};
}

/// The finite double that the whole of `text` spells (std::from_chars, a leading '+' allowed);
/// none for anything else, also for blanks around it, inf and nan.
inline std::optional<double> parse_number(std::string_view text) {
  const char* first = text.data();
  const char* last = text.data() + text.size();
  if (first != last && *first == '+') ++first;
  double v = 0.0;
  const auto [ptr, ec] = std::from_chars(first, last, v);
  if (first == last || ec != std::errc() || ptr != last || !std::isfinite(v)) return std::nullopt;
  return v;
}

}  // namespace rtt::io::detail
