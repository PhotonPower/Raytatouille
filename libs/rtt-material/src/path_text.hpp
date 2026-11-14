#pragma once

/// @file path_text.hpp
/// Paths in messages as UTF-8 (#333): path::string() uses the native narrow encoding, which on
/// Windows garbles names outside ASCII; u8string() gives the UTF-8 bytes on every platform.

#include <filesystem>
#include <string>

namespace rtt::material::detail {

/// The UTF-8 bytes of `path` (or of its part, e.g. path.filename()) as a std::string.
inline std::string path_text(const std::filesystem::path& path) {
  const std::u8string text = path.u8string();
  std::string out;
  out.reserve(text.size());
  for (const char8_t c : text) out.push_back(static_cast<char>(c));
  return out;
}

}  // namespace rtt::material::detail
