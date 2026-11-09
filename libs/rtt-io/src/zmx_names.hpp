#pragma once

/// @file zmx_names.hpp
/// Readings of ZMX names shared by the import (ADR 0044) and the export (ADR 0045).

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <string>

namespace rtt::io::detail {

/// `s` in upper case (ASCII).
inline std::string upper(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(),
                 [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
  return s;
}

/// MIL number as a glass name: six digits or the form 517.642 (optiland, docs/quellen.md). The
/// import reads such a name as a model glass (ADR 0044, point 7).
inline bool mil_number(const std::string& name) {
  const auto digits = [&](std::size_t from, std::size_t n) {
    for (std::size_t i = from; i < from + n; ++i) {
      if (std::isdigit(static_cast<unsigned char>(name[i])) == 0) return false;
    }
    return true;
  };
  if (name.size() == 6) return digits(0, 6);
  if (name.size() == 7) return name[3] == '.' && digits(0, 3) && digits(4, 3);
  return false;
}

}  // namespace rtt::io::detail
