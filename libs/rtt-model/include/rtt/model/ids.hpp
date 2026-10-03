#pragma once

/// @file ids.hpp
/// Stable identifiers. Paths, tolerances and analyses refer to surfaces only via SurfaceId.

#include <compare>
#include <functional>
#include <string>
#include <utility>

namespace rtt::model {

/// Stable, unique, human-readable identifier of a surface (e.g. "L1.S1").
class SurfaceId {
 public:
  SurfaceId() = default;
  explicit SurfaceId(std::string value) : value_(std::move(value)) {}

  [[nodiscard]] const std::string& str() const noexcept { return value_; }
  [[nodiscard]] bool empty() const noexcept { return value_.empty(); }

  auto operator<=>(const SurfaceId&) const = default;

 private:
  std::string value_;
};

}  // namespace rtt::model

template <>
struct std::hash<rtt::model::SurfaceId> {
  std::size_t operator()(const rtt::model::SurfaceId& id) const noexcept {
    return std::hash<std::string>{}(id.str());
  }
};
