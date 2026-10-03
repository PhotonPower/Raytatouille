#pragma once

/// @file real.hpp
/// Scalar abstraction for templated physics code (ADR 06).

#include <type_traits>

namespace rtt::math {

/// Marks scalar types that the templated physics libraries (geom, material, coating, polar)
/// accept. Built-in floating point types qualify. Automatic-differentiation types (dual
/// numbers) are enabled later by specialising this trait; no other code has to change.
template <class T>
struct is_real : std::bool_constant<std::is_floating_point_v<T>> {};

template <class T>
inline constexpr bool is_real_v = is_real<T>::value;

/// Concept used to constrain templated physics functions.
template <class T>
concept Real = is_real_v<T>;

}  // namespace rtt::math
