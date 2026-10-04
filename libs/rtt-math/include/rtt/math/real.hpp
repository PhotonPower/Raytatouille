#pragma once

/// @file real.hpp
/// Scalar abstraction for templated physics code (ADR 0006, ADR 0015).

#include <type_traits>

namespace rtt::math {

/// Marks scalar types that the templated physics libraries (geom, material, coating, polar)
/// accept. Only double qualifies (ADR 0003, ADR 0015): absolute tolerances such as 1e-12 mm
/// are unreachable in float. Automatic-differentiation types (dual numbers over double) are
/// enabled later by specialising this trait; no other code has to change.
template <class T>
struct is_real : std::false_type {};

template <>
struct is_real<double> : std::true_type {};

template <class T>
inline constexpr bool is_real_v = is_real<T>::value;

/// Concept used to constrain templated physics functions.
template <class T>
concept Real = is_real_v<T>;

}  // namespace rtt::math
