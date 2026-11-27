#pragma once

/// @file random.hpp
/// Random perturbations of tolerance items (ADR 0036, point 5; #259): a counter-based generator
/// (Philox-4x32-10), so that the value of an item depends only on the seed, the trial and the
/// item's id, not on threads, on the order of the draws or on the other items.
///
/// - Key (k0, k1) = (lower, upper 32 bits of the seed).
/// - Counter (c0, c1, c2, c3) = (trial, lower and upper 32 bits of h, call j) with
///   h = model::tolerance_id_hash(item id), the FNV-1a 64 hash that validate also uses
///   for the collision check (ADR 0035, point 9).
/// - One call gives four words (w0, w1, w2, w3) and so two uniform numbers in [0, 1):
///   u_a = uniform53(w0, w1) and u_b = uniform53(w2, w3).
///
/// The distributions (limit T, ADR 0035, point 2):
/// - uniform: T (2 u_a - 1) from call 0;
/// - endpoint: -T for u_a < 1/2, else +T (call 0);
/// - normal: sigma = T/2, truncated at +-T by rejection. Draw k = 0, 1, ..., 63 is
///   z_k (Box-Muller, call floor(k/2), z1 for even and z2 for odd k); the first sigma z_k with
///   |z_k| <= 2 is taken. After 64 rejected draws (probability about 1e-86) the value is
///   T (2 u_a - 1) from call 32.
///
/// Values are bitwise the same within a build (threads, repetition). Box-Muller uses log, cos
/// and sin, which are not correctly rounded: normal values may differ between platforms in the
/// last bits, the uniform and endpoint values do not.

#include <array>
#include <cstdint>
#include <span>
#include <vector>

#include "rtt/math/philox.hpp"
#include "rtt/model/tolerance.hpp"

namespace rtt::tolerance {

// The generator lives in rtt-math since ADR 0057, point 1 (shared with rtt-nonseq); the names
// stay available here unchanged.
using math::philox4x32_10;
using math::PhiloxCounter;
using math::PhiloxKey;
using math::uniform53;

/// Number of normal draws before the fallback of the truncated normal distribution.
inline constexpr int kMaxNormalDraws = 64;

/// The perturbation of `item` in trial `trial` (0-based) with seed `seed`, in the unit of the
/// item (that of its limit: mm, degree for the kinds ending in Deg, dimensionless otherwise),
/// in [-T, T]. Depends only on (seed, trial, item.id, item.limit, item.distribution).
/// @throws std::invalid_argument if the limit is not finite or <= 0, or the distribution is not
///         a valid ToleranceDistribution
[[nodiscard]] double perturbation(std::uint64_t seed,
                                  std::uint32_t trial,
                                  const model::ToleranceItem& item);

/// perturbation() of every item, in the order of `items`.
/// @throws as perturbation()
[[nodiscard]] std::vector<double> perturbations(std::uint64_t seed,
                                                std::uint32_t trial,
                                                std::span<const model::ToleranceItem> items);

}  // namespace rtt::tolerance
