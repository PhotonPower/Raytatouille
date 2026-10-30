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

#include "rtt/model/tolerance.hpp"

namespace rtt::tolerance {

/// Counter of Philox-4x32 (four 32-bit words).
using PhiloxCounter = std::array<std::uint32_t, 4>;
/// Key of Philox-4x32 (two 32-bit words).
using PhiloxKey = std::array<std::uint32_t, 2>;

/// Number of normal draws before the fallback of the truncated normal distribution.
inline constexpr int kMaxNormalDraws = 64;

/// Philox-4x32 with 10 rounds (Salmon, Moraes, Dror, Shaw, "Parallel Random Numbers: As Easy
/// as 1, 2, 3", SC11, Sec. 4.3: R' = mulhi(R, M) ^ k ^ L; key schedule after Eq. (6)). Each round
/// maps (c0, c1, c2, c3) to (hi(M1 c2) ^ c1 ^ k0, lo(M1 c2), hi(M0 c0) ^ c3 ^ k1, lo(M0 c0)) with
/// M0 = 0xD2511F53 and M1 = 0xCD9E8D57 (64-bit products), and then bumps the key by the Weyl
/// constants 0x9E3779B9 and 0xBB67AE85 (upper halves of the 64-bit constants of Eq. (6)). The
/// paper leaves the word order open; this one reproduces the known-answer vectors of Random123
/// (tests/kat_vectors, ADR 0036, point 5).
[[nodiscard]] std::array<std::uint32_t, 4> philox4x32_10(const PhiloxCounter& counter,
                                                         const PhiloxKey& key) noexcept;

/// Uniform number in [0, 1) with 53 bits from two words: (w0 << 21 | w1 >> 11) * 2^-53, with w0
/// widened to 64 bits before the shift (32 + 21 bits). 0 for w0 = w1 = 0, 1 - 2^-53 for all bits
/// set.
[[nodiscard]] double uniform53(std::uint32_t w0, std::uint32_t w1) noexcept;

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
