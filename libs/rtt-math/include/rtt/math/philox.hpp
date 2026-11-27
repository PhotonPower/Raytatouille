#pragma once

/// @file philox.hpp
/// Counter-based random numbers (ADR 0036, point 5; moved here from rtt-tolerance by ADR 0057,
/// point 1): Philox-4x32-10 maps a counter and a key to four 32-bit words without state, so a
/// value depends only on (counter, key), not on threads or on the order of the draws.
/// uniform53 turns two words into a double in [0, 1). Values are bitwise the same on every
/// platform (integer arithmetic and one exact scaling).

#include <array>
#include <cstdint>

namespace rtt::math {

/// Counter of Philox-4x32 (four 32-bit words).
using PhiloxCounter = std::array<std::uint32_t, 4>;
/// Key of Philox-4x32 (two 32-bit words).
using PhiloxKey = std::array<std::uint32_t, 2>;

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

}  // namespace rtt::math
