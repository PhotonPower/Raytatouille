#include "rtt/math/philox.hpp"

#include <cmath>

namespace rtt::math {

std::array<std::uint32_t, 4> philox4x32_10(const PhiloxCounter& counter,
                                           const PhiloxKey& key) noexcept {
  constexpr std::uint64_t kM0 = 0xD2511F53U;
  constexpr std::uint64_t kM1 = 0xCD9E8D57U;
  constexpr std::uint32_t kWeyl0 = 0x9E3779B9U;
  constexpr std::uint32_t kWeyl1 = 0xBB67AE85U;
  std::array<std::uint32_t, 4> c = counter;
  std::uint32_t k0 = key[0];
  std::uint32_t k1 = key[1];
  for (int round = 0; round < 10; ++round) {
    const std::uint64_t p0 = kM0 * c[0];
    const std::uint64_t p1 = kM1 * c[2];
    c = {static_cast<std::uint32_t>(p1 >> 32) ^ c[1] ^ k0, static_cast<std::uint32_t>(p1),
         static_cast<std::uint32_t>(p0 >> 32) ^ c[3] ^ k1, static_cast<std::uint32_t>(p0)};
    k0 += kWeyl0;
    k1 += kWeyl1;
  }
  return c;
}

double uniform53(std::uint32_t w0, std::uint32_t w1) noexcept {
  const std::uint64_t bits = (std::uint64_t{w0} << 21) | (w1 >> 11);
  return std::ldexp(static_cast<double>(bits), -53);
}

}  // namespace rtt::math
