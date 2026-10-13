#include "rtt/tolerance/random.hpp"

#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>

#include "random_detail.hpp"

namespace rtt::tolerance {

namespace {

/// Words of the Philox call `call` for (seed, trial, id hash) (ADR 0036, point 5).
std::array<std::uint32_t, 4> draw_call(const PhiloxKey& key,
                                       std::uint32_t trial,
                                       std::uint64_t hash,
                                       std::uint32_t call) {
  return philox4x32_10(
      {trial, static_cast<std::uint32_t>(hash), static_cast<std::uint32_t>(hash >> 32), call}, key);
}

}  // namespace

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

std::uint64_t fnv1a64(std::string_view bytes) noexcept {
  constexpr std::uint64_t kOffsetBasis = 0xCBF29CE484222325ULL;
  constexpr std::uint64_t kPrime = 0x00000100000001B3ULL;
  std::uint64_t hash = kOffsetBasis;
  for (const char ch : bytes) {
    hash ^= static_cast<unsigned char>(ch);
    hash *= kPrime;  // modulo 2^64 by unsigned overflow
  }
  return hash;
}

double uniform53(std::uint32_t w0, std::uint32_t w1) noexcept {
  const std::uint64_t bits = (std::uint64_t{w0} << 21) | (w1 >> 11);
  return std::ldexp(static_cast<double>(bits), -53);
}

double perturbation(std::uint64_t seed, std::uint32_t trial, const model::ToleranceItem& item) {
  const double limit = item.limit;
  if (!std::isfinite(limit) || !(limit > 0.0)) {
    throw std::invalid_argument("tolerance: the limit of item '" + item.id +
                                "' must be finite and > 0");
  }
  const PhiloxKey key{static_cast<std::uint32_t>(seed), static_cast<std::uint32_t>(seed >> 32)};
  const std::uint64_t hash = fnv1a64(item.id);
  const auto uniform_a = [&](std::uint32_t call) {
    const auto w = draw_call(key, trial, hash, call);
    return uniform53(w[0], w[1]);
  };
  switch (item.distribution) {
    case model::ToleranceDistribution::Uniform:
      return limit * (2.0 * uniform_a(0) - 1.0);
    case model::ToleranceDistribution::Endpoint:
      return uniform_a(0) < 0.5 ? -limit : limit;
    case model::ToleranceDistribution::Normal: {
      // Draw k: call floor(k/2), z1 for even and z2 for odd k. truncated_normal asks for
      // k = 0, 1, 2, ... in order, so the pair of a call is computed once, at its even draw.
      std::pair<double, double> pair{};
      const auto normal = [&](int k) {
        if (k % 2 == 0) {
          const auto w = draw_call(key, trial, hash, static_cast<std::uint32_t>(k / 2));
          pair = detail::box_muller(uniform53(w[0], w[1]), uniform53(w[2], w[3]));
          return pair.first;
        }
        return pair.second;
      };
      // After kMaxNormalDraws rejected draws: u_a of call 32 (ADR 0036, point 5).
      const auto fallback = [&] {
        return uniform_a(static_cast<std::uint32_t>(kMaxNormalDraws / 2));
      };
      return detail::truncated_normal(limit, kMaxNormalDraws, normal, fallback);
    }
  }
  throw std::invalid_argument("tolerance: item '" + item.id + "' has an unknown distribution");
}

std::vector<double> perturbations(std::uint64_t seed,
                                  std::uint32_t trial,
                                  std::span<const model::ToleranceItem> items) {
  std::vector<double> values;
  values.reserve(items.size());
  for (const model::ToleranceItem& i : items) values.push_back(perturbation(seed, trial, i));
  return values;
}

}  // namespace rtt::tolerance
