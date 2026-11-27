// Counter-based generator Philox-4x32-10 and 53-bit uniform numbers in rtt-math (ADR 0057,
// point 1): moved from rtt-tolerance unchanged, so that the non-sequential tracer (rtt-nonseq)
// and the tolerancing draw from the same generator. rtt::tolerance names the same functions.

#include <array>
#include <bit>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstdint>

#include "rtt/math/philox.hpp"

namespace {

bool same(double a, double b) {
  return std::bit_cast<std::uint64_t>(a) == std::bit_cast<std::uint64_t>(b);
}

}  // namespace

TEST_CASE("philox: Philox-4x32-10 gives the known-answer vectors", "[math][philox]") {
  // Known-answer vectors "philox4x32 10" from Random123, tests/kat_vectors
  // (github.com/DEShawResearch/random123): counter, key, expected output (docs/quellen.md,
  // Salmon et al., SC11). License of Random123, file LICENSE, verbatim:
  //
  // Copyright 2010-2012, D. E. Shaw Research.
  // All rights reserved.
  //
  // Redistribution and use in source and binary forms, with or without
  // modification, are permitted provided that the following conditions are
  // met:
  //
  // * Redistributions of source code must retain the above copyright
  //   notice, this list of conditions, and the following disclaimer.
  //
  // * Redistributions in binary form must reproduce the above copyright
  //   notice, this list of conditions, and the following disclaimer in the
  //   documentation and/or other materials provided with the distribution.
  //
  // * Neither the name of D. E. Shaw Research nor the names of its
  //   contributors may be used to endorse or promote products derived from
  //   this software without specific prior written permission.
  //
  // THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  // "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  // LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  // A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
  // OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  // SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
  // LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  // DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  // THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  // (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  // OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
  struct Kat {
    rtt::math::PhiloxCounter counter;
    rtt::math::PhiloxKey key;
    std::array<std::uint32_t, 4> expected;
  };
  const std::array<Kat, 3> kats{{
      {{0x00000000, 0x00000000, 0x00000000, 0x00000000},
       {0x00000000, 0x00000000},
       {0x6627e8d5, 0xe169c58d, 0xbc57ac4c, 0x9b00dbd8}},
      {{0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff},
       {0xffffffff, 0xffffffff},
       {0x408f276d, 0x41c83b0e, 0xa20bc7c6, 0x6d5451fd}},
      {{0x243f6a88, 0x85a308d3, 0x13198a2e, 0x03707344},
       {0xa4093822, 0x299f31d0},
       {0xd16cfe09, 0x94fdcceb, 0x5001e420, 0x24126ea1}},
  }};
  for (const Kat& k : kats) CHECK(rtt::math::philox4x32_10(k.counter, k.key) == k.expected);
}

TEST_CASE("philox: 53-bit uniform numbers from two words", "[math][philox]") {
  using rtt::math::uniform53;
  CHECK(same(uniform53(0, 0), 0.0));
  CHECK(same(uniform53(0xffffffffU, 0xffffffffU), 1.0 - std::ldexp(1.0, -53)));
  // w0 is widened before the shift: its top bit becomes 2^52, i.e. 1/2.
  CHECK(same(uniform53(0x80000000U, 0), 0.5));
  CHECK(same(uniform53(1, 0), std::ldexp(1.0, -32)));
  // Only the upper 21 bits of w1 count.
  CHECK(same(uniform53(0, 0x800U), std::ldexp(1.0, -53)));
  CHECK(same(uniform53(0, 0x7ffU), 0.0));
}
