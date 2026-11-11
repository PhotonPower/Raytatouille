// Guard for the build setting "no floating-point contraction" (ADR 0029, addendum #373).
//
// A fused multiply-add rounds a * b + c once, separate operations round twice. With
// a = 1 + 2^-30, b = 1 - 2^-30, c = -1: a * b = 1 - 2^-60 exactly, rounded to double 1 (the
// spacing below 1 is 2^-53), so a * b + c = 0 with two roundings, but -2^-60 with an FMA.
// Clang fuses by default on targets with an FMA instruction (arm64, -ffp-contract=on), GCC even
// across statements (-ffp-contract=fast); MSVC 2022 does not under /fp:precise. The project
// builds with -ffp-contract=off (cmake/RttTargets.cmake), so this test sees no FMA anywhere.
// The test only discriminates on targets with an FMA instruction (arm64, x86-64 with -mfma);
// on x86-64 without -march, as in the CI, it passes with or without the option.

#include <catch2/catch_test_macros.hpp>
#include <cmath>

TEST_CASE("the build does not contract a * b + c into a fused multiply-add", "[fp]") {
  // volatile: the compiler must not fold the expression at compile time.
  volatile double a = 1.0 + std::ldexp(1.0, -30);
  volatile double b = 1.0 - std::ldexp(1.0, -30);
  volatile double c = -1.0;
  const double x = a;
  const double y = b;
  const double z = c;
  // Precondition: the two roundings and the single rounding differ for these inputs.
  REQUIRE(std::fma(x, y, z) == -std::ldexp(1.0, -60));
  REQUIRE(x * y + z == 0.0);
}
