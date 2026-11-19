// Zernike polynomials in rtt-math (ADR 0051, points 2 and 3): the Fringe set against Niu/Tian
// Table 4, the Noll factor, and the radial derivative.

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <stdexcept>

#include "rtt/math/zernike.hpp"

using rtt::math::FringeIndex;

TEST_CASE("Fringe order follows Niu/Tian Table 4 (ADR 0051)", "[zernike]") {
  // (n, signed m) of Table 4, j = 0 ... 36; file term k = j + 1.
  struct Row {
    int n;
    int m;
  };
  constexpr std::array<Row, 37> kTable4 = {{
      {0, 0},  {1, 1},  {1, -1}, {2, 0},  {2, 2},  {2, -2}, {3, 1},  {3, -1}, {4, 0},  {3, 3},
      {3, -3}, {4, 2},  {4, -2}, {5, 1},  {5, -1}, {6, 0},  {4, 4},  {4, -4}, {5, 3},  {5, -3},
      {6, 2},  {6, -2}, {7, 1},  {7, -1}, {8, 0},  {5, 5},  {5, -5}, {6, 4},  {6, -4}, {7, 3},
      {7, -3}, {8, 2},  {8, -2}, {9, 1},  {9, -1}, {10, 0}, {12, 0},
  }};
  for (int k = 1; k <= 37; ++k) {
    INFO("Fringe term " << k);
    const FringeIndex f = rtt::math::fringe_to_nm(k);
    const Row& r = kTable4[static_cast<std::size_t>(k - 1)];
    CHECK(f.n == r.n);
    CHECK(f.m == std::abs(r.m));
    CHECK(f.sine == (r.m < 0));
  }
  CHECK_THROWS_AS(rtt::math::fringe_to_nm(0), std::invalid_argument);
  CHECK_THROWS_AS(rtt::math::fringe_to_nm(38), std::invalid_argument);
}

TEST_CASE("Fringe polynomials are 1 at the rim and Noll = factor * Fringe (ADR 0051)",
          "[zernike]") {
  // Table 4, j = 8: 6 rho^4 - 6 rho^2 + 1; j = 36: 924 rho^12 - ... + 1, both 1 at rho = 1.
  CHECK(rtt::math::zernike_radial(4, 0, 1.0) == 1.0);
  CHECK(rtt::math::zernike_radial(12, 0, 1.0) == 1.0);
  CHECK(rtt::math::zernike_radial(4, 0, 0.5) == 6.0 * 0.0625 - 6.0 * 0.25 + 1.0);
  // Noll j = 11 is (4, 0): sqrt(5) R_4^0; j = 6 is (2, 2) with the cosine.
  const double rho = 0.7;
  const double theta = 0.3;
  CHECK(rtt::math::zernike(11, rho, theta) ==
        rtt::math::zernike_noll_factor(4, 0) * rtt::math::zernike_radial(4, 0, rho));
  CHECK(std::abs(rtt::math::zernike(6, rho, theta) -
                 std::sqrt(6.0) * rtt::math::zernike_radial(2, 2, rho) * std::cos(2.0 * theta)) <
        1e-15);
}

TEST_CASE("radial derivative and R/rho agree with the polynomial (ADR 0051)", "[zernike]") {
  // Independent: R = sum_s c_s rho^(n - 2s), dR/drho = sum_s c_s (n - 2s) rho^(n - 2s - 1), as
  // power sums in long double (Niu/Tian Eq. (4)). On macOS arm64 long double is double, so the
  // reference rounds like the tested code. Bound: the Horner forms round by at most
  // gamma_(3K+m+1) <= gamma_31 <= 31u of the scale sum_s |c_s| rho^(n-2s) (ADR 0032, point 6; with
  // the derivative one more product), the power sums by at most about (K + 1 + n) u of the same
  // scale (pow is faithful); together below 64 u, and 64 * 2^-52 = 128 u leaves a margin of 2.
  for (int n = 0; n <= rtt::math::kMaxZernikeDegree; ++n) {
    for (int m = n % 2; m <= n; m += 2) {
      const auto c = rtt::math::zernike_radial_coefficients(n, m);
      for (const double rho : {0.0, 1e-9, 0.13, 0.5, 0.97, 1.0, 1.3}) {
        INFO("n " << n << ", m " << m << ", rho " << rho);
        long double r = 0.0L;
        long double d = 0.0L;
        long double scale_r = 0.0L;
        long double scale_d = 0.0L;
        for (std::size_t s = 0; s < c.size(); ++s) {
          const int p = n - 2 * static_cast<int>(s);
          const auto cs = static_cast<long double>(c[s]);
          r += cs * std::pow(static_cast<long double>(rho), p);
          scale_r += std::abs(cs) * std::pow(static_cast<long double>(rho), p);
          if (p >= 1) {
            d += cs * p * std::pow(static_cast<long double>(rho), p - 1);
            scale_d += std::abs(cs) * p * std::pow(static_cast<long double>(rho), p - 1);
          }
        }
        const auto [over, deriv] = rtt::math::zernike_radial_over_rho_and_derivative(n, m, rho);
        const long double eps = 0x1.0p-52L;
        CHECK(std::abs(static_cast<long double>(deriv) - d) <= 64.0L * eps * (scale_d + 1e-300L));
        CHECK(std::abs(static_cast<long double>(rtt::math::zernike_radial(n, m, rho)) - r) <=
              64.0L * eps * (scale_r + 1e-300L));
        if (m >= 1 && rho > 0.0) {
          CHECK(std::abs(static_cast<long double>(over) * static_cast<long double>(rho) - r) <=
                64.0L * eps * (scale_r + 1e-300L));
        }
      }
    }
  }
  // Near the vertex R_1^1 / rho = 1 exactly and R_3^1 / rho = 3 rho^2 - 2 without cancellation.
  CHECK(rtt::math::zernike_radial_over_rho_and_derivative(1, 1, 1e-300).first == 1.0);
  CHECK(rtt::math::zernike_radial_over_rho_and_derivative(3, 1, 0.0).first == -2.0);
  CHECK(rtt::math::zernike_radial_over_rho_and_derivative(3, 1, 0.0).second == -2.0);
}
