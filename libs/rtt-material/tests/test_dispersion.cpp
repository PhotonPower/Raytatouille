#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <stdexcept>

#include "rtt/material/dispersion.hpp"

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using rtt::material::CauchyCoefficients;
using rtt::material::ConradyCoefficients;
using rtt::material::DispersionMaterial;
using rtt::material::HerzbergerCoefficients;
using rtt::material::refractive_index;
using rtt::material::SchottCoefficients;
using rtt::material::Sellmeier1Coefficients;
using rtt::material::Sellmeier2Coefficients;
using rtt::material::Sellmeier3Coefficients;
using rtt::material::Sellmeier4Coefficients;
using rtt::material::Sellmeier5Coefficients;
using rtt::material::WavelengthRange;

namespace {

constexpr double kRel = 1e-12;  // issue #23

// Invented coefficients; the expected values below are evaluated by hand at lambda = 1 um
// (lambda^2 = 1) and lambda = 2 um (lambda^2 = 4), so every power of lambda is exact. Formulas:
// docs/quellen.md (OpticStudio User Guide, "The Glass Dispersion Formulas"; Palmer, Cauchy).

const Sellmeier1Coefficients kSell1{{1.0, 0.5, 0.25}, {0.5, 0.75, 2.0}};
const Sellmeier3Coefficients kSell3{{1.0, 0.5, 0.25, 0.125}, {0.5, 0.75, 2.0, 3.0}};
const Sellmeier5Coefficients kSell5{{1.0, 0.5, 0.25, 0.125, 0.0625}, {0.5, 0.75, 2.0, 3.0, 5.0}};

void require_n(double actual, double expected) {
  REQUIRE_THAT(actual, WithinRel(expected, kRel));
}

}  // namespace

TEST_CASE("Schott: n^2 = a0 + a1 l^2 + a2 l^-2 + a3 l^-4 + a4 l^-6 + a5 l^-8", "[dispersion]") {
  const SchottCoefficients c{{2.0, 0.01, 0.02, 0.003, 0.0004, 0.00005}};
  // l = 1: n^2 = 2 + 0.01 + 0.02 + 0.003 + 0.0004 + 0.00005 = 2.03345
  require_n(refractive_index(c, 1.0), std::sqrt(2.03345));
  // l = 2: n^2 = 2 + 0.04 + 0.005 + 0.0001875 + 0.00000625 + 0.0000001953125 = 2.0451939453125
  require_n(refractive_index(c, 2.0), std::sqrt(2.0451939453125));
}

TEST_CASE("Sellmeier 1, 3, 5: n^2 - 1 = sum K l^2 / (l^2 - L)", "[dispersion]") {
  // l = 1: terms 1/0.5 = 2, 0.5/0.25 = 2, 0.25/(-1) = -0.25, 0.125/(-2) = -0.0625,
  //        0.0625/(-4) = -0.015625
  require_n(refractive_index(kSell1, 1.0), std::sqrt(4.75));
  require_n(refractive_index(kSell3, 1.0), std::sqrt(4.6875));
  require_n(refractive_index(kSell5, 1.0), std::sqrt(4.671875));
  // l = 2: terms 4/3.5 = 8/7, 2/3.25 = 8/13, 1/2, 0.5/1 = 0.5, 0.25/(-1) = -0.25
  const double s1 = 1.0 + 8.0 / 7.0 + 8.0 / 13.0 + 0.5;
  require_n(refractive_index(kSell1, 2.0), std::sqrt(s1));
  require_n(refractive_index(kSell3, 2.0), std::sqrt(s1 + 0.5));
  require_n(refractive_index(kSell5, 2.0), std::sqrt(s1 + 0.5 - 0.25));
}

TEST_CASE("Sellmeier 2: n^2 - 1 = A + B1 l^2/(l^2 - l1^2) + B2/(l^2 - l2^2)", "[dispersion]") {
  const Sellmeier2Coefficients c{0.5, 1.0, 0.5, 0.25, 3.0};  // lambda1^2 = 0.25, lambda2^2 = 9
  // l = 1: 1 + 0.5 + 1/0.75 + 0.25/(-8) = 1.5 + 4/3 - 1/32
  require_n(refractive_index(c, 1.0), std::sqrt(1.5 + 4.0 / 3.0 - 1.0 / 32.0));
  // l = 2: 1 + 0.5 + 4/3.75 + 0.25/(-5) = 1.5 + 16/15 - 0.05
  require_n(refractive_index(c, 2.0), std::sqrt(1.5 + 16.0 / 15.0 - 0.05));
}

TEST_CASE("Sellmeier 4: n^2 = A + B l^2/(l^2 - C) + D l^2/(l^2 - E)", "[dispersion]") {
  const Sellmeier4Coefficients c{1.5, 1.0, 0.25, 0.5, 9.0};
  // l = 1: 1.5 + 1/0.75 + 0.5/(-8) = 1.5 + 4/3 - 1/16
  require_n(refractive_index(c, 1.0), std::sqrt(1.5 + 4.0 / 3.0 - 1.0 / 16.0));
  // l = 2: 1.5 + 4/3.75 + 2/(-5) = 1.5 + 16/15 - 0.4
  require_n(refractive_index(c, 2.0), std::sqrt(1.5 + 16.0 / 15.0 - 0.4));
}

TEST_CASE("Herzberger: n = A + B L + C L^2 + D l^2 + E l^4 + F l^6, L = 1/(l^2 - 0.028)",
          "[dispersion]") {
  const HerzbergerCoefficients c{1.5, 0.01, 0.001, -0.002, 0.0001, -0.00001};
  // l = 1: L = 1/0.972; D l^2 + E l^4 + F l^6 = -0.002 + 0.0001 - 0.00001
  require_n(refractive_index(c, 1.0),
            1.5 + 0.01 / 0.972 + 0.001 / (0.972 * 0.972) - 0.002 + 0.0001 - 0.00001);
  // l = 2: L = 1/3.972; -0.002 * 4 + 0.0001 * 16 - 0.00001 * 64 = -0.008 + 0.0016 - 0.00064
  require_n(refractive_index(c, 2.0),
            1.5 + 0.01 / 3.972 + 0.001 / (3.972 * 3.972) - 0.008 + 0.0016 - 0.00064);
}

TEST_CASE("Conrady: n = n0 + A/l + B/l^3.5", "[dispersion]") {
  const ConradyCoefficients c{1.5, 0.01, 0.001};
  require_n(refractive_index(c, 1.0), 1.511);
  // l = 2: 2^3.5 = 8 sqrt(2)
  require_n(refractive_index(c, 2.0), 1.5 + 0.005 + 0.001 / (8.0 * std::sqrt(2.0)));
}

TEST_CASE("Cauchy: n = A + B/l^2 + C/l^4", "[dispersion]") {
  const CauchyCoefficients c{1.5, 0.004, 0.0002};
  require_n(refractive_index(c, 1.0), 1.5042);
  require_n(refractive_index(c, 2.0), 1.5010125);  // 1.5 + 0.001 + 0.0002/16
}

TEST_CASE("DispersionMaterial evaluates its formula and knows its range", "[dispersion]") {
  const DispersionMaterial m(kSell1, WavelengthRange{0.3, 2.5});
  REQUIRE(m.wavelength_range_um() == WavelengthRange{0.3, 2.5});
  REQUIRE(m.formula() == rtt::material::DispersionFormula{kSell1});
  // Real index, no absorption; temperature and pressure are ignored until #25.
  for (const double t : {-20.0, 20.0, 80.0}) {
    for (const double p : {0.0, 1.0, 2.0}) {
      const auto n = m.index(1.0, t, p);
      require_n(n.real(), std::sqrt(4.75));
      REQUIRE(n.imag() == 0.0);
    }
  }
  require_n(m.index(2.0, 20.0, 1.0).real(), std::sqrt(1.0 + 8.0 / 7.0 + 8.0 / 13.0 + 0.5));
}

TEST_CASE("DispersionMaterial rejects invalid ranges", "[dispersion]") {
  const auto make = [](double lo, double hi) { return DispersionMaterial(kSell1, {lo, hi}); };
  REQUIRE_THROWS_AS(make(0.0, 1.0), std::invalid_argument);
  REQUIRE_THROWS_AS(make(-0.5, 1.0), std::invalid_argument);
  REQUIRE_THROWS_AS(make(1.0, 1.0), std::invalid_argument);
  REQUIRE_THROWS_AS(make(2.0, 1.0), std::invalid_argument);
  REQUIRE_THROWS_AS(make(0.3, INFINITY), std::invalid_argument);
  REQUIRE_THROWS_AS(make(NAN, 1.0), std::invalid_argument);
  REQUIRE_NOTHROW(make(0.3, 2.5));
}

TEST_CASE("N-BK7 from the SCHOTT data sheet", "[dispersion]") {
  // SCHOTT N-BK7 data sheet (as of 01-Dec-2023, docs/quellen.md): Sellmeier constants
  // B1..B3, C1..C3 (C in um^2) and the refractive indices n_d (587.6 nm), n_F (486.1 nm),
  // n_C (656.3 nm), given to 5 decimals. The SCHOTT Sellmeier form is Sellmeier 1 with
  // K = B, L = C.
  const DispersionMaterial nbk7(Sellmeier1Coefficients{{1.039612120, 0.231792344, 1.010469450},
                                                       {0.006000699, 0.0200179144, 103.56065300}},
                                WavelengthRange{0.3, 2.5});
  const auto n = [&](double wl) { return nbk7.index(wl, 20.0, 1.0).real(); };
  // Within the rounding of the data sheet (5th decimal): |n - n_sheet| <= 5e-6.
  REQUIRE_THAT(n(0.5876), WithinAbs(1.51680, 5e-6));
  REQUIRE_THAT(n(0.4861), WithinAbs(1.52238, 5e-6));
  REQUIRE_THAT(n(0.6563), WithinAbs(1.51432, 5e-6));
}
