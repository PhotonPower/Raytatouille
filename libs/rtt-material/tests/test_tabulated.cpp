#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

#include "rtt/material/tabulated.hpp"

using Catch::Matchers::WithinRel;
using rtt::material::IndexSample;
using rtt::material::TabulatedMaterial;
using rtt::material::WavelengthRange;
using rtt::math::Complex;

namespace {

const std::vector<IndexSample> kSamples{{0.4, 1.6, 0.01}, {0.5, 1.55, 0.0}, {0.8, 1.5, 0.002}};

void require_index(const Complex& actual, double n, double kappa) {
  REQUIRE_THAT(actual.real(), WithinRel(n, 1e-12));
  if (kappa == 0.0) {
    REQUIRE(actual.imag() == 0.0);
  } else {
    REQUIRE_THAT(actual.imag(), WithinRel(kappa, 1e-12));
  }
}

}  // namespace

TEST_CASE("tabulated data are exact at the samples", "[tabulated]") {
  const TabulatedMaterial m(kSamples);
  for (const auto& s : kSamples) {
    INFO(s.wavelength_um);
    REQUIRE(m.index(s.wavelength_um, 20.0, 1.0) == Complex(s.n, s.kappa));
  }
}

TEST_CASE("tabulated data are linear in the wavelength between samples", "[tabulated]") {
  const TabulatedMaterial m(kSamples);
  // Midpoint of [0.4, 0.5]: averages. At 0.65, halfway in [0.5, 0.8]: averages too.
  require_index(m.index(0.45, 20.0, 1.0), 1.575, 0.005);
  require_index(m.index(0.65, 20.0, 1.0), 1.525, 0.001);
  // A quarter of the way into [0.5, 0.8]: n = 1.55 - 0.25 * 0.05, kappa = 0.25 * 0.002.
  require_index(m.index(0.575, 20.0, 1.0), 1.5375, 0.0005);
  // Temperature and pressure are ignored.
  REQUIRE(m.index(0.65, -40.0, 0.0) == m.index(0.65, 80.0, 3.0));
}

TEST_CASE("tabulated range and values beyond it", "[tabulated]") {
  const TabulatedMaterial m(kSamples);
  REQUIRE(m.wavelength_range_um() == WavelengthRange{0.4, 0.8});
  REQUIRE(m.wavelength_range_um()->contains(0.4));
  REQUIRE(m.wavelength_range_um()->contains(0.8));
  REQUIRE_FALSE(m.wavelength_range_um()->contains(0.39));
  // Outside the range the nearest end value (rtt-compile rejects such wavelengths).
  REQUIRE(m.index(0.3, 20.0, 1.0) == Complex(1.6, 0.01));
  REQUIRE(m.index(1.2, 20.0, 1.0) == Complex(1.5, 0.002));
}

TEST_CASE("tabulated data are validated", "[tabulated]") {
  const auto bad = [](std::vector<IndexSample> s) {
    REQUIRE_THROWS_AS(TabulatedMaterial(std::move(s)), std::invalid_argument);
  };
  bad({});
  bad({{0.5, 1.5, 0.0}});                        // fewer than 2 samples
  bad({{0.5, 1.5, 0.0}, {0.5, 1.6, 0.0}});       // not strictly increasing
  bad({{0.6, 1.5, 0.0}, {0.5, 1.6, 0.0}});       // decreasing
  bad({{0.0, 1.5, 0.0}, {0.5, 1.6, 0.0}});       // wavelength <= 0
  bad({{0.4, 0.0, 0.0}, {0.5, 1.6, 0.0}});       // n <= 0
  bad({{0.4, 1.5, -0.1}, {0.5, 1.6, 0.0}});      // kappa < 0
  bad({{0.4, NAN, 0.0}, {0.5, 1.6, 0.0}});       // not finite
  bad({{0.4, 1.5, 0.0}, {INFINITY, 1.6, 0.0}});  // not finite
  REQUIRE_NOTHROW(TabulatedMaterial({{0.4, 1.5, 0.0}, {0.5, 1.6, 0.0}}));
}
