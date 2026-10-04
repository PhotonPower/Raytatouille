#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <complex>
#include <cstddef>
#include <numbers>
#include <stdexcept>
#include <utility>
#include <vector>

#include "rtt/coating/ideal.hpp"
#include "rtt/coating/tabulated.hpp"
#include "rtt/coating/thickness.hpp"
#include "rtt/coating/transfer_matrix.hpp"

using Catch::Matchers::WithinAbs;
using Complex = std::complex<double>;
using rtt::coating::Amplitudes;
using rtt::coating::IdealCoating;
using rtt::coating::Powers;
using rtt::coating::TabulatedCoating;

namespace {

constexpr double kPi = std::numbers::pi;

void require_close(Complex actual, Complex expected, double tol) {
  REQUIRE_THAT(actual.real(), WithinAbs(expected.real(), tol));
  REQUIRE_THAT(actual.imag(), WithinAbs(expected.imag(), tol));
}

Amplitudes<double> constant(Complex v) {
  return {v, v, v, v};
}

}  // namespace

TEST_CASE("ideal mirror: r_s = -1, r_p = +1, no transmission", "[coating][ideal]") {
  // Decided for #56/#58: an ideal mirror without retardance is r_s = -1, r_p = +1 in the
  // basis [s, p, k] (Convention A), the limit of a perfect conductor (test_transfer_matrix).
  const IdealCoating mirror = IdealCoating::mirror();
  for (const double xi : {0.0, 0.5, 0.9}) {
    const Amplitudes<double> a = mirror.amplitudes(1.0, 1.52, xi);
    REQUIRE(a.rs == Complex(-1.0));
    REQUIRE(a.rp == Complex(1.0));
    REQUIRE(a.ts == Complex(0.0));
    REQUIRE(a.tp == Complex(0.0));
  }
}

TEST_CASE("ideal anti-reflection and beam splitter conserve power", "[coating][ideal]") {
  for (const double reflectance : {0.0, 0.3, 1.0}) {
    const IdealCoating coating(reflectance);
    for (const double theta : {0.0, 30.0, 60.0}) {
      const double xi = std::sin(theta * kPi / 180.0);
      const Amplitudes<double> a = coating.amplitudes(1.0, 1.52, xi);
      const Powers<double> p = rtt::coating::stack_powers(a, Complex(1.0), Complex(1.52), xi);
      INFO("R " << reflectance << ", theta " << theta);
      REQUIRE_THAT(p.reflectance_s, WithinAbs(reflectance, 1e-15));
      REQUIRE_THAT(p.reflectance_p, WithinAbs(reflectance, 1e-15));
      REQUIRE_THAT(p.transmittance_s, WithinAbs(1.0 - reflectance, 1e-15));
      REQUIRE_THAT(p.transmittance_p, WithinAbs(1.0 - reflectance, 1e-15));
      REQUIRE(a.ts.imag() == 0.0);  // no phase on transmission
      REQUIRE(a.rs == -a.rp);       // no retardance on reflection
    }
  }
}

TEST_CASE("ideal coating under total internal reflection reflects like the bare interface",
          "[coating][ideal]") {
  const double xi = 1.52 * std::sin(60.0 * kPi / 180.0);  // glass -> air beyond critical
  const Amplitudes<double> a = IdealCoating::anti_reflection().amplitudes(1.52, 1.0, xi);
  const Amplitudes<double> bare = rtt::coating::interface_amplitudes<double>(1.52, 1.0, xi);
  REQUIRE(a.rs == bare.rs);
  REQUIRE(a.rp == bare.rp);
  REQUIRE_THAT(std::abs(a.rs), WithinAbs(1.0, 1e-15));
}

TEST_CASE("ideal coating rejects a reflectance outside [0, 1]", "[coating][ideal]") {
  REQUIRE_THROWS_AS(IdealCoating(-0.1), std::invalid_argument);
  REQUIRE_THROWS_AS(IdealCoating(1.1), std::invalid_argument);
  REQUIRE_THROWS_AS(IdealCoating(std::nan("")), std::invalid_argument);
  REQUIRE(IdealCoating(0.5).reflectance() == 0.5);
}

TEST_CASE("tabulated coating interpolates bilinearly in real and imaginary part",
          "[coating][tabulated]") {
  // Values v(angle, wavelength) = angle + i wavelength are reproduced exactly by bilinear
  // interpolation (linear in both).
  const std::vector<double> angles = {0.0, 0.5, kPi / 2.0};
  const std::vector<double> wavelengths = {0.4, 0.6, 0.8};
  std::vector<Amplitudes<double>> values;
  for (const double a : angles) {
    for (const double w : wavelengths) values.push_back(constant({a, w}));
  }
  const TabulatedCoating table(angles, wavelengths, values);
  REQUIRE(table.min_wavelength_um() == 0.4);
  REQUIRE(table.max_wavelength_um() == 0.8);
  static_assert(noexcept(table.amplitudes(0.0, 0.5)));
  for (const double a : {0.0, 0.25, 0.5, 1.0, kPi / 2.0}) {
    for (const double w : {0.4, 0.45, 0.6, 0.73, 0.8}) {
      const Amplitudes<double> v = table.amplitudes(a, w);
      require_close(v.rs, {a, w}, 1e-15);
      require_close(v.tp, {a, w}, 1e-15);
    }
  }
  // Outside the wavelength range the nearest column is used (the range is checked at
  // compile time, #61); angles are clamped to [0, pi/2].
  require_close(table.amplitudes(0.25, 0.3).rs, {0.25, 0.4}, 1e-15);
  require_close(table.amplitudes(0.25, 1.0).rs, {0.25, 0.8}, 1e-15);
  require_close(table.amplitudes(-0.1, 0.6).rs, {0.0, 0.6}, 1e-15);
}

TEST_CASE("tabulated coating with one wavelength", "[coating][tabulated]") {
  const TabulatedCoating table({0.0, kPi / 2.0}, {0.55}, {constant(0.0), constant({1.0, -1.0})});
  require_close(table.amplitudes(kPi / 4.0, 0.55).rp, {0.5, -0.5}, 1e-15);
  require_close(table.amplitudes(kPi / 4.0, 0.9).rp, {0.5, -0.5}, 1e-15);
}

TEST_CASE("tabulated coating checks its grid at construction", "[coating][tabulated]") {
  const auto make = [](std::vector<double> angles, std::vector<double> wavelengths,
                       std::size_t values) {
    return TabulatedCoating(std::move(angles), std::move(wavelengths),
                            std::vector<Amplitudes<double>>(values));
  };
  REQUIRE_NOTHROW(make({0.0, kPi / 2.0}, {0.5}, 2));
  REQUIRE_THROWS_AS(make({0.0, 1.0}, {0.5}, 2), std::invalid_argument);        // not to pi/2
  REQUIRE_THROWS_AS(make({0.1, kPi / 2.0}, {0.5}, 2), std::invalid_argument);  // not from 0
  REQUIRE_THROWS_AS(make({0.0}, {0.5}, 1), std::invalid_argument);             // one angle
  REQUIRE_THROWS_AS(make({0.0, 1.0, 1.0, kPi / 2.0}, {0.5}, 4), std::invalid_argument);
  REQUIRE_THROWS_AS(make({0.0, kPi / 2.0}, {}, 0), std::invalid_argument);
  REQUIRE_THROWS_AS(make({0.0, kPi / 2.0}, {0.6, 0.5}, 4), std::invalid_argument);
  REQUIRE_THROWS_AS(make({0.0, kPi / 2.0}, {-0.5, 0.5}, 4), std::invalid_argument);
  REQUIRE_THROWS_AS(make({0.0, kPi / 2.0}, {0.5}, 3), std::invalid_argument);  // value count
  REQUIRE_THROWS_AS(
      TabulatedCoating({0.0, kPi / 2.0}, {0.5}, {constant(0.0), constant({std::nan(""), 0.0})}),
      std::invalid_argument);
}

TEST_CASE("thickness: physical or quarter waves", "[coating][thickness]") {
  using rtt::coating::PhysicalThickness;
  using rtt::coating::QuarterWaves;
  REQUIRE(rtt::coating::physical_thickness_um(PhysicalThickness{0.12}, 1.5) == 0.12);
  REQUIRE_THAT(rtt::coating::physical_thickness_um(QuarterWaves{2.0, 0.6}, 1.5),
               WithinAbs(2.0 * 0.6 / (4.0 * 1.5), 1e-16));
  REQUIRE_THROWS_AS(rtt::coating::physical_thickness_um(PhysicalThickness{-0.1}, 1.5),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(rtt::coating::physical_thickness_um(QuarterWaves{-1.0, 0.6}, 1.5),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(rtt::coating::physical_thickness_um(QuarterWaves{1.0, 0.0}, 1.5),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(rtt::coating::physical_thickness_um(QuarterWaves{1.0, 0.6}, 0.0),
                    std::invalid_argument);
}
