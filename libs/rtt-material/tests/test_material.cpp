#include <array>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <thread>
#include <type_traits>
#include <vector>

#include "rtt/material/material.hpp"

using rtt::material::Material;
using rtt::material::MaterialLibrary;
using rtt::material::UnknownMaterial;
using rtt::math::Complex;

namespace {

// Wavelengths (um) and temperatures (degC) at which constant media must not change.
constexpr std::array<double, 3> kWavelengthsUm{0.4861, 0.5876, 10.6};
constexpr std::array<double, 3> kTemperaturesC{-40.0, 20.0, 80.0};

void require_constant_index(const Material& material, const Complex& expected) {
  for (const double wl : kWavelengthsUm) {
    for (const double t : kTemperaturesC) {
      REQUIRE(material.index(wl, t) == expected);
    }
  }
}

}  // namespace

TEST_CASE("UnknownMaterial is a std::runtime_error", "[material]") {
  STATIC_REQUIRE(std::is_base_of_v<std::runtime_error, UnknownMaterial>);
}

TEST_CASE("VACUUM and AIR have index exactly 1 in M1", "[material]") {
  const MaterialLibrary lib;
  require_constant_index(*lib.resolve("VACUUM"), Complex(1.0, 0.0));
  // Ciddor air follows in M2; until then AIR is defined as n = 1 (issue #2).
  require_constant_index(*lib.resolve("AIR"), Complex(1.0, 0.0));
}

TEST_CASE("CONST:<n> gives a constant real index", "[material]") {
  const MaterialLibrary lib;
  // N-BK7 at the d line, n = 1.5168 (Schott data sheet); the parser must round exactly.
  require_constant_index(*lib.resolve("CONST:1.5168"), Complex(1.5168, 0.0));
  require_constant_index(*lib.resolve("CONST:4"), Complex(4.0, 0.0));
  require_constant_index(*lib.resolve("CONST:1.38e0"), Complex(1.38, 0.0));
  // Indices below 1 are valid (e.g. metals, X-rays).
  require_constant_index(*lib.resolve("CONST:0.5"), Complex(0.5, 0.0));
}

TEST_CASE("CONST:<n>,<kappa> gives a constant complex index", "[material]") {
  const MaterialLibrary lib;
  require_constant_index(*lib.resolve("CONST:1.5,0.1"), Complex(1.5, 0.1));
  require_constant_index(*lib.resolve("CONST:1.5,0"), Complex(1.5, 0.0));
  // Metal-like medium with kappa > n.
  require_constant_index(*lib.resolve("CONST:0.2,3.5"), Complex(0.2, 3.5));
  require_constant_index(*lib.resolve("CONST:1.5,1e-7"), Complex(1.5, 1e-7));
}

TEST_CASE("malformed CONST references are rejected", "[material]") {
  const MaterialLibrary lib;
  const std::vector<std::string> bad{
      "CONST:",           // no value
      "CONST:abc",        // not a number
      "CONST:1.5,-0.1",   // negative kappa means gain, not allowed
      "CONST:1.5,",       // missing kappa
      "CONST:,0.1",       // missing n
      "CONST:1.5,0.1,2",  // too many values
      "CONST:1.5x",       // trailing characters
      "CONST: 1.5",       // whitespace
      "CONST:1.5, 0.1",   // whitespace
      "CONST:+1.5",       // leading plus
      "CONST:0",          // n must be positive
      "CONST:-1.5",       // n must be positive
      "CONST:nan",        // not finite
      "CONST:inf",        // not finite
      "CONST:1.5,inf",    // not finite
      "CONST:1e999",      // out of range
  };
  for (const auto& reference : bad) {
    INFO(reference);
    REQUIRE_THROWS_AS(lib.resolve(reference), UnknownMaterial);
  }
}

TEST_CASE("unknown references are rejected with a clear message", "[material]") {
  const MaterialLibrary lib;
  REQUIRE_THROWS_WITH(lib.resolve("SCHOTT:N-BK7"),
                      Catch::Matchers::ContainsSubstring("M2") &&
                          Catch::Matchers::ContainsSubstring("SCHOTT:N-BK7"));
  for (const char* reference : {"", "vacuum", "Air", "VACUUM ", "const:1.5", "BK7"}) {
    INFO(reference);
    REQUIRE_THROWS_AS(lib.resolve(reference), UnknownMaterial);
  }
}

TEST_CASE("resolve returns the same object for identical references", "[material]") {
  const MaterialLibrary lib;
  const auto a = lib.resolve("CONST:1.5168");
  REQUIRE(a != nullptr);
  REQUIRE(lib.resolve("CONST:1.5168") == a);
  REQUIRE(lib.resolve("VACUUM") == lib.resolve("VACUUM"));
  REQUIRE(lib.resolve("AIR") == lib.resolve("AIR"));
  REQUIRE(lib.resolve("CONST:1.5,0.1") == lib.resolve("CONST:1.5,0.1"));
  REQUIRE(lib.resolve("CONST:1.5") != lib.resolve("CONST:1.6"));
}

TEST_CASE("resolve is thread-safe", "[material]") {
  const MaterialLibrary lib;
  const std::vector<std::string> references{"VACUUM", "AIR", "CONST:1.5168", "CONST:0.2,3.5",
                                            "CONST:1.7"};
  constexpr int kThreads = 8;
  constexpr int kRepetitions = 200;

  // Each thread records the first object it got per reference and whether a later call
  // returned a different one. Catch2 assertions are not thread-safe, so checks run afterwards.
  struct Seen {
    std::vector<const Material*> first;
    bool mismatch = false;
  };
  std::vector<Seen> seen(kThreads, Seen{std::vector<const Material*>(references.size()), false});
  {
    std::vector<std::jthread> threads;
    threads.reserve(kThreads);
    for (auto& s : seen) {
      threads.emplace_back([&lib, &references, &s] {
        for (int rep = 0; rep < kRepetitions; ++rep) {
          for (std::size_t r = 0; r < references.size(); ++r) {
            const Material* material = lib.resolve(references[r]).get();
            if (rep == 0) {
              s.first[r] = material;
            } else if (material != s.first[r]) {
              s.mismatch = true;
            }
          }
        }
      });
    }
  }
  for (std::size_t r = 0; r < references.size(); ++r) {
    INFO(references[r]);
    const Material* expected = lib.resolve(references[r]).get();
    for (const auto& s : seen) {
      REQUIRE_FALSE(s.mismatch);
      REQUIRE(s.first[r] == expected);
    }
  }
}

TEST_CASE("kappa = -0 is stored as +0", "[material]") {
  // The sign of a zero imaginary part selects the branch of std::sqrt on its cut (C++
  // [complex.value.ops], C99 Annex G), e.g. cos(theta_t) under TIR. kappa >= 0 must hold
  // including the sign bit.
  const MaterialLibrary lib;
  const Complex index = lib.resolve("CONST:1.5,-0")->index(0.5876, 20.0);
  REQUIRE(index == Complex(1.5, 0.0));
  REQUIRE_FALSE(std::signbit(index.imag()));
  REQUIRE_FALSE(std::signbit(lib.resolve("CONST:1.5,-0.0e3")->index(0.5876, 20.0).imag()));
}

TEST_CASE("a comma separates n and kappa", "[material]") {
  // There is no decimal comma: "CONST:1,5" is n = 1, kappa = 5.
  const MaterialLibrary lib;
  require_constant_index(*lib.resolve("CONST:1,5"), Complex(1.0, 5.0));
}
