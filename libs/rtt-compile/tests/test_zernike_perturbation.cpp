// The Zernike error of compile::Perturbation (#316, ADR 0051, point 7; ADR 0036, addendum #316):
// appended after the surface's own terms, normalized to the semi-diameter of the aperture (here
// the half diagonal of a rectangle), equal to the term written by hand; all zero is no change;
// the API errors.
//
// GCC 13 at -O3 (CI ci-linux-gcc) reports -Wmaybe-uninitialized in std::string for the optional
// row name of a Param moved into the term variant (the other alternatives inlined on unreachable
// branches); GCC 16 and clang do not. A false positive, suppressed for this test translation unit
// only and only for GCC, as in libs/rtt-trace/tests/test_diffraction.cpp.
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif

#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "rtt/coating/catalog.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/compile/perturbation.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"

using rtt::compile::CompiledSystem;
using rtt::compile::Perturbation;
using rtt::model::Element;
using rtt::model::Param;
using rtt::model::SurfaceId;
using rtt::model::System;

namespace {

/// m1/singlet_const with a rectangular aperture 4 x 3 (half widths) on L1.S1 and, with `own`, a
/// Noll Z_11 term of its own there (radius 9).
System singlet(bool own) {
  System s = rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m1/singlet_const.rtt.json");
  auto& surface = std::get<Element>(s.root.children[1].value).surfaces[0];
  surface.aperture = rtt::model::RectangularAperture{Param(4.0), Param(3.0)};
  if (own) {
    std::vector<Param> a(11, Param(0.0));
    a[10] = Param(2e-4);
    surface.shape.terms.emplace_back(std::in_place_type<rtt::model::ZernikeSag>, Param(9.0),
                                     std::move(a));
  }
  return s;
}

CompiledSystem compiled(const System& s, const Perturbation& p = {}) {
  const rtt::material::MaterialLibrary lib;
  return rtt::compile::compile(s, lib, rtt::coating::CoatingLibrary{}, 0, p);
}

const rtt::compile::CompiledShape& shape(const CompiledSystem& cs) {
  const std::optional<std::uint32_t> i = cs.find_surface(SurfaceId("L1.S1"));
  REQUIRE(i.has_value());
  return cs.surfaces()[i.value_or(0)].shape;
}

}  // namespace

TEST_CASE("a Zernike error goes after the surface's own terms (ADR 0051, point 7)",
          "[compile][zernike]") {
  Perturbation p;
  p.zernike[SurfaceId("L1.S1")] = {0.0, 0.0, 0.0, 0.0, 0.0, 3e-4, 0.0, -1e-4};  // Z_6, Z_8
  const CompiledSystem trial = compiled(singlet(true), p);
  // By hand: the same term after the own one, radius = half diagonal hypot(4, 3) = 5.
  System hand = singlet(true);
  std::vector<Param> a;
  for (const double v : p.zernike.at(SurfaceId("L1.S1"))) a.emplace_back(v);
  std::get<Element>(hand.root.children[1].value)
      .surfaces[0]
      .shape.terms.emplace_back(std::in_place_type<rtt::model::ZernikeSag>, Param(5.0),
                                std::move(a));
  const CompiledSystem by_hand = compiled(hand);
  const auto* t = std::get_if<rtt::geom::Stacked<double>>(&shape(trial));
  const auto* h = std::get_if<rtt::geom::Stacked<double>>(&shape(by_hand));
  REQUIRE(t != nullptr);
  REQUIRE(h != nullptr);
  REQUIRE(t->terms().size() == 2);
  const auto& added = std::get<rtt::geom::ZernikeTerm<double>>(t->terms()[1]);
  CHECK(added.normalization_radius() == std::hypot(4.0, 3.0));
  for (const auto& [x, y] : {std::pair{0.0, 0.0}, std::pair{1.5, -2.0}, std::pair{-3.9, 2.9}}) {
    CHECK(t->sag(x, y) == h->sag(x, y));
    CHECK(t->grad(x, y) == h->grad(x, y));
  }
  const CompiledSystem plain = compiled(singlet(true));
  const auto* nominal = std::get_if<rtt::geom::Stacked<double>>(&shape(plain));
  REQUIRE(nominal != nullptr);
  CHECK(t->sag(1.5, -2.0) != nominal->sag(1.5, -2.0));  // content guard
}

TEST_CASE("an all-zero Zernike error changes nothing", "[compile][zernike]") {
  Perturbation p;
  p.zernike[SurfaceId("L1.S1")] = {0.0, 0.0, 0.0, 0.0};
  CHECK(std::holds_alternative<rtt::geom::Conic<double>>(shape(compiled(singlet(false), p))));
}

TEST_CASE("Zernike error: API errors of compile", "[compile][zernike]") {
  Perturbation p;
  p.zernike[SurfaceId("NOPE")] = {0.0, 1e-4};
  CHECK_THROWS_AS(compiled(singlet(false), p), std::invalid_argument);
  p.zernike.clear();
  p.zernike[SurfaceId("L1.S1")] = std::vector<double>(232, 0.0);
  CHECK_THROWS_AS(compiled(singlet(false), p), std::invalid_argument);
  p.zernike[SurfaceId("L1.S1")] = {0.0, std::numeric_limits<double>::quiet_NaN()};
  CHECK_THROWS_AS(compiled(singlet(false), p), std::invalid_argument);
  p.zernike[SurfaceId("IMG")] = {0.0, 1e-4};  // the detector surface has no aperture
  p.zernike.erase(SurfaceId("L1.S1"));
  CHECK_THROWS_AS(compiled(singlet(false), p), std::invalid_argument);
}
