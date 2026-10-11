// The ideal lens in rtt-paraxial (#178, ADR 0031, point 6; tests of point 9 for PR 3): a y-nu
// step with phi = |n| / f, so rear_focal_length = f and efl = f / |n|; Seidel contributions 0;
// the ideal cylinder lens is not rotationally symmetric (ParaxialError).
//
// Setup of tests/reference/r2/ideal_lens.rtt.json: f = 50 mm at z = 0, object at 75 mm, stop at
// z = -1. Tolerances relative 1e-14: a handful of y-nu steps round by a few u = 1.1e-16.

#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <optional>
#include <string>
#include <variant>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
#include "rtt/paraxial/paraxial.hpp"
#include "rtt/paraxial/seidel.hpp"

using rtt::compile::CompiledSystem;
using rtt::model::Element;
using rtt::model::IdealLens;
using rtt::model::Param;
using rtt::model::System;
using rtt::paraxial::FirstOrder;

namespace {

System reference() {
  return rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/r2/ideal_lens.rtt.json");
}

Element& child(System& s, std::size_t i) {
  return std::get<Element>(s.root.children[i].value);
}

CompiledSystem compiled(const System& s) {
  return rtt::compile::compile(s, rtt::material::MaterialLibrary{});
}

FirstOrder first_order(const CompiledSystem& cs, const char* path) {
  return rtt::paraxial::first_order(cs, *cs.find_path(path), cs.reference_wavelength());
}

bool close(double a, double b, double relative = 1e-14) {
  return std::abs(a - b) <= relative * std::abs(b);
}

}  // namespace

TEST_CASE("ideal lens in VACUUM: efl = f, bfl = f, principal planes in the lens",
          "[ideal][paraxial]") {
  System s = reference();
  s.environment.medium = "VACUUM";
  const FirstOrder fo = first_order(compiled(s), "ideal lens");
  REQUIRE(fo.efl.has_value());
  CHECK(close(*fo.efl, 50.0));
  CHECK(close(fo.rear_focal_length.value(), 50.0));
  CHECK(close(fo.bfl.value(), 50.0));
  CHECK(std::abs(fo.front_principal_z.value()) <= 1e-12);
  CHECK(std::abs(fo.rear_principal_z.value()) <= 1e-12);
  // Object at 75 mm: image at s' = 1 / (1/50 - 1/75) = 150 behind the lens, m = s'/s = -2.
  CHECK(close(fo.image_z.value(), 150.0));
  CHECK(close(fo.lateral_magnification.value(), -2.0));
}

TEST_CASE("ideal lens in AIR: rear_focal_length = f, efl = f / n_air", "[ideal][paraxial]") {
  const CompiledSystem cs = compiled(reference());
  const double n_air = cs.media()[0].index[cs.reference_wavelength()].real();
  REQUIRE(n_air > 1.0);
  const FirstOrder fo = first_order(cs, "ideal lens");
  CHECK(close(fo.rear_focal_length.value(), 50.0));
  CHECK(close(fo.efl.value(), 50.0 / n_air));
  CHECK(close(fo.bfl.value(), 50.0));
}

TEST_CASE("two ideal lenses: efl = f1 f2 / (f1 + f2 - e)", "[ideal][paraxial]") {
  // Two y-nu steps of the y-nu trace (apply() in paraxial.cpp) give Phi = phi1 + phi2 - e phi1
  // phi2: efl = 1 / Phi = f1 f2 / (f1 + f2 - e). f1 = 50 at z = 0, f2 = 80 at z = 20 in VACUUM.
  System s = reference();
  s.environment.medium = "VACUUM";
  Element& second = child(s, 2);
  second.surfaces[0].interaction = IdealLens{Param(80.0), std::nullopt};
  second.pose.position[2] = Param(20.0);
  using rtt::model::EventKind;
  using rtt::model::SurfaceId;
  s.paths.push_back({"two lenses",
                     false,
                     {{SurfaceId{"STO"}, EventKind::Transmit, 0},
                      {SurfaceId{"IL"}, EventKind::Transmit, 0},
                      {SurfaceId{"CL"}, EventKind::Transmit, 0},
                      {SurfaceId{"IMG"}, EventKind::Transmit, 0}}});
  const FirstOrder fo = first_order(compiled(s), "two lenses");
  CHECK(close(fo.efl.value(), 50.0 * 80.0 / (50.0 + 80.0 - 20.0)));
}

TEST_CASE("ideal lens: Seidel contributions exactly 0", "[ideal][paraxial]") {
  const CompiledSystem cs = compiled(reference());
  const rtt::paraxial::Seidel sd = rtt::paraxial::seidel(cs, *cs.find_path("ideal lens"),
                                                         cs.reference_wavelength(), std::nullopt);
  bool found = false;
  for (const rtt::paraxial::SeidelSurface& s : sd.surfaces) {
    if (cs.surfaces()[s.surface].id.str() != "IL") continue;
    found = true;
    CHECK(s.terms.s1 == 0.0);
    CHECK(s.terms.s2 == 0.0);
    CHECK(s.terms.s3 == 0.0);
    CHECK(s.terms.s4 == 0.0);
    CHECK(s.terms.s5 == 0.0);
    CHECK(s.terms.c_l == 0.0);
    CHECK(s.terms.c_t == 0.0);
    CHECK(s.y != 0.0);  // content guard: the marginal ray passes the lens off the axis
  }
  CHECK(found);
}
