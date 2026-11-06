// The anamorphic shapes in the rotationally symmetric paraxial model (#315, ADR 0042, point 1,
// "Scheitelentwicklung", and acceptance 2): a biconic with equal sections and a spherical toroid
// give first_order, prescription and seidel bit for bit as the conic; a real biconic, toroid or
// cylinder has no rotationally symmetric model (ADR 0041, point 2).

#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <variant>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/paraxial/paraxial.hpp"
#include "rtt/paraxial/prescription.hpp"
#include "rtt/paraxial/seidel.hpp"

using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::model::Element;
using rtt::model::Param;
using rtt::model::System;
using rtt::paraxial::ChromaticPair;
using rtt::paraxial::FirstOrder;
using rtt::paraxial::ParaxialError;
using rtt::paraxial::Prescription;
using rtt::paraxial::Seidel;

namespace {

/// m1/singlet_const: L1.S1 (/root/children/1/surfaces/0) is a sphere R = 51.68; the conic
/// constant `k` is set here so that the Seidel cap K enters.
System singlet(double k) {
  System s = rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m1/singlet_const.rtt.json");
  std::get<rtt::model::Conic>(std::get<Element>(s.root.children[1].value).surfaces[0].shape.base)
      .conic = Param(k);
  return s;
}

System with_front(System s, const rtt::model::BaseShape& base) {
  std::get<Element>(s.root.children[1].value).surfaces[0].shape.base = base;
  return s;
}

CompiledSystem compiled(const System& s) {
  const rtt::material::MaterialLibrary lib;
  return rtt::compile::compile(s, lib);
}

void same(const std::optional<rtt::paraxial::Pupil>& a,
          const std::optional<rtt::paraxial::Pupil>& b) {
  REQUIRE(a.has_value() == b.has_value());
  if (!a || !b) return;
  CHECK(a->z == b->z);
  CHECK(a->diameter == b->diameter);
}

void same(const FirstOrder& a, const FirstOrder& b) {
  CHECK(a.object_index == b.object_index);
  CHECK(a.image_index == b.image_index);
  CHECK(a.image_direction == b.image_direction);
  CHECK(a.power == b.power);
  CHECK(a.efl == b.efl);
  CHECK(a.front_focal_length == b.front_focal_length);
  CHECK(a.rear_focal_length == b.rear_focal_length);
  CHECK(a.ffl == b.ffl);
  CHECK(a.bfl == b.bfl);
  CHECK(a.front_focal_z == b.front_focal_z);
  CHECK(a.rear_focal_z == b.rear_focal_z);
  CHECK(a.front_principal_z == b.front_principal_z);
  CHECK(a.rear_principal_z == b.rear_principal_z);
  CHECK(a.image_z == b.image_z);
  CHECK(a.lateral_magnification == b.lateral_magnification);
  CHECK(a.angular_magnification == b.angular_magnification);
  same(a.entrance_pupil, b.entrance_pupil);
  same(a.exit_pupil, b.exit_pupil);
}

void same(const std::optional<rtt::paraxial::PrescriptionRay>& a,
          const std::optional<rtt::paraxial::PrescriptionRay>& b) {
  REQUIRE(a.has_value() == b.has_value());
  if (!a || !b) return;
  CHECK(a->y == b->y);
  CHECK(a->u == b->u);
  CHECK(a->i == b->i);
}

void same(const std::optional<rtt::paraxial::RayStart>& a,
          const std::optional<rtt::paraxial::RayStart>& b) {
  REQUIRE(a.has_value() == b.has_value());
  if (!a || !b) return;
  CHECK(a->z == b->z);
  CHECK(a->y == b->y);
  CHECK(a->u == b->u);
}

void same(const Prescription& a, const Prescription& b) {
  REQUIRE(a.surfaces.size() == b.surfaces.size());
  for (std::size_t i = 0; i < a.surfaces.size(); ++i) {
    INFO("event " << i);
    CHECK(a.surfaces[i].surface == b.surfaces[i].surface);
    CHECK(a.surfaces[i].z == b.surfaces[i].z);
    CHECK(a.surfaces[i].n == b.surfaces[i].n);
    same(a.surfaces[i].marginal, b.surfaces[i].marginal);
    same(a.surfaces[i].chief, b.surfaces[i].chief);
    CHECK(a.surfaces[i].lagrange == b.surfaces[i].lagrange);
  }
  CHECK(a.total_track == b.total_track);
  CHECK(a.object_distance == b.object_distance);
  CHECK(a.paraxial_working_f_number == b.paraxial_working_f_number);
  CHECK(a.paraxial_image_na == b.paraxial_image_na);
  CHECK(a.lagrange_invariant == b.lagrange_invariant);
  same(a.marginal_start, b.marginal_start);
  same(a.chief_start, b.chief_start);
  same(a.first_order, b.first_order);
}

void same(const rtt::paraxial::SeidelTerms& a, const rtt::paraxial::SeidelTerms& b) {
  CHECK(a.s1 == b.s1);
  CHECK(a.s2 == b.s2);
  CHECK(a.s3 == b.s3);
  CHECK(a.s4 == b.s4);
  CHECK(a.s5 == b.s5);
  CHECK(a.c_l == b.c_l);
  CHECK(a.c_t == b.c_t);
}

void same(const Seidel& a, const Seidel& b) {
  REQUIRE(a.surfaces.size() == b.surfaces.size());
  for (std::size_t i = 0; i < a.surfaces.size(); ++i) {
    INFO("event " << i);
    CHECK(a.surfaces[i].surface == b.surfaces[i].surface);
    CHECK(a.surfaces[i].y == b.surfaces[i].y);
    CHECK(a.surfaces[i].y_bar == b.surfaces[i].y_bar);
    CHECK(a.surfaces[i].a == b.surfaces[i].a);
    CHECK(a.surfaces[i].a_bar == b.surfaces[i].a_bar);
    CHECK(a.surfaces[i].lagrange == b.surfaces[i].lagrange);
    same(a.surfaces[i].terms, b.surfaces[i].terms);
  }
  same(a.sum, b.sum);
  CHECK(a.lagrange == b.lagrange);
  same(std::optional{a.marginal}, std::optional{b.marginal});
  same(std::optional{a.chief}, std::optional{b.chief});
  REQUIRE(a.chromatic.has_value() == b.chromatic.has_value());
}

/// Runs f and returns the ParaxialError it throws.
template <typename F>
ParaxialError paraxial_error(F&& f) {
  try {
    f();
  } catch (const ParaxialError& e) {
    return e;
  }
  FAIL("no ParaxialError");
  return ParaxialError("unreachable");
}

}  // namespace

TEST_CASE("a biconic with equal sections is the conic in the paraxial model, bit for bit",
          "[paraxial][anamorphic]") {
  // ADR 0042, acceptance 2: equal parameters give c once for both diagonals and the cap (k, 0),
  // so first_order, prescription and seidel are bitwise those of the conic.
  for (const double k : {0.0, -0.6}) {
    INFO("k = " << k);
    const System conic = singlet(k);
    const System biconic =
        with_front(conic, rtt::model::Biconic{Param(51.68), Param(51.68), Param(k), Param(k)});
    const CompiledSystem a = compiled(conic);
    const CompiledSystem b = compiled(biconic);
    // Content guard: a focal system with a pupil and non-zero aberrations.
    const FirstOrder fa = rtt::paraxial::first_order(a, PathId{0}, 1);
    REQUIRE(fa.efl.has_value());
    REQUIRE(fa.entrance_pupil.has_value());
    const Seidel sa = rtt::paraxial::seidel(a, PathId{0}, 1, ChromaticPair{0, 2});
    REQUIRE(sa.sum.s1 != 0.0);
    REQUIRE(sa.sum.s3 != 0.0);
    REQUIRE(sa.sum.c_l != 0.0);
    for (std::uint16_t w = 0; w < 3; ++w) {
      INFO("wavelength " << w);
      same(rtt::paraxial::first_order(b, PathId{0}, w),
           rtt::paraxial::first_order(a, PathId{0}, w));
      same(rtt::paraxial::prescription(b, PathId{0}, w),
           rtt::paraxial::prescription(a, PathId{0}, w));
      same(rtt::paraxial::seidel(b, PathId{0}, w, ChromaticPair{0, 2}),
           rtt::paraxial::seidel(a, PathId{0}, w, ChromaticPair{0, 2}));
    }
  }
}

TEST_CASE("a spherical toroid is the sphere in the paraxial model, bit for bit",
          "[paraxial][anamorphic]") {
  // ADR 0042, point 1: radius_of_rotation == radius, k = 0 and no coefficients give the
  // Hessian (c, 0, c) with c assigned once and the cap (0, 0).
  const System sphere = singlet(0.0);
  const System toroid =
      with_front(sphere, rtt::model::Toroid{Param(51.68), Param(0.0), {}, Param(51.68)});
  const CompiledSystem a = compiled(sphere);
  const CompiledSystem b = compiled(toroid);
  same(rtt::paraxial::first_order(b, PathId{0}, 1), rtt::paraxial::first_order(a, PathId{0}, 1));
  same(rtt::paraxial::prescription(b, PathId{0}, 1), rtt::paraxial::prescription(a, PathId{0}, 1));
  same(rtt::paraxial::seidel(b, PathId{0}, 1, ChromaticPair{0, 2}),
       rtt::paraxial::seidel(a, PathId{0}, 1, ChromaticPair{0, 2}));
}

TEST_CASE("a real biconic, toroid or cylinder has no rotationally symmetric paraxial model",
          "[paraxial][anamorphic][errors]") {
  // ADR 0041, point 2 (table) and ADR 0042, point 2 ("Unverändert"): unequal diagonals of the
  // Hessian are a ParaxialError at the surface, with its pointer.
  const System base = singlet(0.0);
  for (const rtt::model::BaseShape& shape :
       {rtt::model::BaseShape{
            rtt::model::Biconic{Param(51.68), Param(60.0), Param(0.0), Param(0.0)}},
        rtt::model::BaseShape{
            rtt::model::Biconic{Param(51.68), Param(51.68), Param(0.0), Param(-0.5)}},
        rtt::model::BaseShape{rtt::model::Toroid{Param(51.68), Param(0.0), {}, Param(80.0)}},
        rtt::model::BaseShape{rtt::model::Cylinder{Param(51.68), Param(0.0), {}}}}) {
    const CompiledSystem cs = compiled(with_front(base, shape));
    const auto check = [](const ParaxialError& e, const std::string& part) {
      INFO(e.what());
      CHECK(std::string(e.what()).find(part) != std::string::npos);
      REQUIRE(e.surface() == rtt::model::SurfaceId("L1.S1"));
      CHECK(e.location() == "/root/children/1/surfaces/0");
    };
    if (std::holds_alternative<rtt::model::Biconic>(shape) &&
        std::get<rtt::model::Biconic>(shape).radius_x ==
            std::get<rtt::model::Biconic>(shape).radius_y) {
      // Equal curvatures, unequal conics: first order works, the Seidel cap does not.
      CHECK_NOTHROW((void)rtt::paraxial::first_order(cs, PathId{0}, 1));
      check(paraxial_error([&] { (void)rtt::paraxial::seidel(cs, PathId{0}, 1); }),
            "not rotationally symmetric to fourth order");
      continue;
    }
    check(paraxial_error([&] { (void)rtt::paraxial::first_order(cs, PathId{0}, 1); }),
          "not rotationally symmetric in second order");
    check(paraxial_error([&] { (void)rtt::paraxial::prescription(cs, PathId{0}, 1); }),
          "not rotationally symmetric in second order");
    check(paraxial_error([&] { (void)rtt::paraxial::seidel(cs, PathId{0}, 1); }),
          "not rotationally symmetric in second order");
  }
}
