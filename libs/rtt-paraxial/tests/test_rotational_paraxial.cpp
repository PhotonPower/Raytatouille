// The odd asphere and the axicon in the paraxial model (#314, ADR 0043, point 3 and acceptance
// 6): an odd asphere with only even terms is the even asphere bit for bit; A2 shifts the vertex
// curvature to c + 2 A2; A3 has a first-order model but no Seidel cap; A1 and the axicon have no
// paraxial model at the vertex (F6).

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

/// m1/singlet_const with `base` at L1.S1 (/root/children/1/surfaces/0).
System singlet(const rtt::model::BaseShape& base) {
  System s = rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m1/singlet_const.rtt.json");
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
  CHECK(a.lagrange_invariant == b.lagrange_invariant);
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
    CHECK(a.surfaces[i].a == b.surfaces[i].a);
    CHECK(a.surfaces[i].a_bar == b.surfaces[i].a_bar);
    same(a.surfaces[i].terms, b.surfaces[i].terms);
  }
  same(a.sum, b.sum);
  CHECK(a.lagrange == b.lagrange);
}

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

void at_front_surface(const ParaxialError& e, const std::string& part) {
  INFO(e.what());
  CHECK(std::string(e.what()).find(part) != std::string::npos);
  REQUIRE(e.surface() == rtt::model::SurfaceId("L1.S1"));
  CHECK(e.location() == "/root/children/1/surfaces/0");
}

}  // namespace

TEST_CASE("an odd asphere with only even terms is the even asphere in the paraxial model",
          "[paraxial][rotational]") {
  // ADR 0043, acceptance 6: R = 64 is dyadic, so curvature = 1.0 / 64 is exactly the double
  // that compile forms for the even asphere; vertex expansions bit for bit give first_order,
  // prescription and seidel bit for bit.
  for (const double k : {0.0, -0.8}) {
    INFO("k = " << k);
    const System even = singlet(
        rtt::model::EvenAsphere{Param(64.0), Param(k), {Param(-2e-6), Param(3e-9), Param(-1e-12)}});
    const System odd =
        singlet(rtt::model::OddAsphere{Param(1.0 / 64.0),
                                       Param(k),
                                       {Param(0.0), Param(0.0), Param(0.0), Param(-2e-6),
                                        Param(0.0), Param(3e-9), Param(0.0), Param(-1e-12)}});
    const CompiledSystem a = compiled(even);
    const CompiledSystem b = compiled(odd);
    // Content guard: a focal system with a pupil, the cap A4 enters the spherical aberration.
    const FirstOrder fa = rtt::paraxial::first_order(a, PathId{0}, 1);
    REQUIRE(fa.efl.has_value());
    REQUIRE(fa.entrance_pupil.has_value());
    const Seidel sa = rtt::paraxial::seidel(a, PathId{0}, 1, ChromaticPair{0, 2});
    const Seidel sphere =
        rtt::paraxial::seidel(compiled(singlet(rtt::model::Conic{Param(64.0), Param(k)})),
                              PathId{0}, 1, ChromaticPair{0, 2});
    REQUIRE(sa.sum.s1 != sphere.sum.s1);
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

TEST_CASE("an odd asphere with A2 has the first order of a conic with c + 2 A2",
          "[paraxial][rotational]") {
  // ADR 0043, acceptance 6: c = 1/64 and A2 = 1/256 give c_v = 1/32 exactly, the curvature of
  // the conic R = 32. The first order depends only on c_v, so it is the conic's bit for bit.
  const System conic = singlet(rtt::model::Conic{Param(32.0), Param(-0.5)});
  const System odd = singlet(
      rtt::model::OddAsphere{Param(1.0 / 64.0), Param(-0.5), {Param(0.0), Param(1.0 / 256.0)}});
  const CompiledSystem a = compiled(conic);
  const CompiledSystem b = compiled(odd);
  for (std::uint16_t w = 0; w < 3; ++w) {
    INFO("wavelength " << w);
    same(rtt::paraxial::first_order(b, PathId{0}, w), rtt::paraxial::first_order(a, PathId{0}, w));
  }
  // Content guard: the base curvature alone (R = 64) gives another focal length.
  const FirstOrder base = rtt::paraxial::first_order(
      compiled(singlet(rtt::model::Conic{Param(64.0), Param(-0.5)})), PathId{0}, 1);
  REQUIRE(base.efl.has_value());
  CHECK(*base.efl != *rtt::paraxial::first_order(a, PathId{0}, 1).efl);
  // The Seidel cap differs (native K with A4_eff, ADR 0043, point 3), but it exists.
  CHECK_NOTHROW((void)rtt::paraxial::seidel(b, PathId{0}, 1));
}

TEST_CASE("an odd asphere with A3 has a first-order model but no Seidel cap",
          "[paraxial][rotational][errors]") {
  // ADR 0043, point 3: order 2.
  const System base = singlet(rtt::model::Conic{Param(51.68), Param(0.0)});
  const System odd = singlet(rtt::model::OddAsphere{
      Param(1.0 / 51.68), Param(0.0), {Param(0.0), Param(0.0), Param(2e-5), Param(1e-6)}});
  const CompiledSystem cs = compiled(odd);
  const FirstOrder f = rtt::paraxial::first_order(cs, PathId{0}, 1);
  same(f, rtt::paraxial::first_order(compiled(base), PathId{0}, 1));
  at_front_surface(paraxial_error([&] { (void)rtt::paraxial::seidel(cs, PathId{0}, 1); }),
                   "no expansion to fourth order");
}

TEST_CASE("an odd asphere with A1 and an axicon have no paraxial model at the vertex",
          "[paraxial][rotational][errors]") {
  // ADR 0043, point 3 (F6): order 0, a ParaxialError at the surface with its pointer.
  for (const rtt::model::BaseShape& shape :
       {rtt::model::BaseShape{rtt::model::OddAsphere{Param(0.02), Param(0.0), {Param(1e-3)}}},
        rtt::model::BaseShape{rtt::model::Axicon{Param(-5.0)}}}) {
    const CompiledSystem cs = compiled(singlet(shape));
    at_front_surface(paraxial_error([&] { (void)rtt::paraxial::first_order(cs, PathId{0}, 1); }),
                     "no paraxial model at the vertex");
    at_front_surface(paraxial_error([&] { (void)rtt::paraxial::prescription(cs, PathId{0}, 1); }),
                     "no paraxial model at the vertex");
    at_front_surface(paraxial_error([&] { (void)rtt::paraxial::seidel(cs, PathId{0}, 1); }),
                     "no paraxial model at the vertex");
  }
  // An axicon with alpha = 0 is the plane.
  const System plane = singlet(rtt::model::Plane{});
  const System flat = singlet(rtt::model::Axicon{Param(0.0)});
  same(rtt::paraxial::first_order(compiled(flat), PathId{0}, 1),
       rtt::paraxial::first_order(compiled(plane), PathId{0}, 1));
}
