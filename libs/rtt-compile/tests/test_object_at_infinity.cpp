// with_object_at_infinity (#397, ADR 0030, addendum #397): the variant of a compiled system with
// the object at infinity, one field at 0 degree and the system aperture StopSize, for efl and bfl
// of a finite object. Everything the compiled system holds besides object, fields and aperture
// is the input's, the medium of the object space included (ADR 0037).

#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <string>
#include <utility>
#include <variant>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/compile/errors.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"

using rtt::compile::CompiledSystem;
using rtt::compile::CompileError;
using rtt::compile::with_object_at_infinity;
using rtt::material::MaterialLibrary;
using rtt::model::Element;
using rtt::model::ElementKind;
using rtt::model::FieldType;
using rtt::model::Param;
using rtt::model::Pose;
using rtt::model::Surface;
using rtt::model::SurfaceId;
using rtt::model::System;

namespace {

Surface surface(const std::string& id, double z_mm, double radius) {
  Surface s;
  s.id = SurfaceId(id);
  s.pose = Pose::along_z(z_mm);
  s.shape.base = rtt::model::Conic{Param(radius), Param(0.0)};
  return s;
}

/// A lens CONST:1.6 (R 50 / -200, 5 mm) with the stop on its back surface (r = 8,
/// aperture.stop), a detector at z = 120; finite object at 200 mm in front of the reference
/// surface L.S1, object height fields, object-space NA 0.04 read real.
System finite_system() {
  System s;
  s.name = "finite";
  s.environment.medium = "VACUUM";
  s.wavelengths = {{0.5876, 1.0, true}, {0.6563, 1.0, false}};
  s.object.at_infinity = false;
  s.object.distance = Param(200.0);
  s.object.reference = SurfaceId("L.S1");
  s.aperture = {rtt::model::SystemApertureType::ObjectSpaceNA, Param(0.04)};
  s.aperture.stop = SurfaceId("L.S2");
  s.aperture.reading = rtt::model::ApertureReading::Real;
  s.fields = {FieldType::ObjectHeight, {{0.0, 0.0, 1.0}, {0.0, 3.0, 1.0}}};
  s.root.name = "root";
  s.paths = {{"main", true, {}}};
  Surface back = surface("L.S2", 5.0, -200.0);
  back.aperture = rtt::model::CircularAperture{8.0, 0.0};
  Surface img;
  img.id = SurfaceId("IMG");
  s.root.children = {
      {Element{"L",
               ElementKind::Lens,
               Pose::along_z(10.0),
               "CONST:1.6",
               {surface("L.S1", 0.0, 50.0), back}}},
      {Element{"D", ElementKind::Detector, Pose::along_z(120.0), std::nullopt, {img}}}};
  return s;
}

CompiledSystem compiled(const System& s) {
  const MaterialLibrary lib;
  return rtt::compile::compile(s, lib);
}

/// The model variant of the addendum #256 as rtt-optim built it before #397.
System model_variant(const System& s) {
  System v = s;
  v.object.at_infinity = true;
  v.object.reference.reset();
  v.fields = {FieldType::AngleDeg, {{0.0, 0.0, 1.0}}};
  v.aperture.type = rtt::model::SystemApertureType::StopSize;
  v.aperture.value = Param(0.0);
  v.aperture.reading = rtt::model::ApertureReading::Paraxial;
  return v;
}

/// Same surfaces, media and paths, bit for bit.
void check_same_optics(const CompiledSystem& a, const CompiledSystem& b) {
  REQUIRE(a.surfaces().size() == b.surfaces().size());
  for (std::size_t i = 0; i < a.surfaces().size(); ++i) {
    CHECK(a.surfaces()[i].id == b.surfaces()[i].id);
    CHECK(a.surfaces()[i].to_global.rotation() == b.surfaces()[i].to_global.rotation());
    CHECK(a.surfaces()[i].to_global.translation() == b.surfaces()[i].to_global.translation());
  }
  REQUIRE(a.media().size() == b.media().size());
  for (std::size_t m = 0; m < a.media().size(); ++m)
    CHECK(a.media()[m].index == b.media()[m].index);
  REQUIRE(a.paths().size() == b.paths().size());
  for (std::size_t p = 0; p < a.paths().size(); ++p) {
    const auto& ea = a.paths()[p].events;
    const auto& eb = b.paths()[p].events;
    REQUIRE(ea.size() == eb.size());
    for (std::size_t e = 0; e < ea.size(); ++e) {
      CHECK(ea[e].surface == eb[e].surface);
      CHECK(ea[e].medium_before == eb[e].medium_before);
      CHECK(ea[e].medium_after == eb[e].medium_after);
    }
  }
  CHECK(a.stop_surface() == b.stop_surface());
  CHECK(a.wavelengths_um() == b.wavelengths_um());
  CHECK(a.reference_wavelength() == b.reference_wavelength());
  CHECK(a.configuration() == b.configuration());
}

}  // namespace

TEST_CASE("with_object_at_infinity: object, fields and aperture of the variant (#397)",
          "[compile][variant]") {
  const CompiledSystem nominal = compiled(finite_system());
  REQUIRE_FALSE(nominal.object().at_infinity);
  REQUIRE(nominal.object_reference_z().has_value());
  const CompiledSystem v = with_object_at_infinity(nominal);
  CHECK(v.object().at_infinity);
  CHECK_FALSE(v.object().reference.has_value());
  CHECK_FALSE(v.object_reference_z().has_value());
  CHECK(v.object().medium == nominal.object().medium);
  CHECK(v.fields() == rtt::model::FieldSet{FieldType::AngleDeg, {{0.0, 0.0, 1.0}}});
  CHECK(v.aperture().type == rtt::model::SystemApertureType::StopSize);
  CHECK(v.aperture().value == Param(0.0));
  CHECK(v.aperture().reading == rtt::model::ApertureReading::Paraxial);
  CHECK(v.aperture().stop == SurfaceId("L.S2"));  // a stop on a surface stays (#400)
  check_same_optics(v, nominal);
  // The input is unchanged.
  CHECK_FALSE(nominal.object().at_infinity);
  CHECK(nominal.fields().points.size() == 2);
}

TEST_CASE("with_object_at_infinity equals the compiled model variant (#397)",
          "[compile][variant]") {
  // Without an object medium the model variant of the addendum #256 compiles; both carry the same
  // object, fields and aperture and the same optics.
  const System s = finite_system();
  const CompiledSystem a = with_object_at_infinity(compiled(s));
  const CompiledSystem b = compiled(model_variant(s));
  CHECK(a.object() == b.object());
  CHECK(a.fields() == b.fields());
  CHECK(a.aperture() == b.aperture());
  check_same_optics(a, b);
}

TEST_CASE("with_object_at_infinity keeps the medium of the object space (#397)",
          "[compile][variant]") {
  // ADR 0037: the model variant fails validate (object.medium_infinite); the variant of the
  // compiled system keeps the medium, so the first surface still refracts out of n_obj.
  System s = finite_system();
  s.object.medium = "CONST:1.33";
  CHECK_THROWS_AS(compiled(model_variant(s)), CompileError);
  const CompiledSystem nominal = compiled(s);
  const CompiledSystem v = with_object_at_infinity(nominal);
  CHECK(v.object().medium == std::optional<std::string>("CONST:1.33"));
  check_same_optics(v, nominal);
  const auto& first = v.paths().at(0).events.at(0);
  CHECK(v.media().at(first.medium_before).index.at(0).real() == 1.33);
}

TEST_CASE("with_object_at_infinity of an object at infinity changes only fields and aperture",
          "[compile][variant]") {
  System s = finite_system();
  s.object.at_infinity = true;
  s.object.reference.reset();
  s.aperture = {rtt::model::SystemApertureType::EntrancePupilDiameter, Param(10.0)};
  s.aperture.stop = SurfaceId("L.S2");
  s.fields = {FieldType::AngleDeg, {{0.0, 0.0, 1.0}, {0.0, 5.0, 1.0}}};
  const CompiledSystem nominal = compiled(s);
  const CompiledSystem v = with_object_at_infinity(nominal);
  CHECK(v.object() == nominal.object());
  CHECK(v.fields().points.size() == 1);
  CHECK(v.aperture().type == rtt::model::SystemApertureType::StopSize);
  check_same_optics(v, nominal);
}
