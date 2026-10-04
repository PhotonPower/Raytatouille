#include <catch2/catch_test_macros.hpp>
#include <optional>
#include <string>

#include "test_support.hpp"

using namespace rtt::model;
using rtt::model::test::element;
using rtt::model::test::has_error_at;
using rtt::model::test::make_singlet;

namespace {

/// The singlet of make_singlet() turned into a cemented doublet: three surfaces, the lens
/// material is reset so that each test chooses the material form itself.
System make_doublet() {
  System s = make_singlet();
  Element& lens = element(s, 1);
  Surface s3;
  s3.id = SurfaceId("L1.S3");
  s3.pose = Pose::along_z(6.0);
  s3.shape.base = Conic{Param(-120.0), Param(0.0)};
  lens.surfaces.push_back(s3);
  lens.material.reset();
  return s;
}

}  // namespace

TEST_CASE("cemented lens with one material per segment is valid", "[validate][segments]") {
  System s = make_doublet();
  element(s, 1).segment_materials = {"SCHOTT:N-BK7", "SCHOTT:F2"};
  const auto d = validate(s);
  INFO((d.empty() ? std::string() : to_string(d.front())));
  REQUIRE(d.empty());
}

TEST_CASE("cemented lens with the single-material shorthand is valid", "[validate][segments]") {
  System s = make_doublet();
  element(s, 1).material = "SCHOTT:N-BK7";
  REQUIRE(validate(s).empty());
}

TEST_CASE("a list with one entry is valid for a lens with two surfaces", "[validate][segments]") {
  System s = make_singlet();
  element(s, 1).material.reset();
  element(s, 1).segment_materials = {"SCHOTT:N-BK7"};
  REQUIRE(validate(s).empty());
}

TEST_CASE("plate with one material per segment is valid", "[validate][segments]") {
  System s = make_doublet();
  Element& plate = element(s, 1);
  plate.kind = ElementKind::Plate;
  for (Surface& surface : plate.surfaces) surface.shape.base = Plane{};
  plate.segment_materials = {"SCHOTT:N-BK7", "SCHOTT:F2"};
  REQUIRE(validate(s).empty());
}

TEST_CASE("segment material list must have N-1 entries", "[validate][segments]") {
  System s = make_doublet();
  element(s, 1).segment_materials = {"SCHOTT:N-BK7"};
  REQUIRE(has_error_at(validate(s), "/root/children/1/material"));

  s = make_doublet();
  element(s, 1).segment_materials = {"SCHOTT:N-BK7", "SCHOTT:F2", "SCHOTT:N-SF6"};
  REQUIRE(has_error_at(validate(s), "/root/children/1/material"));

  s = make_doublet();
  element(s, 1).kind = ElementKind::Plate;
  for (Surface& surface : element(s, 1).surfaces) surface.shape.base = Plane{};
  element(s, 1).segment_materials = {"SCHOTT:N-BK7"};
  REQUIRE(has_error_at(validate(s), "/root/children/1/material"));
}

TEST_CASE("empty segment materials are reported at the list entry", "[validate][segments]") {
  System s = make_doublet();
  element(s, 1).segment_materials = {"SCHOTT:N-BK7", ""};
  REQUIRE(has_error_at(validate(s), "/root/children/1/material/1"));

  s = make_doublet();
  element(s, 1).segment_materials = {"", "SCHOTT:F2"};
  REQUIRE(has_error_at(validate(s), "/root/children/1/material/0"));
}

TEST_CASE("shorthand and list must not both be set", "[validate][segments]") {
  System s = make_doublet();
  element(s, 1).material = "SCHOTT:N-BK7";
  element(s, 1).segment_materials = {"SCHOTT:N-BK7", "SCHOTT:F2"};
  REQUIRE(has_error_at(validate(s), "/root/children/1/material"));
}

TEST_CASE("cemented lens without any material is invalid", "[validate][segments]") {
  const System s = make_doublet();
  REQUIRE(has_error_at(validate(s), "/root/children/1/material"));
}

TEST_CASE("other element kinds have no segment material list", "[validate][segments]") {
  System s = make_singlet();
  element(s, 2).segment_materials = {"SCHOTT:N-BK7"};  // detector
  REQUIRE(has_error_at(validate(s), "/root/children/2/material"));

  s = make_singlet();
  element(s, 0).segment_materials = {"SCHOTT:N-BK7"};  // stop
  REQUIRE(has_error_at(validate(s), "/root/children/0/material"));

  s = make_singlet();
  Element& mirror = element(s, 1);
  mirror.kind = ElementKind::Mirror;
  mirror.material.reset();
  mirror.segment_materials = {"SCHOTT:N-BK7"};
  REQUIRE(has_error_at(validate(s), "/root/children/1/material"));
}

TEST_CASE("segment_material resolves shorthand and list", "[model][segments]") {
  System s = make_doublet();
  Element& lens = element(s, 1);
  REQUIRE_FALSE(lens.segment_material(0).has_value());

  lens.material = "SCHOTT:N-BK7";
  REQUIRE(lens.segment_material(0) == std::optional<std::string>("SCHOTT:N-BK7"));
  REQUIRE(lens.segment_material(1) == std::optional<std::string>("SCHOTT:N-BK7"));
  REQUIRE_FALSE(lens.segment_material(2).has_value());  // only 2 segments for 3 surfaces

  lens.material.reset();
  lens.segment_materials = {"SCHOTT:N-BK7", "SCHOTT:F2"};
  REQUIRE(lens.segment_material(0) == std::optional<std::string>("SCHOTT:N-BK7"));
  REQUIRE(lens.segment_material(1) == std::optional<std::string>("SCHOTT:F2"));
  REQUIRE_FALSE(lens.segment_material(2).has_value());
}
