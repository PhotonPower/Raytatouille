// object.medium (#282, ADR 0037): only for a finite object; the NA check for air follows the
// medium of the object space.

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <string>
#include <vector>

#include "test_support.hpp"

using namespace rtt::model;
using rtt::model::test::has_error_at;
using rtt::model::test::make_singlet;

namespace {

System finite_system() {
  System s = make_singlet();
  s.object.at_infinity = false;
  s.object.distance = Param(100.0);
  s.fields = {FieldType::ObjectHeight, {{0.0, 0.0, 1.0}}};
  return s;
}

bool has(const std::vector<Diagnostic>& d, const std::string& code, const std::string& location) {
  return std::any_of(d.begin(), d.end(),
                     [&](const Diagnostic& x) { return x.code == code && x.location == location; });
}

bool has_code(const std::vector<Diagnostic>& d, const std::string& code) {
  return std::any_of(d.begin(), d.end(), [&](const Diagnostic& x) { return x.code == code; });
}

}  // namespace

TEST_CASE("object.medium needs a finite object", "[validate][object-medium]") {
  System s = make_singlet();  // object at infinity
  s.object.medium = "CONST:1.5";
  const auto d = validate(s);
  REQUIRE(has(d, "object.medium_infinite", "/object/medium"));
  REQUIRE(has_error_at(d, "/object/medium"));
}

TEST_CASE("object.medium with a finite object is valid", "[validate][object-medium]") {
  System s = finite_system();
  s.object.medium = "CONST:1.5";
  const auto d = validate(s);
  INFO((d.empty() ? std::string() : to_string(d.front())));
  REQUIRE(d.empty());
}

TEST_CASE("na_not_physical follows the medium of the object space", "[validate][object-medium]") {
  SECTION("air as the environment, no object medium: NA >= 1 is reported, as before") {
    System s = finite_system();
    s.aperture = {SystemApertureType::ObjectSpaceNA, Param(1.2)};
    REQUIRE(has(validate(s), "aperture.na_not_physical", "/aperture/value"));
  }
  SECTION("glass as the object medium: NA >= 1 is not reported by validate") {
    System s = finite_system();
    s.aperture = {SystemApertureType::ObjectSpaceNA, Param(1.2)};
    s.object.medium = "CONST:1.5";
    REQUIRE_FALSE(has_code(validate(s), "aperture.na_not_physical"));
  }
  SECTION("air as the object medium in a glass environment: NA >= 1 is reported") {
    System s = finite_system();
    s.aperture = {SystemApertureType::ObjectSpaceNA, Param(1.2)};
    s.environment.medium = "CONST:1.5";
    s.object.medium = "AIR";
    REQUIRE(has(validate(s), "aperture.na_not_physical", "/aperture/value"));
  }
  SECTION("glass as the environment, no object medium: not reported, as before") {
    System s = finite_system();
    s.aperture = {SystemApertureType::ObjectSpaceNA, Param(1.2)};
    s.environment.medium = "CONST:1.5";
    REQUIRE_FALSE(has_code(validate(s), "aperture.na_not_physical"));
  }
}
