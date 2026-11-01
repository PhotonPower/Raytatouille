// aperture.reading (#296, ADR 0040, point 1): real only for object_na; another aperture type with
// reading real can only come from an edit and is reported at /aperture/reading.

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <string>
#include <vector>

#include "test_support.hpp"

using namespace rtt::model;
using rtt::model::test::make_singlet;

namespace {

bool has(const std::vector<Diagnostic>& d, const std::string& code, const std::string& location) {
  return std::any_of(d.begin(), d.end(),
                     [&](const Diagnostic& x) { return x.code == code && x.location == location; });
}

}  // namespace

TEST_CASE("aperture.reading real needs object_na", "[validate][aperture-reading]") {
  System s = make_singlet();  // epd
  s.aperture.reading = ApertureReading::Real;
  CHECK(has(validate(s), "aperture.reading_not_na", "/aperture/reading"));
  // Content guard: object_na with a finite object takes it.
  s.object.at_infinity = false;
  s.object.distance = Param(100.0);
  s.fields = {FieldType::ObjectHeight, {{0.0, 0.0, 1.0}}};
  s.aperture.type = SystemApertureType::ObjectSpaceNA;
  s.aperture.value = Param(0.1);
  CHECK(!has(validate(s), "aperture.reading_not_na", "/aperture/reading"));
  // Paraxial is the default and valid for every type.
  System p = make_singlet();
  CHECK(!has(validate(p), "aperture.reading_not_na", "/aperture/reading"));
}
