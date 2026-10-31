// object.reference (#288, ADR 0038, point 7): a surface of the system, only for a finite object.
// Any surface is allowed, also the image surface or one on no path.

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <string>
#include <vector>

#include "test_support.hpp"

using namespace rtt::model;
using rtt::model::test::make_singlet;

namespace {

System finite_system(const std::string& reference) {
  System s = make_singlet();
  s.object.at_infinity = false;
  s.object.distance = Param(100.0);
  s.object.reference = SurfaceId(reference);
  s.fields = {FieldType::ObjectHeight, {{0.0, 0.0, 1.0}}};
  return s;
}

bool has(const std::vector<Diagnostic>& d, const std::string& code, const std::string& location) {
  return std::any_of(d.begin(), d.end(),
                     [&](const Diagnostic& x) { return x.code == code && x.location == location; });
}

bool has_prefix(const std::vector<Diagnostic>& d, const std::string& prefix) {
  return std::any_of(d.begin(), d.end(),
                     [&](const Diagnostic& x) { return x.code.starts_with(prefix); });
}

}  // namespace

TEST_CASE("object.reference names a surface of the system", "[validate][object-reference]") {
  CHECK(!has_prefix(validate(finite_system("L1.S1")), "object.reference"));
  // Any surface, also the stop and the image surface.
  CHECK(!has_prefix(validate(finite_system("STO")), "object.reference"));
  CHECK(!has_prefix(validate(finite_system("IMG")), "object.reference"));
  const auto d = validate(finite_system("L9.S1"));
  CHECK(has(d, "object.reference_surface_unknown", "/object/reference"));
  CHECK(!has(d, "object.reference_at_infinity", "/object/reference"));
}

TEST_CASE("object.reference needs a finite object", "[validate][object-reference]") {
  System s = make_singlet();  // object at infinity
  s.object.reference = SurfaceId("L1.S1");
  const auto d = validate(s);
  CHECK(has(d, "object.reference_at_infinity", "/object/reference"));
  CHECK(!has(d, "object.reference_surface_unknown", "/object/reference"));
}

TEST_CASE("without object.reference nothing is reported", "[validate][object-reference]") {
  System s = make_singlet();
  s.object.at_infinity = false;
  s.object.distance = Param(100.0);
  s.fields = {FieldType::ObjectHeight, {{0.0, 0.0, 1.0}}};
  CHECK(!has_prefix(validate(s), "object.reference"));
}
