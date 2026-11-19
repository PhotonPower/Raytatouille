// validate() on the regular polygon and polygon apertures (#320, ADR 0050, point 2): every new
// code at its pointer, the simple-polygon rules, and a bound circumradius skipped before the
// resolution.

#include <algorithm>
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <limits>
#include <string>
#include <vector>

#include "rtt/model/model.hpp"
#include "test_support.hpp"

using namespace rtt::model;
using namespace rtt::model::test;

namespace {

const std::string kAp = "/root/children/1/surfaces/0/aperture";

std::vector<Diagnostic> with(const Aperture& a) {
  System s = make_singlet();
  element(s, 1).surfaces[0].aperture = a;
  return validate(s);
}

/// Exactly one error, with `code` at `location`.
bool only(const std::vector<Diagnostic>& d, const std::string& code, const std::string& location) {
  std::size_t errors = 0;
  bool found = false;
  for (const Diagnostic& x : d) {
    if (x.severity != Severity::Error) continue;
    ++errors;
    found = found || (x.code == code && x.location == location);
  }
  return errors == 1 && found;
}

PolygonAperture poly(std::vector<std::array<double, 2>> v) {
  return PolygonAperture{std::move(v)};
}

}  // namespace

TEST_CASE("regular polygon: valid forms and every code", "[validate][polygon]") {
  CHECK_FALSE(has_errors(with(RegularPolygonAperture{3, Param(5.0), 0.0})));
  CHECK_FALSE(has_errors(with(RegularPolygonAperture{1024, Param(5.0), -720.0})));
  CHECK(only(with(RegularPolygonAperture{2, Param(5.0), 0.0}), "surface_aperture.sides_invalid",
             kAp + "/sides"));
  CHECK(only(with(RegularPolygonAperture{1025, Param(5.0), 0.0}), "surface_aperture.sides_invalid",
             kAp + "/sides"));
  CHECK(only(with(RegularPolygonAperture{6, Param(-1.0), 0.0}),
             "surface_aperture.circumradius_invalid", kAp + "/circumradius"));
  CHECK(only(with(RegularPolygonAperture{6, Param(std::numeric_limits<double>::quiet_NaN()), 0.0}),
             "surface_aperture.circumradius_invalid", kAp + "/circumradius"));
  CHECK(only(with(RegularPolygonAperture{6, Param(5.0), std::numeric_limits<double>::infinity()}),
             "surface_aperture.rotation_invalid", kAp + "/rotation_deg"));
}

TEST_CASE("regular polygon: a bound circumradius is checked after the resolution",
          "[validate][polygon]") {
  System s = make_singlet();
  ParameterRow r;
  r.name = "R";
  r.form = -1.0;  // invalid, but only compile sees the resolved value (#353)
  s.parameters.push_back(r);
  element(s, 1).surfaces[0].aperture = RegularPolygonAperture{6, Param::bound("R"), 0.0};
  CHECK_FALSE(has_errors(validate(s)));
  CHECK(has_error_at(validate(resolve_parameters(s, 0)), kAp + "/circumradius"));
}

TEST_CASE("polygon: simple polygons are valid, convex or concave, either orientation",
          "[validate][polygon]") {
  std::vector<std::array<double, 2>> l = {{0.0, 0.0}, {4.0, 0.0}, {4.0, 2.0},
                                          {2.0, 2.0}, {2.0, 4.0}, {0.0, 4.0}};
  CHECK_FALSE(has_errors(with(poly(l))));
  std::reverse(l.begin(), l.end());
  CHECK_FALSE(has_errors(with(poly(l))));
  CHECK_FALSE(has_errors(with(poly({{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}}))));
  // A straight vertex (collinear, same direction) is allowed.
  CHECK_FALSE(has_errors(with(poly({{0.0, 0.0}, {1.0, 0.0}, {2.0, 0.0}, {2.0, 2.0}}))));
}

TEST_CASE("polygon: vertices_invalid and polygon_self_intersecting", "[validate][polygon]") {
  const std::string v = kAp + "/vertices";
  CHECK(only(with(poly({{0.0, 0.0}, {1.0, 0.0}})), "surface_aperture.vertices_invalid", v));
  CHECK(only(with(poly({{0.0, 0.0}, {1.0, std::numeric_limits<double>::infinity()}, {0.0, 1.0}})),
             "surface_aperture.vertices_invalid", v + "/1"));
  // two equal consecutive vertices, also the last and the first
  CHECK(only(with(poly({{0.0, 0.0}, {1.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}})),
             "surface_aperture.vertices_invalid", v + "/2"));
  CHECK(only(with(poly({{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}, {0.0, 0.0}})),
             "surface_aperture.vertices_invalid", v + "/0"));
  // a U-turn at (2, 0): the edge back to (1, 0) overlaps the edge to it
  CHECK(only(with(poly({{0.0, 0.0}, {2.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}})),
             "surface_aperture.vertices_invalid", v + "/1"));
  // three collinear points: area 0, found as the U-turn at the last of them
  CHECK(only(with(poly({{0.0, 0.0}, {1.0, 1.0}, {2.0, 2.0}})), "surface_aperture.vertices_invalid",
             v + "/2"));
  // a bow tie, and a polygon whose non-adjacent edges only touch at a vertex
  CHECK(only(with(poly({{0.0, 0.0}, {2.0, 2.0}, {2.0, 0.0}, {0.0, 2.0}})),
             "surface_aperture.polygon_self_intersecting", v));
  CHECK(only(with(poly({{0.0, 0.0}, {4.0, 0.0}, {2.0, 2.0}, {4.0, 4.0}, {0.0, 4.0}, {2.0, 2.0}})),
             "surface_aperture.polygon_self_intersecting", v));
  CHECK(only(with(poly({{0.0, 0.0}, {4.0, 0.0}, {4.0, 4.0}, {2.0, 0.0}, {0.0, 4.0}})),
             "surface_aperture.polygon_self_intersecting", v));
}
