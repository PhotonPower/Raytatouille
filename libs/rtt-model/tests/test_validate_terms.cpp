// Validation of shape terms (#316, ADR 0051, point 6, PR 2): a plate stays plane-parallel, and
// the normalization extent warning covers the Zernike terms.

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <string>
#include <vector>

#include "rtt/model/model.hpp"
#include "test_support.hpp"

using namespace rtt::model;
using namespace rtt::model::test;

namespace {

bool has(const std::vector<Diagnostic>& d, const std::string& code, const std::string& location) {
  return std::any_of(d.begin(), d.end(),
                     [&](const Diagnostic& x) { return x.code == code && x.location == location; });
}

bool has_code(const std::vector<Diagnostic>& d, const std::string& code) {
  return std::any_of(d.begin(), d.end(), [&](const Diagnostic& x) { return x.code == code; });
}

ZernikeSag defocus(Param radius) {
  return ZernikeSag{std::move(radius), {Param(0.0), Param(0.0), Param(0.0), Param(1e-4)}};
}

}  // namespace

TEST_CASE("a plate surface with terms is not a plane (ADR 0051, point 6)", "[validate][terms]") {
  System s = make_singlet();
  Element& e = element(s, 1);
  e.kind = ElementKind::Plate;
  for (Surface& f : e.surfaces) f.shape = ShapeStack{};
  REQUIRE_FALSE(has_code(validate(s), "element.plate_surface_not_plane"));  // planes pass
  e.surfaces[0].shape.terms = {defocus(Param(10.0))};
  const std::vector<Diagnostic> d = validate(s);
  CHECK(has(d, "element.plate_surface_not_plane", "/root/children/1/surfaces/0/shape"));
  CHECK_FALSE(has(d, "element.plate_surface_not_plane", "/root/children/1/surfaces/1/shape"));
  // The same term on a lens surface is fine.
  e.kind = ElementKind::Lens;
  CHECK_FALSE(has_code(validate(s), "element.plate_surface_not_plane"));
}

TEST_CASE("the normalization extent warning covers Zernike terms (ADR 0051)", "[validate][terms]") {
  const std::string term = "/root/children/1/surfaces/0/shape/terms/0/normalization_radius";
  const auto with = [](Param radius, double aperture) {
    System s = make_singlet();
    Surface& f = element(s, 1).surfaces[0];
    f.aperture = CircularAperture{Param(aperture), Param(0.0)};
    f.shape.terms = {defocus(std::move(radius)), ZernikeFringeSag{Param(20.0), {Param(1e-4)}}};
    return validate(s);
  };
  // Aperture 12 beyond R = 10: warning at the term; the Fringe term with R = 20 has none.
  const std::vector<Diagnostic> beyond = with(Param(10.0), 12.0);
  CHECK(has(beyond, "shape.normalization_radius_exceeded", term));
  CHECK_FALSE(has(beyond, "shape.normalization_radius_exceeded",
                  "/root/children/1/surfaces/0/shape/terms/1/normalization_radius"));
  // Equal: no warning (1e-9 relative).
  CHECK_FALSE(has_code(with(Param(12.0), 12.0), "shape.normalization_radius_exceeded"));
  // A bound radius has its value only after resolving: no warning here.
  CHECK_FALSE(has(with(Param::bound("R"), 12.0), "shape.normalization_radius_exceeded", term));
}
