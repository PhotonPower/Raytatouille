// aperture.stop (#53, ADR 0039, points 2 and 8): the stop on a surface of a lens, mirror, plate or
// thin element; at most one stop in the system; the stop surface needs an aperture.

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <string>
#include <vector>

#include "test_support.hpp"

using namespace rtt::model;
using rtt::model::test::make_singlet;

namespace {

/// The singlet without its stop element (L1 becomes child 0), the stop on `surface`.
System surface_stop(const std::string& surface) {
  System s = make_singlet();
  s.root.children.erase(s.root.children.begin());
  s.aperture.stop = SurfaceId(surface);
  return s;
}

/// A surface with a circular aperture and no other content.
Surface plane_surface(const std::string& id, double z = 0.0) {
  Surface s;
  s.id = SurfaceId(id);
  s.pose = Pose::along_z(z);
  s.aperture = CircularAperture{12.0, 0.0};
  return s;
}

bool has(const std::vector<Diagnostic>& d, const std::string& code, const std::string& location) {
  return std::any_of(d.begin(), d.end(),
                     [&](const Diagnostic& x) { return x.code == code && x.location == location; });
}

bool has_code(const std::vector<Diagnostic>& d, const std::string& code) {
  return std::any_of(d.begin(), d.end(), [&](const Diagnostic& x) { return x.code == code; });
}

bool has_stop_code(const std::vector<Diagnostic>& d) {
  return std::any_of(d.begin(), d.end(), [](const Diagnostic& x) {
    return x.code.starts_with("stop.") || x.code == "aperture.stop_missing";
  });
}

}  // namespace

TEST_CASE("aperture.stop on a lens surface is valid", "[validate][aperture-stop]") {
  CHECK(!has_stop_code(validate(surface_stop("L1.S1"))));
  CHECK(!has_stop_code(validate(surface_stop("L1.S2"))));
  // stop_size takes the stop surface (ADR 0039, point 1): no aperture.stop_missing.
  System s = surface_stop("L1.S1");
  s.aperture = {SystemApertureType::StopSize, Param(), SurfaceId("L1.S1")};
  CHECK(!has_stop_code(validate(s)));
  // Content guard: stop_size without any stop is still an error.
  s.aperture.stop.reset();
  CHECK(has(validate(s), "aperture.stop_missing", "/aperture/type"));
}

TEST_CASE("aperture.stop on a mirror, plate or thin element surface is valid",
          "[validate][aperture-stop]") {
  System s = surface_stop("L1.S1");
  Element mirror{
      "M", ElementKind::Mirror, Pose::along_z(200.0), std::nullopt, {plane_surface("M.S")}};
  Element plate{"P",
                ElementKind::Plate,
                Pose::along_z(150.0),
                "SCHOTT:N-BK7",
                {plane_surface("P.S1"), plane_surface("P.S2", 2.0)}};
  Element thin{
      "T", ElementKind::ThinElement, Pose::along_z(120.0), std::nullopt, {plane_surface("T.S")}};
  s.root.children.push_back({mirror});
  s.root.children.push_back({plate});
  s.root.children.push_back({thin});
  for (const char* id : {"M.S", "P.S1", "P.S2", "T.S"}) {
    INFO(id);
    s.aperture.stop = SurfaceId(id);
    CHECK(!has_stop_code(validate(s)));
  }
}

TEST_CASE("aperture.stop names a surface of the system", "[validate][aperture-stop]") {
  const auto d = validate(surface_stop("L9.S1"));
  CHECK(has(d, "stop.surface_unknown", "/aperture/stop"));
  CHECK(!has_code(d, "stop.surface_kind"));
}

TEST_CASE("aperture.stop on a stop or detector surface is an error", "[validate][aperture-stop]") {
  const auto image = validate(surface_stop("IMG"));
  CHECK(has(image, "stop.surface_kind", "/aperture/stop"));
  CHECK(!has_code(image, "stop.surface_unknown"));
  // The surface of a stop element: the kind is wrong, and there are two stops.
  System s = make_singlet();
  s.aperture.stop = SurfaceId("STO");
  const auto d = validate(s);
  CHECK(has(d, "stop.surface_kind", "/aperture/stop"));
  CHECK(has(d, "stop.multiple", "/aperture/stop"));
}

TEST_CASE("a stop element and aperture.stop are two stops", "[validate][aperture-stop]") {
  System s = make_singlet();  // with its stop element
  s.aperture.stop = SurfaceId("L1.S1");
  const auto d = validate(s);
  CHECK(has(d, "stop.multiple", "/aperture/stop"));
  CHECK(!has_code(d, "stop.surface_kind"));
  // Content guard: either one alone is valid.
  CHECK(!has_stop_code(validate(make_singlet())));
  CHECK(!has_stop_code(validate(surface_stop("L1.S1"))));
}

TEST_CASE("the stop surface needs an aperture", "[validate][aperture-stop]") {
  System s = surface_stop("L1.S1");
  auto& lens = std::get<Element>(s.root.children[0].value);
  lens.surfaces[0].aperture.reset();
  const auto d = validate(s);
  CHECK(has(d, "stop.aperture_missing", "/root/children/0/surfaces/0/aperture"));
  // Only the stop surface: without aperture on L1.S2 nothing is reported.
  System t = surface_stop("L1.S1");
  std::get<Element>(t.root.children[0].value).surfaces[1].aperture.reset();
  CHECK(!has_stop_code(validate(t)));
}
