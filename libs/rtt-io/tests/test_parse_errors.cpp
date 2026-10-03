#include <catch2/catch_test_macros.hpp>
#include <string>

#include "rtt/io/json_io.hpp"

namespace {

/// Minimal valid file; tests replace one marked fragment to provoke a single error.
std::string minimal(const std::string& extra_top = "",
                    const std::string& units = R"("mm")",
                    const std::string& version = R"("0.1.0")",
                    const std::string& surface = R"({"id": "IMG"})") {
  return R"({"schema_version": )" + version + R"(, "units": {"length": )" + units +
         R"(, "wavelength": "um"}, "wavelengths": [{"um": 0.55, "reference": true}],)"
         R"( "aperture": {"type": "epd", "value": 10.0}, "fields": {"points": [{}]},)"
         R"( "root": {"type": "assembly", "name": "sys", "children": [)"
         R"({"type": "detector", "name": "image", "surfaces": [)" +
         surface + R"(]}]}, "paths": [{"name": "main", "events": "auto"}])" + extra_top + "}";
}

std::string error_pointer(const std::string& text) {
  try {
    (void)rtt::io::parse_system(text);
  } catch (const rtt::io::ParseError& e) {
    return e.pointer();
  }
  return "<no error>";
}

}  // namespace

TEST_CASE("minimal file parses", "[io][errors]") {
  REQUIRE_NOTHROW(rtt::io::parse_system(minimal()));
}

TEST_CASE("structural errors carry a JSON pointer", "[io][errors]") {
  REQUIRE(error_pointer(minimal(R"(, "colour": "red")")) == "/colour");
  REQUIRE(error_pointer(minimal("", R"("inch")")) == "/units/length");
  REQUIRE(error_pointer(minimal("", R"("mm")", R"("1.0.0")")) == "/schema_version");
  REQUIRE(error_pointer(minimal("", R"("mm")", R"("0.1")")) == "/schema_version");
  REQUIRE(error_pointer(minimal("", R"("mm")", R"("0.1.0")", R"({"id": 5})")) ==
          "/root/children/0/surfaces/0/id");
  REQUIRE(error_pointer(minimal("", R"("mm")", R"("0.1.0")",
                                R"({"id": "IMG", "shape": {"base": {"type": "torus"}}})")) ==
          "/root/children/0/surfaces/0/shape/base/type");
  REQUIRE(error_pointer(minimal("", R"("mm")", R"("0.1.0")",
                                R"({"id": "IMG", "shape": {"base": {"type": "conic"}}})")) ==
          "/root/children/0/surfaces/0/shape/base");
  REQUIRE(error_pointer(minimal("", R"("mm")", R"("0.1.0")",
                                R"({"id": "IMG", "pose": {"position": [0, 0]}})")) ==
          "/root/children/0/surfaces/0/pose/position");
  REQUIRE(error_pointer(minimal("", R"("mm")", R"("0.1.0")",
                                R"({"id": "IMG", "pose": {"position": [0, 0, true]}})")) ==
          "/root/children/0/surfaces/0/pose/position/2");
}

TEST_CASE("invalid JSON is reported", "[io][errors]") {
  REQUIRE_THROWS_AS(rtt::io::parse_system("{ not json"), rtt::io::ParseError);
}

TEST_CASE("parameters accept plain numbers and objects", "[io]") {
  const auto s = rtt::io::parse_system(
      minimal("", R"("mm")", R"("0.1.0")",
              R"({"id": "IMG", "pose": {"position": [0, {"value": 1.5, "variable": true}, 0]}})"));
  const auto& e = std::get<rtt::model::Element>(s.root.children[0].value);
  REQUIRE(e.surfaces[0].pose.position[1].value == 1.5);
  REQUIRE(e.surfaces[0].pose.position[1].variable);
}
