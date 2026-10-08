#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
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

TEST_CASE("duplicate keys are errors at the second occurrence (#68)", "[io][errors]") {
  // ADR 0008 addendum: nlohmann-json alone would keep the last value silently.
  REQUIRE(error_pointer(R"({"schema_version": "9.9.9", )" + minimal().substr(1)) ==
          "/schema_version");
  REQUIRE(error_pointer(minimal(R"(, "units": {"length": "mm", "wavelength": "um"})")) == "/units");
  REQUIRE(error_pointer(minimal("", R"("mm")", R"("0.1.0")", R"({"id": "IMG", "id": "X"})")) ==
          "/root/children/0/surfaces/0/id");
  REQUIRE(error_pointer(minimal("", R"("mm")", R"("0.1.0")",
                                R"({"id": "IMG", "aperture": {"type": "circular",)"
                                R"( "radius": 1.0, "radius": 2.0}})")) ==
          "/root/children/0/surfaces/0/aperture/radius");
  // RFC 6901 escaping in the pointer: '/' -> ~1, '~' -> ~0.
  REQUIRE(error_pointer(minimal(R"(, "a/b~c": 1, "a/b~c": 2)")) == "/a~1b~0c");
  try {
    (void)rtt::io::parse_system(minimal(R"(, "name": "a", "name": "b")"));
    FAIL("no ParseError");
  } catch (const rtt::io::ParseError& e) {
    REQUIRE(e.pointer() == "/name");
    REQUIRE(std::string(e.what()) == "/name: duplicate key 'name'");
  }
}

TEST_CASE("a number that overflows double is a ParseError (#68)", "[io][errors]") {
  // nlohmann-json reports it as out_of_range 406, which escaped parse_system before #68.
  const std::string text =
      minimal("", R"("mm")", R"("0.1.0")", R"({"id": "IMG", "pose": {"position": [0, 0, 1e400]}})");
  REQUIRE(error_pointer(text).empty());
  REQUIRE_THROWS_WITH(rtt::io::parse_system(text), Catch::Matchers::ContainsSubstring("1e400"));
}

TEST_CASE("an integer outside the range of int is a ParseError, not cut off (#35)",
          "[io][errors]") {
  // Event::order is an int; before #35 E, get<int>() silently cut 4294967296 to 0.
  const auto with_order = [](const std::string& order) {
    std::string text = minimal();
    const std::string automatic = R"("events": "auto")";
    text.replace(text.find(automatic), automatic.size(),
                 R"("events": [{"surface": "IMG", "kind": "diffract", "order": )" + order + "}]");
    return text;
  };
  for (const char* order : {"2147483648", "4294967296", "-2147483649", "18446744073709551615",
                            "-9223372036854775808"}) {
    INFO(order);
    REQUIRE(error_pointer(with_order(order)) == "/paths/0/events/0/order");
    REQUIRE_THROWS_WITH(rtt::io::parse_system(with_order(order)),
                        Catch::Matchers::ContainsSubstring("out of range"));
  }
  REQUIRE(rtt::io::parse_system(with_order("2147483647")).paths[0].events[0].order == 2147483647);
  REQUIRE(rtt::io::parse_system(with_order("-2147483648")).paths[0].events[0].order ==
          -2147483647 - 1);
  REQUIRE(error_pointer(with_order("1.0")) == "/paths/0/events/0/order");  // not an integer
}

TEST_CASE("parameters accept plain numbers and objects", "[io]") {
  const auto s = rtt::io::parse_system(
      minimal("", R"("mm")", R"("0.1.0")",
              R"({"id": "IMG", "pose": {"position": [0, {"value": 1.5, "variable": true}, 0]}})"));
  const auto& e = std::get<rtt::model::Element>(s.root.children[0].value);
  REQUIRE(e.surfaces[0].pose.position[1].value == 1.5);
  REQUIRE(e.surfaces[0].pose.position[1].variable);
}
