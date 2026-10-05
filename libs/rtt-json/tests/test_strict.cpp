#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <nlohmann/json.hpp>
#include <string>
#include <vector>

#include "rtt/json/strict.hpp"

using Catch::Matchers::ContainsSubstring;
using rtt::json::parse_strict;
using rtt::json::StrictParseError;

namespace {

StrictParseError error_of(const std::string& text) {
  try {
    (void)parse_strict(text);
  } catch (const StrictParseError& e) {
    return e;
  }
  FAIL("no StrictParseError for " << text);
  return StrictParseError("", "");
}

}  // namespace

TEST_CASE("valid documents parse like nlohmann::json::parse", "[json]") {
  for (const std::string text :
       {R"({"a": 1, "b": [1, 2, {"c": "d"}], "e": {"f": null}})", R"([1, "x", true])", "3.5",
        R"({"a": {"a": {"a": 1}}})", R"([{"k": 1}, {"k": 2}])", "{}", "[]"}) {
    INFO(text);
    REQUIRE(parse_strict(text) == nlohmann::json::parse(text));
  }
}

TEST_CASE("a duplicate key is an error with the pointer of its second occurrence", "[json]") {
  struct Case {
    std::string text;
    std::string pointer;
  };
  const std::vector<Case> cases = {
      {R"({"a": 1, "a": 2})", "/a"},
      {R"({"x": {"y": 1, "z": 2, "y": 3}})", "/x/y"},
      {R"({"list": [{"k": 1}, {"k": 2, "k": 3}]})", "/list/1/k"},
      // Array indices count primitive, array and object elements before the duplicate.
      {R"({"a": [1, [2, 3], {"b": 1, "b": 2}]})", "/a/2/b"},
      {R"({"c": {"d": [{"e": 1}, [], {"e": 2, "e": 3}]}})", "/c/d/2/e"},
      {R"([{"a": 1, "a": 1}])", "/0/a"},
      {R"({"o": {"p": [{"q": {"r": 1, "r": 1}}]}})", "/o/p/0/q/r"},
      // RFC 6901 escaping of the reference tokens: '/' -> ~1, '~' -> ~0.
      {R"({"a/b": {"~": 1, "~": 2}})", "/a~1b/~0"},
      // The duplicate of a key whose first value was an object or array.
      {R"({"s": {"t": 1}, "s": 2})", "/s"},
      {R"({"s": [1], "s": [2]})", "/s"},
  };
  for (const Case& c : cases) {
    INFO(c.text);
    const StrictParseError e = error_of(c.text);
    REQUIRE(e.pointer() == c.pointer);
    REQUIRE_THAT(e.message(), ContainsSubstring("duplicate key"));
    REQUIRE_THAT(std::string(e.what()), ContainsSubstring(c.pointer + ": duplicate key"));
  }
}

TEST_CASE("syntax errors and number overflow have an empty pointer", "[json]") {
  for (const std::string text : {"{", R"({"a": })", R"({"a": 1e400})", "[1e999]", ""}) {
    INFO(text);
    const StrictParseError e = error_of(text);
    REQUIRE(e.pointer().empty());
    REQUIRE_THAT(e.message(), ContainsSubstring("invalid JSON"));
    REQUIRE(std::string(e.what()).starts_with("/: invalid JSON: "));
  }
}

TEST_CASE("the message names the duplicate key", "[json]") {
  const StrictParseError e = error_of(R"({"x": {"k~/": 1, "k~/": 2}})");
  REQUIRE(e.pointer() == "/x/k~0~1");
  REQUIRE(e.message() == "duplicate key 'k~/'");
  REQUIRE(std::string(e.what()) == "/x/k~0~1: duplicate key 'k~/'");
}

TEST_CASE("pointer_token escapes per RFC 6901", "[json]") {
  REQUIRE(rtt::json::pointer_token("plain") == "plain");
  REQUIRE(rtt::json::pointer_token("a/b") == "a~1b");
  REQUIRE(rtt::json::pointer_token("~1") == "~01");
  REQUIRE(rtt::json::pointer_token("") == "");
}
