// JSON Patch engine on generic JSON (ADR 0024 point 3). The Appendix A cases are the examples of
// RFC 6902 (IETF, April 2013), Appendix A.1-A.16, with their documents and patches verbatim;
// Code Components of the RFC are under the Simplified BSD License (RFC 6902, "Copyright
// Notice"; docs/quellen.md).

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cstddef>
#include <nlohmann/json.hpp>
#include <optional>
#include <string>

#include "json_patch.hpp"
#include "rtt/json/strict.hpp"

using nlohmann::json;
using rtt::io::detail::apply_json_patch;
using rtt::io::detail::PatchApplied;
using rtt::io::detail::PatchFailure;

namespace {

/// Applies the patch, checks that its inverse restores the document and returns the result.
json applied(const std::string& document, const std::string& patch) {
  const json doc = json::parse(document);
  const PatchApplied r = apply_json_patch(doc, json::parse(patch));
  INFO("inverse " << r.inverse.dump());
  REQUIRE(apply_json_patch(r.document, r.inverse).document == doc);
  return r.document;
}

PatchFailure failure(const std::string& document, const std::string& patch) {
  const json doc = json::parse(document);
  try {
    (void)apply_json_patch(doc, json::parse(patch));
  } catch (const PatchFailure& e) {
    return e;
  }
  FAIL("no PatchFailure for " << patch);
  return PatchFailure("", "", std::nullopt, "");
}

}  // namespace

TEST_CASE("RFC 6902 Appendix A.1-A.7, A.10, A.11, A.14, A.16: results", "[io][patch][rfc6902]") {
  // A.1 Adding an Object Member
  CHECK(applied(R"({ "foo": "bar"})", R"([ { "op": "add", "path": "/baz", "value": "qux" } ])") ==
        json::parse(R"({ "baz": "qux", "foo": "bar" })"));
  // A.2 Adding an Array Element
  CHECK(applied(R"({ "foo": [ "bar", "baz" ] })",
                R"([ { "op": "add", "path": "/foo/1", "value": "qux" } ])") ==
        json::parse(R"({ "foo": [ "bar", "qux", "baz" ] })"));
  // A.3 Removing an Object Member
  CHECK(applied(R"({ "baz": "qux", "foo": "bar" })", R"([ { "op": "remove", "path": "/baz" } ])") ==
        json::parse(R"({ "foo": "bar" })"));
  // A.4 Removing an Array Element
  CHECK(applied(R"({ "foo": [ "bar", "qux", "baz" ] })",
                R"([ { "op": "remove", "path": "/foo/1" } ])") ==
        json::parse(R"({ "foo": [ "bar", "baz" ] })"));
  // A.5 Replacing a Value
  CHECK(applied(R"({ "baz": "qux", "foo": "bar" })",
                R"([ { "op": "replace", "path": "/baz", "value": "boo" } ])") ==
        json::parse(R"({ "baz": "boo", "foo": "bar" })"));
  // A.6 Moving a Value
  CHECK(
      applied(R"({ "foo": { "bar": "baz", "waldo": "fred" }, "qux": { "corge": "grault" } })",
              R"([ { "op": "move", "from": "/foo/waldo", "path": "/qux/thud" } ])") ==
      json::parse(R"({ "foo": { "bar": "baz" }, "qux": { "corge": "grault", "thud": "fred" } })"));
  // A.7 Moving an Array Element
  CHECK(applied(R"({ "foo": [ "all", "grass", "cows", "eat" ] })",
                R"([ { "op": "move", "from": "/foo/1", "path": "/foo/3" } ])") ==
        json::parse(R"({ "foo": [ "all", "cows", "eat", "grass" ] })"));
  // A.10 Adding a Nested Member Object
  CHECK(applied(R"({ "foo": "bar" })",
                R"([ { "op": "add", "path": "/child", "value": { "grandchild": { } } } ])") ==
        json::parse(R"({ "foo": "bar", "child": { "grandchild": { } } })"));
  // A.11 Ignoring Unrecognized Elements
  CHECK(applied(R"({ "foo": "bar" })",
                R"([ { "op": "add", "path": "/baz", "value": "qux", "xyz": 123 } ])") ==
        json::parse(R"({ "foo": "bar", "baz": "qux" })"));
  // A.14 ~ Escape Ordering
  CHECK(applied(R"({ "/": 9, "~1": 10 })", R"([ {"op": "test", "path": "/~01", "value": 10} ])") ==
        json::parse(R"({ "/": 9, "~1": 10 })"));
  // A.16 Adding an Array Value
  CHECK(applied(R"({ "foo": ["bar"] })",
                R"([ { "op": "add", "path": "/foo/-", "value": ["abc", "def"] } ])") ==
        json::parse(R"({ "foo": ["bar", ["abc", "def"]] })"));
}

TEST_CASE("RFC 6902 Appendix A.8, A.9, A.12, A.13, A.15: tests and errors",
          "[io][patch][rfc6902]") {
  // A.8 Testing a Value: Success
  const std::string a8 = R"({ "baz": "qux", "foo": [ "a", 2, "c" ] })";
  CHECK(applied(a8, R"([ { "op": "test", "path": "/baz", "value": "qux" },
                         { "op": "test", "path": "/foo/1", "value": 2 } ])") == json::parse(a8));
  // A.9 Testing a Value: Error
  const PatchFailure a9 =
      failure(R"({ "baz": "qux" })", R"([ { "op": "test", "path": "/baz", "value": "bar" } ])");
  CHECK(a9.code() == "edit.test_failed");
  CHECK((a9.pointer() == "/baz" && a9.op_index() == 0));
  // A.12 Adding to a Nonexistent Target
  const PatchFailure a12 =
      failure(R"({ "foo": "bar" })", R"([ { "op": "add", "path": "/baz/bat", "value": "qux" } ])");
  CHECK(a12.code() == "edit.path_not_found");
  CHECK(a12.pointer() == "/baz/bat");
  // A.13 Invalid JSON Patch Document: a duplicate "op"; the strict parser rejects it (ADR 0020),
  // so it never reaches the engine as an "add" or a "remove".
  CHECK_THROWS_AS(rtt::json::parse_strict(
                      R"([ { "op": "add", "path": "/baz", "value": "qux", "op": "remove" } ])"),
                  rtt::json::StrictParseError);
  // A.15 Comparing Strings and Numbers
  CHECK(failure(R"({ "/": 9, "~1": 10 })", R"([ {"op": "test", "path": "/~01", "value": "10"} ])")
            .code() == "edit.test_failed");
}

TEST_CASE("deviations of nlohmann::json::patch are errors here", "[io][patch]") {
  // add, copy and move into a primitive parent (nlohmann: JSON_ASSERT abort in debug builds)
  CHECK(failure(R"({"a": "s"})", R"([{"op": "add", "path": "/a/b", "value": 1}])").code() ==
        "edit.path_not_found");
  CHECK(failure(R"({"a": "s", "b": 1})", R"([{"op": "copy", "from": "/b", "path": "/a/b"}])")
            .code() == "edit.path_not_found");
  CHECK(failure(R"({"a": "s", "b": 1})", R"([{"op": "move", "from": "/b", "path": "/a/b"}])")
            .code() == "edit.path_not_found");
  // remove under a primitive parent (nlohmann: silently nothing)
  CHECK(failure(R"({"name": "x"})", R"([{"op": "remove", "path": "/name/a"}])").code() ==
        "edit.path_not_found");
  // move into its own child (nlohmann: not rejected)
  const PatchFailure into_child =
      failure(R"({"c": [{"c": []}]})", R"([{"op": "move", "from": "/c/0", "path": "/c/0/c/0"}])");
  CHECK(into_child.code() == "edit.patch_invalid");
  CHECK(into_child.op_index() == 0);
}

TEST_CASE("pointers, indices and the move target", "[io][patch]") {
  // "/a/1" is no prefix of "/a/10": token by token
  CHECK(applied(R"({"a": [0,1,2,3,4,5,6,7,8,9,10,11]})",
                R"([{"op": "move", "from": "/a/1", "path": "/a/10"}])") ==
        json::parse(R"({"a": [0,2,3,4,5,6,7,8,9,10,1,11]})"));
  // from == path is allowed and changes nothing
  CHECK(applied(R"({"a": {"b": 1}})", R"([{"op": "move", "from": "/a/b", "path": "/a/b"}])") ==
        json::parse(R"({"a": {"b": 1}})"));
  // "-" of a move target is the end after removing the source
  CHECK(applied(R"({"a": [1, 2, 3]})", R"([{"op": "move", "from": "/a/0", "path": "/a/-"}])") ==
        json::parse(R"({"a": [2, 3, 1]})"));
  // index == size is allowed for add, not for replace or remove; no leading zero; no "-" there
  CHECK(applied(R"({"a": [1]})", R"([{"op": "add", "path": "/a/1", "value": 2}])") ==
        json::parse(R"({"a": [1, 2]})"));
  CHECK(failure(R"({"a": [1]})", R"([{"op": "replace", "path": "/a/1", "value": 2}])").code() ==
        "edit.path_not_found");
  // A token of the wrong form is an invalid pointer (ADR 0024 point 4), an index past the end
  // does not exist.
  CHECK(failure(R"({"a": [1, 2]})", R"([{"op": "remove", "path": "/a/01"}])").code() ==
        "edit.patch_invalid");
  CHECK(failure(R"({"a": [1, 2]})", R"([{"op": "remove", "path": "/a/x"}])").code() ==
        "edit.patch_invalid");
  CHECK(failure(R"({"a": [1]})", R"([{"op": "remove", "path": "/a/-"}])").code() ==
        "edit.patch_invalid");
  CHECK(failure(R"({"a": [[1]]})", R"([{"op": "add", "path": "/a/-/0", "value": 2}])").code() ==
        "edit.patch_invalid");
  CHECK(failure(R"({"a": [1]})", R"([{"op": "add", "path": "/a/2", "value": 2}])").code() ==
        "edit.path_not_found");
  CHECK(failure(R"({"a": [1]})", R"([{"op": "remove", "path": "/a/99999999999999999999999"}])")
            .code() == "edit.path_not_found");
  // in an object "-" is an ordinary key
  CHECK(applied(R"({"a": {}})", R"([{"op": "add", "path": "/a/-", "value": 1}])") ==
        json::parse(R"({"a": {"-": 1}})"));
  // keys with "~" and "/" in the patch and in its inverse (RFC 6901 Sec. 3)
  CHECK(applied(R"({})", R"([{"op": "add", "path": "/~01", "value": 1},
                             {"op": "add", "path": "/a~1b", "value": 2}])") ==
        json::parse(R"({"~1": 1, "a/b": 2})"));
  CHECK(applied(R"({"~1": 1, "a/b": {"c": 2}})", R"([{"op": "remove", "path": "/~01"},
                  {"op": "move", "from": "/a~1b/c", "path": "/x"}])") ==
        json::parse(R"({"a/b": {}, "x": 2})"));
  // a pointer must start with "/" and "~" must be followed by 0 or 1
  CHECK(failure(R"({"a": 1})", R"([{"op": "remove", "path": "a"}])").code() ==
        "edit.patch_invalid");
  CHECK(failure(R"({"a~2": 1})", R"([{"op": "remove", "path": "/a~2"}])").code() ==
        "edit.patch_invalid");
}

TEST_CASE("test compares by value", "[io][patch]") {
  CHECK(applied(R"({"a": 2})", R"([{"op": "test", "path": "/a", "value": 2.0}])") ==
        json::parse(R"({"a": 2})"));
  CHECK(applied(R"({"a": {"x": 1, "y": [1, 2]}})",
                R"([{"op": "test", "path": "/a", "value": {"y": [1, 2], "x": 1}}])") ==
        json::parse(R"({"a": {"x": 1, "y": [1, 2]}})"));
  CHECK(failure(R"({"a": [1, 2]})", R"([{"op": "test", "path": "/a", "value": [2, 1]}])").code() ==
        "edit.test_failed");
}

TEST_CASE("inverse of move", "[io][patch]") {
  // Onto an existing member whose way back shifts array indices (second review of ADR 0024).
  // After removing /r/0, /r/1 is the former /r/2 and its member b exists.
  CHECK(applied(R"({"r": [{"a": 1}, {"b": 2}, {"b": 3}]})",
                R"([{"op": "move", "from": "/r/0", "path": "/r/1/b"}])") ==
        json::parse(R"({"r": [{"b": 2}, {"b": {"a": 1}}]})"));
  // Into an ancestor of the source (RFC-conformant): a plain "move back" would break the prefix
  // rule.
  CHECK(applied(R"({"s": {"radius": {"value": 50.0, "variable": true}}})",
                R"([{"op": "move", "from": "/s/radius/value", "path": "/s/radius"}])") ==
        json::parse(R"({"s": {"radius": 50.0}})"));
  // Onto an existing object member
  CHECK(applied(R"({"a": {"p": 1}, "b": {"p": 2}})",
                R"([{"op": "move", "from": "/a/p", "path": "/b/p"}])") ==
        json::parse(R"({"a": {}, "b": {"p": 1}})"));
}

TEST_CASE("a patch is applied completely or not at all", "[io][patch]") {
  const json doc = json::parse(R"({"a": 1, "b": [1]})");
  try {
    (void)apply_json_patch(doc, json::parse(R"([{"op": "replace", "path": "/a", "value": 2},
                                               {"op": "remove", "path": "/b/5"}])"));
    FAIL("no PatchFailure");
  } catch (const PatchFailure& e) {
    CHECK(e.op_index() == 1);
    CHECK(e.pointer() == "/b/5");
  }
  CHECK(doc == json::parse(R"({"a": 1, "b": [1]})"));
  // An inverse of several operations restores the original.
  CHECK(applied(R"({"a": 1, "b": [1, 2], "c": {"d": 3}})",
                R"([{"op": "replace", "path": "/a", "value": 5},
                    {"op": "add", "path": "/b/0", "value": 0},
                    {"op": "move", "from": "/c/d", "path": "/b/-"},
                    {"op": "copy", "from": "/b", "path": "/c/e"},
                    {"op": "remove", "path": "/a"}])") ==
        json::parse(R"({"b": [0, 1, 2, 3], "c": {"e": [0, 1, 2, 3]}})"));
}

TEST_CASE("the root and invalid patch documents", "[io][patch]") {
  // RFC 6902 Sec. 4.1: an add at the root replaces the whole document.
  CHECK(applied(R"({"a": 1})", R"([{"op": "add", "path": "", "value": [1]}])") ==
        json::parse("[1]"));
  CHECK(failure(R"({"a": 1})", R"([{"op": "remove", "path": ""}])").code() == "edit.patch_invalid");
  const PatchFailure not_array = failure(R"({"a": 1})", R"({"op": "remove", "path": "/a"})");
  CHECK(not_array.code() == "edit.patch_invalid");
  CHECK(!not_array.op_index().has_value());
  CHECK(failure(R"({"a": 1})", R"([{"op": "delete", "path": "/a"}])").code() ==
        "edit.patch_invalid");
  CHECK(failure(R"({"a": 1})", R"([{"op": "add", "path": "/b"}])").code() ==
        "edit.patch_invalid");  // no value
  CHECK(failure(R"({"a": 1})", R"([{"op": "move", "path": "/b"}])").code() ==
        "edit.patch_invalid");  // no from
  CHECK(failure(R"({"a": 1})", R"([{"path": "/a"}])").code() == "edit.patch_invalid");
  CHECK(failure(R"({"a": 1})", R"([7])").code() == "edit.patch_invalid");
}
