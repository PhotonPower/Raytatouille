#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <string_view>

#include "rtt/diagnostics/codes.hpp"

using rtt::diagnostics::CodeInfo;
using rtt::diagnostics::DiagnosticCode;
using rtt::diagnostics::find_code;
using rtt::diagnostics::kCodes;
using rtt::diagnostics::Severity;

namespace {

/// "<group>.<what>": lower-case letters, digits and '_', exactly one dot, no empty part.
bool well_formed(std::string_view code) {
  const auto dot = code.find('.');
  if (dot == std::string_view::npos || dot == 0 || dot + 1 == code.size()) return false;
  if (code.find('.', dot + 1) != std::string_view::npos) return false;
  return std::all_of(code.begin(), code.end(), [](char c) {
    return (c >= 'a' && c <= 'z') || (c >= '0' && c <= '9') || c == '_' || c == '.';
  });
}

// Checked at compile time: a registered code converts and carries its registry severity. An
// unregistered code does not compile (consteval constructor).
static_assert(DiagnosticCode("material.unknown").severity() == Severity::Error);
static_assert(DiagnosticCode("aperture.na_not_physical").severity() == Severity::Warning);
static_assert(DiagnosticCode("material.unknown").str() == "material.unknown");

}  // namespace

TEST_CASE("the registry is sorted, unique and well formed", "[diagnostics]") {
  const auto by_code = [](const CodeInfo& a, const CodeInfo& b) { return a.code < b.code; };
  REQUIRE(std::is_sorted(kCodes.begin(), kCodes.end(), by_code));
  const auto same = [](const CodeInfo& a, const CodeInfo& b) { return a.code == b.code; };
  REQUIRE(std::adjacent_find(kCodes.begin(), kCodes.end(), same) == kCodes.end());
  for (const CodeInfo& info : kCodes) {
    INFO(info.code);
    REQUIRE(well_formed(info.code));
    REQUIRE_FALSE(info.summary.empty());
  }
}

TEST_CASE("find_code returns the registry entry or nullptr", "[diagnostics]") {
  const CodeInfo* info = find_code("paths.unknown_surface");
  REQUIRE(info != nullptr);
  REQUIRE(info->code == "paths.unknown_surface");
  REQUIRE(info->severity == Severity::Error);
  REQUIRE(find_code("paths.unknown") == nullptr);
  REQUIRE(find_code("") == nullptr);
  REQUIRE(find_code("PATHS.UNKNOWN_SURFACE") == nullptr);  // case-sensitive
}

TEST_CASE("every code names its producer", "[diagnostics]") {
  for (const CodeInfo& info : kCodes) {
    INFO(info.code);
    REQUIRE((info.producer == "validate" || info.producer == "compile" ||
             info.producer == "analysis" || info.producer == "agf" || info.producer == "edit" ||
             info.producer == "io"));
  }
  REQUIRE(find_code("stop.not_on_path")->producer == "compile");
  REQUIRE(find_code("rays.lost")->producer == "analysis");
  REQUIRE(find_code("agf.duplicate_glass")->producer == "agf");
  REQUIRE(find_code("edit.patch_invalid")->producer == "edit");
  REQUIRE(find_code("io.pickup_dropped")->producer == "io");
}
