#include "rtt/io/edit.hpp"

#include <cstddef>
#include <map>
#include <nlohmann/json.hpp>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "json_format.hpp"
#include "json_patch.hpp"
#include "json_tree.hpp"
#include "rtt/diagnostics/codes.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/json/strict.hpp"

namespace rtt::io {

EditError::EditError(std::string code,
                     std::string location,
                     std::optional<std::size_t> op_index,
                     std::vector<model::Diagnostic> diagnostics,
                     const std::string& message)
    : std::invalid_argument("[" + code + "] " + (location.empty() ? std::string("/") : location) +
                            ": " + message),
      code_(std::move(code)),
      location_(std::move(location)),
      op_index_(op_index),
      diagnostics_(std::move(diagnostics)) {}

namespace {

using Json = nlohmann::json;

/// Codes of the edit functions, checked against the registry at compile time (ADR 0022).
constexpr diagnostics::DiagnosticCode kPatchInvalid{"edit.patch_invalid"};
constexpr diagnostics::DiagnosticCode kReadOnly{"edit.read_only"};
constexpr diagnostics::DiagnosticCode kInvalidValue{"edit.invalid_value"};
constexpr diagnostics::DiagnosticCode kBaseNotRepresentable{"edit.base_not_representable"};
// edit.path_not_found and edit.test_failed come from the engine (json_patch.cpp); registered too.
constexpr diagnostics::DiagnosticCode kPathNotFound{"edit.path_not_found"};
constexpr diagnostics::DiagnosticCode kTestFailed{"edit.test_failed"};

[[noreturn]] void fail(diagnostics::DiagnosticCode code,
                       std::string location,
                       std::optional<std::size_t> op_index,
                       const std::string& message) {
  throw EditError(std::string(code.str()), std::move(location), op_index, {}, message);
}

/// True if the pointer is "", or on or below /schema_version or /units (ADR 0024 point 4).
bool read_only(const std::string& pointer) {
  if (pointer.empty()) return true;
  for (const std::string_view top :
       {std::string_view("/schema_version"), std::string_view("/units")}) {
    if (pointer.starts_with(top) && (pointer.size() == top.size() || pointer[top.size()] == '/')) {
      return true;
    }
  }
  return false;
}

/// Rejects writes to read-only places before the patch is applied.
void check_read_only(const Json& patch) {
  if (!patch.is_array()) return;  // reported by the engine
  for (std::size_t i = 0; i < patch.size(); ++i) {
    const Json& op = patch[i];
    if (!op.is_object()) continue;
    const auto name = op.find("op");
    if (name == op.end() || !name->is_string() || *name == "test") continue;
    for (const char* key : {"path", "from"}) {
      if (std::string_view(key) == "from" && *name != "move") continue;  // copy only reads
      const auto it = op.find(key);
      if (it != op.end() && it->is_string() && read_only(it->get<std::string>())) {
        fail(kReadOnly, it->get<std::string>(), i,
             "operation " + std::to_string(i) + ": '" + it->get<std::string>() +
                 "' cannot be changed (schema version, units, whole document)");
      }
    }
  }
}

std::map<std::string, std::size_t> error_counts(const std::vector<model::Diagnostic>& d) {
  std::map<std::string, std::size_t> counts;
  for (const model::Diagnostic& x : d) {
    if (x.severity == model::Severity::Error) ++counts[x.code];
  }
  return counts;
}

}  // namespace

std::string to_edit_json(const model::System& system) {
  try {
    return detail::format_canonical(detail::write_system(system, detail::Form::Edit));
  } catch (const nlohmann::json::exception& e) {  // e.g. a string that is not valid UTF-8
    throw std::invalid_argument(std::string("the system cannot be written: ") + e.what());
  }
}

PatchResult apply_patch(const model::System& system, std::string_view patch_json, EditCheck check) {
  Json patch;
  try {
    patch = rtt::json::parse_strict(patch_json);  // duplicate keys are errors (RFC 6902 A.13)
  } catch (const rtt::json::StrictParseError& e) {
    // The place is inside the patch document, not the system: keep it in the message only.
    fail(kPatchInvalid, "", std::nullopt,
         "invalid patch document" +
             (e.pointer().empty() ? std::string() : " at '" + e.pointer() + "'") + ": " +
             e.message());
  }
  check_read_only(patch);

  Json document;
  try {
    document = detail::plain(detail::write_system(system, detail::Form::Edit));
  } catch (const std::invalid_argument& e) {
    fail(kBaseNotRepresentable, "", std::nullopt,
         std::string("the system cannot be written: ") + e.what());
  }

  detail::PatchApplied applied;
  try {
    applied = detail::apply_json_patch(document, patch);
  } catch (const detail::PatchFailure& e) {
    // The engine uses only edit.patch_invalid, edit.path_not_found and edit.test_failed.
    const diagnostics::DiagnosticCode code = e.code() == kPathNotFound.str() ? kPathNotFound
                                             : e.code() == kTestFailed.str() ? kTestFailed
                                                                             : kPatchInvalid;
    fail(code, e.pointer(), e.op_index(), e.message());
  }

  PatchResult result;
  try {
    result.system = detail::read_system(applied.document);
  } catch (const ParseError& e) {
    // what() is "<pointer>: <message>"; keep the message only.
    const std::string text = e.what();
    const std::string prefix = (e.pointer().empty() ? std::string("/") : e.pointer()) + ": ";
    fail(kInvalidValue, e.pointer(), std::nullopt,
         "the result is not a valid system: " +
             (text.starts_with(prefix) ? text.substr(prefix.size()) : text));
  }
  result.diagnostics = model::validate(result.system);

  if (check == EditCheck::NoNewErrors) {
    const std::map<std::string, std::size_t> before = error_counts(model::validate(system));
    const std::map<std::string, std::size_t> after = error_counts(result.diagnostics);
    std::vector<model::Diagnostic> grown;
    for (const model::Diagnostic& d : result.diagnostics) {
      if (d.severity != model::Severity::Error) continue;
      const auto b = before.find(d.code);
      if (after.at(d.code) > (b == before.end() ? 0 : b->second)) grown.push_back(d);
    }
    if (!grown.empty()) {
      const model::Diagnostic first = grown.front();
      throw EditError(first.code, first.location, std::nullopt, std::move(grown),
                      "the patch adds an error: " + first.message);
    }
  }
  try {
    result.inverse = applied.inverse.dump();
  } catch (const nlohmann::json::exception& e) {  // a string of the system that is not UTF-8
    fail(kBaseNotRepresentable, "", std::nullopt,
         std::string("the inverse cannot be written: ") + e.what());
  }
  return result;
}

}  // namespace rtt::io
