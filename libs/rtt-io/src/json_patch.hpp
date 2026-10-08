#pragma once

/// @file json_patch.hpp
/// JSON Patch (RFC 6902) on generic JSON with inverse patches (ADR 0024 point 3). Internal to
/// rtt-io: the public API is rtt/io/edit.hpp. nlohmann::json::patch is not used because it
/// deviates from RFC 6902 (an add, copy or move into a primitive parent hits an assertion, a
/// remove under a primitive parent does nothing, a move into its own child is not rejected).

#include <cstddef>
#include <nlohmann/json.hpp>
#include <optional>
#include <stdexcept>
#include <string>

namespace rtt::io::detail {

/// Failure of a JSON Patch: a stable code of the group edit.* (edit.patch_invalid,
/// edit.path_not_found, edit.test_failed), the JSON pointer of the place (RFC 6901, "" for the
/// whole patch or document) and the index of the failing operation (none if the patch document
/// itself is not an array).
class PatchFailure : public std::runtime_error {
 public:
  PatchFailure(std::string code,
               std::string pointer,
               std::optional<std::size_t> op_index,
               const std::string& message);
  [[nodiscard]] const std::string& code() const noexcept { return code_; }
  [[nodiscard]] const std::string& pointer() const noexcept { return pointer_; }
  [[nodiscard]] std::optional<std::size_t> op_index() const noexcept { return op_index_; }
  /// Description without the pointer.
  [[nodiscard]] const std::string& message() const noexcept { return message_; }

 private:
  std::string code_;
  std::string pointer_;
  std::optional<std::size_t> op_index_;
  std::string message_;
};

/// Result of a patch: the new document and the inverse patch document. Applying `inverse` to
/// `document` gives back the original document.
struct PatchApplied {
  nlohmann::json document;
  nlohmann::json inverse;
};

/// Applies the RFC 6902 patch document `patch` (an array of operations) to a copy of `document`.
///
/// RFC 6902 Sec. 4.1-4.6 (operations) and Sec. 5 (the patch is applied completely or not at all):
/// - add, remove, replace, move, copy, test; members other than op, path, from and value are
///   ignored (Appendix A.11);
/// - pointers per RFC 6901: "~1" decodes to "/", then "~0" to "~" (so "~01" is "~1", A.14);
/// - an array index is "0" or digits without a leading zero, another token is an invalid pointer
///   (edit.patch_invalid); add and the targets of move and copy accept index <= size and "-"
///   (append), the other places index < size (else edit.path_not_found); the target index of
///   move is checked against the array after removing the source (Sec. 4.4);
/// - the "from" of a move must not be a proper prefix of its "path" (compared token by token;
///   from == path is allowed and changes nothing);
/// - test compares by value: object members in any order, numbers by numeric value (Sec. 4.6).
///
/// The inverse is built per operation against the state before it, the operations of the
/// inverse in reverse order. A move from F to P is decomposed like Sec. 4.4 (remove F, then
/// add P): its inverse is "remove P" (or "replace P" with its value after removing F, if P was an
/// existing object member), then "add F" with the moved value.
///
/// @throws PatchFailure for an invalid patch document, operation, pointer or array index token
///         (edit.patch_invalid), a missing target, source or parent, an index out of range, or a
///         parent that is no object or array (edit.path_not_found), and a failed test
///         (edit.test_failed); `document` is never changed
[[nodiscard]] PatchApplied apply_json_patch(const nlohmann::json& document,
                                            const nlohmann::json& patch);

}  // namespace rtt::io::detail
