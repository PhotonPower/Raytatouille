#pragma once

/// @file edit.hpp
/// Changing a system with JSON Patch (RFC 6902) on its edit form, with inverse patches for undo
/// (ADR 0024). Patches cross this interface as JSON text, so nlohmann stays private (ADR 0020).
///
/// The edit form is the system file with every value written: defaults are not omitted, every
/// unbound Param is an object {"value": ..., "variable": ...} (plus "min"/"max" if set), a bound
/// one {"param": ...} (ADR 0029); every row of the parameter table carries "variable"; optional
/// members without a value (surface aperture, element material, bounds, configurations when
/// there are none) are missing. It has the same
/// schema and schema version as the file; units and signs are those of the file format (mm,
/// um, degree for keys ending in _deg; docs/architecture.md). Every JSON pointer into the
/// canonical file is also valid in the edit form.

#include <cstddef>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "rtt/model/system.hpp"
#include "rtt/model/validate.hpp"

namespace rtt::io {

/// What apply_patch checks after a patch beyond the strict reading of the result.
enum class EditCheck : std::uint8_t {
  /// Reject the patch if, for any diagnostic code, the number of errors of validate() grows
  /// compared with the input system. Systems that already have errors stay editable.
  NoNewErrors,
  /// Only the strict reading; for undo and redo, which return to a state accepted before.
  StructureOnly,
};

/// Result of apply_patch.
struct PatchResult {
  model::System system;  ///< the changed system
  /// RFC 6902 patch document (JSON text) that turns the edit form of `system` back into the
  /// input; apply it with EditCheck::StructureOnly.
  std::string inverse;
  std::vector<model::Diagnostic> diagnostics;  ///< all diagnostics of validate(system)
};

/// A patch that cannot be applied. The input system is never changed.
///
/// code() is a registered code: edit.patch_invalid, edit.path_not_found, edit.read_only,
/// edit.test_failed, edit.invalid_value, edit.base_not_representable, or for EditCheck::NoNewErrors
/// the code of the first new error of validate() (e.g. surface.id_duplicate).
class EditError : public std::invalid_argument {
 public:
  EditError(std::string code,
            std::string location,
            std::optional<std::size_t> op_index,
            std::vector<model::Diagnostic> diagnostics,
            const std::string& message);
  /// Stable code (docs/diagnostics.md).
  [[nodiscard]] const std::string& code() const noexcept { return code_; }
  /// JSON pointer (RFC 6901) into the edit form: the operation's path or from, the place the
  /// strict reader rejected, or the location of the first new error.
  [[nodiscard]] const std::string& location() const noexcept { return location_; }
  /// Index of the failing operation in the patch; empty for errors of the whole patch document
  /// and for new errors of validate().
  [[nodiscard]] std::optional<std::size_t> op_index() const noexcept { return op_index_; }
  /// For new errors of validate(): all error diagnostics of the codes whose number grew.
  [[nodiscard]] const std::vector<model::Diagnostic>& diagnostics() const noexcept {
    return diagnostics_;
  }

 private:
  std::string code_;
  std::string location_;
  std::optional<std::size_t> op_index_;
  std::vector<model::Diagnostic> diagnostics_;
};

/// The edit form of `system` as JSON text (key order and indentation as to_json).
/// @throws std::invalid_argument if the system cannot be written (a non-finite number, a string
///         that is not valid UTF-8, or material and segment_materials both set)
[[nodiscard]] std::string to_edit_json(const model::System& system);

/// Applies an RFC 6902 patch document (JSON text) to the edit form of `system` (ADR 0024).
///
/// Steps, each failing with the error named: parse the patch (edit.patch_invalid); reject
/// operations other than test that write on or below /schema_version or /units, or the whole
/// document, checked for all operations before any is applied (edit.read_only); write the edit
/// form (edit.base_not_representable); apply the operations in order (RFC 6902 Sec. 4.1-4.6,
/// all or nothing per Sec. 5; edit.patch_invalid, edit.path_not_found, edit.test_failed); read
/// the result strictly (edit.invalid_value); run validate() and the check. The inverse is built
/// per operation; a move is decomposed into remove and add (Sec. 4.4).
///
/// @param system     input system, never changed
/// @param patch_json RFC 6902 patch document: a JSON array of operations, duplicate keys are an
///                   error (ADR 0020)
/// @param check      NoNewErrors for commands, StructureOnly for undo and redo
/// @throws EditError if the patch cannot be applied or the check fails. A string of the input
///         that is not valid UTF-8 (only possible from C++) is reported only if it ends up in
///         the inverse (edit.base_not_representable); otherwise writing the result fails later.
[[nodiscard]] PatchResult apply_patch(const model::System& system,
                                      std::string_view patch_json,
                                      EditCheck check = EditCheck::NoNewErrors);

}  // namespace rtt::io
