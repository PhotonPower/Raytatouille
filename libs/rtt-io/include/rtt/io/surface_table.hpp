#pragma once

/// @file surface_table.hpp
/// The CSV surface table of a coaxial system, in both directions (ADR 0045, point 7).
///
/// One row for the object and one per surface of the exported path, in path order, after a
/// header with the column names: row, surface, element, kind, shape, radius, conic, a4 ... a16,
/// thickness, material, aperture, aperture_1, aperture_2, stop, focal_length, object_distance.
/// Wavelengths, fields, system aperture and environment stay in the JSON file: reading needs a
/// base system. Format after RFC 4180, section 2: comma separated, CRLF written (CRLF and LF
/// read), UTF-8, quotes as in rules 6 and 7. Units: lengths in mm (radius, thickness, apertures,
/// focal length, object distance), asphere coefficients a_n in mm^(1-n); radii with the sign of
/// the model; coordinates global, +z the optical axis, thicknesses along z (negative after an odd
/// number of mirrors). Numbers in the shortest form that reads back to the same double;
/// "Infinity" for an object at infinity.

#include <cstddef>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "rtt/io/zmx.hpp"
#include "rtt/model/system.hpp"
#include "rtt/model/validate.hpp"

namespace rtt::material {
class MaterialLibrary;
}

namespace rtt::io {

/// Result of writing a table: the text (UTF-8, CRLF) and the warnings (codes export.*, located
/// as for ZmxExport).
struct TableExport {
  std::string text;
  std::vector<model::Diagnostic> warnings;
};

/// Result of reading a table: the system and the warnings (table.base_dropped at the pointer of
/// the dropped block, e.g. "/optimization", and the warnings of validate and compile).
struct TableImport {
  model::System system;
  std::vector<model::Diagnostic> warnings;
};

/// A table that cannot be read. code() is a registered table.* code; line() is the 1-based line
/// of the file where the field starts (the header is line 1), column() the 1-based number of the
/// field in its record; both 0 for the whole file (a missing column, table.invalid_result).
class TableError : public std::runtime_error {
 public:
  /// @param code    registered table.* code
  /// @param line    1-based line, 0 for the whole file
  /// @param column  1-based field number, 0 for the whole record or file
  /// @param message human-readable text; what() is "line <l>, column <c>: <message> [<code>]",
  ///                without the position for line 0
  TableError(std::string code, std::size_t line, std::size_t column, const std::string& message);

  [[nodiscard]] const std::string& code() const noexcept { return code_; }
  [[nodiscard]] std::size_t line() const noexcept { return line_; }
  [[nodiscard]] std::size_t column() const noexcept { return column_; }

 private:
  std::string code_;
  std::size_t line_;
  std::size_t column_;
};

/// Writes one path of one configuration of a coaxial system as a surface table (ADR 0045,
/// point 7): the subset and the errors of the ZMX export (points 2 and 3, the same export.*
/// codes), shapes plane, conic and even_asphere up to A16 (others: export.surface_unsupported),
/// the ideal lens in focal_length and object_distance. Unlike the ZMX export, CONST: and VACUUM
/// materials, any environment and a design conjugate are allowed. A surface with the axis -e_z
/// is written mirrored (R -> -R, A_i -> -A_i).
/// @throws ExportError, std::invalid_argument and compile::CompileError as export_zmx()
[[nodiscard]] TableExport write_surface_table(const model::System& system,
                                              const material::MaterialLibrary& materials,
                                              const ExportOptions& options = {});

/// Reads a surface table into the tree, object and stop of `base` (ADR 0045, point 7) in the
/// canonical form of point 6: a flat tree, the first element at z = 0, every further element and
/// inner surface relative_to_preceding with (0, 0, t), one automatic path (named as the first
/// path of `base`, else "main"). The blocks of `base` that refer to its tree (further paths,
/// optimization, tolerances, parameter table, configurations) are dropped with
/// table.base_dropped; Params of `base` bound to the table keep their value of configuration 0.
/// @param csv  the bytes of the file; a UTF-8 byte order mark is skipped
/// @throws TableError with table.parse, table.column, table.value or table.duplicate at the
///         field, and table.invalid_result if validate or compile reject the result
[[nodiscard]] TableImport read_surface_table(std::string_view csv,
                                             model::System base,
                                             const material::MaterialLibrary& materials);

}  // namespace rtt::io
