#pragma once

/// @file zmx.hpp
/// Import of the sequential part of an OpticStudio file (.zmx) into the model (ADR 0044), and
/// the export of a coaxial system into such a file (ADR 0045).
///
/// The meanings of the quantities come from the OpticStudio User Guide, the mapping of the
/// keywords from readers with a free licence (docs/quellen.md, rows "ZMX-Import"); readings
/// that only one reader or only the test corpus supports are marked as such in ADR 0044.
/// Units: the lengths of the file are scaled by its UNIT into mm; angles stay in degree,
/// wavelengths in um. OpticStudio indices are relative to air, the model's are absolute:
/// wavelengths are converted from air to vacuum (compile::line_vacuum_wavelength_um), an F/# is
/// divided and an object-space NA multiplied by n_air(lambda_vac,ref, 20 degC, 1 atm) after
/// Ciddor (ADR 0044, point 4). Coordinates: the vertex of Zemax surface 1 is the global origin,
/// +z the optical axis; radii keep the sign of OpticStudio, which is the model's.

#include <cstddef>
#include <filesystem>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "rtt/model/system.hpp"
#include "rtt/model/validate.hpp"

namespace rtt::material {
class MaterialLibrary;
}

namespace rtt::io {

/// Result of an import: the system and the warnings (codes zmx.*, plus the warnings of validate
/// and compile with their own codes). The location of a zmx.* warning is "line <n>" (1-based)
/// or "file" for the whole file.
struct ZmxImport {
  model::System system;
  std::vector<model::Diagnostic> warnings;
};

/// A file that cannot be imported. code() is a registered zmx.* code (docs/diagnostics.md),
/// line() the 1-based line of the file, 0 for the whole file (e.g. zmx.invalid_result).
class ZmxError : public std::runtime_error {
 public:
  /// @param code    registered zmx.* code
  /// @param line    1-based line, 0 for the whole file
  /// @param message human-readable text; what() is "line <n>: <message> [<code>]"
  ZmxError(std::string code, std::size_t line, const std::string& message);

  [[nodiscard]] const std::string& code() const noexcept { return code_; }
  [[nodiscard]] std::size_t line() const noexcept { return line_; }

 private:
  std::string code_;
  std::size_t line_;
};

/// Imports the bytes of a .zmx file (UTF-16 LE with BOM FF FE, otherwise ASCII or UTF-8).
/// Glass names are resolved through the GCAT catalogues in `materials` (ADR 0044, point 7).
/// The result has passed model::validate and compile::compile; it is written with a single
/// automatic path, AIR at 20 degC and 1 atm, and element and surface names after ADR 0044,
/// point 5.
/// @throws ZmxError for every rule of ADR 0044 that rejects the file, and with
///         zmx.invalid_result if validate or compile reject the result
[[nodiscard]] ZmxImport import_zmx(std::string_view bytes,
                                   const material::MaterialLibrary& materials);

/// Reads `file` and imports it as import_zmx().
/// @throws ZmxError as import_zmx(); std::runtime_error if the file cannot be read
[[nodiscard]] ZmxImport load_zmx(const std::filesystem::path& file,
                                 const material::MaterialLibrary& materials);

/// What to export (ADR 0045, point 3; the CSV surface table uses the same options).
struct ExportOptions {
  /// Name of the path to export; none: path 0. The other paths are dropped with the warning
  /// export.paths_dropped.
  std::optional<std::string> path;
  /// Index of the configuration (ADR 0029), exported with the parameters resolved for it; none:
  /// configuration 0, an error (export.multi_config) for a system with more than one.
  std::optional<std::size_t> configuration;
};

/// Result of an export: the file text (UTF-8, CRLF) and the warnings (codes export.*; the
/// location is the JSON pointer of the surface, e.g. "/root/children/1/surfaces/0", or of the
/// part of the system, e.g. "/paths" or "/name").
struct ZmxExport {
  std::string text;
  std::vector<model::Diagnostic> warnings;
};

/// A system that cannot be exported. code() is a registered export.* code (docs/diagnostics.md),
/// surface() the id of the surface concerned, empty if the cause is not one surface.
class ExportError : public std::runtime_error {
 public:
  /// @param code    registered export.* code
  /// @param surface id of the surface concerned, empty for the whole system
  /// @param message human-readable text; what() is "<surface>: <message> [<code>]", without
  ///                "<surface>: " for the whole system
  ExportError(std::string code, std::string surface, const std::string& message);

  [[nodiscard]] const std::string& code() const noexcept { return code_; }
  [[nodiscard]] const std::string& surface() const noexcept { return surface_; }

 private:
  std::string code_;
  std::string surface_;
};

/// Writes one path of a coaxial system as a sequential .zmx file (ADR 0045, points 1 to 4): the
/// inverse of import_zmx, so that the import reads every exported file. Coaxial: every vertex on
/// the z axis and every surface axis +-e_z, with local x and y parallel to +-e_x and +-e_y, to
/// 1e-12 (mm, or on the unit vectors); a surface with the axis -e_z is written mirrored
/// (R -> -R, A_i -> -A_i). Lengths in mm (UNIT MM); radii with the sign of the model, which is
/// OpticStudio's; thicknesses along the global z, negative after an odd number of mirrors (the
/// frame of OpticStudio does not turn at a mirror). OpticStudio indices are relative to air:
/// wavelengths are written as lambda_air = lambda_vac / n_air(lambda_vac), an F/# multiplied and
/// an object-space NA divided by n_air(lambda_vac,ref), n_air after Ciddor at 20 degC and 1 atm.
/// Numbers in the shortest form that reads back to the same double (std::to_chars).
/// @param system    the system; compiled with `materials` and without coatings (a coating is
///                  dropped with a warning)
/// @param materials the glass catalogues; the GCAT line lists the catalogues of the glasses used
/// @throws ExportError for every rule of ADR 0045 that rejects the system
/// @throws std::invalid_argument for an unknown path name or a configuration index beyond the
///         configurations of the system, or a system with an inactive node (ADR 0049, point 6;
///         supported later)
/// @throws compile::CompileError if the system does not compile (validate errors, unknown
///         materials)
[[nodiscard]] ZmxExport export_zmx(const model::System& system,
                                   const material::MaterialLibrary& materials,
                                   const ExportOptions& options = {});

/// Exports as export_zmx() and writes the text to `file` (binary, so the CRLF stay).
/// @return the warnings of the export
/// @throws as export_zmx(); std::runtime_error if the file cannot be written
[[nodiscard]] std::vector<model::Diagnostic> save_zmx(const std::filesystem::path& file,
                                                      const model::System& system,
                                                      const material::MaterialLibrary& materials,
                                                      const ExportOptions& options = {});

}  // namespace rtt::io
