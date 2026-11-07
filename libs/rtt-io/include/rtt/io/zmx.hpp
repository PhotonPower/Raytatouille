#pragma once

/// @file zmx.hpp
/// Import of the sequential part of an OpticStudio file (.zmx) into the model (ADR 0044).
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

}  // namespace rtt::io
