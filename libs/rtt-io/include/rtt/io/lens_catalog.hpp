#pragma once
/// @file lens_catalog.hpp
/// Lens catalogues as JSON files (ADR 0053): stock singlets, cemented doublets and triplets that
/// the design helpers of the Python API (rt.design) turn into systems.
///
/// A catalogue file (schema/raytatouille-lenses.schema.json):
///
///     {
///       "format": "raytatouille-lenses", "schema_version": "0.1.0", "catalog": "DEMO",
///       "lenses": [
///         {"name": "PCX-25-100", "kind": "singlet", "description": "...",
///          "radii": [51.68, null], "thicknesses": [4.0], "materials": ["CONST:1.5168"],
///          "diameter": 25.4}
///       ]
///     }
///
/// - Lengths in mm. Radii with the sign convention of the system format (centre of curvature in
///   +z: positive), null is plano, 0 is an error. Thicknesses > 0, gaps >= 0.
/// - Counts per kind: radii 2 / 3 / 6, thicknesses 1 / 2 / 3, materials 1 / 2 / 3 for singlet /
///   doublet (cemented) / triplet; "gaps" has exactly 2 entries for a triplet and is absent
///   otherwise.
/// - Materials are non-empty references CATALOG:NAME, resolved only when a system is compiled.
/// - "diameter" (optional, > 0) is the clear diameter of every surface.
/// - Catalogue and lens names are not empty and contain no ':'; lens names are unique.
/// - Strict parser as for system and coating files (ADR 0008, 0020): unknown keys, wrong types and
///   wrong counts are errors with the file and a JSON pointer.

#include <cstdint>
#include <filesystem>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace rtt::io {

/// Value of "format" in a lens catalogue.
inline constexpr std::string_view kLensCatalogFormat = "raytatouille-lenses";
/// Version of the lens catalogue format read by this library (ADR 0053).
inline constexpr std::string_view kLensCatalogSchemaVersion = "0.1.0";

/// Kind of a catalogue lens.
enum class LensKind : std::uint8_t {
  Singlet,  ///< one lens, two surfaces
  Doublet,  ///< cemented doublet, three surfaces
  Triplet,  ///< three singlets with two air gaps
};

/// One lens of a catalogue. All lengths in mm.
struct LensEntry {
  std::string name;                          ///< unique within the catalogue, no ':'
  std::string description;                   ///< free text, may be empty
  LensKind kind = LensKind::Singlet;         ///< determines the counts below
  std::vector<std::optional<double>> radii;  ///< vertex radii in order along +z; nullopt = plano
  std::vector<double> thicknesses;           ///< centre thicknesses of the glass segments, > 0
  std::vector<double> gaps;                  ///< triplet only: air gaps between the lenses, >= 0
  std::vector<std::string> materials;        ///< one reference CATALOG:NAME per glass segment
  std::optional<double> diameter;            ///< clear diameter of every surface, > 0

  bool operator==(const LensEntry&) const = default;
};

/// Contents of one lens catalogue file.
struct LensCatalog {
  std::string name;               ///< catalogue name
  std::vector<LensEntry> lenses;  ///< in file order

  bool operator==(const LensCatalog&) const = default;
};

/// Error in a lens catalogue. what() names the file and the JSON pointer.
class LensCatalogError : public std::runtime_error {
 public:
  /// @param file    file name (may be empty for text without file)
  /// @param pointer JSON pointer to the offending value, empty for the whole file
  /// @param message description of the error
  LensCatalogError(std::string file, std::string pointer, const std::string& message);

  /// File name as given to the parser; empty for text without file.
  [[nodiscard]] const std::string& file() const noexcept { return file_; }
  /// JSON pointer to the offending value; empty for the whole file (I/O, syntax).
  [[nodiscard]] const std::string& pointer() const noexcept { return pointer_; }

 private:
  std::string file_;
  std::string pointer_;
};

/// Name of a kind as in the file: "singlet", "doublet" or "triplet".
[[nodiscard]] std::string_view lens_kind_name(LensKind kind) noexcept;

/// Parses a lens catalogue from JSON text.
/// @param json_text contents of a catalogue file
/// @param file      file name for the error messages (may be empty)
/// @throws LensCatalogError at the first error, with the JSON pointer
[[nodiscard]] LensCatalog parse_lens_catalog(std::string_view json_text,
                                             const std::string& file = {});

/// Reads a lens catalogue file.
/// @param file path of the catalogue file (any name; the format is recognised by "format")
/// @throws LensCatalogError if the file cannot be read or is malformed
[[nodiscard]] LensCatalog load_lens_catalog(const std::filesystem::path& file);

}  // namespace rtt::io
