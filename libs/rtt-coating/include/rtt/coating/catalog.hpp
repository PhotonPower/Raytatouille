#pragma once

/// @file catalog.hpp
/// Coating designs from JSON catalogue files and the CoatingLibrary that resolves the model's
/// CoatingRef names (ADR 0019).
///
/// A catalogue file (schema/raytatouille-coatings.schema.json):
///
///     {
///       "format": "raytatouille-coatings", "schema_version": "0.1.0", "catalog": "DEMO",
///       "coatings": [
///         {"name": "AR_MGF2", "description": "...", "design_wavelength_um": 0.55,
///          "layers": [{"material": "CONST:1.38", "qwot": 1.0}]}
///       ]
///     }
///
/// - Layers are listed from the ambient side to the substrate (transfer_matrix.hpp). Each layer
///   has a material reference (resolved by rtt::material::MaterialLibrary when the system is
///   compiled) and exactly one of "thickness_um" (physical thickness, um) or "qwot" (number of
///   quarter waves at the design wavelength, thickness.hpp).
/// - "design_wavelength_um" (vacuum, um) is required if a layer uses "qwot".
/// - A design is referenced as "CATALOG:NAME", e.g. "DEMO:AR_MGF2". Catalogue and coating names
///   are not empty and contain no ':'.
/// - Strict parser as for system files (ADR 0008): unknown keys and wrong types are errors with
///   the file and a JSON pointer.

#include <filesystem>
#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "rtt/coating/thickness.hpp"

namespace rtt::coating {

/// Value of "format" in a coating catalogue.
inline constexpr std::string_view kCatalogFormat = "raytatouille-coatings";
/// Version of the catalogue format read and described by this library.
inline constexpr std::string_view kCatalogSchemaVersion = "0.1.0";

/// One layer of a coating design.
struct LayerSpec {
  std::string material;  ///< material reference, e.g. "CONST:1.38" or "SCHOTT:N-BK7"
  Thickness thickness;   ///< physical (um) or quarter waves at the design wavelength
  bool operator==(const LayerSpec&) const = default;
};

/// A coating design: layers from the ambient side to the substrate.
struct CoatingDesign {
  std::string name;                            ///< name within the catalogue
  std::string description;                     ///< free text, may be empty
  std::optional<double> design_wavelength_um;  ///< vacuum, um; required with QWOT layers
  std::vector<LayerSpec> layers;               ///< at least one
  bool operator==(const CoatingDesign&) const = default;
};

/// Contents of one catalogue file.
struct CoatingCatalog {
  std::string name;  ///< catalogue name, the prefix of references ("DEMO" in "DEMO:AR_MGF2")
  std::vector<CoatingDesign> coatings;
};

/// Error in a coating catalogue. what() names the file and the JSON pointer.
class CoatingCatalogError : public std::runtime_error {
 public:
  /// @param file    file name (may be empty for text without file)
  /// @param pointer JSON pointer to the offending value, empty for the whole file
  /// @param message description of the error
  CoatingCatalogError(std::string file, std::string pointer, const std::string& message);

  [[nodiscard]] const std::string& file() const noexcept { return file_; }
  /// JSON pointer to the offending value; empty for the whole file (I/O, syntax).
  [[nodiscard]] const std::string& pointer() const noexcept { return pointer_; }

 private:
  std::string file_;
  std::string pointer_;
};

/// A reference that the CoatingLibrary cannot resolve.
class UnknownCoating : public std::runtime_error {
 public:
  using std::runtime_error::runtime_error;
};

/// Parses a catalogue from JSON text.
/// @param json_text contents of a catalogue file
/// @param file      file name for error messages
/// @throws CoatingCatalogError for syntax errors, unknown or missing keys, wrong types, a wrong
///         format or schema version, empty or duplicate names, names with ':', an empty layer
///         list, a layer without or with both thickness keys, negative or non-finite values, or
///         a QWOT layer without design wavelength
[[nodiscard]] CoatingCatalog parse_coating_catalog(std::string_view json_text,
                                                   const std::string& file = "");

/// Reads a catalogue file.
/// @throws CoatingCatalogError as parse_coating_catalog(), or if the file cannot be read
[[nodiscard]] CoatingCatalog load_coating_catalog(const std::filesystem::path& file);

/// Resolves coating references "CATALOG:NAME" to designs. Thread-safe.
class CoatingLibrary {
 public:
  CoatingLibrary() = default;

  /// Registers the designs of a catalogue under "CATALOG:NAME". A catalogue built in code is
  /// checked like a parsed file.
  /// @throws std::invalid_argument if the catalogue name is empty, contains ':' or is already
  ///         registered, or a design breaks a rule of the file format (name, duplicate name,
  ///         no layers, empty material, negative or non-finite thickness or count, design
  ///         wavelength not > 0); nothing is registered then
  void add(const CoatingCatalog& catalog);

  /// Loads catalogue files: a .json file, or a directory whose *.json files (case-insensitive,
  /// not recursive) are loaded in sorted order. All or nothing: on an error nothing is added.
  /// @throws CoatingCatalogError for a malformed file
  /// @throws std::invalid_argument if the path has no .json file or a catalogue name is already
  ///         registered (also twice in one directory)
  void add_catalog(const std::filesystem::path& path);

  /// The design for a reference, never null; kept alive independently of the library.
  /// @throws UnknownCoating if the reference is malformed or unknown
  [[nodiscard]] std::shared_ptr<const CoatingDesign> resolve(std::string_view reference) const;

  /// The design for a reference, or null if it is unknown (no exception).
  [[nodiscard]] std::shared_ptr<const CoatingDesign> find(std::string_view reference) const;

 private:
  mutable std::mutex mutex_;
  std::set<std::string, std::less<>> catalogs_;
  std::map<std::string, std::shared_ptr<const CoatingDesign>, std::less<>> designs_;
};

}  // namespace rtt::coating
