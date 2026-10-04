#pragma once

/// @file material.hpp
/// Optical media and the resolution of material references (docs/architecture.md, "Medien").
///
/// Conventions (docs/architecture.md, "Konventionen"):
/// - wavelengths are vacuum wavelengths in micrometre (um),
/// - temperatures are in degree Celsius, pressures in atm,
/// - the complex refractive index is n_c = n + i*kappa with kappa >= 0 for absorption,
///   matching plane waves ~ exp(i(k.r - omega t)),
/// - indices are absolute (against vacuum); catalogue data relative to air are converted with
///   Ciddor air (follows with #25, until then formula values are used unchanged).

#include <filesystem>
#include <functional>
#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>

#include "rtt/math/types.hpp"

namespace rtt::material {

/// Closed interval of vacuum wavelengths in micrometre on which a material is defined.
struct WavelengthRange {
  double min_um = 0.0;
  double max_um = 0.0;

  /// True if min_um <= wavelength_um <= max_um.
  [[nodiscard]] bool contains(double wavelength_um) const noexcept {
    return min_um <= wavelength_um && wavelength_um <= max_um;
  }

  bool operator==(const WavelengthRange&) const = default;
};

/// An optical medium. Implementations are immutable and safe to share between threads.
class Material {
 public:
  Material() = default;
  Material(const Material&) = delete;
  Material& operator=(const Material&) = delete;
  Material(Material&&) = delete;
  Material& operator=(Material&&) = delete;
  virtual ~Material() = default;

  /// Complex refractive index n + i*kappa (kappa >= 0 means absorption), absolute.
  /// Does not throw; called while compiling a system for tracing. Outside
  /// wavelength_range_um() the result is not meaningful (rtt-compile rejects such systems).
  /// @param wavelength_um vacuum wavelength in micrometre
  /// @param temperature_c medium temperature in degree Celsius
  /// @param pressure_atm  medium pressure in atm (ADR 0014, point 3)
  [[nodiscard]] virtual math::Complex index(double wavelength_um,
                                            double temperature_c,
                                            double pressure_atm) const = 0;

  /// Wavelengths on which the material data are valid; std::nullopt means unbounded.
  [[nodiscard]] virtual std::optional<WavelengthRange> wavelength_range_um() const {
    return std::nullopt;
  }
};

/// Thrown by MaterialLibrary::resolve for references that are unknown or malformed.
class UnknownMaterial : public std::runtime_error {
 public:
  using std::runtime_error::runtime_error;
};

/// Resolves material references from the model (e.g. "VACUUM", "CONST:1.5168") to materials.
///
/// Supported references (M1):
/// | Reference           | Medium                                                  |
/// | ------------------- | ------------------------------------------------------- |
/// | `VACUUM`            | n = 1                                                   |
/// | `AIR`               | dry air after Ciddor at T and p (rtt/material/air.hpp)  |
/// | `CONST:<n>`         | constant real index, n > 0, e.g. `CONST:1.5168`         |
/// | `CONST:<n>,<kappa>` | constant complex index n + i*kappa, n > 0, kappa >= 0   |
///
/// Numbers use the C locale ('.' as decimal separator, optional exponent) without
/// whitespace or a leading '+'; they must be finite. kappa = -0 is stored as +0.
/// Matching is case-sensitive. Further materials are registered with add(); glass catalogues
/// in the AGF format are loaded with add_catalog() and resolved as `KATALOG:NAME`.
///
/// resolve, add and add_catalog are thread-safe; resolve returns the same object for identical
/// reference strings for the lifetime of the library. Object identity does not mean "same medium":
/// different strings such as `CONST:1.5` and `CONST:1.50` give distinct objects with equal
/// indices. The cache holds every resolved or added reference until the library is destroyed.
class MaterialLibrary {
 public:
  MaterialLibrary() = default;

  /// Returns the material for a reference string.
  /// @param reference material reference as written in the model, e.g. "CONST:1.5168"
  /// @return the material; never null and kept alive independently of the library
  /// @throws UnknownMaterial if the reference is unknown or malformed
  [[nodiscard]] std::shared_ptr<const Material> resolve(std::string_view reference) const;

  /// Registers a material under a name; resolve(name) then returns exactly this object.
  /// @param name     reference string, e.g. "SCHOTT:N-BK7"; not empty, not VACUUM, AIR or
  ///                 starting with "CONST:", and not yet registered or resolved
  /// @param material the material; not null
  /// @throws std::invalid_argument if the name is empty, reserved or already in use, or the
  ///         material is null
  void add(std::string name, std::shared_ptr<const Material> material);

  /// Loads AGF glass catalogues (rtt/material/agf.hpp). Each file becomes the catalogue named
  /// after the file name without extension in upper case (`schott.agf` -> `SCHOTT`); its glasses
  /// resolve as `SCHOTT:N-BK7`. Glasses with an unsupported dispersion formula are listed and
  /// make resolve() throw UnknownMaterial with formula, glass, file and line (#42).
  /// Manufacturer catalogues are not shipped with Raytatouille; the user provides them.
  /// @param path an .agf file, or a directory whose *.agf files (case-insensitive, not
  ///             recursive) are loaded in sorted order
  /// @throws AgfError for malformed files (with file and line)
  /// @throws std::invalid_argument if the path has no .agf file, a catalogue name is empty,
  ///         reserved or already loaded, or a glass name is already in use; nothing is
  ///         registered in that case
  void add_catalog(const std::filesystem::path& path);

 private:
  mutable std::mutex mutex_;
  mutable std::map<std::string, std::shared_ptr<const Material>, std::less<>> cache_;
  /// Loaded catalogues: name -> file.
  std::map<std::string, std::string, std::less<>> catalogs_;
  /// Catalogue glasses that exist but cannot be evaluated: reference -> message.
  std::map<std::string, std::string, std::less<>> unsupported_;
};

}  // namespace rtt::material
