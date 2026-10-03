#pragma once

/// @file material.hpp
/// Optical media and the resolution of material references (docs/architecture.md, "Medien").
///
/// Conventions (docs/architecture.md, "Konventionen"):
/// - wavelengths are vacuum wavelengths in micrometre (um),
/// - temperatures are in degree Celsius,
/// - the complex refractive index is n_c = n + i*kappa with kappa >= 0 for absorption,
///   matching plane waves ~ exp(i(k.r - omega t)).

#include <functional>
#include <map>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <string_view>

#include "rtt/math/types.hpp"

namespace rtt::material {

/// An optical medium. Implementations are immutable and safe to share between threads.
class Material {
 public:
  Material() = default;
  Material(const Material&) = delete;
  Material& operator=(const Material&) = delete;
  Material(Material&&) = delete;
  Material& operator=(Material&&) = delete;
  virtual ~Material() = default;

  /// Complex refractive index n + i*kappa (kappa >= 0 means absorption).
  /// Does not throw; called while compiling a system for tracing.
  /// @param wavelength_um vacuum wavelength in micrometre
  /// @param temperature_c medium temperature in degree Celsius
  [[nodiscard]] virtual math::Complex index(double wavelength_um, double temperature_c) const = 0;
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
/// | `AIR`               | n = 1 (Ciddor air and relative glass data follow in M2) |
/// | `CONST:<n>`         | constant real index, n > 0, e.g. `CONST:1.5168`         |
/// | `CONST:<n>,<kappa>` | constant complex index n + i*kappa, n > 0, kappa >= 0   |
///
/// Numbers use the C locale ('.' as decimal separator, optional exponent) without
/// whitespace or a leading '+'; they must be finite. kappa = -0 is stored as +0.
/// Matching is case-sensitive. Catalogue references such as `SCHOTT:N-BK7` are not
/// supported yet (M2).
///
/// resolve is thread-safe and returns the same object for identical reference strings for
/// the lifetime of the library. Object identity does not mean "same medium": different
/// strings such as `CONST:1.5` and `CONST:1.50` give distinct objects with equal indices.
/// The cache holds every resolved reference until the library is destroyed.
class MaterialLibrary {
 public:
  MaterialLibrary() = default;

  /// Returns the material for a reference string.
  /// @param reference material reference as written in the model, e.g. "CONST:1.5168"
  /// @return the material; never null and kept alive independently of the library
  /// @throws UnknownMaterial if the reference is unknown or malformed
  [[nodiscard]] std::shared_ptr<const Material> resolve(std::string_view reference) const;

 private:
  mutable std::mutex mutex_;
  mutable std::map<std::string, std::shared_ptr<const Material>, std::less<>> cache_;
};

}  // namespace rtt::material
