#pragma once

/// @file system.hpp
/// Root of the optical model. Pure data, no tracing logic.

#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

#include "rtt/model/element.hpp"
#include "rtt/model/param.hpp"
#include "rtt/model/path.hpp"

namespace rtt::model {

/// Version of the file format written by this library. rtt-io also reads 0.1 files and
/// migrates them (0.2: per-segment materials of Lens and Plate, ADR 0017).
inline constexpr std::string_view kSchemaVersion = "0.2.0";

struct Wavelength {
  double um = 0.0;  ///< vacuum wavelength in micrometre
  double weight = 1.0;
  bool reference = false;  ///< exactly one wavelength is the reference
  bool operator==(const Wavelength&) const = default;
};

enum class FieldType : std::uint8_t { AngleDeg, ObjectHeight, ParaxialImageHeight };

struct Field {
  double x = 0.0;  ///< degree or mm depending on FieldType
  double y = 0.0;
  double weight = 1.0;
  bool operator==(const Field&) const = default;
};

struct FieldSet {
  FieldType type = FieldType::AngleDeg;
  std::vector<Field> points;
  bool operator==(const FieldSet&) const = default;
};

/// What SystemAperture::value means.
enum class SystemApertureType : std::uint8_t {
  EntrancePupilDiameter,  ///< mm
  /// Dimensionless F-number at infinite conjugates: EPD = |EFL| / F#, also for a finite object.
  ImageSpaceFNumber,
  /// Object-space numerical aperture NA, dimensionless, finite objects only, read paraxially:
  /// the paraxial marginal slope from the axial object point is u = NA / n, n the refractive
  /// index of object space (first_order, ray aiming of object-space telecentric systems). With
  /// paraxial aiming the marginal ray has tan U = NA / n; with real aiming it hits the paraxial
  /// stop edge R_s, and its tan U differs from NA / n by the pupil aberration.
  /// Greivenkamp, OPTI-502, Sec. 9, p. 9-34:
  /// NA = n sin U ≈ n u; a real marginal ray with n sin U = NA differs for a large NA
  /// (docs/quellen.md). With an object-space telecentric system only ObjectSpaceNA and
  /// StopSize define the bundle.
  ObjectSpaceNA,
  StopSize,  ///< use the stop aperture as defined; the value is not used
};

struct SystemAperture {
  SystemApertureType type = SystemApertureType::EntrancePupilDiameter;
  Param value;
  bool operator==(const SystemAperture&) const = default;
};

struct ObjectSpace {
  bool at_infinity = true;
  Param distance;  ///< mm from the object to the global origin along -z; used if !at_infinity
  bool operator==(const ObjectSpace&) const = default;
};

struct Environment {
  double temperature_c = 20.0;
  double pressure_atm = 1.0;
  std::string medium = "AIR";
  bool operator==(const Environment&) const = default;
};

struct System {
  std::string schema_version{kSchemaVersion};
  std::string name;
  Environment environment;
  ObjectSpace object;
  std::vector<Wavelength> wavelengths;
  SystemAperture aperture;
  FieldSet fields;
  Assembly root;
  std::vector<Path> paths;
  bool operator==(const System&) const = default;
};

}  // namespace rtt::model
