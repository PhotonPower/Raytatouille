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

/// Version of the file format written by this library.
inline constexpr std::string_view kSchemaVersion = "0.1.0";

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

enum class SystemApertureType : std::uint8_t {
  EntrancePupilDiameter,  ///< mm
  ImageSpaceFNumber,      ///< dimensionless
  ObjectSpaceNA,          ///< dimensionless
  StopSize,               ///< use the stop aperture as defined
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
