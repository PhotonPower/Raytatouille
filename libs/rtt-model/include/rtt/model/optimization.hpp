#pragma once

/// @file optimization.hpp
/// Merit function of a system: operands and generators (ADR 0030, points 2-4; file section
/// "optimization", docs/dateiformat.md). Pure data: the solver and the evaluation are in
/// rtt-optim (#166, #167, #168). The defaults below are the only place of these values.

#include <cstdint>
#include <optional>
#include <string>
#include <variant>
#include <vector>

#include "rtt/model/ids.hpp"

namespace rtt::model {

/// What every operand has (ADR 0030, point 1): residual sqrt(weight) (value - target).
struct OperandCommon {
  /// Target in the unit of the operand's value (mm, waves, dimensionless; see the operand).
  double target = 0.0;
  double weight = 1.0;  ///< >= 0 and finite (validate: merit.weight_invalid)
  /// Name of a configuration (ADR 0029); none: configuration 0 (the first, or the nominal one
  /// without a configurations section).
  std::optional<std::string> configuration;
  bool operator==(const OperandCommon&) const = default;
};

/// Paraxial quantity of a first-order operand (rtt::paraxial::first_order).
enum class FirstOrderQuantity : std::uint8_t {
  Efl,            ///< effective focal length, mm (file: "efl")
  Bfl,            ///< back focal length, mm (file: "bfl")
  ImageFNumber,   ///< image-space F-number, dimensionless (file: "image_fnumber")
  Magnification,  ///< paraxial magnification, dimensionless (file: "magnification")
};

/// A paraxial quantity of a path at a wavelength.
struct FirstOrderOperand {
  OperandCommon common;
  FirstOrderQuantity quantity = FirstOrderQuantity::Efl;
  std::string path;  ///< name of a path
  /// Index into System::wavelengths; none: the reference wavelength.
  std::optional<std::uint16_t> wavelength;
  bool operator==(const FirstOrderOperand&) const = default;
};

/// Coordinate of a real ray (file: "ray_x", "ray_y").
enum class RayCoordinate : std::uint8_t { X, Y };

/// A coordinate of one real ray at a surface, in mm in the local coordinates of the surface.
struct RayOperand {
  OperandCommon common;
  RayCoordinate coordinate = RayCoordinate::Y;
  std::string path;   ///< name of a path
  SurfaceId surface;  ///< a surface on the path
  /// Which event at `surface` (0-based) if the path meets it more than once; none: the path meets
  /// it once (validate: merit.surface_ambiguous otherwise).
  std::optional<std::uint32_t> occurrence;
  std::uint16_t field = 0;  ///< index into System::fields.points
  /// Normalised pupil coordinates as in rtt::trace::make_rays: the unit circle is the rim of the
  /// paraxial entrance pupil, +y meridional; any finite value.
  double px = 0.0;
  double py = 0.0;
  /// Index into System::wavelengths; none: the reference wavelength.
  std::optional<std::uint16_t> wavelength;
  bool operator==(const RayOperand&) const = default;
};

/// Reference point of a spot RMS (rtt::analysis::spot statistics).
enum class SpotReference : std::uint8_t {
  Centroid,  ///< about the centroid of the arrived rays (file: "centroid")
  Chief,     ///< about the chief ray (file: "chief")
};

/// RMS spot radius of one field, mm (rtt::analysis::spot with a hexapolar pupil).
struct SpotRmsOperand {
  OperandCommon common;
  std::string path;         ///< name of a path
  std::uint16_t field = 0;  ///< index into System::fields.points
  /// Index into System::wavelengths; none: the reference wavelength, or all wavelengths with
  /// their weights if `polychromatic` (then no index; a read error otherwise).
  std::optional<std::uint16_t> wavelength;
  bool polychromatic = false;
  SpotReference reference = SpotReference::Centroid;
  int rings = 6;  ///< rings of the hexapolar pupil (as analysis::SpotOptions), >= 1
  bool operator==(const SpotRmsOperand&) const = default;
};

/// RMS wavefront error of one field in waves (rtt::analysis::opd_map).
struct OpdRmsOperand {
  OperandCommon common;
  std::string path;         ///< name of a path
  std::uint16_t field = 0;  ///< index into System::fields.points
  /// Index into System::wavelengths; none: the reference wavelength.
  std::optional<std::uint16_t> wavelength;
  int grid = 33;  ///< points per side of the pupil grid (as analysis::OpdOptions), >= 1
  bool operator==(const OpdRmsOperand&) const = default;
};

/// The value of a row of the parameter table in a configuration (ADR 0029), unit of the row.
struct ParamValueOperand {
  OperandCommon common;
  std::string parameter;  ///< name of a row of System::parameters
  bool operator==(const ParamValueOperand&) const = default;
};

/// One operand (ADR 0030, point 3).
using Operand =
    std::variant<FirstOrderOperand, RayOperand, SpotRmsOperand, OpdRmsOperand, ParamValueOperand>;

/// Spot generator (file: "rms_spot", ADR 0030, point 4): two residuals per ray, about the
/// centroid or the chief ray, with Gaussian quadrature in the pupil (#168); target 0.
struct SpotGenerator {
  std::string path;                          ///< name of a path
  std::optional<std::string> configuration;  ///< as OperandCommon::configuration
  /// Indices into System::fields.points; none: all fields (with their weights). Not empty.
  std::optional<std::vector<std::uint16_t>> fields;
  /// Indices into System::wavelengths; none: all wavelengths (with their weights). Not empty.
  std::optional<std::vector<std::uint16_t>> wavelengths;
  SpotReference reference = SpotReference::Centroid;
  int rings = 3;        ///< Gauss-Legendre rings in rho^2, >= 1
  int arms = 6;         ///< arms per ring, evenly spaced, >= 1
  double weight = 1.0;  ///< >= 0 and finite
  bool operator==(const SpotGenerator&) const = default;
};

/// Wavefront generator (file: "rms_wavefront", ADR 0030, point 4): one residual per ray in
/// waves; sampling as SpotGenerator; target 0.
struct WavefrontGenerator {
  std::string path;                                       ///< name of a path
  std::optional<std::string> configuration;               ///< as OperandCommon::configuration
  std::optional<std::vector<std::uint16_t>> fields;       ///< as SpotGenerator::fields
  std::optional<std::vector<std::uint16_t>> wavelengths;  ///< as SpotGenerator::wavelengths
  int rings = 3;                                          ///< as SpotGenerator::rings
  int arms = 6;                                           ///< as SpotGenerator::arms
  double weight = 1.0;                                    ///< >= 0 and finite
  bool operator==(const WavefrontGenerator&) const = default;
};

using Generator = std::variant<SpotGenerator, WavefrontGenerator>;

/// The merit function: operands and generators in file order (ADR 0030, point 2). Empty: no
/// merit function.
struct Optimization {
  std::vector<Operand> operands;
  std::vector<Generator> generators;
  [[nodiscard]] bool empty() const noexcept { return operands.empty() && generators.empty(); }
  bool operator==(const Optimization&) const = default;
};

}  // namespace rtt::model
