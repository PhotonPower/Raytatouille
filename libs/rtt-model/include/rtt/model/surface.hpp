#pragma once

/// @file surface.hpp
/// Surface with its layer stack: shape, aperture, phases, interaction
/// (docs/architecture.md, "Layer-Stack einer Fläche").

#include <array>
#include <optional>
#include <string>
#include <variant>
#include <vector>

#include "rtt/model/ids.hpp"
#include "rtt/model/param.hpp"
#include "rtt/model/pose.hpp"

namespace rtt::model {

// ---------------------------------------------------------------- shape -----

/// Flat surface z = 0.
struct Plane {
  bool operator==(const Plane&) const = default;
};

/// Conic of revolution. Sphere for conic = 0. Radius in mm, positive if the centre of
/// curvature lies on the +z side of the vertex; must be finite and non-zero (use Plane).
struct Conic {
  Param radius;
  Param conic;
  bool operator==(const Conic&) const = default;
};

/// Conic plus even polynomial: z = conic(r) + A4 r^4 + A6 r^6 + ... (coefficients[0] = A4).
struct EvenAsphere {
  Param radius;
  Param conic;
  std::vector<Param> coefficients;  ///< A4, A6, A8, ... in mm^(1-2k)
  bool operator==(const EvenAsphere&) const = default;
};

using BaseShape = std::variant<Plane, Conic, EvenAsphere>;

/// Additive Zernike sag term, Noll ordering (coefficients[0] = j = 1), coefficients in mm.
struct ZernikeSag {
  Param normalization_radius;  ///< mm
  std::vector<Param> coefficients;
  bool operator==(const ZernikeSag&) const = default;
};

using ShapeTerm = std::variant<ZernikeSag>;

/// Total sag = base + sum of terms (terms evaluated in the local xy plane).
struct ShapeStack {
  BaseShape base = Plane{};
  std::vector<ShapeTerm> terms;
  bool operator==(const ShapeStack&) const = default;
};

// ------------------------------------------------------------- aperture -----

struct CircularAperture {
  double radius = 0.0;        ///< mm
  double inner_radius = 0.0;  ///< mm, > 0 makes an annulus
  bool operator==(const CircularAperture&) const = default;
};

struct RectangularAperture {
  double half_width_x = 0.0;  ///< mm
  double half_width_y = 0.0;  ///< mm
  bool operator==(const RectangularAperture&) const = default;
};

struct EllipticalAperture {
  double semi_axis_x = 0.0;  ///< mm
  double semi_axis_y = 0.0;  ///< mm
  bool operator==(const EllipticalAperture&) const = default;
};

using Aperture = std::variant<CircularAperture, RectangularAperture, EllipticalAperture>;

// --------------------------------------------------------------- phases -----

/// Straight-line grating. orientation_deg = 0 means grooves parallel to the local y axis.
struct LinearGrating {
  Param lines_per_mm;
  double orientation_deg = 0.0;
  bool operator==(const LinearGrating&) const = default;
};

/// Rotationally symmetric phase phi = sum c_k rho^(2k), rho = r / normalization_radius,
/// coefficients in radian, coefficients[0] belongs to rho^2.
struct RadialPhase {
  Param normalization_radius;  ///< mm
  std::vector<Param> coefficients;
  bool operator==(const RadialPhase&) const = default;
};

using PhaseLayer = std::variant<LinearGrating, RadialPhase>;

// ---------------------------------------------------------- interaction -----

/// Uncoated interface; Fresnel coefficients from the two media (default).
struct Fresnel {
  bool operator==(const Fresnel&) const = default;
};
/// Lossless mirror, R = 1 for s and p, no phase difference.
struct IdealMirror {
  bool operator==(const IdealMirror&) const = default;
};
/// Lossless anti-reflection interface, T = 1.
struct IdealAntiReflection {
  bool operator==(const IdealAntiReflection&) const = default;
};
/// Fully absorbing surface.
struct Absorber {
  bool operator==(const Absorber&) const = default;
};
/// Lossless beam splitter, T = 1 - R per polarization.
struct IdealBeamSplitter {
  double reflectance_s = 0.5;
  double reflectance_p = 0.5;
  bool operator==(const IdealBeamSplitter&) const = default;
};
/// Thin-film coating defined in a coating catalog.
struct CoatingRef {
  std::string name;
  bool operator==(const CoatingRef&) const = default;
};
/// Ideal linear polarizer; axis in element coordinates, projected perpendicular to k.
struct IdealPolarizer {
  std::array<double, 3> transmission_axis{1.0, 0.0, 0.0};
  double extinction_ratio = 0.0;  ///< T_min / T_max, 0 = perfect
  bool operator==(const IdealPolarizer&) const = default;
};
/// Ideal linear retarder; positive retardance delays the slow axis.
struct IdealRetarder {
  std::array<double, 3> fast_axis{1.0, 0.0, 0.0};
  double retardance_waves = 0.25;
  bool operator==(const IdealRetarder&) const = default;
};

using Interaction = std::variant<Fresnel,
                                 IdealMirror,
                                 IdealAntiReflection,
                                 IdealBeamSplitter,
                                 CoatingRef,
                                 IdealPolarizer,
                                 IdealRetarder,
                                 Absorber>;

// -------------------------------------------------------------- surface -----

struct Surface {
  SurfaceId id;
  Pose pose;  ///< placement in the element
  ShapeStack shape;
  std::optional<Aperture> aperture;
  std::vector<PhaseLayer> phases;
  Interaction interaction = Fresnel{};
  bool operator==(const Surface&) const = default;
};

}  // namespace rtt::model
