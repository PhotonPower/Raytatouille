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
  std::vector<Param> coefficients;  ///< A4, A6, A8, ...: A_n in mm^(1 - n) for the term A_n r^n
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

/// Ideal lens without thickness (ADR 0031): a collinear map that images the plane at s onto the
/// plane at s' with 1/s' = 1/s + 1/f for every ray, not only paraxially. Only at the single plane
/// surface of a thin_element without shape terms or phases (validate:
/// interaction.ideal_lens_not_allowed); only `transmit` acts there (compile:
/// paths.ideal_lens_event).
struct IdealLens {
  /// f in mm: the physical distance from the lens to the rear focal point in the surrounding
  /// medium, measured along the propagation; f > 0 converges. Independent of the wavelength.
  Param focal_length;
  /// Design conjugate for the optical path (ADR 0031, point 3), mm: positive for a real object in
  /// front of the lens (as System::object), negative for a virtual object; s = -object_distance.
  /// None: object at infinity. The OPD is exact only for rays from this plane.
  std::optional<Param> object_distance;
  bool operator==(const IdealLens&) const = default;
};
/// Ideal cylinder lens (ADR 0031, point 4): the ideal lens in the direction of power
/// b = (-sin psi, cos psi, 0) only; the slope along the axis a = (cos psi, sin psi, 0) stays.
struct IdealCylinderLens {
  Param focal_length;     ///< as IdealLens::focal_length, mm
  double axis_deg = 0.0;  ///< psi, local surface frame, from the x to the y axis (as a grating)
  std::optional<Param> object_distance;  ///< as IdealLens::object_distance, mm
  bool operator==(const IdealCylinderLens&) const = default;
};

using Interaction = std::variant<Fresnel,
                                 IdealMirror,
                                 IdealAntiReflection,
                                 IdealBeamSplitter,
                                 CoatingRef,
                                 IdealPolarizer,
                                 IdealRetarder,
                                 Absorber,
                                 IdealLens,
                                 IdealCylinderLens>;

// -------------------------------------------------------------- surface -----

/// Power fraction of one diffraction order at a surface with phase layers (ADR 0025, point 5):
/// independent of polarization and wavelength in M4; plain numbers, not Param.
struct DiffractionEfficiency {
  int order = 0;            ///< diffraction order m, sign as in ADR 0025
  double efficiency = 1.0;  ///< power fraction, 0 <= efficiency <= 1
  bool operator==(const DiffractionEfficiency&) const = default;
};

struct Surface {
  SurfaceId id;
  Pose pose;  ///< placement in the element
  ShapeStack shape;
  std::optional<Aperture> aperture;
  std::vector<PhaseLayer> phases;
  Interaction interaction = Fresnel{};
  /// Efficiency per diffraction order (ADR 0025, point 5). Absent: every order has efficiency 1
  /// (geometry only). Present: orders not listed have efficiency 0. Only at a surface with a
  /// phase layer, never empty. Last member with an explicit initializer, so that aggregate
  /// initialisations stay free of -Wmissing-field-initializers.
  std::optional<std::vector<DiffractionEfficiency>>
      diffraction_efficiency{};  // NOLINT(readability-redundant-member-init)
  bool operator==(const Surface&) const = default;
};

}  // namespace rtt::model
