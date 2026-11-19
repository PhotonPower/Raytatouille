#pragma once

/// @file surface.hpp
/// Surface with its layer stack: shape, aperture, phases, interaction
/// (docs/architecture.md, "Layer-Stack einer Fläche").

#include <array>
#include <cstdint>
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

/// Biconic (ADR 0042, point 1): z = (c_x x^2 + c_y y^2) / (1 + sqrt(1 - (1 + k_x) c_x^2 x^2 -
/// (1 + k_y) c_y^2 y^2)) with c_x = 1/radius_x, c_y = 1/radius_y, in the local frame of the
/// surface (its rotation about z orients the sections). Radii in mm with the sign of Conic::radius
/// in their section (x-z, y-z); both finite and non-zero. No asphere terms.
struct Biconic {
  Param radius_x;  ///< R_x, mm, section x-z
  Param radius_y;  ///< R_y, mm, section y-z
  Param conic_x;   ///< k_x, dimensionless
  Param conic_y;   ///< k_y, dimensionless
  bool operator==(const Biconic&) const = default;
};

/// Toroid (ADR 0042, point 1): the profile z_y(y) = c y^2 / (1 + sqrt(1 - (1 + k) c^2 y^2)) +
/// A4 y^4 + A6 y^6 + ... in the local y-z plane, c = 1/radius, rotated about the axis parallel to
/// y through (x, z) = (0, radius_of_rotation). Radii in mm, sign as Conic::radius (positive: the
/// centre of curvature, or the axis of rotation, on the +z side of the vertex); both finite and
/// non-zero. radius_of_rotation == radius with conic 0 and no coefficients is a sphere.
struct Toroid {
  Param radius;                     ///< R of the y-z profile, mm
  Param conic;                      ///< k of the profile, dimensionless
  std::vector<Param> coefficients;  ///< A4, A6, ... of the profile: A_n in mm^(1 - n)
  Param radius_of_rotation;         ///< R, mm: the axis of rotation lies at z = R
  bool operator==(const Toroid&) const = default;
};

/// Cylinder (ADR 0042, point 1): the profile z = z_y(y) of Toroid, flat in x. Radius in mm,
/// sign as Conic::radius; finite and non-zero.
struct Cylinder {
  Param radius;                     ///< R of the y-z profile, mm
  Param conic;                      ///< k of the profile, dimensionless
  std::vector<Param> coefficients;  ///< A4, A6, ... of the profile: A_n in mm^(1 - n)
  bool operator==(const Cylinder&) const = default;
};

/// Odd asphere (ADR 0043, point 1; OpticStudio "Odd Asphere"): z = c r^2 / (1 + sqrt(1 - (1 + k)
/// c^2 r^2)) + A1 r + A2 r^2 + A3 r^3 + ... with r^2 = x^2 + y^2. The base carries its curvature
/// c (finite, 0 allowed: a plane base) instead of a radius (decision F10), so that c stays a
/// variable also through 0; a conic constant != 0 has no effect at c = 0. The vertex curvature is
/// c + 2 A2, not c. A1 != 0 makes a cone at the vertex (no paraxial model there, ADR 0043,
/// point 3).
struct OddAsphere {
  Param curvature;                  ///< c of the base conic, 1/mm, positive as 1/Conic::radius
  Param conic;                      ///< k of the base conic, dimensionless
  std::vector<Param> coefficients;  ///< A1, A2, A3, ...: A_n in mm^(1 - n) for the term A_n r^n
  bool operator==(const OddAsphere&) const = default;
};

/// Axicon (ADR 0043, point 1, derived there): the cone z = r tan(alpha) about the local z axis,
/// apex at the vertex. alpha is the base angle between the mantle and the plane z = 0, with the
/// sign of a curvature: alpha > 0 lets the mantle rise towards +z; alpha = 0 is the plane;
/// |alpha| < 90 deg (decision F5). No paraxial model at the apex (F6).
struct Axicon {
  Param angle_deg;  ///< base angle alpha, degree
  bool operator==(const Axicon&) const = default;
};

/// Forbes Q-con asphere (ADR 0043, point 1; Forbes, Opt. Express 15(8), 5218 (2007), Eqs. (2),
/// (3), (8)): z = c r^2 / (1 + sqrt(1 - (1 + k) c^2 r^2)) + u^4 sum_m a_m Q^con_m(u^2) with
/// u = r / normalization_radius. The base carries its curvature c (finite, 0 allowed) as the odd
/// asphere; Forbes writes the conic as epsilon = 1 + k. Beyond the normalization radius the
/// polynomials are extrapolated (a warning if the aperture reaches beyond it).
struct QCon {
  Param curvature;                  ///< c of the base conic, 1/mm, positive as 1/Conic::radius
  Param conic;                      ///< k of the base conic, dimensionless
  Param normalization_radius;       ///< rho, mm, finite and > 0 (half the clear aperture)
  std::vector<Param> coefficients;  ///< a_0, a_1, ...: mm
  bool operator==(const QCon&) const = default;
};

/// Forbes Q-bfs asphere with a conic base (ADR 0043, point 1; Forbes, Opt. Express 19(10), 9923
/// (2011), Eq. (2.4)): z = c r^2 / (1 + phi) + u^2 (1 - u^2) / sigma sum_m a_m Q^bfs_m(u^2) with
/// phi = sqrt(1 - (1 + k) c^2 r^2), 1/sigma = sqrt(1 - k c^2 r^2) / phi and u = r /
/// normalization_radius. Defined for phi > 0; (1 + k) c^2 rho^2 < 1 is required. With k = 0, c is
/// the curvature of the best-fit sphere through the vertex and the rim u = 1. The vertex
/// curvature is c + 2 S(0) / rho^2, not c.
struct QBfs {
  Param curvature;                  ///< c of the base conic, 1/mm, positive as 1/Conic::radius
  Param conic;                      ///< k of the base conic, dimensionless
  Param normalization_radius;       ///< rho, mm, finite and > 0 (half the clear aperture)
  std::vector<Param> coefficients;  ///< a_0, a_1, ...: mm
  bool operator==(const QBfs&) const = default;
};

using BaseShape = std::
    variant<Plane, Conic, EvenAsphere, Biconic, Toroid, Cylinder, OddAsphere, Axicon, QCon, QBfs>;

/// Additive Zernike sag term, Noll ordering (coefficients[0] = j = 1), orthonormal polynomials
/// (Niu/Tian Eq. (3), as rtt::math::zernike()), coefficients in mm; at most 231 (ADR 0051,
/// point 2). rho = r / normalization_radius, theta from +x counter-clockwise.
struct ZernikeSag {
  Param normalization_radius;  ///< mm
  std::vector<Param> coefficients;
  bool operator==(const ZernikeSag&) const = default;
};

/// Additive Zernike sag term in the Fringe set: coefficients[k] belongs to Fringe term k + 1
/// (Niu/Tian Table 4, j = k), polynomials without normalisation (Eq. (22), 1 at the rim),
/// coefficients in mm; at most 37 (ADR 0051, point 3).
struct ZernikeFringeSag {
  Param normalization_radius;  ///< mm
  std::vector<Param> coefficients;
  bool operator==(const ZernikeFringeSag&) const = default;
};

/// Additive XY polynomial term: z = sum_k c_k (x/R)^i_k (y/R)^j_k without a constant term,
/// coefficients[0] = x, then y, x^2, x y, y^2, x^3, ... (degree ascending, within a degree the
/// power of x descending; [UG] "Extended Polynomial"), coefficients in mm, at most 230 (degree 20;
/// ADR 0051, point 4). R = normalization_radius in mm.
struct XyPolynomial {
  Param normalization_radius;  ///< mm
  std::vector<Param> coefficients;
  bool operator==(const XyPolynomial&) const = default;
};

using ShapeTerm = std::variant<ZernikeSag, ZernikeFringeSag, XyPolynomial>;

/// Total sag = base + sum of terms (terms evaluated in the local xy plane).
struct ShapeStack {
  BaseShape base = Plane{};
  std::vector<ShapeTerm> terms;
  bool operator==(const ShapeStack&) const = default;
};

// ------------------------------------------------------------- aperture -----

// Aperture dimensions are Params (ADR 0050, point 1): variable for the optimizer and bindable to
// the parameter table, but never toleranced. Local surface coordinates, mm.

struct CircularAperture {
  Param radius;        ///< mm, > 0
  Param inner_radius;  ///< mm, > 0 makes an annulus (central obscuration, ADR 0032)
  bool operator==(const CircularAperture&) const = default;
};

struct RectangularAperture {
  Param half_width_x;  ///< mm, > 0
  Param half_width_y;  ///< mm, > 0
  bool operator==(const RectangularAperture&) const = default;
};

struct EllipticalAperture {
  Param semi_axis_x;  ///< mm, > 0
  Param semi_axis_y;  ///< mm, > 0
  bool operator==(const EllipticalAperture&) const = default;
};

/// Regular polygon (ADR 0050, point 2): vertex k at circumradius R and angle
/// phi_k = rotation_deg + 360 k / sides degree, k = 0 ... sides - 1, measured from the local +x
/// axis towards +y; rotation 0 puts the first vertex on +x.
struct RegularPolygonAperture {
  std::uint16_t sides = 0;    ///< 3 ... 1024
  Param circumradius;         ///< R in mm, > 0
  double rotation_deg = 0.0;  ///< theta_0 in degree, finite
  bool operator==(const RegularPolygonAperture&) const = default;
};

/// Simple polygon (ADR 0050, point 2): at least 3 vertices (x, y) in mm in the local frame, in
/// either orientation, convex or concave, without self-intersection; the last vertex joins the
/// first and is not repeated. The vertices are plain numbers, not Params.
struct PolygonAperture {
  std::vector<std::array<double, 2>> vertices;
  bool operator==(const PolygonAperture&) const = default;
};

using Aperture = std::variant<CircularAperture,
                              RectangularAperture,
                              EllipticalAperture,
                              RegularPolygonAperture,
                              PolygonAperture>;

/// The vertices of a polygon aperture in local x, y in mm (ADR 0050, point 2): those of
/// PolygonAperture as given, those of RegularPolygonAperture at R (cos phi_k, sin phi_k) with
/// phi_k = (rotation_deg + 360 k / sides) degree; empty for the other apertures and for fewer
/// than 3 sides. Uses the current value of the circumradius (resolve bound Params first).
[[nodiscard]] std::vector<std::array<double, 2>> aperture_vertices(const Aperture& aperture);

// ---------------------------------------------------------- obscurations -----

// Obscurations (ADR 0050, point 3) block light inside a shape at a surface; the surface is still
// present there. Local surface coordinates, mm. The centre lies at (decenter_x, decenter_y); the
// shape's own axes u, v follow from rotating x, y by rotation_deg (positive turns +x towards +y):
//   u = cos(theta) (x - x0) + sin(theta) (y - y0),  v = -sin(theta) (x - x0) + cos(theta) (y - y0).
// A ray is blocked strictly inside the shape shrunk by kApertureTolerance; the edge passes.

/// Disc of radius `radius` around (decenter_x, decenter_y).
struct CircularObscuration {
  Param radius;      ///< mm, > 0
  Param decenter_x;  ///< x of the centre, mm
  Param decenter_y;  ///< y of the centre, mm
  bool operator==(const CircularObscuration&) const = default;
};

/// Rectangle |u| <= half_width_x, |v| <= half_width_y (a spider arm, for instance).
struct RectangularObscuration {
  Param half_width_x;         ///< mm, > 0, along u
  Param half_width_y;         ///< mm, > 0, along v
  Param decenter_x;           ///< x of the centre, mm
  Param decenter_y;           ///< y of the centre, mm
  double rotation_deg = 0.0;  ///< theta in degree, finite
  bool operator==(const RectangularObscuration&) const = default;
};

/// Ellipse (u / semi_axis_x)^2 + (v / semi_axis_y)^2 <= 1.
struct EllipticalObscuration {
  Param semi_axis_x;          ///< mm, > 0, along u
  Param semi_axis_y;          ///< mm, > 0, along v
  Param decenter_x;           ///< x of the centre, mm
  Param decenter_y;           ///< y of the centre, mm
  double rotation_deg = 0.0;  ///< theta in degree, finite
  bool operator==(const EllipticalObscuration&) const = default;
};

using Obscuration =
    std::variant<CircularObscuration, RectangularObscuration, EllipticalObscuration>;

/// True for a CircularObscuration with decenter (0, 0) exactly: a centred circle, treated like
/// CircularAperture::inner_radius by the virtual rays and the ring cells (ADR 0050, point 3).
[[nodiscard]] bool is_centred_circle(const Obscuration& obscuration) noexcept;

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
  /// phase layer, never empty.
  std::optional<std::vector<DiffractionEfficiency>>
      diffraction_efficiency{};  // NOLINT(readability-redundant-member-init)
  /// Obscurations at this surface (ADR 0050, point 3), checked after the aperture; also allowed
  /// without an aperture. The last members carry explicit initializers, so that aggregate
  /// initialisations stay free of -Wmissing-field-initializers.
  std::vector<Obscuration> obscurations{};  // NOLINT(readability-redundant-member-init)
  bool operator==(const Surface&) const = default;
};

}  // namespace rtt::model
