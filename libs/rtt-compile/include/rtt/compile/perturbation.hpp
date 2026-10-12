#pragma once

/// @file perturbation.hpp
/// Overlay of tolerance perturbations on compile() (ADR 0036, points 1b and 1c; #254): rigid
/// motions of nodes and surfaces and changes of n_d and nu_d of element materials. The nominal
/// model stays unchanged; dimension tolerances are Param offsets on a copy of the model and not
/// part of the overlay (ADR 0036, point 1a).
///
/// Conventions (docs/architecture.md, "Konventionen"): lengths in mm, angles in rad, right-handed
/// coordinates; rotations intrinsic X -> Y -> Z as `rotation_deg` (ADR 0028, point 3), so a
/// positive angle about x turns +y towards +z.

#include <compare>
#include <cstddef>
#include <functional>
#include <map>
#include <optional>
#include <string>

#include "rtt/math/types.hpp"
#include "rtt/model/ids.hpp"

namespace rtt::compile {

/// Rigid motion T of one node or surface X in the coordinates of its nominal frame G_X (ADR 0036,
/// point 1b): T = Tr(t) * Tr(P) * R * Tr(-P), R = Rx(a) * Ry(b) * Rz(c).
struct RigidMotion {
  /// t = (decenter_x, decenter_y, shift_z) in mm, along the axes of the nominal frame of X.
  math::Vec3 translation_mm = math::Vec3::Zero();
  /// (a, b, c) = (tilt_x, tilt_y, tilt_z) in rad, intrinsic X -> Y -> Z like `rotation_deg`.
  math::Vec3 rotation_rad = math::Vec3::Zero();
  /// Pivot P in mm in the coordinates of X. None: the pivot of the pose of X in its own
  /// coordinates (ADR 0035, point 2): the pose pivot for translate_first, Q - t for
  /// rotate_first (the point that the pose maps onto its pivot Q).
  std::optional<math::Vec3> pivot_mm;

  /// True if t and the angles are exactly zero: the motion is the identity and is skipped,
  /// whatever the pivot.
  [[nodiscard]] bool is_identity() const;
};

/// Change of the material of an element segment (ADR 0036, point 1c): Delta n_d and Delta nu_d.
/// The perturbed index keeps every relative partial dispersion:
///   D = n_F - n_C, nu_d = (n_d - 1) / D, D' = (n_d + Delta n_d - 1) / (nu_d + Delta nu_d),
///   n'(lambda) = n(lambda) + Delta n_d + (D'/D - 1) (n(lambda) - n_d)
/// on the real part; kappa is unchanged. With D = 0 (CONST) n' = n + Delta n_d. The lines d, F
/// and C are the air wavelengths kLineDAirUm, kLineFAirUm and kLineCAirUm (TIE-29, p. 1), turned
/// into vacuum wavelengths with line_vacuum_wavelength_um(); the material is evaluated there at
/// the environment temperature and pressure, like its medium. n_d, n_F and n_C are therefore
/// absolute indices: Delta n_d is n_air (about 1.0003) times the same change of the catalogue
/// index relative to air.
struct MaterialChange {
  double delta_nd = 0.0;  ///< Delta n_d, dimensionless
  double delta_vd = 0.0;  ///< Delta nu_d, dimensionless

  /// True if both changes are exactly zero: the entry is skipped.
  [[nodiscard]] bool is_zero() const noexcept { return delta_nd == 0.0 && delta_vd == 0.0; }
};

/// The material of an element segment: `segment` (0-based, between surfaces `segment` and
/// `segment + 1`, ADR 0017) is required for an element with one material per segment
/// (`segment_materials`) and must be empty for an element with one material (shorthand).
struct MaterialTarget {
  std::string element;                 ///< element name (unique, node.name_duplicate)
  std::optional<std::size_t> segment;  ///< see above

  auto operator<=>(const MaterialTarget&) const = default;
};

/// Overlay for compile() (ADR 0036, point 1; #254). Each target has at most one entry: all
/// position and tilt tolerances of a target form one motion (ADR 0035, point 2).
///
/// A perturbed surface s gets G'_s = [prod over X in ancestors(s) and s, outer to inner, of
/// G_X T_X G_X^-1] * G_s, all G the nominal frames (ADR 0036, point 1b). The motion of a node
/// moves its whole subtree rigidly, also children placed relative to a surface outside it; a
/// node outside the subtree stays at its nominal place, also if it is placed relative to a
/// surface inside it. Identity motions and zero material changes are skipped, so an empty
/// overlay gives bitwise the compiled system without overlay.
struct Perturbation {
  /// Motions of assemblies and elements by node name; the root cannot be perturbed.
  std::map<std::string, RigidMotion, std::less<>> nodes;
  /// Motions of single surfaces by SurfaceId.
  std::map<model::SurfaceId, RigidMotion> surfaces;
  /// Material changes of element segments.
  std::map<MaterialTarget, MaterialChange> materials;

  /// True if there are no entries at all.
  [[nodiscard]] bool empty() const noexcept {
    return nodes.empty() && surfaces.empty() && materials.empty();
  }
};

/// Air wavelengths of the Fraunhofer lines d, F and C in um (SCHOTT TIE-29, p. 1 and Eq.
/// (2.1-1), rounded there to 0.01 nm; these values as in test_m2_acceptance.cpp).
inline constexpr double kLineDAirUm = 0.5875618;
inline constexpr double kLineFAirUm = 0.4861327;
inline constexpr double kLineCAirUm = 0.6562725;

/// Vacuum wavelength in um of a line with wavelength `air_um` (in um) in standard air of 20 degC
/// and 1 atm: lambda_vac = lambda_air * n_air(lambda_vac), n_air after Ciddor (air.hpp), solved
/// by fixed-point iteration from lambda_air (at most 10 steps, stops when a step changes
/// nothing). This is the wavelength scale of the catalogues (CatalogMaterial in agf.hpp, #25).
[[nodiscard]] double line_vacuum_wavelength_um(double air_um);

}  // namespace rtt::compile
