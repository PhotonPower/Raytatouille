#pragma once

/// @file compiled_system.hpp
/// Immutable, flattened form of a model::System on which the tracers work
/// (docs/architecture.md, "Kompilierung").
///
/// Conventions (docs/architecture.md, "Konventionen"): lengths in mm, wavelengths in um
/// (vacuum), temperature in degC, right-handed global coordinates with the optical axis +z.
/// Complex indices are n + i*kappa with kappa >= 0 for absorption.

#include <compare>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <variant>
#include <vector>

#include "rtt/coating/transfer_matrix.hpp"
#include "rtt/compile/compiled_element.hpp"
#include "rtt/compile/perturbation.hpp"
#include "rtt/geom/asphere.hpp"
#include "rtt/geom/conic.hpp"
#include "rtt/geom/phase.hpp"
#include "rtt/geom/plane.hpp"
#include "rtt/material/material.hpp"
#include "rtt/math/isometry.hpp"
#include "rtt/math/types.hpp"
#include "rtt/model/model.hpp"

namespace rtt::coating {
class CoatingLibrary;  // rtt/coating/catalog.hpp
}  // namespace rtt::coating

namespace rtt::compile {

/// Resolved base shape of a surface in its local coordinates. Supported in M1: plane, conic
/// and even asphere; further shapes are added to this variant as rtt-geom provides them.
using CompiledShape =
    std::variant<geom::Plane<double>, geom::Conic<double>, geom::EvenAsphere<double>>;

/// Resolved data of an ideal lens or ideal cylinder lens (ADR 0031): Params from the parameter
/// table already resolved for the compiled configuration.
struct CompiledIdealLens {
  double focal_length = 0.0;  ///< f in mm along the propagation, f > 0 converges; never 0
  /// t_o = 1/s of the design conjugate in 1/mm, with s = -object_distance (ADR 0031, point 3);
  /// 0 for an object at infinity (no object_distance).
  double t_o = 0.0;
  /// Cylinder lens only: the direction of power b = (-sin psi, cos psi, 0) in local surface
  /// coordinates, psi = axis_deg in rad; none for the ideal lens (power in every direction).
  std::optional<math::Vec3> power_axis;
};

/// Coating of a surface (ADR 0019): which compiled coating, and which side is the substrate.
struct SurfaceCoating {
  std::uint32_t coating = 0;  ///< index into CompiledSystem::coatings()
  /// Index into CompiledSystem::media() of the substrate: the inside of the element the surface
  /// belongs to. Light that arrives from this medium sees the layers in reverse order.
  std::uint32_t substrate_medium = 0;
};

/// A coating design evaluated at all system wavelengths (ADR 0019).
struct CompiledCoating {
  std::string reference;  ///< coating reference as in the model, e.g. "DEMO:AR_MGF2"
  /// Layers per system wavelength (same order as wavelengths_um()), each from the ambient side
  /// to the substrate: complex index at the environment temperature and pressure, physical
  /// thickness in um (QWOT converted at the design wavelength, same conditions).
  std::vector<std::vector<coating::Layer<double>>> layers;
};

/// One surface with everything the tracer needs, independent of the model.
struct CompiledSurface {
  model::SurfaceId id;  ///< stable surface id from the model, e.g. "L1.S1"
  model::ElementKind element_kind = model::ElementKind::Lens;  ///< kind of the owning element
  std::string element_name;                                    ///< name of the owning element
  std::uint32_t element = 0;  ///< index of the owning element in CompiledSystem::elements()
  /// Media before and after the surface in the surface order of its element (#81), indices into
  /// CompiledSystem::media(): what a refraction through the element's surfaces in order, coming
  /// from the environment, would see (rules of ADR 0017: segments for lenses and segmented
  /// plates, inside/environment toggle otherwise; environment on both sides for an element
  /// without material). Independent of the paths: the media of a path are in its events, and a
  /// reflection or a path running backwards may see the surface differently. For plates and
  /// prisms of one material with more than two surfaces this is only the convention of the
  /// surface order (each refraction toggles inside/environment); which side is glass follows
  /// from CompiledElement::media and CompiledElement::segmented.
  std::uint32_t medium_front = 0;
  std::uint32_t medium_back = 0;
  /// JSON pointer of the surface in the system file, e.g. "/root/children/1/surfaces/0"
  /// (ADR 0022: errors and the GUI point at the surface with it).
  std::string location;
  /// Local surface coordinates -> global coordinates (assembly -> element -> surface poses
  /// chained, see rtt::math::Isometry3::from_pose for the convention).
  math::Isometry3 to_global;
  /// Global coordinates -> local surface coordinates (inverse of to_global).
  math::Isometry3 to_local;
  CompiledShape shape = geom::Plane<double>{};
  std::optional<model::Aperture> aperture;  ///< in local coordinates, mm; none = unbounded
  std::vector<model::PhaseLayer> phases;    ///< as in the model (export, display)
  /// The phase layers resolved for the tracer (ADR 0025, point 1, #127), in the order of
  /// `phases`: parameter values, orientation_deg in rad, coefficients in rad unchanged; phase in
  /// rad over the local (x, y) in mm.
  std::vector<geom::PhaseFunction<double>> phase_functions;
  model::Interaction interaction = model::Fresnel{};
  /// Efficiency per diffraction order as in the model (ADR 0025, point 5): absent = every order
  /// has efficiency 1; present = orders not listed have efficiency 0.
  std::optional<std::vector<model::DiffractionEfficiency>> diffraction_efficiency;
  /// Set for a CoatingRef interaction: the resolved coating and its substrate side.
  std::optional<SurfaceCoating> coating;
  /// Set for IdealPolarizer (transmission axis) and IdealRetarder (fast axis): the axis of the
  /// model, given in element coordinates, rotated into global coordinates (assembly and element
  /// poses, not the surface pose); not normalised. The tracer projects it perpendicular to the
  /// ray (rtt/polar/ideal.hpp; ADR 0021).
  std::optional<math::Vec3> ideal_axis;
  /// True for every surface of a Mirror element without material (ADR 0034, point 7): a Fresnel
  /// reflection there is that of an ideal conductor, r_s = -1, r_p = +1 (docs/architecture.md,
  /// "Polarisation"). A mirror with material reflects with Fresnel against its material, also if
  /// its index equals that of the medium in front.
  bool ideal_conductor = false;
  /// Set for IdealLens and IdealCylinderLens: the resolved lens data (ADR 0031).
  std::optional<CompiledIdealLens> ideal_lens;
};

/// A homogeneous medium evaluated at all system wavelengths: isotropic, or a uniaxial crystal
/// (ADR 0026, point 3). Isotropic media with equal references are shared; every crystal element
/// has a medium of its own, because the optic axis belongs to it, and so has every element
/// segment with a material change of a Perturbation (ADR 0036, point 1c), under the reference of
/// its nominal material.
struct CompiledMedium {
  /// Material reference as in the model, e.g. "AIR", "CONST:1.5168"; for a crystal the
  /// reference of the ordinary index n_O.
  std::string reference;
  /// Complex index n + i*kappa per system wavelength (same order as wavelengths_um()),
  /// evaluated at the environment temperature and pressure; for a crystal n_O (kappa = 0, since
  /// compile rejects absorbing crystals with crystal.absorbing).
  std::vector<math::Complex> index;
  /// Crystal only (empty otherwise): reference of the extraordinary index n_E.
  std::string reference_extraordinary{};  // NOLINT(readability-redundant-member-init)
  /// Crystal only (empty otherwise): n_E per system wavelength, real, same order as `index`.
  std::vector<double> index_extraordinary{};  // NOLINT(readability-redundant-member-init)
  /// Crystal only: optic axis in global coordinates, unit vector; the model's axis in element
  /// coordinates rotated with the assembly and element poses (ADR 0026, point 2). Its sign has
  /// no meaning.
  std::optional<math::Vec3> optic_axis{};  // NOLINT(readability-redundant-member-init)

  /// True for a uniaxial crystal (optic_axis set).
  [[nodiscard]] bool is_crystal() const noexcept { return optic_axis.has_value(); }
};

/// Mode of the ray inside a crystal (ADR 0026, point 3): None outside crystals.
enum class CrystalMode : std::uint8_t {
  None,           ///< the medium before the event is not a crystal
  Ordinary,       ///< the ray entered the crystal with an Ordinary event
  Extraordinary,  ///< the ray entered the crystal with an Extraordinary event
};

/// One step of a compiled path.
struct CompiledEvent {
  std::uint32_t surface = 0;  ///< index into CompiledSystem::surfaces()
  model::EventKind kind = model::EventKind::Refract;
  int order = 0;  ///< diffraction order (ADR 0025); 0 at a surface without phase layer
  std::uint32_t medium_before = 0;  ///< index into CompiledSystem::media()
  std::uint32_t medium_after = 0;   ///< index into CompiledSystem::media()
  /// Index into CompiledSystem::media() of the medium on the other side of the surface: the
  /// medium a Refract at this surface would enter from medium_before (equal to medium_after for
  /// Refract). Fresnel and coatings need it for Reflect (ADR 0021). For an element without
  /// material, and for an inner surface whose side is ambiguous, it equals medium_before.
  std::uint32_t medium_beyond = 0;
  /// True if the ray is inside the element of this surface before the event (in one of its
  /// segments). For a coated surface this is the substrate side (ADR 0019): the tracer then uses
  /// the reversed layer stack (ADR 0021).
  bool from_inside = false;
  /// For an event whose medium before is a crystal: the mode with which the path entered this
  /// crystal (ADR 0026, point 3), so that the exit knows it; None otherwise.
  CrystalMode crystal_mode = CrystalMode::None;
};

/// Named, ordered list of events.
struct CompiledPath {
  std::string name;  ///< path name from the model, unique within the system
  std::vector<CompiledEvent> events;
};

/// Index of a path in CompiledSystem::paths().
struct PathId {
  std::uint32_t index = 0;
  auto operator<=>(const PathId&) const = default;
};

/// Thrown by compile() for invalid models, unknown materials and model features that are not
/// supported yet. diagnostics() holds the errors with JSON pointers into the system file.
class CompileError : public std::runtime_error {
 public:
  explicit CompileError(std::vector<model::Diagnostic> diagnostics);

  [[nodiscard]] const std::vector<model::Diagnostic>& diagnostics() const noexcept {
    return diagnostics_;
  }

 private:
  std::vector<model::Diagnostic> diagnostics_;
};

class CompiledSystem;

/// Compiles a model into an immutable CompiledSystem.
///
/// - Runs model::validate; any error diagnostic throws CompileError. Its warnings are kept in
///   CompiledSystem::diagnostics() (ADR 0022).
/// - Flattens the tree: the global transform of each surface is root pose * assembly poses *
///   element pose * surface pose. Pickups are not evaluated before M5; their value is used.
/// - Resolves every material (environment medium, element materials in both forms: shorthand
///   and one material per segment, ADR 0017) with `materials` and evaluates it once per system
///   wavelength at environment.temperature_c (ADR 0014). Equal references share one medium.
/// - Builds every path. An automatic path visits all surfaces in tree order with Refract for
///   Lens and Plate, Reflect for Mirror and Transmit for Stop, Detector and ThinElement.
/// - Media along a path (docs/architecture.md, "Medien entlang eines Pfads"): the ray starts in
///   the environment medium. Reflect and Transmit keep the medium, with any order (ADR 0025),
///   and so does every event at an element without material. The rule for Refract, Ordinary
///   and Extraordinary depends only on the element, never on the kind of path (decided for #27):
///   | element                                  | automatic path         | explicit path |
///   | Lens, one material or several            | segment rule           | segment rule  |
///   | Plate, one material (shorthand or equal) | 2 surfaces: toggle,    | toggle        |
///   |                                          | > 2: CompileError      |               |
///   | Plate, different segment materials       | segment rule           | segment rule  |
///   | Mirror with substrate                    | 1 surface: toggle,     | toggle        |
///   |                                          | > 1: CompileError (#6) |               |
///   Toggle: between the inside and the environment at every surface (prisms, cubes; the rule
///   of #5). For 2 surfaces it equals the segment rule. Segment rule at surface i of an element
///   with N surfaces (segment j between surfaces j and j + 1, zero-based): from segment i - 1
///   into segment i, or into the environment if i = N - 1; from segment i into segment i - 1,
///   or into the environment if i = 0; from another segment of the same element into the
///   environment; from outside the element through the first surface into segment 0, through
///   the last surface into segment N - 2, and through an inner surface i into the material of
///   segments i - 1 and i if both are the same, otherwise CompileError (ambiguous side).
///   Elements do not nest: entering element B while inside element A leaves A (A's medium
///   before, B's medium after), and leaving B then
///   goes to the environment. A cemented group is one element with one material per segment;
///   two separate elements always meet through the environment.
///
/// Not supported yet (CompileError): Zernike sag terms (M8).
/// Crystals (ADR 0026): every crystal element gets a medium of its own (CompiledMedium with
/// n_O, n_E and the global optic axis; parts resolved at .../material/ordinary and
/// .../extraordinary, crystal.absorbing for kappa != 0), and every event inside a crystal its
/// CompiledEvent::crystal_mode. Path rules of ADR 0026, point 4 (docs/architecture.md): a
/// crystal is entered with Ordinary or Extraordinary (paths.crystal_mode_required), modes only
/// enter crystals (paths.mode_without_crystal), crystal to crystal, reflection inside or at a
/// crystal from outside (except ideal_anti_reflection), an order inside a crystal and the
/// automatic path through a crystal are crystal.unsupported; interactions other than fresnel and
/// ideal_anti_reflection at a crystal surface are crystal.interaction_unsupported.
/// Also a CompileError: a Mirror with substrate material and more than one surface on an
/// automatic path (Mangin mirror; its front surface refracts, so it needs an explicit path
/// Refract, Reflect, Refract); a Plate of one material with more than 2 surfaces on an
/// automatic path (prism or cube, needs an explicit path); and, under the segment rule, a
/// Refract, Ordinary or Extraordinary event at an inner surface reached from outside the
/// element when the segments on its two sides have different materials (ambiguous side).
///
/// Coatings (ADR 0019): every CoatingRef is resolved with `coatings`, its layer materials with
/// `materials`, and the layers are evaluated at every system wavelength at the environment
/// temperature and pressure (QWOT thickness at the design wavelength under the same
/// conditions). The substrate is the inside of the element the surface belongs to: the
/// adjacent segment of a Lens or Plate (the plate's material for a Plate of one material), the
/// substrate of a Mirror. CompileError at .../interaction or .../interaction/name for an
/// unknown coating, an unknown layer material, a QWOT whose design wavelength lies outside the
/// material's range or has Re n <= 0, a system wavelength outside the range of a layer material
/// (only for coatings on surfaces that a path uses), and a surface without unambiguous
/// substrate: an inner surface between two segments, a Mirror without material, a ThinElement,
/// Stop or Detector.
///
/// Configurations (ADR 0029, point 5; #165): `configuration` selects a column of the parameter
/// table. Before the poses are composed, every Param bound to a row gets the row's value in that
/// configuration (model::resolve_parameters); without bound Params the system is compiled as it
/// is, so the result does not depend on `configuration` beyond configuration() and
/// configuration_name(). Errors of the table are diagnostics of validate() and come as
/// CompileError; an index >= model::configuration_count(system) is a CompileError with
/// config.unknown at "/configurations" ("" without a configurations section).
///
/// The result holds no references or pointers into `system`, `materials` or `coatings`.
/// @throws CompileError as described above
/// @throws std::invalid_argument for a phase layer with a non-finite coefficient or orientation
///         (only reachable through the API; validate() does not check them, like the
///         coefficients of aspheres; rtt::geom phase constructors, #127)
[[nodiscard]] CompiledSystem compile(const model::System& system,
                                     const material::MaterialLibrary& materials,
                                     const coating::CoatingLibrary& coatings,
                                     std::size_t configuration = 0);

/// compile() with the overlay `perturbation` (ADR 0036, point 1; #254), applied after the
/// parameter table. Without entries, or with identity motions and zero material changes only,
/// the result is bitwise that of compile() without overlay.
///
/// - Motions: every placement reference (ADR 0028, point 2) uses the nominal frames G; a node or
///   surface X gets G'_X = O * G_X * T_X, with O the product of G_Y T_Y G_Y^-1 over its perturbed
///   ancestors Y, outer to inner (see Perturbation). CompiledSurface::to_global, the optic axes
///   of crystals and the axes of ideal polarizers and retarders follow the perturbed frames of
///   their element (the surface motion does not turn the axes, as the surface pose does not).
/// - Node frames: NodeFrame::reference is O * (nominal reference) and NodeFrame::to_global is
///   G'_X, so to_global = reference * to_isometry(pose) * T_X; the to_global of a surface stays
///   bitwise CompiledSurface::to_global. For a node outside a perturbed subtree that is placed
///   relative to a surface inside it, `reference` is the nominal frame of that surface and then
///   differs from the surface's (perturbed) to_global.
/// - Materials: the perturbed segment gets a medium of its own (MaterialChange). The structure
///   comes from the nominal media: whether a Plate is segmented, and whether an inner surface
///   reached from outside is ambiguous, so the path rules do not change.
/// @throws std::invalid_argument for an entry the model has no target for (unknown node name,
///         the root, unknown SurfaceId or element name), a `segment` that is missing for an
///         element with segment materials, given for one with a single material, out of
///         range, or other than 0 for a list of equal materials, and non-finite values (only
///         reachable through the API; rtt-tolerance builds the overlay from a validated tolerance
///         set); and for `object_reference_z` given for a system without object.reference or with
///         the object at infinity, with a length other than the number of configurations, or with
///         a non-finite value (ADR 0038, point 5)
/// @throws CompileError as compile(), and for a material change (at .../material or
///         .../material/<i> of the element) with tolerances.material_unsupported (crystal,
///         element without material, kappa != 0 at d, F or C, a line outside the range of the
///         material), tolerances.abbe_on_const (Delta nu_d != 0 with n_F = n_C) and
///         tolerances.limit_invalid (n_d + Delta n_d <= 1, or nu_d + Delta nu_d <= 0 with
///         n_F != n_C)
[[nodiscard]] CompiledSystem compile(const model::System& system,
                                     const material::MaterialLibrary& materials,
                                     const coating::CoatingLibrary& coatings,
                                     std::size_t configuration,
                                     const Perturbation& perturbation);

/// compile() with an empty CoatingLibrary: a system with a CoatingRef gives a CompileError.
[[nodiscard]] CompiledSystem compile(const model::System& system,
                                     const material::MaterialLibrary& materials,
                                     std::size_t configuration = 0);

/// Index of the configuration called `name` (ADR 0029, point 5; #169).
/// @throws CompileError with config.unknown at "/configurations" ("" without a configurations
///         section), the name in the message, if the system has no configuration of that name
[[nodiscard]] std::size_t configuration_index(const model::System& system, std::string_view name);

/// compile() for the configuration called `configuration` (configuration_index()). A literal 0
/// or any integer chooses the index overloads above: a standard conversion beats the
/// conversion to std::string_view.
/// @throws CompileError as configuration_index() and as compile()
[[nodiscard]] CompiledSystem compile(const model::System& system,
                                     const material::MaterialLibrary& materials,
                                     const coating::CoatingLibrary& coatings,
                                     std::string_view configuration);

/// compile() with an empty CoatingLibrary for the configuration called `configuration`.
[[nodiscard]] CompiledSystem compile(const model::System& system,
                                     const material::MaterialLibrary& materials,
                                     std::string_view configuration);

/// Frames of one node of the model tree (assembly, element or surface) in a compiled system
/// (ADR 0028, points 2 and 6; #169): the global transform of the frame its pose is given in
/// (the parent for an absolute pose, the last surface before it in tree order for
/// relative_to_preceding, the preceding sibling for relative_to_sibling, the identity for the
/// root) and its own global transform, to_global = reference * to_isometry(pose). A GUI needs
/// both to move a node in 3D and to turn an absolute pose into a relative one. With a
/// Perturbation both are the perturbed frames (compile() with overlay, #254).
struct NodeFrame {
  std::string location;       ///< JSON pointer of the node, e.g. "/root/children/1/surfaces/0"
  math::Isometry3 reference;  ///< frame of the pose -> global
  math::Isometry3 to_global;  ///< node coordinates -> global
};

/// Immutable compiled system. Only const access; safe to read from many threads.
class CompiledSystem {
 public:
  /// System wavelengths in um (vacuum), in model order.
  [[nodiscard]] const std::vector<double>& wavelengths_um() const noexcept {
    return wavelengths_um_;
  }
  /// Index of the reference wavelength in wavelengths_um().
  [[nodiscard]] std::uint16_t reference_wavelength() const noexcept { return reference_wl_; }
  /// Weights of the system wavelengths as in the model (model::Wavelength::weight, same order as
  /// wavelengths_um(), not normalised; dimensionless).
  [[nodiscard]] const std::vector<double>& wavelength_weights() const noexcept {
    return wavelength_weights_;
  }

  /// Environment temperature in degC at which the media were evaluated.
  [[nodiscard]] double temperature_c() const noexcept { return temperature_c_; }

  /// System aperture, field points and object space, copied from the model.
  [[nodiscard]] const model::SystemAperture& aperture() const noexcept { return aperture_; }
  [[nodiscard]] const model::FieldSet& fields() const noexcept { return fields_; }
  [[nodiscard]] const model::ObjectSpace& object() const noexcept { return object_; }

  /// Global z of the axial object point in mm, for a finite object (ADR 0038, points 3 and 4):
  /// -object().distance from the global origin (ADR 0031, addendum #279), or with
  /// object.reference z_S - distance, z_S the vertex z of that surface in the nominal frame of
  /// this configuration (or Perturbation::object_reference_z). The one place that turns the
  /// object distance into a position; for an object at infinity the value has no meaning
  /// (callers use at_infinity first).
  [[nodiscard]] double object_z() const noexcept { return object_z_; }

  /// The z_S that object_z() used (mm, global): the nominal vertex z of the surface that
  /// object.reference names, or Perturbation::object_reference_z; none without a reference or
  /// for an object at infinity (ADR 0038). rtt-tolerance fills the overlay of its trials with
  /// the nominal values from here.
  [[nodiscard]] const std::optional<double>& object_reference_z() const noexcept {
    return object_reference_z_;
  }

  /// All surfaces in tree order.
  [[nodiscard]] const std::vector<CompiledSurface>& surfaces() const noexcept { return surfaces_; }
  /// Elements in tree order; CompiledSurface::element indexes this list (#81).
  [[nodiscard]] const std::vector<CompiledElement>& elements() const noexcept { return elements_; }
  /// Distinct media; media()[environment_medium()] is the environment.
  [[nodiscard]] const std::vector<CompiledMedium>& media() const noexcept { return media_; }
  [[nodiscard]] std::uint32_t environment_medium() const noexcept { return 0; }

  /// Coatings used by the surfaces (SurfaceCoating::coating indexes this list), evaluated at
  /// all system wavelengths.
  [[nodiscard]] const std::vector<CompiledCoating>& coatings() const noexcept { return coatings_; }

  /// All paths in model order.
  [[nodiscard]] const std::vector<CompiledPath>& paths() const noexcept { return paths_; }
  /// Path by id.
  /// @throws std::out_of_range if `id` does not belong to this system
  [[nodiscard]] const CompiledPath& path(PathId id) const { return paths_.at(id.index); }

  /// Index of the surface of the system stop (ADR 0039, point 3): the surface named by
  /// aperture.stop, or else the surface of the Stop element; none without either. The stop of a
  /// path is its first event at this surface (stop_event()), whatever its kind.
  [[nodiscard]] std::optional<std::uint32_t> stop_surface() const noexcept { return stop_surface_; }

  /// Index of the surface with this id, if any.
  [[nodiscard]] std::optional<std::uint32_t> find_surface(const model::SurfaceId& id) const;
  /// Id of the path with this name, if any.
  [[nodiscard]] std::optional<PathId> find_path(std::string_view name) const;

  /// Index of the configuration this system was compiled for (ADR 0029, point 5); 0 for the
  /// nominal configuration of a system without configurations.
  [[nodiscard]] std::size_t configuration() const noexcept { return configuration_; }
  /// Name of that configuration; empty for the nominal configuration (no configurations section).
  [[nodiscard]] const std::string& configuration_name() const noexcept {
    return configuration_name_;
  }

  /// Frames of every node of the model tree in pre-order (an element before its surfaces), see
  /// NodeFrame; the to_global of a surface equals CompiledSurface::to_global bit for bit.
  [[nodiscard]] const std::vector<NodeFrame>& node_frames() const noexcept { return node_frames_; }

  /// Warnings found while compiling (model::validate), with code and JSON pointer (ADR 0022).
  [[nodiscard]] const std::vector<model::Diagnostic>& diagnostics() const noexcept {
    return diagnostics_;
  }

 private:
  friend CompiledSystem compile(const model::System&,
                                const material::MaterialLibrary&,
                                const coating::CoatingLibrary&,
                                std::size_t,
                                const Perturbation&);
  CompiledSystem() = default;

  std::vector<double> wavelengths_um_;
  std::vector<double> wavelength_weights_;
  std::uint16_t reference_wl_ = 0;
  double temperature_c_ = 20.0;
  model::SystemAperture aperture_;
  model::FieldSet fields_;
  model::ObjectSpace object_;
  std::vector<CompiledSurface> surfaces_;
  std::optional<std::uint32_t> stop_surface_;
  double object_z_ = 0.0;                     ///< see object_z()
  std::optional<double> object_reference_z_;  ///< see object_reference_z()
  std::vector<CompiledElement> elements_;
  std::vector<CompiledMedium> media_;
  std::vector<CompiledCoating> coatings_;
  std::vector<CompiledPath> paths_;
  std::vector<model::Diagnostic> diagnostics_;
  std::size_t configuration_ = 0;
  std::string configuration_name_;
  std::vector<NodeFrame> node_frames_;
};

/// The stop event of `path` (ADR 0039, point 3): the index of its first event at
/// CompiledSystem::stop_surface(), whatever the event kind (a double pass or a ghost path may
/// meet the stop surface more than once); none if the system has no stop or the path does not
/// visit it. require_stop() throws for none.
/// @throws std::out_of_range if `path` does not belong to `system`
[[nodiscard]] std::optional<std::size_t> stop_event(const CompiledSystem& system, PathId path);

}  // namespace rtt::compile
