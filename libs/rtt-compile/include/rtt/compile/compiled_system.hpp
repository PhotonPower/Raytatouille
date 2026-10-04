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

#include "rtt/geom/asphere.hpp"
#include "rtt/geom/conic.hpp"
#include "rtt/geom/plane.hpp"
#include "rtt/material/material.hpp"
#include "rtt/math/isometry.hpp"
#include "rtt/math/types.hpp"
#include "rtt/model/model.hpp"

namespace rtt::compile {

/// Resolved base shape of a surface in its local coordinates. Supported in M1: plane, conic
/// and even asphere; further shapes are added to this variant as rtt-geom provides them.
using CompiledShape =
    std::variant<geom::Plane<double>, geom::Conic<double>, geom::EvenAsphere<double>>;

/// One surface with everything the tracer needs, independent of the model.
struct CompiledSurface {
  model::SurfaceId id;  ///< stable surface id from the model, e.g. "L1.S1"
  model::ElementKind element_kind = model::ElementKind::Lens;  ///< kind of the owning element
  std::string element_name;                                    ///< name of the owning element
  /// Local surface coordinates -> global coordinates (assembly -> element -> surface poses
  /// chained, see rtt::math::Isometry3::from_pose for the convention).
  math::Isometry3 to_global;
  /// Global coordinates -> local surface coordinates (inverse of to_global).
  math::Isometry3 to_local;
  CompiledShape shape = geom::Plane<double>{};
  std::optional<model::Aperture> aperture;  ///< in local coordinates, mm; none = unbounded
  std::vector<model::PhaseLayer> phases;    ///< copied as values; evaluated from M4 on
  model::Interaction interaction = model::Fresnel{};
};

/// A homogeneous medium evaluated at all system wavelengths.
struct CompiledMedium {
  std::string reference;  ///< material reference as in the model, e.g. "AIR", "CONST:1.5168"
  /// Complex index n + i*kappa per system wavelength (same order as wavelengths_um()),
  /// evaluated at the environment temperature.
  std::vector<math::Complex> index;
};

/// One step of a compiled path.
struct CompiledEvent {
  std::uint32_t surface = 0;  ///< index into CompiledSystem::surfaces()
  model::EventKind kind = model::EventKind::Refract;
  int order = 0;                    ///< diffraction order (Diffract only)
  std::uint32_t medium_before = 0;  ///< index into CompiledSystem::media()
  std::uint32_t medium_after = 0;   ///< index into CompiledSystem::media()
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
/// - Runs model::validate; any error diagnostic throws CompileError. Warnings are not
///   reported by compile(); call model::validate directly to see them.
/// - Flattens the tree: the global transform of each surface is root pose * assembly poses *
///   element pose * surface pose. Pickups are not evaluated before M5; their value is used.
/// - Resolves every material (environment medium, element materials) with `materials` and
///   evaluates it once per system wavelength at environment.temperature_c (ADR 0014).
/// - Builds every path. An automatic path visits all surfaces in tree order with Refract for
///   Lens and Plate, Reflect for Mirror and Transmit for Stop, Detector and ThinElement.
/// - Media along a path: the ray starts in the environment medium. Refract, Ordinary and
///   Extraordinary at a surface of an element with a material switch between the inside of
///   that element and the environment; all other events, and events at elements without
///   material, keep the medium. Entering element B while inside element A goes from A's
///   medium to B's, and leaving B then goes to the environment; that case is not settled
///   yet and is part of the question on cemented groups (#14).
///
/// Not supported yet (CompileError): Zernike sag terms (M8), and Lens or Plate elements with
/// more than 2 surfaces on an automatic path (cemented groups, see #14).
/// Also a CompileError: a Mirror with substrate material and more than one surface on an
/// automatic path (Mangin mirror; its front surface refracts, so it needs an explicit path
/// Refract, Reflect, Refract).
///
/// The result holds no references or pointers into `system` or `materials`.
/// @throws CompileError as described above
[[nodiscard]] CompiledSystem compile(const model::System& system,
                                     const material::MaterialLibrary& materials);

/// Immutable compiled system. Only const access; safe to read from many threads.
class CompiledSystem {
 public:
  /// System wavelengths in um (vacuum), in model order.
  [[nodiscard]] const std::vector<double>& wavelengths_um() const noexcept {
    return wavelengths_um_;
  }
  /// Index of the reference wavelength in wavelengths_um().
  [[nodiscard]] std::uint16_t reference_wavelength() const noexcept { return reference_wl_; }

  /// Environment temperature in degC at which the media were evaluated.
  [[nodiscard]] double temperature_c() const noexcept { return temperature_c_; }

  /// System aperture, field points and object space, copied from the model.
  [[nodiscard]] const model::SystemAperture& aperture() const noexcept { return aperture_; }
  [[nodiscard]] const model::FieldSet& fields() const noexcept { return fields_; }
  [[nodiscard]] const model::ObjectSpace& object() const noexcept { return object_; }

  /// All surfaces in tree order.
  [[nodiscard]] const std::vector<CompiledSurface>& surfaces() const noexcept { return surfaces_; }
  /// Distinct media; media()[environment_medium()] is the environment.
  [[nodiscard]] const std::vector<CompiledMedium>& media() const noexcept { return media_; }
  [[nodiscard]] std::uint32_t environment_medium() const noexcept { return 0; }

  /// All paths in model order.
  [[nodiscard]] const std::vector<CompiledPath>& paths() const noexcept { return paths_; }
  /// Path by id.
  /// @throws std::out_of_range if `id` does not belong to this system
  [[nodiscard]] const CompiledPath& path(PathId id) const { return paths_.at(id.index); }

  /// Index of the surface with this id, if any.
  [[nodiscard]] std::optional<std::uint32_t> find_surface(const model::SurfaceId& id) const;
  /// Id of the path with this name, if any.
  [[nodiscard]] std::optional<PathId> find_path(std::string_view name) const;

 private:
  friend CompiledSystem compile(const model::System&, const material::MaterialLibrary&);
  CompiledSystem() = default;

  std::vector<double> wavelengths_um_;
  std::uint16_t reference_wl_ = 0;
  double temperature_c_ = 20.0;
  model::SystemAperture aperture_;
  model::FieldSet fields_;
  model::ObjectSpace object_;
  std::vector<CompiledSurface> surfaces_;
  std::vector<CompiledMedium> media_;
  std::vector<CompiledPath> paths_;
};

}  // namespace rtt::compile
