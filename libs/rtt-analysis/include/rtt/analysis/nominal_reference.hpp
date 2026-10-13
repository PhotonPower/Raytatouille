#pragma once

/// @file nominal_reference.hpp
/// Nominal reference of perturbed trials (ADR 0036, point 2, maintainer's decision; ADR 0032,
/// addendum #265): the paraxial data that spot() and the OPD analyses need, taken from the
/// nominal system, while the rays run through a perturbed trial that is not rotationally
/// symmetric (rtt::paraxial::first_order rejects it).
///
/// Conventions: positions in mm in global coordinates (right-handed, optical axis +z).

#include <cstdint>
#include <memory>
#include <optional>
#include <utility>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/math/types.hpp"
#include "rtt/trace/sources.hpp"

namespace rtt::analysis {

namespace detail {
struct NominalData;  // nominal_reference.cpp
}  // namespace detail

/// Paraxial data of one path of a compiled (nominal) system for the analyses: per system
/// wavelength the aiming data (trace::AimReference, ADR 0007, addendum #252) and the centre of
/// the paraxial exit pupil at the reference wavelength. An immutable value with shared content:
/// cheap to copy (the options of spot() and the OPD analyses hold it by value) and safe to read
/// from many threads. Created with nominal_reference(); for a system with configurations, from
/// the system compiled for the configuration of the analysis.
class NominalReference {
 public:
  /// Path the data belong to.
  [[nodiscard]] compile::PathId path() const noexcept;
  /// Aiming data at `wavelength` (index into the wavelengths of the nominal system).
  /// @throws std::invalid_argument for an index out of range; otherwise the exception that
  ///         trace::aim_reference() gave for this wavelength (same type), if it failed
  [[nodiscard]] const trace::AimReference& aim(std::uint16_t wavelength) const;
  /// Centre E of the paraxial exit pupil of the nominal system at its reference wavelength,
  /// global, E = (0, 0, z_XP) mm; none for an exit pupil at infinity (image-space telecentric).
  /// @throws the exception of the paraxial computation (rtt::compile::NoStopError,
  ///         rtt::paraxial::ParaxialError), if it failed
  [[nodiscard]] std::optional<math::Vec3> exit_pupil() const;

 private:
  friend NominalReference nominal_reference(const compile::CompiledSystem& nominal,
                                            compile::PathId path);
  explicit NominalReference(std::shared_ptr<const detail::NominalData> data)
      : data_(std::move(data)) {}
  std::shared_ptr<const detail::NominalData> data_;
};

/// The nominal reference of `path` of `nominal`: trace::aim_reference() at every system
/// wavelength and the paraxial exit pupil at the reference wavelength, computed once. A
/// computation that fails is kept and its exception thrown where the analysis needs it, so an
/// analysis with the reference of the same system fails exactly as without reference.
/// @throws std::invalid_argument if the path does not exist
[[nodiscard]] NominalReference nominal_reference(const compile::CompiledSystem& nominal,
                                                 compile::PathId path);

}  // namespace rtt::analysis
