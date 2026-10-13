#pragma once

/// @file overlay_basis.hpp
/// Content of a NominalBasis (overlay.hpp; ADR 0030, addendum #256), private to rtt-optim.

#include <cstddef>
#include <cstdint>
#include <exception>
#include <optional>
#include <string>
#include <vector>

#include "rtt/analysis/nominal_reference.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/model/system.hpp"
#include "rtt/optim/overlay.hpp"
#include "rtt/trace/sources.hpp"

namespace rtt::optim::detail {

/// What the paraxial operands need of the nominal system at one wavelength of one path.
struct WavelengthData {
  /// The error of rtt::paraxial::first_order of the nominal path, empty if it succeeded; the
  /// path is then not supported by any operand but param_value (ADR 0030, point 3).
  std::string unsupported;
  std::optional<double> efl;            ///< nominal EFL, mm (sign of efl; none: afocal)
  std::optional<double> magnification;  ///< nominal lateral magnification (sign of magnification)
  /// Event index of the last vertex of the path (the last event that changes the ray, as
  /// rtt::paraxial::first_order: reflection, change of index or an ideal lens); none if no event
  /// changes the ray. Start of bfl (addendum #256).
  std::optional<std::size_t> last_vertex;
  /// Aiming data of the axis-parallel bundle of efl and bfl: those of the variant "object at
  /// infinity" (variant_of()) for a finite object, else none (the path's own NominalReference).
  std::optional<trace::AimReference> infinite_aim;
  /// Why the axis-parallel bundle cannot be aimed (e.g. an entrance pupil at infinity with the
  /// object at infinity: object-space telecentric); empty if it can.
  std::string infinite_unsupported;
};

/// One path in one configuration.
struct PathData {
  std::string name;
  analysis::NominalReference reference;
  std::vector<WavelengthData> wavelengths;
};

/// One configuration.
struct ConfigurationData {
  bool object_at_infinity = true;
  std::vector<PathData> paths;  ///< path order of the compiled system
};

struct BasisData {
  std::vector<ConfigurationData> configurations;
  std::size_t wavelengths = 0;
};

/// Access to the content of a NominalBasis.
struct BasisAccess {
  [[nodiscard]] static const BasisData& data(const NominalBasis& basis) noexcept {
    return *basis.data_;
  }
};

/// The variant of `system` for the axis-parallel bundle of efl and bfl (addendum #256): object
/// at infinity, one field at 0 degree, the system aperture StopSize (the stop of each path as
/// defined). Everything else, the tree and the parameter table included, is the input.
[[nodiscard]] model::System variant_of(const model::System& system);

}  // namespace rtt::optim::detail
