#pragma once

/// @file overlay.hpp
/// The merit function of a perturbed trial (ADR 0030, addendum #256; ADR 0036, point 2): the
/// overlay of a tolerance trial (compile::Perturbation, #254) applied at every compile, and the
/// nominal basis, the paraxial data of the nominal system that the trial cannot have itself (it
/// is not rotationally symmetric, so rtt::paraxial::first_order rejects it).
///
/// Conventions: lengths in mm, global coordinates (right-handed, optical axis +z).

#include <cstddef>
#include <memory>
#include <utility>

#include "rtt/coating/catalog.hpp"
#include "rtt/compile/perturbation.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/system.hpp"

namespace rtt::optim {

namespace detail {
struct BasisData;    // overlay_basis.hpp (private)
struct BasisAccess;  // overlay_basis.hpp (private)
}  // namespace detail

/// The nominal side of a tolerance trial (ADR 0036, point 2), computed once from the nominal
/// system for every configuration and path:
/// - the analysis::NominalReference (aiming data per wavelength, exit pupil; #265);
/// - per wavelength the first-order data of the nominal system (signs of efl and magnification,
///   whether a path is supported) and the last vertex of the path (bfl);
/// - for efl and bfl the aiming data of the variant "object at infinity" (addendum #256).
/// An immutable value with shared content: cheap to copy, safe to read from many threads.
/// Computations that fail are kept and reported where an operand needs them.
class NominalBasis {
 public:
  /// Number of configurations of the nominal system (ADR 0029, point 5).
  [[nodiscard]] std::size_t configurations() const noexcept;

 private:
  friend NominalBasis nominal_basis(const model::System& nominal,
                                    const material::MaterialLibrary& materials,
                                    const coating::CoatingLibrary* coatings);
  friend struct detail::BasisAccess;
  explicit NominalBasis(std::shared_ptr<const detail::BasisData> data) : data_(std::move(data)) {}
  std::shared_ptr<const detail::BasisData> data_;
};

/// The nominal basis of `nominal`: every configuration compiled once (without overlay), and for
/// each path and wavelength the data above.
/// @param nominal   a valid system; not kept
/// @param materials material library for compile
/// @param coatings  coating library for compile, or nullptr
/// @throws compile::CompileError if a configuration does not compile
[[nodiscard]] NominalBasis nominal_basis(const model::System& nominal,
                                         const material::MaterialLibrary& materials,
                                         const coating::CoatingLibrary* coatings);

/// The overlay of one trial (ADR 0030, addendum #256): the perturbation applied at every compile
/// of the merit function, the nominal basis for aiming, reference spheres and the paraxial
/// operands.
struct TrialOverlay {
  compile::Perturbation perturbation;
  NominalBasis nominal;
  /// delta of the parabasal rays of the paraxial operands, in normalised pupil coordinates
  /// (ADR 0036, point 2: 1e-3, as the mapping M of ADR 0032, point 4). Changed only by tests
  /// (Richardson bound, addendum #256).
  double parabasal_step = 1e-3;
};

}  // namespace rtt::optim
