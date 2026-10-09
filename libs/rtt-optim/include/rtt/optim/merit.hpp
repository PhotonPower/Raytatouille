#pragma once

/// @file merit.hpp
/// Evaluation of the merit function of a system (ADR 0030, points 1-5, 10; #167, #168): the
/// operands and generators of System::optimization at given values of the variables. One
/// evaluation sets the values in a copy of the system, compiles every configuration the operands
/// and generators use (ADR 0029, point 5) and calls exactly the existing analysis of each
/// operand and the generators of generators.hpp.

#include <cstddef>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include "rtt/coating/catalog.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/system.hpp"
#include "rtt/model/validate.hpp"
#include "rtt/optim/generators.hpp"
#include "rtt/optim/variables.hpp"

namespace rtt::optim {

/// Thrown at the start of an optimization for an input the run cannot use, with diagnostics
/// carrying stable codes (ADR 0022): `optim.no_variables` (no variable in the system),
/// `optim.no_operands` (no operand and no generator) and `merit.operand_unsupported` (an
/// operand on a path its analysis rejects, ADR 0030, points 3 and 10), with JSON pointers into
/// the system file.
class OptimError : public std::invalid_argument {
 public:
  explicit OptimError(std::vector<model::Diagnostic> diagnostics);

  [[nodiscard]] const std::vector<model::Diagnostic>& diagnostics() const noexcept {
    return diagnostics_;
  }

 private:
  std::vector<model::Diagnostic> diagnostics_;
};

/// One evaluation of the merit function.
struct MeritEvaluation {
  /// Value of each operand in file order, in the unit of the operand (mm, waves or
  /// dimensionless; see rtt/model/optimization.hpp); NaN if it is not defined at this point.
  std::vector<double> values;
  /// Residual sqrt(weight) (value - target) of each operand, same order (ADR 0030, point 1).
  /// An undefined value gives a NaN residual also for weight 0: an operand that is only watched
  /// must still be defined, or the evaluation is invalid (ADR 0030, point 10). Then the
  /// residuals of the generators in file order, each block as generator_residuals() (ADR 0030,
  /// addendum #168, point 13); MeritFunction::generator_sizes() gives the block lengths.
  std::vector<double> residuals;
  /// Why an operand has no value (empty string if it has one), same order: an afocal path for
  /// efl/bfl, no paraxial working F-number, a ray that does not arrive at its surface, ...
  std::vector<std::string> undefined;
  /// Warnings of the analyses ("rays.lost", "stop.clips_beam"), in operand order, each once.
  std::vector<model::Diagnostic> warnings;
  /// Rays and mean square of each generator, file order (ADR 0030, addendum #168).
  std::vector<GeneratorStats> generators;

  /// True if every residual is finite.
  [[nodiscard]] bool valid() const noexcept;
};

/// The merit function of a system: its variables (collect_variables) and its operands.
///
/// Operands (ADR 0030, point 3), each in configuration `common.configuration` (none:
/// configuration 0):
/// - efl, bfl, magnification: paraxial::first_order(cs, path, wavelength); efl and bfl are not
///   defined for an afocal path, magnification not for an object at infinity (a start error).
/// - image_fnumber: paraxial::prescription(cs, path, wavelength).paraxial_working_f_number.
/// - ray_x, ray_y: one ray from trace::make_rays (SinglePupilPoint, real aiming) traced with
///   recording; the coordinate of its hit point at the event of `surface` (occurrence) in the
///   local coordinates of that surface, mm. Not defined if the ray is not Alive after that
///   event (lost before or at it); a loss after it does not matter.
/// - spot_rms: analysis::spot with HexapolarPupil{rings} and real aiming; stats.rms_centroid or
///   stats.rms_chief, mm. Polychromatic: all system wavelengths with their weights.
/// - opd_rms: analysis::opd_map with OpdOptions{grid}; rms, waves at the reference wavelength.
/// - param_value: the value of the table row in the configuration (model::evaluate_parameters).
/// The wavelength of an operand is an index into System::wavelengths; none: the reference.
/// Generators rms_spot and rms_wavefront (ADR 0030, point 4; generators.hpp), each in
/// configuration `configuration` (none: configuration 0), append their residuals after those
/// of the operands.
class MeritFunction {
 public:
  /// Prepares the evaluation and checks the operands at the start (ADR 0030, point 10).
  /// @param system    a system that passes model::validate (also its section optimization); kept
  ///                  by reference, so it must outlive this object (a temporary is rejected)
  /// @param materials material library for compile; must outlive this object
  /// @param coatings  coating library for compile, or nullptr; must outlive this object
  /// @throws compile::CompileError if the start system does not compile in a used configuration
  /// @throws OptimError with merit.operand_unsupported for a magnification operand with the
  ///         object at infinity, an operand other than param_value or a generator on a path
  ///         that paraxial::first_order rejects (crystal, order != 0, not rotationally
  ///         symmetric, stop not circular; for a generator at each chosen wavelength), or a
  ///         generator whose chosen field or wavelength weights sum to 0 (ADR 0030, addendum
  ///         #168), pointer /optimization/operands/i or /optimization/generators/g
  MeritFunction(const model::System& system,
                const material::MaterialLibrary& materials,
                const coating::CoatingLibrary* coatings);
  /// No temporary system: the object keeps a reference to it.
  MeritFunction(model::System&& system,
                const material::MaterialLibrary& materials,
                const coating::CoatingLibrary* coatings) = delete;

  /// The variables, in the order of ADR 0030, point 5.
  [[nodiscard]] const std::vector<Variable>& variables() const noexcept { return variables_; }
  /// Start values of the variables (Variable::start), same order.
  [[nodiscard]] std::vector<double> start() const;
  /// Number of residuals m (one per operand, then generator_sizes()), fixed for the run.
  [[nodiscard]] std::size_t size() const noexcept;
  /// Number of residuals of each generator (generator_size()), file order.
  [[nodiscard]] const std::vector<std::size_t>& generator_sizes() const noexcept {
    return generator_sizes_;
  }

  /// Evaluates the operands with the variables at `values`. Thread-safe and deterministic:
  /// nothing is shared between calls except the read-only input and the libraries.
  /// @param values one value per variable, in the unit of the field or row
  /// @throws as compile() and the analyses of the operands (ParaxialError, AnalysisError,
  ///         NoStopError, std::invalid_argument); an undefined value is NaN, not an exception
  [[nodiscard]] MeritEvaluation evaluate(std::span<const double> values) const;

 private:
  const model::System* system_;
  const material::MaterialLibrary* materials_;
  const coating::CoatingLibrary* coatings_;
  std::vector<Variable> variables_;
  std::vector<std::size_t> configurations_;            ///< configuration of each operand
  std::vector<std::size_t> compiled_;                  ///< configurations to compile, ascending
  std::vector<std::optional<std::size_t>> rows_;       ///< table row of a param_value operand
  std::vector<std::size_t> generator_configurations_;  ///< configuration of each generator
  std::vector<std::size_t> generator_sizes_;           ///< residuals of each generator
  std::size_t size_ = 0;                               ///< operands plus generator residuals
};

}  // namespace rtt::optim
