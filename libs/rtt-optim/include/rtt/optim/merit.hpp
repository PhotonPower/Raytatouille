#pragma once

/// @file merit.hpp
/// Evaluation of the merit function of a system (ADR 0030, points 1-3, 5, 10; #167): the
/// operands of System::optimization at given values of the variables. One evaluation sets the
/// values in a copy of the system, compiles every configuration the operands use (ADR 0029,
/// point 5) and calls exactly the existing analysis of each operand.

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
#include "rtt/optim/variables.hpp"

namespace rtt::optim {

/// Thrown at the start of an optimization for an input the run cannot use, with diagnostics
/// carrying stable codes (ADR 0022): `optim.no_variables` (no variable in the system) and
/// `merit.operand_unsupported` (an operand on a path its analysis rejects, ADR 0030, points 3
/// and 10), with JSON pointers into the system file.
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
  std::vector<double> residuals;
  /// Why an operand has no value (empty string if it has one), same order: an afocal path for
  /// efl/bfl, no paraxial working F-number, a ray that does not arrive at its surface, ...
  std::vector<std::string> undefined;
  /// Warnings of the analyses ("rays.lost", "stop.clips_beam"), in operand order, each once.
  std::vector<model::Diagnostic> warnings;

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
/// Generators (ADR 0030, point 4) come with #168; until then the constructor rejects them.
class MeritFunction {
 public:
  /// Prepares the evaluation and checks the operands at the start (ADR 0030, point 10).
  /// @param system    a system that passes model::validate (also its section optimization)
  /// @param materials material library for compile; must outlive this object
  /// @param coatings  coating library for compile, or nullptr; must outlive this object
  /// @throws compile::CompileError if the start system does not compile in a used configuration
  /// @throws OptimError with merit.operand_unsupported for a generator (until #168), a
  ///         magnification operand with the object at infinity, or an operand other than
  ///         param_value on a path that paraxial::first_order rejects (crystal, order != 0, not
  ///         rotationally symmetric, stop not circular), pointer /optimization/operands/i
  MeritFunction(const model::System& system,
                const material::MaterialLibrary& materials,
                const coating::CoatingLibrary* coatings);

  /// The variables, in the order of ADR 0030, point 5.
  [[nodiscard]] const std::vector<Variable>& variables() const noexcept { return variables_; }
  /// Start values of the variables (Variable::start), same order.
  [[nodiscard]] std::vector<double> start() const;
  /// Number of residuals m (one per operand), fixed for the run.
  [[nodiscard]] std::size_t size() const noexcept;

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
  std::vector<std::size_t> configurations_;       ///< configuration of each operand
  std::vector<std::size_t> compiled_;             ///< configurations to compile, ascending
  std::vector<std::optional<std::size_t>> rows_;  ///< table row of a param_value operand
};

}  // namespace rtt::optim
