#include "rtt/optim/merit.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "rtt/analysis/opd.hpp"
#include "rtt/analysis/spot.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/diagnostics/codes.hpp"
#include "rtt/math/types.hpp"
#include "rtt/model/optimization.hpp"
#include "rtt/model/parameters.hpp"
#include "rtt/paraxial/paraxial.hpp"
#include "rtt/paraxial/prescription.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/ray_paths.hpp"
#include "rtt/trace/sequential.hpp"
#include "rtt/trace/sources.hpp"

namespace rtt::optim {

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

model::Diagnostic diagnostic(diagnostics::DiagnosticCode code,
                             std::string location,
                             std::string message) {
  return {code.severity(), std::move(location), std::move(message), std::string(code.str())};
}

std::string what_of(const std::vector<model::Diagnostic>& diagnostics) {
  std::string out = "optimize:";
  for (const model::Diagnostic& d : diagnostics) {
    out += " " + d.code + " at '" + d.location + "': " + d.message + ";";
  }
  return out;
}

std::string operand_pointer(std::size_t i) {
  return "/optimization/operands/" + std::to_string(i);
}

// Accessors over the operand variant with get_if (no std::visit returning a reference: GCC 13
// warns -Wdangling-reference for that at the call site).
const model::OperandCommon& common_of(const model::Operand& op) {
  if (const auto* o = std::get_if<model::FirstOrderOperand>(&op)) return o->common;
  if (const auto* o = std::get_if<model::RayOperand>(&op)) return o->common;
  if (const auto* o = std::get_if<model::SpotRmsOperand>(&op)) return o->common;
  if (const auto* o = std::get_if<model::OpdRmsOperand>(&op)) return o->common;
  return std::get<model::ParamValueOperand>(op).common;
}

/// Name of the path of an operand; nullptr for param_value.
const std::string* path_of(const model::Operand& op) {
  if (const auto* o = std::get_if<model::FirstOrderOperand>(&op)) return &o->path;
  if (const auto* o = std::get_if<model::RayOperand>(&op)) return &o->path;
  if (const auto* o = std::get_if<model::SpotRmsOperand>(&op)) return &o->path;
  if (const auto* o = std::get_if<model::OpdRmsOperand>(&op)) return &o->path;
  return nullptr;
}

/// Wavelength index of an operand for the support check (none: the reference wavelength).
std::optional<std::uint16_t> wavelength_of(const model::Operand& op) {
  if (const auto* o = std::get_if<model::FirstOrderOperand>(&op)) return o->wavelength;
  if (const auto* o = std::get_if<model::RayOperand>(&op)) return o->wavelength;
  if (const auto* o = std::get_if<model::SpotRmsOperand>(&op)) return o->wavelength;
  if (const auto* o = std::get_if<model::OpdRmsOperand>(&op)) return o->wavelength;
  return std::nullopt;
}

compile::CompiledSystem compile_one(const model::System& system,
                                    const material::MaterialLibrary& materials,
                                    const coating::CoatingLibrary* coatings,
                                    std::size_t configuration) {
  return coatings != nullptr ? compile::compile(system, materials, *coatings, configuration)
                             : compile::compile(system, materials, configuration);
}

compile::PathId path_id(const compile::CompiledSystem& cs, const std::string& name) {
  const std::optional<compile::PathId> id = cs.find_path(name);
  if (!id) throw std::invalid_argument("optimize: unknown path '" + name + "'");
  return *id;
}

/// Value and, if undefined, the reason.
struct Value {
  double value = kNaN;
  std::string undefined;
};

Value defined_or(const std::optional<double>& v, const char* reason) {
  if (v) return {*v, {}};
  return {kNaN, reason};
}

Value first_order_value(const compile::CompiledSystem& cs, const model::FirstOrderOperand& op) {
  const compile::PathId path = path_id(cs, op.path);
  const std::uint16_t wl = op.wavelength.value_or(cs.reference_wavelength());
  switch (op.quantity) {
    case model::FirstOrderQuantity::Efl:
      return defined_or(paraxial::first_order(cs, path, wl).efl, "the path is afocal (no EFL)");
    case model::FirstOrderQuantity::Bfl:
      return defined_or(paraxial::first_order(cs, path, wl).bfl, "the path is afocal (no BFL)");
    case model::FirstOrderQuantity::Magnification:
      return defined_or(paraxial::first_order(cs, path, wl).lateral_magnification,
                        "no paraxial magnification");
    case model::FirstOrderQuantity::ImageFNumber:
      return defined_or(paraxial::prescription(cs, path, wl).paraxial_working_f_number,
                        "no paraxial working F-number (no marginal ray, or collimated output)");
  }
  throw std::invalid_argument("optimize: invalid first-order quantity");
}

/// ADR 0030, point 3: one real ray, its hit point at the event of `surface` in local
/// coordinates of that surface. The ray is traced along the whole path with recording
/// (rtt/trace/ray_paths.hpp): slot e + 1 is the state after event e.
Value ray_value(const compile::CompiledSystem& cs, const model::RayOperand& op) {
  const compile::PathId path = path_id(cs, op.path);
  const std::uint16_t wl = op.wavelength.value_or(cs.reference_wavelength());
  const std::optional<std::uint32_t> surface = cs.find_surface(op.surface);
  if (!surface) throw std::invalid_argument("optimize: unknown surface of a ray operand");
  const std::vector<compile::CompiledEvent>& events = cs.path(path).events;
  const std::uint32_t occurrence = op.occurrence.value_or(0);
  std::optional<std::size_t> event;
  std::uint32_t seen = 0;
  for (std::size_t e = 0; e < events.size() && !event; ++e) {
    if (events[e].surface == *surface && seen++ == occurrence) event = e;
  }
  if (!event) throw std::invalid_argument("optimize: the path does not meet the ray surface");
  const std::uint16_t field = op.field;
  trace::RayBatch rays = trace::make_rays(cs, path, std::span<const std::uint16_t>(&field, 1), wl,
                                          trace::SinglePupilPoint{op.px, op.py});
  trace::RayPaths recorded;
  static_cast<void>(trace::SequentialTracer().trace(cs, path, rays, recorded));
  const std::size_t slot = *event + 1;
  if (recorded.count[0] <= slot || recorded.status[slot] != trace::RayStatus::Alive) {
    return {kNaN, "the ray does not pass the surface (lost before or at it, or not aimed)"};
  }
  const math::Vec3 local =
      cs.surfaces()[*surface].to_local.apply_point(recorded.position_at(0, slot));
  return {op.coordinate == model::RayCoordinate::X ? local.x() : local.y(), {}};
}

void add_warnings(std::vector<model::Diagnostic>& out, const std::vector<model::Diagnostic>& in) {
  for (const model::Diagnostic& d : in) {
    const bool known = std::any_of(out.begin(), out.end(), [&d](const model::Diagnostic& o) {
      return o.code == d.code && o.location == d.location && o.message == d.message;
    });
    if (!known) out.push_back(d);
  }
}

}  // namespace

OptimError::OptimError(std::vector<model::Diagnostic> diagnostics)
    : std::invalid_argument(what_of(diagnostics)), diagnostics_(std::move(diagnostics)) {}

bool MeritEvaluation::valid() const noexcept {
  return std::all_of(residuals.begin(), residuals.end(), [](double r) { return std::isfinite(r); });
}

MeritFunction::MeritFunction(const model::System& system,
                             const material::MaterialLibrary& materials,
                             const coating::CoatingLibrary* coatings)
    : system_(&system),
      materials_(&materials),
      coatings_(coatings),
      variables_(collect_variables(system)) {
  // The model first, also when no operand needs compile (only param_value operands).
  std::vector<model::Diagnostic> errors = model::validate(system);
  std::erase_if(errors,
                [](const model::Diagnostic& d) { return d.severity != model::Severity::Error; });
  if (!errors.empty()) throw compile::CompileError(std::move(errors));

  const model::Optimization& merit = system.optimization;
  std::vector<model::Diagnostic> unsupported;
  // Generators come with #168; until then they would be ignored, so they are rejected.
  for (std::size_t g = 0; g < merit.generators.size(); ++g) {
    unsupported.push_back(diagnostic("merit.operand_unsupported",
                                     "/optimization/generators/" + std::to_string(g),
                                     "generators are not evaluated yet (#168)"));
  }
  for (std::size_t i = 0; i < merit.operands.size(); ++i) {
    const model::Operand& op = merit.operands[i];
    const model::OperandCommon& common = common_of(op);
    std::size_t k = 0;
    if (common.configuration) {
      const std::optional<std::size_t> found =
          model::find_configuration(system, *common.configuration);
      if (!found) throw std::invalid_argument("optimize: unknown configuration");
      k = *found;
    }
    configurations_.push_back(k);
    if (const auto* p = std::get_if<model::ParamValueOperand>(&op)) {
      const std::optional<std::size_t> row = model::find_parameter(system, p->parameter);
      if (!row) throw std::invalid_argument("optimize: unknown parameter row");
      rows_.push_back(row);
      continue;
    }
    rows_.push_back(std::nullopt);
    compiled_.push_back(k);
    const auto* f = std::get_if<model::FirstOrderOperand>(&op);
    if (f != nullptr && f->quantity == model::FirstOrderQuantity::Magnification &&
        system.object.at_infinity) {
      // ADR 0030, point 3 (#167, Q3): the lateral magnification needs a finite object.
      unsupported.push_back(diagnostic("merit.operand_unsupported", operand_pointer(i),
                                       "magnification with the object at infinity"));
    }
  }
  std::sort(compiled_.begin(), compiled_.end());
  compiled_.erase(std::unique(compiled_.begin(), compiled_.end()), compiled_.end());

  // Paths that the analyses reject (ADR 0030, point 3): first_order at the start, which every
  // operand other than param_value needs (paraxial values, aiming).
  for (const std::size_t k : compiled_) {
    const compile::CompiledSystem cs = compile_one(system, materials, coatings, k);
    for (std::size_t i = 0; i < merit.operands.size(); ++i) {
      if (rows_[i] || configurations_[i] != k) continue;
      const model::Operand& op = merit.operands[i];
      const std::uint16_t wl = wavelength_of(op).value_or(cs.reference_wavelength());
      try {
        static_cast<void>(paraxial::first_order(cs, path_id(cs, *path_of(op)), wl));
      } catch (const paraxial::ParaxialError& e) {
        unsupported.push_back(diagnostic("merit.operand_unsupported", operand_pointer(i),
                                         std::string("the path is not supported: ") + e.what()));
      }
    }
  }
  if (!unsupported.empty()) throw OptimError(std::move(unsupported));
}

std::vector<double> MeritFunction::start() const {
  std::vector<double> out;
  out.reserve(variables_.size());
  for (const Variable& v : variables_) out.push_back(v.start);
  return out;
}

std::size_t MeritFunction::size() const noexcept {
  return system_->optimization.operands.size();
}

MeritEvaluation MeritFunction::evaluate(std::span<const double> values) const {
  const model::System system = with_values(*system_, variables_, values);
  // Every used configuration once, in ascending order (ADR 0030, point 13).
  std::vector<std::optional<compile::CompiledSystem>> compiled(model::configuration_count(system));
  for (const std::size_t k : compiled_) {
    compiled[k].emplace(compile_one(system, *materials_, coatings_, k));
  }
  std::optional<model::ParameterValues> table;

  const std::vector<model::Operand>& operands = system_->optimization.operands;
  MeritEvaluation out;
  out.values.reserve(operands.size());
  out.residuals.reserve(operands.size());
  out.undefined.reserve(operands.size());
  for (std::size_t i = 0; i < operands.size(); ++i) {
    const model::Operand& op = operands[i];
    Value v;
    if (rows_[i]) {
      if (!table) table = model::evaluate_parameters(system);
      v.value = table->at(*rows_[i], configurations_[i]);
      if (!std::isfinite(v.value)) v.undefined = "the parameter row has no finite value";
    } else {
      const compile::CompiledSystem& cs = *compiled[configurations_[i]];
      if (const auto* f = std::get_if<model::FirstOrderOperand>(&op)) {
        v = first_order_value(cs, *f);
      } else if (const auto* r = std::get_if<model::RayOperand>(&op)) {
        v = ray_value(cs, *r);
      } else if (const auto* s = std::get_if<model::SpotRmsOperand>(&op)) {
        analysis::SpotOptions options;
        options.sampling = trace::HexapolarPupil{s->rings};
        const std::optional<std::uint16_t> wl =
            s->polychromatic ? std::nullopt
                             : std::optional(s->wavelength.value_or(cs.reference_wavelength()));
        const analysis::SpotDiagram spot =
            analysis::spot(cs, path_id(cs, s->path), s->field, wl, options);
        v.value = s->reference == model::SpotReference::Centroid ? spot.stats.rms_centroid
                                                                 : spot.stats.rms_chief;
        add_warnings(out.warnings, spot.warnings);
      } else {
        const auto& o = std::get<model::OpdRmsOperand>(op);
        analysis::OpdOptions options;
        options.grid = o.grid;
        const analysis::OpdMap map =
            analysis::opd_map(cs, path_id(cs, o.path), o.field,
                              o.wavelength.value_or(cs.reference_wavelength()), options);
        v.value = map.rms;
        add_warnings(out.warnings, map.warnings);
      }
    }
    // ADR 0030, point 1: r = sqrt(w) (v - t).
    const model::OperandCommon& common = common_of(op);
    out.values.push_back(v.value);
    out.residuals.push_back(std::sqrt(common.weight) * (v.value - common.target));
    out.undefined.push_back(std::move(v.undefined));
  }
  return out;
}

}  // namespace rtt::optim
