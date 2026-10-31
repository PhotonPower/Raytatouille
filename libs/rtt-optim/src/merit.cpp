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

#include "generators_detail.hpp"
#include "overlay_basis.hpp"
#include "parabasal.hpp"
#include "rtt/analysis/opd.hpp"
#include "rtt/analysis/spot.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/diagnostics/codes.hpp"
#include "rtt/math/types.hpp"
#include "rtt/model/optimization.hpp"
#include "rtt/model/parameters.hpp"
#include "rtt/optim/generators.hpp"
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

std::string generator_pointer(std::size_t g) {
  return "/optimization/generators/" + std::to_string(g);
}

/// What both generator types have, for the start check.
struct GeneratorCommon {
  const std::string* path = nullptr;
  const std::optional<std::string>* configuration = nullptr;
  const std::optional<std::vector<std::uint16_t>>* fields = nullptr;
  const std::optional<std::vector<std::uint16_t>>* wavelengths = nullptr;
};

GeneratorCommon generator_common(const model::Generator& g) {
  if (const auto* s = std::get_if<model::SpotGenerator>(&g)) {
    return {&s->path, &s->configuration, &s->fields, &s->wavelengths};
  }
  const auto& w = std::get<model::WavefrontGenerator>(g);
  return {&w.path, &w.configuration, &w.fields, &w.wavelengths};
}

/// Sum of the model weights of the chosen entries (none: all); a generator needs a positive sum
/// to normalise them (ADR 0030, addendum #168).
template <typename Entries>
double chosen_weight(const std::optional<std::vector<std::uint16_t>>& list,
                     const Entries& entries) {
  double total = 0.0;
  if (list) {
    for (const std::uint16_t i : *list) total += entries[i].weight;
  } else {
    for (const auto& e : entries) total += e.weight;
  }
  return total;
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

/// compile() of one configuration, with the overlay's perturbation if there is one (#254).
compile::CompiledSystem compile_one(const model::System& system,
                                    const material::MaterialLibrary& materials,
                                    const coating::CoatingLibrary* coatings,
                                    std::size_t configuration,
                                    const TrialOverlay* overlay) {
  if (overlay != nullptr) {
    const coating::CoatingLibrary none;
    return compile::compile(system, materials, coatings != nullptr ? *coatings : none,
                            configuration, overlay->perturbation);
  }
  return coatings != nullptr ? compile::compile(system, materials, *coatings, configuration)
                             : compile::compile(system, materials, configuration);
}

/// The nominal data of `path` (by name) in configuration `k` of the overlay's basis.
/// @throws std::invalid_argument if the basis has no such path (another system)
const detail::PathData& nominal_path(const TrialOverlay& overlay,
                                     std::size_t k,
                                     const compile::CompiledSystem& cs,
                                     compile::PathId path) {
  const detail::BasisData& basis = detail::BasisAccess::data(overlay.nominal);
  const std::vector<detail::PathData>& paths = basis.configurations.at(k).paths;
  if (path.index >= paths.size() || paths[path.index].name != cs.path(path).name) {
    throw std::invalid_argument("optimize: the nominal basis has no path '" + cs.path(path).name +
                                "' (it belongs to another system)");
  }
  return paths[path.index];
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
Value ray_value(const compile::CompiledSystem& cs,
                const model::RayOperand& op,
                const analysis::NominalReference* nominal) {
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
  const std::span<const std::uint16_t> fields(&field, 1);
  // With an overlay aimed with the nominal aiming data (#255; ADR 0030, addendum #256).
  trace::RayBatch rays =
      nominal != nullptr
          ? trace::make_rays(cs, path, fields, wl, trace::SinglePupilPoint{op.px, op.py},
                             trace::Aiming::Real, nominal->aim(wl))
          : trace::make_rays(cs, path, fields, wl, trace::SinglePupilPoint{op.px, op.py});
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
                             const coating::CoatingLibrary* coatings,
                             const TrialOverlay* overlay)
    : system_(&system),
      materials_(&materials),
      coatings_(coatings),
      overlay_(overlay),
      variables_(collect_variables(system)) {
  // The model first, also when no operand needs compile (only param_value operands).
  std::vector<model::Diagnostic> errors = model::validate(system);
  std::erase_if(errors,
                [](const model::Diagnostic& d) { return d.severity != model::Severity::Error; });
  if (!errors.empty()) throw compile::CompileError(std::move(errors));
  if (overlay != nullptr &&
      overlay->nominal.configurations() != model::configuration_count(system)) {
    throw std::invalid_argument(
        "optimize: the nominal basis has another number of configurations than the system");
  }
  if (overlay != nullptr &&
      detail::BasisAccess::data(overlay->nominal).wavelengths != system.wavelengths.size()) {
    throw std::invalid_argument(
        "optimize: the nominal basis has another number of wavelengths than the system");
  }

  const model::Optimization& merit = system.optimization;
  std::vector<model::Diagnostic> unsupported;
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
  // Generators (ADR 0030, point 4; addendum #168): configuration and the weights of the choice;
  // validate has checked the indices (merit.index_out_of_range).
  for (std::size_t g = 0; g < merit.generators.size(); ++g) {
    const GeneratorCommon gc = generator_common(merit.generators[g]);
    std::size_t k = 0;
    if (*gc.configuration) {
      const std::optional<std::size_t> found =
          model::find_configuration(system, **gc.configuration);
      if (!found) throw std::invalid_argument("optimize: unknown configuration");
      k = *found;
    }
    generator_configurations_.push_back(k);
    compiled_.push_back(k);
    if (!(chosen_weight(*gc.fields, system.fields.points) > 0.0)) {
      unsupported.push_back(diagnostic("merit.operand_unsupported", generator_pointer(g),
                                       "the weights of the chosen fields sum to 0"));
    }
    if (!(chosen_weight(*gc.wavelengths, system.wavelengths) > 0.0)) {
      unsupported.push_back(diagnostic("merit.operand_unsupported", generator_pointer(g),
                                       "the weights of the chosen wavelengths sum to 0"));
    }
  }
  generator_sizes_.assign(merit.generators.size(), 0);
  std::sort(compiled_.begin(), compiled_.end());
  compiled_.erase(std::unique(compiled_.begin(), compiled_.end()), compiled_.end());

  // Paths that the analyses reject (ADR 0030, point 3): first_order at the start, which every
  // operand other than param_value needs (paraxial values, aiming). With an overlay that of the
  // nominal system, kept in its basis (addendum #256).
  for (const std::size_t k : compiled_) {
    const compile::CompiledSystem cs = compile_one(system, materials, coatings, k, overlay);
    for (std::size_t i = 0; i < merit.operands.size(); ++i) {
      if (rows_[i] || configurations_[i] != k) continue;
      const model::Operand& op = merit.operands[i];
      const std::uint16_t wl = wavelength_of(op).value_or(cs.reference_wavelength());
      if (overlay != nullptr) {
        const detail::PathData& nominal = nominal_path(*overlay, k, cs, path_id(cs, *path_of(op)));
        const detail::WavelengthData& w = nominal.wavelengths.at(wl);
        if (!w.unsupported.empty()) {
          unsupported.push_back(diagnostic("merit.operand_unsupported", operand_pointer(i),
                                           "the path is not supported: " + w.unsupported));
          continue;
        }
        const auto* f = std::get_if<model::FirstOrderOperand>(&op);
        if (f != nullptr && (f->quantity == model::FirstOrderQuantity::Efl ||
                             f->quantity == model::FirstOrderQuantity::Bfl)) {
          if (!w.infinite_unsupported.empty()) {
            unsupported.push_back(diagnostic(
                "merit.operand_unsupported", operand_pointer(i),
                "efl and bfl of a trial need an axis-parallel bundle aimed at the stop, which the "
                "nominal system cannot aim (addendum #256): " +
                    w.infinite_unsupported));
          } else if (!system.object.at_infinity) {
            variants_.push_back(k);
          }
        }
        continue;
      }
      try {
        static_cast<void>(paraxial::first_order(cs, path_id(cs, *path_of(op)), wl));
      } catch (const paraxial::ParaxialError& e) {
        unsupported.push_back(diagnostic("merit.operand_unsupported", operand_pointer(i),
                                         std::string("the path is not supported: ") + e.what()));
      }
    }
    for (std::size_t g = 0; g < merit.generators.size(); ++g) {
      if (generator_configurations_[g] != k) continue;
      const model::Generator& gen = merit.generators[g];
      generator_sizes_[g] = generator_size(gen, cs);
      const GeneratorCommon gc = generator_common(gen);
      std::vector<std::uint16_t> wavelengths;
      if (*gc.wavelengths) {
        wavelengths = **gc.wavelengths;
      } else {
        for (std::size_t l = 0; l < cs.wavelengths_um().size(); ++l) {
          wavelengths.push_back(static_cast<std::uint16_t>(l));
        }
      }
      for (const std::uint16_t wl : wavelengths) {
        if (overlay != nullptr) {
          const detail::PathData& nominal = nominal_path(*overlay, k, cs, path_id(cs, *gc.path));
          const std::string& why = nominal.wavelengths.at(wl).unsupported;
          if (why.empty()) continue;
          unsupported.push_back(diagnostic("merit.operand_unsupported", generator_pointer(g),
                                           "the path is not supported: " + why));
          break;
        }
        try {
          static_cast<void>(paraxial::first_order(cs, path_id(cs, *gc.path), wl));
        } catch (const paraxial::ParaxialError& e) {
          unsupported.push_back(diagnostic("merit.operand_unsupported", generator_pointer(g),
                                           std::string("the path is not supported: ") + e.what()));
          break;
        }
      }
    }
  }
  if (!unsupported.empty()) throw OptimError(std::move(unsupported));
  std::sort(variants_.begin(), variants_.end());
  variants_.erase(std::unique(variants_.begin(), variants_.end()), variants_.end());
  size_ = merit.operands.size();
  for (const std::size_t n : generator_sizes_) size_ += n;
}

std::vector<double> MeritFunction::start() const {
  std::vector<double> out;
  out.reserve(variables_.size());
  for (const Variable& v : variables_) out.push_back(v.start);
  return out;
}

std::size_t MeritFunction::size() const noexcept {
  return size_;
}

MeritEvaluation MeritFunction::evaluate(std::span<const double> values) const {
  const model::System system = with_values(*system_, variables_, values);
  // Every used configuration once, in ascending order (ADR 0030, point 13).
  std::vector<std::optional<compile::CompiledSystem>> compiled(model::configuration_count(system));
  for (const std::size_t k : compiled_) {
    compiled[k].emplace(compile_one(system, *materials_, coatings_, k, overlay_));
  }
  // The variant with the object at infinity for efl and bfl with a finite object, compiled with
  // the same perturbation (addendum #256).
  std::vector<std::optional<compile::CompiledSystem>> variants(compiled.size());
  if (!variants_.empty()) {
    const model::System variant = detail::variant_of(system);
    // The variant has no object.reference, so the nominal reference z of the overlay does not
    // apply to it (ADR 0038, point 5); the rest of the perturbation does.
    std::optional<TrialOverlay> variant_overlay;
    if (overlay_ != nullptr) {
      variant_overlay.emplace(*overlay_);
      variant_overlay->perturbation.object_reference_z.clear();
    }
    const TrialOverlay* const vo = variant_overlay ? &*variant_overlay : nullptr;
    for (const std::size_t k : variants_) {
      variants[k].emplace(compile_one(variant, *materials_, coatings_, k, vo));
    }
  }
  std::optional<model::ParameterValues> table;

  const std::vector<model::Operand>& operands = system_->optimization.operands;
  MeritEvaluation out;
  out.values.reserve(operands.size());
  out.residuals.reserve(size_);
  out.undefined.reserve(operands.size());
  out.operand_arrived.reserve(operands.size());
  for (std::size_t i = 0; i < operands.size(); ++i) {
    const model::Operand& op = operands[i];
    Value v;
    std::vector<std::uint8_t> arrived;  // spot_rms and opd_rms only (addendum #299)
    if (const std::optional<std::size_t> row = rows_[i]; row.has_value()) {
      if (!table) table.emplace(model::evaluate_parameters(system));
      v.value = table.value().at(row.value(), configurations_[i]);
      if (!std::isfinite(v.value)) v.undefined = "the parameter row has no finite value";
    } else {
      // Compiled in the loop above: configurations_[i] is in compiled_ for this operand.
      const std::optional<compile::CompiledSystem>& slot = compiled[configurations_[i]];
      if (!slot) throw std::logic_error("optimize: configuration of an operand not compiled");
      const compile::CompiledSystem& cs = *slot;
      const std::size_t k = configurations_[i];
      const compile::PathId path = path_id(cs, *path_of(op));
      const detail::PathData* nominal =
          overlay_ != nullptr ? &nominal_path(*overlay_, k, cs, path) : nullptr;
      if (const auto* f = std::get_if<model::FirstOrderOperand>(&op)) {
        if (nominal != nullptr) {
          const std::uint16_t wl = f->wavelength.value_or(cs.reference_wavelength());
          const std::optional<compile::CompiledSystem>& inf = variants[k];
          const detail::ParabasalValue p = detail::parabasal_value(
              *f, cs, inf ? &inf.value() : nullptr, nominal->wavelengths.at(wl), nominal->reference,
              system.object.at_infinity, path, wl, overlay_->parabasal_step);
          v.value = p.value;
          v.undefined = p.undefined;
        } else {
          v = first_order_value(cs, *f);
        }
      } else if (const auto* r = std::get_if<model::RayOperand>(&op)) {
        v = ray_value(cs, *r, nominal != nullptr ? &nominal->reference : nullptr);
      } else if (const auto* s = std::get_if<model::SpotRmsOperand>(&op)) {
        analysis::SpotOptions options;
        options.sampling = trace::HexapolarPupil{s->rings};
        if (nominal != nullptr) options.reference = nominal->reference;
        const std::optional<std::uint16_t> wl =
            s->polychromatic ? std::nullopt
                             : std::optional(s->wavelength.value_or(cs.reference_wavelength()));
        analysis::SpotDiagram spot =
            analysis::spot(cs, path_id(cs, s->path), s->field, wl, options);
        v.value = s->reference == model::SpotReference::Centroid ? spot.stats.rms_centroid
                                                                 : spot.stats.rms_chief;
        add_warnings(out.warnings, spot.warnings);
        arrived = std::move(spot.arrived);
      } else {
        const auto& o = std::get<model::OpdRmsOperand>(op);
        analysis::OpdOptions options;
        options.grid = o.grid;
        if (nominal != nullptr) options.reference = nominal->reference;
        const analysis::OpdMap map =
            analysis::opd_map(cs, path_id(cs, o.path), o.field,
                              o.wavelength.value_or(cs.reference_wavelength()), options);
        v.value = map.rms;
        add_warnings(out.warnings, map.warnings);
        arrived.reserve(map.points.size());
        for (const analysis::OpdPoint& q : map.points) {
          arrived.push_back(q.status == trace::RayStatus::Alive ? std::uint8_t{1}
                                                                : std::uint8_t{0});
        }
      }
    }
    // ADR 0030, point 1: r = sqrt(w) (v - t).
    const model::OperandCommon& common = common_of(op);
    out.values.push_back(v.value);
    out.residuals.push_back(std::sqrt(common.weight) * (v.value - common.target));
    out.undefined.push_back(std::move(v.undefined));
    out.operand_arrived.push_back(std::move(arrived));
  }

  // Generators after the operands, each in its configuration (ADR 0030, addendum #168).
  const std::vector<model::Generator>& generators = system_->optimization.generators;
  out.residuals.resize(size_);
  out.generators.resize(generators.size());
  std::size_t offset = operands.size();
  for (std::size_t g = 0; g < generators.size(); ++g) {
    const std::size_t k = generator_configurations_[g];
    const std::optional<compile::CompiledSystem>& slot = compiled[k];
    if (!slot) throw std::logic_error("optimize: configuration of a generator not compiled");
    const analysis::NominalReference* nominal = nullptr;
    if (overlay_ != nullptr) {
      const GeneratorCommon gc = generator_common(generators[g]);
      nominal = &nominal_path(*overlay_, k, *slot, path_id(*slot, *gc.path)).reference;
    }
    detail::generator_residuals(
        generators[g], *slot, nominal,
        std::span<double>(out.residuals).subspan(offset, generator_sizes_[g]), out.generators[g]);
    offset += generator_sizes_[g];
  }
  return out;
}

}  // namespace rtt::optim
