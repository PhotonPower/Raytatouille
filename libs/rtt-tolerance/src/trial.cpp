#include "rtt/tolerance/trial.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <functional>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <numbers>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <variant>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/diagnostics/codes.hpp"
#include "rtt/math/types.hpp"
#include "rtt/model/parameters.hpp"
#include "rtt/optim/merit.hpp"
#include "rtt/optim/overlay.hpp"
#include "rtt/optim/variables.hpp"

namespace rtt::tolerance {

namespace detail {

/// What one item does to a trial (ADR 0036, point 1).
struct ItemAction {
  enum class Kind : std::uint8_t {
    ParamAdd,    ///< (a) offset on a Param: radius, conic, thickness, param
    ParamScale,  ///< (a) curvature_relative: R' = R / (1 + delta)
    RowAdd,      ///< (a) parameter: the row in every configuration (K6)
    Rigid,       ///< (b) one component of the rigid motion of a node or surface
    Material,    ///< (c) Delta n_d or Delta nu_d of an element segment
  };
  Kind kind = Kind::ParamAdd;
  std::size_t index = 0;  ///< Param (for_each_param order) or row
  bool surface = false;   ///< Rigid: a surface (SurfaceId) or a node (name)
  std::string node;
  model::SurfaceId surface_id;
  int component = 0;  ///< Rigid: 0..2 translation x, y, z; 3..5 rotation x, y, z
  compile::MaterialTarget material;
  bool abbe = false;  ///< Material: Delta nu_d, else Delta n_d
  /// ParamAdd, ParamScale: the Param is the base radius of a surface (radius,
  /// curvature_relative, param at /shape/base/radius); its nominal value in mm.
  std::optional<double> radius_mm;
};

/// The target of a compensator: a Param or a row of the form value.
struct CompensatorTarget {
  bool row = false;
  std::size_t index = 0;  ///< Param (for_each_param order) or row
  std::string pointer;    ///< optim::Variable::pointer of the variable it becomes
};

/// A base radius of the nominal system in one configuration (ADR 0036, addendum #293).
struct NominalRadius {
  std::string pointer;               ///< .../shape/base/radius (edit form, ADR 0024)
  std::string surface;               ///< its SurfaceId
  std::optional<std::string> row;    ///< the row it is bound to, if any
  std::optional<std::size_t> index;  ///< that row's index in System::parameters
  std::vector<double> value;         ///< mm per configuration, finite and non-zero (model::Conic)
};

struct ContextData {
  model::System nominal;
  const material::MaterialLibrary* materials = nullptr;
  const coating::CoatingLibrary* coatings = nullptr;
  std::vector<model::ToleranceItem> items;
  std::vector<ItemAction> actions;
  std::map<std::string, math::Vec3, std::less<>> node_pivots;
  std::map<model::SurfaceId, math::Vec3> surface_pivots;
  std::vector<CompensatorTarget> compensators;
  model::Optimization compensator_merit;
  model::Optimization criteria;  ///< the criteria's operands, for the merit function
  TrialOptions options;
  std::optional<optim::NominalBasis> basis;
  TrialResult nominal_result;
  /// Per configuration the nominal vertex z of the reference of the object distance, empty
  /// without one (ADR 0038, point 5): the object of a trial stays at its nominal reference.
  std::vector<double> object_reference_z;
  /// Per configuration the nominal base radii (ADR 0036, addendum #293).
  std::vector<NominalRadius> radii;
};

struct ContextAccess {
  static const ContextData& data(const TrialContext& context) { return *context.data_; }
};

}  // namespace detail

namespace {

using detail::ContextData;
using detail::ItemAction;

model::Diagnostic diagnostic(diagnostics::DiagnosticCode code,
                             std::string location,
                             std::string message) {
  return {code.severity(), std::move(location), std::move(message), std::string(code.str())};
}

/// The JSON pointers (edit form, ADR 0024) of the nodes and surfaces of a system, by name.
struct Anchors {
  std::unordered_map<std::string, std::string> nodes;
  std::unordered_map<std::string, std::string> surfaces;
  std::unordered_map<std::string, const model::Element*> elements;

  explicit Anchors(const model::System& s) { assembly(s.root, "/root"); }

  void assembly(const model::Assembly& a, const std::string& loc) {
    for (std::size_t i = 0; i < a.children.size(); ++i) {
      const std::string child = loc + "/children/" + std::to_string(i);
      if (const auto* sub = std::get_if<model::Assembly>(&a.children[i].value)) {
        nodes.emplace(sub->name, child);
        assembly(*sub, child);
      } else {
        const auto& e = std::get<model::Element>(a.children[i].value);
        nodes.emplace(e.name, child);
        elements.emplace(e.name, &e);
        for (std::size_t k = 0; k < e.surfaces.size(); ++k) {
          surfaces.emplace(e.surfaces[k].id.str(), child + "/surfaces/" + std::to_string(k));
        }
      }
    }
  }

  /// The pointer of the anchor of `t` (node or surface).
  /// @throws std::invalid_argument if it does not exist
  [[nodiscard]] const std::string& of(const model::ToleranceTarget& t) const {
    const auto& map = t.anchor == model::ToleranceAnchor::Surface ? surfaces : nodes;
    const auto it = t.anchor == model::ToleranceAnchor::Parameter ? map.end() : map.find(t.name);
    if (it == map.end()) {
      throw std::invalid_argument("tolerance: the target '" + t.name + "' does not exist");
    }
    return it->second;
  }
};

/// The Params of a system by pointer, in the order of model::for_each_param.
std::unordered_map<std::string, std::size_t> param_indices(const model::System& s) {
  std::unordered_map<std::string, std::size_t> out;
  model::for_each_param(s, [&](std::string_view pointer, const model::Param& /*param*/) {
    out.emplace(std::string(pointer), out.size());
  });
  return out;
}

/// The Params of a system in the order of model::for_each_param, for writing.
std::vector<model::Param*> params_of(model::System& s) {
  std::vector<model::Param*> out;
  model::for_each_param(s,
                        [&](std::string_view /*pointer*/, model::Param& p) { out.push_back(&p); });
  return out;
}

std::size_t param_index(const std::unordered_map<std::string, std::size_t>& params,
                        const model::System& s,
                        const std::string& pointer,
                        const std::string& what) {
  const auto it = params.find(pointer);
  if (it == params.end()) {
    throw std::invalid_argument("tolerance: " + what + " has no Param at '" + pointer + "'");
  }
  bool bound = false;
  model::for_each_param(s, [&](std::string_view p, const model::Param& param) {
    if (p == pointer) bound = param.is_bound();
  });
  if (bound) {
    throw std::invalid_argument("tolerance: " + what + " at the bound Param '" + pointer +
                                "'; tolerate the row (ADR 0035, point 3)");
  }
  return it->second;
}

std::size_t row_of(const model::System& s, const std::string& name, const std::string& what) {
  const std::optional<std::size_t> row = model::find_parameter(s, name);
  if (!row) throw std::invalid_argument("tolerance: " + what + ": no row '" + name + "'");
  if (std::holds_alternative<model::ParameterExpression>(s.parameters[*row].form)) {
    throw std::invalid_argument("tolerance: " + what + " at the expression row '" + name + "'");
  }
  return *row;
}

/// The action of one item (ADR 0036, point 1). Material items are checked in check_materials().
ItemAction action_of(const model::ToleranceItem& item,
                     const model::System& s,
                     const Anchors& anchors,
                     const std::unordered_map<std::string, std::size_t>& params) {
  using K = model::ToleranceKind;
  const std::string what = "item '" + item.id + "'";
  ItemAction a;
  switch (item.kind) {
    case K::Radius:
    case K::CurvatureRelative:
    case K::Conic: {
      if (item.target.anchor != model::ToleranceAnchor::Surface) {
        throw std::invalid_argument("tolerance: " + what + " needs a surface");
      }
      const std::string field = item.kind == K::Conic ? "/shape/base/conic" : "/shape/base/radius";
      a.kind = item.kind == K::CurvatureRelative ? ItemAction::Kind::ParamScale
                                                 : ItemAction::Kind::ParamAdd;
      a.index = param_index(params, s, anchors.of(item.target) + field, what);
      if (item.kind != K::Conic) a.radius_mm = 0.0;  // the nominal value: prepare()
      return a;
    }
    case K::Thickness:
      a.index = param_index(params, s, anchors.of(item.target) + "/pose/position/2", what);
      return a;
    case K::Param:
      if (!item.target.param) throw std::invalid_argument("tolerance: " + what + " needs param");
      a.index = param_index(params, s, anchors.of(item.target) + *item.target.param, what);
      if (item.target.anchor == model::ToleranceAnchor::Surface &&
          *item.target.param == "/shape/base/radius") {
        a.radius_mm = 0.0;  // the nominal value: prepare()
      }
      return a;
    case K::Parameter:
      if (item.target.anchor != model::ToleranceAnchor::Parameter) {
        throw std::invalid_argument("tolerance: " + what + " needs a parameter row");
      }
      a.kind = ItemAction::Kind::RowAdd;
      a.index = row_of(s, item.target.name, what);
      return a;
    case K::DecenterX:
    case K::DecenterY:
    case K::ShiftZ:
    case K::TiltXDeg:
    case K::TiltYDeg:
    case K::TiltZDeg: {
      static_cast<void>(anchors.of(item.target));  // exists
      a.kind = ItemAction::Kind::Rigid;
      a.surface = item.target.anchor == model::ToleranceAnchor::Surface;
      if (a.surface) {
        a.surface_id = model::SurfaceId(item.target.name);
      } else {
        a.node = item.target.name;
      }
      a.component = item.kind == K::DecenterX   ? 0
                    : item.kind == K::DecenterY ? 1
                    : item.kind == K::ShiftZ    ? 2
                    : item.kind == K::TiltXDeg  ? 3
                    : item.kind == K::TiltYDeg  ? 4
                                                : 5;
      return a;
    }
    case K::Index:
    case K::Abbe: {
      if (item.target.anchor != model::ToleranceAnchor::Node ||
          !anchors.elements.contains(item.target.name)) {
        throw std::invalid_argument("tolerance: " + what + " needs an element");
      }
      a.kind = ItemAction::Kind::Material;
      a.material = compile::MaterialTarget{item.target.name, item.target.segment};
      a.abbe = item.kind == K::Abbe;
      return a;
    }
  }
  throw std::invalid_argument("tolerance: " + what + " has an unknown kind");
}

/// The glass of each material item gives n_d and, for abbe, nu_d at the lines d, F, C
/// (compile::fraunhofer_indices, ADR 0036, point 1c; #260): checked once here, before any trial,
/// as compile would report it in the first trial (tolerances.material_unsupported,
/// material.unknown, tolerances.abbe_on_const at the material of the element).
void check_materials(const ContextData& d, const Anchors& anchors) {
  std::vector<model::Diagnostic> errors;
  for (std::size_t i = 0; i < d.items.size(); ++i) {
    const ItemAction& a = d.actions[i];
    if (a.kind != ItemAction::Kind::Material) continue;
    const model::Element& e = *anchors.elements.at(a.material.element);
    std::string location = anchors.nodes.at(a.material.element) + "/material";
    const std::string* shorthand = e.material ? &*e.material : nullptr;
    if (e.crystal || (shorthand == nullptr && e.segment_materials.empty())) {
      errors.push_back(diagnostic("tolerances.material_unsupported", location,
                                  "element '" + e.name + "' has no isotropic material"));
      continue;
    }
    // The segment rules of compile (ADR 0035, point 2), here before any trial: the shorthand
    // takes no segment, a list needs one, a list of equal entries is one material (segment 0).
    const std::string what = "tolerance: item '" + d.items[i].id + "': element '" + e.name + "'";
    const std::vector<std::string>& list = e.segment_materials;
    std::string reference;
    if (list.empty()) {
      if (a.material.segment.has_value()) {
        throw std::invalid_argument(what + " has one material and takes no segment");
      }
      reference = *shorthand;
    } else {
      if (!a.material.segment.has_value()) {
        throw std::invalid_argument(what + " needs a segment");
      }
      const std::size_t k = *a.material.segment;
      const bool one_material =
          std::adjacent_find(list.begin(), list.end(), std::not_equal_to<>()) == list.end();
      if (one_material && k != 0) {
        throw std::invalid_argument(what +
                                    " has one material in all segments and takes only segment 0");
      }
      if (k >= list.size()) {
        throw std::invalid_argument(what + " needs a segment below " + std::to_string(list.size()));
      }
      reference = list[k];
      location += "/" + std::to_string(k);
    }
    std::shared_ptr<const material::Material> material;
    try {
      material = d.materials->resolve(reference);
    } catch (const material::UnknownMaterial& x) {
      errors.push_back(diagnostic("material.unknown", location, x.what()));
      continue;
    }
    const compile::FraunhoferIndices lines = compile::fraunhofer_indices(
        *material, d.nominal.environment.temperature_c, d.nominal.environment.pressure_atm);
    if (lines.status != compile::FraunhoferIndices::Status::Ok) {
      errors.push_back(diagnostic("tolerances.material_unsupported", location,
                                  "material '" + reference +
                                      "' has no real index at the lines d, F and C (item '" +
                                      d.items[i].id + "')"));
    } else if (a.abbe && !lines.abbe_number()) {
      errors.push_back(
          diagnostic("tolerances.abbe_on_const", location,
                     "material '" + reference + "' has n_F = n_C (item '" + d.items[i].id + "')"));
    }
  }
  if (!errors.empty()) throw compile::CompileError(std::move(errors));
}

/// The values applied to the nominal system (ADR 0036, point 1).
Applied applied(const ContextData& d, std::span<const double> values) {
  if (values.size() != d.items.size()) {
    throw std::invalid_argument("tolerance: " + std::to_string(values.size()) + " values for " +
                                std::to_string(d.items.size()) + " items");
  }
  for (const double v : values) {
    if (!std::isfinite(v)) throw std::invalid_argument("tolerance: a value is not finite");
  }
  Applied out{d.nominal, {}};
  const std::vector<model::Param*> params = params_of(out.system);
  for (std::size_t i = 0; i < values.size(); ++i) {
    const ItemAction& a = d.actions[i];
    const double v = values[i];
    switch (a.kind) {
      case ItemAction::Kind::ParamAdd:
        params[a.index]->value += v;
        break;
      case ItemAction::Kind::ParamScale:
        // C' = C (1 + delta), R' = R / (1 + delta) (ADR 0035, point 2), for any finite delta:
        // R' is infinite for delta = -1 and of the other sign for delta < -1. T < 1 does not
        // exclude that, since the sensitivity confirms at T' > T (ADR 0036, point 6); run()
        // makes such a trial Invalid (radius_invalid, ADR 0036, addendum #257).
        params[a.index]->value /= 1.0 + v;
        break;
      case ItemAction::Kind::RowAdd: {
        model::ParameterForm& form = out.system.parameters[a.index].form;
        if (auto* value = std::get_if<double>(&form)) {
          *value += v;
        } else if (auto* list = std::get_if<std::vector<double>>(&form)) {
          for (double& x : *list) x += v;  // every configuration (K6)
        }
        break;
      }
      case ItemAction::Kind::Rigid: {
        const auto pivot = [&]() -> std::optional<math::Vec3> {
          if (a.surface) {
            const auto it = d.surface_pivots.find(a.surface_id);
            if (it != d.surface_pivots.end()) return it->second;
          } else {
            const auto it = d.node_pivots.find(a.node);
            if (it != d.node_pivots.end()) return it->second;
          }
          return std::nullopt;
        }();
        compile::RigidMotion& m =
            a.surface ? out.perturbation.surfaces[a.surface_id] : out.perturbation.nodes[a.node];
        m.pivot_mm = pivot;
        if (a.component < 3) {
          m.translation_mm[a.component] += v;
        } else {
          m.rotation_rad[a.component - 3] += v * std::numbers::pi / 180.0;
        }
        break;
      }
      case ItemAction::Kind::Material: {
        compile::MaterialChange& c = out.perturbation.materials[a.material];
        (a.abbe ? c.delta_vd : c.delta_nd) += v;
        break;
      }
    }
  }
  out.perturbation.object_reference_z = d.object_reference_z;
  return out;
}

bool cancelled(const trace::RunControl* control) {
  return control != nullptr && control->cancel.has_value() && control->cancel->cancelled();
}

TrialResult cancelled_result(const ContextData& d) {
  TrialResult r;
  r.status = TrialStatus::Cancelled;
  r.reason = "cancelled";
  r.criteria.resize(d.nominal.tolerances.criteria.size());
  for (CriterionValue& c : r.criteria) c.reason = "cancelled";
  return r;
}

/// A trial that is Invalid before its criteria: every criterion without value, with `reason`.
TrialResult invalid_result(const ContextData& d, TrialResult r, const std::string& reason) {
  r.status = TrialStatus::Invalid;
  r.reason = reason;
  r.criteria.assign(d.criteria.operands.size(), CriterionValue{std::nullopt, reason, std::nullopt});
  return r;
}

/// Why the radii of `trial` make no valid trial, or none (ADR 0036, addendum #257): the Param of
/// a radius item that ends at 0, at +-infinity or with the other sign than its nominal value
/// describes another lens (a plane is the shape Plane, a change of sign another bending).
std::optional<std::string> radius_invalid(const ContextData& d, model::System& trial) {
  std::vector<model::Param*> params;
  for (std::size_t i = 0; i < d.actions.size(); ++i) {
    const ItemAction& a = d.actions[i];
    if (!a.radius_mm) continue;
    if (params.empty()) params = params_of(trial);
    const double r = params[a.index]->value;
    if (std::isfinite(r) && r != 0.0 && std::signbit(r) == std::signbit(*a.radius_mm)) continue;
    std::ostringstream why;
    why << "item '" << d.items[i].id << "': the radius " << r << " mm (nominal " << *a.radius_mm
        << " mm) is zero, infinite or of the other sign";
    return why.str();
  }
  return std::nullopt;
}

constexpr std::string_view kRadiusSuffix = "/shape/base/radius";

/// The base radius Params of `system` in tree order, with their pointers (edit form).
std::vector<std::pair<std::string, const model::Param*>> base_radii(const model::System& system) {
  std::vector<std::pair<std::string, const model::Param*>> out;
  model::for_each_param(system, [&](std::string_view pointer, const model::Param& p) {
    if (pointer.ends_with(kRadiusSuffix)) out.emplace_back(std::string(pointer), &p);
  });
  return out;
}

/// The nominal base radii with surface, row and value per configuration (ADR 0036, addendum
/// #293). A radius bound to a row takes the row's value in each configuration.
std::vector<detail::NominalRadius> nominal_radii(const model::System& nominal,
                                                 const Anchors& anchors) {
  std::unordered_map<std::string, std::string> surfaces;  // pointer -> SurfaceId
  for (const auto& [name, pointer] : anchors.surfaces) surfaces.emplace(pointer, name);
  const model::ParameterValues table = model::evaluate_parameters(nominal);
  std::vector<detail::NominalRadius> out;
  for (const auto& [pointer, param] : base_radii(nominal)) {
    detail::NominalRadius r;
    r.pointer = pointer;
    const std::string owner = pointer.substr(0, pointer.size() - kRadiusSuffix.size());
    if (const auto it = surfaces.find(owner); it != surfaces.end()) r.surface = it->second;
    if (param->param) {
      r.row = *param->param;
      r.index = model::find_parameter(nominal, *param->param);
    }
    for (std::size_t c = 0; c < table.configurations; ++c) {
      r.value.push_back(r.index ? table.at(*r.index, c) : param->value);
    }
    out.push_back(std::move(r));
  }
  return out;
}

/// Why a base radius of `system` makes no valid trial, or none (ADR 0036, addendum #293): in
/// any configuration, with the parameter table evaluated, a base radius that is 0, infinite or
/// of the other sign than its nominal value. This covers what radius_invalid does not see: rows
/// (also through expression rows and per configuration) and compensators. A row that is not
/// defined (NaN: an error in the table) is left to compile, which reports it with the criteria.
/// One walk over the Params and one evaluation of the table; validate does not run again.
std::optional<std::string> radius_sign_invalid(const ContextData& d,
                                               const model::System& system,
                                               std::string_view stage) {
  if (d.radii.empty()) return std::nullopt;
  const std::vector<std::pair<std::string, const model::Param*>> radii = base_radii(system);
  // Invariant: a trial has the tree of the nominal (values change, never the structure), so its
  // base radii come in the same order with the same pointers.
  const auto mismatch = [] {
    return std::logic_error("tolerance: the base radii of a trial differ from the nominal");
  };
  if (radii.size() != d.radii.size()) throw mismatch();
  for (std::size_t k = 0; k < radii.size(); ++k) {
    if (radii[k].first != d.radii[k].pointer) throw mismatch();
  }
  const bool bound =
      std::any_of(d.radii.begin(), d.radii.end(),
                  [](const detail::NominalRadius& n) { return n.index.has_value(); });
  const model::ParameterValues table =
      bound ? model::evaluate_parameters(system) : model::ParameterValues{};
  const std::size_t configurations = d.radii.front().value.size();
  for (std::size_t c = 0; c < configurations; ++c) {
    for (std::size_t k = 0; k < d.radii.size(); ++k) {
      const detail::NominalRadius& n = d.radii[k];
      const double r = n.index ? table.at(*n.index, c) : radii[k].second->value;
      if (std::isnan(r)) continue;  // an undefined row: compile reports the table error
      if (std::isfinite(r) && r != 0.0 && std::signbit(r) == std::signbit(n.value[c])) continue;
      std::ostringstream why;
      why << std::setprecision(std::numeric_limits<double>::max_digits10);  // 50.00000012 != 50
      why << stage << "surface '" << n.surface << "'";
      if (configurations > 1) why << ", configuration '" << system.configurations[c].name << "'";
      why << ": the radius " << r << " mm (nominal " << n.value[c] << " mm";
      if (n.row) why << ", row '" << *n.row << "'";
      why << ") is zero, infinite or of the other sign";
      return why.str();
    }
  }
  return std::nullopt;
}

/// The value of every criterion of `system` with the overlay (ADR 0036, point 2). One merit
/// function for all; if it throws (an analysis error of one criterion), each criterion alone, so
/// that the others keep their values.
std::vector<CriterionValue> criteria_values(const ContextData& d,
                                            const model::System& compensated,
                                            const optim::TrialOverlay& overlay,
                                            std::vector<model::Diagnostic>& warnings) {
  const std::size_t n = d.criteria.operands.size();
  std::vector<CriterionValue> out(n);
  if (n == 0) return out;
  const auto fill = [&](std::size_t k, const optim::MeritEvaluation& e, std::size_t j) {
    if (std::isfinite(e.values[j])) {
      out[k].value = e.values[j];
    } else {
      out[k].reason = e.undefined[j].empty() ? "not defined" : e.undefined[j];
    }
  };
  model::System s = compensated;
  s.optimization = d.criteria;
  std::optional<optim::MeritEvaluation> all;
  try {
    const optim::MeritFunction merit(s, *d.materials, d.coatings, &overlay);
    all = merit.evaluate(merit.start());
  } catch (const std::exception&) {
    // Any error of one criterion (OptimError, CompileError, an analysis error): below, criterion
    // by criterion, each with its own reason.
    all.reset();
  }
  if (all.has_value()) {
    const optim::MeritEvaluation& e = all.value();
    for (std::size_t k = 0; k < n; ++k) fill(k, e, k);
    warnings.insert(warnings.end(), e.warnings.begin(), e.warnings.end());
    return out;
  }
  for (std::size_t k = 0; k < n; ++k) {
    s.optimization = model::Optimization{{d.criteria.operands[k]}, {}};
    try {
      const optim::MeritFunction merit(s, *d.materials, d.coatings, &overlay);
      const optim::MeritEvaluation e = merit.evaluate(merit.start());
      fill(k, e, 0);
      warnings.insert(warnings.end(), e.warnings.begin(), e.warnings.end());
    } catch (const std::exception& x) {
      out[k].reason = x.what();
    }
  }
  return out;
}

/// One trial (ADR 0036, point 4). `nominal`: the criteria values of Phi_0, none while Phi_0 itself
/// is computed (then each value is its own nominal value).
TrialResult run(const ContextData& d,
                std::span<const double> values,
                const trace::RunControl* control,
                const std::vector<CriterionValue>* nominal) {
  Applied a = applied(d, values);
  if (cancelled(control)) return cancelled_result(d);
  if (!d.basis.has_value()) throw std::logic_error("tolerance: the context has no nominal basis");
  if (const std::optional<std::string> why = radius_invalid(d, a.system)) {
    return invalid_result(d, {}, *why);
  }
  if (const std::optional<std::string> why = radius_sign_invalid(d, a.system, "")) {
    return invalid_result(d, {}, *why);
  }
  const optim::TrialOverlay overlay{std::move(a.perturbation), d.basis.value()};
  TrialResult r;
  model::System trial = std::move(a.system);

  // Compensators (ADR 0036, point 3): all variable flags cleared, then the compensators variable
  // with the travel around the perturbed value, intersected with Param.min/max.
  if (!d.compensators.empty()) {
    const std::vector<model::Param*> params = params_of(trial);
    for (model::Param* p : params) p->variable = false;
    for (model::ParameterRow& row : trial.parameters) row.variable = false;
    bool any = false;
    for (std::size_t c = 0; c < d.compensators.size(); ++c) {
      const detail::CompensatorTarget& t = d.compensators[c];
      const model::Compensator& k = d.nominal.tolerances.compensators[c];
      double x = 0.0;
      std::optional<double>* min = nullptr;
      std::optional<double>* max = nullptr;
      bool* variable = nullptr;
      if (t.row) {
        model::ParameterRow& row = trial.parameters[t.index];
        x = std::get<double>(row.form);
        min = &row.min;
        max = &row.max;
        variable = &row.variable;
      } else {
        model::Param& p = *params[t.index];
        x = p.value;
        min = &p.min;
        max = &p.max;
        variable = &p.variable;
      }
      CompensatorValue v{x, x, false, false};
      // No inner start: x on or outside a bound of the Param stays (ADR 0036, point 3).
      if ((min->has_value() && !(x > min->value())) || (max->has_value() && !(x < max->value()))) {
        v.infeasible = true;
        r.warnings.push_back(diagnostic(
            "tolerances.compensator_infeasible", "/tolerances/compensators/" + std::to_string(c),
            "the perturbed value lies on or outside a bound of the compensator's "
            "target; it stays at that value"));
      } else {
        const double lo = x + k.travel_min;
        const double hi = x + k.travel_max;
        *min = min->has_value() ? std::max(min->value(), lo) : lo;
        *max = max->has_value() ? std::min(max->value(), hi) : hi;
        *variable = true;
        any = true;
      }
      r.compensators.push_back(v);
    }
    if (any) {
      trial.optimization = d.compensator_merit;
      trace::RunControl cancel_only;
      if (control != nullptr) cancel_only.cancel = control->cancel;
      optim::OptimResult result;
      try {
        result = optim::optimize(trial, *d.materials, d.coatings, d.options.compensators,
                                 cancel_only, overlay);
      } catch (const std::exception& x) {
        return invalid_result(d, std::move(r),
                              std::string("the compensator run is not possible: ") + x.what());
      }
      if (result.status == optim::LmStatus::Cancelled) return cancelled_result(d);
      r.lm_status = result.status;
      if (!result.history.empty()) r.mu = result.history.back().mu;
      r.lm_iterations = result.iterations;
      for (std::size_t c = 0; c < d.compensators.size(); ++c) {
        for (const optim::VariableValue& v : result.variables) {
          if (v.pointer != d.compensators[c].pointer) continue;
          r.compensators[c].end = v.end;
          r.compensators[c].at_bound = v.at_bound;
        }
      }
      r.warnings.insert(r.warnings.end(), result.diagnostics.begin(), result.diagnostics.end());
      if (result.status == optim::LmStatus::Failed) {
        return invalid_result(d, std::move(r),
                              "the compensator run failed (Levenberg-Marquardt status failed)");
      }
      trial = std::move(result.system);
      // Only the end of the run counts: LM-internal states do not reach the result.
      if (const std::optional<std::string> why =
              radius_sign_invalid(d, trial, "after the compensator run: ")) {
        return invalid_result(d, std::move(r), *why);
      }
    }
  }

  // Criteria with the nominal reference (ADR 0036, point 2).
  if (cancelled(control)) return cancelled_result(d);
  r.criteria = criteria_values(d, trial, overlay, r.warnings);
  const std::vector<model::Criterion>& criteria = d.nominal.tolerances.criteria;
  for (std::size_t k = 0; k < r.criteria.size(); ++k) {
    CriterionValue& c = r.criteria[k];
    const std::optional<double> phi0 = nominal != nullptr ? (*nominal)[k].value : c.value;
    c.satisfied = satisfies(criteria[k], c.value, phi0);
    if (!c.value && r.status == TrialStatus::Valid) {
      r.status = TrialStatus::Invalid;
      r.reason = "criterion '" + criteria[k].id + "': " + c.reason;
    }
  }
  return r;
}

}  // namespace

std::span<const model::ToleranceItem> TrialContext::items() const noexcept {
  return data_->items;
}

std::span<const model::Criterion> TrialContext::criteria() const noexcept {
  return data_->nominal.tolerances.criteria;
}

std::span<const model::Compensator> TrialContext::compensators() const noexcept {
  return data_->nominal.tolerances.compensators;
}

const std::optional<model::MonteCarloDefaults>& TrialContext::monte_carlo_defaults()
    const noexcept {
  return data_->nominal.tolerances.monte_carlo;
}

const TrialResult& TrialContext::nominal() const noexcept {
  return data_->nominal_result;
}

TrialContext prepare(const model::System& nominal,
                     const material::MaterialLibrary& materials,
                     const coating::CoatingLibrary* coatings,
                     std::span<const model::ToleranceItem> items,
                     const TrialOptions& options,
                     const trace::RunControl* control) {
  std::vector<model::Diagnostic> errors = model::validate(nominal);
  std::erase_if(errors,
                [](const model::Diagnostic& x) { return x.severity != model::Severity::Error; });
  if (!errors.empty()) throw compile::CompileError(std::move(errors));

  auto d = std::make_shared<ContextData>();
  d->nominal = nominal;
  d->materials = &materials;
  d->coatings = coatings;
  d->items.assign(items.begin(), items.end());
  d->options = options;
  const model::ToleranceSet& set = nominal.tolerances;

  const Anchors anchors(nominal);
  const std::unordered_map<std::string, std::size_t> params = param_indices(nominal);
  for (const model::ToleranceItem& item : d->items) {
    d->actions.push_back(action_of(item, nominal, anchors, params));
  }
  {
    std::vector<double> nominal_values;
    model::for_each_param(nominal, [&](std::string_view /*pointer*/, const model::Param& p) {
      nominal_values.push_back(p.value);
    });
    for (ItemAction& a : d->actions) {
      if (a.radius_mm) a.radius_mm = nominal_values[a.index];
    }
  }
  check_materials(*d, anchors);
  // ADR 0038, point 5: the nominal vertex z of the object reference, one per configuration.
  if (nominal.object.reference && !nominal.object.at_infinity) {
    const coating::CoatingLibrary none;
    for (std::size_t c = 0; c < model::configuration_count(nominal); ++c) {
      const compile::CompiledSystem cs =
          compile::compile(nominal, materials, coatings != nullptr ? *coatings : none, c);
      const std::optional<double> z = cs.object_reference_z();
      if (!z) throw std::logic_error("prepare: no nominal vertex z for object.reference");
      d->object_reference_z.push_back(*z);
    }
  }
  d->radii = nominal_radii(nominal, anchors);
  // Pivots of the targets (ADR 0035, point 2): one per target, for explicit and generated items.
  for (const model::ToleranceOverride& o : set.overrides) {
    const auto* p = std::get_if<model::TolerancePivotOverride>(&o);
    if (p == nullptr) continue;
    const math::Vec3 pivot(p->pivot[0], p->pivot[1], p->pivot[2]);
    if (p->target.anchor == model::ToleranceAnchor::Surface) {
      d->surface_pivots.emplace(model::SurfaceId(p->target.name), pivot);
    } else {
      d->node_pivots.emplace(p->target.name, pivot);
    }
  }
  // Compensators (ADR 0035, point 6): a Param or a row of the form value.
  for (const model::Compensator& c : set.compensators) {
    detail::CompensatorTarget t;
    if (c.target.anchor == model::ToleranceAnchor::Parameter) {
      t.row = true;
      t.index = row_of(nominal, c.target.name, "compensator");
      if (!std::holds_alternative<double>(nominal.parameters[t.index].form)) {
        throw std::invalid_argument(
            "tolerance: compensator at a row that is not of the form value");
      }
      t.pointer = "/parameters/" + std::to_string(t.index) + "/value";
    } else {
      if (!c.target.param) throw std::invalid_argument("tolerance: compensator without param");
      const std::string pointer = anchors.of(c.target) + *c.target.param;
      t.index = param_index(params, nominal, pointer, "compensator");
      t.pointer = pointer + "/value";
    }
    d->compensators.push_back(std::move(t));
  }
  if (!d->compensators.empty()) {
    d->compensator_merit = set.compensator_merit ? *set.compensator_merit : nominal.optimization;
  }
  for (const model::Criterion& c : set.criteria) d->criteria.operands.push_back(c.operand);

  const optim::NominalBasis& basis =
      d->basis.emplace(optim::nominal_basis(nominal, materials, coatings));
  // Start checks with the nominal basis (ADR 0030, addendum #256): an operand the nominal system
  // does not support is an input error of the run, before any trial.
  {
    const optim::TrialOverlay overlay{{}, basis};
    model::System s = nominal;
    s.optimization = d->criteria;
    if (!s.optimization.empty()) {
      static_cast<void>(optim::MeritFunction(s, materials, coatings, &overlay));
    }
    if (!d->compensators.empty()) {
      s.optimization = d->compensator_merit;
      static_cast<void>(optim::MeritFunction(s, materials, coatings, &overlay));
    }
  }
  // Phi_0 (ADR 0036, points 2 and 6): the trial at zero, compensators included.
  const std::vector<double> zeros(d->items.size(), 0.0);
  d->nominal_result = cancelled(control) ? cancelled_result(*d) : run(*d, zeros, control, nullptr);
  return TrialContext(std::move(d));
}

TrialResult evaluate(const TrialContext& context,
                     std::span<const double> values,
                     const trace::RunControl* control) {
  const ContextData& d = detail::ContextAccess::data(context);
  if (cancelled(control)) {
    static_cast<void>(applied(d, values));  // input errors first
    return cancelled_result(d);
  }
  return run(d, values, control, &d.nominal_result.criteria);
}

std::optional<bool> satisfies(const model::Criterion& criterion,
                              std::optional<double> value,
                              std::optional<double> nominal_value) {
  if (!value) return std::nullopt;
  const double v = *value;
  if (criterion.max_change && !nominal_value) return std::nullopt;
  bool ok = true;
  if (criterion.min) ok = ok && *criterion.min <= v;
  if (criterion.max) ok = ok && v <= *criterion.max;
  if (criterion.max_change && nominal_value) {
    ok = ok && std::abs(v - *nominal_value) <= *criterion.max_change;
  }
  return ok;
}

Applied apply(const TrialContext& context, std::span<const double> values) {
  return applied(detail::ContextAccess::data(context), values);
}

}  // namespace rtt::tolerance
