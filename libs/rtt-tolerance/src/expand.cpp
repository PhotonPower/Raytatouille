#include "rtt/tolerance/expand.hpp"

#include <array>
#include <cstddef>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <variant>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/compile/perturbation.hpp"
#include "rtt/model/validate.hpp"
#include "rtt/tolerance/classes.hpp"

namespace rtt::tolerance {

namespace {

/// An element and its JSON pointer in the system file.
struct ElementAt {
  const model::Element* element = nullptr;
  std::string location;
};

void collect(const model::Assembly& a,
             const std::string& location,
             std::unordered_map<std::string, ElementAt>& out) {
  for (std::size_t i = 0; i < a.children.size(); ++i) {
    const std::string child = location + "/children/" + std::to_string(i);
    if (const auto* sub = std::get_if<model::Assembly>(&a.children[i].value)) {
      collect(*sub, child, out);
    } else {
      const auto& e = std::get<model::Element>(a.children[i].value);
      out.emplace(e.name, ElementAt{&e, child});
    }
  }
}

/// The limits p * nu_d of the generated abbe items (ADR 0035, addendum #260). Reports the glass
/// of an item that cannot give nu_d into `errors`; skips abbe at a glass without dispersion.
class AbbeLimits {
 public:
  AbbeLimits(const model::System& system, const material::MaterialLibrary& materials)
      : system_(system), materials_(materials) {
    collect(system.root, "/root", elements_);
  }

  /// nu_d of the glass of `c`, or none (reported, or skipped with ConstAbbe in `c`).
  std::optional<double> nu_d(model::ClassTolerance& c) {
    const auto it = elements_.find(c.target.name);
    if (it == elements_.end()) return std::nullopt;  // validate reports unknown targets
    const model::Element& e = *it->second.element;
    std::string reference;
    std::string location = it->second.location + "/material";
    if (!e.segment_materials.empty()) {
      const std::size_t segment = c.target.segment.value_or(0);
      reference = e.segment_materials.at(segment);
      location += "/" + std::to_string(segment);
    } else if (e.material) {
      reference = *e.material;
    } else {
      return std::nullopt;  // the class skips abbe there (NoMaterial)
    }
    std::shared_ptr<const material::Material> material;
    try {
      material = materials_.resolve(reference);
    } catch (const material::UnknownMaterial& x) {
      report("material.unknown", location, x.what());
      return std::nullopt;
    }
    const compile::FraunhoferIndices lines = compile::fraunhofer_indices(
        *material, system_.environment.temperature_c, system_.environment.pressure_atm);
    if (lines.status != compile::FraunhoferIndices::Status::Ok) {
      constexpr std::array<const char*, 3> kNames{"d", "F", "C"};
      const std::string why = lines.status == compile::FraunhoferIndices::Status::Complex
                                  ? "has a complex index at the line "
                                  : "is not defined at the line ";
      report("tolerances.material_unsupported", location,
             "material '" + reference + "' " + why + kNames.at(lines.line) +
                 ": no class limit for abbe (ADR 0035, addendum #260)");
      return std::nullopt;
    }
    const std::optional<double> vd = lines.abbe_number();
    if (!vd) c.skip = model::ClassSkip::ConstAbbe;  // n_F = n_C: no nu_d
    return vd;
  }

  std::vector<model::Diagnostic> errors;

 private:
  void report(const std::string& code, const std::string& location, const std::string& message) {
    model::Diagnostic d;
    d.severity = model::Severity::Error;
    d.location = location;
    d.message = message;
    d.code = code;
    errors.push_back(std::move(d));
  }

  const model::System& system_;
  const material::MaterialLibrary& materials_;
  std::unordered_map<std::string, ElementAt> elements_;
};

}  // namespace

ToleranceItems expand(const model::System& system, const material::MaterialLibrary& materials) {
  // The input is checked at the API boundary, as optimize does (rtt-optim, merit.cpp).
  std::vector<model::Diagnostic> errors = model::validate(system);
  std::erase_if(errors,
                [](const model::Diagnostic& d) { return d.severity != model::Severity::Error; });
  if (!errors.empty()) throw compile::CompileError(std::move(errors));

  ToleranceItems out;
  out.classes = model::class_tolerances(system);
  const model::ToleranceSet& set = system.tolerances;
  out.items.reserve(out.classes.items.size() + set.items.size());
  out.origins.reserve(out.classes.items.size() + set.items.size());

  // Generated items in tree order (ADR 0035, point 5): the single value, else the class table
  // and the distribution normal; abbe relative to nu_d of the glass (addendum #260).
  AbbeLimits abbe(system, materials);
  for (std::size_t i = 0; i < out.classes.items.size(); ++i) {
    model::ClassTolerance& c = out.classes.items[i];
    if (!c.generated()) continue;
    model::ToleranceItem item;
    item.id = c.id;
    item.target = c.target;
    item.kind = c.kind;
    if (const std::optional<std::size_t> o = c.value_override) {
      const auto& v = std::get<model::ToleranceValueOverride>(set.overrides[*o]);
      item.limit = v.limit;
      item.distribution = v.distribution;
    } else if (c.kind == model::ToleranceKind::Abbe) {
      const std::optional<double> vd = abbe.nu_d(c);
      if (!vd) continue;  // reported, or skipped without dispersion
      item.limit = abbe_share(c.tolerance_class) * *vd;
      item.distribution = model::ToleranceDistribution::Normal;
    } else {
      item.limit = class_limit(c.tolerance_class, c.kind, c.level);
      item.distribution = model::ToleranceDistribution::Normal;
    }
    out.items.push_back(std::move(item));
    out.origins.push_back({ItemOrigin::Source::Class, i});
  }
  if (!abbe.errors.empty()) throw compile::CompileError(std::move(abbe.errors));
  // Then the explicit items in file order.
  for (std::size_t i = 0; i < set.items.size(); ++i) {
    out.items.push_back(set.items[i]);
    out.origins.push_back({ItemOrigin::Source::Explicit, i});
  }
  return out;
}

}  // namespace rtt::tolerance
