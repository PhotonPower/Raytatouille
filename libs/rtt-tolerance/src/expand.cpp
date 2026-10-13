#include "rtt/tolerance/expand.hpp"

#include <cstddef>
#include <optional>
#include <utility>
#include <variant>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/model/validate.hpp"
#include "rtt/tolerance/classes.hpp"

namespace rtt::tolerance {

ToleranceItems expand(const model::System& system) {
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
  // and the distribution normal.
  for (std::size_t i = 0; i < out.classes.items.size(); ++i) {
    const model::ClassTolerance& c = out.classes.items[i];
    if (!c.generated()) continue;
    model::ToleranceItem item;
    item.id = c.id;
    item.target = c.target;
    item.kind = c.kind;
    if (const std::optional<std::size_t> o = c.value_override) {
      const auto& v = std::get<model::ToleranceValueOverride>(set.overrides[*o]);
      item.limit = v.limit;
      item.distribution = v.distribution;
    } else {
      item.limit = class_limit(c.tolerance_class, c.kind, c.level);
      item.distribution = model::ToleranceDistribution::Normal;
    }
    out.items.push_back(std::move(item));
    out.origins.push_back({ItemOrigin::Source::Class, i});
  }
  // Then the explicit items in file order.
  for (std::size_t i = 0; i < set.items.size(); ++i) {
    out.items.push_back(set.items[i]);
    out.origins.push_back({ItemOrigin::Source::Explicit, i});
  }
  return out;
}

}  // namespace rtt::tolerance
