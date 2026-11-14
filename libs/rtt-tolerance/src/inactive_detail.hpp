#pragma once

/// @file inactive_detail.hpp
/// Targets of the tolerance set in inactive subtrees (ADR 0049, point 5): shared by expand() and
/// prepare(). One rule with rtt-optim (collect_variables): a target is inactive if the JSON
/// pointer it acts on lies in an inactive node (model::inactive_nodes, in_inactive_node).

#include <cstddef>
#include <optional>
#include <string>
#include <variant>
#include <vector>

#include "rtt/model/parameters.hpp"
#include "rtt/model/system.hpp"
#include "rtt/model/tolerance.hpp"

namespace rtt::tolerance::detail {

namespace inactive {

/// The JSON pointer of the node or surface `target` names below `a` (at `location`), or none.
inline std::optional<std::string> anchor(const model::Assembly& a,
                                         const std::string& location,
                                         const model::ToleranceTarget& target) {
  for (std::size_t i = 0; i < a.children.size(); ++i) {
    std::string child = location + "/children/" + std::to_string(i);
    if (const auto* sub = std::get_if<model::Assembly>(&a.children[i].value)) {
      if (target.anchor == model::ToleranceAnchor::Node && sub->name == target.name) return child;
      if (std::optional<std::string> found = anchor(*sub, child, target)) return found;
      continue;
    }
    const auto& e = std::get<model::Element>(a.children[i].value);
    if (target.anchor == model::ToleranceAnchor::Node && e.name == target.name) return child;
    if (target.anchor == model::ToleranceAnchor::Surface) {
      for (std::size_t k = 0; k < e.surfaces.size(); ++k) {
        if (e.surfaces[k].id.str() == target.name) {
          return child + "/surfaces/" + std::to_string(k);
        }
      }
    }
  }
  return std::nullopt;
}

}  // namespace inactive

/// True if `target` acts in an inactive node of `system` (active false at the node, its
/// element or an ancestor): the pointer of its node or surface, followed by `param` if it has
/// one (a param may point into a child of an assembly). A row of the parameter table is never
/// inactive; an unknown target is not inactive (validate reports it).
inline bool target_inactive(const model::System& system, const model::ToleranceTarget& target) {
  if (target.anchor == model::ToleranceAnchor::Parameter) return false;
  const std::vector<std::string> nodes = model::inactive_nodes(system);
  if (nodes.empty()) return false;
  const std::optional<std::string> at = inactive::anchor(system.root, "/root", target);
  if (!at) return false;
  return model::in_inactive_node(nodes, *at + target.param.value_or(std::string()));
}

}  // namespace rtt::tolerance::detail
