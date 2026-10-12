#include "tolerance_targets.hpp"

#include <charconv>
#include <cmath>
#include <system_error>
#include <utility>
#include <variant>

#include "rtt/model/parameters.hpp"

namespace rtt::model::detail {

TargetIndex::TargetIndex(const System& system) : root_name_(system.root.name) {
  std::vector<std::string> ancestors;
  walk(system.root, "/root", ancestors);
  for (const ParameterRow& r : system.parameters) rows_.emplace(r.name, &r);
  for_each_param(system, [this](std::string_view pointer, const Param& p) {
    params_.emplace(std::string(pointer), &p);
  });
}

void TargetIndex::walk(const Assembly& a,
                       const std::string& pointer,
                       std::vector<std::string>& ancestors) {
  for (std::size_t i = 0; i < a.children.size(); ++i) {
    const std::string child = pointer + "/children/" + std::to_string(i);
    if (const auto* sub = std::get_if<Assembly>(&a.children[i].value)) {
      nodes_.emplace(sub->name, NodeEntry{child, sub, nullptr, ancestors});
      ancestors.push_back(sub->name);
      walk(*sub, child, ancestors);
      ancestors.pop_back();
    } else {
      const Element& e = std::get<Element>(a.children[i].value);
      nodes_.emplace(e.name, NodeEntry{child, nullptr, &e, ancestors});
      for (std::size_t j = 0; j < e.surfaces.size(); ++j) {
        surfaces_.emplace(
            e.surfaces[j].id.str(),
            SurfaceEntry{child + "/surfaces/" + std::to_string(j), &e.surfaces[j], &e, j});
      }
    }
  }
}

const NodeEntry* TargetIndex::node(const std::string& name) const {
  const auto it = nodes_.find(name);
  return it == nodes_.end() ? nullptr : &it->second;
}

const SurfaceEntry* TargetIndex::surface(const std::string& id) const {
  const auto it = surfaces_.find(id);
  return it == surfaces_.end() ? nullptr : &it->second;
}

const ParameterRow* TargetIndex::row(const std::string& name) const {
  const auto it = rows_.find(name);
  return it == rows_.end() ? nullptr : it->second;
}

const Param* TargetIndex::param(const std::string& pointer) const {
  const auto it = params_.find(pointer);
  return it == params_.end() ? nullptr : it->second;
}

std::optional<std::string> TargetIndex::anchor_pointer(const ToleranceTarget& target) const {
  if (target.anchor == ToleranceAnchor::Surface) {
    if (const SurfaceEntry* s = surface(target.name)) return s->pointer;
  } else if (target.anchor == ToleranceAnchor::Node) {
    if (const NodeEntry* n = node(target.name)) return n->pointer;
  }
  return std::nullopt;
}

bool is_rigid(ToleranceKind kind) noexcept {
  switch (kind) {
    case ToleranceKind::DecenterX:
    case ToleranceKind::DecenterY:
    case ToleranceKind::ShiftZ:
    case ToleranceKind::TiltXDeg:
    case ToleranceKind::TiltYDeg:
    case ToleranceKind::TiltZDeg:
      return true;
    default:
      return false;
  }
}

std::optional<std::string> quantity(const TargetIndex& index,
                                    const ToleranceTarget& target,
                                    ToleranceKind kind) {
  if (kind == ToleranceKind::Parameter) {
    if (target.anchor != ToleranceAnchor::Parameter || index.row(target.name) == nullptr) {
      return std::nullopt;
    }
    return "row:" + target.name;
  }
  const std::optional<std::string> anchor = index.anchor_pointer(target);
  if (!anchor) return std::nullopt;
  switch (kind) {
    case ToleranceKind::Radius:
    case ToleranceKind::CurvatureRelative:
      if (target.anchor != ToleranceAnchor::Surface) return std::nullopt;
      return *anchor + "/shape/base/radius";
    case ToleranceKind::Conic:
      if (target.anchor != ToleranceAnchor::Surface) return std::nullopt;
      return *anchor + "/shape/base/conic";
    case ToleranceKind::Thickness:
      return *anchor + "/pose/position/2";
    case ToleranceKind::Param:
      if (!target.param) return std::nullopt;
      return *anchor + *target.param;
    case ToleranceKind::Index:
    case ToleranceKind::Abbe:
      return *anchor + ":" + std::string(tolerance_kind_name(kind)) + ":" +
             (target.segment ? std::to_string(*target.segment) : std::string("-"));
    default:  // rigid
      return *anchor + ":" + std::string(tolerance_kind_name(kind));
  }
}

MaterialClass classify_material(std::string_view reference) {
  if (reference == "VACUUM") return MaterialClass::RealConst;
  constexpr std::string_view kConst = "CONST:";
  if (!reference.starts_with(kConst)) return MaterialClass::Other;
  const std::string_view args = reference.substr(kConst.size());
  const std::size_t comma = args.find(',');
  if (comma == std::string_view::npos) return MaterialClass::RealConst;
  // "CONST:<n>,<kappa>" (rtt-material, material.cpp): a metal if kappa > 0; a malformed
  // reference is left to compile (material.unknown).
  const std::string_view kappa_text = args.substr(comma + 1);
  double kappa = 0.0;
  const auto [end, ec] =
      std::from_chars(kappa_text.data(), kappa_text.data() + kappa_text.size(), kappa);
  if (ec != std::errc{} || end != kappa_text.data() + kappa_text.size()) {
    return MaterialClass::Other;
  }
  return kappa > 0.0 ? MaterialClass::Metal : MaterialClass::RealConst;
}

}  // namespace rtt::model::detail
