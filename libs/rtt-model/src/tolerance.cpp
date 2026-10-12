#include <algorithm>
#include <array>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <variant>
#include <vector>

#include "rtt/model/tolerance_classes.hpp"
#include "tolerance_targets.hpp"

namespace rtt::model {

std::string_view tolerance_kind_name(ToleranceKind kind) noexcept {
  switch (kind) {
    case ToleranceKind::Radius:
      return "radius";
    case ToleranceKind::CurvatureRelative:
      return "curvature_relative";
    case ToleranceKind::Conic:
      return "conic";
    case ToleranceKind::Thickness:
      return "thickness";
    case ToleranceKind::Param:
      return "param";
    case ToleranceKind::Parameter:
      return "parameter";
    case ToleranceKind::DecenterX:
      return "decenter_x";
    case ToleranceKind::DecenterY:
      return "decenter_y";
    case ToleranceKind::ShiftZ:
      return "shift_z";
    case ToleranceKind::TiltXDeg:
      return "tilt_x_deg";
    case ToleranceKind::TiltYDeg:
      return "tilt_y_deg";
    case ToleranceKind::TiltZDeg:
      return "tilt_z_deg";
    case ToleranceKind::Index:
      return "index";
    case ToleranceKind::Abbe:
      return "abbe";
  }
  return "?";
}

std::uint64_t tolerance_id_hash(std::string_view id) noexcept {
  // FNV-1a, 64 bit (RFC 9923, Sec. 2: hash = offset_basis; per octet hash = hash xor octet,
  // then hash = hash * FNV_Prime mod 2^64; Sec. 5: the 64-bit constants).
  constexpr std::uint64_t kOffsetBasis = 0xCBF29CE484222325ULL;
  constexpr std::uint64_t kPrime = 0x00000100000001B3ULL;
  std::uint64_t hash = kOffsetBasis;
  for (const char c : id) {
    hash ^= static_cast<std::uint64_t>(static_cast<unsigned char>(c));
    hash *= kPrime;  // unsigned arithmetic: mod 2^64
  }
  return hash;
}

namespace {

using detail::MaterialClass;

/// The kinds a class generates per level, in the order of ADR 0035, point 5.
constexpr std::array kSurfaceKinds = {ToleranceKind::CurvatureRelative, ToleranceKind::Conic,
                                      ToleranceKind::Thickness};
constexpr std::array kElementKinds = {ToleranceKind::Thickness, ToleranceKind::DecenterX,
                                      ToleranceKind::DecenterY, ToleranceKind::TiltXDeg,
                                      ToleranceKind::TiltYDeg};
constexpr std::array kMaterialKinds = {ToleranceKind::Index, ToleranceKind::Abbe};
constexpr std::array kAssemblyKinds = {ToleranceKind::DecenterX, ToleranceKind::DecenterY,
                                       ToleranceKind::ShiftZ, ToleranceKind::TiltXDeg,
                                       ToleranceKind::TiltYDeg};

ClassSkip material_skip(const std::string& reference, ToleranceKind kind) {
  switch (detail::classify_material(reference)) {
    case MaterialClass::Metal:
      return ClassSkip::Metal;
    case MaterialClass::RealConst:
      return kind == ToleranceKind::Abbe ? ClassSkip::ConstAbbe : ClassSkip::None;
    case MaterialClass::Other:
      break;
  }
  return ClassSkip::None;
}

/// True if every entry of a segment material list is the same string (ADR 0035, addendum #254).
bool uniform(const std::vector<std::string>& materials) {
  return std::all_of(materials.begin(), materials.end(),
                     [&](const std::string& m) { return m == materials.front(); });
}

ClassSkip bound_skip(const Param& p) {
  return p.is_bound() ? ClassSkip::Bound : ClassSkip::None;
}

class ClassGenerator {
 public:
  explicit ClassGenerator(const System& system) : system_(system), index_(system) {
    for (const ToleranceOverride& o : system.tolerances.overrides) {
      if (const auto* c = std::get_if<ToleranceClassOverride>(&o)) {
        classes_.emplace(c->node, c->tolerance_class);  // the first one counts
      }
    }
  }

  ClassTolerances run() {
    walk(system_.root, system_.tolerances.tolerance_class);
    apply_items();
    apply_overrides();
    return std::move(out_);
  }

 private:
  ToleranceClass class_of(const std::string& name, ToleranceClass inherited) const {
    const auto it = classes_.find(name);
    return it == classes_.end() ? inherited : it->second;
  }

  void walk(const Assembly& a, ToleranceClass inherited) {
    for (const Node& n : a.children) {
      if (const auto* sub = std::get_if<Assembly>(&n.value)) {
        const ToleranceClass c = class_of(sub->name, inherited);
        effective_.emplace(sub->name, c);
        if (c != ToleranceClass::None) {
          for (const ToleranceKind k : kAssemblyKinds) {
            add(node_target(sub->name), k, ToleranceLevel::Assembly, c, ClassSkip::None);
          }
        }
        walk(*sub, c);
      } else {
        const Element& e = std::get<Element>(n.value);
        const ToleranceClass c = class_of(e.name, inherited);
        effective_.emplace(e.name, c);
        if (c != ToleranceClass::None) element(e, c);
      }
    }
  }

  void element(const Element& e, ToleranceClass c) {
    for (const ToleranceKind k : kElementKinds) {
      const ClassSkip skip =
          k == ToleranceKind::Thickness ? bound_skip(e.pose.position[2]) : ClassSkip::None;
      add(node_target(e.name), k, ToleranceLevel::Element, c, skip);
    }
    const bool in_beam = e.kind == ElementKind::Lens || e.kind == ElementKind::Plate;
    for (const ToleranceKind k : kMaterialKinds) {
      if (!in_beam) {
        add(node_target(e.name), k, ToleranceLevel::Element, c, ClassSkip::NoMaterialInBeam);
      } else if (e.crystal) {
        add(node_target(e.name), k, ToleranceLevel::Element, c, ClassSkip::Crystal);
      } else if (!e.segment_materials.empty()) {
        // A uniform list is one material and one body (ADR 0035, addendum #254): segment 0.
        const std::size_t n = uniform(e.segment_materials) ? 1 : e.segment_materials.size();
        for (std::size_t i = 0; i < n; ++i) {
          ToleranceTarget t = node_target(e.name);
          t.segment = static_cast<std::uint32_t>(i);
          add(t, k, ToleranceLevel::Element, c, material_skip(e.segment_materials[i], k));
        }
      } else if (const std::optional<std::string>& m = e.material) {
        add(node_target(e.name), k, ToleranceLevel::Element, c, material_skip(*m, k));
      } else {
        add(node_target(e.name), k, ToleranceLevel::Element, c, ClassSkip::NoMaterial);
      }
    }
    for (std::size_t j = 0; j < e.surfaces.size(); ++j) surface(e.surfaces[j], j, c);
  }

  void surface(const Surface& s, std::size_t j, ToleranceClass c) {
    ToleranceTarget t;
    t.anchor = ToleranceAnchor::Surface;
    t.name = s.id.str();
    ClassSkip radius = ClassSkip::Plane;
    ClassSkip conic = ClassSkip::Plane;
    if (const auto* q = std::get_if<Conic>(&s.shape.base)) {
      radius = bound_skip(q->radius);
      conic = q->conic.is_bound()     ? ClassSkip::Bound
              : q->conic.value == 0.0 ? ClassSkip::Sphere
                                      : ClassSkip::None;
    } else if (const auto* a = std::get_if<EvenAsphere>(&s.shape.base)) {
      radius = bound_skip(a->radius);
      conic = bound_skip(a->conic);
    }
    add(t, ToleranceKind::CurvatureRelative, ToleranceLevel::Surface, c, radius);
    add(t, ToleranceKind::Conic, ToleranceLevel::Surface, c, conic);
    // The thickness of the segment before surface j >= 1 (ADR 0035, point 5).
    if (j >= 1) {
      add(t, ToleranceKind::Thickness, ToleranceLevel::Surface, c, bound_skip(s.pose.position[2]));
    }
  }

  static ToleranceTarget node_target(const std::string& name) {
    ToleranceTarget t;
    t.anchor = ToleranceAnchor::Node;
    t.name = name;
    return t;
  }

  void add(const ToleranceTarget& t,
           ToleranceKind kind,
           ToleranceLevel level,
           ToleranceClass c,
           ClassSkip skip) {
    ClassTolerance item;
    item.id = std::string(t.anchor == ToleranceAnchor::Surface ? "class:surface:" : "class:node:") +
              t.name + ":" + std::string(tolerance_kind_name(kind));
    if (t.segment) item.id += ":" + std::to_string(*t.segment);
    item.target = t;
    item.kind = kind;
    item.level = level;
    item.tolerance_class = c;
    item.skip = skip;
    out_.items.push_back(std::move(item));
  }

  /// Explicit items take precedence by perturbed quantity (ADR 0035, points 2 and 5).
  void apply_items() {
    std::unordered_map<std::string, std::size_t> by_quantity;
    const std::vector<ToleranceItem>& items = system_.tolerances.items;
    for (std::size_t i = 0; i < items.size(); ++i) {
      if (const std::optional<std::string> q =
              detail::quantity(index_, items[i].target, items[i].kind)) {
        by_quantity.emplace(*q, i);  // the first one counts
      }
    }
    for (ClassTolerance& c : out_.items) {
      if (const std::optional<std::string> q = detail::quantity(index_, c.target, c.kind)) {
        if (const auto it = by_quantity.find(*q); it != by_quantity.end()) {
          c.replaced_by = it->second;
        }
      }
    }
  }

  void apply_overrides() {
    const std::vector<ToleranceOverride>& overrides = system_.tolerances.overrides;
    for (std::size_t i = 0; i < overrides.size(); ++i) {
      const auto* v = std::get_if<ToleranceValueOverride>(&overrides[i]);
      if (v == nullptr) continue;
      out_.value_overrides.push_back(status(i, *v));
    }
  }

  OverrideStatus status(std::size_t i, const ToleranceValueOverride& v) {
    for (ClassTolerance& c : out_.items) {
      if (c.kind != v.kind || c.target.anchor != v.target.anchor ||
          c.target.name != v.target.name || c.target.segment != v.target.segment) {
        continue;
      }
      if (c.skip != ClassSkip::None) return {i, OverrideEffect::Skipped, c.skip};
      if (c.replaced_by) return {i, OverrideEffect::Replaced, ClassSkip::None};
      if (!c.value_override) c.value_override = i;  // the first one counts
      return {i, OverrideEffect::Applied, ClassSkip::None};
    }
    // No such item: the class none at the target, or a kind or target no class generates.
    std::string node = v.target.name;
    if (v.target.anchor == ToleranceAnchor::Surface) {
      const detail::SurfaceEntry* s = index_.surface(v.target.name);
      if (s == nullptr) return {i, OverrideEffect::NotGenerated, ClassSkip::None};
      node = s->element->name;
    }
    const auto it = effective_.find(node);
    if (it != effective_.end() && it->second == ToleranceClass::None && generatable(v)) {
      return {i, OverrideEffect::ClassNone, ClassSkip::None};
    }
    return {i, OverrideEffect::NotGenerated, ClassSkip::None};
  }

  /// True if a class generates the kind of `v` at the level of its target.
  bool generatable(const ToleranceValueOverride& v) const {
    const auto contains = [&](const auto& kinds) {
      for (const ToleranceKind k : kinds) {
        if (k == v.kind) return true;
      }
      return false;
    };
    if (v.target.anchor == ToleranceAnchor::Surface) return contains(kSurfaceKinds);
    const detail::NodeEntry* n = index_.node(v.target.name);
    if (n == nullptr) return false;
    if (n->assembly != nullptr) return contains(kAssemblyKinds);
    return contains(kElementKinds) || contains(kMaterialKinds);
  }

  const System& system_;
  detail::TargetIndex index_;
  std::unordered_map<std::string, ToleranceClass> classes_;    ///< class overrides by node
  std::unordered_map<std::string, ToleranceClass> effective_;  ///< effective class by node
  ClassTolerances out_;
};

}  // namespace

ClassTolerances class_tolerances(const System& system) {
  return ClassGenerator(system).run();
}

}  // namespace rtt::model
