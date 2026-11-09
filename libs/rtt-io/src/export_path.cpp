#include "export_path.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/material/material.hpp"
#include "rtt/math/isometry.hpp"
#include "rtt/model/parameters.hpp"
#include "rtt/model/path.hpp"
#include "rtt/model/pose.hpp"

namespace rtt::io {

ExportError::ExportError(std::string code, std::string surface, const std::string& message)
    : std::runtime_error((surface.empty() ? std::string() : surface + ": ") + message + " [" +
                         code + "]"),
      code_(std::move(code)),
      surface_(std::move(surface)) {}

namespace detail {

void export_fail(diagnostics::DiagnosticCode code,
                 const std::string& surface,
                 const std::string& message) {
  throw ExportError(std::string(code.str()), surface, message);
}

model::Diagnostic export_warning(diagnostics::DiagnosticCode code,
                                 std::string location,
                                 std::string message) {
  return {model::Severity::Warning, std::move(location), std::move(message),
          std::string(code.str())};
}

namespace {

using model::Element;
using model::ElementKind;
using model::EventKind;

/// Threshold of the coaxial subset (ADR 0045, point 2): mm for the vertex, dimensionless on the
/// unit vectors of the frame. Chosen like the symmetry check of paraxial::first_order, whose
/// constant lies in a private header of rtt-paraxial.
constexpr double kCoaxial = 1e-12;

/// A surface of the model tree in tree order (pre-order), as CompiledSystem::surfaces().
struct TreeSurface {
  const model::Surface* surface = nullptr;
  const Element* element = nullptr;
  /// t of the pose that places this surface relative to the surface before it, for a canonical
  /// system: the element pose for the first surface of an element, else the surface pose.
  double t = 0.0;
};

void collect(const model::Assembly& a, std::vector<TreeSurface>& out) {
  for (const model::Node& n : a.children) {
    if (const auto* e = std::get_if<Element>(&n.value)) {
      for (std::size_t i = 0; i < e->surfaces.size(); ++i) {
        const model::Pose& p = i == 0 ? e->pose : e->surfaces[i].pose;
        out.push_back({&e->surfaces[i], e, p.position[2].value});
      }
    } else {
      collect(std::get<model::Assembly>(n.value), out);
    }
  }
}

/// True for a pose (0, 0, t) without rotation, relative to `reference`.
bool plain_shift(const model::Pose& p, model::PoseReference reference) {
  return p.reference == reference && p.position[0].value == 0.0 && p.position[1].value == 0.0 &&
         p.rotation_deg[0].value == 0.0 && p.rotation_deg[1].value == 0.0 &&
         p.rotation_deg[2].value == 0.0;
}

/// The canonical form of ADR 0045, point 6 (on the resolved system).
bool canonical(const model::System& s) {
  if (!s.root.pose.is_identity() || s.object.reference) return false;
  if (s.paths.size() != 1 || !s.paths.front().automatic) return false;
  bool first = true;
  for (const model::Node& n : s.root.children) {
    const auto* e = std::get_if<Element>(&n.value);
    if (e == nullptr) return false;
    const auto& m = e->segment_materials;
    if (!m.empty() && std::all_of(m.begin(), m.end(), [&](const auto& x) { return x == m[0]; })) {
      return false;
    }
    if (first) {
      if (!plain_shift(e->pose, model::PoseReference::Absolute) ||
          e->pose.position[2].value != 0.0) {
        return false;
      }
    } else if (!plain_shift(e->pose, model::PoseReference::RelativeToPreceding)) {
      return false;
    }
    first = false;
    for (std::size_t i = 0; i < e->surfaces.size(); ++i) {
      const model::Pose& p = e->surfaces[i].pose;
      if (i == 0 ? !p.is_identity() : !plain_shift(p, model::PoseReference::RelativeToPreceding)) {
        return false;
      }
    }
  }
  return true;
}

/// The kind of the event the automatic path has at an element of this kind (compile()).
EventKind automatic_kind(ElementKind k) {
  switch (k) {
    case ElementKind::Lens:
    case ElementKind::Plate:
      return EventKind::Refract;
    case ElementKind::Mirror:
      return EventKind::Reflect;
    case ElementKind::ThinElement:
    case ElementKind::Stop:
    case ElementKind::Detector:
      return EventKind::Transmit;
  }
  return EventKind::Refract;
}

/// The shape seen from the global frame of a surface with the axis -e_z: z_global = -f(r), the
/// root of the conic holds only c^2, so -f has the curvature -c, the same conic and -A_i
/// (ADR 0045, point 2). + 0.0 keeps a zero coefficient +0.
void mirror_shape(model::ShapeStack& shape) {
  const auto negate = [](model::Param& p) { p.value = -p.value + 0.0; };
  if (auto* c = std::get_if<model::Conic>(&shape.base)) {
    negate(c->radius);
  } else if (auto* a = std::get_if<model::EvenAsphere>(&shape.base)) {
    negate(a->radius);
    for (model::Param& x : a->coefficients) negate(x);
  } else if (auto* o = std::get_if<model::OddAsphere>(&shape.base)) {
    negate(o->curvature);
    for (model::Param& x : o->coefficients) negate(x);
  }
  // other shapes have no export mapping (export.surface_unsupported in the writers)
}

/// Checks the interaction of one surface on the path and drops a coating (ADR 0045, point 4).
void check_interaction(ExportRow& row, const Element& e, std::vector<model::Diagnostic>& w) {
  const std::string& id = row.surface.id.str();
  if (!row.surface.phases.empty()) {
    export_fail("export.interaction_unsupported", id, "a phase layer or grating has no mapping");
  }
  if (e.crystal) {
    export_fail("export.interaction_unsupported", id, "a crystal has no mapping");
  }
  model::Interaction& i = row.surface.interaction;
  if (std::holds_alternative<model::CoatingRef>(i) ||
      std::holds_alternative<model::IdealAntiReflection>(i)) {
    w.push_back(export_warning("export.coating_dropped", row.location,
                               "coating of " + id + " not exported"));
    i = model::Fresnel{};
  } else if (std::holds_alternative<model::IdealMirror>(i)) {
    if (e.kind != ElementKind::Mirror) {
      export_fail("export.interaction_unsupported", id,
                  "an ideal mirror outside a mirror element has no mapping");
    }
  } else if (std::holds_alternative<model::IdealLens>(i)) {
    if (e.kind != ElementKind::ThinElement) {
      export_fail("export.interaction_unsupported", id,
                  "an ideal lens outside a thin element has no mapping");
    }
  } else if (!std::holds_alternative<model::Fresnel>(i)) {
    export_fail("export.interaction_unsupported", id, "the interaction has no mapping");
  }
}

}  // namespace

ExportPath walk_export_path(const model::System& system,
                            const material::MaterialLibrary& materials,
                            const ExportOptions& options) {
  ExportPath out;

  // Configuration (ADR 0045, point 3; F6).
  const std::size_t count = model::configuration_count(system);
  if (!options.configuration && count > 1) {
    export_fail("export.multi_config", "",
                "the system has " + std::to_string(count) + " configurations; choose one");
  }
  const std::size_t configuration = options.configuration.value_or(0);
  if (configuration >= count) {
    throw std::invalid_argument("configuration " + std::to_string(configuration) +
                                " is beyond the " + std::to_string(count) + " of the system");
  }

  // Compile without coatings: they are dropped anyway, and no coating library is needed.
  model::System bare = system;
  const auto strip = [](const auto& self, model::Assembly& a) -> void {
    for (model::Node& n : a.children) {
      if (auto* e = std::get_if<Element>(&n.value)) {
        for (model::Surface& s : e->surfaces) {
          if (std::holds_alternative<model::CoatingRef>(s.interaction)) {
            s.interaction = model::Fresnel{};
          }
        }
      } else {
        self(self, std::get<model::Assembly>(n.value));
      }
    }
  };
  strip(strip, bare.root);
  const compile::CompiledSystem cs = compile::compile(bare, materials, configuration);
  out.system = model::resolve_parameters(system, configuration);

  // Path (F6).
  compile::PathId path_id{0};
  if (options.path) {
    const auto found = cs.find_path(*options.path);
    if (!found) throw std::invalid_argument("the system has no path '" + *options.path + "'");
    path_id = *found;
  }
  if (cs.paths().size() > 1) {
    out.warnings.push_back(
        export_warning("export.paths_dropped", "/paths",
                       "only the path '" + cs.path(path_id).name + "' is exported, " +
                           std::to_string(cs.paths().size() - 1) + " other paths are dropped"));
  }
  const compile::CompiledPath& path = cs.path(path_id);

  std::vector<TreeSurface> tree;
  collect(out.system.root, tree);
  if (tree.size() != cs.surfaces().size()) {
    throw std::logic_error("export: the model tree and the compiled surfaces differ");
  }

  const std::uint32_t env = cs.environment_medium();
  const std::string& env_name = cs.media()[env].reference;

  // Object space before the path (ADR 0045, point 3).
  if (out.system.object.medium && *out.system.object.medium != env_name) {
    const compile::CompiledEvent* first = path.events.empty() ? nullptr : &path.events.front();
    const bool lens =
        first != nullptr && first->kind == EventKind::Refract && first->medium_after != env;
    if (!lens) {
      export_fail("export.object_medium_unsupported",
                  first != nullptr ? cs.surfaces()[first->surface].id.str() : std::string(),
                  "glass in object space needs a path that starts in a lens");
    }
    out.object_medium = out.system.object.medium;
  }

  std::set<std::uint32_t> visited;
  for (std::size_t k = 0; k < path.events.size(); ++k) {
    const compile::CompiledEvent& e = path.events[k];
    const compile::CompiledSurface& cs_surface = cs.surfaces()[e.surface];
    const TreeSurface& t = tree[e.surface];
    const Element& element = *t.element;
    const std::string& id = cs_surface.id.str();

    // refract, the default kind of an event, keeps the medium at a stop or detector (elements
    // without material) as transmit does (ADR 0045, addendum PR B)
    const bool same_as_transmit =
        e.kind == EventKind::Refract &&
        (element.kind == ElementKind::Stop || element.kind == ElementKind::Detector);
    if ((e.kind != automatic_kind(element.kind) && !same_as_transmit) || e.order != 0) {
      export_fail("export.event_unsupported", id,
                  "a ZMX surface holds only the event of the automatic path of its element");
    }
    if (!visited.insert(e.surface).second) {
      export_fail("export.event_unsupported", id, "the path visits the surface twice");
    }
    if (element.kind == ElementKind::Detector && k + 1 != path.events.size()) {
      export_fail("export.event_unsupported", id, "a detector before the last event");
    }

    // Coaxial (ADR 0045, point 2).
    const math::Vec3& p = cs_surface.to_global.translation();
    const math::Mat3& r = cs_surface.to_global.rotation();
    bool coaxial = std::abs(p.x()) <= kCoaxial && std::abs(p.y()) <= kCoaxial;
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        if (i != j && std::abs(r(i, j)) > kCoaxial) coaxial = false;
      }
    }
    if (!coaxial) {
      export_fail("export.not_coaxial", id,
                  "only coaxial systems are exported: the vertex must lie on the z axis and the "
                  "surface axes along +-x, +-y, +-z");
    }

    const bool glass_before = e.medium_before != env;
    const bool glass_after = e.medium_after != env;
    if (e.kind == EventKind::Reflect && glass_before) {
      export_fail("export.mirror_in_glass", id, "a mirror inside glass (Mangin)");
    }
    if ((element.kind == ElementKind::Stop || element.kind == ElementKind::ThinElement ||
         element.kind == ElementKind::Detector) &&
        (glass_before || glass_after)) {
      export_fail("export.element_in_glass", id,
                  "a stop, thin element or detector needs the environment on both sides");
    }

    ExportRow row;
    row.surface = *t.surface;
    row.element = element.name;
    row.kind = element.kind;
    row.location = cs_surface.location;
    row.z = p.z();
    row.mirror = e.kind == EventKind::Reflect;
    row.stop = cs.stop_surface() && *cs.stop_surface() == e.surface;
    if (glass_after && !row.mirror) row.medium_after = cs.media()[e.medium_after].reference;
    if (r(2, 2) < 0.0) mirror_shape(row.surface.shape);
    check_interaction(row, element, out.warnings);
    if (element.kind == ElementKind::Mirror &&
        (element.material || !element.segment_materials.empty())) {
      out.warnings.push_back(export_warning("export.mirror_material_dropped", row.location,
                                            "substrate of " + element.name + " not exported"));
    }
    out.rows.push_back(std::move(row));
  }

  if (out.rows.empty() || out.rows.back().kind != ElementKind::Detector) {
    export_fail("export.image_missing",
                out.rows.empty() ? std::string() : out.rows.back().surface.id.str(),
                "the last event of the path must be a detector (the image surface)");
  }
  const auto stop =
      std::find_if(out.rows.begin(), out.rows.end(), [](const ExportRow& x) { return x.stop; });
  if (stop == out.rows.end()) {
    export_fail("export.stop_missing", "", "the exported path has no stop");
  }
  if (stop->surface.aperture &&
      !std::holds_alternative<model::CircularAperture>(*stop->surface.aperture)) {
    export_fail("export.aperture_unsupported", stop->surface.id.str(),
                "a rectangular or elliptical stop has no mapping");
  }

  for (std::size_t i = 0; i < cs.surfaces().size(); ++i) {
    if (!visited.contains(static_cast<std::uint32_t>(i))) {
      out.warnings.push_back(
          export_warning("export.surface_off_path", cs.surfaces()[i].location,
                         "surface " + cs.surfaces()[i].id.str() + " is not on the exported path"));
    }
  }

  // Thicknesses (ADR 0045, point 4): the t of the poses if the system is canonical (then the
  // path is automatic and visits every surface in tree order), else differences of z.
  const bool exact = canonical(out.system);
  for (std::size_t k = 0; k + 1 < out.rows.size(); ++k) {
    out.rows[k].thickness = exact ? tree[k + 1].t : out.rows[k + 1].z - out.rows[k].z;
  }
  if (!out.system.object.at_infinity) {
    out.object_distance =
        exact ? out.system.object.distance.value : out.rows.front().z - cs.object_z();
  }
  return out;
}

}  // namespace detail
}  // namespace rtt::io
