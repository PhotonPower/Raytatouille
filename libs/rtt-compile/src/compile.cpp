#include <cstddef>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "rtt/compile/compiled_system.hpp"

namespace rtt::compile {
namespace {

std::string idx(const std::string& base, std::size_t i) {
  return base + "/" + std::to_string(i);
}

std::string join_errors(const std::vector<model::Diagnostic>& diagnostics) {
  std::string text = "cannot compile system";
  for (const auto& d : diagnostics) text += "\n  " + model::to_string(d);
  return text;
}

/// What the path builder needs to know about an element.
struct ElementInfo {
  model::ElementKind kind = model::ElementKind::Lens;
  std::optional<std::uint32_t> medium;  ///< medium inside; none if the element has no material
  std::string location;                 ///< JSON pointer of the element
  std::uint32_t first_surface = 0;      ///< surfaces [first_surface, first_surface + count)
  std::uint32_t surface_count = 0;
};

/// Collects the pieces of a CompiledSystem and all errors found on the way.
class Compiler {
 public:
  Compiler(const model::System& system, const material::MaterialLibrary& materials)
      : system_(system), materials_(materials) {
    for (const auto& w : system.wavelengths) wavelengths_um_.push_back(w.um);
  }

  void run() {
    if (wavelengths_um_.size() > std::numeric_limits<std::uint16_t>::max()) {
      error("/wavelengths", "at most 65535 wavelengths are supported (RayBatch::wl is 16 bit)");
    }
    // The environment is always media_[0] (CompiledSystem::environment_medium()).
    environment_ = medium(system_.environment.medium, "/environment/medium").value_or(0);
    add_assembly(system_.root, model::to_isometry(system_.root.pose), "/root");
    for (const model::Path& path : system_.paths) paths_.push_back(build_path(path));
  }

  std::vector<model::Diagnostic> errors_;
  std::vector<double> wavelengths_um_;
  std::vector<CompiledSurface> surfaces_;
  std::vector<CompiledMedium> media_;
  std::vector<CompiledPath> paths_;

 private:
  void error(std::string location, std::string message) {
    errors_.push_back({model::Severity::Error, std::move(location), std::move(message)});
  }

  /// Index of the medium for `reference`, resolving and evaluating it on first use.
  std::optional<std::uint32_t> medium(const std::string& reference, const std::string& location) {
    if (const auto it = medium_index_.find(reference); it != medium_index_.end()) {
      return it->second;
    }
    std::shared_ptr<const material::Material> material;
    try {
      material = materials_.resolve(reference);
    } catch (const material::UnknownMaterial& e) {
      error(location, e.what());
      return std::nullopt;
    }
    // Materials are evaluated once per wavelength into constants (ADR 0014).
    CompiledMedium m{reference, {}};
    m.index.reserve(wavelengths_um_.size());
    for (const double wl : wavelengths_um_) {
      m.index.push_back(material->index(wl, system_.environment.temperature_c));
    }
    const auto index = static_cast<std::uint32_t>(media_.size());
    media_.push_back(std::move(m));
    medium_index_.emplace(reference, index);
    return index;
  }

  void add_assembly(const model::Assembly& assembly,
                    const math::Isometry3& to_global,
                    const std::string& location) {
    for (std::size_t i = 0; i < assembly.children.size(); ++i) {
      const std::string child_location = idx(location + "/children", i);
      const auto& value = assembly.children[i].value;
      if (const auto* sub = std::get_if<model::Assembly>(&value)) {
        add_assembly(*sub, to_global * model::to_isometry(sub->pose), child_location);
      } else {
        const auto& element = std::get<model::Element>(value);
        add_element(element, to_global * model::to_isometry(element.pose), child_location);
      }
    }
  }

  void add_element(const model::Element& element,
                   const math::Isometry3& to_global,
                   const std::string& location) {
    ElementInfo info;
    info.kind = element.kind;
    info.location = location;
    info.first_surface = static_cast<std::uint32_t>(surfaces_.size());
    info.surface_count = static_cast<std::uint32_t>(element.surfaces.size());
    if (element.material) info.medium = medium(*element.material, location + "/material");
    const auto element_index = static_cast<std::uint32_t>(elements_.size());
    elements_.push_back(info);

    for (std::size_t j = 0; j < element.surfaces.size(); ++j) {
      const model::Surface& s = element.surfaces[j];
      const std::string surface_location = idx(location + "/surfaces", j);
      CompiledSurface c;
      c.id = s.id;
      c.element_kind = element.kind;
      c.element_name = element.name;
      c.to_global = to_global * model::to_isometry(s.pose);
      c.to_local = c.to_global.inverse();
      c.shape = compile_shape(s.shape, surface_location + "/shape");
      c.aperture = s.aperture;
      c.phases = s.phases;
      c.interaction = s.interaction;
      surface_index_.emplace(s.id, static_cast<std::uint32_t>(surfaces_.size()));
      surface_element_.push_back(element_index);
      surfaces_.push_back(std::move(c));
    }
  }

  CompiledShape compile_shape(const model::ShapeStack& shape, const std::string& location) {
    for (std::size_t k = 0; k < shape.terms.size(); ++k) {
      error(idx(location + "/terms", k), "Zernike sag terms are not supported before M8");
    }
    if (const auto* conic = std::get_if<model::Conic>(&shape.base)) {
      // c = 1/R (docs/architecture.md, Konventionen); validate() guarantees R finite, != 0.
      return geom::Conic<double>(1.0 / conic->radius.value, conic->conic.value);
    }
    if (std::holds_alternative<model::EvenAsphere>(shape.base)) {
      error(location + "/base", "even asphere: not yet supported by the tracer, see #6");
    }
    return geom::Plane<double>{};
  }

  CompiledPath build_path(const model::Path& path) {
    CompiledPath compiled{path.name, {}};
    if (path.automatic) {
      for (const ElementInfo& e : elements_) {
        const bool lens_or_plate =
            e.kind == model::ElementKind::Lens || e.kind == model::ElementKind::Plate;
        if (lens_or_plate && e.surface_count > 2) {
          if (cemented_reported_.emplace(e.location).second) {
            error(e.location + "/surfaces",
                  "cemented groups: not supported before M2, see #14 (automatic path through an "
                  "element with more than 2 surfaces)");
          }
          continue;
        }
        const model::EventKind kind = automatic_event(e.kind);
        for (std::uint32_t s = 0; s < e.surface_count; ++s) {
          compiled.events.push_back({e.first_surface + s, kind, 0, 0, 0});
        }
      }
    } else {
      for (const model::Event& event : path.events) {
        // validate() guarantees that every referenced surface exists.
        compiled.events.push_back(
            {surface_index_.at(event.surface), event.kind, event.order, 0, 0});
      }
    }
    assign_media(compiled.events);
    return compiled;
  }

  /// Lens and Plate refract, Mirror reflects, Stop, Detector and ThinElement transmit (#5).
  static model::EventKind automatic_event(model::ElementKind kind) {
    switch (kind) {
      case model::ElementKind::Lens:
      case model::ElementKind::Plate:
        return model::EventKind::Refract;
      case model::ElementKind::Mirror:
        return model::EventKind::Reflect;
      case model::ElementKind::ThinElement:
      case model::ElementKind::Stop:
      case model::ElementKind::Detector:
        return model::EventKind::Transmit;
    }
    return model::EventKind::Transmit;
  }

  /// Media before and after each event, rule decided for #5: the ray starts in the environment;
  /// Refract, Ordinary and Extraordinary at an element with a material switch between the
  /// inside of that element and the environment; every other event, and every event at an
  /// element without material, keeps the medium.
  void assign_media(std::vector<CompiledEvent>& events) const {
    std::optional<std::uint32_t> inside;  // element the ray is in, none = environment
    std::uint32_t current = environment_;
    for (CompiledEvent& event : events) {
      event.medium_before = current;
      const bool crosses = event.kind == model::EventKind::Refract ||
                           event.kind == model::EventKind::Ordinary ||
                           event.kind == model::EventKind::Extraordinary;
      const std::uint32_t element = surface_element_[event.surface];
      const auto& element_medium = elements_[element].medium;
      if (crosses && element_medium) {
        if (inside == element) {
          inside.reset();
          current = environment_;
        } else {
          inside = element;
          current = *element_medium;
        }
      }
      event.medium_after = current;
    }
  }

  const model::System& system_;
  const material::MaterialLibrary& materials_;
  std::uint32_t environment_ = 0;
  std::vector<ElementInfo> elements_;
  std::vector<std::uint32_t> surface_element_;  // owning element per surface
  std::map<model::SurfaceId, std::uint32_t> surface_index_;
  std::map<std::string, std::uint32_t, std::less<>> medium_index_;
  std::set<std::string> cemented_reported_;  // elements already reported as cemented groups
};

}  // namespace

CompileError::CompileError(std::vector<model::Diagnostic> diagnostics)
    : std::runtime_error(join_errors(diagnostics)), diagnostics_(std::move(diagnostics)) {}

CompiledSystem compile(const model::System& system, const material::MaterialLibrary& materials) {
  std::vector<model::Diagnostic> diagnostics = model::validate(system);
  if (model::has_errors(diagnostics)) {
    std::erase_if(diagnostics,
                  [](const model::Diagnostic& d) { return d.severity != model::Severity::Error; });
    throw CompileError(std::move(diagnostics));
  }

  Compiler compiler(system, materials);
  compiler.run();
  if (!compiler.errors_.empty()) throw CompileError(std::move(compiler.errors_));

  CompiledSystem cs;
  cs.wavelengths_um_ = std::move(compiler.wavelengths_um_);
  for (std::size_t i = 0; i < system.wavelengths.size(); ++i) {
    if (system.wavelengths[i].reference) cs.reference_wl_ = static_cast<std::uint16_t>(i);
  }
  cs.temperature_c_ = system.environment.temperature_c;
  cs.aperture_ = system.aperture;
  cs.fields_ = system.fields;
  cs.object_ = system.object;
  cs.surfaces_ = std::move(compiler.surfaces_);
  cs.media_ = std::move(compiler.media_);
  cs.paths_ = std::move(compiler.paths_);
  return cs;
}

std::optional<std::uint32_t> CompiledSystem::find_surface(const model::SurfaceId& id) const {
  for (std::size_t i = 0; i < surfaces_.size(); ++i) {
    if (surfaces_[i].id == id) return static_cast<std::uint32_t>(i);
  }
  return std::nullopt;
}

std::optional<PathId> CompiledSystem::find_path(std::string_view name) const {
  for (std::size_t i = 0; i < paths_.size(); ++i) {
    if (paths_[i].name == name) return PathId{static_cast<std::uint32_t>(i)};
  }
  return std::nullopt;
}

}  // namespace rtt::compile
