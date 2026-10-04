#include <algorithm>
#include <array>
#include <charconv>
#include <cstddef>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <system_error>
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
  /// Lens and Plate: one medium per segment (segment i between surfaces i and i + 1, ADR 0017).
  /// Mirror: the substrate, if any, as one body. Empty if the element has no material.
  std::vector<std::uint32_t> media;
  /// JSON pointer of the material of each entry of `media`: .../material for the shorthand,
  /// .../material/<i> for a list entry.
  std::vector<std::string> media_locations;
  bool segmented = false;           ///< Lens or Plate: the media follow the segments
  std::string location;             ///< JSON pointer of the element
  std::uint32_t first_surface = 0;  ///< surfaces [first_surface, first_surface + count)
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
    for (std::size_t p = 0; p < system_.paths.size(); ++p) {
      paths_.push_back(build_path(system_.paths[p], idx("/paths", p)));
    }
    check_wavelength_ranges();
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
      m.index.push_back(
          material->index(wl, system_.environment.temperature_c, system_.environment.pressure_atm));
    }
    const auto index = static_cast<std::uint32_t>(media_.size());
    media_.push_back(std::move(m));
    media_checks_.push_back({material->wavelength_range_um(), location});
    medium_index_.emplace(reference, index);
    return index;
  }

  /// Every system wavelength must lie in the valid range of every medium that a path uses
  /// (decided for #23); unused media are not checked. One error per medium, at the first place
  /// on a path where the ray meets it: /environment/medium or the material of an element.
  /// Runs only if nothing failed before: after an unresolved material the medium indices of
  /// the events are placeholders that may not even exist in media_.
  void check_wavelength_ranges() {
    if (!errors_.empty()) return;
    std::vector<std::optional<std::string>> where(media_.size());
    for (const CompiledPath& path : paths_) {
      for (const CompiledEvent& event : path.events) {
        const ElementInfo& element = elements_[surface_element_[event.surface]];
        for (const std::uint32_t m : {event.medium_before, event.medium_after}) {
          if (where[m]) continue;
          if (m == environment_) {
            where[m] = "/environment/medium";
          } else if (const auto j = std::find(element.media.begin(), element.media.end(), m);
                     j != element.media.end()) {
            // The material of this element: shorthand or list entry (ADR 0017, #27).
            where[m] = element.media_locations[static_cast<std::size_t>(j - element.media.begin())];
          } else {
            where[m] = media_checks_[m].location;  // medium of an element entered earlier
          }
        }
      }
    }
    for (std::size_t m = 0; m < media_.size(); ++m) {
      const std::optional<material::WavelengthRange>& range = media_checks_[m].range;
      const std::optional<std::string>& location = where[m];
      if (!location.has_value() || !range.has_value()) continue;
      const material::WavelengthRange& valid = range.value();
      for (const double wl : wavelengths_um_) {
        if (valid.contains(wl)) continue;
        error(location.value(), "wavelength " + number(wl) + " um is outside the valid range [" +
                                    number(valid.min_um) + ", " + number(valid.max_um) +
                                    "] um of material '" + media_[m].reference + "'");
        break;
      }
    }
  }

  /// Shortest text that reads back to the same double, independent of the global locale.
  static std::string number(double value) {
    std::array<char, 32> buffer{};
    const auto [end, ec] = std::to_chars(buffer.data(), buffer.data() + buffer.size(), value);
    return ec == std::errc{} ? std::string(buffer.data(), end) : std::string("?");
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
    // Unresolved materials are reported by medium(); index 0 is only a placeholder then.
    info.segmented =
        element.kind == model::ElementKind::Lens || element.kind == model::ElementKind::Plate;
    if (info.segmented && !element.segment_materials.empty()) {
      for (std::size_t i = 0; i < element.segment_materials.size(); ++i) {
        info.media_locations.push_back(idx(location + "/material", i));
        info.media.push_back(
            medium(element.segment_materials[i], info.media_locations.back()).value_or(0));
      }
    } else if (element.material) {
      // The shorthand holds for all segments; resolved once, reported once.
      const std::uint32_t m = medium(*element.material, location + "/material").value_or(0);
      const std::size_t count = info.segmented ? element.surfaces.size() - 1 : 1;
      info.media.assign(count, m);
      info.media_locations.assign(count, location + "/material");
    }
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
    if (const auto* asphere = std::get_if<model::EvenAsphere>(&shape.base)) {
      // Same c = 1/R as the conic; coefficients[0] = A4 in model and rtt-geom.
      std::vector<double> coefficients;
      coefficients.reserve(asphere->coefficients.size());
      for (const model::Param& a : asphere->coefficients) coefficients.push_back(a.value);
      return geom::EvenAsphere<double>(1.0 / asphere->radius.value, asphere->conic.value,
                                       std::move(coefficients));
    }
    return geom::Plane<double>{};
  }

  /// `location` is the JSON pointer of the path in the system file.
  CompiledPath build_path(const model::Path& path, const std::string& location) {
    CompiledPath compiled{path.name, {}};
    if (path.automatic) {
      for (const ElementInfo& e : elements_) {
        // A mirror on a substrate with several surfaces (Mangin mirror) refracts at its front
        // surface; reflecting at every surface would be wrong (#6).
        if (e.kind == model::ElementKind::Mirror && !e.media.empty() && e.surface_count > 1) {
          if (mangin_reported_.emplace(e.location).second) {
            error(e.location + "/surfaces",
                  "mirror with substrate material on the automatic path: use an explicit path "
                  "(Refract, Reflect, Refract)");
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
    assign_media(compiled.events, location + "/events");
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

  /// Media before and after each event (rules decided for #5 and #27, docs/architecture.md,
  /// "Medien entlang eines Pfads"). The ray starts in the environment. Reflect, Transmit and
  /// Diffract keep the medium, and so does every event at an element without material.
  /// Refract, Ordinary and Extraordinary
  /// - at surface i of a Lens or Plate with N surfaces (segment j between surfaces j and j + 1):
  ///   from segment i - 1 into segment i (the environment if i = N - 1), from segment i into
  ///   segment i - 1 (the environment if i = 0), from any other segment of the same element into
  ///   the environment; from outside the element through the first surface into segment 0,
  ///   through the last surface into segment N - 2, and through an inner surface into its two
  ///   neighbouring segments if they have the same material, otherwise it is an error (the
  ///   side is ambiguous);
  /// - at a Mirror with substrate: toggle between the substrate and the environment (#6).
  /// Outside an element means in the environment or in another element: elements do not nest,
  /// so entering B while in A leaves A, and leaving B goes to the environment.
  /// `location` is the JSON pointer of the events; errors add the event index.
  void assign_media(std::vector<CompiledEvent>& events, const std::string& location) {
    std::optional<std::uint32_t> inside;  // element the ray is in, none = environment
    std::uint32_t segment = 0;            // segment of `inside` the ray is in
    std::uint32_t current = environment_;
    for (std::size_t k = 0; k < events.size(); ++k) {
      CompiledEvent& event = events[k];
      event.medium_before = current;
      const bool crosses = event.kind == model::EventKind::Refract ||
                           event.kind == model::EventKind::Ordinary ||
                           event.kind == model::EventKind::Extraordinary;
      const std::uint32_t element = surface_element_[event.surface];
      const ElementInfo& info = elements_[element];
      if (crosses && !info.media.empty()) {
        const std::uint32_t i = event.surface - info.first_surface;  // surface in the element
        const std::uint32_t last = info.surface_count - 1;
        bool leave = false;
        if (!info.segmented) {
          leave = inside == element;
          segment = 0;
        } else if (inside == element && i == segment + 1) {  // forward through surface i
          leave = i == last;
          segment = i;
        } else if (inside == element && i == segment) {  // backward through surface i
          leave = i == 0;
          if (!leave) segment = i - 1;
        } else if (inside == element) {
          leave = true;  // the surface does not bound the segment (prism on an explicit path)
        } else if (i == 0 || i == last) {
          segment = i == 0 ? 0 : last - 1;
        } else if (info.media[i - 1] == info.media[i]) {
          segment = i;  // both sides are the same material (e.g. a prism of one glass)
        } else {
          error(idx(location, k),
                "inner surface " + surfaces_[event.surface].id.str() + " of element '" +
                    surfaces_[event.surface].element_name +
                    "' reached from outside the element is ambiguous: the segments on its two "
                    "sides have different materials ('" +
                    media_[info.media[i - 1]].reference + "', '" + media_[info.media[i]].reference +
                    "'); enter a cemented group through its first or last surface (ADR 0017)");
          leave = true;
        }
        if (leave) {
          inside.reset();
          current = environment_;
        } else {
          inside = element;
          current = info.media[segment];
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
  /// Valid wavelength range and first referencing JSON pointer, per entry of media_.
  struct MediumCheck {
    std::optional<material::WavelengthRange> range;
    std::string location;
  };
  std::vector<MediumCheck> media_checks_;
  std::set<std::string> mangin_reported_;  // mirrors already reported as Mangin mirrors
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
    cs.wavelength_weights_.push_back(system.wavelengths[i].weight);
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
