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
#include <stdexcept>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

#include "rtt/coating/catalog.hpp"
#include "rtt/coating/thickness.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/compile/errors.hpp"
#include "rtt/diagnostics/codes.hpp"

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
  /// The media follow the segments: every Lens, and a Plate with different segment materials.
  /// Otherwise (Mirror, Plate of one material) Refract toggles inside <-> environment.
  bool segmented = false;
  std::string location;             ///< JSON pointer of the element
  std::uint32_t first_surface = 0;  ///< surfaces [first_surface, first_surface + count)
  std::uint32_t surface_count = 0;
};

/// Collects the pieces of a CompiledSystem and all errors found on the way.
class Compiler {
 public:
  Compiler(const model::System& system,
           const material::MaterialLibrary& materials,
           const coating::CoatingLibrary& coatings)
      : system_(system), materials_(materials), coatings_(coatings) {
    for (const auto& w : system.wavelengths) wavelengths_um_.push_back(w.um);
  }

  void run() {
    if (wavelengths_um_.size() > std::numeric_limits<std::uint16_t>::max()) {
      report("wavelengths.too_many", "/wavelengths",
             "at most 65535 wavelengths are supported (RayBatch::wl is 16 bit)");
    }
    // The environment is always media_[0] (CompiledSystem::environment_medium()).
    environment_ = medium(system_.environment.medium, "/environment/medium").value_or(0);
    add_assembly(system_.root, model::to_isometry(system_.root.pose), "/root");
    for (std::size_t p = 0; p < system_.paths.size(); ++p) {
      paths_.push_back(build_path(system_.paths[p], idx("/paths", p)));
    }
    check_stop_on_paths();
    check_wavelength_ranges();
  }

  std::vector<model::Diagnostic> errors_;
  std::vector<model::Diagnostic> warnings_;
  std::vector<double> wavelengths_um_;
  std::vector<CompiledSurface> surfaces_;
  std::vector<CompiledElement>
      compiled_elements_;  ///< in tree order, see CompiledSystem::elements()
  std::vector<CompiledMedium> media_;
  std::vector<CompiledCoating> compiled_coatings_;
  std::vector<CompiledPath> paths_;

 private:
  /// Adds a diagnostic with the severity of `code` in the registry: errors make compile() throw,
  /// warnings end up in CompiledSystem::diagnostics().
  void report(diagnostics::DiagnosticCode code, std::string location, std::string message) {
    auto& list = code.severity() == model::Severity::Error ? errors_ : warnings_;
    list.push_back(
        {code.severity(), std::move(location), std::move(message), std::string(code.str())});
  }

  /// Warning for every path that does not visit the stop element of a system that has one
  /// (ADR 0023): aiming, pupils, Seidel sums and OPD on such a path throw NoStopError.
  void check_stop_on_paths() {
    const auto is_stop = [&](std::uint32_t surface) {
      return surfaces_[surface].element_kind == model::ElementKind::Stop;
    };
    const auto stop = std::find_if(
        surfaces_.begin(), surfaces_.end(),
        [](const CompiledSurface& s) { return s.element_kind == model::ElementKind::Stop; });
    if (stop == surfaces_.end()) return;
    for (std::size_t p = 0; p < paths_.size(); ++p) {
      const auto& events = paths_[p].events;
      if (std::any_of(events.begin(), events.end(),
                      [&](const CompiledEvent& e) { return is_stop(e.surface); })) {
        continue;
      }
      report("stop.not_on_path", idx("/paths", p),
             "path '" + paths_[p].name + "' does not visit the stop '" + stop->element_name +
                 "': aiming, pupils, Seidel sums and OPD on it need a stop (NoStopError)");
    }
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
      report("material.unknown", location, e.what());
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
    check_coating_ranges();
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
        report("material.wavelength_out_of_range", location.value(),
               "wavelength " + number(wl) + " um is outside the valid range [" +
                   number(valid.min_um) + ", " + number(valid.max_um) + "] um of material '" +
                   media_[m].reference + "'");
        break;
      }
    }
  }

  /// Every system wavelength must lie in the valid range of every layer material of a coating
  /// on a surface that a path uses (as for media, #23). One error per coating and layer, at
  /// the interaction of the first such surface.
  void check_coating_ranges() {
    std::vector<std::optional<std::string>> where(compiled_coatings_.size());
    for (const CompiledPath& path : paths_) {
      for (const CompiledEvent& event : path.events) {
        const auto& coating = surfaces_[event.surface].coating;
        if (coating && !where[coating->coating]) {
          where[coating->coating] = surface_locations_[event.surface] + "/interaction/name";
        }
      }
    }
    for (const CoatingCheck& check : coating_checks_) {
      const std::optional<std::string>& location = where[check.coating];
      if (!location.has_value() || !check.range.has_value()) continue;
      const material::WavelengthRange& valid = check.range.value();
      for (const double wl : wavelengths_um_) {
        if (valid.contains(wl)) continue;
        report("coating.wavelength_out_of_range", location.value(),
               "wavelength " + number(wl) + " um is outside the valid range [" +
                   number(valid.min_um) + ", " + number(valid.max_um) + "] um of material '" +
                   check.material + "' in layer " + std::to_string(check.layer) + " of coating '" +
                   compiled_coatings_[check.coating].reference + "'");
        break;
      }
    }
  }

  /// Substrate side and compiled design of a CoatingRef at surface j of an element (ADR 0019).
  /// `location` is the JSON pointer of the interaction.
  std::optional<SurfaceCoating> surface_coating(const model::CoatingRef& ref,
                                                const ElementInfo& info,
                                                std::size_t j,
                                                const std::string& location) {
    std::optional<std::uint32_t> substrate;
    const std::size_t last = info.surface_count - 1;
    switch (info.kind) {
      case model::ElementKind::Lens:
      case model::ElementKind::Plate:
        if (info.media.empty()) {
          report("coating.no_substrate", location,
                 "coating on an element without material: no substrate (ADR 0019)");
        } else if (!info.segmented || j == 0) {
          // Every face of a plate of one material bounds its inside; the first surface of a
          // segmented element bounds segment 0.
          substrate = info.media[0];
        } else if (j == last) {
          substrate = info.media[last - 1];
        } else {
          report("coating.substrate_ambiguous", location,
                 "coating on an inner surface between two segments: the substrate side is "
                 "ambiguous (ADR 0019)");
        }
        break;
      case model::ElementKind::Mirror:
        if (info.media.empty()) {
          report("coating.no_substrate", location,
                 "coating on a mirror without substrate material: no substrate "
                 "(ADR 0019; use ideal_mirror for a mirror without substrate)");
        } else {
          substrate = info.media[0];
        }
        break;
      case model::ElementKind::ThinElement:
      case model::ElementKind::Stop:
      case model::ElementKind::Detector:
        report("coating.not_allowed", location,
               "coating on a surface of a thin element, stop or detector: only lens, plate and "
               "mirror surfaces have a substrate (ADR 0019)");
        break;
    }
    const std::optional<std::uint32_t> coating = compiled_coating(ref.name, location + "/name");
    if (!substrate || !coating) return std::nullopt;
    return SurfaceCoating{*coating, *substrate};
  }

  /// Index of the compiled coating for `reference`, resolving and evaluating it on first use.
  std::optional<std::uint32_t> compiled_coating(const std::string& reference,
                                                const std::string& location) {
    if (const auto it = coating_index_.find(reference); it != coating_index_.end()) {
      return it->second;
    }
    std::shared_ptr<const coating::CoatingDesign> design;
    try {
      design = coatings_.resolve(reference);
    } catch (const coating::UnknownCoating& e) {
      report("coating.unknown", location, e.what());
      coating_index_.emplace(reference, std::nullopt);  // report once
      return std::nullopt;
    }
    const auto index = static_cast<std::uint32_t>(compiled_coatings_.size());
    const double t = system_.environment.temperature_c;
    const double p = system_.environment.pressure_atm;
    CompiledCoating compiled{
        reference, std::vector<std::vector<coating::Layer<double>>>(wavelengths_um_.size())};
    bool ok = true;
    std::vector<CoatingCheck> checks;  // appended to coating_checks_ only on success
    for (std::size_t k = 0; k < design->layers.size(); ++k) {
      const coating::LayerSpec& spec = design->layers[k];
      const std::string what = "layer " + std::to_string(k) + " of coating '" + reference + "'";
      std::shared_ptr<const material::Material> material;
      try {
        material = materials_.resolve(spec.material);
      } catch (const material::UnknownMaterial& e) {
        report("coating.layer_material_unknown", location, what + ": " + e.what());
        ok = false;
        continue;
      }
      // QWOT: d = count lambda0 / (4 Re n(lambda0)) at the environment temperature and pressure
      // (thickness.hpp, Byrnes Eq. (8) with (2)); the design wavelength must lie in the range of
      // the material. physical_thickness_um() also rejects Re n <= 0 or not finite.
      double index_at_design = 1.0;  // unused for a physical thickness
      if (const auto* qwot = std::get_if<coating::QuarterWaves>(&spec.thickness)) {
        const double l0 = qwot->design_wavelength_um;
        const auto range = material->wavelength_range_um();
        if (range && !range->contains(l0)) {
          report("coating.design_wavelength_out_of_range", location,
                 what + ": design wavelength " + number(l0) +
                     " um is outside the valid range of material '" + spec.material + "'");
          ok = false;
          continue;
        }
        index_at_design = material->index(l0, t, p).real();
      }
      double thickness_um = 0.0;
      try {
        thickness_um = coating::physical_thickness_um(spec.thickness, index_at_design);
      } catch (const std::invalid_argument& e) {
        report("coating.thickness_invalid", location, what + ": " + e.what());
        ok = false;
        continue;
      }
      for (std::size_t w = 0; w < wavelengths_um_.size(); ++w) {
        compiled.layers[w].push_back({material->index(wavelengths_um_[w], t, p), thickness_um});
      }
      checks.push_back({index, k, spec.material, material->wavelength_range_um()});
    }
    if (!ok) {
      coating_index_.emplace(reference, std::nullopt);
      return std::nullopt;
    }
    compiled_coatings_.push_back(std::move(compiled));
    coating_checks_.insert(coating_checks_.end(), checks.begin(), checks.end());
    coating_index_.emplace(reference, index);
    return index;
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
    const bool lens_or_plate =
        element.kind == model::ElementKind::Lens || element.kind == model::ElementKind::Plate;
    if (lens_or_plate && !element.segment_materials.empty()) {
      for (std::size_t i = 0; i < element.segment_materials.size(); ++i) {
        info.media_locations.push_back(idx(location + "/material", i));
        info.media.push_back(
            medium(element.segment_materials[i], info.media_locations.back()).value_or(0));
      }
    } else if (element.material) {
      // The shorthand holds for all segments; resolved once, reported once.
      const std::uint32_t m = medium(*element.material, location + "/material").value_or(0);
      const std::size_t count = lens_or_plate ? element.surfaces.size() - 1 : 1;
      info.media.assign(count, m);
      info.media_locations.assign(count, location + "/material");
    }
    // Equal references share one medium index, so equal indices mean one material (#27).
    const bool uniform = std::adjacent_find(info.media.begin(), info.media.end(),
                                            std::not_equal_to<>()) == info.media.end();
    info.segmented = element.kind == model::ElementKind::Lens ||
                     (element.kind == model::ElementKind::Plate && !uniform);
    const auto element_index = static_cast<std::uint32_t>(elements_.size());
    elements_.push_back(info);

    for (std::size_t j = 0; j < element.surfaces.size(); ++j) {
      const model::Surface& s = element.surfaces[j];
      const std::string surface_location = idx(location + "/surfaces", j);
      CompiledSurface c;
      c.id = s.id;
      c.element_kind = element.kind;
      c.element_name = element.name;
      c.element = element_index;
      c.location = surface_location;
      c.to_global = to_global * model::to_isometry(s.pose);
      c.to_local = c.to_global.inverse();
      c.shape = compile_shape(s.shape, surface_location + "/shape");
      c.aperture = s.aperture;
      c.phases = s.phases;
      c.interaction = s.interaction;
      // Axes of ideal elements are given in element coordinates (ADR 0021).
      if (const auto* polarizer = std::get_if<model::IdealPolarizer>(&s.interaction)) {
        const auto& a = polarizer->transmission_axis;
        c.ideal_axis = to_global.apply_vector(math::Vec3(a[0], a[1], a[2]));
      } else if (const auto* retarder = std::get_if<model::IdealRetarder>(&s.interaction)) {
        const auto& a = retarder->fast_axis;
        c.ideal_axis = to_global.apply_vector(math::Vec3(a[0], a[1], a[2]));
      }
      if (const auto* ref = std::get_if<model::CoatingRef>(&s.interaction)) {
        c.coating = surface_coating(*ref, info, j, surface_location + "/interaction");
      }
      surface_index_.emplace(s.id, static_cast<std::uint32_t>(surfaces_.size()));
      surface_locations_.push_back(surface_location);
      surface_element_.push_back(element_index);
      surfaces_.push_back(std::move(c));
    }

    // Media in front of and behind each surface in the element's surface order (#81): a
    // refraction through the surfaces in order, from the environment, with the rules of cross()
    // (ADR 0017). Independent of the paths.
    MediumState state{std::nullopt, 0, environment_};
    for (std::uint32_t j = 0; j < info.surface_count; ++j) {
      CompiledSurface& c = surfaces_[info.first_surface + j];
      c.medium_front = state.current;
      // In order from outside through the first surface, cross() is never ambiguous.
      if (const std::optional<MediumState> crossed = cross(info.first_surface + j, state)) {
        state = *crossed;
      }
      c.medium_back = state.current;
    }
    compiled_elements_.push_back(CompiledElement{element.name, element.kind, info.first_surface,
                                                 info.surface_count, info.media, info.segmented});
  }

  CompiledShape compile_shape(const model::ShapeStack& shape, const std::string& location) {
    for (std::size_t k = 0; k < shape.terms.size(); ++k) {
      report("shape.zernike_unsupported", idx(location + "/terms", k),
             "Zernike sag terms are not supported before M8");
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
            report("paths.mangin_mirror_automatic", e.location + "/surfaces",
                   "mirror with substrate material on the automatic path: use an explicit path "
                   "(Refract, Reflect, Refract)");
          }
          continue;
        }
        // A plate of one material with more than 2 surfaces is a prism or cube: which faces
        // the ray uses depends on the design, so the tree order is no path (#27).
        if (e.kind == model::ElementKind::Plate && !e.media.empty() && !e.segmented &&
            e.surface_count > 2) {
          if (uniform_plate_reported_.emplace(e.location).second) {
            report("paths.uniform_plate_automatic", e.location + "/surfaces",
                   "plate with more than 2 surfaces and uniform material needs an explicit path "
                   "(prism or cube: the faces used depend on the design)");
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

  /// Where a ray is along a path: element it is in (none = environment), segment of that
  /// element, and the current medium.
  struct MediumState {
    std::optional<std::uint32_t> inside;
    std::uint32_t segment = 0;
    std::uint32_t current = 0;
  };

  /// State after crossing `surface` (Refract, Ordinary, Extraordinary) from `s` by the rules
  /// above; none if the side is ambiguous (inner surface between different materials reached
  /// from outside). An element without material does not change the state.
  [[nodiscard]] std::optional<MediumState> cross(std::uint32_t surface,
                                                 const MediumState& s) const {
    const std::uint32_t element = surface_element_[surface];
    const ElementInfo& info = elements_[element];
    if (info.media.empty()) return s;
    const std::uint32_t i = surface - info.first_surface;  // surface in the element
    const std::uint32_t last = info.surface_count - 1;
    std::uint32_t segment = s.segment;
    bool leave = false;
    if (!info.segmented) {
      leave = s.inside == element;
      segment = 0;
    } else if (s.inside == element && i == segment + 1) {  // forward through surface i
      leave = i == last;
      segment = i;
    } else if (s.inside == element && i == segment) {  // backward through surface i
      leave = i == 0;
      if (!leave) segment = i - 1;
    } else if (s.inside == element) {
      leave = true;  // the surface does not bound the segment
    } else if (i == 0 || i == last) {
      segment = i == 0 ? 0 : last - 1;
    } else if (info.media[i - 1] == info.media[i]) {
      segment = i;  // both sides are the same material
    } else {
      return std::nullopt;
    }
    if (leave) return MediumState{std::nullopt, segment, environment_};
    return MediumState{element, segment, info.media[segment]};
  }

  /// Media before and after each event (rules decided for #5 and #27, docs/architecture.md,
  /// "Medien entlang eines Pfads"). The ray starts in the environment. Reflect, Transmit and
  /// Diffract keep the medium, and so does every event at an element without material.
  /// Refract, Ordinary and Extraordinary
  /// - at surface i of a Lens, or of a Plate with different segment materials, with N surfaces
  ///   (segment j between surfaces j and j + 1):
  ///   from segment i - 1 into segment i (the environment if i = N - 1), from segment i into
  ///   segment i - 1 (the environment if i = 0), from any other segment of the same element into
  ///   the environment; from outside the element through the first surface into segment 0,
  ///   through the last surface into segment N - 2, and through an inner surface into its two
  ///   neighbouring segments if they have the same material, otherwise it is an error (the
  ///   side is ambiguous);
  /// - at a Plate of one material (shorthand, or a list of equal entries; prisms, cubes) and
  ///   at a Mirror with substrate (#6): toggle between the inside and the environment at every
  ///   surface. For a Plate with 2 surfaces this equals the segment rule.
  /// Outside an element means in the environment or in another element: elements do not nest,
  /// so entering B while in A leaves A, and leaving B goes to the environment.
  /// medium_beyond is the medium a crossing would reach, also for events that do not cross
  /// (#61); from_inside tells whether the ray is inside the surface's element before the event.
  /// `location` is the JSON pointer of the events; errors add the event index.
  void assign_media(std::vector<CompiledEvent>& events, const std::string& location) {
    MediumState state{std::nullopt, 0, environment_};
    for (std::size_t k = 0; k < events.size(); ++k) {
      CompiledEvent& event = events[k];
      event.medium_before = state.current;
      event.from_inside = state.inside == surface_element_[event.surface];
      const bool crosses = event.kind == model::EventKind::Refract ||
                           event.kind == model::EventKind::Ordinary ||
                           event.kind == model::EventKind::Extraordinary;
      const std::optional<MediumState> crossed = cross(event.surface, state);
      // The medium on the other side of the surface, also for events that do not cross it
      // (Fresnel and coatings on Reflect need it, #61); an ambiguous side counts as none.
      event.medium_beyond = crossed ? crossed->current : state.current;
      if (crosses) {
        if (crossed) {
          state = *crossed;
        } else {
          const std::uint32_t element = surface_element_[event.surface];
          const ElementInfo& info = elements_[element];
          const std::uint32_t i = event.surface - info.first_surface;
          report("paths.inner_surface_ambiguous", idx(location, k),
                 "inner surface " + surfaces_[event.surface].id.str() + " of element '" +
                     surfaces_[event.surface].element_name +
                     "' reached from outside the element is ambiguous: the segments on its two "
                     "sides have different materials ('" +
                     media_[info.media[i - 1]].reference + "', '" +
                     media_[info.media[i]].reference +
                     "'); enter a cemented group through its first or last surface (ADR 0017)");
          state = MediumState{std::nullopt, state.segment, environment_};
        }
      }
      event.medium_after = state.current;
    }
  }

  const model::System& system_;
  const material::MaterialLibrary& materials_;
  const coating::CoatingLibrary& coatings_;
  std::vector<std::string> surface_locations_;  // JSON pointer per surface
  /// Compiled coating per reference; none once it failed (reported once).
  std::map<std::string, std::optional<std::uint32_t>, std::less<>> coating_index_;
  /// Valid wavelength range of each layer material, checked for coatings on a path.
  struct CoatingCheck {
    std::uint32_t coating;
    std::size_t layer;
    std::string material;
    std::optional<material::WavelengthRange> range;
  };
  std::vector<CoatingCheck> coating_checks_;
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
  std::set<std::string> mangin_reported_;         // mirrors already reported as Mangin mirrors
  std::set<std::string> uniform_plate_reported_;  // plates already reported, see build_path
};

}  // namespace

CompileError::CompileError(std::vector<model::Diagnostic> diagnostics)
    : std::runtime_error(join_errors(diagnostics)), diagnostics_(std::move(diagnostics)) {}

NoStopError::NoStopError(std::string path_name, std::string location)
    : std::invalid_argument("the path '" + path_name + "' (" + location +
                            ") has no stop: aiming at the stop, the pupils, Seidel sums and OPD "
                            "need a Stop element on the path"),
      path_name_(std::move(path_name)),
      location_(std::move(location)) {}

void require_stop(const CompiledSystem& system, PathId path) {
  const CompiledPath& p = system.path(path);
  const bool has_stop = std::any_of(p.events.begin(), p.events.end(), [&](const CompiledEvent& e) {
    return system.surfaces()[e.surface].element_kind == model::ElementKind::Stop;
  });
  if (!has_stop) throw NoStopError(p.name, "/paths/" + std::to_string(path.index));
}

CompiledSystem compile(const model::System& system, const material::MaterialLibrary& materials) {
  const coating::CoatingLibrary none;
  return compile(system, materials, none);
}

CompiledSystem compile(const model::System& system,
                       const material::MaterialLibrary& materials,
                       const coating::CoatingLibrary& coatings) {
  std::vector<model::Diagnostic> diagnostics = model::validate(system);
  if (model::has_errors(diagnostics)) {
    std::erase_if(diagnostics,
                  [](const model::Diagnostic& d) { return d.severity != model::Severity::Error; });
    throw CompileError(std::move(diagnostics));
  }

  Compiler compiler(system, materials, coatings);
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
  cs.elements_ = std::move(compiler.compiled_elements_);
  cs.media_ = std::move(compiler.media_);
  cs.coatings_ = std::move(compiler.compiled_coatings_);
  cs.paths_ = std::move(compiler.paths_);
  cs.diagnostics_ = std::move(diagnostics);  // only warnings are left
  cs.diagnostics_.insert(cs.diagnostics_.end(), compiler.warnings_.begin(),
                         compiler.warnings_.end());
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
