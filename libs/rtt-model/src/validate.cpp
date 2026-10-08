#include "rtt/model/validate.hpp"

#include <algorithm>
#include <cmath>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>

namespace rtt::model {

namespace {

std::string idx(std::string_view base, std::size_t i) {
  return std::string(base) + "/" + std::to_string(i);
}

bool finite_positive(double v) {
  return std::isfinite(v) && v > 0.0;
}

bool nonzero_axis(const std::array<double, 3>& a) {
  return std::isfinite(a[0]) && std::isfinite(a[1]) && std::isfinite(a[2]) &&
         (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]) > 1e-24;
}

class Validator {
 public:
  explicit Validator(const System& system) : system_(system) {}

  std::vector<Diagnostic> run() {
    check_wavelengths();
    check_aperture();
    check_fields();
    check_object_and_environment();
    check_assembly(system_.root, "/root");
    check_stops();
    check_paths();
    return std::move(out_);
  }

 private:
  /// Adds a diagnostic with the severity of `code` in the registry.
  void report(diagnostics::DiagnosticCode code, std::string location, std::string message) {
    out_.push_back(
        {code.severity(), std::move(location), std::move(message), std::string(code.str())});
  }

  void check_wavelengths() {
    if (system_.wavelengths.empty()) {
      report("wavelengths.empty", "/wavelengths", "at least one wavelength is required");
      return;
    }
    std::size_t references = 0;
    for (std::size_t i = 0; i < system_.wavelengths.size(); ++i) {
      const Wavelength& w = system_.wavelengths[i];
      const std::string loc = idx("/wavelengths", i);
      if (!finite_positive(w.um))
        report("wavelengths.value_invalid", loc + "/um", "wavelength must be finite and > 0 um");
      if (!std::isfinite(w.weight) || w.weight < 0.0)
        report("wavelengths.weight_invalid", loc + "/weight", "weight must be >= 0");
      if (w.reference) ++references;
    }
    if (references != 1) {
      report("wavelengths.reference_count", "/wavelengths",
             "exactly one wavelength must be the reference, found " + std::to_string(references));
    }
  }

  void check_aperture() {
    if (system_.aperture.type != SystemApertureType::StopSize &&
        !finite_positive(system_.aperture.value.value)) {
      report("aperture.value_invalid", "/aperture/value", "aperture value must be finite and > 0");
    }
    if (system_.aperture.type == SystemApertureType::ObjectSpaceNA &&
        system_.aperture.value.value >= 1.0 && system_.environment.medium == "AIR") {
      report("aperture.na_not_physical", "/aperture/value",
             "object-space NA >= 1 in air is not physical");
    }
  }

  void check_fields() {
    if (system_.fields.points.empty()) {
      report("fields.empty", "/fields/points", "at least one field point is required");
    }
    for (std::size_t i = 0; i < system_.fields.points.size(); ++i) {
      const Field& f = system_.fields.points[i];
      const std::string loc = idx("/fields/points", i);
      if (!std::isfinite(f.x) || !std::isfinite(f.y))
        report("fields.coordinate_invalid", loc, "field coordinates must be finite");
      if (!std::isfinite(f.weight) || f.weight < 0.0)
        report("fields.weight_invalid", loc + "/weight", "weight must be >= 0");
    }
  }

  void check_object_and_environment() {
    if (!system_.object.at_infinity && !finite_positive(system_.object.distance.value)) {
      report("object.distance_invalid", "/object/distance",
             "finite object distance must be > 0 mm");
    }
    const Environment& e = system_.environment;
    if (!std::isfinite(e.temperature_c) || e.temperature_c <= -273.15) {
      report("environment.temperature_invalid", "/environment/temperature_c",
             "temperature must be above absolute zero");
    }
    if (!std::isfinite(e.pressure_atm) || e.pressure_atm < 0.0) {
      report("environment.pressure_invalid", "/environment/pressure_atm",
             "pressure must be >= 0 atm");
    }
    if (e.medium.empty())
      report("environment.medium_empty", "/environment/medium", "medium must not be empty");
  }

  void check_assembly(const Assembly& a, const std::string& loc) {
    if (a.name.empty()) report("node.name_empty", loc + "/name", "assembly name must not be empty");
    check_name(a.name, loc);
    for (std::size_t i = 0; i < a.children.size(); ++i) {
      const std::string child_loc = idx(loc + "/children", i);
      const auto& v = a.children[i].value;
      if (const auto* sub = std::get_if<Assembly>(&v)) {
        check_assembly(*sub, child_loc);
      } else {
        check_element(std::get<Element>(v), child_loc);
      }
    }
  }

  void check_name(const std::string& name, const std::string& loc) {
    if (name.empty()) return;
    auto [it, inserted] = node_names_.emplace(name, loc);
    if (!inserted)
      report("node.name_duplicate", loc + "/name",
             "duplicate node name '" + name + "' (first at " + it->second + ")");
  }

  void check_element(const Element& e, const std::string& loc) {
    if (e.name.empty()) report("node.name_empty", loc + "/name", "element name must not be empty");
    check_name(e.name, loc);
    const std::size_t n = e.surfaces.size();
    switch (e.kind) {
      case ElementKind::Lens:
        if (n < 2)
          report("element.surface_count", loc + "/surfaces", "a lens needs at least 2 surfaces");
        check_segment_materials(e, loc, "a lens");
        break;
      case ElementKind::Plate:
        if (n < 2)
          report("element.surface_count", loc + "/surfaces", "a plate needs at least 2 surfaces");
        check_segment_materials(e, loc, "a plate");
        for (std::size_t i = 0; i < n; ++i) {
          if (!std::holds_alternative<Plane>(e.surfaces[i].shape.base)) {
            report("element.plate_surface_not_plane", idx(loc + "/surfaces", i) + "/shape",
                   "plate surfaces must be planes");
          }
        }
        break;
      case ElementKind::Mirror:
        if (n < 1)
          report("element.surface_count", loc + "/surfaces", "a mirror needs at least 1 surface");
        if (!e.segment_materials.empty()) {
          report("element.material_list_not_allowed", loc + "/material",
                 "a mirror takes a single substrate material, not a list");
        }
        break;
      case ElementKind::ThinElement:
      case ElementKind::Stop:
      case ElementKind::Detector:
        if (n != 1)
          report("element.surface_count", loc + "/surfaces",
                 "this element kind needs exactly 1 surface");
        if (e.material || !e.segment_materials.empty()) {
          report("element.material_not_allowed", loc + "/material",
                 "this element kind has no material");
        }
        break;
    }
    if (e.kind == ElementKind::Stop) {
      stops_.push_back(loc);
      if (n == 1 && !e.surfaces[0].aperture) {
        report("stop.aperture_missing", loc + "/surfaces/0/aperture", "a stop needs an aperture");
      }
    }
    if (e.material && e.material->empty())
      report("element.material_empty", loc + "/material", "material must not be empty");
    check_crystal(e, loc);
    for (std::size_t i = 0; i < n; ++i) check_surface(e.surfaces[i], idx(loc + "/surfaces", i));
  }

  /// Crystal material and optic axis (ADR 0026, points 1 and 2).
  void check_crystal(const Element& e, const std::string& loc) {
    const std::string mloc = loc + "/material";
    const std::string aloc = loc + "/optic_axis";
    if (e.crystal) {
      if (e.kind != ElementKind::Lens && e.kind != ElementKind::Plate)
        report("crystal.kind_not_allowed", mloc, "only a lens or a plate can be a crystal");
      if (e.material || !e.segment_materials.empty()) {
        report("crystal.material_conflict", mloc,
               "crystal and isotropic material are both set; use one of them");
      }
      if (e.crystal->ordinary.empty())
        report("element.material_empty", mloc + "/ordinary", "material must not be empty");
      if (e.crystal->extraordinary.empty())
        report("element.material_empty", mloc + "/extraordinary", "material must not be empty");
      if (!e.optic_axis)
        report("crystal.optic_axis_missing", aloc, "a crystal needs an optic axis");
    } else if (e.optic_axis) {
      report("crystal.optic_axis_not_allowed", aloc, "optic axis at an element without crystal");
    }
    if (e.optic_axis && !nonzero_axis(*e.optic_axis))
      report("crystal.optic_axis_invalid", aloc, "optic axis must be finite and non-zero");
  }

  /// Lens and Plate: N surfaces have N - 1 segments, given either as one shorthand material
  /// for all segments or as a list with one entry per segment (ADR 0017).
  void check_segment_materials(const Element& e, const std::string& loc, std::string_view what) {
    if (e.crystal) return;  // check_crystal()
    const std::string mloc = loc + "/material";
    const std::size_t n = e.surfaces.size();
    const auto& list = e.segment_materials;
    if (e.material && !list.empty()) {
      report("element.material_both", mloc,
             "use either one material for all segments or a list, not both");
      return;
    }
    if (!e.material && list.empty()) {
      report("element.material_missing", mloc, std::string(what) + " needs a material");
      return;
    }
    if (list.empty()) return;
    const std::size_t segments = n < 2 ? 0 : n - 1;
    if (list.size() != segments) {
      report("element.segment_count", mloc,
             std::string(what) + " with " + std::to_string(n) + " surfaces needs " +
                 std::to_string(segments) + " segment materials, found " +
                 std::to_string(list.size()));
    }
    for (std::size_t i = 0; i < list.size(); ++i) {
      if (list[i].empty())
        report("element.material_empty", idx(mloc, i), "material must not be empty");
    }
  }

  void check_surface(const Surface& s, const std::string& loc) {
    if (s.id.empty()) {
      report("surface.id_empty", loc + "/id", "surface id must not be empty");
    } else {
      auto [it, inserted] = surface_ids_.emplace(s.id.str(), loc);
      if (!inserted) {
        report("surface.id_duplicate", loc + "/id",
               "duplicate surface id '" + s.id.str() + "' (first at " + it->second + ")");
      }
    }
    check_shape(s.shape, loc + "/shape");
    if (s.aperture) check_aperture(*s.aperture, loc + "/aperture");
    for (std::size_t i = 0; i < s.phases.size(); ++i) {
      check_phase(s.phases[i], idx(loc + "/phases", i));
    }
    if (!s.phases.empty()) phase_surfaces_.insert(s.id.str());
    if (s.diffraction_efficiency) check_efficiency(s, loc + "/diffraction_efficiency");
    check_interaction(s.interaction, loc + "/interaction");
  }

  /// Efficiency per diffraction order (ADR 0025, point 5).
  void check_efficiency(const Surface& s, const std::string& loc) {
    const auto& list = *s.diffraction_efficiency;
    if (s.phases.empty()) {
      report("surface.efficiency_invalid", loc,
             "diffraction efficiency at a surface without phase layer");
    }
    if (list.empty()) {
      report("surface.efficiency_invalid", loc,
             "empty list of diffraction efficiencies (it would block every order)");
    }
    std::unordered_set<int> orders;
    for (std::size_t i = 0; i < list.size(); ++i) {
      const DiffractionEfficiency& d = list[i];
      if (!std::isfinite(d.efficiency) || d.efficiency < 0.0 || d.efficiency > 1.0) {
        report("surface.efficiency_invalid", idx(loc, i) + "/efficiency",
               "efficiency must be in [0, 1]");
      }
      if (!orders.insert(d.order).second) {
        report("surface.efficiency_invalid", idx(loc, i) + "/order",
               "order " + std::to_string(d.order) + " is listed twice");
      }
    }
  }

  void check_radius(const Param& r, const std::string& loc) {
    if (!std::isfinite(r.value) || r.value == 0.0) {
      report("shape.radius_invalid", loc,
             "radius must be finite and non-zero (use a plane for flat surfaces)");
    }
  }

  void check_shape(const ShapeStack& shape, const std::string& loc) {
    const std::string base = loc + "/base";
    if (const auto* c = std::get_if<Conic>(&shape.base)) {
      check_radius(c->radius, base + "/radius");
    } else if (const auto* a = std::get_if<EvenAsphere>(&shape.base)) {
      check_radius(a->radius, base + "/radius");
      if (a->coefficients.empty())
        report("shape.asphere_without_coefficients", base + "/coefficients",
               "asphere without coefficients");
    }
    for (std::size_t i = 0; i < shape.terms.size(); ++i) {
      const auto& z = std::get<ZernikeSag>(shape.terms[i]);
      const std::string t = idx(loc + "/terms", i);
      if (!finite_positive(z.normalization_radius.value)) {
        report("shape.zernike_radius_invalid", t + "/normalization_radius",
               "normalization radius must be > 0 mm");
      }
      if (z.coefficients.empty())
        report("shape.zernike_without_coefficients", t + "/coefficients",
               "Zernike term without coefficients");
    }
  }

  void check_aperture(const Aperture& a, const std::string& loc) {
    if (const auto* c = std::get_if<CircularAperture>(&a)) {
      if (!finite_positive(c->radius))
        report("surface_aperture.radius_invalid", loc + "/radius", "radius must be > 0 mm");
      if (!std::isfinite(c->inner_radius) || c->inner_radius < 0.0 ||
          c->inner_radius >= c->radius) {
        report("surface_aperture.inner_radius_invalid", loc + "/inner_radius",
               "inner radius must be >= 0 and < radius");
      }
    } else if (const auto* r = std::get_if<RectangularAperture>(&a)) {
      if (!finite_positive(r->half_width_x) || !finite_positive(r->half_width_y)) {
        report("surface_aperture.half_width_invalid", loc, "half widths must be > 0 mm");
      }
    } else if (const auto* e = std::get_if<EllipticalAperture>(&a)) {
      if (!finite_positive(e->semi_axis_x) || !finite_positive(e->semi_axis_y)) {
        report("surface_aperture.semi_axis_invalid", loc, "semi axes must be > 0 mm");
      }
    }
  }

  void check_phase(const PhaseLayer& p, const std::string& loc) {
    if (const auto* g = std::get_if<LinearGrating>(&p)) {
      if (!finite_positive(g->lines_per_mm.value))
        report("phase.lines_per_mm_invalid", loc + "/lines_per_mm", "must be > 0");
    } else if (const auto* r = std::get_if<RadialPhase>(&p)) {
      if (!finite_positive(r->normalization_radius.value)) {
        report("phase.radius_invalid", loc + "/normalization_radius",
               "normalization radius must be > 0 mm");
      }
    }
  }

  void check_interaction(const Interaction& i, const std::string& loc) {
    if (const auto* bs = std::get_if<IdealBeamSplitter>(&i)) {
      for (const double r : {bs->reflectance_s, bs->reflectance_p}) {
        if (!std::isfinite(r) || r < 0.0 || r > 1.0)
          report("interaction.reflectance_invalid", loc, "reflectance must be in [0, 1]");
      }
    } else if (const auto* c = std::get_if<CoatingRef>(&i)) {
      if (c->name.empty())
        report("interaction.coating_name_empty", loc + "/name", "coating name must not be empty");
    } else if (const auto* pol = std::get_if<IdealPolarizer>(&i)) {
      if (!nonzero_axis(pol->transmission_axis))
        report("interaction.axis_invalid", loc + "/transmission_axis", "axis must be non-zero");
      if (!std::isfinite(pol->extinction_ratio) || pol->extinction_ratio < 0.0 ||
          pol->extinction_ratio > 1.0) {
        report("interaction.extinction_ratio_invalid", loc + "/extinction_ratio",
               "extinction ratio must be in [0, 1]");
      }
    } else if (const auto* ret = std::get_if<IdealRetarder>(&i)) {
      if (!nonzero_axis(ret->fast_axis))
        report("interaction.axis_invalid", loc + "/fast_axis", "axis must be non-zero");
      if (!std::isfinite(ret->retardance_waves))
        report("interaction.retardance_invalid", loc + "/retardance_waves", "must be finite");
    }
  }

  void check_stops() {
    if (stops_.size() > 1) report("stop.multiple", stops_[1], "only one stop element is allowed");
    if (stops_.empty() && system_.aperture.type == SystemApertureType::StopSize) {
      report("aperture.stop_missing", "/aperture/type",
             "aperture type 'stop_size' requires a stop element");
    }
  }

  void check_paths() {
    if (system_.paths.empty()) {
      report("paths.empty", "/paths", "at least one path is required");
      return;
    }
    std::unordered_set<std::string> names;
    for (std::size_t i = 0; i < system_.paths.size(); ++i) {
      const Path& p = system_.paths[i];
      const std::string loc = idx("/paths", i);
      if (p.name.empty()) report("paths.name_empty", loc + "/name", "path name must not be empty");
      if (!names.insert(p.name).second)
        report("paths.name_duplicate", loc + "/name", "duplicate path name '" + p.name + "'");
      if (!p.automatic && p.events.empty())
        report("paths.events_empty", loc + "/events", "explicit path without events");
      for (std::size_t k = 0; k < p.events.size(); ++k) {
        const Event& e = p.events[k];
        const std::string eloc = idx(loc + "/events", k);
        if (!surface_ids_.contains(e.surface.str())) {
          report("paths.unknown_surface", eloc + "/surface",
                 "unknown surface id '" + e.surface.str() + "'");
        }
        // ADR 0025: an event diffracts exactly when its surface has a phase layer.
        if (e.order != 0 && surface_ids_.contains(e.surface.str()) &&
            !phase_surfaces_.contains(e.surface.str())) {
          report("paths.order_not_allowed", eloc + "/order",
                 "diffraction order at a surface without phase layer");
        }
      }
    }
  }

  const System& system_;
  std::vector<Diagnostic> out_;
  std::unordered_map<std::string, std::string> surface_ids_;
  std::unordered_set<std::string> phase_surfaces_;  // ids of surfaces with a phase layer
  std::unordered_map<std::string, std::string> node_names_;
  std::vector<std::string> stops_;
};

}  // namespace

std::vector<Diagnostic> validate(const System& system) {
  return Validator(system).run();
}

bool has_errors(const std::vector<Diagnostic>& diagnostics) {
  return std::any_of(diagnostics.begin(), diagnostics.end(),
                     [](const Diagnostic& d) { return d.severity == Severity::Error; });
}

std::string to_string(const Diagnostic& diagnostic) {
  std::string text = diagnostic.severity == Severity::Error ? "error " : "warning ";
  if (!diagnostic.code.empty()) text += "[" + diagnostic.code + "] ";
  return text + diagnostic.location + ": " + diagnostic.message;
}

}  // namespace rtt::model
