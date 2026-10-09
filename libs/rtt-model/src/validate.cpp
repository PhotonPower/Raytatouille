#include "rtt/model/validate.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <variant>
#include <vector>

#include "rtt/model/parameters.hpp"

namespace rtt::model {

namespace {

std::string idx(std::string_view base, std::size_t i) {
  return std::string(base) + "/" + std::to_string(i);
}

bool finite_positive(double v) {
  return std::isfinite(v) && v > 0.0;
}

/// [A-Za-z_][A-Za-z0-9_]* in ASCII (ADR 0029, point 1).
bool valid_row_name(const std::string& name) {
  const auto letter = [](char c) {
    return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || c == '_';
  };
  const auto digit = [](char c) { return c >= '0' && c <= '9'; };
  if (name.empty() || !letter(name[0])) return false;
  return std::all_of(name.begin() + 1, name.end(), [&](char c) { return letter(c) || digit(c); });
}

/// Empty or only ASCII white space.
bool blank(const std::string& text) {
  return std::all_of(text.begin(), text.end(), [](char c) {
    return c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\f' || c == '\v';
  });
}

bool nonzero_axis(const std::array<double, 3>& a) {
  return std::isfinite(a[0]) && std::isfinite(a[1]) && std::isfinite(a[2]) &&
         (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]) > 1e-24;
}

class Validator {
 public:
  explicit Validator(const System& system) : system_(system) {}

  std::vector<Diagnostic> run() {
    check_configurations();
    check_parameters();
    // The evaluating codes of the table (#164): expressions and non-finite results.
    static_cast<void>(evaluate_parameters(system_, &out_));
    check_wavelengths();
    check_aperture();
    check_fields();
    check_object_and_environment();
    check_assembly(system_.root, "/root", false);
    check_stops();
    check_paths();
    check_optimization();
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
    const Param& v = system_.aperture.value;
    const bool checked = system_.aperture.type != SystemApertureType::StopSize;
    check_param(v, "/aperture/value", checked);
    if (v.is_bound()) return;  // its value comes from the table (ADR 0029, point 3)
    if (checked && !finite_positive(v.value)) {
      report("aperture.value_invalid", "/aperture/value", "aperture value must be finite and > 0");
    }
    if (system_.aperture.type == SystemApertureType::ObjectSpaceNA && v.value >= 1.0 &&
        system_.environment.medium == "AIR") {
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
    const Param& d = system_.object.distance;
    check_param(d, "/object/distance", !system_.object.at_infinity);
    if (!system_.object.at_infinity && !d.is_bound() && !finite_positive(d.value)) {
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

  /// `has_sibling`: the node is not the first child of its assembly (ADR 0028, point 5).
  void check_assembly(const Assembly& a, const std::string& loc, bool has_sibling) {
    if (a.name.empty()) report("node.name_empty", loc + "/name", "assembly name must not be empty");
    check_name(a.name, loc);
    check_node_pose(a.pose, loc, has_sibling);
    for (std::size_t i = 0; i < a.children.size(); ++i) {
      const std::string child_loc = idx(loc + "/children", i);
      const auto& v = a.children[i].value;
      if (const auto* sub = std::get_if<Assembly>(&v)) {
        check_assembly(*sub, child_loc, i > 0);
      } else {
        check_element(std::get<Element>(v), child_loc, i > 0);
      }
    }
  }

  /// Reference of an assembly or element pose (ADR 0028, point 5): the preceding surface in tree
  /// order (pre-order, so surfaces_seen_ counts the surfaces before the node) or the preceding
  /// sibling must exist.
  void check_node_pose(const Pose& pose, const std::string& loc, bool has_sibling) {
    const std::string rloc = loc + "/pose/reference";
    if (pose.reference == PoseReference::RelativeToPreceding && surfaces_seen_ == 0) {
      report("pose.no_preceding", rloc, "no surface precedes this node in tree order");
    }
    if (pose.reference == PoseReference::RelativeToSibling && !has_sibling) {
      report("pose.no_sibling", rloc, "the first node of an assembly has no preceding sibling");
    }
    check_pose_values(pose, loc + "/pose");
  }

  /// Params and pivot of a pose.
  void check_pose_values(const Pose& pose, const std::string& loc) {
    for (std::size_t k = 0; k < 3; ++k) {
      check_param(pose.position[k], idx(loc + "/position", k), false);
      check_param(pose.rotation_deg[k], idx(loc + "/rotation_deg", k), false);
      if (!std::isfinite(pose.pivot[k])) {
        report("value.not_finite", idx(loc + "/pivot", k), "pivot must be finite");
      }
    }
  }

  // ------------------------------------------------- parameter table (ADR 0029) -----

  void check_configurations() {
    std::unordered_set<std::string> names;
    for (std::size_t k = 0; k < system_.configurations.size(); ++k) {
      const std::string& name = system_.configurations[k].name;
      const std::string loc = idx("/configurations", k) + "/name";
      if (blank(name)) {
        report("configurations.name_invalid", loc, "configuration name must not be blank");
      } else if (!names.insert(name).second) {
        report("configurations.name_duplicate", loc, "duplicate configuration name '" + name + "'");
      }
    }
  }

  void check_parameters() {
    const std::size_t columns = std::max<std::size_t>(1, system_.configurations.size());
    for (std::size_t i = 0; i < system_.parameters.size(); ++i) {
      const ParameterRow& row = system_.parameters[i];
      const std::string loc = idx("/parameters", i);
      if (!valid_row_name(row.name)) {
        report("parameters.name_invalid", loc + "/name",
               "row name '" + row.name + "' is not of the form [A-Za-z_][A-Za-z0-9_]*");
      } else if (!row_names_.insert(row.name).second) {
        report("parameters.name_duplicate", loc + "/name", "duplicate row name '" + row.name + "'");
      }
      if (const auto* value = std::get_if<double>(&row.form)) {
        check_value(*value, loc + "/value");
        check_bounds(row.min, row.max, loc, {{*value, loc + "/value"}});
      } else if (const auto* values = std::get_if<std::vector<double>>(&row.form)) {
        if (values->size() != columns) {
          report("parameters.values_count", loc + "/values",
                 "expected " + std::to_string(columns) + " values (one per configuration), found " +
                     std::to_string(values->size()));
        }
        std::vector<std::pair<double, std::string>> located;
        for (std::size_t k = 0; k < values->size(); ++k) {
          check_value((*values)[k], idx(loc + "/values", k));
          located.emplace_back((*values)[k], idx(loc + "/values", k));
        }
        check_bounds(row.min, row.max, loc, located);
      } else {
        // A derived row follows its inputs (ADR 0029, point 4).
        if (row.variable) {
          report("parameters.variable_expression", loc + "/variable",
                 "a row with an expression cannot be variable");
        }
        if (row.min || row.max) {
          report("bounds.invalid", loc + (row.min ? "/min" : "/max"),
                 "bounds are only allowed at rows with value or values");
        }
      }
    }
  }

  void check_value(double v, const std::string& loc) {
    if (!std::isfinite(v)) report("value.not_finite", loc, "value must be finite");
  }

  /// min < max if both are set; values outside are a warning (a start value may lie outside).
  void check_bounds(const std::optional<double>& min,
                    const std::optional<double>& max,
                    const std::string& loc,
                    const std::vector<std::pair<double, std::string>>& values) {
    bool finite = true;
    if (min && !std::isfinite(*min)) {
      report("value.not_finite", loc + "/min", "bound must be finite");
      finite = false;
    }
    if (max && !std::isfinite(*max)) {
      report("value.not_finite", loc + "/max", "bound must be finite");
      finite = false;
    }
    if (!finite) return;
    if (min && max && !(*min < *max)) {
      report("bounds.invalid", loc + "/min", "min must be below max");
      return;
    }
    for (const auto& [v, at] : values) {
      if (std::isfinite(v) && ((min && v < *min) || (max && v > *max))) {
        report("bounds.value_outside", at, "value lies outside its bounds");
      }
    }
  }

  /// Every Param (ADR 0029, point 3): a bound one names an existing row and has neither
  /// variable nor bounds; an unbound one has a finite value (`checked`: the field has its own
  /// check of the value, which keeps its code) and valid bounds.
  void check_param(const Param& p, const std::string& loc, bool checked) {
    if (const std::optional<std::string>& name = p.param) {
      if (!row_names_.contains(*name)) {
        report("param.unknown_parameter", loc + "/param", "no parameter row '" + *name + "'");
      }
      if (p.variable || p.min || p.max) {
        report("param.bound_conflict", loc, "a bound Param has neither variable nor bounds");
      }
      return;
    }
    if (!checked) check_value(p.value, loc + "/value");
    check_bounds(p.min, p.max, loc, {{p.value, loc + "/value"}});
  }

  void check_name(const std::string& name, const std::string& loc) {
    if (name.empty()) return;
    auto [it, inserted] = node_names_.emplace(name, loc);
    if (!inserted)
      report("node.name_duplicate", loc + "/name",
             "duplicate node name '" + name + "' (first at " + it->second + ")");
  }

  void check_element(const Element& e, const std::string& loc, bool has_sibling) {
    if (e.name.empty()) report("node.name_empty", loc + "/name", "element name must not be empty");
    check_name(e.name, loc);
    check_node_pose(e.pose, loc, has_sibling);
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
    for (std::size_t i = 0; i < n; ++i) {
      const std::string sloc = idx(loc + "/surfaces", i);
      // The first surface defines the element's frame (ADR 0028, point 5); its reference would
      // lie outside the element. One diagnostic for both kinds.
      if (i == 0 && e.surfaces[0].pose.reference != PoseReference::Absolute) {
        report("pose.relative_first_surface", sloc + "/pose/reference",
               "the first surface of an element must be placed absolutely");
      }
      check_pose_values(e.surfaces[i].pose, sloc + "/pose");
      check_surface(e.surfaces[i], sloc);
    }
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
    ++surfaces_seen_;
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
    if (const auto& efficiency = s.diffraction_efficiency) {
      check_efficiency(s, *efficiency, loc + "/diffraction_efficiency");
    }
    check_interaction(s.interaction, loc + "/interaction");
  }

  /// Efficiency per diffraction order (ADR 0025, point 5); `list` is the surface's list.
  void check_efficiency(const Surface& s,
                        const std::vector<DiffractionEfficiency>& list,
                        const std::string& loc) {
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
    check_param(r, loc, true);
    if (r.is_bound()) return;  // its value comes from the table
    if (!std::isfinite(r.value) || r.value == 0.0) {
      report("shape.radius_invalid", loc,
             "radius must be finite and non-zero (use a plane for flat surfaces)");
    }
  }

  void check_shape(const ShapeStack& shape, const std::string& loc) {
    const std::string base = loc + "/base";
    if (const auto* c = std::get_if<Conic>(&shape.base)) {
      check_radius(c->radius, base + "/radius");
      check_param(c->conic, base + "/conic", false);
    } else if (const auto* a = std::get_if<EvenAsphere>(&shape.base)) {
      check_radius(a->radius, base + "/radius");
      check_param(a->conic, base + "/conic", false);
      for (std::size_t k = 0; k < a->coefficients.size(); ++k) {
        check_param(a->coefficients[k], idx(base + "/coefficients", k), false);
      }
      if (a->coefficients.empty())
        report("shape.asphere_without_coefficients", base + "/coefficients",
               "asphere without coefficients");
    }
    for (std::size_t i = 0; i < shape.terms.size(); ++i) {
      const auto& z = std::get<ZernikeSag>(shape.terms[i]);
      const std::string t = idx(loc + "/terms", i);
      check_param(z.normalization_radius, t + "/normalization_radius", true);
      for (std::size_t k = 0; k < z.coefficients.size(); ++k) {
        check_param(z.coefficients[k], idx(t + "/coefficients", k), false);
      }
      if (!z.normalization_radius.is_bound() && !finite_positive(z.normalization_radius.value)) {
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
      check_param(g->lines_per_mm, loc + "/lines_per_mm", true);
      if (!std::isfinite(g->orientation_deg))
        report("value.not_finite", loc + "/orientation_deg", "orientation must be finite");
      if (!g->lines_per_mm.is_bound() && !finite_positive(g->lines_per_mm.value))
        report("phase.lines_per_mm_invalid", loc + "/lines_per_mm", "must be > 0");
    } else if (const auto* r = std::get_if<RadialPhase>(&p)) {
      check_param(r->normalization_radius, loc + "/normalization_radius", true);
      for (std::size_t k = 0; k < r->coefficients.size(); ++k) {
        check_param(r->coefficients[k], idx(loc + "/coefficients", k), false);
      }
      if (!r->normalization_radius.is_bound() && !finite_positive(r->normalization_radius.value)) {
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

  // ------------------------------------------------ merit function (ADR 0030) -----

  /// The section "optimization" (ADR 0030, points 2-4; #162 part B): names, indices, surfaces on
  /// the path, weights, sampling and finite numbers. Whether an operand can be evaluated on its
  /// path (merit.operand_unsupported) is decided when the optimization starts (#167).
  void check_optimization() {
    const Optimization& o = system_.optimization;
    for (std::size_t i = 0; i < o.operands.size(); ++i) {
      const std::string loc = idx("/optimization/operands", i);
      std::visit([&](const auto& op) { check_operand(op, loc); }, o.operands[i]);
    }
    for (std::size_t i = 0; i < o.generators.size(); ++i) {
      const std::string loc = idx("/optimization/generators", i);
      std::visit([&](const auto& g) { check_generator(g, loc); }, o.generators[i]);
    }
  }

  void check_common(const OperandCommon& c, const std::string& loc) {
    if (!std::isfinite(c.target)) report("value.not_finite", loc + "/target", "must be finite");
    check_weight(c.weight, loc);
    check_configuration(c.configuration, loc);
  }

  void check_weight(double weight, const std::string& loc) {
    if (!std::isfinite(weight) || weight < 0.0) {
      report("merit.weight_invalid", loc + "/weight", "weight must be finite and >= 0");
    }
  }

  void check_configuration(const std::optional<std::string>& name, const std::string& loc) {
    if (!name) return;
    const auto& cs = system_.configurations;
    if (std::none_of(cs.begin(), cs.end(), [&](const Configuration& c) { return c.name == *name; }))
      report("merit.unknown_configuration", loc + "/configuration",
             "unknown configuration '" + *name + "'");
  }

  /// The path of that name, or none (reported as merit.unknown_path).
  const Path* check_path(const std::string& name, const std::string& loc) {
    const auto& ps = system_.paths;
    const auto it =
        std::find_if(ps.begin(), ps.end(), [&](const Path& p) { return p.name == name; });
    if (it != ps.end()) return &*it;
    report("merit.unknown_path", loc + "/path", "unknown path '" + name + "'");
    return nullptr;
  }

  void check_field(std::uint16_t field, const std::string& loc) {
    if (field >= system_.fields.points.size()) {
      report("merit.index_out_of_range", loc,
             "field " + std::to_string(field) + " does not exist (" +
                 std::to_string(system_.fields.points.size()) + " fields)");
    }
  }

  void check_wavelength(std::uint16_t wavelength, const std::string& loc) {
    if (wavelength >= system_.wavelengths.size()) {
      report("merit.index_out_of_range", loc,
             "wavelength " + std::to_string(wavelength) + " does not exist (" +
                 std::to_string(system_.wavelengths.size()) + " wavelengths)");
    }
  }

  void check_sampling(int n, const std::string& loc) {
    if (n < 1) report("merit.sampling_invalid", loc, "must be at least 1");
  }

  void check_operand(const FirstOrderOperand& op, const std::string& loc) {
    check_common(op.common, loc);
    check_path(op.path, loc);
    if (op.wavelength) check_wavelength(*op.wavelength, loc + "/wavelength");
  }

  void check_operand(const RayOperand& op, const std::string& loc) {
    check_common(op.common, loc);
    if (const Path* p = check_path(op.path, loc)) check_ray_surface(*p, op, loc);
    check_field(op.field, loc + "/field");
    if (!std::isfinite(op.px)) report("value.not_finite", loc + "/px", "must be finite");
    if (!std::isfinite(op.py)) report("value.not_finite", loc + "/py", "must be finite");
    if (op.wavelength) check_wavelength(*op.wavelength, loc + "/wavelength");
  }

  /// The surface of a ray operand must be on its path; `occurrence` picks one of several events
  /// at it. An automatic path meets every surface of the system once.
  void check_ray_surface(const Path& p, const RayOperand& op, const std::string& loc) {
    const std::string& id = op.surface.str();
    std::size_t count = 0;
    if (p.automatic) {
      count = surface_ids_.contains(id) ? 1 : 0;
    } else {
      count = static_cast<std::size_t>(
          std::count_if(p.events.begin(), p.events.end(),
                        [&](const Event& e) { return e.surface == op.surface; }));
    }
    if (count == 0) {
      report("merit.surface_not_on_path", loc + "/surface",
             "path '" + p.name + "' does not meet surface '" + id + "'");
    } else if (op.occurrence && *op.occurrence >= count) {
      report("merit.surface_not_on_path", loc + "/occurrence",
             "path '" + p.name + "' meets surface '" + id + "' " + std::to_string(count) +
                 " time(s); occurrence is 0-based");
    } else if (!op.occurrence && count > 1) {
      report("merit.surface_ambiguous", loc + "/surface",
             "path '" + p.name + "' meets surface '" + id + "' " + std::to_string(count) +
                 " times; set occurrence");
    }
  }

  void check_operand(const SpotRmsOperand& op, const std::string& loc) {
    check_common(op.common, loc);
    check_path(op.path, loc);
    check_field(op.field, loc + "/field");
    if (op.wavelength) {
      if (op.polychromatic) {
        report("merit.polychromatic_wavelength", loc + "/wavelength",
               "a polychromatic spot uses all wavelengths; it has no wavelength");
      } else {
        check_wavelength(*op.wavelength, loc + "/wavelength");
      }
    }
    check_sampling(op.rings, loc + "/rings");
  }

  void check_operand(const OpdRmsOperand& op, const std::string& loc) {
    check_common(op.common, loc);
    check_path(op.path, loc);
    check_field(op.field, loc + "/field");
    if (op.wavelength) check_wavelength(*op.wavelength, loc + "/wavelength");
    check_sampling(op.grid, loc + "/grid");
  }

  void check_operand(const ParamValueOperand& op, const std::string& loc) {
    check_common(op.common, loc);
    if (!row_names_.contains(op.parameter)) {
      report("merit.unknown_parameter", loc + "/parameter",
             "unknown parameter row '" + op.parameter + "'");
    }
  }

  template <class G>
  void check_generator(const G& g, const std::string& loc) {
    check_path(g.path, loc);
    check_configuration(g.configuration, loc);
    check_selection(g.fields, loc + "/fields",
                    [&](std::uint16_t k, const std::string& at) { check_field(k, at); });
    check_selection(g.wavelengths, loc + "/wavelengths",
                    [&](std::uint16_t k, const std::string& at) { check_wavelength(k, at); });
    check_sampling(g.rings, loc + "/rings");
    check_sampling(g.arms, loc + "/arms");
    check_weight(g.weight, loc);
  }

  template <class F>
  void check_selection(const std::optional<std::vector<std::uint16_t>>& list,
                       const std::string& loc,
                       F check) {
    if (!list) return;
    if (list->empty()) {
      report("merit.selection_empty", loc, "a selection is not empty; omit it for all");
      return;
    }
    for (std::size_t k = 0; k < list->size(); ++k) check((*list)[k], idx(loc, k));
  }

  const System& system_;
  std::vector<Diagnostic> out_;
  std::unordered_map<std::string, std::string> surface_ids_;
  std::unordered_set<std::string> phase_surfaces_;  // ids of surfaces with a phase layer
  std::unordered_map<std::string, std::string> node_names_;
  std::vector<std::string> stops_;
  std::unordered_set<std::string> row_names_;  // names of the parameter rows (ADR 0029)
  std::size_t surfaces_seen_ = 0;              // surfaces before the node in tree order
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
