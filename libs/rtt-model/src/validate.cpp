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
  void error(std::string location, std::string message) {
    out_.push_back({Severity::Error, std::move(location), std::move(message)});
  }

  void warning(std::string location, std::string message) {
    out_.push_back({Severity::Warning, std::move(location), std::move(message)});
  }

  void check_wavelengths() {
    if (system_.wavelengths.empty()) {
      error("/wavelengths", "at least one wavelength is required");
      return;
    }
    std::size_t references = 0;
    for (std::size_t i = 0; i < system_.wavelengths.size(); ++i) {
      const Wavelength& w = system_.wavelengths[i];
      const std::string loc = idx("/wavelengths", i);
      if (!finite_positive(w.um)) error(loc + "/um", "wavelength must be finite and > 0 um");
      if (!std::isfinite(w.weight) || w.weight < 0.0) error(loc + "/weight", "weight must be >= 0");
      if (w.reference) ++references;
    }
    if (references != 1) {
      error("/wavelengths",
            "exactly one wavelength must be the reference, found " + std::to_string(references));
    }
  }

  void check_aperture() {
    if (system_.aperture.type != SystemApertureType::StopSize &&
        !finite_positive(system_.aperture.value.value)) {
      error("/aperture/value", "aperture value must be finite and > 0");
    }
    if (system_.aperture.type == SystemApertureType::ObjectSpaceNA &&
        system_.aperture.value.value >= 1.0 && system_.environment.medium == "AIR") {
      warning("/aperture/value", "object-space NA >= 1 in air is not physical");
    }
  }

  void check_fields() {
    if (system_.fields.points.empty()) {
      error("/fields/points", "at least one field point is required");
    }
    for (std::size_t i = 0; i < system_.fields.points.size(); ++i) {
      const Field& f = system_.fields.points[i];
      const std::string loc = idx("/fields/points", i);
      if (!std::isfinite(f.x) || !std::isfinite(f.y))
        error(loc, "field coordinates must be finite");
      if (!std::isfinite(f.weight) || f.weight < 0.0) error(loc + "/weight", "weight must be >= 0");
    }
  }

  void check_object_and_environment() {
    if (!system_.object.at_infinity && !finite_positive(system_.object.distance.value)) {
      error("/object/distance", "finite object distance must be > 0 mm");
    }
    const Environment& e = system_.environment;
    if (!std::isfinite(e.temperature_c) || e.temperature_c <= -273.15) {
      error("/environment/temperature_c", "temperature must be above absolute zero");
    }
    if (!std::isfinite(e.pressure_atm) || e.pressure_atm < 0.0) {
      error("/environment/pressure_atm", "pressure must be >= 0 atm");
    }
    if (e.medium.empty()) error("/environment/medium", "medium must not be empty");
  }

  void check_assembly(const Assembly& a, const std::string& loc) {
    if (a.name.empty()) error(loc + "/name", "assembly name must not be empty");
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
      error(loc + "/name", "duplicate node name '" + name + "' (first at " + it->second + ")");
  }

  void check_element(const Element& e, const std::string& loc) {
    if (e.name.empty()) error(loc + "/name", "element name must not be empty");
    check_name(e.name, loc);
    const std::size_t n = e.surfaces.size();
    switch (e.kind) {
      case ElementKind::Lens:
        if (n < 2) error(loc + "/surfaces", "a lens needs at least 2 surfaces");
        if (!e.material) error(loc + "/material", "a lens needs a material");
        break;
      case ElementKind::Plate:
        if (n < 2) error(loc + "/surfaces", "a plate needs at least 2 surfaces");
        if (!e.material) error(loc + "/material", "a plate needs a material");
        for (std::size_t i = 0; i < n; ++i) {
          if (!std::holds_alternative<Plane>(e.surfaces[i].shape.base)) {
            error(idx(loc + "/surfaces", i) + "/shape", "plate surfaces must be planes");
          }
        }
        break;
      case ElementKind::Mirror:
        if (n < 1) error(loc + "/surfaces", "a mirror needs at least 1 surface");
        break;
      case ElementKind::ThinElement:
      case ElementKind::Stop:
      case ElementKind::Detector:
        if (n != 1) error(loc + "/surfaces", "this element kind needs exactly 1 surface");
        if (e.material) error(loc + "/material", "this element kind has no material");
        break;
    }
    if (e.kind == ElementKind::Stop) {
      stops_.push_back(loc);
      if (n == 1 && !e.surfaces[0].aperture) {
        error(loc + "/surfaces/0/aperture", "a stop needs an aperture");
      }
    }
    if (e.material && e.material->empty()) error(loc + "/material", "material must not be empty");
    for (std::size_t i = 0; i < n; ++i) check_surface(e.surfaces[i], idx(loc + "/surfaces", i));
  }

  void check_surface(const Surface& s, const std::string& loc) {
    if (s.id.empty()) {
      error(loc + "/id", "surface id must not be empty");
    } else {
      auto [it, inserted] = surface_ids_.emplace(s.id.str(), loc);
      if (!inserted) {
        error(loc + "/id",
              "duplicate surface id '" + s.id.str() + "' (first at " + it->second + ")");
      }
    }
    check_shape(s.shape, loc + "/shape");
    if (s.aperture) check_aperture(*s.aperture, loc + "/aperture");
    for (std::size_t i = 0; i < s.phases.size(); ++i) {
      check_phase(s.phases[i], idx(loc + "/phases", i));
    }
    check_interaction(s.interaction, loc + "/interaction");
  }

  void check_radius(const Param& r, const std::string& loc) {
    if (!std::isfinite(r.value) || r.value == 0.0) {
      error(loc, "radius must be finite and non-zero (use a plane for flat surfaces)");
    }
  }

  void check_shape(const ShapeStack& shape, const std::string& loc) {
    const std::string base = loc + "/base";
    if (const auto* c = std::get_if<Conic>(&shape.base)) {
      check_radius(c->radius, base + "/radius");
    } else if (const auto* a = std::get_if<EvenAsphere>(&shape.base)) {
      check_radius(a->radius, base + "/radius");
      if (a->coefficients.empty()) warning(base + "/coefficients", "asphere without coefficients");
    }
    for (std::size_t i = 0; i < shape.terms.size(); ++i) {
      const auto& z = std::get<ZernikeSag>(shape.terms[i]);
      const std::string t = idx(loc + "/terms", i);
      if (!finite_positive(z.normalization_radius.value)) {
        error(t + "/normalization_radius", "normalization radius must be > 0 mm");
      }
      if (z.coefficients.empty()) warning(t + "/coefficients", "Zernike term without coefficients");
    }
  }

  void check_aperture(const Aperture& a, const std::string& loc) {
    if (const auto* c = std::get_if<CircularAperture>(&a)) {
      if (!finite_positive(c->radius)) error(loc + "/radius", "radius must be > 0 mm");
      if (!std::isfinite(c->inner_radius) || c->inner_radius < 0.0 ||
          c->inner_radius >= c->radius) {
        error(loc + "/inner_radius", "inner radius must be >= 0 and < radius");
      }
    } else if (const auto* r = std::get_if<RectangularAperture>(&a)) {
      if (!finite_positive(r->half_width_x) || !finite_positive(r->half_width_y)) {
        error(loc, "half widths must be > 0 mm");
      }
    } else if (const auto* e = std::get_if<EllipticalAperture>(&a)) {
      if (!finite_positive(e->semi_axis_x) || !finite_positive(e->semi_axis_y)) {
        error(loc, "semi axes must be > 0 mm");
      }
    }
  }

  void check_phase(const PhaseLayer& p, const std::string& loc) {
    if (const auto* g = std::get_if<LinearGrating>(&p)) {
      if (!finite_positive(g->lines_per_mm.value)) error(loc + "/lines_per_mm", "must be > 0");
    } else if (const auto* r = std::get_if<RadialPhase>(&p)) {
      if (!finite_positive(r->normalization_radius.value)) {
        error(loc + "/normalization_radius", "normalization radius must be > 0 mm");
      }
    }
  }

  void check_interaction(const Interaction& i, const std::string& loc) {
    if (const auto* bs = std::get_if<IdealBeamSplitter>(&i)) {
      for (const double r : {bs->reflectance_s, bs->reflectance_p}) {
        if (!std::isfinite(r) || r < 0.0 || r > 1.0) error(loc, "reflectance must be in [0, 1]");
      }
    } else if (const auto* c = std::get_if<CoatingRef>(&i)) {
      if (c->name.empty()) error(loc + "/name", "coating name must not be empty");
    } else if (const auto* pol = std::get_if<IdealPolarizer>(&i)) {
      if (!nonzero_axis(pol->transmission_axis))
        error(loc + "/transmission_axis", "axis must be non-zero");
      if (!std::isfinite(pol->extinction_ratio) || pol->extinction_ratio < 0.0 ||
          pol->extinction_ratio > 1.0) {
        error(loc + "/extinction_ratio", "extinction ratio must be in [0, 1]");
      }
    } else if (const auto* ret = std::get_if<IdealRetarder>(&i)) {
      if (!nonzero_axis(ret->fast_axis)) error(loc + "/fast_axis", "axis must be non-zero");
      if (!std::isfinite(ret->retardance_waves)) error(loc + "/retardance_waves", "must be finite");
    }
  }

  void check_stops() {
    if (stops_.size() > 1) error(stops_[1], "only one stop element is allowed");
    if (stops_.empty() && system_.aperture.type == SystemApertureType::StopSize) {
      error("/aperture/type", "aperture type 'stop_size' requires a stop element");
    }
  }

  void check_paths() {
    if (system_.paths.empty()) {
      error("/paths", "at least one path is required");
      return;
    }
    std::unordered_set<std::string> names;
    for (std::size_t i = 0; i < system_.paths.size(); ++i) {
      const Path& p = system_.paths[i];
      const std::string loc = idx("/paths", i);
      if (p.name.empty()) error(loc + "/name", "path name must not be empty");
      if (!names.insert(p.name).second)
        error(loc + "/name", "duplicate path name '" + p.name + "'");
      if (!p.automatic && p.events.empty()) error(loc + "/events", "explicit path without events");
      for (std::size_t k = 0; k < p.events.size(); ++k) {
        const Event& e = p.events[k];
        const std::string eloc = idx(loc + "/events", k);
        if (!surface_ids_.contains(e.surface.str())) {
          error(eloc + "/surface", "unknown surface id '" + e.surface.str() + "'");
        }
        if (e.kind != EventKind::Diffract && e.order != 0) {
          error(eloc + "/order", "order is only allowed for diffract events");
        }
      }
    }
  }

  const System& system_;
  std::vector<Diagnostic> out_;
  std::unordered_map<std::string, std::string> surface_ids_;
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
  return std::string(diagnostic.severity == Severity::Error ? "error " : "warning ") +
         diagnostic.location + ": " + diagnostic.message;
}

}  // namespace rtt::model
