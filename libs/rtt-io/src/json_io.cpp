#include "rtt/io/json_io.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <initializer_list>
#include <limits>
#include <nlohmann/json.hpp>
#include <sstream>
#include <stdexcept>
#include <utility>

#include "json_format.hpp"
#include "json_tree.hpp"
#include "rtt/json/strict.hpp"

namespace rtt::io {

namespace {

using Json = nlohmann::json;
using OJson = nlohmann::ordered_json;
using namespace rtt::model;

// ===================================================================== tables =====

template <class E>
using EnumTable = std::initializer_list<std::pair<std::string_view, E>>;

const EnumTable<ElementKind> kElementKinds = {
    {"lens", ElementKind::Lens},   {"mirror", ElementKind::Mirror},
    {"plate", ElementKind::Plate}, {"thin_element", ElementKind::ThinElement},
    {"stop", ElementKind::Stop},   {"detector", ElementKind::Detector},
};

const EnumTable<EventKind> kEventKinds = {
    {"refract", EventKind::Refract},
    {"reflect", EventKind::Reflect},
    {"transmit", EventKind::Transmit},
    {"ordinary", EventKind::Ordinary},
    {"extraordinary", EventKind::Extraordinary},
    {"diffract", EventKind::Diffract},
};

const EnumTable<FieldType> kFieldTypes = {
    {"angle_deg", FieldType::AngleDeg},
    {"object_height", FieldType::ObjectHeight},
    {"paraxial_image_height", FieldType::ParaxialImageHeight},
};

const EnumTable<SystemApertureType> kApertureTypes = {
    {"epd", SystemApertureType::EntrancePupilDiameter},
    {"image_fnumber", SystemApertureType::ImageSpaceFNumber},
    {"object_na", SystemApertureType::ObjectSpaceNA},
    {"stop_size", SystemApertureType::StopSize},
};

template <class E>
std::string_view enum_name(E value, const EnumTable<E>& table) {
  for (const auto& [name, v] : table) {
    if (v == value) return name;
  }
  return "?";
}

// ===================================================================== reader =====

/// Location of the value being read, as a JSON pointer.
class Ctx {
 public:
  Ctx() = default;
  explicit Ctx(std::string pointer) : pointer_(std::move(pointer)) {}
  /// The key as a reference token: "~" and "/" are masked (RFC 6901, Sec. 3).
  [[nodiscard]] Ctx at(std::string_view key) const {
    return Ctx(pointer_ + "/" + rtt::json::pointer_token(key));
  }
  [[nodiscard]] Ctx at(std::size_t i) const { return Ctx(pointer_ + "/" + std::to_string(i)); }
  [[noreturn]] void fail(const std::string& message) const { throw ParseError(pointer_, message); }

 private:
  std::string pointer_;
};

std::string_view type_name(const Json& j) {
  return j.type_name();
}

/// Checks that `j` is an object containing only allowed keys.
void expect_object(const Json& j, const Ctx& c, std::initializer_list<std::string_view> allowed) {
  if (!j.is_object()) c.fail("expected an object, got " + std::string(type_name(j)));
  for (const auto& item : j.items()) {
    bool ok = false;
    for (const std::string_view key : allowed) ok = ok || item.key() == key;
    if (!ok) c.at(item.key()).fail("unknown key '" + item.key() + "'");
  }
}

const Json* find(const Json& j, std::string_view key) {
  const auto it = j.find(std::string(key));
  return it == j.end() ? nullptr : &*it;
}

const Json& require(const Json& j, std::string_view key, const Ctx& c) {
  const Json* v = find(j, key);
  if (v == nullptr) c.fail("missing required key '" + std::string(key) + "'");
  return *v;
}

double read_number(const Json& j, const Ctx& c) {
  if (!j.is_number()) c.fail("expected a number, got " + std::string(type_name(j)));
  return j.get<double>();
}

int read_int(const Json& j, const Ctx& c) {
  if (!j.is_number_integer()) c.fail("expected an integer, got " + std::string(type_name(j)));
  // Strict (ADR 0008): a value outside int is an error, not cut off by get<int>() (#35).
  constexpr auto kMin = std::numeric_limits<int>::min();
  constexpr auto kMax = std::numeric_limits<int>::max();
  const bool in_range = j.is_number_unsigned()
                            ? j.get<std::uint64_t>() <= static_cast<std::uint64_t>(kMax)
                            : j.get<std::int64_t>() >= kMin && j.get<std::int64_t>() <= kMax;
  if (!in_range) {
    c.fail("integer " + j.dump() + " out of range (" + std::to_string(kMin) + " ... " +
           std::to_string(kMax) + ")");
  }
  return j.get<int>();
}

bool read_bool(const Json& j, const Ctx& c) {
  if (!j.is_boolean()) c.fail("expected true or false, got " + std::string(type_name(j)));
  return j.get<bool>();
}

std::string read_string(const Json& j, const Ctx& c) {
  if (!j.is_string()) c.fail("expected a string, got " + std::string(type_name(j)));
  return j.get<std::string>();
}

const Json& read_array(const Json& j, const Ctx& c) {
  if (!j.is_array()) c.fail("expected an array, got " + std::string(type_name(j)));
  return j;
}

/// Required key whose value must be an array.
const Json& require_array(const Json& j, std::string_view key, const Ctx& c) {
  const Json& v = require(j, key, c);
  if (!v.is_array()) c.at(key).fail("expected an array, got " + std::string(type_name(v)));
  return v;
}

template <class E>
E read_enum(const Json& j, const Ctx& c, const EnumTable<E>& table) {
  const std::string s = read_string(j, c);
  for (const auto& [name, v] : table) {
    if (name == s) return v;
  }
  std::string allowed;
  for (const auto& entry : table) {
    allowed += (allowed.empty() ? "" : ", ") + std::string(entry.first);
  }
  c.fail("unknown value '" + s + "', expected one of: " + allowed);
}

/// Param: either a plain number or {"value": n, "variable": b, "pickup": "expr"}.
Param read_param(const Json& j, const Ctx& c) {
  if (j.is_number()) return Param(j.get<double>());
  if (!j.is_object()) c.fail("expected a number or a parameter object");
  expect_object(j, c, {"value", "variable", "pickup"});
  Param p(read_number(require(j, "value", c), c.at("value")));
  if (const Json* v = find(j, "variable")) p.variable = read_bool(*v, c.at("variable"));
  if (const Json* v = find(j, "pickup")) p.pickup = read_string(*v, c.at("pickup"));
  return p;
}

std::vector<Param> read_param_list(const Json& j, const Ctx& c) {
  std::vector<Param> out;
  const Json& a = read_array(j, c);
  out.reserve(a.size());
  for (std::size_t i = 0; i < a.size(); ++i) out.push_back(read_param(a[i], c.at(i)));
  return out;
}

std::array<Param, 3> read_param3(const Json& j, const Ctx& c) {
  const Json& a = read_array(j, c);
  if (a.size() != 3) c.fail("expected exactly 3 values");
  return {read_param(a[0], c.at(0)), read_param(a[1], c.at(1)), read_param(a[2], c.at(2))};
}

std::array<double, 3> read_vec3(const Json& j, const Ctx& c) {
  const Json& a = read_array(j, c);
  if (a.size() != 3) c.fail("expected exactly 3 numbers");
  return {read_number(a[0], c.at(0)), read_number(a[1], c.at(1)), read_number(a[2], c.at(2))};
}

template <class T, class F>
void read_opt(const Json& j, std::string_view key, const Ctx& c, T& target, F reader) {
  if (const Json* v = find(j, key)) target = reader(*v, c.at(key));
}

Pose read_pose(const Json& j, const Ctx& c) {
  expect_object(j, c, {"position", "rotation_deg", "pivot"});
  Pose p;
  read_opt(j, "position", c, p.position, read_param3);
  read_opt(j, "rotation_deg", c, p.rotation_deg, read_param3);
  read_opt(j, "pivot", c, p.pivot, read_vec3);
  return p;
}

std::string read_type(const Json& j, const Ctx& c) {
  if (!j.is_object()) c.fail("expected an object, got " + std::string(type_name(j)));
  return read_string(require(j, "type", c), c.at("type"));
}

BaseShape read_base_shape(const Json& j, const Ctx& c) {
  const std::string type = read_type(j, c);
  if (type == "plane") {
    expect_object(j, c, {"type"});
    return Plane{};
  }
  if (type == "conic") {
    expect_object(j, c, {"type", "radius", "conic"});
    Conic s;
    s.radius = read_param(require(j, "radius", c), c.at("radius"));
    read_opt(j, "conic", c, s.conic, read_param);
    return s;
  }
  if (type == "even_asphere") {
    expect_object(j, c, {"type", "radius", "conic", "coefficients"});
    EvenAsphere s;
    s.radius = read_param(require(j, "radius", c), c.at("radius"));
    read_opt(j, "conic", c, s.conic, read_param);
    read_opt(j, "coefficients", c, s.coefficients, read_param_list);
    return s;
  }
  c.at("type").fail("unknown shape type '" + type + "', expected plane, conic or even_asphere");
}

ShapeTerm read_shape_term(const Json& j, const Ctx& c) {
  const std::string type = read_type(j, c);
  if (type == "zernike_sag") {
    expect_object(j, c, {"type", "normalization_radius", "coefficients"});
    ZernikeSag z;
    z.normalization_radius =
        read_param(require(j, "normalization_radius", c), c.at("normalization_radius"));
    read_opt(j, "coefficients", c, z.coefficients, read_param_list);
    return z;
  }
  c.at("type").fail("unknown shape term type '" + type + "', expected zernike_sag");
}

ShapeStack read_shape(const Json& j, const Ctx& c) {
  expect_object(j, c, {"base", "terms"});
  ShapeStack s;
  read_opt(j, "base", c, s.base, read_base_shape);
  if (const Json* t = find(j, "terms")) {
    const Ctx tc = c.at("terms");
    const Json& a = read_array(*t, tc);
    for (std::size_t i = 0; i < a.size(); ++i) {
      s.terms.push_back(read_shape_term(a[i], c.at("terms").at(i)));
    }
  }
  return s;
}

Aperture read_aperture(const Json& j, const Ctx& c) {
  const std::string type = read_type(j, c);
  if (type == "circular") {
    expect_object(j, c, {"type", "radius", "inner_radius"});
    CircularAperture a;
    a.radius = read_number(require(j, "radius", c), c.at("radius"));
    read_opt(j, "inner_radius", c, a.inner_radius, read_number);
    return a;
  }
  if (type == "rectangular") {
    expect_object(j, c, {"type", "half_width_x", "half_width_y"});
    return RectangularAperture{read_number(require(j, "half_width_x", c), c.at("half_width_x")),
                               read_number(require(j, "half_width_y", c), c.at("half_width_y"))};
  }
  if (type == "elliptical") {
    expect_object(j, c, {"type", "semi_axis_x", "semi_axis_y"});
    return EllipticalAperture{read_number(require(j, "semi_axis_x", c), c.at("semi_axis_x")),
                              read_number(require(j, "semi_axis_y", c), c.at("semi_axis_y"))};
  }
  c.at("type").fail("unknown aperture type '" + type +
                    "', expected circular, rectangular or elliptical");
}

PhaseLayer read_phase(const Json& j, const Ctx& c) {
  const std::string type = read_type(j, c);
  if (type == "linear_grating") {
    expect_object(j, c, {"type", "lines_per_mm", "orientation_deg"});
    LinearGrating g;
    g.lines_per_mm = read_param(require(j, "lines_per_mm", c), c.at("lines_per_mm"));
    read_opt(j, "orientation_deg", c, g.orientation_deg, read_number);
    return g;
  }
  if (type == "radial_phase") {
    expect_object(j, c, {"type", "normalization_radius", "coefficients"});
    RadialPhase r;
    r.normalization_radius =
        read_param(require(j, "normalization_radius", c), c.at("normalization_radius"));
    read_opt(j, "coefficients", c, r.coefficients, read_param_list);
    return r;
  }
  c.at("type").fail("unknown phase type '" + type + "', expected linear_grating or radial_phase");
}

Interaction read_interaction(const Json& j, const Ctx& c) {
  const std::string type = read_type(j, c);
  if (type == "fresnel") {
    expect_object(j, c, {"type"});
    return Fresnel{};
  }
  if (type == "ideal_mirror") {
    expect_object(j, c, {"type"});
    return IdealMirror{};
  }
  if (type == "ideal_anti_reflection") {
    expect_object(j, c, {"type"});
    return IdealAntiReflection{};
  }
  if (type == "absorber") {
    expect_object(j, c, {"type"});
    return Absorber{};
  }
  if (type == "ideal_beam_splitter") {
    expect_object(j, c, {"type", "reflectance_s", "reflectance_p"});
    IdealBeamSplitter b;
    read_opt(j, "reflectance_s", c, b.reflectance_s, read_number);
    read_opt(j, "reflectance_p", c, b.reflectance_p, read_number);
    return b;
  }
  if (type == "coating") {
    expect_object(j, c, {"type", "name"});
    return CoatingRef{read_string(require(j, "name", c), c.at("name"))};
  }
  if (type == "ideal_polarizer") {
    expect_object(j, c, {"type", "transmission_axis", "extinction_ratio"});
    IdealPolarizer p;
    read_opt(j, "transmission_axis", c, p.transmission_axis, read_vec3);
    read_opt(j, "extinction_ratio", c, p.extinction_ratio, read_number);
    return p;
  }
  if (type == "ideal_retarder") {
    expect_object(j, c, {"type", "fast_axis", "retardance_waves"});
    IdealRetarder r;
    read_opt(j, "fast_axis", c, r.fast_axis, read_vec3);
    read_opt(j, "retardance_waves", c, r.retardance_waves, read_number);
    return r;
  }
  c.at("type").fail("unknown interaction type '" + type + "'");
}

Surface read_surface(const Json& j, const Ctx& c) {
  expect_object(j, c, {"id", "pose", "shape", "aperture", "phases", "interaction"});
  Surface s;
  s.id = SurfaceId(read_string(require(j, "id", c), c.at("id")));
  read_opt(j, "pose", c, s.pose, read_pose);
  read_opt(j, "shape", c, s.shape, read_shape);
  if (const Json* a = find(j, "aperture")) s.aperture = read_aperture(*a, c.at("aperture"));
  if (const Json* p = find(j, "phases")) {
    const Ctx pc = c.at("phases");
    const Json& a = read_array(*p, pc);
    for (std::size_t i = 0; i < a.size(); ++i)
      s.phases.push_back(read_phase(a[i], c.at("phases").at(i)));
  }
  read_opt(j, "interaction", c, s.interaction, read_interaction);
  return s;
}

/// Element material: a string (shorthand, one material for all segments) or a non-empty
/// array of strings (one material per segment, ADR 0017).
void read_material(const Json& j, const Ctx& c, Element& e) {
  if (j.is_string()) {
    e.material = j.get<std::string>();
    return;
  }
  if (!j.is_array()) {
    c.fail("expected a string or an array of strings, got " + std::string(type_name(j)));
  }
  if (j.empty()) c.fail("expected at least one material");
  e.segment_materials.reserve(j.size());
  for (std::size_t i = 0; i < j.size(); ++i) {
    e.segment_materials.push_back(read_string(j[i], c.at(i)));
  }
}

Node read_node(const Json& j, const Ctx& c);

Assembly read_assembly_body(const Json& j, const Ctx& c) {
  expect_object(j, c, {"type", "name", "pose", "children"});
  Assembly a;
  a.name = read_string(require(j, "name", c), c.at("name"));
  read_opt(j, "pose", c, a.pose, read_pose);
  const Json& children = require_array(j, "children", c);
  a.children.reserve(children.size());
  for (std::size_t i = 0; i < children.size(); ++i) {
    a.children.push_back(read_node(children[i], c.at("children").at(i)));
  }
  return a;
}

Node read_node(const Json& j, const Ctx& c) {
  const std::string type = read_type(j, c);
  if (type == "assembly") return Node{read_assembly_body(j, c)};
  expect_object(j, c, {"type", "name", "pose", "material", "surfaces"});
  Element e;
  e.kind = read_enum(require(j, "type", c), c.at("type"), kElementKinds);
  e.name = read_string(require(j, "name", c), c.at("name"));
  read_opt(j, "pose", c, e.pose, read_pose);
  if (const Json* m = find(j, "material")) read_material(*m, c.at("material"), e);
  const Json& surfaces = require_array(j, "surfaces", c);
  e.surfaces.reserve(surfaces.size());
  for (std::size_t i = 0; i < surfaces.size(); ++i) {
    e.surfaces.push_back(read_surface(surfaces[i], c.at("surfaces").at(i)));
  }
  return Node{std::move(e)};
}

Path read_path(const Json& j, const Ctx& c) {
  expect_object(j, c, {"name", "events"});
  Path p;
  p.name = read_string(require(j, "name", c), c.at("name"));
  const Json& events = require(j, "events", c);
  const Ctx ec = c.at("events");
  if (events.is_string()) {
    if (events.get<std::string>() != "auto") ec.fail("expected \"auto\" or a list of events");
    p.automatic = true;
    return p;
  }
  const Json& a = read_array(events, ec);
  for (std::size_t i = 0; i < a.size(); ++i) {
    const Ctx c2 = ec.at(i);
    expect_object(a[i], c2, {"surface", "kind", "order"});
    Event e;
    e.surface = SurfaceId(read_string(require(a[i], "surface", c2), c2.at("surface")));
    if (const Json* k = find(a[i], "kind")) e.kind = read_enum(*k, c2.at("kind"), kEventKinds);
    read_opt(a[i], "order", c2, e.order, read_int);
    p.events.push_back(std::move(e));
  }
  return p;
}

// ================================================================== migration =====

/// "major.minor" of a version string with exactly two dots, otherwise empty.
std::string major_minor(const std::string& version) {
  if (std::count(version.begin(), version.end(), '.') != 2) return {};
  return version.substr(0, version.rfind('.'));
}

/// Migration 0.1 -> 0.2 (ADR 0017). 0.2 only adds the array form of "material"; every 0.1 file
/// is a 0.2 file with the same meaning. An array material in a file that claims 0.1 is an error.
void check_no_material_lists(const Json& node, const Ctx& c) {
  if (!node.is_object()) return;  // structural errors are reported by the reader
  if (const Json* m = find(node, "material"); m != nullptr && m->is_array()) {
    c.at("material").fail("material lists need schema_version 0.2 or later");
  }
  if (const Json* children = find(node, "children"); children != nullptr && children->is_array()) {
    for (std::size_t i = 0; i < children->size(); ++i) {
      check_no_material_lists((*children)[i], c.at("children").at(i));
    }
  }
}

/// Checks the schema version of `j` and migrates supported older versions to the current one.
/// Pre-1.0 rule: major and minor must match a supported version, patch may differ. The reader
/// then reads `j` as a current file and the writer writes kSchemaVersion.
void migrate(const Json& j, const std::string& version, const Ctx& c) {
  const std::string ours(kSchemaVersion);
  const std::string mm = major_minor(version);
  if (!mm.empty() && mm == major_minor(ours)) return;
  if (mm == "0.1") {
    // 0.1 -> 0.2: content unchanged, only lists are new.
    if (const Json* root = find(j, "root")) check_no_material_lists(*root, Ctx().at("root"));
    return;
  }
  c.fail("incompatible schema_version '" + version + "', this build reads " + ours +
         " and migrates 0.1");
}

System read_system_tree(const Json& j) {
  const Ctx c;
  expect_object(j, c,
                {"schema_version", "name", "units", "environment", "object", "wavelengths",
                 "aperture", "fields", "root", "paths"});
  System s;
  s.schema_version = read_string(require(j, "schema_version", c), c.at("schema_version"));
  migrate(j, s.schema_version, c.at("schema_version"));
  s.schema_version = std::string(kSchemaVersion);
  read_opt(j, "name", c, s.name, read_string);

  {
    const Ctx u = c.at("units");
    const Json& units = require(j, "units", c);
    expect_object(units, u, {"length", "wavelength"});
    if (read_string(require(units, "length", u), u.at("length")) != "mm") {
      u.at("length").fail("length unit must be \"mm\"");
    }
    if (read_string(require(units, "wavelength", u), u.at("wavelength")) != "um") {
      u.at("wavelength").fail("wavelength unit must be \"um\"");
    }
  }

  if (const Json* e = find(j, "environment")) {
    const Ctx ec = c.at("environment");
    expect_object(*e, ec, {"temperature_c", "pressure_atm", "medium"});
    read_opt(*e, "temperature_c", ec, s.environment.temperature_c, read_number);
    read_opt(*e, "pressure_atm", ec, s.environment.pressure_atm, read_number);
    read_opt(*e, "medium", ec, s.environment.medium, read_string);
  }

  if (const Json* o = find(j, "object")) {
    const Ctx oc = c.at("object");
    expect_object(*o, oc, {"at_infinity", "distance"});
    s.object.at_infinity = read_bool(require(*o, "at_infinity", oc), oc.at("at_infinity"));
    read_opt(*o, "distance", oc, s.object.distance, read_param);
  }

  {
    const Ctx wc = c.at("wavelengths");
    const Json& a = require_array(j, "wavelengths", c);
    for (std::size_t i = 0; i < a.size(); ++i) {
      const Ctx c2 = wc.at(i);
      expect_object(a[i], c2, {"um", "weight", "reference"});
      Wavelength w;
      w.um = read_number(require(a[i], "um", c2), c2.at("um"));
      read_opt(a[i], "weight", c2, w.weight, read_number);
      read_opt(a[i], "reference", c2, w.reference, read_bool);
      s.wavelengths.push_back(w);
    }
  }

  {
    const Ctx ac = c.at("aperture");
    const Json& a = require(j, "aperture", c);
    expect_object(a, ac, {"type", "value"});
    s.aperture.type = read_enum(require(a, "type", ac), ac.at("type"), kApertureTypes);
    if (s.aperture.type == SystemApertureType::StopSize) {
      read_opt(a, "value", ac, s.aperture.value, read_param);
    } else {
      s.aperture.value = read_param(require(a, "value", ac), ac.at("value"));
    }
  }

  {
    const Ctx fc = c.at("fields");
    const Json& f = require(j, "fields", c);
    expect_object(f, fc, {"type", "points"});
    if (const Json* t = find(f, "type")) s.fields.type = read_enum(*t, fc.at("type"), kFieldTypes);
    const Json& pts = require_array(f, "points", fc);
    for (std::size_t i = 0; i < pts.size(); ++i) {
      const Ctx pc = fc.at("points").at(i);
      expect_object(pts[i], pc, {"x", "y", "weight"});
      Field p;
      read_opt(pts[i], "x", pc, p.x, read_number);
      read_opt(pts[i], "y", pc, p.y, read_number);
      read_opt(pts[i], "weight", pc, p.weight, read_number);
      s.fields.points.push_back(p);
    }
  }

  {
    const Ctx rc = c.at("root");
    const Json& r = require(j, "root", c);
    if (read_type(r, rc) != "assembly") rc.at("type").fail("root must be an assembly");
    s.root = read_assembly_body(r, rc);
  }

  {
    const Ctx pc = c.at("paths");
    const Json& a = require_array(j, "paths", c);
    for (std::size_t i = 0; i < a.size(); ++i) s.paths.push_back(read_path(a[i], pc.at(i)));
  }
  return s;
}

// ===================================================================== writer =====
//
// Two forms (json_tree.hpp): the canonical form omits values equal to their defaults and writes
// a plain Param as a number; the edit form (ADR 0024) writes every value and every Param as an
// object. Optional members without a value (surface aperture, element material, pickup) are
// missing in both.

OJson num(double v) {
  if (!std::isfinite(v)) throw std::invalid_argument("cannot write non-finite number to JSON");
  return OJson(v);
}

class Writer {
 public:
  explicit Writer(detail::Form form) : edit_(form == detail::Form::Edit) {}

  OJson system(const System& s) const {
    OJson o;
    o["schema_version"] = std::string(kSchemaVersion);
    if (edit_ || !s.name.empty()) o["name"] = s.name;
    o["units"] = OJson{{"length", "mm"}, {"wavelength", "um"}};

    if (edit_ || s.environment != Environment{}) {
      const Environment d;
      OJson e;
      if (edit_ || s.environment.temperature_c != d.temperature_c)
        e["temperature_c"] = num(s.environment.temperature_c);
      if (edit_ || s.environment.pressure_atm != d.pressure_atm)
        e["pressure_atm"] = num(s.environment.pressure_atm);
      if (edit_ || s.environment.medium != d.medium) e["medium"] = s.environment.medium;
      o["environment"] = std::move(e);
    }

    if (edit_ || s.object != ObjectSpace{}) {
      OJson ob;
      ob["at_infinity"] = s.object.at_infinity;
      if (edit_ || s.object.distance != Param{}) ob["distance"] = param(s.object.distance);
      o["object"] = std::move(ob);
    }

    OJson wl = OJson::array();
    for (const Wavelength& w : s.wavelengths) {
      OJson wo;
      wo["um"] = num(w.um);
      if (edit_ || w.weight != 1.0) wo["weight"] = num(w.weight);
      if (edit_ || w.reference) wo["reference"] = w.reference;
      wl.push_back(std::move(wo));
    }
    o["wavelengths"] = std::move(wl);

    OJson ap;
    ap["type"] = std::string(enum_name(s.aperture.type, kApertureTypes));
    if (edit_ || s.aperture.type != SystemApertureType::StopSize || s.aperture.value != Param{}) {
      ap["value"] = param(s.aperture.value);
    }
    o["aperture"] = std::move(ap);

    OJson fields;
    if (edit_ || s.fields.type != FieldType::AngleDeg)
      fields["type"] = std::string(enum_name(s.fields.type, kFieldTypes));
    OJson pts = OJson::array();
    for (const Field& f : s.fields.points) {
      OJson fo = OJson::object();
      if (edit_ || f.x != 0.0) fo["x"] = num(f.x);
      if (edit_ || f.y != 0.0) fo["y"] = num(f.y);
      if (edit_ || f.weight != 1.0) fo["weight"] = num(f.weight);
      pts.push_back(std::move(fo));
    }
    fields["points"] = std::move(pts);
    o["fields"] = std::move(fields);

    o["root"] = assembly(s.root);

    OJson paths = OJson::array();
    for (const Path& p : s.paths) paths.push_back(path(p));
    o["paths"] = std::move(paths);
    return o;
  }

 private:
  OJson param(const Param& p) const {
    if (!edit_ && p.is_plain()) return num(p.value);
    OJson o;
    o["value"] = num(p.value);
    if (edit_ || p.variable) o["variable"] = p.variable;
    if (p.pickup) o["pickup"] = *p.pickup;
    return o;
  }

  OJson param_list(const std::vector<Param>& v) const {
    OJson a = OJson::array();
    for (const Param& p : v) a.push_back(param(p));
    return a;
  }

  OJson param3(const std::array<Param, 3>& v) const {
    return OJson::array({param(v[0]), param(v[1]), param(v[2])});
  }

  static OJson vec3(const std::array<double, 3>& v) {
    return OJson::array({num(v[0]), num(v[1]), num(v[2])});
  }

  void put_pose(OJson& o, const Pose& p) const {
    if (!edit_ && p.is_identity()) return;
    const Pose d;
    OJson j;
    if (edit_ || p.position != d.position) j["position"] = param3(p.position);
    if (edit_ || p.rotation_deg != d.rotation_deg) j["rotation_deg"] = param3(p.rotation_deg);
    if (edit_ || p.pivot != d.pivot) j["pivot"] = vec3(p.pivot);
    o["pose"] = std::move(j);
  }

  OJson base_shape(const BaseShape& b) const {
    OJson o;
    if (std::holds_alternative<Plane>(b)) {
      o["type"] = "plane";
    } else if (const auto* c = std::get_if<Conic>(&b)) {
      o["type"] = "conic";
      o["radius"] = param(c->radius);
      if (edit_ || c->conic != Param{}) o["conic"] = param(c->conic);
    } else if (const auto* a = std::get_if<EvenAsphere>(&b)) {
      o["type"] = "even_asphere";
      o["radius"] = param(a->radius);
      if (edit_ || a->conic != Param{}) o["conic"] = param(a->conic);
      if (edit_ || !a->coefficients.empty()) o["coefficients"] = param_list(a->coefficients);
    }
    return o;
  }

  OJson shape(const ShapeStack& s) const {
    OJson o;
    if (edit_ || !std::holds_alternative<Plane>(s.base)) o["base"] = base_shape(s.base);
    if (edit_ || !s.terms.empty()) {
      OJson a = OJson::array();
      for (const ShapeTerm& t : s.terms) {
        const auto& z = std::get<ZernikeSag>(t);
        OJson zo;
        zo["type"] = "zernike_sag";
        zo["normalization_radius"] = param(z.normalization_radius);
        if (edit_ || !z.coefficients.empty()) zo["coefficients"] = param_list(z.coefficients);
        a.push_back(std::move(zo));
      }
      o["terms"] = std::move(a);
    }
    return o;
  }

  OJson aperture(const Aperture& a) const {
    OJson o;
    if (const auto* c = std::get_if<CircularAperture>(&a)) {
      o["type"] = "circular";
      o["radius"] = num(c->radius);
      if (edit_ || c->inner_radius != 0.0) o["inner_radius"] = num(c->inner_radius);
    } else if (const auto* r = std::get_if<RectangularAperture>(&a)) {
      o["type"] = "rectangular";
      o["half_width_x"] = num(r->half_width_x);
      o["half_width_y"] = num(r->half_width_y);
    } else if (const auto* e = std::get_if<EllipticalAperture>(&a)) {
      o["type"] = "elliptical";
      o["semi_axis_x"] = num(e->semi_axis_x);
      o["semi_axis_y"] = num(e->semi_axis_y);
    }
    return o;
  }

  OJson phase(const PhaseLayer& p) const {
    OJson o;
    if (const auto* g = std::get_if<LinearGrating>(&p)) {
      o["type"] = "linear_grating";
      o["lines_per_mm"] = param(g->lines_per_mm);
      if (edit_ || g->orientation_deg != 0.0) o["orientation_deg"] = num(g->orientation_deg);
    } else if (const auto* r = std::get_if<RadialPhase>(&p)) {
      o["type"] = "radial_phase";
      o["normalization_radius"] = param(r->normalization_radius);
      if (edit_ || !r->coefficients.empty()) o["coefficients"] = param_list(r->coefficients);
    }
    return o;
  }

  OJson interaction(const Interaction& i) const {
    OJson o;
    if (std::holds_alternative<Fresnel>(i)) {
      o["type"] = "fresnel";
    } else if (std::holds_alternative<IdealMirror>(i)) {
      o["type"] = "ideal_mirror";
    } else if (std::holds_alternative<IdealAntiReflection>(i)) {
      o["type"] = "ideal_anti_reflection";
    } else if (std::holds_alternative<Absorber>(i)) {
      o["type"] = "absorber";
    } else if (const auto* b = std::get_if<IdealBeamSplitter>(&i)) {
      const IdealBeamSplitter d;
      o["type"] = "ideal_beam_splitter";
      if (edit_ || b->reflectance_s != d.reflectance_s) o["reflectance_s"] = num(b->reflectance_s);
      if (edit_ || b->reflectance_p != d.reflectance_p) o["reflectance_p"] = num(b->reflectance_p);
    } else if (const auto* c = std::get_if<CoatingRef>(&i)) {
      o["type"] = "coating";
      o["name"] = c->name;
    } else if (const auto* p = std::get_if<IdealPolarizer>(&i)) {
      const IdealPolarizer d;
      o["type"] = "ideal_polarizer";
      if (edit_ || p->transmission_axis != d.transmission_axis)
        o["transmission_axis"] = vec3(p->transmission_axis);
      if (edit_ || p->extinction_ratio != d.extinction_ratio)
        o["extinction_ratio"] = num(p->extinction_ratio);
    } else if (const auto* r = std::get_if<IdealRetarder>(&i)) {
      const IdealRetarder d;
      o["type"] = "ideal_retarder";
      if (edit_ || r->fast_axis != d.fast_axis) o["fast_axis"] = vec3(r->fast_axis);
      if (edit_ || r->retardance_waves != d.retardance_waves)
        o["retardance_waves"] = num(r->retardance_waves);
    }
    return o;
  }

  OJson surface(const Surface& s) const {
    OJson o;
    o["id"] = s.id.str();
    put_pose(o, s.pose);
    if (edit_ || s.shape != ShapeStack{}) o["shape"] = shape(s.shape);
    if (s.aperture) o["aperture"] = aperture(*s.aperture);
    if (edit_ || !s.phases.empty()) {
      OJson a = OJson::array();
      for (const PhaseLayer& p : s.phases) a.push_back(phase(p));
      o["phases"] = std::move(a);
    }
    if (edit_ || !std::holds_alternative<Fresnel>(s.interaction))
      o["interaction"] = interaction(s.interaction);
    return o;
  }

  OJson assembly(const Assembly& a) const {
    OJson o;
    o["type"] = "assembly";
    o["name"] = a.name;
    put_pose(o, a.pose);
    OJson children = OJson::array();
    for (const Node& c : a.children) children.push_back(node(c));
    o["children"] = std::move(children);
    return o;
  }

  OJson node(const Node& n) const {
    if (const auto* a = std::get_if<Assembly>(&n.value)) return assembly(*a);
    const auto& e = std::get<Element>(n.value);
    OJson o;
    o["type"] = std::string(enum_name(e.kind, kElementKinds));
    o["name"] = e.name;
    put_pose(o, e.pose);
    if (e.material && !e.segment_materials.empty()) {
      throw std::invalid_argument("element '" + e.name +
                                  "': material and segment_materials are both set");
    }
    if (e.material) {
      o["material"] = *e.material;
    } else if (!e.segment_materials.empty()) {
      OJson a = OJson::array();
      for (const std::string& m : e.segment_materials) a.push_back(m);
      o["material"] = std::move(a);
    }
    OJson surfaces = OJson::array();
    for (const Surface& s : e.surfaces) surfaces.push_back(surface(s));
    o["surfaces"] = std::move(surfaces);
    return o;
  }

  OJson path(const Path& p) const {
    OJson o;
    o["name"] = p.name;
    if (p.automatic) {
      o["events"] = "auto";
      return o;
    }
    OJson events = OJson::array();
    for (const Event& e : p.events) {
      OJson eo;
      eo["surface"] = e.surface.str();
      if (edit_ || e.kind != EventKind::Refract)
        eo["kind"] = std::string(enum_name(e.kind, kEventKinds));
      if (edit_ || e.order != 0) eo["order"] = e.order;
      events.push_back(std::move(eo));
    }
    o["events"] = std::move(events);
    return o;
  }

  bool edit_;
};

OJson write_tree(const System& s, detail::Form form) {
  return Writer(form).system(s);
}

}  // namespace

model::System parse_system(std::string_view json_text) {
  Json j;
  try {
    // Duplicate keys, syntax errors and number overflow (ADR 0008 addendum, ADR 0020).
    j = rtt::json::parse_strict(json_text);
  } catch (const rtt::json::StrictParseError& e) {
    throw ParseError(e.pointer(), e.message());
  }
  return read_system_tree(j);
}

model::System load_system(const std::filesystem::path& file) {
  const std::ifstream in(file, std::ios::binary);
  if (!in) throw std::runtime_error("cannot open '" + file.string() + "'");
  std::ostringstream buffer;
  buffer << in.rdbuf();
  return parse_system(buffer.str());
}

std::string to_json(const model::System& system) {
  return detail::format_canonical(write_tree(system, detail::Form::Canonical));
}

model::System detail::read_system(const nlohmann::json& j) {
  return read_system_tree(j);
}

nlohmann::ordered_json detail::write_system(const model::System& s, Form form) {
  return write_tree(s, form);
}

nlohmann::json detail::plain(const nlohmann::ordered_json& j) {
  switch (j.type()) {
    case nlohmann::ordered_json::value_t::object: {
      Json o = Json::object();
      for (const auto& item : j.items()) o[item.key()] = plain(item.value());
      return o;
    }
    case nlohmann::ordered_json::value_t::array: {
      Json a = Json::array();
      for (const auto& v : j) a.push_back(plain(v));
      return a;
    }
    case nlohmann::ordered_json::value_t::string:
      return j.get<std::string>();
    case nlohmann::ordered_json::value_t::boolean:
      return j.get<bool>();
    case nlohmann::ordered_json::value_t::number_integer:
      return j.get<std::int64_t>();
    case nlohmann::ordered_json::value_t::number_unsigned:
      return j.get<std::uint64_t>();
    case nlohmann::ordered_json::value_t::number_float:
      return j.get<double>();
    default:
      return nullptr;
  }
}

void save_system(const model::System& system, const std::filesystem::path& file) {
  const std::string text = to_json(system);
  std::ofstream out(file, std::ios::binary | std::ios::trunc);
  if (!out) throw std::runtime_error("cannot write '" + file.string() + "'");
  out << text;
  if (!out) throw std::runtime_error("write failed for '" + file.string() + "'");
}

}  // namespace rtt::io
