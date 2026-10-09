#include "rtt/io/json_io.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <initializer_list>
#include <limits>
#include <nlohmann/json.hpp>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "json_format.hpp"
#include "json_tree.hpp"
#include "rtt/diagnostics/codes.hpp"
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

const EnumTable<PoseReference> kPoseReferences = {
    {"absolute", PoseReference::Absolute},
    {"relative_to_preceding", PoseReference::RelativeToPreceding},
    {"relative_to_sibling", PoseReference::RelativeToSibling},
};

const EnumTable<PoseOrder> kPoseOrders = {
    {"translate_first", PoseOrder::TranslateFirst},
    {"rotate_first", PoseOrder::RotateFirst},
};

/// Operand types of the paraxial quantities (ADR 0030, point 3).
const EnumTable<FirstOrderQuantity> kFirstOrderTypes = {
    {"efl", FirstOrderQuantity::Efl},
    {"bfl", FirstOrderQuantity::Bfl},
    {"image_fnumber", FirstOrderQuantity::ImageFNumber},
    {"magnification", FirstOrderQuantity::Magnification},
};

const EnumTable<SpotReference> kSpotReferences = {
    {"centroid", SpotReference::Centroid},
    {"chief", SpotReference::Chief},
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
  [[nodiscard]] const std::string& pointer() const noexcept { return pointer_; }

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

/// Checks that `j` is an array (a check without return value, so that no reference to the
/// parameter is handed back; the callers keep using `j`).
void expect_array(const Json& j, const Ctx& c) {
  if (!j.is_array()) c.fail("expected an array, got " + std::string(type_name(j)));
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

/// Param (ADR 0029, point 3): a plain number, {"value": n, "variable": b, "min": n, "max": n}
/// (only value required) or a reference {"param": "NAME"} without anything else.
Param read_param(const Json& j, const Ctx& c) {
  if (j.is_number()) return Param(j.get<double>());
  if (!j.is_object()) c.fail("expected a number or a parameter object");
  if (const Json* name = find(j, "param")) {
    expect_object(j, c, {"param"});
    return Param::bound(read_string(*name, c.at("param")));
  }
  expect_object(j, c, {"value", "variable", "min", "max"});
  Param p(read_number(require(j, "value", c), c.at("value")));
  if (const Json* v = find(j, "variable")) p.variable = read_bool(*v, c.at("variable"));
  if (const Json* v = find(j, "min")) p.min = read_number(*v, c.at("min"));
  if (const Json* v = find(j, "max")) p.max = read_number(*v, c.at("max"));
  return p;
}

std::vector<Param> read_param_list(const Json& j, const Ctx& c) {
  std::vector<Param> out;
  expect_array(j, c);
  const Json& a = j;
  out.reserve(a.size());
  for (std::size_t i = 0; i < a.size(); ++i) out.push_back(read_param(a[i], c.at(i)));
  return out;
}

std::array<Param, 3> read_param3(const Json& j, const Ctx& c) {
  expect_array(j, c);
  const Json& a = j;
  if (a.size() != 3) c.fail("expected exactly 3 values");
  return {read_param(a[0], c.at(0)), read_param(a[1], c.at(1)), read_param(a[2], c.at(2))};
}

std::array<double, 3> read_vec3(const Json& j, const Ctx& c) {
  expect_array(j, c);
  const Json& a = j;
  if (a.size() != 3) c.fail("expected exactly 3 numbers");
  return {read_number(a[0], c.at(0)), read_number(a[1], c.at(1)), read_number(a[2], c.at(2))};
}

template <class T, class F>
void read_opt(const Json& j, std::string_view key, const Ctx& c, T& target, F reader) {
  if (const Json* v = find(j, key)) target = reader(*v, c.at(key));
}

Pose read_pose(const Json& j, const Ctx& c) {
  expect_object(j, c, {"reference", "order", "position", "rotation_deg", "pivot"});
  Pose p;
  if (const Json* r = find(j, "reference")) {
    p.reference = read_enum(*r, c.at("reference"), kPoseReferences);
  }
  if (const Json* o = find(j, "order")) p.order = read_enum(*o, c.at("order"), kPoseOrders);
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
    expect_array(*t, tc);
    const Json& a = *t;
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

/// Efficiency per diffraction order (ADR 0025): [{"order": m, "efficiency": eta}, ...]. Values
/// and duplicates are checked by validate.
std::vector<DiffractionEfficiency> read_efficiency(const Json& j, const Ctx& c) {
  expect_array(j, c);
  const Json& a = j;
  std::vector<DiffractionEfficiency> out;
  out.reserve(a.size());
  for (std::size_t i = 0; i < a.size(); ++i) {
    const Ctx ci = c.at(i);
    expect_object(a[i], ci, {"order", "efficiency"});
    out.push_back({read_int(require(a[i], "order", ci), ci.at("order")),
                   read_number(require(a[i], "efficiency", ci), ci.at("efficiency"))});
  }
  return out;
}

Surface read_surface(const Json& j, const Ctx& c) {
  expect_object(
      j, c, {"id", "pose", "shape", "aperture", "phases", "diffraction_efficiency", "interaction"});
  Surface s;
  s.id = SurfaceId(read_string(require(j, "id", c), c.at("id")));
  read_opt(j, "pose", c, s.pose, read_pose);
  read_opt(j, "shape", c, s.shape, read_shape);
  if (const Json* a = find(j, "aperture")) s.aperture = read_aperture(*a, c.at("aperture"));
  if (const Json* p = find(j, "phases")) {
    const Ctx pc = c.at("phases");
    expect_array(*p, pc);
    const Json& a = *p;
    for (std::size_t i = 0; i < a.size(); ++i)
      s.phases.push_back(read_phase(a[i], c.at("phases").at(i)));
  }
  if (const Json* d = find(j, "diffraction_efficiency")) {
    s.diffraction_efficiency = read_efficiency(*d, c.at("diffraction_efficiency"));
  }
  read_opt(j, "interaction", c, s.interaction, read_interaction);
  return s;
}

/// Element material: a string (shorthand, one material for all segments), a non-empty
/// array of strings (one material per segment, ADR 0017) or a crystal object
/// {"ordinary": ..., "extraordinary": ...} (ADR 0026).
void read_material(const Json& j, const Ctx& c, Element& e) {
  if (j.is_string()) {
    e.material = j.get<std::string>();
    return;
  }
  if (j.is_object()) {
    expect_object(j, c, {"ordinary", "extraordinary"});
    e.crystal = CrystalMaterial{read_string(require(j, "ordinary", c), c.at("ordinary")),
                                read_string(require(j, "extraordinary", c), c.at("extraordinary"))};
    return;
  }
  if (!j.is_array()) {
    c.fail("expected a string, an array of strings or a crystal object, got " +
           std::string(type_name(j)));
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
  expect_object(j, c, {"type", "name", "pose", "material", "optic_axis", "surfaces"});
  Element e;
  e.kind = read_enum(require(j, "type", c), c.at("type"), kElementKinds);
  e.name = read_string(require(j, "name", c), c.at("name"));
  read_opt(j, "pose", c, e.pose, read_pose);
  if (const Json* m = find(j, "material")) read_material(*m, c.at("material"), e);
  if (const Json* a = find(j, "optic_axis")) e.optic_axis = read_vec3(*a, c.at("optic_axis"));
  const Json& surfaces = require_array(j, "surfaces", c);
  e.surfaces.reserve(surfaces.size());
  for (std::size_t i = 0; i < surfaces.size(); ++i) {
    e.surfaces.push_back(read_surface(surfaces[i], c.at("surfaces").at(i)));
  }
  return Node{std::move(e)};
}

/// `before_0_3`: the file predates schema 0.3, where the event kind "diffract" existed.
Path read_path(const Json& j, const Ctx& c, bool before_0_3) {
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
  expect_array(events, ec);
  const Json& a = events;
  for (std::size_t i = 0; i < a.size(); ++i) {
    const Ctx c2 = ec.at(i);
    expect_object(a[i], c2, {"surface", "kind", "order"});
    Event e;
    e.surface = SurfaceId(read_string(require(a[i], "surface", c2), c2.at("surface")));
    if (const Json* k = find(a[i], "kind")) {
      // Migration to 0.3 (ADR 0025, point 4): "diffract" kept the medium like "transmit" and is
      // read as "transmit" with the same order; the media of the path stay the same.
      if (before_0_3 && k->is_string() && k->get<std::string>() == "diffract") {
        e.kind = EventKind::Transmit;
      } else {
        e.kind = read_enum(*k, c2.at("kind"), kEventKinds);
      }
    }
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

/// Material forms newer than the file: 0.2 adds the array form of "material" (ADR 0017), 0.3
/// the crystal object (ADR 0026). Such a form in a file that claims an older version is an
/// error, so that a version is never silently upgraded.
void check_material_forms(const Json& node, const Ctx& c, bool lists_allowed) {
  if (!node.is_object()) return;  // structural errors are reported by the reader
  if (const Json* m = find(node, "material")) {
    if (m->is_array() && !lists_allowed) {
      c.at("material").fail("material lists need schema_version 0.2 or later");
    }
    if (m->is_object()) c.at("material").fail("crystal materials need schema_version 0.3 or later");
  }
  if (const Json* children = find(node, "children"); children != nullptr && children->is_array()) {
    for (std::size_t i = 0; i < children->size(); ++i) {
      check_material_forms((*children)[i], c.at("children").at(i), lists_allowed);
    }
  }
}

/// Which migrations a file needs.
struct Migration {
  bool before_0_3 = false;  ///< the event kind "diffract" existed
  bool before_0_4 = false;  ///< Params could carry a pickup; no 0.4 forms
};

/// Checks the schema version of `j` and tells which migrations to the current version it needs.
/// Pre-1.0 rule: major and minor must match a supported version, patch may differ. The reader
/// then reads `j` (after migrate_0_4 for files before 0.4) as a current file and the writer
/// writes kSchemaVersion.
/// - 0.1 -> 0.2 (ADR 0017): content unchanged, only material lists are new.
/// - 0.2 -> 0.3 (ADR 0025, 0026): the event kind "diffract" is read as "transmit" (read_path);
///   crystals, the optic axis and diffraction efficiencies are new.
/// - 0.3 -> 0.4 (ADR 0028, 0029): the pickup of a Param is dropped with io.pickup_dropped;
///   Pose.reference/order, the parameter table, configurations, bounds and the merit function
///   (ADR 0030) are new.
Migration migrate(const Json& j, const std::string& version, const Ctx& c) {
  const std::string ours(kSchemaVersion);
  const std::string mm = major_minor(version);
  if (!mm.empty() && mm == major_minor(ours)) return {};
  if (mm == "0.1" || mm == "0.2") {
    if (const Json* root = find(j, "root")) {
      check_material_forms(*root, Ctx().at("root"), mm == "0.2");
    }
    return {true, true};
  }
  if (mm == "0.3") return {false, true};
  c.fail("incompatible schema_version '" + version + "', this build reads " + ours +
         " and migrates 0.1, 0.2 and 0.3");
}

/// A Param object: no "type" (the system aperture {"type", "value"} has one) and "value" or
/// "param".
bool param_object(const Json& j) {
  return j.is_object() && !j.contains("type") && (j.contains("value") || j.contains("param"));
}

/// Files before 0.4 (ADR 0029, point 6): drops every pickup, keeping value and variable, with
/// the warning io.pickup_dropped at its pointer in the file; a 0.4 form (Pose.reference or
/// .order, a bound Param or bounds, the parameter table, configurations or the merit function) is
/// an error, so that a version is never silently upgraded.
void migrate_0_4(Json& j, const Ctx& c, std::vector<Diagnostic>* warnings, bool top) {
  const auto newer = [](const Ctx& at, std::string_view what) {
    at.fail(std::string(what) + " needs schema_version 0.4 or later");
  };
  if (j.is_array()) {
    for (std::size_t i = 0; i < j.size(); ++i) migrate_0_4(j[i], c.at(i), warnings, false);
    return;
  }
  if (!j.is_object()) return;
  if (top) {
    if (j.contains("configurations")) newer(c.at("configurations"), "configurations");
    if (j.contains("parameters")) newer(c.at("parameters"), "a parameter table");
    if (j.contains("optimization")) newer(c.at("optimization"), "a merit function");
  }
  if (param_object(j)) {
    for (const char* key : {"param", "min", "max"}) {
      if (j.contains(key)) newer(c.at(key), "'" + std::string(key) + "' at a parameter");
    }
    if (const auto it = j.find("pickup"); it != j.end()) {
      const std::string text = it->is_string() ? it->get<std::string>() : it->dump();
      if (warnings != nullptr) {
        constexpr diagnostics::DiagnosticCode kCode = "io.pickup_dropped";
        warnings->push_back(
            {kCode.severity(), c.at("pickup").pointer(),
             "pickup '" + text + "' dropped: schema 0.4 has a parameter table instead (ADR 0029)",
             std::string(kCode.str())});
      }
      j.erase(it);
    }
    return;
  }
  if (const auto pose = j.find("pose"); pose != j.end() && pose->is_object()) {
    for (const char* key : {"reference", "order"}) {
      if (pose->contains(key)) newer(c.at("pose").at(key), "'" + std::string(key) + "' of a pose");
    }
  }
  for (auto& [key, value] : j.items()) migrate_0_4(value, c.at(key), warnings, false);
}

/// Rows of the parameter table (ADR 0029, point 1): a name and exactly one of value, values,
/// expression; variable at every row; min and max only at value and values rows.
ParameterRow read_parameter_row(const Json& j, const Ctx& c) {
  expect_object(j, c, {"name", "value", "values", "expression", "variable", "min", "max"});
  ParameterRow r;
  r.name = read_string(require(j, "name", c), c.at("name"));
  const Json* value = find(j, "value");
  const Json* values = find(j, "values");
  const Json* expression = find(j, "expression");
  const int forms =
      (value != nullptr ? 1 : 0) + (values != nullptr ? 1 : 0) + (expression != nullptr ? 1 : 0);
  if (forms != 1) c.fail("a parameter row needs exactly one of value, values, expression");
  if (value != nullptr) {
    r.form = read_number(*value, c.at("value"));
  } else if (values != nullptr) {
    const Ctx vc = c.at("values");
    expect_array(*values, vc);
    std::vector<double> v;
    v.reserve(values->size());
    for (std::size_t k = 0; k < values->size(); ++k)
      v.push_back(read_number((*values)[k], vc.at(k)));
    r.form = std::move(v);
  } else {
    for (const char* key : {"min", "max"}) {
      if (find(j, key) != nullptr) {
        c.at(key).fail("bounds are only allowed at rows with value or values");
      }
    }
    r.form = ParameterExpression{read_string(*expression, c.at("expression"))};
  }
  if (const Json* v = find(j, "variable")) r.variable = read_bool(*v, c.at("variable"));
  if (const Json* v = find(j, "min")) r.min = read_number(*v, c.at("min"));
  if (const Json* v = find(j, "max")) r.max = read_number(*v, c.at("max"));
  return r;
}

/// Index into the fields or wavelengths of the system: an integer 0 ... 65535 (the range against
/// the system is validate's, merit.index_out_of_range).
std::uint16_t read_index(const Json& j, const Ctx& c) {
  const int v = read_int(j, c);
  if (v < 0 || v > std::numeric_limits<std::uint16_t>::max()) {
    c.fail("index " + std::to_string(v) + " out of range (0 ... 65535)");
  }
  return static_cast<std::uint16_t>(v);
}

/// Fields or wavelengths of a generator: a non-empty list of indices (missing means all).
std::vector<std::uint16_t> read_selection(const Json& j, const Ctx& c) {
  expect_array(j, c);
  if (j.empty()) c.fail("a selection is not empty; omit it for all");
  std::vector<std::uint16_t> out;
  out.reserve(j.size());
  for (std::size_t k = 0; k < j.size(); ++k) out.push_back(read_index(j[k], c.at(k)));
  return out;
}

std::optional<std::string> read_configuration(const Json& j, const Ctx& c) {
  if (const Json* v = find(j, "configuration")) return read_string(*v, c.at("configuration"));
  return std::nullopt;
}

/// target (required), weight and configuration of an operand (ADR 0030, point 2).
OperandCommon read_common(const Json& j, const Ctx& c) {
  OperandCommon o;
  o.target = read_number(require(j, "target", c), c.at("target"));
  read_opt(j, "weight", c, o.weight, read_number);
  o.configuration = read_configuration(j, c);
  return o;
}

std::string read_path_name(const Json& j, const Ctx& c) {
  return read_string(require(j, "path", c), c.at("path"));
}

/// One operand (ADR 0030, point 3); keys that do not belong to its type are errors.
Operand read_operand(const Json& j, const Ctx& c) {
  const std::string type = read_type(j, c);
  for (const auto& [name, quantity] : kFirstOrderTypes) {
    if (type != name) continue;
    expect_object(j, c, {"type", "path", "configuration", "wavelength", "target", "weight"});
    FirstOrderOperand o;
    o.common = read_common(j, c);
    o.quantity = quantity;
    o.path = read_path_name(j, c);
    read_opt(j, "wavelength", c, o.wavelength, read_index);
    return o;
  }
  if (type == "ray_x" || type == "ray_y") {
    expect_object(j, c,
                  {"type", "path", "configuration", "surface", "occurrence", "field", "px", "py",
                   "wavelength", "target", "weight"});
    RayOperand o;
    o.common = read_common(j, c);
    o.coordinate = type == "ray_x" ? RayCoordinate::X : RayCoordinate::Y;
    o.path = read_path_name(j, c);
    o.surface = SurfaceId(read_string(require(j, "surface", c), c.at("surface")));
    if (const Json* v = find(j, "occurrence")) {
      const int n = read_int(*v, c.at("occurrence"));
      if (n < 0) c.at("occurrence").fail("occurrence is 0-based and not negative");
      o.occurrence = static_cast<std::uint32_t>(n);
    }
    read_opt(j, "field", c, o.field, read_index);
    read_opt(j, "px", c, o.px, read_number);
    read_opt(j, "py", c, o.py, read_number);
    read_opt(j, "wavelength", c, o.wavelength, read_index);
    return o;
  }
  if (type == "spot_rms") {
    expect_object(j, c,
                  {"type", "path", "configuration", "field", "wavelength", "polychromatic",
                   "reference", "rings", "target", "weight"});
    SpotRmsOperand o;
    o.common = read_common(j, c);
    o.path = read_path_name(j, c);
    read_opt(j, "field", c, o.field, read_index);
    read_opt(j, "polychromatic", c, o.polychromatic, read_bool);
    if (o.polychromatic && find(j, "wavelength") != nullptr) {
      c.at("wavelength").fail("a polychromatic spot uses all wavelengths and has no wavelength");
    }
    read_opt(j, "wavelength", c, o.wavelength, read_index);
    if (const Json* r = find(j, "reference")) {
      o.reference = read_enum(*r, c.at("reference"), kSpotReferences);
    }
    read_opt(j, "rings", c, o.rings, read_int);
    return o;
  }
  if (type == "opd_rms") {
    expect_object(
        j, c, {"type", "path", "configuration", "field", "wavelength", "grid", "target", "weight"});
    OpdRmsOperand o;
    o.common = read_common(j, c);
    o.path = read_path_name(j, c);
    read_opt(j, "field", c, o.field, read_index);
    read_opt(j, "wavelength", c, o.wavelength, read_index);
    read_opt(j, "grid", c, o.grid, read_int);
    return o;
  }
  if (type == "param_value") {
    expect_object(j, c, {"type", "parameter", "configuration", "target", "weight"});
    ParamValueOperand o;
    o.common = read_common(j, c);
    o.parameter = read_string(require(j, "parameter", c), c.at("parameter"));
    return o;
  }
  c.at("type").fail("unknown operand type '" + type +
                    "', expected efl, bfl, image_fnumber, magnification, ray_x, ray_y, spot_rms, "
                    "opd_rms or param_value");
}

/// The keys every generator has (ADR 0030, point 4); no target (always 0).
template <class G>
void read_generator_common(const Json& j, const Ctx& c, G& g) {
  g.path = read_path_name(j, c);
  g.configuration = read_configuration(j, c);
  read_opt(j, "fields", c, g.fields, read_selection);
  read_opt(j, "wavelengths", c, g.wavelengths, read_selection);
  read_opt(j, "rings", c, g.rings, read_int);
  read_opt(j, "arms", c, g.arms, read_int);
  read_opt(j, "weight", c, g.weight, read_number);
}

Generator read_generator(const Json& j, const Ctx& c) {
  const std::string type = read_type(j, c);
  if (type == "rms_spot") {
    expect_object(j, c,
                  {"type", "path", "configuration", "fields", "wavelengths", "reference", "rings",
                   "arms", "weight"});
    SpotGenerator g;
    read_generator_common(j, c, g);
    if (const Json* r = find(j, "reference")) {
      g.reference = read_enum(*r, c.at("reference"), kSpotReferences);
    }
    return g;
  }
  if (type == "rms_wavefront") {
    expect_object(
        j, c,
        {"type", "path", "configuration", "fields", "wavelengths", "rings", "arms", "weight"});
    WavefrontGenerator g;
    read_generator_common(j, c, g);
    return g;
  }
  c.at("type").fail("unknown generator type '" + type + "', expected rms_spot or rms_wavefront");
}

/// The section "optimization" (ADR 0030, point 2): both lists optional.
Optimization read_optimization(const Json& j, const Ctx& c) {
  expect_object(j, c, {"operands", "generators"});
  Optimization o;
  if (const Json* a = find(j, "operands")) {
    const Ctx oc = c.at("operands");
    expect_array(*a, oc);
    for (std::size_t i = 0; i < a->size(); ++i)
      o.operands.push_back(read_operand((*a)[i], oc.at(i)));
  }
  if (const Json* a = find(j, "generators")) {
    const Ctx gc = c.at("generators");
    expect_array(*a, gc);
    for (std::size_t i = 0; i < a->size(); ++i) {
      o.generators.push_back(read_generator((*a)[i], gc.at(i)));
    }
  }
  return o;
}

System read_migrated(const Json& j, const Migration& m);

/// `warnings` may be null: the warnings of the migration are dropped.
System read_system_tree(const Json& j, std::vector<Diagnostic>* warnings) {
  const Ctx c;
  if (!j.is_object()) c.fail("expected an object, got " + std::string(type_name(j)));
  const std::string version = read_string(require(j, "schema_version", c), c.at("schema_version"));
  const Migration m = migrate(j, version, c.at("schema_version"));
  if (!m.before_0_4) return read_migrated(j, m);
  Json current = j;
  migrate_0_4(current, c, warnings, true);
  return read_migrated(current, m);
}

/// Reads `j` as a current file; `m` tells the reader which older forms to accept.
System read_migrated(const Json& j, const Migration& m) {
  const Ctx c;
  const bool before_0_3 = m.before_0_3;
  expect_object(
      j, c,
      {"schema_version", "name", "units", "environment", "object", "wavelengths", "aperture",
       "fields", "configurations", "parameters", "root", "paths", "optimization"});
  System s;
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

  if (const Json* a = find(j, "configurations")) {
    const Ctx cc = c.at("configurations");
    expect_array(*a, cc);
    if (a->empty()) cc.fail("a configurations section is not empty; omit it for one configuration");
    for (std::size_t k = 0; k < a->size(); ++k) {
      const Ctx kc = cc.at(k);
      expect_object((*a)[k], kc, {"name"});
      s.configurations.push_back({read_string(require((*a)[k], "name", kc), kc.at("name"))});
    }
  }

  if (const Json* a = find(j, "parameters")) {
    const Ctx pc = c.at("parameters");
    expect_array(*a, pc);
    s.parameters.reserve(a->size());
    for (std::size_t i = 0; i < a->size(); ++i) {
      s.parameters.push_back(read_parameter_row((*a)[i], pc.at(i)));
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
    for (std::size_t i = 0; i < a.size(); ++i) {
      s.paths.push_back(read_path(a[i], pc.at(i), before_0_3));
    }
  }

  if (const Json* o = find(j, "optimization")) {
    s.optimization = read_optimization(*o, c.at("optimization"));
  }
  return s;
}

// ===================================================================== writer =====
//
// Two forms (json_tree.hpp): the canonical form omits values equal to their defaults and writes
// a plain Param as a number; the edit form (ADR 0024) writes every value and every Param as an
// object. Optional members without a value (surface aperture, diffraction efficiency, element
// material, optic axis, bounds) are missing in both. A bound Param is {"param": …} in both, without
// a value; the edit form writes configurations only if there are any (an empty section is not
// readable, ADR 0029) and parameters and the merit function always (both lists, so that a patch
// can add the first entry). Operands always carry their target (ADR 0030, point 2).

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

    if (!s.configurations.empty()) {
      OJson cs = OJson::array();
      for (const Configuration& c : s.configurations) cs.push_back(OJson{{"name", c.name}});
      o["configurations"] = std::move(cs);
    }
    if (edit_ || !s.parameters.empty()) {
      OJson rows = OJson::array();
      for (const ParameterRow& r : s.parameters) rows.push_back(parameter_row(r));
      o["parameters"] = std::move(rows);
    }

    o["root"] = assembly(s.root);

    OJson paths = OJson::array();
    for (const Path& p : s.paths) paths.push_back(path(p));
    o["paths"] = std::move(paths);

    if (edit_ || !s.optimization.empty()) o["optimization"] = optimization(s.optimization);
    return o;
  }

 private:
  OJson param(const Param& p) const {
    // A bound Param has no value of its own (ADR 0029, point 3): only the reference.
    if (const std::optional<std::string>& name = p.param) return OJson{{"param", *name}};
    if (!edit_ && p.is_plain()) return num(p.value);
    OJson o;
    o["value"] = num(p.value);
    if (edit_ || p.variable) o["variable"] = p.variable;
    if (p.min) o["min"] = num(*p.min);
    if (p.max) o["max"] = num(*p.max);
    return o;
  }

  OJson parameter_row(const ParameterRow& r) const {
    OJson o;
    o["name"] = r.name;
    if (const auto* value = std::get_if<double>(&r.form)) {
      o["value"] = num(*value);
    } else if (const auto* values = std::get_if<std::vector<double>>(&r.form)) {
      OJson a = OJson::array();
      for (const double v : *values) a.push_back(num(v));
      o["values"] = std::move(a);
    } else {
      o["expression"] = std::get<ParameterExpression>(r.form).text;
    }
    if (edit_ || r.variable) o["variable"] = r.variable;
    if (r.min) o["min"] = num(*r.min);
    if (r.max) o["max"] = num(*r.max);
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
    if (edit_ || p.reference != d.reference)
      j["reference"] = std::string(enum_name(p.reference, kPoseReferences));
    if (edit_ || p.order != d.order) j["order"] = std::string(enum_name(p.order, kPoseOrders));
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
    if (s.diffraction_efficiency) {
      OJson a = OJson::array();
      for (const DiffractionEfficiency& d : *s.diffraction_efficiency) {
        OJson eo;
        eo["order"] = d.order;
        eo["efficiency"] = num(d.efficiency);
        a.push_back(std::move(eo));
      }
      o["diffraction_efficiency"] = std::move(a);
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
    if (e.crystal && (e.material || !e.segment_materials.empty())) {
      throw std::invalid_argument("element '" + e.name +
                                  "': crystal and isotropic material are both set");
    }
    if (e.material) {
      o["material"] = *e.material;
    } else if (!e.segment_materials.empty()) {
      OJson a = OJson::array();
      for (const std::string& m : e.segment_materials) a.push_back(m);
      o["material"] = std::move(a);
    } else if (e.crystal) {
      OJson m;
      m["ordinary"] = e.crystal->ordinary;
      m["extraordinary"] = e.crystal->extraordinary;
      o["material"] = std::move(m);
    }
    if (e.optic_axis) o["optic_axis"] = vec3(*e.optic_axis);
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

  // ---- merit function (ADR 0030): keys in the order type, path or parameter, configuration,
  // the keys of the type, target, weight.

  OJson optimization(const Optimization& m) const {
    OJson o;
    if (edit_ || !m.operands.empty()) {
      OJson a = OJson::array();
      for (const Operand& op : m.operands) {
        a.push_back(std::visit([this](const auto& x) { return operand(x); }, op));
      }
      o["operands"] = std::move(a);
    }
    if (edit_ || !m.generators.empty()) {
      OJson a = OJson::array();
      for (const Generator& g : m.generators) {
        a.push_back(std::visit([this](const auto& x) { return generator(x); }, g));
      }
      o["generators"] = std::move(a);
    }
    return o;
  }

  static void put_configuration(OJson& o, const std::optional<std::string>& configuration) {
    if (configuration) o["configuration"] = *configuration;
  }

  void put_target_weight(OJson& o, const OperandCommon& c) const {
    o["target"] = num(c.target);
    if (edit_ || c.weight != OperandCommon{}.weight) o["weight"] = num(c.weight);
  }

  static void put_index(OJson& o, const char* key, const std::optional<std::uint16_t>& index) {
    if (index) o[key] = *index;
  }

  static OJson selection(const std::vector<std::uint16_t>& list) {
    OJson a = OJson::array();
    for (const std::uint16_t k : list) a.push_back(k);
    return a;
  }

  OJson operand(const FirstOrderOperand& op) const {
    OJson o;
    o["type"] = std::string(enum_name(op.quantity, kFirstOrderTypes));
    o["path"] = op.path;
    put_configuration(o, op.common.configuration);
    put_index(o, "wavelength", op.wavelength);
    put_target_weight(o, op.common);
    return o;
  }

  OJson operand(const RayOperand& op) const {
    const RayOperand d;
    OJson o;
    o["type"] = op.coordinate == RayCoordinate::X ? "ray_x" : "ray_y";
    o["path"] = op.path;
    put_configuration(o, op.common.configuration);
    o["surface"] = op.surface.str();
    if (op.occurrence) o["occurrence"] = *op.occurrence;
    if (edit_ || op.field != d.field) o["field"] = op.field;
    if (edit_ || op.px != d.px) o["px"] = num(op.px);
    if (edit_ || op.py != d.py) o["py"] = num(op.py);
    put_index(o, "wavelength", op.wavelength);
    put_target_weight(o, op.common);
    return o;
  }

  OJson operand(const SpotRmsOperand& op) const {
    const SpotRmsOperand d;
    OJson o;
    o["type"] = "spot_rms";
    o["path"] = op.path;
    put_configuration(o, op.common.configuration);
    if (edit_ || op.field != d.field) o["field"] = op.field;
    put_index(o, "wavelength", op.wavelength);
    if (edit_ || op.polychromatic != d.polychromatic) o["polychromatic"] = op.polychromatic;
    if (edit_ || op.reference != d.reference) {
      o["reference"] = std::string(enum_name(op.reference, kSpotReferences));
    }
    if (edit_ || op.rings != d.rings) o["rings"] = op.rings;
    put_target_weight(o, op.common);
    return o;
  }

  OJson operand(const OpdRmsOperand& op) const {
    const OpdRmsOperand d;
    OJson o;
    o["type"] = "opd_rms";
    o["path"] = op.path;
    put_configuration(o, op.common.configuration);
    if (edit_ || op.field != d.field) o["field"] = op.field;
    put_index(o, "wavelength", op.wavelength);
    if (edit_ || op.grid != d.grid) o["grid"] = op.grid;
    put_target_weight(o, op.common);
    return o;
  }

  OJson operand(const ParamValueOperand& op) const {
    OJson o;
    o["type"] = "param_value";
    o["parameter"] = op.parameter;
    put_configuration(o, op.common.configuration);
    put_target_weight(o, op.common);
    return o;
  }

  template <class G>
  void put_generator_sampling(OJson& o, const G& g) const {
    const G d;
    if (edit_ || g.rings != d.rings) o["rings"] = g.rings;
    if (edit_ || g.arms != d.arms) o["arms"] = g.arms;
    if (edit_ || g.weight != d.weight) o["weight"] = num(g.weight);
  }

  template <class G>
  static void put_generator_head(OJson& o, const char* type, const G& g) {
    o["type"] = type;
    o["path"] = g.path;
    put_configuration(o, g.configuration);
    if (g.fields) o["fields"] = selection(*g.fields);
    if (g.wavelengths) o["wavelengths"] = selection(*g.wavelengths);
  }

  OJson generator(const SpotGenerator& g) const {
    OJson o;
    put_generator_head(o, "rms_spot", g);
    if (edit_ || g.reference != SpotGenerator{}.reference) {
      o["reference"] = std::string(enum_name(g.reference, kSpotReferences));
    }
    put_generator_sampling(o, g);
    return o;
  }

  OJson generator(const WavefrontGenerator& g) const {
    OJson o;
    put_generator_head(o, "rms_wavefront", g);
    put_generator_sampling(o, g);
    return o;
  }

  bool edit_;
};

OJson write_tree(const System& s, detail::Form form) {
  return Writer(form).system(s);
}

}  // namespace

namespace {

System parse_with(std::string_view json_text, std::vector<Diagnostic>* warnings) {
  Json j;
  try {
    // Duplicate keys, syntax errors and number overflow (ADR 0008 addendum, ADR 0020).
    j = rtt::json::parse_strict(json_text);
  } catch (const rtt::json::StrictParseError& e) {
    throw ParseError(e.pointer(), e.message());
  }
  return read_system_tree(j, warnings);
}

std::string read_file(const std::filesystem::path& file) {
  const std::ifstream in(file, std::ios::binary);
  if (!in) throw std::runtime_error("cannot open '" + file.string() + "'");
  std::ostringstream buffer;
  buffer << in.rdbuf();
  return buffer.str();
}

}  // namespace

model::System parse_system(std::string_view json_text, std::vector<model::Diagnostic>& warnings) {
  // Appended only when reading succeeds, so that `warnings` is unchanged if it throws.
  std::vector<Diagnostic> found;
  System s = parse_with(json_text, &found);
  warnings.insert(warnings.end(), found.begin(), found.end());
  return s;
}

model::System parse_system(std::string_view json_text) {
  return parse_with(json_text, nullptr);
}

model::System load_system(const std::filesystem::path& file,
                          std::vector<model::Diagnostic>& warnings) {
  return parse_system(read_file(file), warnings);
}

model::System load_system(const std::filesystem::path& file) {
  return parse_system(read_file(file));
}

std::string to_json(const model::System& system) {
  return detail::format_canonical(write_tree(system, detail::Form::Canonical));
}

model::System detail::read_system(const nlohmann::json& j) {
  return read_system_tree(j, nullptr);
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
