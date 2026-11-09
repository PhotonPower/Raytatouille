#include "rtt/io/surface_table.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

#include "export_path.hpp"
#include "numbers.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/element.hpp"
#include "rtt/model/parameters.hpp"
#include "rtt/model/path.hpp"
#include "rtt/model/surface.hpp"

namespace rtt::io {

TableError::TableError(std::string code,
                       std::size_t line,
                       std::size_t column,
                       const std::string& message)
    : std::runtime_error((line == 0 ? std::string()
                                    : "line " + std::to_string(line) + ", column " +
                                          std::to_string(column) + ": ") +
                         message + " [" + code + "]"),
      code_(std::move(code)),
      line_(line),
      column_(column) {}

namespace {

using detail::export_fail;
using detail::ExportRow;
using model::ElementKind;

/// The columns of ADR 0045, point 7, in the order they are written.
constexpr std::array<std::string_view, 22> kColumns{"row",          "surface",
                                                    "element",      "kind",
                                                    "shape",        "radius",
                                                    "conic",        "a4",
                                                    "a6",           "a8",
                                                    "a10",          "a12",
                                                    "a14",          "a16",
                                                    "thickness",    "material",
                                                    "aperture",     "aperture_1",
                                                    "aperture_2",   "stop",
                                                    "focal_length", "object_distance"};

enum Column : std::uint8_t {
  kRow,
  kSurface,
  kElement,
  kKind,
  kShape,
  kRadius,
  kConic,
  kA4,  // a4 ... a16 follow
  kThickness = kA4 + 7,
  kMaterial,
  kAperture,
  kAperture1,
  kAperture2,
  kStop,
  kFocalLength,
  kObjectDistance,
};

constexpr std::array<std::pair<ElementKind, std::string_view>, 6> kKinds{{
    {ElementKind::Lens, "lens"},
    {ElementKind::Plate, "plate"},
    {ElementKind::Mirror, "mirror"},
    {ElementKind::ThinElement, "thin_element"},
    {ElementKind::Stop, "stop"},
    {ElementKind::Detector, "detector"},
}};

std::string_view kind_name(ElementKind k) {
  for (const auto& [kind, name] : kKinds) {
    if (kind == k) return name;
  }
  return "lens";
}

// ------------------------------------------------------------------------------- writing -----

/// A field after RFC 4180, rules 6 and 7: quoted if it holds a comma, a quote or a line break,
/// a quote inside doubled.
std::string quoted(const std::string& field) {
  if (field.find_first_of(",\"\r\n") == std::string::npos) return field;
  std::string out = "\"";
  for (const char c : field) {
    out += c;
    if (c == '"') out += '"';
  }
  out += '"';
  return out;
}

/// A number that is written empty when it is 0 (conic, a4 ... a16, inner radius; point 7).
std::string optional_number(double v) {
  return v == 0.0 ? std::string() : detail::format_number(v);
}

void write_shape(const ExportRow& row, std::array<std::string, kColumns.size()>& f) {
  const std::string& id = row.surface.id.str();
  const model::ShapeStack& shape = row.surface.shape;
  if (!shape.terms.empty()) {
    export_fail("export.surface_unsupported", id, "shape terms have no column in the table");
  }
  if (std::holds_alternative<model::Plane>(shape.base)) {
    f[kShape] = "plane";
  } else if (const auto* c = std::get_if<model::Conic>(&shape.base)) {
    f[kShape] = "conic";
    f[kRadius] = detail::format_number(c->radius.value);
    f[kConic] = optional_number(c->conic.value);
  } else if (const auto* e = std::get_if<model::EvenAsphere>(&shape.base)) {
    f[kShape] = "even_asphere";
    f[kRadius] = detail::format_number(e->radius.value);
    f[kConic] = optional_number(e->conic.value);
    for (std::size_t i = 0; i < e->coefficients.size(); ++i) {
      if (i < 7) {
        f[kA4 + i] = optional_number(e->coefficients[i].value);
      } else if (e->coefficients[i].value != 0.0) {
        export_fail("export.surface_unsupported", id, "the table holds the terms up to A16 only");
      }
    }
  } else {
    // an odd asphere too, also with c = 0 (decision G1 on #318)
    export_fail("export.surface_unsupported", id,
                "the table holds plane, conic and even_asphere surfaces only");
  }
}

void write_aperture(const ExportRow& row, std::array<std::string, kColumns.size()>& f) {
  const std::optional<model::Aperture>& a = row.surface.aperture;
  if (!a.has_value()) return;
  if (const auto* c = std::get_if<model::CircularAperture>(&a.value())) {
    f[kAperture] = "circular";
    f[kAperture1] = detail::format_number(c->radius);
    f[kAperture2] = optional_number(c->inner_radius);
  } else if (const auto* r = std::get_if<model::RectangularAperture>(&a.value())) {
    f[kAperture] = "rectangular";
    f[kAperture1] = detail::format_number(r->half_width_x);
    f[kAperture2] = detail::format_number(r->half_width_y);
  } else if (const auto* e = std::get_if<model::EllipticalAperture>(&a.value())) {
    f[kAperture] = "elliptical";
    f[kAperture1] = detail::format_number(e->semi_axis_x);
    f[kAperture2] = detail::format_number(e->semi_axis_y);
  }
}

void append_record(std::string& text, const std::array<std::string, kColumns.size()>& f) {
  for (std::size_t i = 0; i < f.size(); ++i) {
    if (i > 0) text += ',';
    text += quoted(f[i]);
  }
  text += "\r\n";  // RFC 4180, rule 1
}

// ------------------------------------------------------------------------------- reading -----

[[noreturn]] void table_fail(diagnostics::DiagnosticCode code,
                             std::size_t line,
                             std::size_t column,
                             const std::string& message) {
  throw TableError(std::string(code.str()), line, column, message);
}

struct Field {
  std::string text;
  std::size_t line = 0;    ///< 1-based line where the field starts
  std::size_t column = 0;  ///< 1-based number of the field in its record
};

using Record = std::vector<Field>;

/// True if `s` is well-formed UTF-8 (no overlong forms, no surrogates, at most U+10FFFF).
bool valid_utf8(std::string_view s) {
  std::size_t i = 0;
  while (i < s.size()) {
    const auto b = static_cast<unsigned char>(s[i]);
    std::size_t n = 0;
    std::uint32_t cp = 0;
    if (b < 0x80) {
      ++i;
      continue;
    }
    if (b >= 0xC2 && b <= 0xDF) {
      n = 1;
      cp = b & 0x1Fu;
    } else if (b >= 0xE0 && b <= 0xEF) {
      n = 2;
      cp = b & 0x0Fu;
    } else if (b >= 0xF0 && b <= 0xF4) {
      n = 3;
      cp = b & 0x07u;
    } else {
      return false;
    }
    for (std::size_t k = 1; k <= n; ++k) {
      if (i + k >= s.size()) return false;
      const auto c = static_cast<unsigned char>(s[i + k]);
      if ((c & 0xC0u) != 0x80u) return false;
      cp = (cp << 6u) | (c & 0x3Fu);
    }
    if ((n == 2 && cp < 0x800) || (n == 3 && cp < 0x10000) || cp > 0x10FFFF ||
        (cp >= 0xD800 && cp <= 0xDFFF)) {
      return false;
    }
    i += n + 1;
  }
  return true;
}

/// The records of a CSV text after RFC 4180, section 2 (CRLF and LF read, the last record may
/// end without a line end, a byte order mark is skipped).
std::vector<Record> parse_csv(std::string_view s) {
  if (s.starts_with("\xEF\xBB\xBF")) s.remove_prefix(3);
  std::vector<Record> records;
  std::size_t i = 0;
  std::size_t line = 1;
  while (i < s.size()) {
    Record record;
    bool end_of_record = false;
    while (!end_of_record) {
      Field f;
      f.line = line;
      f.column = record.size() + 1;
      if (i < s.size() && s[i] == '"') {
        ++i;  // rule 5 to 7: a quoted field, "" for a quote inside
        bool closed = false;
        while (i < s.size()) {
          if (s[i] == '"') {
            if (i + 1 < s.size() && s[i + 1] == '"') {
              f.text += '"';
              i += 2;
              continue;
            }
            ++i;
            closed = true;
            break;
          }
          if (s[i] == '\n') ++line;
          f.text += s[i++];
        }
        if (!closed) table_fail("table.parse", f.line, f.column, "a quote is not closed");
      } else {
        while (i < s.size() && s[i] != ',' && s[i] != '\n' && s[i] != '\r') {
          if (s[i] == '"') {
            table_fail("table.parse", f.line, f.column,
                       "a quote inside a field that does not start with one (RFC 4180, rule 5)");
          }
          f.text += s[i++];
        }
      }
      if (!valid_utf8(f.text)) table_fail("table.parse", f.line, f.column, "not UTF-8");
      record.push_back(std::move(f));
      if (i >= s.size()) {
        end_of_record = true;  // rule 2: the last record may end without a line end
      } else if (s[i] == ',') {
        ++i;
      } else if (s[i] == '\n' || (s[i] == '\r' && i + 1 < s.size() && s[i + 1] == '\n')) {
        i += s[i] == '\r' ? 2 : 1;
        ++line;
        end_of_record = true;
      } else {
        table_fail("table.parse", record.back().line, record.back().column,
                   s[i] == '\r' ? "a carriage return without a line feed"
                                : "text after the closing quote of a field");
      }
    }
    records.push_back(std::move(record));
  }
  // an empty line at the end, as text editors add one, is no record
  if (records.size() > 1 && records.back().size() == 1 && records.back().front().text.empty()) {
    records.pop_back();
  }
  return records;
}

/// One surface row, read.
struct Row {
  const Record* record = nullptr;
  ElementKind kind = ElementKind::Lens;
  model::Surface surface;
  std::optional<double> thickness;
  std::string material;
  bool stop = false;
};

class Reader {
 public:
  explicit Reader(const std::array<std::size_t, 22>& index) : index_(index) {}

  const Field& at(const Record& r, Column c) const { return r[index_[c]]; }
  const std::string& text(const Record& r, Column c) const { return at(r, c).text; }

  [[noreturn]] void bad(const Record& r, Column c, const std::string& message) const {
    const Field& f = at(r, c);
    table_fail("table.value", f.line, f.column, message);
  }

  double number(const Record& r, Column c) const {
    const std::optional<double> v = detail::parse_number(text(r, c));
    if (!v) bad(r, c, "'" + text(r, c) + "' is not a finite number");
    return *v;
  }

  std::optional<double> optional_number(const Record& r, Column c) const {
    if (text(r, c).empty()) return std::nullopt;
    return number(r, c);
  }

  void require_empty(const Record& r, Column c, const std::string& why) const {
    if (!text(r, c).empty()) bad(r, c, why);
  }

 private:
  const std::array<std::size_t, 22>& index_;
};

model::ShapeStack read_shape(const Reader& in, const Record& r) {
  model::ShapeStack shape;
  const std::string& name = in.text(r, kShape);
  const auto no_terms = [&] {
    for (std::size_t i = 0; i < 7; ++i) {
      in.require_empty(r, static_cast<Column>(kA4 + i), "only an even_asphere has a4 to a16");
    }
  };
  if (name == "plane") {
    in.require_empty(r, kRadius, "a plane has no radius");
    in.require_empty(r, kConic, "a plane has no conic");
    no_terms();
    shape.base = model::Plane{};
  } else if (name == "conic" || name == "even_asphere") {
    if (in.text(r, kRadius).empty()) in.bad(r, kRadius, name + " needs a radius");
    const double radius = in.number(r, kRadius);
    const double conic = in.optional_number(r, kConic).value_or(0.0);
    if (name == "conic") {
      no_terms();
      shape.base = model::Conic{radius, conic};
    } else {
      model::EvenAsphere a{radius, conic, {}};
      std::size_t last = 0;  // number of coefficients up to the last field that is set
      for (std::size_t i = 0; i < 7; ++i) {
        if (!in.text(r, static_cast<Column>(kA4 + i)).empty()) last = i + 1;
      }
      for (std::size_t i = 0; i < last; ++i) {
        a.coefficients.emplace_back(
            in.optional_number(r, static_cast<Column>(kA4 + i)).value_or(0.0));
      }
      shape.base = a;
    }
  } else {
    in.bad(r, kShape, "unknown shape '" + name + "' (plane, conic, even_asphere)");
  }
  return shape;
}

std::optional<model::Aperture> read_aperture(const Reader& in, const Record& r) {
  const std::string& kind = in.text(r, kAperture);
  if (kind.empty()) {
    in.require_empty(r, kAperture1, "aperture_1 without an aperture");
    in.require_empty(r, kAperture2, "aperture_2 without an aperture");
    return std::nullopt;
  }
  if (kind != "circular" && kind != "rectangular" && kind != "elliptical") {
    in.bad(r, kAperture, "unknown aperture '" + kind + "' (circular, rectangular, elliptical)");
  }
  if (in.text(r, kAperture1).empty()) in.bad(r, kAperture1, kind + " needs aperture_1");
  const double a1 = in.number(r, kAperture1);
  if (kind == "circular") {
    return model::CircularAperture{a1, in.optional_number(r, kAperture2).value_or(0.0)};
  }
  if (in.text(r, kAperture2).empty()) in.bad(r, kAperture2, kind + " needs aperture_2");
  const double a2 = in.number(r, kAperture2);
  if (kind == "rectangular") return model::RectangularAperture{a1, a2};
  return model::EllipticalAperture{a1, a2};
}

}  // namespace

TableExport write_surface_table(const model::System& system,
                                const material::MaterialLibrary& materials,
                                const ExportOptions& options) {
  detail::ExportPath path = detail::walk_export_path(system, materials, options);
  TableExport out;
  out.warnings = std::move(path.warnings);

  std::array<std::string, kColumns.size()> f;
  for (std::size_t i = 0; i < kColumns.size(); ++i) f[i] = kColumns[i];
  append_record(out.text, f);

  f.fill({});
  f[kRow] = "object";
  f[kThickness] = path.system.object.at_infinity ? std::string("Infinity")
                                                 : detail::format_number(path.object_distance);
  f[kMaterial] = path.object_medium.value_or(std::string());
  append_record(out.text, f);

  // consecutive rows of one element form it (point 7): an element left and entered again would
  // come back as two elements of the same name (decision G2 on #318)
  std::set<std::string> closed;
  for (std::size_t k = 0; k < path.rows.size(); ++k) {
    const ExportRow& row = path.rows[k];
    const std::string& id = row.surface.id.str();
    if (k > 0 && path.rows[k - 1].element != row.element) {
      closed.insert(path.rows[k - 1].element);
      if (closed.contains(row.element)) {
        export_fail("export.event_unsupported", id,
                    "the path leaves the element " + row.element +
                        " and enters it again; the table holds an element in one piece");
      }
    }
    f.fill({});
    f[kRow] = "surface";
    f[kSurface] = id;
    f[kElement] = row.element;
    f[kKind] = kind_name(row.kind);
    write_shape(row, f);
    if (k + 1 < path.rows.size()) f[kThickness] = detail::format_number(row.thickness);
    // An inner surface of a lens or plate names its segment, also one of the environment
    // medium (an air gap in a cemented group); the reader needs a material there.
    const bool inner = (row.kind == ElementKind::Lens || row.kind == ElementKind::Plate) &&
                       k + 1 < path.rows.size() && path.rows[k + 1].element == row.element;
    if (inner) {
      f[kMaterial] = row.medium_reference;
    } else if (!row.mirror) {
      f[kMaterial] = row.medium_after.value_or(std::string());
    }
    write_aperture(row, f);
    if (row.stop) f[kStop] = "true";
    if (const auto* lens = std::get_if<model::IdealLens>(&row.surface.interaction)) {
      f[kFocalLength] = detail::format_number(lens->focal_length.value);
      if (lens->object_distance) {
        f[kObjectDistance] = detail::format_number(lens->object_distance->value);
      }
    }
    append_record(out.text, f);
  }
  return out;
}

TableImport read_surface_table(std::string_view csv,
                               model::System base,
                               const material::MaterialLibrary& materials) {
  const std::vector<Record> records = parse_csv(csv);

  // Header (rule 3): the columns by name, in any order.
  if (records.empty()) table_fail("table.column", 0, 0, "the table has no header");
  const Record& header = records.front();
  if (header.size() == 1 && header.front().text.find(';') != std::string::npos) {
    table_fail("table.parse", 1, 1,
               "the separator is ';', not ','. A spreadsheet with German settings writes ';' "
               "and decimal commas; save the table as CSV with commas and decimal points");
  }
  std::array<std::size_t, 22> index{};
  std::array<bool, 22> found{};
  for (const Field& h : header) {
    const auto it = std::find(kColumns.begin(), kColumns.end(), h.text);
    if (it == kColumns.end()) {
      table_fail("table.column", h.line, h.column, "unknown column '" + h.text + "'");
    }
    const auto c = static_cast<std::size_t>(it - kColumns.begin());
    if (found[c]) table_fail("table.column", h.line, h.column, "column '" + h.text + "' twice");
    found[c] = true;
    index[c] = h.column - 1;
  }
  std::string missing;
  for (std::size_t c = 0; c < kColumns.size(); ++c) {
    if (!found[c]) missing += (missing.empty() ? "" : ", ") + std::string(kColumns[c]);
  }
  if (!missing.empty()) table_fail("table.column", 0, 0, "missing columns: " + missing);
  for (const Record& r : records) {
    if (r.size() != header.size()) {
      const std::size_t column = std::min(r.size(), header.size()) + 1;
      const std::size_t line = r.size() > header.size() ? r[header.size()].line : r.back().line;
      table_fail("table.parse", line, column,
                 "the record has " + std::to_string(r.size()) + " fields, the header " +
                     std::to_string(header.size()) + " (RFC 4180, rule 4)");
    }
  }
  const Reader in(index);

  // Object row (decision G4: only thickness and material).
  if (records.size() < 2 || in.text(records[1], kRow) != "object") {
    if (records.size() < 2) table_fail("table.value", 0, 0, "the table has no object row");
    in.bad(records[1], kRow, "the first record after the header is the object row");
  }
  const Record& object = records[1];
  for (std::size_t c = 0; c < kColumns.size(); ++c) {
    if (c != kRow && c != kThickness && c != kMaterial) {
      in.require_empty(object, static_cast<Column>(c),
                       "the object row holds only thickness and material");
    }
  }
  model::ObjectSpace object_space;
  const std::string& distance = in.text(object, kThickness);
  if (distance == "Infinity") {
    object_space.at_infinity = true;
  } else {
    if (distance.empty()) in.bad(object, kThickness, "the object needs a distance or Infinity");
    object_space.at_infinity = false;
    object_space.distance = in.number(object, kThickness);
  }
  if (!in.text(object, kMaterial).empty()) object_space.medium = in.text(object, kMaterial);

  // Surface rows.
  std::vector<Row> rows;
  std::set<std::string> ids;
  const Record* stop_field = nullptr;
  for (std::size_t k = 2; k < records.size(); ++k) {
    const Record& r = records[k];
    if (in.text(r, kRow) != "surface") {
      in.bad(r, kRow,
             in.text(r, kRow) == "object" ? "a second object row"
                                          : "the row is 'object' or 'surface'");
    }
    Row row;
    row.record = &r;
    const std::string& kind = in.text(r, kKind);
    const auto it =
        std::find_if(kKinds.begin(), kKinds.end(), [&](const auto& x) { return x.second == kind; });
    if (it == kKinds.end()) in.bad(r, kKind, "unknown kind '" + kind + "'");
    row.kind = it->first;
    if (!ids.insert(in.text(r, kSurface)).second) {
      const Field& f = in.at(r, kSurface);
      table_fail("table.duplicate", f.line, f.column, "surface id '" + f.text + "' used twice");
    }
    row.surface.id = model::SurfaceId(in.text(r, kSurface));
    row.surface.shape = read_shape(in, r);
    row.surface.aperture = read_aperture(in, r);
    if (k + 1 < records.size() || !in.text(r, kThickness).empty()) {
      if (in.text(r, kThickness).empty()) in.bad(r, kThickness, "a thickness is needed");
      row.thickness = in.number(r, kThickness);
    }
    row.material = in.text(r, kMaterial);
    const std::string& stop = in.text(r, kStop);
    if (!stop.empty()) {
      if (stop != "true") in.bad(r, kStop, "stop is 'true' or empty");
      if (stop_field != nullptr) in.bad(r, kStop, "a second stop");
      stop_field = &r;
      row.stop = true;
    }
    if (!in.text(r, kFocalLength).empty()) {
      if (row.kind != ElementKind::ThinElement) {
        in.bad(r, kFocalLength, "only a thin_element has a focal length (ideal lens)");
      }
      model::IdealLens lens{in.number(r, kFocalLength), std::nullopt};
      if (const std::optional<double> s = in.optional_number(r, kObjectDistance)) {
        lens.object_distance = model::Param(*s);
      }
      row.surface.interaction = lens;
    } else {
      in.require_empty(r, kObjectDistance, "object_distance belongs to a focal length");
    }
    if (row.kind == ElementKind::Mirror) row.surface.interaction = model::IdealMirror{};
    rows.push_back(std::move(row));
  }

  // Elements: consecutive rows of one name (point 7), in the canonical form of point 6.
  std::vector<model::Node> children;
  std::set<std::string> names;
  std::optional<model::SurfaceId> aperture_stop;
  std::size_t k = 0;
  while (k < rows.size()) {
    const std::string& name = in.text(*rows[k].record, kElement);
    if (!names.insert(name).second) {
      const Field& f = in.at(*rows[k].record, kElement);
      table_fail("table.duplicate", f.line, f.column,
                 "element '" + name + "' again after other rows");
    }
    std::size_t end = k;
    while (end < rows.size() && in.text(*rows[end].record, kElement) == name) {
      if (rows[end].kind != rows[k].kind) {
        in.bad(*rows[end].record, kKind, "the rows of element '" + name + "' differ in kind");
      }
      ++end;
    }
    model::Element e{name, rows[k].kind, {}, std::nullopt, {}};
    if (k > 0) {
      e.pose.reference = model::PoseReference::RelativeToPreceding;
      e.pose.position[2] = model::Param(rows[k - 1].thickness.value_or(0.0));
    }
    const bool glass = e.kind == ElementKind::Lens || e.kind == ElementKind::Plate;
    std::vector<std::string> media;
    for (std::size_t i = k; i < end; ++i) {
      Row& row = rows[i];
      const bool last = i + 1 == end;
      if (!glass || last) {
        in.require_empty(*row.record, kMaterial,
                         glass ? "the last surface of a lens or plate leads into the environment"
                               : "only lens and plate surfaces have a material after them");
      } else if (row.material.empty()) {
        in.bad(*row.record, kMaterial, "an inner surface of a lens or plate needs its material");
      } else {
        media.push_back(row.material);
      }
      if (i > k) {
        row.surface.pose.reference = model::PoseReference::RelativeToPreceding;
        row.surface.pose.position[2] = model::Param(rows[i - 1].thickness.value_or(0.0));
      }
      if (row.stop && e.kind != ElementKind::Stop) aperture_stop = row.surface.id;
      e.surfaces.push_back(std::move(row.surface));
    }
    if (!media.empty()) {
      if (std::all_of(media.begin(), media.end(),
                      [&](const std::string& m) { return m == media.front(); })) {
        e.material = media.front();
      } else {
        e.segment_materials = media;
      }
    }
    children.push_back(model::Node{std::move(e)});
    k = end;
  }

  // The base system (point 7): tree, object and stop replaced, the blocks that refer to the old
  // tree dropped with a warning, bound Params resolved for configuration 0.
  TableImport out;
  try {
    if (!base.parameters.empty()) base = model::resolve_parameters(base, 0);
  } catch (const std::exception& e) {
    table_fail("table.invalid_result", 0, 0,
               std::string("the parameter table of the base system is invalid: ") + e.what());
  }
  const auto drop = [&](const char* where, const char* what) {
    out.warnings.push_back(detail::export_warning("table.base_dropped", where,
                                                  std::string(what) +
                                                      " of the base system dropped: it refers "
                                                      "to the tree the table replaces"));
  };
  const std::string path_name = base.paths.empty() ? std::string("main") : base.paths[0].name;
  if (base.paths.size() > 1 || (base.paths.size() == 1 && !base.paths[0].automatic)) {
    drop("/paths", "the paths");
  }
  if (!base.optimization.empty()) drop("/optimization", "the optimization");
  if (!base.tolerances.empty()) drop("/tolerances", "the tolerances");
  if (!base.parameters.empty()) drop("/parameters", "the parameter table");
  if (!base.configurations.empty()) drop("/configurations", "the configurations");
  base.optimization = {};
  base.tolerances = {};
  base.parameters.clear();
  base.configurations.clear();
  base.paths = {model::Path{path_name, true, {}}};
  base.object = object_space;
  base.aperture.stop = aperture_stop;
  base.root = model::Assembly{base.root.name, {}, std::move(children)};

  try {
    const compile::CompiledSystem cs = compile::compile(base, materials);
    for (const model::Diagnostic& d : cs.diagnostics()) out.warnings.push_back(d);
  } catch (const compile::CompileError& e) {
    std::string message;
    for (const model::Diagnostic& d : e.diagnostics()) {
      message += (message.empty() ? "" : "; ") + model::to_string(d);
    }
    table_fail("table.invalid_result", 0, 0, "the system of the table is invalid: " + message);
  } catch (const std::exception& e) {
    table_fail("table.invalid_result", 0, 0,
               std::string("the system of the table is invalid: ") + e.what());
  }
  out.system = std::move(base);
  return out;
}

}  // namespace rtt::io
