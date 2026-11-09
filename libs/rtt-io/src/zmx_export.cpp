#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "export_path.hpp"
#include "rtt/io/zmx.hpp"
#include "rtt/material/air.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/element.hpp"
#include "rtt/model/surface.hpp"
#include "rtt/model/system.hpp"
#include "zmx_names.hpp"

namespace rtt::io {

namespace {

using detail::export_fail;
using detail::export_warning;
using detail::ExportRow;

/// Most fields and wavelengths of a file: a reading of the corpus, the test files and
/// ray-optics (zmxread.py); no page of the User Guide names them (ADR 0045, point 4).
constexpr std::size_t kMaxFields = 12;
constexpr std::size_t kMaxWavelengths = 24;

/// The shortest text that reads back to the same double ([charconv.to.chars], paragraphs 2, 7
/// and 8; ADR 0045, F3).
std::string num(double v) {
  std::array<char, 32> buf{};
  const auto [end, ec] = std::to_chars(buf.data(), buf.data() + buf.size(), v);
  if (ec != std::errc()) throw std::logic_error("to_chars failed");
  return {buf.data(), end};
}

/// The lines of the file; each ends with CRLF as in the files of OpticStudio.
class Writer {
 public:
  void line(const std::string& s) { text_ += s + "\r\n"; }
  /// A line of a surface block, indented: the import gives indented lines to the last SURF.
  void sub(const std::string& s) { line("  " + s); }
  [[nodiscard]] std::string take() { return std::move(text_); }

 private:
  std::string text_;
};

/// The glass of a material reference "CATALOG:NAME", checked so that the import finds it as a
/// catalogue glass (ADR 0045, point 4; ADR 0044, point 7).
struct Glass {
  std::string catalog;
  std::string name;
};

Glass glass_of(const std::string& reference,
               const std::string& surface,
               const std::vector<std::string>& catalogs) {
  const std::size_t colon = reference.find(':');
  const std::string catalog =
      colon == std::string::npos ? std::string() : reference.substr(0, colon);
  const std::string name = colon == std::string::npos ? reference : reference.substr(colon + 1);
  if (catalog.empty() || std::find(catalogs.begin(), catalogs.end(), catalog) == catalogs.end()) {
    export_fail("export.material_unsupported", surface,
                reference +
                    " is no glass of a loaded catalogue (CONST: and VACUUM have no "
                    "mapping, a model glass would be another glass)");
  }
  // The import splits a line at white space: a blank in a GCAT or GLAS word splits the name.
  const auto blank = [](const std::string& s) {
    return s.empty() ||
           std::any_of(s.begin(), s.end(), [](unsigned char c) { return c <= 0x20 || c == 0x7F; });
  };
  if (blank(catalog)) {
    export_fail("export.material_unsupported", surface,
                "the catalogue name '" + catalog + "' holds a blank, the GCAT line would split it");
  }
  const std::string u = detail::upper(name);
  if (blank(name) || u == "AIR" || u == "MIRROR" || name == "___BLANK" ||
      detail::mil_number(name)) {
    export_fail("export.material_unsupported", surface,
                "the import would not read the glass name '" + name + "' as a catalogue glass");
  }
  return {catalog, name};
}

/// GCAT: the catalogues so ordered that the import, which takes the first catalogue holding a
/// name, finds each glass in its own catalogue. A glass of C whose name also resolves in D
/// needs C before D; Kahn's algorithm, ties in name order. No order: export.glass_ambiguous.
std::vector<std::string> gcat_order(const std::map<std::string, std::set<std::string>>& used,
                                    const std::map<std::string, std::string>& first_surface,
                                    const material::MaterialLibrary& lib) {
  const auto resolves = [&](const std::string& reference) {
    try {
      (void)lib.resolve(reference);
      return true;
    } catch (const material::UnknownMaterial&) {
      return false;
    }
  };
  // The import matches catalogue names without case: a second loaded catalogue of the same name
  // in another case could be the one it finds.
  for (const auto& [catalog, glasses] : used) {
    for (const std::string& other : lib.catalogs()) {
      if (other != catalog && detail::upper(other) == detail::upper(catalog)) {
        std::string message = "catalogues ";
        message.append(catalog).append(" and ").append(other).append(" differ only in case");
        export_fail("export.glass_ambiguous", first_surface.at(catalog), message);
      }
    }
  }
  std::map<std::string, std::set<std::string>> after;  // C -> catalogues that must follow C
  std::map<std::string, std::size_t> before;           // number of catalogues that must precede
  for (const auto& [catalog, glasses] : used) before.emplace(catalog, 0);
  for (const auto& [catalog, glasses] : used) {
    for (const std::string& name : glasses) {
      for (const auto& [other, unused] : used) {
        std::string reference = other;
        reference.append(":").append(name);
        if (other != catalog && resolves(reference) && after[catalog].insert(other).second) {
          ++before[other];
        }
      }
    }
  }
  std::vector<std::string> order;
  std::set<std::string> ready;
  for (const auto& [catalog, n] : before) {
    if (n == 0) ready.insert(catalog);
  }
  while (!ready.empty()) {
    const std::string c = *ready.begin();
    ready.erase(ready.begin());
    order.push_back(c);
    for (const std::string& d : after[c]) {
      if (--before[d] == 0) ready.insert(d);
    }
  }
  if (order.size() != used.size()) {
    std::string left;
    for (const auto& [catalog, n] : before) {
      if (n > 0) left += (left.empty() ? "" : ", ") + catalog;
    }
    export_fail("export.glass_ambiguous", "",
                "no order of the catalogues " + left +
                    " finds every glass in its own catalogue (a name is in several of them)");
  }
  return order;
}

/// 1/R in 1/mm; a radius so small that 1/R overflows has no CURV the import reads.
double curvature(double radius, const std::string& surface) {
  const double c = 1.0 / radius;
  if (!std::isfinite(c)) {
    export_fail("export.surface_unsupported", surface, "the curvature 1/R is not finite");
  }
  return c;
}

/// TYPE, CURV, CONI and PARM of a surface in the global frame (ADR 0045, point 4).
void write_shape(Writer& w, const ExportRow& row) {
  const std::string& id = row.surface.id.str();
  const model::ShapeStack& shape = row.surface.shape;
  if (!shape.terms.empty()) {
    export_fail("export.surface_unsupported", id, "shape terms have no mapping");
  }
  if (const auto* lens = std::get_if<model::IdealLens>(&row.surface.interaction)) {
    // UG "Paraxial": PARM 1 focal length, PARM 2 OPD mode; 2 is the infinite conjugate
    w.sub("TYPE PARAXIAL");
    w.sub("CURV 0");
    w.sub("PARM 1 " + num(lens->focal_length.value));
    w.sub(std::string("PARM 2 ") + (lens->object_distance ? "0" : "2"));
    return;
  }
  const auto even = [&](double c, double k, const std::array<double, 7>& a) {
    w.sub("TYPE EVENASPH");
    w.sub("CURV " + num(c));
    w.sub("CONI " + num(k));
    w.sub("PARM 1 0");
    for (std::size_t i = 0; i < a.size(); ++i)
      w.sub("PARM " + std::to_string(i + 2) + " " + num(a[i]));
  };
  if (std::holds_alternative<model::Plane>(shape.base)) {
    w.sub("TYPE STANDARD");
    w.sub("CURV 0");
  } else if (const auto* c = std::get_if<model::Conic>(&shape.base)) {
    w.sub("TYPE STANDARD");
    w.sub("CURV " + num(curvature(c->radius.value, id)));
    w.sub("CONI " + num(c->conic.value));
  } else if (const auto* e = std::get_if<model::EvenAsphere>(&shape.base)) {
    // PARM 2 ... 8 = A4 ... A16 (ADR 0044, point 6)
    std::array<double, 7> a{};
    for (std::size_t i = 0; i < e->coefficients.size(); ++i) {
      if (i < a.size()) {
        a[i] = e->coefficients[i].value;
      } else if (e->coefficients[i].value != 0.0) {
        export_fail("export.surface_unsupported", id, "EVENASPH holds the terms up to A16 only");
      }
    }
    even(curvature(e->radius.value, id), e->conic.value, a);
  } else if (const auto* o = std::get_if<model::OddAsphere>(&shape.base)) {
    // c = 0 and only even terms A4 ... A16: EVENASPH with c = 0, as the import reads it back
    std::array<double, 7> a{};
    bool ok = o->curvature.value == 0.0;
    for (std::size_t i = 0; i < o->coefficients.size(); ++i) {
      const std::size_t n = i + 1;  // A_n
      const double v = o->coefficients[i].value;
      if (n >= 4 && n <= 16 && n % 2 == 0) {
        a[n / 2 - 2] = v;
      } else if (v != 0.0) {
        ok = false;
      }
    }
    if (!ok) {
      export_fail("export.surface_unsupported", id,
                  "an odd asphere is exported only with c = 0 and the even terms A4 to A16");
    }
    even(0.0, o->conic.value, a);
  } else {
    export_fail("export.surface_unsupported", id, "the shape has no ZMX surface type");
  }
}

/// DIAM and the aperture lines (ADR 0045, point 4; ADR 0044, point 8).
void write_aperture(Writer& w, const ExportRow& row) {
  const std::optional<model::Aperture>& a = row.surface.aperture;
  if (row.stop && a) {
    // the stop: fixed semi-diameter (rule 2), an annulus also as CLAP (rule 1); validate gives
    // every stop an aperture, the walk only lets a circle through
    const auto& c = std::get<model::CircularAperture>(*a);
    w.sub("DIAM " + num(c.radius) + " 1");
    if (c.inner_radius > 0.0) w.sub("CLAP " + num(c.inner_radius) + " " + num(c.radius));
    return;
  }
  w.sub("DIAM 0 0");
  if (!a) return;
  if (const auto* c = std::get_if<model::CircularAperture>(&*a)) {
    w.sub("CLAP " + num(c->inner_radius) + " " + num(c->radius));
  } else if (const auto* r = std::get_if<model::RectangularAperture>(&*a)) {
    w.sub("SQAP " + num(r->half_width_x) + " " + num(r->half_width_y));
  } else if (const auto* e = std::get_if<model::EllipticalAperture>(&*a)) {
    w.sub("ELAP " + num(e->semi_axis_x) + " " + num(e->semi_axis_y));
  }
}

}  // namespace

ZmxExport export_zmx(const model::System& system,
                     const material::MaterialLibrary& materials,
                     const ExportOptions& options) {
  detail::ExportPath path = detail::walk_export_path(system, materials, options);
  const model::System& s = path.system;
  ZmxExport out;
  out.warnings = std::move(path.warnings);

  // System data (ADR 0045, point 4).
  if (s.environment.medium != "AIR" || s.environment.temperature_c != 20.0 ||
      s.environment.pressure_atm != 1.0) {
    export_fail("export.environment_unsupported", "",
                "the indices of a ZMX file are relative to AIR at 20 degC and 1 atm");
  }
  if (s.fields.points.size() > kMaxFields || s.wavelengths.size() > kMaxWavelengths) {
    export_fail("export.limit_exceeded", "",
                "a ZMX file holds at most 12 fields and 24 wavelengths");
  }
  std::size_t ref = 0;
  for (std::size_t i = 0; i < s.wavelengths.size(); ++i) {
    if (s.wavelengths[i].reference) ref = i;
  }
  const double n_air = material::ciddor_air_index(s.wavelengths[ref].um, 20.0, 1.0);
  std::string aperture;
  using AT = model::SystemApertureType;
  switch (s.aperture.type) {
    case AT::EntrancePupilDiameter:
      aperture = "ENPD " + num(s.aperture.value.value);
      break;
    case AT::ImageSpaceFNumber:
      aperture = "FNUM " + num(s.aperture.value.value * n_air) + " 0";
      break;
    case AT::ObjectSpaceNA:
      if (s.aperture.reading != model::ApertureReading::Real) {
        export_fail("export.aperture_unsupported", "",
                    "object_na read paraxially has no mapping (OBNA is the real NA)");
      }
      aperture = "OBNA " + num(s.aperture.value.value / n_air) + " 0";
      break;
    case AT::StopSize:
      aperture = "FLOA";
      break;
  }

  // Glasses and their catalogues (ADR 0045, point 4; F4).
  const std::vector<std::string> catalogs = materials.catalogs();
  std::map<std::string, std::set<std::string>> used;
  std::map<std::string, std::string> first_surface;
  std::optional<Glass> object_glass;
  const auto use = [&](const std::string& reference, const std::string& surface) {
    Glass g = glass_of(reference, surface, catalogs);
    used[g.catalog].insert(g.name);
    first_surface.emplace(g.catalog, surface);
    return g;
  };
  if (path.object_medium) object_glass = use(*path.object_medium, "");
  std::vector<std::optional<Glass>> glasses;
  glasses.reserve(path.rows.size());
  for (const ExportRow& row : path.rows) {
    glasses.push_back(row.medium_after
                          ? std::optional<Glass>(use(*row.medium_after, row.surface.id.str()))
                          : std::nullopt);
  }
  const std::vector<std::string> gcat = gcat_order(used, first_surface, materials);

  // NAME: a control character would break the line.
  std::string name = s.name;
  if (std::any_of(name.begin(), name.end(), [](unsigned char c) { return c < 0x20; })) {
    std::replace_if(name.begin(), name.end(), [](unsigned char c) { return c < 0x20; }, ' ');
    out.warnings.push_back(export_warning("export.name_changed", "/name",
                                          "control characters of the name replaced by spaces"));
  }

  Writer w;
  w.line("VERS 60801 67 21410");
  w.line("MODE SEQ");
  w.line(name.empty() ? std::string("NAME") : "NAME " + name);
  w.line("UNIT MM");
  w.line(aperture);
  w.line("ENVD 20 1 0");
  if (!gcat.empty()) {
    std::string g = "GCAT";
    for (const std::string& c : gcat) g += " " + c;
    w.line(g);
  }
  static constexpr std::array<int, 3> kFieldTypes{0, 1, 2};  // angle, object height, image height
  w.line("FTYP " + std::to_string(kFieldTypes[static_cast<std::size_t>(s.fields.type)]) + " 0 " +
         std::to_string(s.fields.points.size()) + " " + std::to_string(s.wavelengths.size()) +
         " 0 0 0");
  std::string xs = "XFLN";
  std::string ys = "YFLN";
  std::string ws = "FWGN";
  for (const model::Field& f : s.fields.points) {
    xs += " " + num(f.x);
    ys += " " + num(f.y);
    ws += " " + num(f.weight);
  }
  w.line(xs);
  w.line(ys);
  w.line(ws);
  for (std::size_t i = 0; i < s.wavelengths.size(); ++i) {
    // lambda_air = lambda_vac / n_air(lambda_vac): the equation the import solves as a fixed
    // point (compile::line_vacuum_wavelength_um), solved for lambda_air (ADR 0045, point 4)
    const double um = s.wavelengths[i].um;
    w.line("WAVM " + std::to_string(i + 1) + " " +
           num(um / material::ciddor_air_index(um, 20.0, 1.0)) + " " +
           num(s.wavelengths[i].weight));
  }
  w.line("PWAV " + std::to_string(ref + 1));

  w.line("SURF 0");
  w.sub("TYPE STANDARD");
  w.sub("CURV 0");
  w.sub(s.object.at_infinity ? std::string("DISZ INFINITY") : "DISZ " + num(path.object_distance));
  if (object_glass) w.sub("GLAS " + object_glass->name + " 0");
  w.sub("DIAM 0 0");
  for (std::size_t k = 0; k < path.rows.size(); ++k) {
    const ExportRow& row = path.rows[k];
    if (row.kind == model::ElementKind::ThinElement && !row.surface.aperture &&
        !std::holds_alternative<model::IdealLens>(row.surface.interaction)) {
      // an air-air surface without aperture: the import drops it (ADR 0044, point 5)
      export_fail("export.surface_unsupported", row.surface.id.str(),
                  "a thin element without aperture and ideal lens has no ZMX surface the import "
                  "keeps");
    }
    w.line("SURF " + std::to_string(k + 1));
    if (row.stop) w.sub("STOP");
    write_shape(w, row);
    w.sub("DISZ " + num(row.thickness));
    if (row.mirror) {
      w.sub("GLAS MIRROR 0");
    } else if (const std::optional<Glass>& g = glasses[k]; g.has_value()) {
      w.sub("GLAS " + g->name + " 0");
    }
    write_aperture(w, row);
    if (const auto* lens = std::get_if<model::IdealLens>(&row.surface.interaction);
        lens != nullptr && lens->object_distance) {
      out.warnings.push_back(export_warning(
          "export.paraxial_opd", row.location,
          "ideal lens " + row.surface.id.str() +
              " written with OPD mode 0: the rays are the same, the design conjugate is not kept"));
    }
  }
  w.line("MNUM 1");
  out.text = w.take();
  return out;
}

std::vector<model::Diagnostic> save_zmx(const std::filesystem::path& file,
                                        const model::System& system,
                                        const material::MaterialLibrary& materials,
                                        const ExportOptions& options) {
  ZmxExport out = export_zmx(system, materials, options);
  std::ofstream f(file, std::ios::binary | std::ios::trunc);
  if (!f) throw std::runtime_error("cannot write " + file.string());
  f << out.text;
  if (!f) throw std::runtime_error("cannot write " + file.string());
  return std::move(out.warnings);
}

}  // namespace rtt::io
