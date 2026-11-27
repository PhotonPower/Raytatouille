#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include "export_path.hpp"
#include "numbers.hpp"
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

std::string num(double v) {
  return detail::format_number(v);
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

/// A surface with one shape term (ADR 0045, addendum #430): zernike_sag as SZERNSAG,
/// zernike_fringe_sag as FZERNSAG, both on a plane, conic or even asphere up to A16 (UG "Zernike
/// Standard Sag", "Zernike Fringe Sag": the base of "Even Asphere", PARM 2 ... 8 = A4 ... A16;
/// an odd asphere with c = 0 and only these even terms as well, addendum #436);
/// xy_polynomial as XPOLYNOM on a plane or conic (UG "Extended Polynomial"). Extrapolate = 1
/// (PARM 0), since the model evaluates the terms everywhere (ADR 0051, point 1); the terms as
/// extra data, XDAT 1 = number of terms, 2 = normalization radius, 2 + j = A_j (XDAT n =
/// parameter 12 + n, ADR 0044, addendum #430), lengths in mm (UNIT MM). The ZMX frame must be
/// that of the surface: a turned or mirrored frame would change the signs of the terms
/// (mirror_shape rejects the mirrored one).
void write_term_shape(Writer& w, const ExportRow& row) {
  const std::string& id = row.surface.id.str();
  const model::ShapeStack& shape = row.surface.shape;
  if (shape.terms.size() != 1) {
    export_fail("export.surface_unsupported", id, "more than one shape term has no mapping");
  }
  if (row.local_xy != std::array<double, 4>{1.0, 0.0, 0.0, 1.0}) {
    export_fail("export.surface_unsupported", id,
                "a shape term on a surface turned in its ZMX frame has no mapping");
  }
  const model::ShapeTerm& term = shape.terms.front();
  const auto* noll = std::get_if<model::ZernikeSag>(&term);
  const auto* fringe = std::get_if<model::ZernikeFringeSag>(&term);
  const auto* xy = std::get_if<model::XyPolynomial>(&term);
  if (noll == nullptr && fringe == nullptr && xy == nullptr) {
    export_fail("export.surface_unsupported", id, "the shape term has no ZMX surface type");
  }
  double c = 0.0;
  double k = 0.0;
  std::array<double, 7> a{};  // A4 ... A16
  if (const auto* conic = std::get_if<model::Conic>(&shape.base)) {
    c = curvature(conic->radius.value, id);
    k = conic->conic.value;
  } else if (const auto* e = std::get_if<model::EvenAsphere>(&shape.base);
             e != nullptr && xy == nullptr) {
    c = curvature(e->radius.value, id);
    k = e->conic.value;
    for (std::size_t i = 0; i < e->coefficients.size(); ++i) {
      if (i < a.size()) {
        a[i] = e->coefficients[i].value;
      } else if (e->coefficients[i].value != 0.0) {
        export_fail("export.surface_unsupported", id,
                    "the base of a Zernike surface holds the terms up to A16 only");
      }
    }
  } else if (const auto* o = std::get_if<model::OddAsphere>(&shape.base);
             o != nullptr && xy == nullptr) {
    // the base the import makes of CURV 0 (ADR 0044, point 6; ADR 0045, addendum #436): c = 0
    // and only the even terms A4 ... A16, as EVENASPH with c = 0 (addendum #368); k has no
    // effect at c = 0 and is written as 0
    if (o->curvature.value != 0.0) {
      export_fail("export.surface_unsupported", id,
                  "the base of a Zernike surface is an odd asphere with c = 0 only");
    }
    for (std::size_t i = 0; i < o->coefficients.size(); ++i) {
      const std::size_t n = i + 1;  // A_n
      const double v = o->coefficients[i].value;
      if (n >= 4 && n <= 16 && n % 2 == 0) {
        a[n / 2 - 2] = v;
      } else if (v != 0.0) {
        export_fail("export.surface_unsupported", id,
                    "the base of a Zernike surface holds the even terms A4 to A16 only");
      }
    }
  } else if (!std::holds_alternative<model::Plane>(shape.base)) {
    export_fail("export.surface_unsupported", id,
                xy != nullptr ? "XPOLYNOM has a plane or conic base only"
                              : "the base of a Zernike surface is a plane, conic, even asphere or "
                                "odd asphere with c = 0");
  }
  const double radius = noll != nullptr     ? noll->normalization_radius.value
                        : fringe != nullptr ? fringe->normalization_radius.value
                                            : xy->normalization_radius.value;
  const std::vector<model::Param>& coefficients = noll != nullptr     ? noll->coefficients
                                                  : fringe != nullptr ? fringe->coefficients
                                                                      : xy->coefficients;
  w.sub(std::string("TYPE ") + (noll != nullptr     ? "SZERNSAG"
                                : fringe != nullptr ? "FZERNSAG"
                                                    : "XPOLYNOM"));
  w.sub("CURV " + num(c));
  w.sub("CONI " + num(k));
  if (xy == nullptr) {
    w.sub("PARM 0 1");
    w.sub("PARM 1 0");
    for (std::size_t i = 0; i < a.size(); ++i)
      w.sub("PARM " + std::to_string(i + 2) + " " + num(a[i]));
  }
  // as OpticStudio writes extra data (ray-optics test file HoO-V2C18Ex46.zmx)
  const auto xdat = [&](std::size_t n, double v) {
    w.sub("XDAT " + std::to_string(n) + " " + num(v) + " 0 0 1 0 0 \"\"");
  };
  xdat(1, static_cast<double>(coefficients.size()));
  xdat(2, radius);
  for (std::size_t j = 0; j < coefficients.size(); ++j) xdat(j + 3, coefficients[j].value);
}

/// TYPE, CURV, CONI and PARM of a surface in the global frame (ADR 0045, point 4).
void write_shape(Writer& w, const ExportRow& row) {
  const std::string& id = row.surface.id.str();
  const model::ShapeStack& shape = row.surface.shape;
  if (!shape.terms.empty()) {
    write_term_shape(w, row);
    return;
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
    // EVENASPH first (ADR 0045, addendum #368): c = 0 and only even terms A4 ... A16 stay
    // EVENASPH with c = 0, which the import reads back as the same odd asphere; the others with
    // terms up to A8 are ODDASPHE. UG "Odd Asphere": z = c r^2 / (1 + sqrt(1 - (1 + k) c^2 r^2))
    // + beta_1 r + ... + beta_8 r^8, parameters 1 to 8 = beta_1 ... beta_8, in lens units; the
    // file is in mm (UNIT MM), so beta_n = A_n (ADR 0043, point 1; ADR 0044, addendum #360).
    std::array<double, 7> even_terms{};
    std::array<double, 8> beta{};
    bool as_even = o->curvature.value == 0.0;
    bool as_odd = true;
    for (std::size_t i = 0; i < o->coefficients.size(); ++i) {
      const std::size_t n = i + 1;  // A_n
      const double v = o->coefficients[i].value;
      if (n >= 4 && n <= 16 && n % 2 == 0) {
        even_terms[n / 2 - 2] = v;
      } else if (v != 0.0) {
        as_even = false;
      }
      if (i < beta.size()) {
        beta[i] = v;
      } else if (v != 0.0) {
        as_odd = false;
      }
    }
    if (as_even) {
      even(0.0, o->conic.value, even_terms);
    } else if (as_odd) {
      w.sub("TYPE ODDASPHE");
      w.sub("CURV " + num(o->curvature.value));
      w.sub("CONI " + num(o->conic.value));
      for (std::size_t i = 0; i < beta.size(); ++i)
        w.sub("PARM " + std::to_string(i + 1) + " " + num(beta[i]));
    } else {
      export_fail("export.surface_unsupported", id,
                  "an odd asphere is exported as EVENASPH (c = 0 and only the even terms A4 to "
                  "A16) or as ODDASPHE (the terms up to A8)");
    }
  } else {
    export_fail("export.surface_unsupported", id, "the shape has no ZMX surface type");
  }
}

template <class>
inline constexpr bool kUnhandledExportAperture = false;

/// DIAM and the aperture lines (ADR 0045, point 4; ADR 0044, point 8).
void write_aperture(Writer& w, const ExportRow& row) {
  const std::optional<model::Aperture>& a = row.surface.aperture;
  if (row.stop && a) {
    // the stop: fixed semi-diameter (rule 2), an annulus also as CLAP (rule 1); validate gives
    // every stop an aperture, the walk only lets a circle through
    const auto& c = std::get<model::CircularAperture>(*a);
    w.sub("DIAM " + num(c.radius.value) + " 1");
    if (c.inner_radius.value > 0.0)
      w.sub("CLAP " + num(c.inner_radius.value) + " " + num(c.radius.value));
    return;
  }
  w.sub("DIAM 0 0");
  // An obscuration (ADR 0050, point 3): on a surface without aperture, unrotated (checked in
  // export_zmx); OBSC rmin rmax with rmin = 0, SQOB and ELOB with the half widths, OBDC with the
  // decenter of the centre (UG "Aperture Type Codes and Parameters"), all in the x, y of the ZMX
  // frame. ExportRow::local_xy is diagonal (+-1): it flips the signs of the decenter, never swaps
  // x and y, so the half widths stay.
  if (!row.surface.obscurations.empty()) {
    const std::array<double, 4>& m = row.local_xy;
    std::visit(
        [&](const auto& o) {
          using O = std::decay_t<decltype(o)>;
          if constexpr (std::is_same_v<O, model::CircularObscuration>) {
            w.sub("OBSC 0 " + num(o.radius.value));
          } else if constexpr (std::is_same_v<O, model::RectangularObscuration>) {
            w.sub("SQOB " + num(o.half_width_x.value) + " " + num(o.half_width_y.value));
          } else {
            static_assert(std::is_same_v<O, model::EllipticalObscuration>,
                          "write_aperture: unhandled obscuration");
            w.sub("ELOB " + num(o.semi_axis_x.value) + " " + num(o.semi_axis_y.value));
          }
          const double dx = o.decenter_x.value;
          const double dy = o.decenter_y.value;
          if (dx != 0.0 || dy != 0.0) {
            // + 0.0 turns a -0 from a sign flip of 0 into 0
            w.sub("OBDC " + num(m[0] * dx + 0.0) + " " + num(m[3] * dy + 0.0));
          }
        },
        row.surface.obscurations.front());
  }
  if (!a) return;
  std::visit(
      [&w](const auto& ap) {
        using A = std::decay_t<decltype(ap)>;
        if constexpr (std::is_same_v<A, model::CircularAperture>) {
          w.sub("CLAP " + num(ap.inner_radius.value) + " " + num(ap.radius.value));
        } else if constexpr (std::is_same_v<A, model::RectangularAperture>) {
          w.sub("SQAP " + num(ap.half_width_x.value) + " " + num(ap.half_width_y.value));
        } else if constexpr (std::is_same_v<A, model::EllipticalAperture>) {
          w.sub("ELAP " + num(ap.semi_axis_x.value) + " " + num(ap.semi_axis_y.value));
        } else if constexpr (std::is_same_v<A, model::RegularPolygonAperture> ||
                             std::is_same_v<A, model::PolygonAperture>) {
          // Unreachable: walk_export_path refuses them (export.aperture_unsupported, ADR 0050).
        } else {
          static_assert(kUnhandledExportAperture<A>, "write_aperture: unhandled aperture");
        }
      },
      *a);
}

}  // namespace

ZmxExport export_zmx(const model::System& system,
                     const material::MaterialLibrary& materials,
                     const ExportOptions& options) {
  detail::ExportPath path = detail::walk_export_path(system, materials, options, true);
  // OpticStudio has one aperture type per surface and no rotation of an aperture (UG "Aperture
  // Type Codes and Parameters", SURP codes 20 to 24): one obscuration on a surface without
  // aperture, not rotated, or export.aperture_unsupported before anything is written
  // (ADR 0050, point 3).
  for (const detail::ExportRow& row : path.rows) {
    const model::Surface& surface = row.surface;
    if (surface.obscurations.empty()) continue;
    const std::string id = surface.id.str();
    if (surface.aperture) {
      export_fail("export.aperture_unsupported", id,
                  "a surface with an aperture and an obscuration has no mapping (one aperture "
                  "type per ZMX surface)");
    }
    if (surface.obscurations.size() > 1) {
      export_fail("export.aperture_unsupported", id,
                  "more than one obscuration on a surface has no mapping");
    }
    const bool rotated = std::visit(
        [](const auto& o) {
          if constexpr (std::is_same_v<std::decay_t<decltype(o)>, model::CircularObscuration>) {
            return false;
          } else {
            return o.rotation_deg != 0.0;
          }
        },
        surface.obscurations.front());
    if (rotated) {
      export_fail("export.aperture_unsupported", id, "a rotated obscuration has no mapping");
    }
  }
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

  // Coordinate breaks of a folded path (ADR 0045, addendum #366, point 3): flag 0, the medium
  // of the surface before them (ADR 0044, addendum PR C).
  std::size_t number = 0;
  const auto breaks = [&](const std::vector<detail::CoordinateBreak>& list,
                          const std::optional<Glass>& glass) {
    for (const detail::CoordinateBreak& b : list) {
      w.line("SURF " + std::to_string(++number));
      w.sub("TYPE COORDBRK");
      w.sub("PARM 1 " + num(b.dx));
      w.sub("PARM 2 " + num(b.dy));
      w.sub("PARM 3 " + num(b.tilt_deg[0]));
      w.sub("PARM 4 " + num(b.tilt_deg[1]));
      w.sub("PARM 5 " + num(b.tilt_deg[2]));
      w.sub("PARM 6 0");
      w.sub("DISZ " + num(b.thickness));
      if (glass) w.sub("GLAS " + glass->name + " 0");
    }
  };

  w.line("SURF 0");
  // A tilted object plane (ADR 0055): TILTSURF with PARM 1 = sx, PARM 2 = sy (UG "Tilted",
  // z = x tan(theta_x) + y tan(theta_y)), as the import reads it.
  const double slope_x = s.object.slope[0].value;
  const double slope_y = s.object.slope[1].value;
  const bool tilted_object = !s.object.at_infinity && (slope_x != 0.0 || slope_y != 0.0);
  w.sub(tilted_object ? "TYPE TILTSURF" : "TYPE STANDARD");
  w.sub("CURV 0");
  if (tilted_object) {
    w.sub("PARM 1 " + num(slope_x));
    w.sub("PARM 2 " + num(slope_y));
  }
  w.sub(s.object.at_infinity ? std::string("DISZ INFINITY") : "DISZ " + num(path.object_distance));
  if (object_glass) w.sub("GLAS " + object_glass->name + " 0");
  w.sub("DIAM 0 0");
  breaks(path.object_breaks, object_glass);
  for (std::size_t k = 0; k < path.rows.size(); ++k) {
    const ExportRow& row = path.rows[k];
    if (row.kind == model::ElementKind::ThinElement && !row.surface.aperture &&
        !std::holds_alternative<model::IdealLens>(row.surface.interaction)) {
      // an air-air surface without aperture: the import drops it (ADR 0044, point 5)
      export_fail("export.surface_unsupported", row.surface.id.str(),
                  "a thin element without aperture and ideal lens has no ZMX surface the import "
                  "keeps");
    }
    w.line("SURF " + std::to_string(++number));
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
    breaks(row.breaks, row.mirror ? std::nullopt : glasses[k]);
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
