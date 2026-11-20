// Writer of the result format in C++ (ADR 0023; ADR 0052, point 1): byte for byte the text of
// raytatouille.results.to_json with indent=None. The keys, their order, the Python kinds, dtypes
// and shapes follow the Python binding (libs/rtt-py/src/bind_*.cpp) and the attribute lists of
// results.py (_FIELDS); each function names its sources.

#include "rtt/io/result_json.hpp"

#include <array>
#include <charconv>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <type_traits>
#include <utility>
#include <vector>

#include "rtt/io/result_names.hpp"

namespace rtt::io {

std::string python_float_repr(double x) {
  // The shortest digits of std::to_chars in scientific notation (the general form would choose
  // the shorter of fixed and scientific, not Python's rule), then Python's layout
  // (float_repr_style "short": repr uses scientific notation for decpt <= -4 or decpt > 16).
  if (!std::isfinite(x)) {
    throw std::invalid_argument(
        "python_float_repr: only finite values (NaN and infinities are "
        "written as strings, ADR 0023)");
  }
  std::array<char, 64> buffer{};
  const std::to_chars_result r =
      std::to_chars(buffer.data(), buffer.data() + buffer.size(), x, std::chars_format::scientific);
  if (r.ec != std::errc{}) throw std::logic_error("python_float_repr: to_chars failed");
  std::string_view s(buffer.data(), static_cast<std::size_t>(r.ptr - buffer.data()));
  const bool negative = !s.empty() && s.front() == '-';
  if (negative) s.remove_prefix(1);
  const std::size_t e = s.find('e');
  const std::string_view mantissa = s.substr(0, e);
  int exponent = 0;
  const std::string_view exp_text = s.substr(e + 1);
  const char* first = exp_text.data();
  if (*first == '+') ++first;
  if (std::from_chars(first, exp_text.data() + exp_text.size(), exponent).ec != std::errc{}) {
    throw std::logic_error("python_float_repr: from_chars failed");
  }
  std::string digits;
  for (const char c : mantissa) {
    if (c != '.') digits += c;
  }
  const int decpt = exponent + 1;  // value = 0.digits * 10^decpt
  std::string out = negative ? "-" : "";
  const auto n = static_cast<int>(digits.size());
  if (decpt <= -4 || decpt > 16) {
    out += digits.front();
    if (n > 1) {
      out += '.';
      out += digits.substr(1);
    }
    out += exponent < 0 ? "e-" : "e+";
    const int magnitude = exponent < 0 ? -exponent : exponent;
    if (magnitude < 10) out += '0';
    out += std::to_string(magnitude);
  } else if (decpt <= 0) {
    out += "0.";
    out.append(static_cast<std::size_t>(-decpt), '0');
    out += digits;
  } else if (n <= decpt) {
    out += digits;
    out.append(static_cast<std::size_t>(decpt - n), '0');
    out += ".0";
  } else {
    out += digits.substr(0, static_cast<std::size_t>(decpt));
    out += '.';
    out += digits.substr(static_cast<std::size_t>(decpt));
  }
  return out;
}

namespace {

/// JSON text as json.dumps(..., indent=None, ensure_ascii=False) writes it.
class Writer {
 public:
  void begin_object() {
    separate();
    out_ += '{';
    first_.push_back(true);
  }
  void end_object() {
    out_ += '}';
    first_.pop_back();
  }
  void begin_array() {
    separate();
    out_ += '[';
    first_.push_back(true);
  }
  void end_array() {
    out_ += ']';
    first_.pop_back();
  }
  void key(std::string_view k) {
    separate();
    quoted(k);
    out_ += ": ";
    after_key_ = true;
  }
  void null() {
    separate();
    out_ += "null";
  }
  void boolean(bool b) {
    separate();
    out_ += b ? "true" : "false";
  }
  template <class I>
  void integer(I v) {
    separate();
    out_ += std::to_string(v);
  }
  /// A float scalar: repr, or {"float": "NaN" | "Infinity" | "-Infinity"}.
  void number(double v) {
    if (std::isfinite(v)) {
      separate();
      out_ += python_float_repr(v);
      return;
    }
    begin_object();
    key("float");
    string(special(v));
    end_object();
  }
  /// A float array element: repr, or the string "NaN", "Infinity", "-Infinity".
  void element(double v) {
    if (std::isfinite(v)) {
      separate();
      out_ += python_float_repr(v);
    } else {
      string(special(v));
    }
  }
  void string(std::string_view s) {
    separate();
    quoted(s);
  }
  template <class T>
  void optional_number(const std::optional<T>& v) {
    if (!v) {
      null();
    } else if constexpr (std::is_floating_point_v<T>) {
      number(*v);
    } else {
      integer(*v);
    }
  }

  std::string take() { return std::move(out_); }

 private:
  static std::string_view special(double v) {
    if (std::isnan(v)) return "NaN";
    return v > 0 ? "Infinity" : "-Infinity";
  }

  void separate() {
    if (after_key_) {
      after_key_ = false;
      return;
    }
    if (first_.empty()) return;
    if (!first_.back()) out_ += ", ";
    first_.back() = false;
  }

  /// A string as json.dumps with ensure_ascii=False: \" \\ \b \f \n \r \t, other characters below
  /// U+0020 as \u00xx with lowercase hex digits, everything else raw (UTF-8 passes through).
  void quoted(std::string_view s) {
    out_ += '"';
    for (const char c : s) {
      switch (c) {
        case '"':
          out_ += "\\\"";
          break;
        case '\\':
          out_ += "\\\\";
          break;
        case '\b':
          out_ += "\\b";
          break;
        case '\f':
          out_ += "\\f";
          break;
        case '\n':
          out_ += "\\n";
          break;
        case '\r':
          out_ += "\\r";
          break;
        case '\t':
          out_ += "\\t";
          break;
        default:
          if (static_cast<unsigned char>(c) < 0x20) {
            constexpr std::string_view kHex = "0123456789abcdef";
            out_ += "\\u00";
            out_ += kHex[static_cast<unsigned char>(c) >> 4U];
            out_ += kHex[static_cast<unsigned char>(c) & 0xFU];
          } else {
            out_ += c;
          }
      }
    }
    out_ += '"';
  }

  std::string out_;
  std::vector<bool> first_;
  bool after_key_ = false;
};

/// The envelope around `data` (results.encode_json).
template <class Data>
std::string envelope(std::string_view type, Data&& data) {
  Writer w;
  w.begin_object();
  w.key("format");
  w.string("raytatouille-result");
  w.key("schema_version");
  w.string(kResultSchemaVersion);
  w.key("type");
  w.string(type);
  w.key("data");
  data(w);
  w.end_object();
  return w.take();
}

/// A NumPy array {"dtype", "shape", "values"} with `count` values from get(i), in C order.
template <class Get>
void array(Writer& w,
           std::string_view dtype,
           std::initializer_list<std::size_t> shape,
           std::size_t count,
           Get&& get) {
  w.begin_object();
  w.key("dtype");
  w.string(dtype);
  w.key("shape");
  w.begin_array();
  for (const std::size_t n : shape) w.integer(n);
  w.end_array();
  w.key("values");
  w.begin_array();
  for (std::size_t i = 0; i < count; ++i) {
    using V = decltype(get(i));
    const V v = get(i);
    if constexpr (std::is_same_v<V, double>) {
      w.element(v);
    } else if constexpr (std::is_same_v<V, std::complex<double>>) {
      w.begin_array();
      w.element(v.real());
      w.element(v.imag());
      w.end_array();
    } else if constexpr (std::is_same_v<V, bool>) {
      w.boolean(v);
    } else {
      static_assert(std::is_integral_v<V>, "array: double, complex, bool or integer elements");
      w.integer(v);
    }
  }
  w.end_array();
  w.end_object();
}

/// A float64 array of shape (n,) from a span.
void f64(Writer& w, std::span<const double> v) {
  array(w, "float64", {v.size()}, v.size(), [&](std::size_t i) { return v[i]; });
}

/// A float64 array (3,) of a vector.
void vec3(Writer& w, const math::Vec3& v) {
  array(w, "float64", {3}, 3, [&](std::size_t i) { return v[static_cast<Eigen::Index>(i)]; });
}

/// A float64 array (2, 2), row-major, of a 2x2 matrix stored as std::array<double, 4>.
void matrix2(Writer& w, const std::array<double, 4>& m) {
  array(w, "float64", {2, 2}, 4, [&](std::size_t i) { return m[i]; });
}

/// A float64 (or complex128) square map (side, side) as the binding shapes it (bind_psf.cpp,
/// square()): `side` from the result (Psf::length, 2 Mtf::grid - 1), (0, 0) if empty.
template <class T>
void square(Writer& w, std::string_view dtype, const std::vector<T>& v, std::size_t side) {
  if (v.empty()) side = 0;
  if (side * side != v.size()) throw std::logic_error("result_json: a map is not side x side");
  array(w, dtype, {side, side}, v.size(), [&](std::size_t i) { return v[i]; });
}

void diagnostics(Writer& w, const std::vector<model::Diagnostic>& list) {
  // Diagnostic (results.py _FIELDS): severity, code, location, message.
  w.begin_array();
  for (const model::Diagnostic& d : list) {
    w.begin_object();
    w.key("severity");
    w.string(result_name(d.severity));
    w.key("code");
    w.string(d.code);
    w.key("location");
    w.string(d.location);
    w.key("message");
    w.string(d.message);
    w.end_object();
  }
  w.end_array();
}

void point2(Writer& w, const analysis::Point2& p) {
  w.begin_object();
  w.key("x");
  w.number(p.x);
  w.key("y");
  w.number(p.y);
  w.end_object();
}

void losses(Writer& w, const analysis::RayLosses& l) {
  // RayLosses: launched, by_status (a list of ints), worst_surface, worst_surface_count.
  w.begin_object();
  w.key("launched");
  w.integer(l.launched);
  w.key("by_status");
  w.begin_array();
  for (const std::size_t n : l.by_status) w.integer(n);
  w.end_array();
  w.key("worst_surface");
  w.optional_number(l.worst_surface);
  w.key("worst_surface_count");
  w.integer(l.worst_surface_count);
  w.end_object();
}

void psf_wavelengths(Writer& w, const std::vector<analysis::PsfWavelength>& v) {
  // PsfWavelengths (bind_psf.cpp): columns wavelength u16, weight, energy f64, length i32,
  // outside_field f64.
  const std::size_t n = v.size();
  w.begin_object();
  w.key("wavelength");
  array(w, "uint16", {n}, n, [&](std::size_t i) { return v[i].wavelength; });
  w.key("weight");
  array(w, "float64", {n}, n, [&](std::size_t i) { return v[i].weight; });
  w.key("energy");
  array(w, "float64", {n}, n, [&](std::size_t i) { return v[i].energy; });
  w.key("length");
  array(w, "int32", {n}, n, [&](std::size_t i) { return static_cast<std::int32_t>(v[i].length); });
  w.key("outside_field");
  array(w, "float64", {n}, n, [&](std::size_t i) { return v[i].outside_field; });
  w.end_object();
}

void seidel_terms(Writer& w, const paraxial::SeidelTerms& t) {
  w.begin_object();
  w.key("s1");
  w.number(t.s1);
  w.key("s2");
  w.number(t.s2);
  w.key("s3");
  w.number(t.s3);
  w.key("s4");
  w.number(t.s4);
  w.key("s5");
  w.number(t.s5);
  w.key("c_l");
  w.number(t.c_l);
  w.key("c_t");
  w.number(t.c_t);
  w.end_object();
}

void ray_start(Writer& w, const paraxial::RayStart& r) {
  w.begin_object();
  w.key("z");
  w.number(r.z);
  w.key("y");
  w.number(r.y);
  w.key("u");
  w.number(r.u);
  w.end_object();
}

void pupil(Writer& w, const std::optional<paraxial::Pupil>& p) {
  if (!p) {
    w.null();
    return;
  }
  w.begin_object();
  w.key("z");
  w.optional_number(p->z);
  w.key("diameter");
  w.optional_number(p->diameter);
  w.end_object();
}

void mtf_curve(Writer& w, const analysis::MtfCurve& c) {
  // MtfCurve: direction f64 (2,), frequency, modulation f64 (m,), otf complex128 (m,) or (0,).
  w.begin_object();
  w.key("direction");
  array(w, "float64", {2}, 2, [&](std::size_t i) { return c.direction[i]; });
  w.key("frequency");
  f64(w, c.frequency);
  w.key("modulation");
  f64(w, c.modulation);
  w.key("otf");
  array(w, "complex128", {c.otf.size()}, c.otf.size(), [&](std::size_t i) { return c.otf[i]; });
  w.end_object();
}

}  // namespace

std::string to_result_json(const trace::RayBatch& r) {
  // bind_trace.cpp; results.py _FIELDS["RayBatch"]: the columns, status as uint8, prt always
  // (complex128 (n, 3, 3), data[9 i + 3 row + col] = prt(row, col)[i]), q likewise or null
  // without Q, then wave_x/y/z and mode_index.
  return envelope("RayBatch", [&](Writer& w) {
    const std::size_t n = r.size();
    const auto column = [&](std::string_view name, std::span<const double> v) {
      w.key(name);
      f64(w, v);
    };
    w.begin_object();
    column("pos_x", r.pos_x());
    column("pos_y", r.pos_y());
    column("pos_z", r.pos_z());
    column("dir_x", r.dir_x());
    column("dir_y", r.dir_y());
    column("dir_z", r.dir_z());
    w.key("wl");
    array(w, "uint16", {n}, n, [&](std::size_t i) { return r.wl()[i]; });
    column("opl", r.opl());
    column("weight", r.weight());
    w.key("field");
    array(w, "uint16", {n}, n, [&](std::size_t i) { return r.field()[i]; });
    column("pupil_x", r.pupil_x());
    column("pupil_y", r.pupil_y());
    w.key("last_surface");
    array(w, "uint32", {n}, n, [&](std::size_t i) { return r.last_surface()[i]; });
    w.key("status");
    array(w, "uint8", {n}, n,
          [&](std::size_t i) { return static_cast<std::uint8_t>(r.status()[i]); });
    w.key("prt");
    array(w, "complex128", {n, 3, 3}, 9 * n, [&](std::size_t k) {
      return r.prt(static_cast<int>(k % 9 / 3), static_cast<int>(k % 3))[k / 9];
    });
    w.key("q");
    if (r.has_q()) {
      array(w, "float64", {n, 3, 3}, 9 * n, [&](std::size_t k) {
        return r.q(static_cast<int>(k % 9 / 3), static_cast<int>(k % 3))[k / 9];
      });
    } else {
      w.null();
    }
    column("wave_x", r.wave_x());
    column("wave_y", r.wave_y());
    column("wave_z", r.wave_z());
    column("mode_index", r.mode_index());
    w.end_object();
  });
}

std::string to_result_json(const trace::RayPaths& p) {
  // bind_trace.cpp: slots, n_events, ray_count, then the arrays (N rays, S slots); prt and q
  // null when empty.
  return envelope("RayPaths", [&](Writer& w) {
    const std::size_t rays = p.ray_count();
    const std::size_t s = p.slots;
    w.begin_object();
    w.key("slots");
    w.integer(p.slots);
    w.key("n_events");
    w.integer(p.event_surfaces.size());
    w.key("ray_count");
    w.integer(rays);
    w.key("ray_indices");
    array(w, "uint64", {rays}, rays,
          [&](std::size_t i) { return static_cast<std::uint64_t>(p.ray_indices[i]); });
    w.key("event_surfaces");
    array(w, "uint32", {p.event_surfaces.size()}, p.event_surfaces.size(),
          [&](std::size_t i) { return p.event_surfaces[i]; });
    w.key("position");
    array(w, "float64", {rays, s, 3}, p.position.size(),
          [&](std::size_t i) { return p.position[i]; });
    w.key("direction");
    array(w, "float64", {rays, s, 3}, p.direction.size(),
          [&](std::size_t i) { return p.direction[i]; });
    w.key("opl");
    array(w, "float64", {rays, s}, p.opl.size(), [&](std::size_t i) { return p.opl[i]; });
    w.key("weight");
    array(w, "float64", {rays, s}, p.weight.size(), [&](std::size_t i) { return p.weight[i]; });
    w.key("status");
    array(w, "uint8", {rays, s}, p.status.size(),
          [&](std::size_t i) { return static_cast<std::uint8_t>(p.status[i]); });
    w.key("count");
    array(w, "uint32", {rays}, p.count.size(), [&](std::size_t i) { return p.count[i]; });
    w.key("lost_at");
    array(w, "int32", {rays}, p.lost_at.size(), [&](std::size_t i) { return p.lost_at[i]; });
    w.key("polarization");
    w.boolean(p.polarization);
    w.key("prt");
    if (p.prt.empty()) {
      w.null();
    } else {
      array(w, "complex128", {rays, s, 3, 3}, p.prt.size(),
            [&](std::size_t i) { return p.prt[i]; });
    }
    w.key("q");
    if (p.q.empty()) {
      w.null();
    } else {
      array(w, "float64", {rays, s, 3, 3}, p.q.size(), [&](std::size_t i) { return p.q[i]; });
    }
    w.end_object();
  });
}

std::string to_result_json(const trace::TraceStats& t) {
  // TraceStats: rays, a list of the counts by status.
  return envelope("TraceStats", [&](Writer& w) {
    w.begin_object();
    w.key("rays");
    w.begin_array();
    for (const std::size_t n : t.rays) w.integer(n);
    w.end_array();
    w.end_object();
  });
}

std::string to_result_json(const paraxial::FirstOrder& f) {
  // bind_paraxial.cpp; results.py _FIELDS["FirstOrder"].
  return envelope("FirstOrder", [&](Writer& w) {
    w.begin_object();
    w.key("object_index");
    w.number(f.object_index);
    w.key("image_index");
    w.number(f.image_index);
    w.key("image_direction");
    w.integer(f.image_direction);
    w.key("power");
    w.number(f.power);
    const std::array<std::pair<std::string_view, const std::optional<double>*>, 12> optional{{
        {"efl", &f.efl},
        {"front_focal_length", &f.front_focal_length},
        {"rear_focal_length", &f.rear_focal_length},
        {"ffl", &f.ffl},
        {"bfl", &f.bfl},
        {"front_focal_z", &f.front_focal_z},
        {"rear_focal_z", &f.rear_focal_z},
        {"front_principal_z", &f.front_principal_z},
        {"rear_principal_z", &f.rear_principal_z},
        {"image_z", &f.image_z},
        {"lateral_magnification", &f.lateral_magnification},
        {"angular_magnification", &f.angular_magnification},
    }};
    for (const auto& [name, value] : optional) {
      w.key(name);
      w.optional_number(*value);
    }
    w.key("entrance_pupil");
    pupil(w, f.entrance_pupil);
    w.key("exit_pupil");
    pupil(w, f.exit_pupil);
    w.end_object();
  });
}

std::string to_result_json(const paraxial::Seidel& s) {
  // bind_analysis.cpp: surfaces (columns over the events), sum, lagrange, marginal, chief,
  // chromatic.
  return envelope("Seidel", [&](Writer& w) {
    const std::vector<paraxial::SeidelSurface>& v = s.surfaces;
    const std::size_t n = v.size();
    const auto column = [&](std::string_view name, auto get) {
      w.key(name);
      array(w, "float64", {n}, n, [&](std::size_t i) { return get(v[i]); });
    };
    w.begin_object();
    w.key("surfaces");
    w.begin_object();
    w.key("surface");
    array(w, "uint32", {n}, n, [&](std::size_t i) { return v[i].surface; });
    column("y", [](const paraxial::SeidelSurface& q) { return q.y; });
    column("y_bar", [](const paraxial::SeidelSurface& q) { return q.y_bar; });
    column("a", [](const paraxial::SeidelSurface& q) { return q.a; });
    column("a_bar", [](const paraxial::SeidelSurface& q) { return q.a_bar; });
    column("lagrange", [](const paraxial::SeidelSurface& q) { return q.lagrange; });
    column("s1", [](const paraxial::SeidelSurface& q) { return q.terms.s1; });
    column("s2", [](const paraxial::SeidelSurface& q) { return q.terms.s2; });
    column("s3", [](const paraxial::SeidelSurface& q) { return q.terms.s3; });
    column("s4", [](const paraxial::SeidelSurface& q) { return q.terms.s4; });
    column("s5", [](const paraxial::SeidelSurface& q) { return q.terms.s5; });
    column("c_l", [](const paraxial::SeidelSurface& q) { return q.terms.c_l; });
    column("c_t", [](const paraxial::SeidelSurface& q) { return q.terms.c_t; });
    w.end_object();
    w.key("sum");
    seidel_terms(w, s.sum);
    w.key("lagrange");
    w.number(s.lagrange);
    w.key("marginal");
    ray_start(w, s.marginal);
    w.key("chief");
    ray_start(w, s.chief);
    w.key("chromatic");
    if (s.chromatic) {
      w.begin_object();
      w.key("first");
      w.integer(s.chromatic->first);
      w.key("second");
      w.integer(s.chromatic->second);
      w.end_object();
    } else {
      w.null();
    }
    w.end_object();
  });
}

std::string to_result_json(const analysis::SpotDiagram& d) {
  // bind_analysis.cpp; _FIELDS["SpotDiagram"]: x, y, wavelengths, weight from the points.
  return envelope("SpotDiagram", [&](Writer& w) {
    const std::size_t n = d.points.size();
    w.begin_object();
    w.key("field");
    w.integer(d.field);
    w.key("wavelength");
    w.optional_number(d.wavelength);
    w.key("image_surface");
    w.integer(d.image_surface);
    w.key("chief");
    point2(w, d.chief);
    w.key("stats");
    w.begin_object();
    w.key("centroid");
    point2(w, d.stats.centroid);
    w.key("rms_centroid");
    w.number(d.stats.rms_centroid);
    w.key("rms_chief");
    w.number(d.stats.rms_chief);
    w.key("geo_centroid");
    w.number(d.stats.geo_centroid);
    w.key("geo_chief");
    w.number(d.stats.geo_chief);
    w.end_object();
    w.key("rays_launched");
    w.integer(d.rays_launched);
    w.key("rays_arrived");
    w.integer(d.rays_arrived);
    w.key("vignetted_fraction");
    w.number(d.vignetted_fraction);
    w.key("losses");
    losses(w, d.losses);
    w.key("warnings");
    diagnostics(w, d.warnings);
    w.key("x");
    array(w, "float64", {n}, n, [&](std::size_t i) { return d.points[i].x; });
    w.key("y");
    array(w, "float64", {n}, n, [&](std::size_t i) { return d.points[i].y; });
    w.key("wavelengths");
    array(w, "uint16", {n}, n, [&](std::size_t i) { return d.points[i].wavelength; });
    w.key("weight");
    array(w, "float64", {n}, n, [&](std::size_t i) { return d.points[i].weight; });
    w.end_object();
  });
}

std::string to_result_json(const analysis::OpdMap& m) {
  // bind_analysis.cpp; _FIELDS["OpdMap"], ReferenceSphere, OpdPoints (status uint8).
  return envelope("OpdMap", [&](Writer& w) {
    const std::size_t n = m.points.size();
    w.begin_object();
    w.key("field");
    w.integer(m.field);
    w.key("wavelength");
    w.integer(m.wavelength);
    w.key("sphere");
    w.begin_object();
    w.key("centre");
    vec3(w, m.sphere.centre);
    w.key("radius");
    w.number(m.sphere.radius);
    w.end_object();
    w.key("points");
    w.begin_object();
    w.key("px");
    array(w, "float64", {n}, n, [&](std::size_t i) { return m.points[i].px; });
    w.key("py");
    array(w, "float64", {n}, n, [&](std::size_t i) { return m.points[i].py; });
    w.key("w");
    array(w, "float64", {n}, n, [&](std::size_t i) { return m.points[i].w; });
    w.key("status");
    array(w, "uint8", {n}, n,
          [&](std::size_t i) { return static_cast<std::uint8_t>(m.points[i].status); });
    w.end_object();
    w.key("rms");
    w.number(m.rms);
    w.key("pv");
    w.number(m.pv);
    w.key("arrived");
    w.integer(m.arrived);
    w.key("vignetted");
    w.integer(m.vignetted);
    w.key("losses");
    losses(w, m.losses);
    w.key("warnings");
    diagnostics(w, m.warnings);
    w.end_object();
  });
}

std::string to_result_json(const analysis::ZernikeFit& f) {
  // bind_pupil.cpp; _FIELDS["ZernikeFit"]: path as its index.
  return envelope("ZernikeFit", [&](Writer& w) {
    w.begin_object();
    w.key("path");
    w.integer(f.path.index);
    w.key("field");
    w.integer(f.field);
    w.key("wavelength");
    w.integer(f.wavelength);
    w.key("wavelength_um");
    w.number(f.wavelength_um);
    w.key("reference_wavelength_um");
    w.number(f.reference_wavelength_um);
    w.key("coefficients");
    f64(w, f.coefficients);
    w.key("residual_rms");
    w.number(f.residual_rms);
    w.key("sigma");
    w.number(f.sigma);
    w.key("sigma_bf");
    w.number(f.sigma_bf);
    w.key("strehl");
    w.number(f.strehl);
    w.key("cells");
    w.integer(f.cells);
    w.key("arrived_fraction");
    w.number(f.arrived_fraction);
    w.key("warnings");
    diagnostics(w, f.warnings);
    w.end_object();
  });
}

std::string to_result_json(const analysis::Mtf& m) {
  // bind_psf.cpp; _FIELDS["Mtf"]: map and otf_map square or (0, 0).
  return envelope("Mtf", [&](Writer& w) {
    w.begin_object();
    w.key("field");
    w.integer(m.field);
    w.key("grid");
    w.integer(m.grid);
    w.key("centre");
    vec3(w, m.centre);
    w.key("e_x");
    vec3(w, m.e_x);
    w.key("e_y");
    vec3(w, m.e_y);
    w.key("tangential");
    mtf_curve(w, m.tangential);
    w.key("sagittal");
    mtf_curve(w, m.sagittal);
    w.key("frequency_matrix");
    matrix2(w, m.frequency_matrix);
    w.key("map");
    const auto side = static_cast<std::size_t>(2 * m.grid - 1);
    square(w, "float64", m.map, side);
    w.key("otf_map");
    square(w, "complex128", m.otf_map, side);
    w.key("edge_energy_deficit");
    w.number(m.edge_energy_deficit);
    w.key("wavelengths");
    psf_wavelengths(w, m.wavelengths);
    w.key("losses");
    losses(w, m.losses);
    w.key("warnings");
    diagnostics(w, m.warnings);
    w.end_object();
  });
}

std::string to_result_json(const analysis::Psf& p) {
  // bind_psf.cpp; _FIELDS["Psf"]: intensity (length, length) or (0, 0).
  return envelope("Psf", [&](Writer& w) {
    w.begin_object();
    w.key("field");
    w.integer(p.field);
    w.key("grid");
    w.integer(p.grid);
    w.key("length");
    w.integer(p.length);
    w.key("centre");
    vec3(w, p.centre);
    w.key("e_x");
    vec3(w, p.e_x);
    w.key("e_y");
    vec3(w, p.e_y);
    w.key("lattice");
    matrix2(w, p.lattice);
    w.key("intensity");
    square(w, "float64", p.intensity, static_cast<std::size_t>(p.length));
    w.key("strehl");
    w.number(p.strehl);
    w.key("peak");
    w.begin_object();
    w.key("value");
    w.number(p.peak.value);
    w.key("k1");
    w.integer(p.peak.k1);
    w.key("k2");
    w.integer(p.peak.k2);
    w.key("x");
    w.number(p.peak.x);
    w.key("y");
    w.number(p.peak.y);
    w.end_object();
    w.key("edge_energy_deficit");
    w.number(p.edge_energy_deficit);
    w.key("polychromatic_scale_error");
    w.number(p.polychromatic_scale_error);
    w.key("wavelengths");
    psf_wavelengths(w, p.wavelengths);
    w.key("losses");
    losses(w, p.losses);
    w.key("warnings");
    diagnostics(w, p.warnings);
    w.end_object();
  });
}

std::string to_result_json(const optim::OptimResult& r) {
  // bind_optim.cpp; _FIELDS["OptimResult"] (without the system): history as OptimIterations
  // columns, then operands, generators, variables, diagnostics, counts.
  return envelope("OptimResult", [&](Writer& w) {
    const std::vector<optim::OptimIteration>& h = r.history;
    const std::size_t n = h.size();
    w.begin_object();
    w.key("status");
    w.string(result_name(r.status));
    w.key("patch");
    w.string(r.patch);
    w.key("history");
    w.begin_object();
    w.key("k");
    array(w, "int32", {n}, n, [&](std::size_t i) { return static_cast<std::int32_t>(h[i].k); });
    w.key("phi");
    array(w, "float64", {n}, n, [&](std::size_t i) { return h[i].phi; });
    w.key("mu");
    array(w, "float64", {n}, n, [&](std::size_t i) { return h[i].mu; });
    w.key("rho");
    array(w, "float64", {n}, n, [&](std::size_t i) { return h[i].rho; });
    w.key("step_norm");
    array(w, "float64", {n}, n, [&](std::size_t i) { return h[i].step_norm; });
    w.key("accepted");
    array(w, "bool", {n}, n, [&](std::size_t i) { return h[i].accepted; });
    w.key("evaluations");
    array(w, "uint64", {n}, n,
          [&](std::size_t i) { return static_cast<std::uint64_t>(h[i].evaluations); });
    w.end_object();
    w.key("operands");
    w.begin_array();
    for (const optim::OperandValue& o : r.operands) {
      w.begin_object();
      w.key("pointer");
      w.string(o.pointer);
      w.key("value");
      w.number(o.value);
      w.key("target");
      w.number(o.target);
      w.key("weight");
      w.number(o.weight);
      w.key("contribution");
      w.number(o.contribution);
      w.end_object();
    }
    w.end_array();
    w.key("generators");
    w.begin_array();
    for (const optim::GeneratorValue& g : r.generators) {
      w.begin_object();
      w.key("pointer");
      w.string(g.pointer);
      w.key("rms");
      w.number(g.rms);
      w.key("weight");
      w.number(g.weight);
      w.key("contribution");
      w.number(g.contribution);
      w.key("rays_launched");
      w.integer(g.rays_launched);
      w.key("rays_lost");
      w.integer(g.rays_lost);
      w.end_object();
    }
    w.end_array();
    w.key("variables");
    w.begin_array();
    for (const optim::VariableValue& v : r.variables) {
      w.begin_object();
      w.key("pointer");
      w.string(v.pointer);
      w.key("row");
      w.string(v.row);
      w.key("configuration");
      w.optional_number(v.configuration);
      w.key("start");
      w.number(v.start);
      w.key("end");
      w.number(v.end);
      w.key("changed");
      w.boolean(v.changed);
      w.key("at_bound");
      w.boolean(v.at_bound);
      w.end_object();
    }
    w.end_array();
    w.key("diagnostics");
    diagnostics(w, r.diagnostics);
    w.key("iterations");
    w.integer(r.iterations);
    w.key("evaluations");
    w.integer(r.evaluations);
    w.key("failed_evaluations");
    w.integer(r.failed_evaluations);
    w.end_object();
  });
}

}  // namespace rtt::io
