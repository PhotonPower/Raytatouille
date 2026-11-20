// rtt trace and rtt analyze (ADR 0052, points 2 to 5): load, validate and compile the system,
// run one computation with the options of the Python API, write the result in the format
// raytatouille-result (rtt::io::to_result_json, byte for byte as to_json in Python) or as CSV.

#include "compute.hpp"

#include <oneapi/tbb/task_arena.h>

#include <charconv>
#include <climits>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <fstream>
#include <functional>
#include <map>
#include <numeric>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <variant>
#include <vector>

#include "rtt/analysis/nominal_reference.hpp"
#include "rtt/analysis/opd.hpp"
#include "rtt/analysis/psf.hpp"
#include "rtt/analysis/pupil.hpp"
#include "rtt/analysis/spot.hpp"
#include "rtt/analysis/zernike.hpp"
#include "rtt/coating/catalog.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/io/result_json.hpp"
#include "rtt/io/result_names.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/parameters.hpp"
#include "rtt/model/validate.hpp"
#include "rtt/paraxial/paraxial.hpp"
#include "rtt/paraxial/seidel.hpp"
#include "rtt/trace/ray_paths.hpp"
#include "rtt/trace/sequential.hpp"
#include "rtt/trace/sources.hpp"

namespace rtt::cli {

namespace {

constexpr int kOk = 0;
constexpr int kFailed = 1;
constexpr int kUsage = 2;

/// A wrong call (exit 2): command, kind, option, missing or unreadable value, a range the
/// command line checks itself, or --format csv without CSV.
class UsageError : public std::runtime_error {
 public:
  using std::runtime_error::runtime_error;
};

/// The system cannot be computed (exit 1); the diagnostics are printed already.
class Failed : public std::runtime_error {
 public:
  using std::runtime_error::runtime_error;
};

// ------------------------------------------------------------------------------- options ---

/// One option of a command: its name, whether it takes a value, whether it may repeat, help.
struct OptionSpec {
  std::string_view name;
  std::string_view value;  ///< name of the value in the help ("" for a flag)
  bool repeat = false;
  std::string_view help;
};

/// The parsed arguments: operands and options in the order given.
struct Arguments {
  std::vector<std::string> operands;
  std::map<std::string, std::vector<std::string>, std::less<>> values;
  std::set<std::string, std::less<>> flags;

  [[nodiscard]] const std::string* value(std::string_view name) const {
    const auto it = values.find(name);
    return it == values.end() ? nullptr : &it->second.back();
  }
  [[nodiscard]] bool flag(std::string_view name) const { return flags.contains(name); }
};

Arguments parse(const std::vector<std::string>& args,
                std::size_t first,
                const std::vector<OptionSpec>& specs) {
  Arguments a;
  for (std::size_t i = first; i < args.size(); ++i) {
    const std::string& arg = args[i];
    if (!arg.starts_with("-") || arg == "-") {
      a.operands.push_back(arg);
      continue;
    }
    const std::string name = arg == "-o" ? "--output" : arg;
    const OptionSpec* spec = nullptr;
    for (const OptionSpec& s : specs) {
      if (s.name == name) spec = &s;
    }
    if (spec == nullptr) throw UsageError("unknown option " + arg);
    if (spec->value.empty()) {
      a.flags.insert(name);
      continue;
    }
    if (i + 1 >= args.size()) throw UsageError("option " + arg + " needs a value");
    std::vector<std::string>& v = a.values[name];
    if (!v.empty() && !spec->repeat) throw UsageError("option " + arg + " given twice");
    v.push_back(args[++i]);
  }
  return a;
}

template <class T>
T integer(const std::string& text, std::string_view option, T min, T max) {
  T v{};
  const auto [end, ec] = std::from_chars(text.data(), text.data() + text.size(), v);
  if (ec != std::errc{} || end != text.data() + text.size() || v < min || v > max) {
    throw UsageError("option " + std::string(option) + ": '" + text + "' is not an integer in [" +
                     std::to_string(min) + ", " + std::to_string(max) + "]");
  }
  return v;
}

double number(const std::string& text, std::string_view option) {
  double v = 0.0;
  const auto [end, ec] = std::from_chars(text.data(), text.data() + text.size(), v);
  if (ec != std::errc{} || end != text.data() + text.size() || !std::isfinite(v)) {
    throw UsageError("option " + std::string(option) + ": '" + text + "' is not a finite number");
  }
  return v;
}

std::vector<std::string> split(const std::string& text, char sep) {
  std::vector<std::string> parts;
  std::string part;
  std::istringstream in(text);
  while (std::getline(in, part, sep)) parts.push_back(part);
  if (!text.empty() && text.back() == sep) parts.emplace_back();
  return parts;
}

/// The pupil sampling of the shorthand of the Python API (raytatouille.analysis.sampling):
/// hexapolar:N, grid:N, fan_x:N, fan_y:N, random:N[:SEED], single:PX,PY.
trace::PupilSampling sampling(const std::string& text) {
  const std::vector<std::string> parts = split(text, ':');
  const std::string kind = parts.empty() ? std::string() : parts[0];
  constexpr std::string_view kOption = "--rays";
  constexpr int kIntMax = 2147483647;
  if (kind == "hexapolar" && parts.size() == 2) {
    return trace::HexapolarPupil{integer<int>(parts[1], kOption, 0, kIntMax)};
  }
  if ((kind == "grid" || kind == "fan_x" || kind == "fan_y") && parts.size() == 2) {
    const int n = integer<int>(parts[1], kOption, 1, kIntMax);
    if (kind == "grid") return trace::GridPupil{n};
    if (kind == "fan_x") return trace::FanXPupil{n};
    return trace::FanYPupil{n};
  }
  if (kind == "random" && (parts.size() == 2 || parts.size() == 3)) {
    constexpr std::uint64_t kMax = UINT64_MAX;
    const auto count = integer<std::uint64_t>(parts[1], kOption, 0, kMax);
    const std::uint64_t seed =
        parts.size() == 3 ? integer<std::uint64_t>(parts[2], kOption, 0, kMax) : 0;
    return trace::RandomPupil{count, seed};
  }
  if (kind == "single" && parts.size() == 2) {
    const std::vector<std::string> xy = split(parts[1], ',');
    if (xy.size() == 2)
      return trace::SinglePupilPoint{number(xy[0], kOption), number(xy[1], kOption)};
  }
  throw UsageError("option --rays: '" + text +
                   "' is not hexapolar:N, grid:N, fan_x:N, fan_y:N, random:N[:SEED] or "
                   "single:PX,PY");
}

trace::Aiming aiming(const Arguments& a) {
  const std::string* v = a.value("--aiming");
  if (v == nullptr || *v == "real") return trace::Aiming::Real;
  if (*v == "paraxial") return trace::Aiming::Paraxial;
  throw UsageError("option --aiming: '" + *v + "' is not real or paraxial");
}

// ---------------------------------------------------------------------------- the system ---

/// The loaded, validated and compiled system with the libraries and the chosen path.
struct Prepared {
  std::string file;
  model::System system;
  material::MaterialLibrary materials;
  coating::CoatingLibrary coatings;
  std::size_t configuration = 0;
  std::optional<compile::CompiledSystem> compiled;
  compile::PathId path{};
};

void print(std::ostream& err, const std::string& file, const std::vector<model::Diagnostic>& list) {
  for (const model::Diagnostic& d : list) err << file << ": " << model::to_string(d) << "\n";
}

/// Loads, validates (errors: exit 1 before computing) and compiles `file` with the libraries
/// and the configuration of the options; prints the warnings of reading and compiling.
void prepare(Prepared& p, const std::string& file, const Arguments& a, std::ostream& err) {
  p.file = file;
  std::vector<model::Diagnostic> read_warnings;
  try {
    p.system = io::load_system(file, read_warnings);
  } catch (const std::exception& e) {
    err << file << ": error " << e.what() << "\n";
    throw Failed("read");
  }
  print(err, file, read_warnings);
  const std::vector<model::Diagnostic> diagnostics = model::validate(p.system);
  print(err, file, diagnostics);
  if (model::has_errors(diagnostics)) throw Failed("validate");
  if (const auto it = a.values.find("--catalog"); it != a.values.end()) {
    for (const std::string& c : it->second) p.materials.add_catalog(c);
  }
  if (const auto it = a.values.find("--coatings"); it != a.values.end()) {
    for (const std::string& c : it->second) p.coatings.add_catalog(c);
  }
  if (const std::string* c = a.value("--config")) {
    if (const std::optional<std::size_t> i = model::find_configuration(p.system, *c)) {
      p.configuration = *i;
    } else {
      const std::size_t n = std::max<std::size_t>(p.system.configurations.size(), 1);
      p.configuration = integer<std::size_t>(*c, "--config", 0, n - 1);
    }
  }
  try {
    p.compiled.emplace(compile::compile(p.system, p.materials, p.coatings, p.configuration));
  } catch (const compile::CompileError& e) {
    print(err, file, e.diagnostics());
    throw Failed("compile");
  }
  print(err, file, p.compiled->diagnostics());
  const compile::CompiledSystem& cs = *p.compiled;
  if (const std::string* name = a.value("--path")) {
    if (const std::optional<compile::PathId> id = cs.find_path(*name)) {
      p.path = *id;
    } else {
      p.path = compile::PathId{integer<std::uint32_t>(
          *name, "--path", 0, static_cast<std::uint32_t>(cs.paths().size() - 1))};
    }
  }
}

std::uint16_t field(const Prepared& p, const Arguments& a) {
  const std::string* v = a.value("--field");
  if (v == nullptr) return 0;
  const auto n = static_cast<std::uint16_t>(p.compiled->fields().points.size());
  return integer<std::uint16_t>(*v, "--field", 0, static_cast<std::uint16_t>(n - 1));
}

std::optional<std::uint16_t> wavelength(const Prepared& p, const Arguments& a) {
  const std::string* v = a.value("--wavelength");
  if (v == nullptr) return std::nullopt;
  const auto n = static_cast<std::uint16_t>(p.compiled->wavelengths_um().size());
  return integer<std::uint16_t>(*v, "--wavelength", 0, static_cast<std::uint16_t>(n - 1));
}

std::optional<std::vector<std::uint16_t>> wavelengths(const Prepared& p, const Arguments& a) {
  const std::string* v = a.value("--wavelengths");
  if (v == nullptr) return std::nullopt;
  const auto n = static_cast<std::uint16_t>(p.compiled->wavelengths_um().size());
  std::vector<std::uint16_t> list;
  for (const std::string& s : split(*v, ',')) {
    list.push_back(
        integer<std::uint16_t>(s, "--wavelengths", 0, static_cast<std::uint16_t>(n - 1)));
  }
  return list;
}

/// The nominal reference of --reference (ADR 0041, point 3): the nominal system compiled with
/// the same libraries and configuration, on the same path.
std::optional<analysis::NominalReference> reference(const Prepared& p,
                                                    const Arguments& a,
                                                    std::ostream& err,
                                                    std::optional<compile::CompiledSystem>& keep) {
  const std::string* file = a.value("--reference");
  if (file == nullptr) return std::nullopt;
  Prepared nominal;
  Arguments same = a;
  same.values.erase("--reference");
  prepare(nominal, *file, same, err);
  keep = std::move(nominal.compiled);
  return analysis::nominal_reference(*keep, p.path);
}

// -------------------------------------------------------------------------------- output ---

/// Writes `text` to `file` through a temporary file in the same directory and a rename, so that
/// no half file remains (ADR 0052, point 2), or to `out` without a file.
void emit(const std::string& text, const std::string* file, std::ostream& out) {
  if (file == nullptr) {
    out << text;
    out.flush();
    return;
  }
  const std::filesystem::path target(*file);
  std::filesystem::path temporary = target;
  temporary += ".rtt-partial";
  {
    std::ofstream f(temporary, std::ios::binary | std::ios::trunc);
    f << text;
    f.close();
    if (!f) throw std::runtime_error("cannot write " + *file);
  }
  std::error_code ec;
  std::filesystem::rename(temporary, target, ec);
  if (ec) {
    std::filesystem::remove(temporary, ec);
    throw std::runtime_error("cannot write " + *file);
  }
}

/// A CSV field: a number as Python's repr, NaN as an empty field (ADR 0052, point 5).
std::string csv(double v) {
  return std::isfinite(v) ? io::python_float_repr(v) : std::string();
}

std::string trace_csv(const trace::RayBatch& r) {
  std::ostringstream s;
  s << "pos_x,pos_y,pos_z,dir_x,dir_y,dir_z,wl,opl,weight,field,pupil_x,pupil_y,last_surface,"
       "status,wave_x,wave_y,wave_z,mode_index\n";
  for (std::size_t i = 0; i < r.size(); ++i) {
    s << csv(r.pos_x()[i]) << ',' << csv(r.pos_y()[i]) << ',' << csv(r.pos_z()[i]) << ','
      << csv(r.dir_x()[i]) << ',' << csv(r.dir_y()[i]) << ',' << csv(r.dir_z()[i]) << ','
      << r.wl()[i] << ',' << csv(r.opl()[i]) << ',' << csv(r.weight()[i]) << ',' << r.field()[i]
      << ',' << csv(r.pupil_x()[i]) << ',' << csv(r.pupil_y()[i]) << ',' << r.last_surface()[i]
      << ',' << io::result_name(r.status()[i]) << ',' << csv(r.wave_x()[i]) << ','
      << csv(r.wave_y()[i]) << ',' << csv(r.wave_z()[i]) << ',' << csv(r.mode_index()[i]) << '\n';
  }
  return s.str();
}

std::string paths_csv(const trace::RayPaths& p) {
  std::ostringstream s;
  s << "ray,slot,status,x,y,z,dx,dy,dz,opl\n";
  for (std::size_t r = 0; r < p.ray_count(); ++r) {
    for (std::size_t k = 0; k < p.slots; ++k) {
      const std::size_t at = r * p.slots + k;
      s << p.ray_indices[r] << ',' << k << ',' << io::result_name(p.status[at]);
      for (std::size_t c = 0; c < 3; ++c) s << ',' << csv(p.position[3 * at + c]);
      for (std::size_t c = 0; c < 3; ++c) s << ',' << csv(p.direction[3 * at + c]);
      s << ',' << csv(p.opl[at]) << '\n';
    }
  }
  return s.str();
}

std::string spot_csv(const analysis::SpotDiagram& d) {
  std::ostringstream s;
  const bool poly = !d.wavelength.has_value();
  s << (poly ? "x,y,wavelength,weight\n" : "x,y\n");
  for (const analysis::SpotPoint& q : d.points) {
    s << csv(q.x) << ',' << csv(q.y);
    if (poly) s << ',' << q.wavelength << ',' << csv(q.weight);
    s << '\n';
  }
  return s.str();
}

std::string mtf_csv(const analysis::Mtf& m) {
  std::ostringstream s;
  const bool phase = !m.tangential.otf.empty();
  s << (phase ? "direction,frequency,modulation,phase\n" : "direction,frequency,modulation\n");
  const auto curve = [&](std::string_view name, const analysis::MtfCurve& c) {
    for (std::size_t i = 0; i < c.frequency.size(); ++i) {
      s << name << ',' << csv(c.frequency[i]) << ',' << csv(c.modulation[i]);
      if (phase) s << ',' << csv(std::arg(c.otf[i]));
      s << '\n';
    }
  };
  curve("tangential", m.tangential);
  curve("sagittal", m.sagittal);
  return s.str();
}

// ---------------------------------------------------------------------------- the commands ---

const std::vector<OptionSpec>& common_options() {
  static const std::vector<OptionSpec> specs = {
      {"--catalog", "file.agf", true, "glass catalogue (AGF), repeatable"},
      {"--coatings", "dir", true, "coating catalogue, repeatable"},
      {"--config", "index|name", false, "configuration (default 0)"},
      {"--path", "name|index", false, "path; a path name, otherwise a 0-based index (default 0)"},
      {"--threads", "n", false, "worker threads (default: all)"},
      {"--output", "file", false, "write the result to file (-o), otherwise to stdout"},
      {"--format", "json|csv", false, "json (default) or csv where the command has CSV"},
      {"--help", "", false, "this help"},
  };
  return specs;
}

struct Kind {
  std::string_view name;
  std::string_view help;
  std::vector<OptionSpec> options;
  bool csv = false;
};

const std::vector<Kind>& analyze_kinds() {
  static const OptionSpec field{"--field", "index", false, "field index (default 0)"};
  static const OptionSpec wl{"--wavelength", "index", false,
                             "wavelength index (default: the reference; spot: all, weighted)"};
  static const OptionSpec wls{"--wavelengths", "i,j,...", false,
                              "wavelength indices (default: all with their weights)"};
  static const OptionSpec aim{"--aiming", "real|paraxial", false, "ray aiming (default real)"};
  static const OptionSpec ref{"--reference", "nominal.rtt.json", false,
                              "aim with the reference of a nominal system (ADR 0041, point 3)"};
  static const OptionSpec focus{"--focus-shift", "mm", false, "shift of the image plane, mm"};
  static const OptionSpec apod{"--apodization", "value", false, "Gaussian apodization (default 0)"};
  static const std::vector<Kind> kinds = {
      {"first-order", "first-order properties (rt.paraxial.first_order)", {wl}},
      {"seidel",
       "Seidel sums (rt.paraxial.seidel)",
       {wl, {"--pair", "i,j", false, "wavelength pair for the colour sums"}}},
      {"spot",
       "spot diagram (rt.analysis.spot); CSV: x,y (polychromatic also wavelength,weight)",
       {field, wl, aim, ref, {"--rays", "sampling", false, "pupil sampling (default hexapolar:6)"}},
       true},
      {"opd",
       "OPD map (rt.analysis.opd_map)",
       {field, wl, aim, ref, focus, {"--grid", "n", false, "pupil grid (default 33)"}}},
      {"zernike",
       "Zernike fit of the pupil function (rt.analysis.pupil_function, zernike_fit)",
       {field,
        wl,
        aim,
        focus,
        apod,
        {"--grid", "n", false, "pupil grid (default 64)"},
        {"--terms", "n", false, "Noll terms (default 37)"}}},
      {"psf",
       "diffraction PSF (rt.analysis.psf)",
       {field,
        wls,
        aim,
        focus,
        apod,
        {"--grid", "n", false, "pupil grid (default 64)"},
        {"--padding", "q", false, "FFT padding (default 4)"}}},
      {"mtf",
       "MTF (rt.analysis.mtf); CSV: direction,frequency,modulation (with --complex-otf phase)",
       {field,
        wls,
        aim,
        focus,
        apod,
        {"--grid", "n", false, "pupil grid (default 64)"},
        {"--map", "", false, "the 2D MTF map (one wavelength)"},
        {"--complex-otf", "", false, "the complex OTF next to the modulation"}},
       true},
  };
  return kinds;
}

const std::vector<OptionSpec>& trace_options() {
  static const std::vector<OptionSpec> specs = {
      {"--field", "index", false, "field index (default: all fields)"},
      {"--wavelength", "index", false, "wavelength index (default: the reference)"},
      {"--rays", "sampling", false, "pupil sampling (default hexapolar:6)"},
      {"--aiming", "real|paraxial", false, "ray aiming (default real)"},
      {"--paths", "", false, "write RayPaths (all events) instead of the RayBatch"},
      {"--max-recorded-rays", "n", false, "limit of rays recorded with --paths"},
      {"--polarization", "", false, "record P and Q per slot with --paths"},
      {"--stats", "file", false, "write TraceStats to file"},
  };
  return specs;
}

void help(std::ostream& os,
          std::string_view usage,
          std::string_view what,
          const std::vector<OptionSpec>& own) {
  os << "usage: " << usage << "\n" << what << "\n\noptions:\n";
  const auto line = [&](const OptionSpec& s) {
    std::string left = "  " + std::string(s.name == "--output" ? "-o, --output" : s.name);
    if (!s.value.empty()) left += " <" + std::string(s.value) + ">";
    os << left << std::string(left.size() < 34 ? 34 - left.size() : 1, ' ') << s.help << "\n";
  };
  for (const OptionSpec& s : own) line(s);
  for (const OptionSpec& s : common_options()) line(s);
}

std::vector<OptionSpec> with_common(std::vector<OptionSpec> own) {
  for (const OptionSpec& s : common_options()) own.push_back(s);
  return own;
}

/// Runs `body` in a TBB arena of --threads workers, or as is.
template <class F>
void in_arena(const Arguments& a, F&& body) {
  if (const std::string* t = a.value("--threads")) {
    oneapi::tbb::task_arena arena(integer<int>(*t, "--threads", 1, 65535));
    arena.execute(body);
  } else {
    body();
  }
}

std::string format_of(const Arguments& a, bool csv_allowed) {
  const std::string* f = a.value("--format");
  if (f == nullptr || *f == "json") return "json";
  if (*f != "csv") throw UsageError("option --format: '" + *f + "' is not json or csv");
  if (!csv_allowed) throw UsageError("option --format csv: this command has no CSV");
  return "csv";
}

int trace_command(const std::vector<std::string>& args, std::ostream& out, std::ostream& err) {
  const std::vector<OptionSpec> specs = with_common(trace_options());
  const Arguments a = parse(args, 1, specs);
  if (a.flag("--help")) {
    help(out, "rtt trace <system.rtt.json> [options]",
         "Rays of a pupil sampling (rt.trace.make_rays) traced through a path (rt.trace.trace).\n"
         "Output: RayBatch, or RayPaths with --paths; CSV: one row per ray (per ray and event).",
         trace_options());
    return kOk;
  }
  if (a.operands.size() != 1) throw UsageError("rtt trace needs one system file");
  const std::string format = format_of(a, true);
  Prepared p;
  prepare(p, a.operands[0], a, err);
  const compile::CompiledSystem& cs = *p.compiled;
  std::vector<std::uint16_t> fields;
  if (a.value("--field") != nullptr) {
    fields = {field(p, a)};
  } else {
    // One by one: resize and iota trip a false -Wnull-dereference in GCC 13 at -O3.
    const std::size_t n = cs.fields().points.size();
    fields.reserve(n);
    for (std::size_t i = 0; i < n; ++i) fields.push_back(static_cast<std::uint16_t>(i));
  }
  const std::uint16_t wl = wavelength(p, a).value_or(cs.reference_wavelength());
  const trace::PupilSampling s =
      a.value("--rays") != nullptr ? sampling(*a.value("--rays")) : trace::HexapolarPupil{6};
  const bool record = a.flag("--paths");
  if (!record && (a.flag("--polarization") || a.value("--max-recorded-rays") != nullptr)) {
    throw UsageError("--polarization and --max-recorded-rays need --paths");
  }
  std::string text;
  std::string stats_text;
  in_arena(a, [&] {
    trace::RayBatch rays = trace::make_rays(cs, p.path, fields, wl, s, aiming(a));
    trace::TraceStats stats;
    if (record) {
      trace::RayPaths paths;
      paths.polarization = a.flag("--polarization");
      const std::size_t limit = a.value("--max-recorded-rays") != nullptr
                                    ? integer<std::size_t>(*a.value("--max-recorded-rays"),
                                                           "--max-recorded-rays", 1, SIZE_MAX)
                                    : trace::kDefaultMaxRecordedRays;
      stats = trace::SequentialTracer().trace(cs, p.path, rays, paths, std::nullopt, limit);
      text = format == "csv" ? paths_csv(paths) : io::to_result_json(paths);
    } else {
      stats = trace::SequentialTracer().trace(cs, p.path, rays);
      text = format == "csv" ? trace_csv(rays) : io::to_result_json(rays);
    }
    stats_text = io::to_result_json(stats);
  });
  emit(text, a.value("--output"), out);
  if (const std::string* f = a.value("--stats")) emit(stats_text, f, out);
  return kOk;
}

int analyze_command(const std::vector<std::string>& args, std::ostream& out, std::ostream& err) {
  const std::vector<Kind>& kinds = analyze_kinds();
  if (args.size() < 2 || args[1] == "--help" || args[1] == "-h") {
    if (args.size() < 2) throw UsageError("rtt analyze needs a kind");
    out << "usage: rtt analyze <kind> <system.rtt.json> [options]\n\nkinds:\n";
    for (const Kind& k : kinds) {
      out << "  " << k.name << std::string(14 - k.name.size(), ' ') << k.help << "\n";
    }
    out << "\nrtt analyze <kind> --help lists the options of a kind.\n";
    return kOk;
  }
  const Kind* kind = nullptr;
  for (const Kind& k : kinds) {
    if (k.name == args[1]) kind = &k;
  }
  if (kind == nullptr) throw UsageError("unknown kind '" + args[1] + "' for rtt analyze");
  const Arguments a = parse(args, 2, with_common(kind->options));
  if (a.flag("--help")) {
    help(out, "rtt analyze " + std::string(kind->name) + " <system.rtt.json> [options]", kind->help,
         kind->options);
    return kOk;
  }
  if (a.operands.size() != 1) throw UsageError("rtt analyze needs one system file");
  const std::string format = format_of(a, kind->csv);
  Prepared p;
  prepare(p, a.operands[0], a, err);
  const compile::CompiledSystem& cs = *p.compiled;
  const compile::PathId path = p.path;
  std::optional<compile::CompiledSystem> nominal;
  const std::optional<analysis::NominalReference> ref = reference(p, a, err, nominal);
  const auto grid = [&](int fallback) {
    const std::string* g = a.value("--grid");
    return g != nullptr ? integer<int>(*g, "--grid", 1, 65535) : fallback;
  };
  const auto real = [&](std::string_view name, double fallback) {
    const std::string* v = a.value(name);
    return v != nullptr ? number(*v, name) : fallback;
  };
  const auto pupil = [&] {
    analysis::PupilFunctionOptions o;
    o.grid = grid(o.grid);
    o.apodization = real("--apodization", o.apodization);
    o.focus_shift = real("--focus-shift", o.focus_shift);
    o.aiming = aiming(a);
    return o;
  };
  const std::uint16_t wl = wavelength(p, a).value_or(cs.reference_wavelength());
  std::string text;
  in_arena(a, [&] {
    if (kind->name == "first-order") {
      text = io::to_result_json(paraxial::first_order(cs, path, wl));
    } else if (kind->name == "seidel") {
      std::optional<paraxial::ChromaticPair> pair;
      if (const std::string* v = a.value("--pair")) {
        const std::vector<std::string> ij = split(*v, ',');
        if (ij.size() != 2) throw UsageError("option --pair: '" + *v + "' is not i,j");
        const auto n = static_cast<std::uint16_t>(cs.wavelengths_um().size() - 1);
        pair = paraxial::ChromaticPair{integer<std::uint16_t>(ij[0], "--pair", 0, n),
                                       integer<std::uint16_t>(ij[1], "--pair", 0, n)};
      }
      text = io::to_result_json(paraxial::seidel(cs, path, wl, pair));
    } else if (kind->name == "spot") {
      analysis::SpotOptions o;
      if (const std::string* r = a.value("--rays")) o.sampling = sampling(*r);
      o.aiming = aiming(a);
      o.reference = ref;
      const analysis::SpotDiagram d = analysis::spot(cs, path, field(p, a), wavelength(p, a), o);
      text = format == "csv" ? spot_csv(d) : io::to_result_json(d);
    } else if (kind->name == "opd") {
      analysis::OpdOptions o;
      o.grid = grid(o.grid);
      o.focus_shift = real("--focus-shift", o.focus_shift);
      o.aiming = aiming(a);
      o.reference = ref;
      text = io::to_result_json(analysis::opd_map(cs, path, field(p, a), wl, o));
    } else if (kind->name == "zernike") {
      analysis::ZernikeOptions z;
      if (const std::string* t = a.value("--terms")) {
        z.terms = integer<int>(*t, "--terms", 1, 231);
      }
      text = io::to_result_json(
          analysis::zernike_fit(analysis::pupil_function(cs, path, field(p, a), wl, pupil()), z));
    } else if (kind->name == "psf") {
      analysis::PsfOptions o;
      o.pupil = pupil();
      if (const std::string* q = a.value("--padding")) {
        o.padding = integer<int>(*q, "--padding", 2, 65535);
      }
      o.wavelengths = wavelengths(p, a);
      text = io::to_result_json(analysis::psf(cs, path, field(p, a), o));
    } else {  // mtf
      analysis::MtfOptions o;
      o.pupil = pupil();
      o.wavelengths = wavelengths(p, a);
      o.outputs.map = a.flag("--map");
      o.outputs.complex_otf = a.flag("--complex-otf");
      const analysis::Mtf m = analysis::mtf(cs, path, field(p, a), o);
      text = format == "csv" ? mtf_csv(m) : io::to_result_json(m);
    }
  });
  emit(text, a.value("--output"), out);
  return kOk;
}

}  // namespace

int run_compute(const std::vector<std::string>& args, std::ostream& out, std::ostream& err) {
  // The system file for the messages: the first operand after the command (and the kind).
  std::string file;
  for (std::size_t i = args[0] == "analyze" ? 2 : 1; i < args.size(); ++i) {
    if (!args[i].starts_with("-")) {
      file = args[i];
      break;
    }
    if (args[i] != "--help" && args[i] != "--paths" && args[i] != "--polarization" &&
        args[i] != "--map" && args[i] != "--complex-otf") {
      ++i;  // the value of an option
    }
  }
  try {
    return args[0] == "trace" ? trace_command(args, out, err) : analyze_command(args, out, err);
  } catch (const UsageError& e) {
    err << "rtt " << args[0] << ": " << e.what() << " (see rtt " << args[0] << " --help)\n";
    return kUsage;
  } catch (const Failed&) {
    return kFailed;  // the diagnostics are printed
  } catch (const std::exception& e) {
    // ParaxialError, AnalysisError, NoStopError, std::invalid_argument of the API, files.
    err << (file.empty() ? std::string("rtt") : file) << ": error " << e.what() << "\n";
    return kFailed;
  }
}

}  // namespace rtt::cli
