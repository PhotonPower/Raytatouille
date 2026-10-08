#include "rtt/material/material.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <optional>
#include <stdexcept>
#include <system_error>
#include <utility>
#include <vector>

#include "rtt/material/agf.hpp"
#include "rtt/material/air.hpp"

namespace rtt::material {
namespace {

constexpr std::string_view kConstPrefix = "CONST:";

/// Medium with a fixed complex index n + i*kappa, independent of wavelength and temperature.
/// Sign convention n + i*kappa, kappa >= 0 for absorption: docs/architecture.md, "Konventionen".
class ConstantMaterial final : public Material {
 public:
  explicit ConstantMaterial(math::Complex index) : index_(index) {}

  [[nodiscard]] math::Complex index(double /*wavelength_um*/,
                                    double /*temperature_c*/,
                                    double /*pressure_atm*/) const override {
    return index_;
  }

 private:
  math::Complex index_;
};

/// Parses the whole of `text` as a finite double in C-locale syntax (no whitespace, no '+').
std::optional<double> parse_finite(std::string_view text) {
  double value = 0.0;
  const char* const first = text.data();
  const char* const last = first + text.size();
  const auto [ptr, ec] = std::from_chars(first, last, value);
  if (ec != std::errc{} || ptr != last || !std::isfinite(value)) return std::nullopt;
  return value;
}

[[noreturn]] void reject(std::string_view reference, std::string_view reason) {
  throw UnknownMaterial("invalid material reference '" + std::string(reference) +
                        "': " + std::string(reason));
}

/// Parses "<n>" or "<n>,<kappa>" (the part after "CONST:").
math::Complex parse_constant(std::string_view reference, std::string_view args) {
  const auto comma = args.find(',');
  const std::string_view n_text = args.substr(0, comma);
  const std::optional<double> n = parse_finite(n_text);
  if (!n) reject(reference, "expected CONST:<n> or CONST:<n>,<kappa> with finite numbers");
  if (*n <= 0.0) reject(reference, "refractive index n must be positive");
  if (comma == std::string_view::npos) return {*n, 0.0};

  const std::optional<double> kappa = parse_finite(args.substr(comma + 1));
  if (!kappa) reject(reference, "expected CONST:<n> or CONST:<n>,<kappa> with finite numbers");
  if (*kappa < 0.0) reject(reference, "extinction coefficient kappa must be >= 0");
  // Adding +0.0 maps -0.0 to +0.0 (IEEE 754-2019, 6.3): the sign of a zero imaginary part
  // selects the branch of complex square roots downstream, so kappa >= 0 includes the sign.
  return {*n, *kappa + 0.0};
}

std::shared_ptr<const Material> make_material(std::string_view reference) {
  if (reference == "VACUUM") return std::make_shared<const ConstantMaterial>(1.0);
  // Dry air after Ciddor at the temperature and pressure of the medium (#25).
  if (reference == "AIR") return std::make_shared<const AirMaterial>();
  if (reference.starts_with(kConstPrefix)) {
    return std::make_shared<const ConstantMaterial>(
        parse_constant(reference, reference.substr(kConstPrefix.size())));
  }
  throw UnknownMaterial(
      "unknown material reference '" + std::string(reference) +
      "': supported are VACUUM, AIR, CONST:<n>, CONST:<n>,<kappa>, catalogue glasses as "
      "KATALOG:NAME after MaterialLibrary::add_catalog and materials registered with "
      "MaterialLibrary::add");
}

}  // namespace

std::shared_ptr<const Material> MaterialLibrary::resolve(std::string_view reference) const {
  const std::scoped_lock lock(mutex_);
  if (const auto it = cache_.find(reference); it != cache_.end()) return it->second;
  if (const auto it = unsupported_.find(reference); it != unsupported_.end()) {
    throw UnknownMaterial(it->second);
  }
  const auto colon = reference.find(':');
  if (colon != std::string_view::npos && !reference.starts_with(kConstPrefix)) {
    const std::string_view catalog = reference.substr(0, colon);
    const std::string text = "unknown material reference '" + std::string(reference) + "': ";
    if (const auto cat = catalogs_.find(catalog); cat != catalogs_.end()) {
      throw UnknownMaterial(text + "glass '" + std::string(reference.substr(colon + 1)) +
                            "' not found in catalog " + std::string(catalog) + " (" +
                            cat->second->file + ")");
    }
    throw UnknownMaterial(text + "catalog " + std::string(catalog) +
                          " not loaded (use MaterialLibrary::add_catalog)");
  }
  auto material = make_material(reference);
  cache_.emplace(std::string(reference), material);
  return material;
}

void MaterialLibrary::add(std::string name, std::shared_ptr<const Material> material) {
  if (name.empty()) throw std::invalid_argument("MaterialLibrary::add: empty material name");
  if (name == "VACUUM" || name == "AIR" || name.starts_with(kConstPrefix)) {
    throw std::invalid_argument("MaterialLibrary::add: '" + name + "' is a reserved name");
  }
  if (material == nullptr) {
    throw std::invalid_argument("MaterialLibrary::add: material '" + name + "' is null");
  }
  const std::scoped_lock lock(mutex_);
  if (cache_.contains(name) || unsupported_.contains(name)) {
    throw std::invalid_argument("MaterialLibrary::add: '" + name + "' is already in use");
  }
  cache_.emplace(std::move(name), std::move(material));
}

namespace {

/// Checks an explicit catalogue name (alias or in-memory catalogue, #85).
void check_catalog_name(const std::string& name, const char* function) {
  const bool bad =
      name.empty() || name == "CONST" || name.find_first_of(": \t\n\r\v\f") != std::string::npos;
  if (bad) {
    throw std::invalid_argument(std::string("MaterialLibrary::") + function + ": catalog name '" +
                                name +
                                "' must not be empty, contain ':' or whitespace, or be CONST");
  }
}

}  // namespace

void MaterialLibrary::add_catalog(const std::filesystem::path& path,
                                  const std::optional<std::string>& name) {
  namespace fs = std::filesystem;
  if (name) check_catalog_name(*name, "add_catalog");
  // Files to load: the path itself, or the *.agf files of a directory in sorted order.
  std::vector<fs::path> files;
  if (fs::is_directory(path)) {
    if (name) {
      throw std::invalid_argument("MaterialLibrary::add_catalog: a catalog name (" + *name +
                                  ") needs a single .agf file, not the directory " + path.string());
    }
    for (const auto& entry : fs::directory_iterator(path)) {
      std::string ext = entry.path().extension().string();
      for (char& c : ext) c = static_cast<char>(c >= 'A' && c <= 'Z' ? c - 'A' + 'a' : c);
      if (entry.is_regular_file() && ext == ".agf") files.push_back(entry.path());
    }
    std::sort(files.begin(), files.end());
    if (files.empty()) {
      throw std::invalid_argument("MaterialLibrary::add_catalog: no .agf file in " + path.string());
    }
  } else {
    files.push_back(path);
  }
  std::vector<AgfCatalog> catalogs;
  for (const fs::path& file : files) {
    AgfCatalog cat = load_agf(file);
    if (name) cat.name = *name;
    if (cat.name.empty() || cat.name == "CONST") {
      throw std::invalid_argument("MaterialLibrary::add_catalog: catalog name '" + cat.name +
                                  "' of " + file.string() + " is empty or reserved");
    }
    catalogs.push_back(std::move(cat));
  }
  register_catalogs(std::move(catalogs), "add_catalog");
}

void MaterialLibrary::add_catalog_text(std::string_view bytes,
                                       const std::string& name,
                                       const std::string& source) {
  check_catalog_name(name, "add_catalog_text");
  std::vector<AgfCatalog> catalogs;
  catalogs.push_back(parse_agf(decode_agf_text(bytes, source), name, source));
  register_catalogs(std::move(catalogs), "add_catalog_text");
}

void MaterialLibrary::register_catalogs(std::vector<AgfCatalog> catalogs, const char* function) {
  const std::string prefix = std::string("MaterialLibrary::") + function + ": ";
  // Build the materials before touching the library (all or nothing).
  struct Entry {
    std::string reference;
    std::shared_ptr<const Material> material;  // null: unsupported, see message
    std::string message;
  };
  std::vector<Entry> entries;
  for (std::size_t i = 0; i < catalogs.size(); ++i) {
    const AgfCatalog& cat = catalogs[i];
    for (std::size_t j = 0; j < i; ++j) {
      if (catalogs[j].name == cat.name) {
        throw std::invalid_argument(prefix + "catalog " + cat.name + " appears twice (" +
                                    catalogs[j].file + ", " + cat.file + ")");
      }
    }
    // Lines of the blocks per name; a name with several blocks is ambiguous (R2, #71).
    std::map<std::string, std::vector<std::size_t>, std::less<>> lines;
    for (const AgfGlass& glass : cat.glasses) lines[glass.name].push_back(glass.line);
    for (const AgfGlass& glass : cat.glasses) {
      const std::vector<std::size_t>& at = lines[glass.name];
      if (at.size() > 1) {
        if (glass.line != at.front()) continue;  // one entry per name
        std::string where;
        for (std::size_t k = 0; k < at.size(); ++k) {
          where += (k == 0 ? "" : (k + 1 == at.size() ? " and " : ", ")) + std::to_string(at[k]);
        }
        const std::string reference = cat.name + ":" + glass.name;
        std::string message = "ambiguous material reference '";
        message += reference;
        message += "': glass ";
        message += glass.name;
        message += " is defined with different data at lines ";
        message += where;
        message += " of ";
        message += cat.file;
        message +=
            " (agf.duplicate_glass_conflict). To use one of the blocks, load it as a catalogue of "
            "its own under an alias, e.g. add_catalog_text(block, \"";
        message += cat.name;
        message += "_ALT\")";
        entries.push_back({reference, nullptr, std::move(message)});
        continue;
      }
      Entry entry{cat.name + ":" + glass.name, nullptr, {}};
      try {
        (void)agf_formula(glass, cat.file + ":" + std::to_string(glass.line));
        entry.material = std::make_shared<const CatalogMaterial>(glass, cat.name);
      } catch (const std::invalid_argument& e) {
        entry.message = "unsupported material reference '" + entry.reference + "': " + e.what();
      }
      entries.push_back(std::move(entry));
    }
  }

  const std::scoped_lock lock(mutex_);
  for (const AgfCatalog& cat : catalogs) {
    if (const auto it = catalogs_.find(cat.name); it != catalogs_.end()) {
      throw std::invalid_argument(prefix + "catalog " + cat.name + " is already loaded from " +
                                  it->second->file);
    }
  }
  for (const Entry& entry : entries) {
    if (cache_.contains(entry.reference) || unsupported_.contains(entry.reference)) {
      throw std::invalid_argument(prefix + "'" + entry.reference + "' is already in use");
    }
  }
  for (AgfCatalog& cat : catalogs) {
    // Warnings in loading order (#71); only after all checks, so a failed call adds none.
    load_warnings_.insert(load_warnings_.end(), cat.warnings.begin(), cat.warnings.end());
    std::string key = cat.name;
    catalogs_.emplace(std::move(key), std::make_shared<const AgfCatalog>(std::move(cat)));
  }
  for (Entry& entry : entries) {
    if (entry.material) {
      cache_.emplace(std::move(entry.reference), std::move(entry.material));
    } else {
      unsupported_.emplace(std::move(entry.reference), std::move(entry.message));
    }
  }
}

std::vector<std::string> MaterialLibrary::catalogs() const {
  const std::scoped_lock lock(mutex_);
  std::vector<std::string> names;
  names.reserve(catalogs_.size());
  for (const auto& [name, cat] : catalogs_) names.push_back(name);
  return names;  // std::map keeps them in ascending order
}

std::vector<LoadWarning> MaterialLibrary::load_warnings() const {
  const std::scoped_lock lock(mutex_);
  return load_warnings_;
}

std::shared_ptr<const AgfCatalog> MaterialLibrary::catalog(std::string_view name) const {
  const std::scoped_lock lock(mutex_);
  if (const auto it = catalogs_.find(name); it != catalogs_.end()) return it->second;
  throw UnknownMaterial("catalog " + std::string(name) +
                        " not loaded (use MaterialLibrary::add_catalog)");
}

std::vector<math::Complex> index_many(const Material& material,
                                      std::span<const double> wavelength_um,
                                      double temperature_c,
                                      double pressure_atm) {
  std::vector<math::Complex> out;
  out.reserve(wavelength_um.size());
  for (const double wl : wavelength_um) {
    out.push_back(material.index(wl, temperature_c, pressure_atm));
  }
  return out;
}

}  // namespace rtt::material
