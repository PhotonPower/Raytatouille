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
  // M1 placeholder: air as n = 1. Ciddor air and glass data relative to air follow in M2.
  if (reference == "AIR") return std::make_shared<const ConstantMaterial>(1.0);
  if (reference.starts_with(kConstPrefix)) {
    return std::make_shared<const ConstantMaterial>(
        parse_constant(reference, reference.substr(kConstPrefix.size())));
  }
  throw UnknownMaterial(
      "unknown material reference '" + std::string(reference) +
      "': supported are VACUUM, AIR, CONST:<n>, CONST:<n>,<kappa> and "
      "materials registered with MaterialLibrary::add; glass catalogs follow in M2");
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
                            "' not found in catalog " + std::string(catalog) + " (" + cat->second +
                            ")");
    }
    throw UnknownMaterial(text + "catalog " + std::string(catalog) +
                          " is not loaded (glass catalogs, M2: MaterialLibrary::add_catalog)");
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

void MaterialLibrary::add_catalog(const std::filesystem::path& path) {
  namespace fs = std::filesystem;
  // Files to load: the path itself, or the *.agf files of a directory in sorted order.
  std::vector<fs::path> files;
  if (fs::is_directory(path)) {
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

  // Parse everything and build the materials before touching the library (all or nothing).
  struct Entry {
    std::string reference;
    std::shared_ptr<const Material> material;  // null: unsupported, see message
    std::string message;
  };
  std::vector<AgfCatalog> catalogs;
  std::vector<Entry> entries;
  for (const fs::path& file : files) {
    AgfCatalog cat = load_agf(file);
    if (cat.name.empty() || cat.name == "CONST") {
      throw std::invalid_argument("MaterialLibrary::add_catalog: catalog name '" + cat.name +
                                  "' of " + file.string() + " is empty or reserved");
    }
    for (const AgfCatalog& other : catalogs) {
      if (other.name == cat.name) {
        throw std::invalid_argument("MaterialLibrary::add_catalog: catalog " + cat.name +
                                    " appears twice (" + other.file + ", " + cat.file + ")");
      }
    }
    for (const AgfGlass& glass : cat.glasses) {
      Entry entry{cat.name + ":" + glass.name, nullptr, {}};
      try {
        (void)agf_formula(glass, cat.file + ":" + std::to_string(glass.line));
        entry.material = std::make_shared<const CatalogMaterial>(glass, cat.name);
      } catch (const std::invalid_argument& e) {
        entry.message = "unknown material reference '" + entry.reference + "': " + e.what();
      }
      entries.push_back(std::move(entry));
    }
    catalogs.push_back(std::move(cat));
  }

  const std::scoped_lock lock(mutex_);
  for (const AgfCatalog& cat : catalogs) {
    if (const auto it = catalogs_.find(cat.name); it != catalogs_.end()) {
      throw std::invalid_argument("MaterialLibrary::add_catalog: catalog " + cat.name +
                                  " is already loaded from " + it->second);
    }
  }
  for (const Entry& entry : entries) {
    if (cache_.contains(entry.reference) || unsupported_.contains(entry.reference)) {
      throw std::invalid_argument("MaterialLibrary::add_catalog: '" + entry.reference +
                                  "' is already in use");
    }
  }
  for (const AgfCatalog& cat : catalogs) catalogs_.emplace(cat.name, cat.file);
  for (Entry& entry : entries) {
    if (entry.material) {
      cache_.emplace(std::move(entry.reference), std::move(entry.material));
    } else {
      unsupported_.emplace(std::move(entry.reference), std::move(entry.message));
    }
  }
}

}  // namespace rtt::material
