#include "rtt/coating/catalog.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <initializer_list>
#include <iterator>
#include <memory>
#include <mutex>
#include <nlohmann/json.hpp>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

namespace rtt::coating {
namespace {

using Json = nlohmann::json;

/// Strict reader of one catalogue (ADR 0008): every error names the file and a JSON pointer.
class Reader {
 public:
  explicit Reader(std::string file) : file_(std::move(file)) {}

  [[noreturn]] void fail(const std::string& pointer, const std::string& message) const {
    throw CoatingCatalogError(file_, pointer, message);
  }

  /// Object with only the keys in `allowed`.
  void object(const Json& j,
              const std::string& pointer,
              std::initializer_list<std::string_view> allowed) const {
    if (!j.is_object()) fail(pointer, "expected an object");
    for (const auto& item : j.items()) {
      if (std::find(allowed.begin(), allowed.end(), item.key()) == allowed.end()) {
        fail(pointer + "/" + item.key(), "unknown key '" + item.key() + "'");
      }
    }
  }

  [[nodiscard]] const Json& required(const Json& j,
                                     const std::string& pointer,
                                     const std::string& key) const {
    const auto it = j.find(key);
    if (it == j.end()) fail(pointer, "missing key '" + key + "'");
    return *it;
  }

  [[nodiscard]] std::string string(const Json& j, const std::string& pointer) const {
    if (!j.is_string()) fail(pointer, "expected a string");
    return j.get<std::string>();
  }

  /// Name of a catalogue or coating: not empty, no ':'.
  [[nodiscard]] std::string name(const Json& j, const std::string& pointer) const {
    std::string text = string(j, pointer);
    if (text.empty()) fail(pointer, "empty name");
    if (text.find(':') != std::string::npos) {
      fail(pointer, "'" + text + "' contains ':' (references are CATALOG:NAME)");
    }
    return text;
  }

  /// Finite number >= 0, or > 0 if `positive`.
  [[nodiscard]] double number(const Json& j, const std::string& pointer, bool positive) const {
    if (!j.is_number()) fail(pointer, "expected a number");
    const double value = j.get<double>();
    if (!std::isfinite(value) || value < 0.0 || (positive && value == 0.0)) {
      fail(pointer, positive ? "expected a finite number > 0" : "expected a finite number >= 0");
    }
    return value;
  }

 private:
  std::string file_;
};

std::string at(const std::string& base, std::size_t i) {
  return base + "/" + std::to_string(i);
}

LayerSpec read_layer(const Reader& r, const Json& j, const std::string& pointer) {
  r.object(j, pointer, {"material", "thickness_um", "qwot"});
  LayerSpec layer;
  layer.material = r.string(r.required(j, pointer, "material"), pointer + "/material");
  if (layer.material.empty()) r.fail(pointer + "/material", "empty material reference");
  const bool physical = j.contains("thickness_um");
  const bool qwot = j.contains("qwot");
  if (physical == qwot) {
    r.fail(pointer, "a layer needs exactly one of 'thickness_um' and 'qwot'");
  }
  if (physical) {
    layer.thickness =
        PhysicalThickness{r.number(j.at("thickness_um"), pointer + "/thickness_um", false)};
  } else {
    // The design wavelength is filled in by read_coating().
    layer.thickness = QuarterWaves{r.number(j.at("qwot"), pointer + "/qwot", false), 0.0};
  }
  return layer;
}

CoatingDesign read_coating(const Reader& r, const Json& j, const std::string& pointer) {
  r.object(j, pointer, {"name", "description", "design_wavelength_um", "layers"});
  CoatingDesign design;
  design.name = r.name(r.required(j, pointer, "name"), pointer + "/name");
  if (j.contains("description")) {
    design.description = r.string(j.at("description"), pointer + "/description");
  }
  if (j.contains("design_wavelength_um")) {
    design.design_wavelength_um =
        r.number(j.at("design_wavelength_um"), pointer + "/design_wavelength_um", true);
  }
  const Json& layers = r.required(j, pointer, "layers");
  const std::string layers_pointer = pointer + "/layers";
  if (!layers.is_array()) r.fail(layers_pointer, "expected an array");
  if (layers.empty()) r.fail(layers_pointer, "a coating needs at least one layer");
  for (std::size_t i = 0; i < layers.size(); ++i) {
    LayerSpec layer = read_layer(r, layers[i], at(layers_pointer, i));
    if (auto* q = std::get_if<QuarterWaves>(&layer.thickness)) {
      if (!design.design_wavelength_um) {
        r.fail(at(layers_pointer, i) + "/qwot",
               "a QWOT layer needs 'design_wavelength_um' in its coating");
      }
      q->design_wavelength_um = *design.design_wavelength_um;
    }
    design.layers.push_back(std::move(layer));
  }
  return design;
}

}  // namespace

CoatingCatalogError::CoatingCatalogError(std::string file,
                                         std::string pointer,
                                         const std::string& message)
    : std::runtime_error((file.empty() ? std::string("coating catalog") : file) + ": " +
                         (pointer.empty() ? std::string("/") : pointer) + ": " + message),
      file_(std::move(file)),
      pointer_(std::move(pointer)) {}

CoatingCatalog parse_coating_catalog(std::string_view json_text, const std::string& file) {
  const Reader r(file);
  Json j;
  try {
    j = Json::parse(json_text);
  } catch (const Json::parse_error& e) {
    r.fail("", std::string("JSON syntax error: ") + e.what());
  }
  r.object(j, "", {"format", "schema_version", "catalog", "coatings"});
  if (r.string(r.required(j, "", "format"), "/format") != kCatalogFormat) {
    r.fail("/format", "expected '" + std::string(kCatalogFormat) + "'");
  }
  if (r.string(r.required(j, "", "schema_version"), "/schema_version") != kCatalogSchemaVersion) {
    r.fail("/schema_version",
           "unsupported version, this library reads " + std::string(kCatalogSchemaVersion));
  }
  CoatingCatalog catalog;
  catalog.name = r.name(r.required(j, "", "catalog"), "/catalog");
  const Json& coatings = r.required(j, "", "coatings");
  if (!coatings.is_array()) r.fail("/coatings", "expected an array");
  for (std::size_t i = 0; i < coatings.size(); ++i) {
    CoatingDesign design = read_coating(r, coatings[i], at("/coatings", i));
    for (const CoatingDesign& other : catalog.coatings) {
      if (other.name == design.name) {
        r.fail(at("/coatings", i) + "/name", "duplicate coating name '" + design.name + "'");
      }
    }
    catalog.coatings.push_back(std::move(design));
  }
  return catalog;
}

CoatingCatalog load_coating_catalog(const std::filesystem::path& file) {
  const std::ifstream in(file, std::ios::binary);
  if (!in) throw CoatingCatalogError(file.string(), "", "cannot open the file");
  std::ostringstream text;
  text << in.rdbuf();
  return parse_coating_catalog(text.str(), file.string());
}

void CoatingLibrary::add(const CoatingCatalog& catalog) {
  if (catalog.name.empty() || catalog.name.find(':') != std::string::npos) {
    throw std::invalid_argument("CoatingLibrary::add: invalid catalog name '" + catalog.name + "'");
  }
  const std::scoped_lock lock(mutex_);
  if (catalogs_.contains(catalog.name)) {
    throw std::invalid_argument("CoatingLibrary::add: catalog " + catalog.name +
                                " is already loaded");
  }
  catalogs_.insert(catalog.name);
  for (const CoatingDesign& design : catalog.coatings) {
    designs_.emplace(catalog.name + ":" + design.name,
                     std::make_shared<const CoatingDesign>(design));
  }
}

void CoatingLibrary::add_catalog(const std::filesystem::path& path) {
  namespace fs = std::filesystem;
  std::vector<fs::path> files;
  if (fs::is_directory(path)) {
    for (const auto& entry : fs::directory_iterator(path)) {
      std::string ext = entry.path().extension().string();
      for (char& c : ext) c = static_cast<char>(c >= 'A' && c <= 'Z' ? c - 'A' + 'a' : c);
      if (entry.is_regular_file() && ext == ".json") files.push_back(entry.path());
    }
    std::sort(files.begin(), files.end());
    if (files.empty()) {
      throw std::invalid_argument("CoatingLibrary::add_catalog: no .json file in " + path.string());
    }
  } else {
    files.push_back(path);
  }
  // Parse everything first (all or nothing).
  std::vector<CoatingCatalog> catalogs;
  for (const fs::path& file : files) {
    CoatingCatalog catalog = load_coating_catalog(file);
    for (const CoatingCatalog& other : catalogs) {
      if (other.name == catalog.name) {
        throw std::invalid_argument("CoatingLibrary::add_catalog: catalog " + catalog.name +
                                    " is defined twice (" + file.string() + ")");
      }
    }
    catalogs.push_back(std::move(catalog));
  }
  const std::scoped_lock lock(mutex_);
  for (const CoatingCatalog& catalog : catalogs) {
    if (catalogs_.contains(catalog.name)) {
      throw std::invalid_argument("CoatingLibrary::add_catalog: catalog " + catalog.name +
                                  " is already loaded");
    }
  }
  for (const CoatingCatalog& catalog : catalogs) {
    catalogs_.insert(catalog.name);
    for (const CoatingDesign& design : catalog.coatings) {
      designs_.emplace(catalog.name + ":" + design.name,
                       std::make_shared<const CoatingDesign>(design));
    }
  }
}

std::shared_ptr<const CoatingDesign> CoatingLibrary::find(std::string_view reference) const {
  const std::scoped_lock lock(mutex_);
  const auto it = designs_.find(reference);
  return it == designs_.end() ? nullptr : it->second;
}

std::shared_ptr<const CoatingDesign> CoatingLibrary::resolve(std::string_view reference) const {
  if (auto design = find(reference)) return design;
  const std::string text = "unknown coating reference '" + std::string(reference) + "': ";
  const auto colon = reference.find(':');
  if (colon == std::string_view::npos) {
    throw UnknownCoating(text + "expected CATALOG:NAME");
  }
  const std::string_view catalog = reference.substr(0, colon);
  const std::scoped_lock lock(mutex_);
  if (catalogs_.contains(catalog)) {
    throw UnknownCoating(text + "coating '" + std::string(reference.substr(colon + 1)) +
                         "' not found in catalog " + std::string(catalog));
  }
  throw UnknownCoating(text + "catalog " + std::string(catalog) +
                       " not loaded (use CoatingLibrary::add_catalog)");
}

}  // namespace rtt::coating
