#include "rtt/io/lens_catalog.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <initializer_list>
#include <nlohmann/json.hpp>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "rtt/json/strict.hpp"

namespace rtt::io {
namespace {

using Json = nlohmann::json;

/// Counts of radii, thicknesses, gaps and materials per kind (ADR 0053, point 2).
struct Counts {
  std::size_t radii;
  std::size_t thicknesses;
  std::size_t gaps;
  std::size_t materials;
};

constexpr Counts counts(LensKind kind) noexcept {
  switch (kind) {
    case LensKind::Singlet:
      return {2, 1, 0, 1};
    case LensKind::Doublet:
      return {3, 2, 0, 2};
    case LensKind::Triplet:
      return {6, 3, 2, 3};
  }
  return {0, 0, 0, 0};
}

std::string at(const std::string& base, std::size_t i) {
  return base + "/" + std::to_string(i);
}

/// Strict reader of one catalogue (ADR 0008): every error names the file and a JSON pointer.
class Reader {
 public:
  explicit Reader(std::string file) : file_(std::move(file)) {}

  [[noreturn]] void fail(const std::string& pointer, const std::string& message) const {
    throw LensCatalogError(file_, pointer, message);
  }

  /// Object with only the keys in `allowed`.
  void object(const Json& j,
              const std::string& pointer,
              std::initializer_list<std::string_view> allowed) const {
    if (!j.is_object()) fail(pointer, "expected an object");
    for (const auto& item : j.items()) {
      if (std::find(allowed.begin(), allowed.end(), item.key()) == allowed.end()) {
        fail(pointer + "/" + rtt::json::pointer_token(item.key()),
             "unknown key '" + item.key() + "'");
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

  /// Name of a catalogue or lens: not empty, no ':'.
  [[nodiscard]] std::string name(const Json& j, const std::string& pointer) const {
    std::string text = string(j, pointer);
    if (text.empty()) fail(pointer, "empty name");
    if (text.find(':') != std::string::npos) fail(pointer, "'" + text + "' contains ':'");
    return text;
  }

  /// Finite number, > 0 if `positive`, else >= 0.
  [[nodiscard]] double length(const Json& j, const std::string& pointer, bool positive) const {
    if (!j.is_number()) fail(pointer, "expected a number");
    const double value = j.get<double>();
    if (!std::isfinite(value) || value < 0.0 || (positive && value == 0.0)) {
      fail(pointer, positive ? "expected a finite number > 0" : "expected a finite number >= 0");
    }
    return value;
  }

  /// Array of exactly `count` entries.
  const Json& array(const Json& j, const std::string& pointer, std::size_t count) const {
    if (!j.is_array()) fail(pointer, "expected an array");
    if (j.size() != count) {
      fail(pointer,
           "expected " + std::to_string(count) + " entries, got " + std::to_string(j.size()));
    }
    return j;
  }

 private:
  std::string file_;
};

LensKind read_kind(const Reader& r, const Json& j, const std::string& pointer) {
  const std::string text = r.string(j, pointer);
  for (const LensKind kind : {LensKind::Singlet, LensKind::Doublet, LensKind::Triplet}) {
    if (text == lens_kind_name(kind)) return kind;
  }
  r.fail(pointer, "unknown kind '" + text + "' (singlet, doublet or triplet)");
}

LensEntry read_lens(const Reader& r, const Json& j, const std::string& pointer) {
  r.object(
      j, pointer,
      {"name", "kind", "description", "radii", "thicknesses", "gaps", "materials", "diameter"});
  LensEntry lens;
  lens.name = r.name(r.required(j, pointer, "name"), pointer + "/name");
  lens.kind = read_kind(r, r.required(j, pointer, "kind"), pointer + "/kind");
  const Counts n = counts(lens.kind);
  if (j.contains("description")) {
    lens.description = r.string(j.at("description"), pointer + "/description");
  }
  const std::string radii_pointer = pointer + "/radii";
  const Json& radii = r.array(r.required(j, pointer, "radii"), radii_pointer, n.radii);
  for (std::size_t i = 0; i < radii.size(); ++i) {
    if (radii[i].is_null()) {
      lens.radii.emplace_back(std::nullopt);
      continue;
    }
    if (!radii[i].is_number()) r.fail(at(radii_pointer, i), "expected a number or null");
    const double radius = radii[i].get<double>();
    if (!std::isfinite(radius) || radius == 0.0) {
      r.fail(at(radii_pointer, i), "expected a finite radius != 0 (null for plano)");
    }
    lens.radii.emplace_back(radius);
  }
  const std::string t_pointer = pointer + "/thicknesses";
  const Json& thicknesses =
      r.array(r.required(j, pointer, "thicknesses"), t_pointer, n.thicknesses);
  for (std::size_t i = 0; i < thicknesses.size(); ++i) {
    lens.thicknesses.push_back(r.length(thicknesses[i], at(t_pointer, i), true));
  }
  if (n.gaps == 0) {
    if (j.contains("gaps")) r.fail(pointer + "/gaps", "only a triplet has gaps");
  } else {
    const std::string g_pointer = pointer + "/gaps";
    const Json& gaps = r.array(r.required(j, pointer, "gaps"), g_pointer, n.gaps);
    for (std::size_t i = 0; i < gaps.size(); ++i) {
      lens.gaps.push_back(r.length(gaps[i], at(g_pointer, i), false));
    }
  }
  const std::string m_pointer = pointer + "/materials";
  const Json& materials = r.array(r.required(j, pointer, "materials"), m_pointer, n.materials);
  for (std::size_t i = 0; i < materials.size(); ++i) {
    std::string material = r.string(materials[i], at(m_pointer, i));
    if (material.empty()) r.fail(at(m_pointer, i), "empty material reference");
    lens.materials.push_back(std::move(material));
  }
  if (j.contains("diameter")) {
    lens.diameter = r.length(j.at("diameter"), pointer + "/diameter", true);
  }
  return lens;
}

}  // namespace

LensCatalogError::LensCatalogError(std::string file,
                                   std::string pointer,
                                   const std::string& message)
    : std::runtime_error((file.empty() ? std::string("lens catalog") : file) + ": " +
                         (pointer.empty() ? std::string("/") : pointer) + ": " + message),
      file_(std::move(file)),
      pointer_(std::move(pointer)) {}

std::string_view lens_kind_name(LensKind kind) noexcept {
  switch (kind) {
    case LensKind::Singlet:
      return "singlet";
    case LensKind::Doublet:
      return "doublet";
    case LensKind::Triplet:
      return "triplet";
  }
  return "";
}

LensCatalog parse_lens_catalog(std::string_view json_text, const std::string& file) {
  const Reader r(file);
  Json j;
  try {
    // Duplicate keys, syntax errors and number overflow (ADR 0008 addendum, ADR 0020).
    j = rtt::json::parse_strict(json_text);
  } catch (const rtt::json::StrictParseError& e) {
    r.fail(e.pointer(), e.message());
  }
  r.object(j, "", {"format", "schema_version", "catalog", "lenses"});
  if (r.string(r.required(j, "", "format"), "/format") != kLensCatalogFormat) {
    r.fail("/format", "expected '" + std::string(kLensCatalogFormat) + "'");
  }
  if (r.string(r.required(j, "", "schema_version"), "/schema_version") !=
      kLensCatalogSchemaVersion) {
    r.fail("/schema_version",
           "unsupported version, this library reads " + std::string(kLensCatalogSchemaVersion));
  }
  LensCatalog catalog;
  catalog.name = r.name(r.required(j, "", "catalog"), "/catalog");
  const Json& lenses = r.required(j, "", "lenses");
  if (!lenses.is_array()) r.fail("/lenses", "expected an array");
  for (std::size_t i = 0; i < lenses.size(); ++i) {
    LensEntry lens = read_lens(r, lenses[i], at("/lenses", i));
    for (const LensEntry& other : catalog.lenses) {
      if (other.name == lens.name) {
        r.fail(at("/lenses", i) + "/name", "duplicate lens name '" + lens.name + "'");
      }
    }
    catalog.lenses.push_back(std::move(lens));
  }
  return catalog;
}

LensCatalog load_lens_catalog(const std::filesystem::path& file) {
  const std::ifstream in(file, std::ios::binary);
  if (!in) throw LensCatalogError(file.string(), "", "cannot open the file");
  std::ostringstream text;
  text << in.rdbuf();
  return parse_lens_catalog(text.str(), file.string());
}

}  // namespace rtt::io
