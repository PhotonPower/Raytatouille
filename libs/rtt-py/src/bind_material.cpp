// Material library bindings (#85): catalogue alias, loading from memory, listing of AGF
// catalogues and a vectorized index. The listing is raw data (_glass_records); the Python types
// (GlassInfo, ClassRange, glass map) are in raytatouille/materials.py.

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/filesystem.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <cstddef>
#include <memory>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "bindings.hpp"
#include "rtt/material/agf.hpp"
#include "rtt/material/material.hpp"

namespace nb = nanobind;
using namespace nb::literals;

namespace rtt::py {
namespace {

using material::AgfClassRange;
using material::AgfGlass;
using material::MaterialLibrary;

using WavelengthsIn = nb::ndarray<const double, nb::c_contig, nb::device::cpu>;

/// NumPy array of `shape` that owns `values`.
template <typename T>
nb::ndarray<nb::numpy, T> owned(std::vector<T>&& values, const std::vector<std::size_t>& shape) {
  auto data = std::make_unique<std::vector<T>>(std::move(values));
  T* pointer = data->data();
  // The capsule owns the copy from here on and frees it with the array.
  nb::capsule owner(data.get(), [](void* p) noexcept {
    const std::unique_ptr<std::vector<T>> owned_data(static_cast<std::vector<T>*>(p));
  });
  static_cast<void>(data.release());
  return nb::ndarray<nb::numpy, T>(pointer, shape.size(), shape.data(), owner);
}

nb::object optional_value(const std::optional<double>& v) {
  return v ? nb::cast(*v) : nb::none();
}

nb::object optional_list(const std::optional<std::vector<std::optional<double>>>& values) {
  if (!values) return nb::none();
  nb::list list;
  for (const auto& v : *values) list.append(optional_value(v));
  return list;
}

nb::object class_range(const std::optional<AgfClassRange>& r) {
  return r ? nb::cast(std::pair{r->low, r->high}) : nb::none();
}

/// Raw data of one glass as written in the AGF file (see AgfGlass); the meaning of the values is
/// applied in raytatouille/materials.py.
nb::dict glass_record(const MaterialLibrary& lib, const std::string& catalog, const AgfGlass& g) {
  nb::dict d;
  d["catalog"] = catalog;
  d["name"] = g.name;
  d["line"] = g.line;
  d["formula"] = g.formula;
  // One source: a glass is supported exactly if the library resolves it (verified formula, LD
  // range, TD with 7 values); otherwise the library's message is the reason.
  try {
    (void)lib.resolve(catalog + ":" + g.name);
    d["unsupported_reason"] = nb::none();
  } catch (const material::UnknownMaterial& e) {
    d["unsupported_reason"] = std::string(e.what());
  }
  d["coefficients"] = g.coefficients;
  d["nd"] = g.nd;
  d["vd"] = g.vd;
  d["exclude_substitution"] = g.exclude_substitution;
  d["status"] = g.status;
  d["melt_frequency"] = g.melt_frequency;
  d["comment"] = g.comment;
  d["wavelength_range_um"] =
      g.range ? nb::cast(std::pair{g.range->min_um, g.range->max_um}) : nb::none();
  d["thermal"] = g.thermal;
  d["extra"] = g.extra;
  d["mechanical"] = optional_list(g.mechanical);
  if (g.other) {
    nb::dict od;
    od["relative_cost"] = optional_value(g.other->relative_cost);
    od["cr"] = class_range(g.other->cr);
    od["fr"] = class_range(g.other->fr);
    od["sr"] = class_range(g.other->sr);
    od["ar"] = class_range(g.other->ar);
    od["pr"] = class_range(g.other->pr);
    d["other"] = od;
  } else {
    d["other"] = nb::none();
  }
  std::vector<double> it;
  it.reserve(3 * g.transmission.size());
  for (const auto& t : g.transmission) {
    it.insert(it.end(), {t.wavelength_um, t.transmittance, t.thickness_mm});
  }
  d["transmission"] = owned(std::move(it), {g.transmission.size(), 3});
  return d;
}

}  // namespace

void bind_material(nb::module_& m) {
  nb::class_<MaterialLibrary>(
      m, "MaterialLibrary",
      "Resolves material references (\"VACUUM\", \"AIR\", \"CONST:<n>\", catalogue glasses).")
      .def(nb::init<>())
      .def(
          "add_catalog",
          [](MaterialLibrary& lib, const std::filesystem::path& path,
             const std::optional<std::string>& name) { lib.add_catalog(path, name); },
          "path"_a, "name"_a = nb::none(),
          "Loads AGF glass catalogues: an .agf file, or a directory whose *.agf files are "
          "loaded in sorted order. The catalogue name is the file name without extension in "
          "upper case (schott.agf -> SCHOTT; glasses resolve as SCHOTT:N-BK7), or `name` for a "
          "single file (alias, e.g. \"SCHOTT_M2\" for a second schott.agf; case-sensitive, not "
          "empty, no ':' or whitespace, not CONST).\n\nRaises AgfError for malformed files and "
          "ValueError for invalid or already used names.")
      .def(
          "add_catalog_text",
          [](MaterialLibrary& lib, const nb::bytes& data, const std::string& name,
             const std::string& source) {
            lib.add_catalog_text(std::string_view(data.c_str(), data.size()), name, source);
          },
          "data"_a, "name"_a, "source"_a = "<memory>",
          "Loads an AGF catalogue from memory under `name` (same rules as the alias of "
          "add_catalog). `data` is the file content (UTF-16LE with byte order mark, UTF-8 or "
          "ANSI); `source` names it in error messages.\n\nRaises AgfError and ValueError like "
          "add_catalog.")
      .def(
          "add_catalog_text",
          [](MaterialLibrary& lib, const std::string& text, const std::string& name,
             const std::string& source) { lib.add_catalog_text(text, name, source); },
          "data"_a, "name"_a, "source"_a = "<memory>",
          "As above, with the catalogue as text (str, encoded as UTF-8).")
      .def("catalogs", &MaterialLibrary::catalogs,
           "Names of the loaded catalogues in ascending order.")
      .def(
          "_glass_records",
          [](const MaterialLibrary& lib, const std::string& catalog) {
            const auto cat = lib.catalog(catalog);
            nb::list records;
            for (const AgfGlass& g : cat->glasses) records.append(glass_record(lib, cat->name, g));
            return records;
          },
          "catalog"_a, "Raw AGF data of the glasses of `catalog` (use MaterialLibrary.glasses).")
      .def(
          "index",
          [](const MaterialLibrary& lib, const std::string& reference, double wavelength_um,
             double temperature_c, double pressure_atm) {
            return lib.resolve(reference)->index(wavelength_um, temperature_c, pressure_atm);
          },
          "reference"_a, "wavelength_um"_a, "temperature_c"_a = 20.0, "pressure_atm"_a = 1.0,
          "Absolute complex index n + i*kappa (kappa >= 0 absorbs) of the material `reference` "
          "at the vacuum wavelength `wavelength_um` in um, temperature in degC and pressure in "
          "atm.\n\nRaises UnknownMaterial if the reference cannot be resolved.")
      .def(
          "index",
          [](const MaterialLibrary& lib, const std::string& reference,
             const WavelengthsIn& wavelength_um, double temperature_c, double pressure_atm) {
            const auto material = lib.resolve(reference);
            const std::span<const double> wl(wavelength_um.data(), wavelength_um.size());
            std::vector<std::size_t> shape;
            for (std::size_t i = 0; i < wavelength_um.ndim(); ++i) {
              shape.push_back(wavelength_um.shape(i));
            }
            return owned(material::index_many(*material, wl, temperature_c, pressure_atm), shape);
          },
          "reference"_a, "wavelength_um"_a, "temperature_c"_a = 20.0, "pressure_atm"_a = 1.0,
          "Vectorized form for a NumPy array of vacuum wavelengths in um: complex128 array of the "
          "same shape, element by element identical (bit for bit) with the scalar call.\n\nRaises "
          "UnknownMaterial if the reference cannot be resolved.");
}

}  // namespace rtt::py
