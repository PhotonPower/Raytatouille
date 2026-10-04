#include <nanobind/nanobind.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/filesystem.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <variant>
#include <vector>

#include "bindings.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/system.hpp"

namespace nb = nanobind;
using namespace nb::literals;

namespace rtt::py {

compile::PathId path_id(const compile::CompiledSystem& system, const PathArg& path) {
  return std::visit(
      [&](const auto& p) -> compile::PathId {
        if constexpr (std::is_same_v<std::decay_t<decltype(p)>, std::string>) {
          const auto id = system.find_path(p);
          if (!id) throw std::invalid_argument("unknown path '" + p + "'");
          return *id;
        } else {
          return compile::PathId{p};
        }
      },
      path);
}

std::uint16_t wavelength_index(const compile::CompiledSystem& system,
                               std::optional<std::uint16_t> wavelength) {
  return wavelength.value_or(system.reference_wavelength());
}

void bind_compile(nb::module_& m) {
  nb::class_<material::MaterialLibrary>(
      m, "MaterialLibrary",
      "Resolves material references (\"VACUUM\", \"AIR\", \"CONST:<n>\", catalogue glasses).")
      .def(nb::init<>())
      .def("add_catalog", &material::MaterialLibrary::add_catalog, "path"_a,
           "Loads AGF glass catalogues: an .agf file, or a directory whose *.agf files are "
           "loaded in sorted order. The catalogue name is the file name without extension in "
           "upper case (schott.agf -> SCHOTT; glasses resolve as SCHOTT:N-BK7).\n\nRaises "
           "AgfError for malformed files and ValueError for names already in use.")
      .def(
          "index",
          [](const material::MaterialLibrary& lib, const std::string& reference,
             double wavelength_um, double temperature_c, double pressure_atm) {
            return lib.resolve(reference)->index(wavelength_um, temperature_c, pressure_atm);
          },
          "reference"_a, "wavelength_um"_a, "temperature_c"_a = 20.0, "pressure_atm"_a = 1.0,
          "Absolute complex index n + i*kappa (kappa >= 0 absorbs) of the material `reference` "
          "at the vacuum wavelength `wavelength_um` in um, temperature in degC and pressure in "
          "atm.\n\nRaises UnknownMaterial if the reference cannot be resolved.");

  nb::class_<compile::CompiledMedium>(m, "CompiledMedium", "Medium evaluated at all wavelengths.")
      .def_ro("reference", &compile::CompiledMedium::reference,
              "Material reference as in the model, e.g. \"AIR\".")
      .def_ro("index", &compile::CompiledMedium::index,
              "Complex index n + i*kappa per system wavelength (kappa >= 0 absorbs).");

  nb::class_<compile::CompiledSystem>(
      m, "CompiledSystem",
      "Immutable compiled form of a System on which paraxial data and tracing work.")
      .def_prop_ro("wavelengths_um", &compile::CompiledSystem::wavelengths_um,
                   "System wavelengths in um (vacuum), in model order.")
      .def_prop_ro("reference_wavelength", &compile::CompiledSystem::reference_wavelength,
                   "Index of the reference wavelength.")
      .def_prop_ro("wavelength_weights", &compile::CompiledSystem::wavelength_weights,
                   "Wavelength weights as in the model (not normalised).")
      .def_prop_ro("temperature_c", &compile::CompiledSystem::temperature_c,
                   "Environment temperature in degC at which the media were evaluated.")
      .def_prop_ro(
          "field_count", [](const compile::CompiledSystem& s) { return s.fields().points.size(); },
          "Number of field points.")
      .def_prop_ro("media", &compile::CompiledSystem::media,
                   "Distinct media; media[0] is the environment.")
      .def_prop_ro(
          "path_names",
          [](const compile::CompiledSystem& s) {
            std::vector<std::string> names;
            for (const auto& p : s.paths()) names.push_back(p.name);
            return names;
          },
          "Path names in model order; a path is addressed by name or by this index.")
      .def_prop_ro(
          "surface_ids",
          [](const compile::CompiledSystem& s) {
            std::vector<std::string> ids;
            for (const auto& surface : s.surfaces()) ids.push_back(surface.id.str());
            return ids;
          },
          "Surface ids in tree order (index = surface index, e.g. RayBatch.last_surface).")
      .def(
          "find_path",
          [](const compile::CompiledSystem& s, const std::string& name) -> std::optional<int> {
            const auto id = s.find_path(name);
            if (!id) return std::nullopt;
            return static_cast<int>(id->index);
          },
          "name"_a, "Index of the path with this name, or None.");

  m.def(
      "compile",
      [](const model::System& system, const material::MaterialLibrary* materials) {
        if (materials != nullptr) return compile::compile(system, *materials);
        const material::MaterialLibrary defaults;
        return compile::compile(system, defaults);
      },
      "system"_a, "materials"_a.none() = nb::none(),
      "Compiles a System. Without `materials` only VACUUM, AIR and CONST: references resolve."
      "\n\nRaises CompileError with the diagnostics for invalid models and unknown materials.");
}

}  // namespace rtt::py
