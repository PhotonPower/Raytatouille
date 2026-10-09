#include <nanobind/nanobind.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/filesystem.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/variant.h>
#include <nanobind/stl/vector.h>

#include <cstddef>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <variant>
#include <vector>

#include "bindings.hpp"
#include "rtt/coating/catalog.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/compile/ghosts.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/path.hpp"
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

void warn(const std::vector<model::Diagnostic>& diagnostics) {
  for (const model::Diagnostic& d : diagnostics) {
    if (d.severity != model::Severity::Warning) continue;
    const nb::object category =
        nb::module_::import_("raytatouille.errors").attr("RaytatouilleWarning");
    nb::module_::import_("warnings")
        .attr("warn")(category(d.message, d.code, d.location), "stacklevel"_a = 1);
  }
}

namespace {

/// The ghost paths of a GhostSystem, bound with columns() (#123, #133).
struct GhostPaths {
  std::vector<compile::GhostPath> points;
};

/// Name of the model path `base` (a name, or an index into system.paths) for
/// compile_with_ghosts(), which takes a name.
/// @throws std::invalid_argument (ValueError) for an index that is not a path of `system`
std::string base_name(const model::System& system, const PathArg& base) {
  if (const auto* name = std::get_if<std::string>(&base)) return *name;
  const std::uint32_t index = std::get<std::uint32_t>(base);
  if (index >= system.paths.size()) {
    throw std::invalid_argument("ghosts: path index " + std::to_string(index) + " does not exist");
  }
  return system.paths[index].name;
}

}  // namespace

void bind_compile(nb::module_& m) {
  bind_material(m);  // before compile(), which takes a MaterialLibrary

  nb::class_<coating::CoatingLibrary>(
      m, "CoatingLibrary",
      "Resolves coating references \"CATALOG:NAME\" to thin-film designs from JSON coating "
      "catalogues (ADR 0019; format raytatouille-coatings). Pass it to compile().")
      .def(nb::init<>())
      .def("add_catalog", &coating::CoatingLibrary::add_catalog, "path"_a,
           "Loads coating catalogues: a .json file, or a directory whose *.json files are "
           "loaded in sorted order. All or nothing: on an error nothing is added.\n\nRaises "
           "CoatingCatalogError for a malformed file (with file and JSON pointer) and ValueError "
           "if the path has no .json file or a catalogue name is already registered.")
      .def(
          "__contains__",
          [](const coating::CoatingLibrary& lib, const std::string& reference) {
            return lib.find(reference) != nullptr;
          },
          "reference"_a, "True if the reference \"CATALOG:NAME\" resolves.");

  nb::class_<compile::CompiledMedium>(m, "CompiledMedium", "Medium evaluated at all wavelengths.")
      .def_ro("reference", &compile::CompiledMedium::reference,
              "Material reference as in the model, e.g. \"AIR\".")
      .def_ro("index", &compile::CompiledMedium::index,
              "Complex index n + i*kappa per system wavelength (kappa >= 0 absorbs); for a "
              "crystal the ordinary index n_O (kappa = 0).")
      .def_prop_ro("is_crystal", &compile::CompiledMedium::is_crystal,
                   "True for a uniaxial crystal (ADR 0026).")
      .def_ro("reference_extraordinary", &compile::CompiledMedium::reference_extraordinary,
              "Crystal only (empty otherwise): reference of the extraordinary index n_E.")
      .def_prop_ro(
          "index_extraordinary",
          [](const compile::CompiledMedium& c) {
            return read_only_array<double>(c.index_extraordinary, [](double v) { return v; });
          },
          nb::rv_policy::reference,
          "Crystal only (empty otherwise): n_E per system wavelength, real, float64 (copy).")
      .def_prop_ro(
          "optic_axis",
          [](const compile::CompiledMedium& c) -> std::optional<ReadOnlyArray<double>> {
            if (!c.optic_axis) return std::nullopt;
            const std::vector<double> xyz = {c.optic_axis->x(), c.optic_axis->y(),
                                             c.optic_axis->z()};
            return read_only_array<double>(xyz, [](double v) { return v; });
          },
          nb::rv_policy::reference,
          "Crystal only: optic axis (x, y, z) in global coordinates, unit vector, float64 "
          "(copy); its sign has no meaning (ADR 0026, point 2). None for an isotropic medium.");

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
      .def_prop_ro(
          "surface_locations",
          [](const compile::CompiledSystem& s) {
            std::vector<std::string> locations;
            for (const auto& surface : s.surfaces()) locations.push_back(surface.location);
            return locations;
          },
          "JSON pointers of the surfaces in the system file, same order as surface_ids, e.g. "
          "\"/root/children/1/surfaces/0\" (ADR 0022).")
      .def_prop_ro("diagnostics", &compile::CompiledSystem::diagnostics,
                   "Warnings found while compiling, with code and JSON pointer (ADR 0022).")
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
      [](const model::System& system, const material::MaterialLibrary* materials,
         const coating::CoatingLibrary* coatings) {
        const material::MaterialLibrary default_materials;
        const coating::CoatingLibrary no_coatings;
        compile::CompiledSystem compiled =
            compile::compile(system, materials != nullptr ? *materials : default_materials,
                             coatings != nullptr ? *coatings : no_coatings);
        warn(compiled.diagnostics());
        return compiled;
      },
      "system"_a, "materials"_a.none() = nb::none(), "coatings"_a.none() = nb::none(),
      "Compiles a System. Without `materials` only VACUUM, AIR and CONST: references resolve; "
      "without `coatings` a surface with a coating reference is an error (ADR 0019).\n\nRaises "
      "CompileError with the diagnostics for invalid models, unknown materials and unknown "
      "coatings.\n\nEvery warning (CompiledSystem.diagnostics) is also issued as a "
      "raytatouille.errors.RaytatouilleWarning with its code and location.");

  // Ghost generator (#123, ADR 0027; Python #133).
  auto ghost_paths = columns<GhostPaths>(
      m, "GhostPaths",
      "Which path of GhostSystem.system is which ghost (ADR 0027), one entry per ghost in the "
      "order of ghost_paths(); read-only NumPy copies.");
  column<std::uint32_t>(
      ghost_paths, "path", [](const compile::GhostPath& g) { return g.path.index; },
      "Index of the ghost path in GhostSystem.system (copy).");
  column<std::uint32_t>(
      ghost_paths, "base", [](const compile::GhostPath& g) { return g.base.index; },
      "Index of the base path it was derived from (copy).");
  column<std::uint32_t>(
      ghost_paths, "surface_j", [](const compile::GhostPath& g) { return g.surface_j; },
      "Surface of the first ghost reflection (back), index into surface_ids (copy).");
  column<std::uint32_t>(
      ghost_paths, "surface_i", [](const compile::GhostPath& g) { return g.surface_i; },
      "Surface of the second ghost reflection (forward), index into surface_ids (copy).");
  column<std::uint64_t>(
      ghost_paths, "event_j", [](const compile::GhostPath& g) { return g.event_j; },
      "Index of the reflecting event j in the base path (copy).");
  column<std::uint64_t>(
      ghost_paths, "event_i", [](const compile::GhostPath& g) { return g.event_i; },
      "Index of the reflecting event i in the base path (copy).");

  nb::class_<compile::GhostSystem>(
      m, "GhostSystem",
      "A system compiled with the ghosts of one base path (ADR 0027): `system` has all paths of "
      "the model, then the ghosts; `ghosts` tells which path is which ghost. Pass it to "
      "raytatouille.analysis.ghost_ranking().")
      .def_ro("system", &compile::GhostSystem::system,
              "The compiled copy: all paths of the model, then the ghost paths.")
      .def_prop_ro(
          "ghosts", [](const compile::GhostSystem& g) { return GhostPaths{g.ghosts}; },
          "One entry per ghost, in the order of ghost_paths().");

  m.def(
      "ghost_paths",
      [](const compile::CompiledSystem& system, const PathArg& base, std::size_t max_paths) {
        return compile::ghost_paths(system, path_id(system, base),
                                    compile::GhostOptions{max_paths});
      },
      "system"_a, "base"_a = 0, nb::kw_only(), "max_paths"_a = compile::GhostOptions{}.max_paths,
      "Two-reflection ghost paths of the path `base` (index or name) of a CompiledSystem, as "
      "explicit model paths (raytatouille.model.Path; ADR 0027). A ghost reflects at refracting "
      "surface j back and at refracting surface i < j forward again, then continues to the "
      "image surface: base[0..j-1], Reflect at j, base[i+1..j-1] reversed, Reflect at i, "
      "base[i+1..end]. N Refract events give N (N - 1) / 2 ghosts, j ascending, then i "
      "ascending; names \"<base> ghost <surface j>/<surface i>\" (\"#k\" with the event index "
      "for a surface that occurs more than once). The ghosts are derived data, never written "
      "into a file by raytatouille.\n\nRaises ValueError for an unknown path, a base path with "
      "a diffraction order or a crystal mode (not supported in M4), one that enters an "
      "element from outside through an inner surface, more than `max_paths` ghosts, or a "
      "ghost name that is already a path name or occurs twice.");
  m.def(
      "compile_with_ghosts",
      [](const model::System& system, const PathArg& base,
         const material::MaterialLibrary* materials, const coating::CoatingLibrary* coatings,
         std::size_t max_paths) {
        const material::MaterialLibrary default_materials;
        const coating::CoatingLibrary no_coatings;
        compile::GhostSystem ghosts = compile::compile_with_ghosts(
            system, base_name(system, base), materials != nullptr ? *materials : default_materials,
            coatings != nullptr ? *coatings : no_coatings, compile::GhostOptions{max_paths});
        warn(ghosts.system.diagnostics());
        return ghosts;
      },
      "system"_a, "base"_a = 0, "materials"_a.none() = nb::none(), "coatings"_a.none() = nb::none(),
      nb::kw_only(), "max_paths"_a = compile::GhostOptions{}.max_paths,
      "Compiles a System with the two-reflection ghosts of its path `base` (index into "
      "System paths or name; ADR 0027): compiles the system, derives the ghosts with "
      "ghost_paths(), appends them to a copy of the model and compiles the copy. The System "
      "itself is not changed. `materials` and `coatings` as for compile(). Pass the result to "
      "raytatouille.analysis.ghost_ranking().\n\nRaises CompileError as compile(), ValueError "
      "for an unknown path and as ghost_paths().\n\nEvery compile warning is also issued as a "
      "raytatouille.errors.RaytatouilleWarning, as by compile().");
}

}  // namespace rtt::py
