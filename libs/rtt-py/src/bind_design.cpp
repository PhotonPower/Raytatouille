// Lens catalogues (ADR 0053, point 3): the loader of rtt-io for rt.design.

#include <nanobind/nanobind.h>
#include <nanobind/stl/filesystem.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/string_view.h>
#include <nanobind/stl/vector.h>

#include <string>

#include "bindings.hpp"
#include "rtt/io/lens_catalog.hpp"
#include "rtt/model/system.hpp"

namespace nb = nanobind;
using namespace nb::literals;

namespace rtt::py {

void bind_design(nb::module_& m) {
  m.attr("SCHEMA_VERSION") = std::string(model::kSchemaVersion);

  nb::class_<io::LensEntry>(
      m, "LensEntry",
      "One lens of a lens catalogue (ADR 0053), read-only. Lengths in mm; radii in order along "
      "+z with the sign convention of the system format (centre of curvature in +z: positive), "
      "None is plano.")
      .def_ro("name", &io::LensEntry::name, "Name, unique within the catalogue.")
      .def_ro("description", &io::LensEntry::description, "Free text, may be empty.")
      .def_prop_ro(
          "kind", [](const io::LensEntry& e) { return std::string(io::lens_kind_name(e.kind)); },
          "\"singlet\", \"doublet\" (cemented) or \"triplet\".")
      .def_ro("radii", &io::LensEntry::radii, "Vertex radii in mm; None is plano.")
      .def_ro("thicknesses", &io::LensEntry::thicknesses,
              "Centre thicknesses of the glass segments in mm.")
      .def_ro("gaps", &io::LensEntry::gaps,
              "Triplet: the two air gaps in mm (vertex to vertex); empty otherwise.")
      .def_ro("materials", &io::LensEntry::materials,
              "One material reference CATALOG:NAME per glass segment.")
      .def_ro("diameter", &io::LensEntry::diameter,
              "Clear diameter of every surface in mm, or None.")
      .def("__repr__", [](const io::LensEntry& e) {
        return "LensEntry(" + e.name + ", " + std::string(io::lens_kind_name(e.kind)) + ")";
      });

  nb::class_<io::LensCatalog>(m, "LensCatalog",
                              "Contents of a lens catalogue file (ADR 0053), read-only.")
      .def_ro("name", &io::LensCatalog::name, "Catalogue name.")
      .def_ro("lenses", &io::LensCatalog::lenses, "The lenses in file order.")
      .def(
          "__getitem__",
          [](const io::LensCatalog& c, const std::string& name) {
            for (const io::LensEntry& e : c.lenses) {
              if (e.name == name) return e;
            }
            throw nb::key_error(("no lens '" + name + "' in catalog '" + c.name + "'").c_str());
          },
          "name"_a, "The lens with this name; KeyError if there is none.")
      .def("__len__", [](const io::LensCatalog& c) { return c.lenses.size(); });

  m.def("load_lens_catalog", &io::load_lens_catalog, "file"_a,
        "Reads a lens catalogue file (format raytatouille-lenses, ADR 0053). Raises "
        "LensCatalogError (with file and JSON pointer) if it cannot be read or is malformed.");
  m.def("parse_lens_catalog", &io::parse_lens_catalog, "text"_a, "file"_a = std::string(),
        "Parses a lens catalogue from JSON text; `file` names it in error messages. Raises "
        "LensCatalogError.");
}

}  // namespace rtt::py
