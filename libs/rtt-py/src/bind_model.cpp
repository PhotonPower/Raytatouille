#include <nanobind/nanobind.h>
#include <nanobind/stl/filesystem.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/string_view.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

#include <filesystem>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "bindings.hpp"
#include "rtt/diagnostics/codes.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/model/system.hpp"
#include "rtt/model/validate.hpp"

namespace nb = nanobind;
using namespace nb::literals;

namespace rtt::py {
namespace {

/// Runs a file operation of rtt-io; I/O failures (std::runtime_error that is not a ParseError)
/// become OSError, ParseError is translated in errors.cpp.
template <typename F>
auto with_os_error(F&& f) {
  try {
    return f();
  } catch (const io::ParseError&) {
    throw;
  } catch (const std::runtime_error& e) {
    PyErr_SetString(PyExc_OSError, e.what());
    throw nb::python_error();
  }
}

}  // namespace

void bind_model(nb::module_& m) {
  nb::enum_<model::Severity>(m, "Severity", "Severity of a Diagnostic.")
      .value("ERROR", model::Severity::Error)
      .value("WARNING", model::Severity::Warning);

  nb::class_<model::Diagnostic>(m, "Diagnostic", "Result of a semantic check of a System.")
      .def_ro("severity", &model::Diagnostic::severity)
      .def_ro("code", &model::Diagnostic::code,
              "Stable code, e.g. \"material.unknown\" (docs/diagnostics.md, ADR 0022).")
      .def_ro("location", &model::Diagnostic::location,
              "JSON pointer into the system file, e.g. \"/root/children/0\".")
      .def_ro("message", &model::Diagnostic::message)
      .def("__str__", [](const model::Diagnostic& d) { return model::to_string(d); })
      .def("__repr__",
           [](const model::Diagnostic& d) { return "Diagnostic(" + model::to_string(d) + ")"; });

  m.def(
      "diagnostic_codes",
      [] {
        std::vector<std::tuple<std::string, model::Severity, std::string, std::string>> codes;
        for (const diagnostics::CodeInfo& info : diagnostics::kCodes) {
          codes.emplace_back(std::string(info.code), info.severity, std::string(info.producer),
                             std::string(info.summary));
        }
        return codes;
      },
      "Registry of the diagnostic codes as (code, severity, producer, summary), sorted by code "
      "(use raytatouille.diagnostics.CODES).");

  nb::class_<model::Environment>(m, "Environment", "Surroundings of the system.")
      .def(nb::init<>())
      .def_rw("temperature_c", &model::Environment::temperature_c, "Temperature in degC.")
      .def_rw("pressure_atm", &model::Environment::pressure_atm, "Air pressure in atm.")
      .def_rw("medium", &model::Environment::medium,
              "Material reference of the surrounding medium, e.g. \"AIR\" or \"VACUUM\".")
      .def("__eq__", &equal<model::Environment>, "other"_a)
      .attr("__hash__") = nb::none();  // mutable with value equality: not hashable

  nb::class_<model::Field>(
      m, "Field",
      "Field value: degree for field angles, mm for object and paraxial image heights (the "
      "system's field type decides).")
      .def(
          "__init__",
          [](model::Field* f, double x, double y, double weight) {
            new (f) model::Field{x, y, weight};
          },
          "x"_a = 0.0, "y"_a = 0.0, "weight"_a = 1.0)
      .def_rw("x", &model::Field::x)
      .def_rw("y", &model::Field::y)
      .def_rw("weight", &model::Field::weight, "Field weight, dimensionless.")
      .def("__eq__", &equal<model::Field>, "other"_a)
      .def("__repr__",
           [](const model::Field& f) {
             return nb::str("Field(x={!r}, y={!r}, weight={!r})").format(f.x, f.y, f.weight);
           })
      .attr("__hash__") = nb::none();  // mutable with value equality: not hashable

  nb::class_<model::Wavelength>(m, "Wavelength", "One system wavelength (read-only copy).")
      .def_ro("um", &model::Wavelength::um, "Vacuum wavelength in um.")
      .def_ro("weight", &model::Wavelength::weight, "Weight, dimensionless.")
      .def_ro("reference", &model::Wavelength::reference,
              "True for the reference wavelength (exactly one per system).")
      .def("__eq__", &equal<model::Wavelength>, "other"_a)
      .def("__hash__", [](const model::Wavelength& w) {
        // Read-only, so hashable by value, consistent with __eq__ (ADR 0024).
        return nb::hash(nb::make_tuple(w.um, w.weight, w.reference));
      });

  bind_model_tree(m);  // the tree types are needed by the properties of System

  nb::class_<model::System>(
      m, "System",
      "Optical system as described by a system file (*.rtt.json). The whole tree is readable "
      "as immutable copies (raytatouille.model); only the name and the environment can be "
      "changed from Python.")
      .def(nb::init<>())
      .def_rw("name", &model::System::name)
      .def_ro("schema_version", &model::System::schema_version,
              "File format version of the model, e.g. \"0.2.0\".")
      .def_rw("environment", &model::System::environment,
              "Surroundings; changes act on this System.")
      .def_prop_ro(
          "wavelengths", [](const model::System& s) { return s.wavelengths; },
          "System wavelengths in model order (copies).")
      .def_prop_ro(
          "object_space", [](const model::System& s) { return s.object; },
          "Position of the object (copy; file key \"object\").")
      .def_prop_ro(
          "aperture", [](const model::System& s) { return s.aperture; }, "System aperture (copy).")
      .def_prop_ro(
          "fields", [](const model::System& s) { return s.fields; }, "Field points (copy).")
      .def_prop_ro(
          "root", [](const model::System& s) { return s.root; },
          "Root assembly of the element tree (copy). Each access copies the whole tree; keep "
          "the result in a variable instead of reading it again in a loop.")
      .def_prop_ro(
          "paths", [](const model::System& s) { return s.paths; }, "Ray paths in order (copies).")
      .def("to_json", &io::to_json,
           "Canonical JSON text (2-space indent, LF, trailing newline, defaults omitted).")
      .def_static("from_json", &io::parse_system, "text"_a,
                  "Parses a system from JSON text.\n\nRaises ParseError for structural errors.")
      .def(
          "save",
          [](const model::System& s, const std::filesystem::path& file) {
            with_os_error([&] {
              io::save_system(s, file);
              return 0;
            });
          },
          "file"_a, "Writes the canonical JSON text to `file`.\n\nRaises OSError on failure.")
      .def("__eq__", &equal<model::System>, "other"_a)
      .attr("__hash__") = nb::none();  // mutable with value equality: not hashable

  m.def(
      "load",
      [](const std::filesystem::path& file) {
        return with_os_error([&] { return io::load_system(file); });
      },
      "file"_a,
      "Reads a system file (*.rtt.json).\n\nRaises ParseError for structural errors and OSError "
      "if the file cannot be read.");
  m.def(
      "save",
      [](const model::System& s, const std::filesystem::path& file) {
        with_os_error([&] {
          io::save_system(s, file);
          return 0;
        });
      },
      "system"_a, "file"_a,
      "Writes the canonical JSON text of `system` to `file`.\n\nRaises OSError on failure.");
  m.def("validate", &model::validate, "system"_a,
        "Runs all semantic checks and returns the diagnostics (errors and warnings); an empty "
        "list means the model is consistent. Does not raise for an inconsistent model.");
}

}  // namespace rtt::py
