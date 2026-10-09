// to_dict() and to_json() on the bound result classes (ADR 0023), to_csv() on the reports
// (#177). The conversions are in raytatouille/results.py and raytatouille/reports.py; these
// methods only forward to them, so that they appear in the stubs.

#include <nanobind/nanobind.h>

#include <initializer_list>

#include "bindings.hpp"

namespace nb = nanobind;
using namespace nb::literals;

namespace rtt::py {

void bind_result_methods(nb::module_& m) {
  // Top-level result types and the nested types that appear inside their data.
  for (const char* name : {"SpotDiagram",
                           "RayFan",
                           "OpdMap",
                           "OpdFan",
                           "LongitudinalColour",
                           "LateralColour",
                           "DistortionSweep",
                           "DistortionPoint",
                           "FieldCurvatureSweep",
                           "FieldCurvaturePoint",
                           "PathTransmission",
                           "PathOplDifference",
                           "GhostRanking",
                           "RaytraceReport",
                           "SystemReport",
                           "DimensionReport",
                           "FirstOrder",
                           "Seidel",
                           "Prescription",
                           "TraceStats",
                           "RayBatch",
                           "RayPaths",
                           "CompiledElement",
                           "Diagnostic",
                           "LoadWarning"}) {
    const nb::object cls = m.attr(name);
    cls.attr("to_dict") = nb::cpp_function(
        [](nb::handle self) {
          return nb::module_::import_("raytatouille.results").attr("to_dict")(self);
        },
        nb::is_method(), nb::scope(cls), nb::name("to_dict"),
        nb::sig("def to_dict(self) -> dict[str, typing.Any]"),
        "The result as plain Python values (ADR 0023): dicts, lists, numbers, strings, None, "
        "NumPy arrays; enums as their names. See raytatouille.results.");
    cls.attr("to_json") = nb::cpp_function(
        [](nb::handle self, nb::handle indent) {
          return nb::module_::import_("raytatouille.results").attr("to_json")(self, indent);
        },
        nb::is_method(), nb::scope(cls), nb::name("to_json"), "indent"_a = nb::none(),
        nb::sig("def to_json(self, indent: int | None = None) -> str"),
        "The result as JSON text in the format raytatouille-result (ADR 0023); read it back "
        "with raytatouille.results.load_json.");
  }
  // CSV of the reports (#177); the conversion is in raytatouille/reports.py.
  for (const char* name : {"RaytraceReport", "SystemReport", "DimensionReport"}) {
    const nb::object cls = m.attr(name);
    cls.attr("to_csv") = nb::cpp_function(
        [](nb::handle self) {
          return nb::module_::import_("raytatouille.reports").attr("to_csv")(self);
        },
        nb::is_method(), nb::scope(cls), nb::name("to_csv"), nb::sig("def to_csv(self) -> str"),
        "The report as CSV text (#177): comma-separated with a header line; floats in the "
        "shortest form that reads back to the same double, NaN and None as empty fields. See "
        "raytatouille.reports.");
  }
}

}  // namespace rtt::py
