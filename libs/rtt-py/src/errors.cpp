#include "rtt/compile/errors.hpp"

#include <nanobind/nanobind.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <exception>
#include <optional>
#include <utility>

#include "bindings.hpp"
#include "rtt/analysis/spot.hpp"
#include "rtt/coating/catalog.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/edit.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/agf.hpp"
#include "rtt/material/material.hpp"
#include "rtt/paraxial/paraxial.hpp"
#include "rtt/trace/run_control.hpp"

namespace nb = nanobind;

namespace rtt::py {
namespace {

/// Raises raytatouille.errors.<name>(*args) as the active Python exception. The classes are
/// pure Python (multiple inheritance, e.g. ParseError(RaytatouilleError, ValueError)); the module
/// is imported on first use and then cached in sys.modules.
template <typename... Args>
void raise(const char* name, Args&&... args) {
  const nb::object cls = nb::module_::import_("raytatouille.errors").attr(name);
  const nb::object error = cls(std::forward<Args>(args)...);
  PyErr_SetObject(cls.ptr(), error.ptr());
}

/// None or the value as a Python object.
template <typename T>
nb::object optional(const std::optional<T>& value) {
  return value ? nb::cast(*value) : nb::none();
}

nb::object optional(const std::optional<model::SurfaceId>& id) {
  return id ? nb::cast(id->str()) : nb::none();
}

}  // namespace

void register_errors(nb::module_& /*m*/) {
  // Exceptions not caught here propagate to nanobind's default translation
  // (std::invalid_argument -> ValueError, std::out_of_range -> IndexError, ...).
  nb::register_exception_translator([](const std::exception_ptr& p, void* /*payload*/) {
    try {
      std::rethrow_exception(p);
    } catch (const io::EditError& e) {
      // Before nanobind's std::invalid_argument -> ValueError (ADR 0024).
      raise("EditError", e.what(), e.code(), e.location(), optional(e.op_index()), e.diagnostics());
    } catch (const io::ParseError& e) {
      raise("ParseError", e.what(), e.pointer());
    } catch (const compile::CompileError& e) {
      raise("CompileError", e.what(), e.diagnostics());
    } catch (const compile::NoStopError& e) {
      // Before nanobind's std::invalid_argument -> ValueError; the Python class also derives
      // from ParaxialError and AnalysisError (ADR 0022).
      raise("NoStopError", e.what(), e.path_name(), e.location());
    } catch (const paraxial::ParaxialError& e) {
      raise("ParaxialError", e.what(), optional(e.surface()), optional(e.location()));
    } catch (const material::UnknownMaterial& e) {
      raise("UnknownMaterial", e.what());
    } catch (const material::AgfError& e) {
      raise("AgfError", e.what(), e.file(), e.line());
    } catch (const analysis::AnalysisError& e) {
      raise("AnalysisError", e.what(), optional(e.surface()), optional(e.location()),
            optional(e.ray_status()), optional(e.field()), optional(e.wavelength()));
    } catch (const coating::CoatingCatalogError& e) {
      raise("CoatingCatalogError", e.what(), e.file(), e.pointer());
    } catch (const trace::Cancelled& e) {
      raise("Cancelled", e.what());
    }
  });
}

}  // namespace rtt::py
