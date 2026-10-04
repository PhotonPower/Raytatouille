#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <exception>
#include <utility>

#include "bindings.hpp"
#include "rtt/analysis/spot.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/agf.hpp"
#include "rtt/material/material.hpp"
#include "rtt/paraxial/paraxial.hpp"

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

}  // namespace

void register_errors(nb::module_& /*m*/) {
  // Exceptions not caught here propagate to nanobind's default translation
  // (std::invalid_argument -> ValueError, std::out_of_range -> IndexError, ...).
  nb::register_exception_translator([](const std::exception_ptr& p, void* /*payload*/) {
    try {
      std::rethrow_exception(p);
    } catch (const io::ParseError& e) {
      raise("ParseError", e.what(), e.pointer());
    } catch (const compile::CompileError& e) {
      raise("CompileError", e.what(), e.diagnostics());
    } catch (const paraxial::ParaxialError& e) {
      raise("ParaxialError", e.what());
    } catch (const material::UnknownMaterial& e) {
      raise("UnknownMaterial", e.what());
    } catch (const material::AgfError& e) {
      raise("AgfError", e.what(), e.file(), e.line());
    } catch (const analysis::AnalysisError& e) {
      raise("AnalysisError", e.what());
    }
  });
}

}  // namespace rtt::py
