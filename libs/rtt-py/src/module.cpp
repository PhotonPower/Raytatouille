#include <nanobind/nanobind.h>

#include "bindings.hpp"

NB_MODULE(_core, m) {
  m.doc() =
      "C++ core of raytatouille (ADR 0002). Use the package raytatouille; this module is "
      "private. Units: lengths in mm, wavelengths in um (vacuum), angles in degree where the "
      "name ends in _deg, temperature in degC, pressure in atm.";
  rtt::py::register_errors(m);
  rtt::py::bind_model(m);
  rtt::py::bind_compile(m);
  rtt::py::bind_paraxial(m);
  rtt::py::bind_trace(m);
}
