#include <nanobind/nanobind.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/variant.h>

#include <cstdint>
#include <optional>

#include "bindings.hpp"
#include "rtt/paraxial/paraxial.hpp"

namespace nb = nanobind;
using namespace nb::literals;

namespace rtt::py {

void bind_paraxial(nb::module_& m) {
  using paraxial::FirstOrder;
  using paraxial::Pupil;

  nb::class_<Pupil>(m, "Pupil", "Paraxial image of the aperture stop.")
      .def_ro("z", &Pupil::z, "Global z in mm; None if the pupil is at infinity.")
      .def_ro("diameter", &Pupil::diameter,
              "Diameter in mm; None if the system aperture does not define it.");

  nb::class_<FirstOrder>(
      m, "FirstOrder",
      "First-order data of a path at one wavelength. Indices are absolute (AIR is Ciddor air), "
      "EFL = 1/power. Positions are global z coordinates in mm; None where a value is not "
      "defined (afocal system, object at infinity, no stop, ...).")
      .def_ro("object_index", &FirstOrder::object_index, "|n| in object space.")
      .def_ro("image_index", &FirstOrder::image_index, "|n| in image space.")
      .def_ro("image_direction", &FirstOrder::image_direction,
              "+1 if light leaves the system towards +z, -1 towards -z.")
      .def_ro("power", &FirstOrder::power, "Power in 1/mm; 0 for afocal systems.")
      .def_ro("efl", &FirstOrder::efl,
              "Effective focal length 1/power in mm, positive for converging systems.")
      .def_ro("front_focal_length", &FirstOrder::front_focal_length,
              "|n| / power in object space (H -> F), mm.")
      .def_ro("rear_focal_length", &FirstOrder::rear_focal_length,
              "|n'| / power in image space (H' -> F'), mm.")
      .def_ro("ffl", &FirstOrder::ffl, "First vertex -> F in mm, positive if F lies upstream.")
      .def_ro("bfl", &FirstOrder::bfl, "Last vertex -> F' in mm, positive if F' lies downstream.")
      .def_ro("front_focal_z", &FirstOrder::front_focal_z, "Global z of F, mm.")
      .def_ro("rear_focal_z", &FirstOrder::rear_focal_z, "Global z of F', mm.")
      .def_ro("front_principal_z", &FirstOrder::front_principal_z, "Global z of H, mm.")
      .def_ro("rear_principal_z", &FirstOrder::rear_principal_z, "Global z of H', mm.")
      .def_ro("image_z", &FirstOrder::image_z,
              "Global z of the paraxial image of the axial object point, mm.")
      .def_ro("lateral_magnification", &FirstOrder::lateral_magnification,
              "m = n u / (n' u') of the axial marginal ray; finite object only.")
      .def_ro("angular_magnification", &FirstOrder::angular_magnification,
              "Ratio u'/u of the chief-ray slopes along the propagation direction.")
      .def_ro("entrance_pupil", &FirstOrder::entrance_pupil, "None if the path has no stop.")
      .def_ro("exit_pupil", &FirstOrder::exit_pupil, "None if the path has no stop.");

  m.def(
      "first_order",
      [](const compile::CompiledSystem& system, const PathArg& path,
         std::optional<std::uint16_t> wavelength) {
        return paraxial::first_order(system, path_id(system, path),
                                     wavelength_index(system, wavelength));
      },
      "system"_a, "path"_a = 0, "wavelength"_a.none() = nb::none(),
      "First-order data (EFL, focal and principal points, pupils, magnification) of `path` "
      "(index or name) at wavelength index `wavelength` (None: reference wavelength).\n\n"
      "Raises ParaxialError if the path is not rotationally symmetric or the path or "
      "wavelength does not exist, ValueError for an unknown path name.");
}

}  // namespace rtt::py
