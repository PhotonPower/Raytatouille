#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/variant.h>
#include <nanobind/stl/vector.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <vector>

#include "bindings.hpp"
#include "rtt/paraxial/paraxial.hpp"

namespace nb = nanobind;
using namespace nb::literals;

namespace rtt::py {

namespace {

/// Read-only NumPy copy (x, y, z) of a vector.
ReadOnlyArray<double> vec3_copy(const math::Vec3& v) {
  const std::array<double, 3> xyz = {v.x(), v.y(), v.z()};
  return read_only_array<double>(xyz, [](double x) { return x; });
}

/// Read-only NumPy copy (3, 3) of a matrix, row-major.
nb::ndarray<nb::numpy, const double, nb::ndim<2>, nb::c_contig> mat3_copy(const math::Mat3& m) {
  auto data = std::make_unique<std::vector<double>>();
  data->reserve(9);
  for (Eigen::Index r = 0; r < 3; ++r) {
    for (Eigen::Index c = 0; c < 3; ++c) data->push_back(m(r, c));
  }
  const std::size_t shape[2] = {3, 3};
  const double* values = data->data();
  const nb::capsule owner = owning_capsule(std::move(data));
  return {values, 2, shape, owner};
}

}  // namespace

void bind_paraxial(nb::module_& m) {
  using paraxial::FirstOrder;
  using paraxial::Pupil;
  using paraxial::RayAtEvent;
  using paraxial::Section;
  using paraxial::UnfoldedFrame;

  nb::class_<UnfoldedFrame>(
      m, "UnfoldedFrame",
      "Map from the unfolded frame of a path section to global coordinates behind plane fold "
      "mirrors (ADR 0054): global = linear @ p + offset, mm. linear is orthogonal with "
      "determinant -1 after an odd number of folds; identity and zero offset without folds.")
      .def_prop_ro(
          "linear", [](const UnfoldedFrame& f) { return mat3_copy(f.linear); },
          nb::rv_policy::reference, "linear (3 x 3), float64, row-major (copy).")
      .def_prop_ro(
          "offset", [](const UnfoldedFrame& f) { return vec3_copy(f.offset); },
          nb::rv_policy::reference, "offset (x, y, z) in mm, float64 (copy).")
      .def(
          "to_global",
          [](const UnfoldedFrame& f, const std::array<double, 3>& p) {
            return vec3_copy(f.to_global(math::Vec3(p[0], p[1], p[2])));
          },
          "p"_a, "Global point (x, y, z) in mm of the unfolded point p (mm).")
      .def(
          "to_unfolded",
          [](const UnfoldedFrame& f, const std::array<double, 3>& g) {
            return vec3_copy(f.to_unfolded(math::Vec3(g[0], g[1], g[2])));
          },
          "g"_a, "Unfolded point (x, y, z) in mm of the global point g (mm).");

  nb::class_<Pupil>(m, "Pupil", "Paraxial image of the aperture stop.")
      .def_ro("z", &Pupil::z,
              "Global z in mm (the exit pupil behind a plane fold mirror: along the unfolded "
              "axis, ADR 0054); None if the pupil is at infinity.")
      .def_ro("diameter", &Pupil::diameter,
              "Diameter in mm; None if the system aperture does not define it.");

  nb::class_<FirstOrder>(
      m, "FirstOrder",
      "First-order data of a path at one wavelength. Indices are absolute (AIR is Ciddor air), "
      "EFL = 1/power. Positions are global z coordinates in mm; behind a plane fold mirror "
      "(ADR 0054) along the unfolded axis, mapped to global by image_space. None where a value "
      "is not defined (afocal system, object at infinity, no stop, ...).")
      .def_ro("object_index", &FirstOrder::object_index, "|n| in object space.")
      .def_ro("image_index", &FirstOrder::image_index, "|n| in image space.")
      .def_ro("image_direction", &FirstOrder::image_direction,
              "+1 if light leaves the system towards +z, -1 towards -z.")
      .def_ro("power", &FirstOrder::power,
              "Power in 1/mm, as computed; for an afocal system 0 up to rounding (below the "
              "afocal threshold of first_order, see the C++ documentation).")
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
      .def_ro("exit_pupil", &FirstOrder::exit_pupil, "None if the path has no stop.")
      .def_ro("folds", &FirstOrder::folds,
              "Number of plane fold mirrors on the path (ADR 0054); 0 along the global z axis.")
      .def_ro("image_space", &FirstOrder::image_space,
              "UnfoldedFrame of image space: maps the unfolded positions behind the folds to "
              "global coordinates (identity without folds).");

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

  // Paraxial computation per section for anamorphic paths (ADR 0042, point 2).
  nb::enum_<Section>(m, "Section",
                     "Section of a path for the paraxial computation of anamorphic systems "
                     "(ADR 0042): the global x-z or y-z plane through the z axis.")
      .value("X", Section::x, "the x-z plane; heights are global x, slopes dx/dz")
      .value("Y", Section::y, "the y-z plane; heights are global y, slopes dy/dz");

  nb::class_<RayAtEvent>(m, "RayAtEvent", "Paraxial ray just after one event of a path.")
      .def_ro("z", &RayAtEvent::z,
              "Global z of the surface vertex, mm; behind a plane fold mirror along the "
              "unfolded axis (ADR 0054).")
      .def_ro("y", &RayAtEvent::y,
              "Ray height at the vertex plane in mm: global y, or global x in the section X.")
      .def_ro("u", &RayAtEvent::u, "Slope dh/dz after the event in the global frame.")
      .def_ro("n", &RayAtEvent::n,
              "Signed index after the event: positive while light travels towards +z, negative "
              "towards -z.");

  m.def(
      "first_order_section",
      [](const compile::CompiledSystem& system, Section section, const PathArg& path,
         std::optional<std::uint16_t> wavelength) {
        return paraxial::first_order_section(system, path_id(system, path),
                                             wavelength_index(system, wavelength), section);
      },
      "system"_a, "section"_a, "path"_a = 0, "wavelength"_a.none() = nb::none(),
      "First-order data of `path` (index or name) in the section `section` (Section.X or "
      "Section.Y) at wavelength index `wavelength` (None: reference wavelength), ADR 0042. "
      "The same quantities as first_order(); heights and pupil diameters are measured along x "
      "or y. Anamorphic surfaces (biconic, toroid, cylinder) and the ideal cylinder lens enter "
      "with their curvature or power in the section; a section without power is afocal (efl "
      "None). For a rotationally symmetric path without folds both sections equal "
      "first_order() bit for bit.\n\n"
      "Raises ParaxialError if a surface is decentred or tilted, for a path folded by plane "
      "mirrors (ADR 0054, point 3), if the principal axes of an "
      "anamorphic surface or the direction of power of an ideal cylinder lens are not aligned "
      "with x and y, or if the path or wavelength does not exist; ValueError for an unknown path "
      "name.");

  m.def(
      "trace_ray_section",
      [](const compile::CompiledSystem& system, Section section, double z0, double h0, double u0,
         const PathArg& path, std::optional<std::uint16_t> wavelength) {
        return paraxial::trace_ray_section(system, path_id(system, path),
                                           wavelength_index(system, wavelength), z0, h0, u0,
                                           section);
      },
      "system"_a, "section"_a, "z0"_a, "h0"_a, "u0"_a, "path"_a = 0,
      "wavelength"_a.none() = nb::none(),
      "Traces one paraxial ray along `path` in the section `section` (ADR 0042): start at the "
      "global z0 in mm in the medium before the first event, travelling towards +z, with height "
      "h0 in mm (global x for Section.X, y for Section.Y) and slope u0 = dh/dz. Returns the ray "
      "after each event of the path, in path order (list of RayAtEvent).\n\n"
      "Raises ParaxialError as first_order_section().");

  m.def(
      "is_anamorphic",
      [](const compile::CompiledSystem& system, const PathArg& path) {
        return paraxial::is_anamorphic(system, path_id(system, path));
      },
      "system"_a, "path"_a = 0,
      "True if an event of `path` (index or name) is at a surface with unequal curvatures "
      "d2z/dx2 != d2z/dy2 or a mixed term d2z/dxdy != 0 at its vertex in its local frame "
      "(biconic, toroid, cylinder), or at an ideal cylinder lens (ADR 0042). Such a path has "
      "no first_order(); use first_order_section().\n\n"
      "Raises ParaxialError if the path does not exist, ValueError for an unknown path name.");
}

}  // namespace rtt::py
