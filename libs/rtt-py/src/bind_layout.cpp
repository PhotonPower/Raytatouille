// Geometry export for the layout (#81): rtt/compile/layout.hpp for raytatouille.layout. The
// functions here take exact dtypes and C-contiguous 1-D arrays; python/raytatouille/layout.py
// broadcasts the user's input and resolves surface ids and plane shorthands.

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <map>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include "bindings.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/compile/layout.hpp"

namespace nb = nanobind;
using namespace nb::literals;

namespace rtt::py {
namespace {

using compile::CompiledSystem;
using math::Vec3;

using RealIn = nb::ndarray<const double, nb::ndim<1>, nb::c_contig, nb::device::cpu>;
template <typename T>
using Out = nb::ndarray<nb::numpy, T>;

/// NumPy array of the given shape that owns `values`.
template <typename T>
Out<T> owned(std::vector<T>&& values, std::initializer_list<std::size_t> shape) {
  auto data = std::make_unique<std::vector<T>>(std::move(values));
  T* pointer = data->data();
  nb::capsule owner(data.get(), [](void* p) noexcept {
    const std::unique_ptr<std::vector<T>> owned_data(static_cast<std::vector<T>*>(p));
  });
  static_cast<void>(data.release());
  return Out<T>(pointer, shape.size(), shape.begin(), owner);
}

/// Model name of an element kind, as in the file format.
std::string kind_name(model::ElementKind kind) {
  switch (kind) {
    case model::ElementKind::Lens:
      return "lens";
    case model::ElementKind::Mirror:
      return "mirror";
    case model::ElementKind::Plate:
      return "plate";
    case model::ElementKind::ThinElement:
      return "thin_element";
    case model::ElementKind::Stop:
      return "stop";
    case model::ElementKind::Detector:
      return "detector";
  }
  return "unknown";
}

Vec3 vec3(const RealIn& a, const char* name) {
  if (a.shape(0) != 3) throw std::invalid_argument(std::string(name) + " must have 3 elements");
  return {a.data()[0], a.data()[1], a.data()[2]};
}

Out<double> polylines_entry(const compile::Polyline& line) {
  std::vector<double> data;
  data.reserve(3 * line.size());
  for (const Vec3& p : line) data.insert(data.end(), {p.x(), p.y(), p.z()});
  return owned(std::move(data), {line.size(), 3});
}

std::vector<Out<double>> polylines(const std::vector<compile::Polyline>& lines) {
  std::vector<Out<double>> out;
  out.reserve(lines.size());
  for (const auto& line : lines) out.push_back(polylines_entry(line));
  return out;
}

/// (shape type, curvature 1/mm, conic constant, polynomial coefficients A4, A6, ...).
std::tuple<std::string, double, double, std::vector<double>> shape_data(
    const compile::CompiledShape& shape) {
  return std::visit(
      [](const auto& s) -> std::tuple<std::string, double, double, std::vector<double>> {
        using S = std::decay_t<decltype(s)>;
        if constexpr (std::is_same_v<S, geom::Plane<double>>) {
          return {"plane", 0.0, 0.0, {}};
        } else if constexpr (std::is_same_v<S, geom::Conic<double>>) {
          return {"conic", s.curvature(), s.conic_constant(), {}};
        } else {
          const auto [c, k] = s.base_conic();
          return {"even_asphere", c, k, s.coefficients()};
        }
      },
      shape);
}

/// (aperture type, parameters in mm) or none.
std::optional<std::tuple<std::string, std::map<std::string, double>>> aperture_data(
    const std::optional<model::Aperture>& aperture) {
  if (!aperture) return std::nullopt;
  return std::visit(
      [](const auto& a) -> std::tuple<std::string, std::map<std::string, double>> {
        using A = std::decay_t<decltype(a)>;
        if constexpr (std::is_same_v<A, model::CircularAperture>) {
          return {"circular", {{"radius", a.radius}, {"inner_radius", a.inner_radius}}};
        } else if constexpr (std::is_same_v<A, model::RectangularAperture>) {
          return {"rectangular",
                  {{"half_width_x", a.half_width_x}, {"half_width_y", a.half_width_y}}};
        } else {
          return {"elliptical", {{"semi_axis_x", a.semi_axis_x}, {"semi_axis_y", a.semi_axis_y}}};
        }
      },
      *aperture);
}

}  // namespace

void bind_layout(nb::module_& m) {
  nb::class_<compile::CompiledElement>(
      m, "CompiledElement",
      "An element of a compiled system in tree order (#81): its surfaces "
      "[first_surface, first_surface + surface_count) and the media of its body.")
      .def_ro("name", &compile::CompiledElement::name, "Element name from the model.")
      .def_prop_ro(
          "kind", [](const compile::CompiledElement& e) { return kind_name(e.kind); },
          "Kind as in the file format: lens, mirror, plate, thin_element, stop, "
          "detector.")
      .def_ro("first_surface", &compile::CompiledElement::first_surface,
              "Index of the first surface of the element.")
      .def_ro("surface_count", &compile::CompiledElement::surface_count,
              "Number of surfaces of the element.")
      .def_ro("media", &compile::CompiledElement::media,
              "Media of the body, indices into CompiledSystem.media: one per segment for lenses "
              "and plates (segment i between surfaces i and i + 1), the substrate for a mirror; "
              "empty without material.")
      .def_ro("segmented", &compile::CompiledElement::segmented,
              "True if the media follow the segments (lenses, plates of different materials); "
              "only these have outlines.");

  m.def(
      "layout_elements", [](const CompiledSystem& s) { return s.elements(); }, "system"_a,
      "Elements of the system in tree order (raytatouille.layout.elements).");

  m.def(
      "layout_surfaces",
      [](const CompiledSystem& s) {
        nb::list out;
        for (const compile::CompiledSurface& c : s.surfaces()) {
          const math::Mat3& r = c.to_global.rotation();
          const Vec3 t = c.to_global.apply_point(Vec3::Zero());
          std::vector<double> rotation(9);
          for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) rotation[static_cast<std::size_t>(3 * i + j)] = r(i, j);
          }
          const auto max_radius =
              std::visit([](const auto& sh) { return sh.max_radius(); }, c.shape);
          out.append(nb::make_tuple(
              c.id.str(), c.element, kind_name(c.element_kind), owned(std::move(rotation), {3, 3}),
              owned(std::vector<double>{t.x(), t.y(), t.z()}, {3}), shape_data(c.shape), max_radius,
              aperture_data(c.aperture), c.medium_front, c.medium_back));
        }
        return out;
      },
      "system"_a,
      "Surface data of the layout as tuples (raytatouille.layout.surfaces builds SurfaceLayout).");

  m.def(
      "layout_sag",
      [](const CompiledSystem& s, std::uint32_t surface, const RealIn& x, const RealIn& y) {
        std::vector<double> z(x.shape(0));
        compile::surface_sag(s, surface, {x.data(), x.shape(0)}, {y.data(), y.shape(0)}, z);
        return owned(std::move(z), {x.shape(0)});
      },
      "system"_a, "surface"_a, "x"_a, "y"_a, "Sag at local points, mm (raytatouille.layout.sag).");

  m.def(
      "layout_normal",
      [](const CompiledSystem& s, std::uint32_t surface, const RealIn& x, const RealIn& y,
         bool global) {
        std::vector<double> n(3 * x.shape(0));
        compile::surface_normal(s, surface, {x.data(), x.shape(0)}, {y.data(), y.shape(0)},
                                global ? compile::Frame::Global : compile::Frame::Local, n);
        return owned(std::move(n), {x.shape(0), 3});
      },
      "system"_a, "surface"_a, "x"_a, "y"_a, "global_frame"_a,
      "Unit normals at local points (raytatouille.layout.normal).");

  m.def(
      "layout_profile",
      [](const CompiledSystem& s, std::uint32_t surface, const RealIn& point, const RealIn& normal,
         std::size_t samples) {
        return polylines(compile::surface_profile(
            s, surface, {vec3(point, "plane point"), vec3(normal, "plane normal")}, samples));
      },
      "system"_a, "surface"_a, "point"_a, "normal"_a, "samples"_a,
      "Profile polylines in a section plane (raytatouille.layout.profile).");

  m.def(
      "layout_outlines",
      [](const CompiledSystem& s, std::uint32_t element, const RealIn& point, const RealIn& normal,
         std::size_t samples) {
        return polylines(compile::element_outlines(
            s, element, {vec3(point, "plane point"), vec3(normal, "plane normal")}, samples));
      },
      "system"_a, "element"_a, "point"_a, "normal"_a, "samples"_a,
      "Closed outlines of the glass segments (raytatouille.layout.outlines).");
}

}  // namespace rtt::py
