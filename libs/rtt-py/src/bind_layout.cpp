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

#include <array>
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
#include "rtt/math/units.hpp"

namespace nb = nanobind;
using namespace nb::literals;

namespace rtt::py {
namespace {

using compile::CompiledSystem;
using math::Vec3;

/// Dependent false for the static_assert over the compile::CompiledShape alternatives in
/// shape_data() (ADR 0041, point 6).
template <class>
inline constexpr bool kUnhandledShape = false;
template <class>
inline constexpr bool kUnhandledAperture = false;

using RealIn = nb::ndarray<const double, nb::ndim<1>, nb::c_contig, nb::device::cpu>;
template <typename T>
using Out = nb::ndarray<nb::numpy, T>;

/// NumPy array of the given shape that owns `values`.
template <typename T>
Out<T> owned(std::vector<T>&& values, std::initializer_list<std::size_t> shape) {
  auto data = std::make_unique<std::vector<T>>(std::move(values));
  T* pointer = data->data();
  // The capsule owns the copy from here on and frees it with the array.
  const nb::capsule owner = owning_capsule(std::move(data));
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

/// (shape type, curvature 1/mm, conic constant, polynomial coefficients, curvature_x 1/mm,
/// conic_constant_x, normalization_radius mm, axicon_angle_deg). The first three describe the y-z
/// section or profile, curvature_x and conic_constant_x the x-z section through the vertex (ADR
/// 0042, point 5): for the rotationally symmetric shapes they equal curvature and conic constant.
/// The coefficients are A4, A6, ... of even asphere, toroid and cylinder, A1, A2, ... of the odd
/// asphere. normalization_radius belongs to Q-con and Q-bfs, axicon_angle_deg to the axicon
/// (ADR 0043, point 4); none for every other shape.
using ShapeData = std::tuple<std::string,
                             double,
                             double,
                             std::vector<double>,
                             double,
                             double,
                             std::optional<double>,
                             std::optional<double>>;

/// shape_data() of one base shape (not a stack).
template <class Shape>
ShapeData base_shape_data(const Shape& shape_in) {
  return [](const auto& s) -> ShapeData {
    using S = std::decay_t<decltype(s)>;
    if constexpr (std::is_same_v<S, geom::Plane<double>>) {
      return {"plane", 0.0, 0.0, {}, 0.0, 0.0, std::nullopt, std::nullopt};
    } else if constexpr (std::is_same_v<S, geom::Conic<double>>) {
      const double c = s.curvature();
      const double k = s.conic_constant();
      return {"conic", c, k, {}, c, k, std::nullopt, std::nullopt};
    } else if constexpr (std::is_same_v<S, geom::EvenAsphere<double>>) {
      const auto [c, k] = s.base_conic();
      return {"even_asphere", c, k, s.coefficients(), c, k, std::nullopt, std::nullopt};
    } else if constexpr (std::is_same_v<S, geom::Biconic<double>>) {
      return {"biconic",       s.curvature_y(), s.conic_y(),  {},
              s.curvature_x(), s.conic_x(),     std::nullopt, std::nullopt};
    } else if constexpr (std::is_same_v<S, geom::Toroid<double>>) {
      // The x section through the vertex is the circle of rotation: c_R, conic 0.
      const auto [c, k] = s.profile().base_conic();
      return {
          "toroid",     c,           k, s.profile().coefficients(), s.curvature_of_rotation(), 0.0,
          std::nullopt, std::nullopt};
    } else if constexpr (std::is_same_v<S, geom::Cylinder<double>>) {
      // Flat in x: the x section through the vertex is a straight line, curvature and
      // conic constant 0.
      const auto [c, k] = s.profile().base_conic();
      return {"cylinder", c, k, s.profile().coefficients(), 0.0, 0.0, std::nullopt, std::nullopt};
    } else if constexpr (std::is_same_v<S, geom::OddAsphere<double>>) {
      // The curvature and conic of the base, not the vertex curvature c + 2 A2 (ADR 0043,
      // point 4); coefficients A1, A2, ...
      const auto [c, k] = s.base_conic();
      return {"odd_asphere", c, k, s.coefficients(), c, k, std::nullopt, std::nullopt};
    } else if constexpr (std::is_same_v<S, geom::Axicon<double>>) {
      // A cone without a base conic; the angle back in degree from the compiled rad (exact
      // up to the rounding of the two conversions).
      return {"axicon", 0.0, 0.0, {}, 0.0, 0.0, std::nullopt, math::rad_to_deg(s.angle())};
    } else if constexpr (std::is_same_v<S, geom::QCon<double>> ||
                         std::is_same_v<S, geom::QBfs<double>>) {
      // The curvature and conic of the base (for Q-bfs not the vertex curvature, ADR 0043,
      // point 4), the coefficients a_m (not the auxiliary b_m of Q-bfs), rho.
      const auto [c, k] = s.base_conic();
      return {std::is_same_v<S, geom::QCon<double>> ? "qcon" : "qbfs",
              c,
              k,
              s.coefficients(),
              c,
              k,
              s.normalization_radius(),
              std::nullopt};
    } else {
      // Every shape has its own branch; a new one must not be labelled even_asphere
      // (ADR 0041, point 6).
      static_assert(kUnhandledShape<S>, "shape_data: unhandled shape");
    }
  }(shape_in);
}

ShapeData shape_data(const compile::CompiledShape& shape) {
  return std::visit(
      [](const auto& s) -> ShapeData {
        using S = std::decay_t<decltype(s)>;
        if constexpr (std::is_same_v<S, geom::Stacked<double>>) {
          // A stack reports its base; the terms are in SurfaceLayout.terms (ADR 0051, point 6).
          return std::visit([](const auto& b) { return base_shape_data(b); }, s.base());
        } else {
          return base_shape_data(s);
        }
      },
      shape);
}

/// The terms of a stack as (type, scalars) (ADR 0051, point 6): Zernike terms and the XY
/// polynomial with normalization_radius and terms (the number of coefficients), the grid sag
/// with nx, ny, spacing_x, spacing_y, decenter_x and decenter_y. Empty for a shape
/// without terms.
nb::list terms_data(const compile::CompiledShape& shape) {
  nb::list out;
  const auto* stacked = std::get_if<geom::Stacked<double>>(&shape);
  if (stacked == nullptr) return out;
  for (const geom::TermVariant<double>& term : stacked->terms()) {
    std::visit(
        [&out](const auto& t) {
          using T = std::decay_t<decltype(t)>;
          if constexpr (std::is_same_v<T, geom::ZernikeTerm<double>>) {
            std::map<std::string, double> p{{"normalization_radius", t.normalization_radius()},
                                            {"terms", static_cast<double>(t.entries().size())}};
            out.append(nb::make_tuple(t.fringe() ? "zernike_fringe_sag" : "zernike_sag", p));
          } else if constexpr (std::is_same_v<T, geom::XyPolynomialTerm<double>>) {
            std::map<std::string, double> p{
                {"normalization_radius", t.normalization_radius()},
                {"terms", static_cast<double>(t.coefficients().size())}};
            out.append(nb::make_tuple("xy_polynomial", p));
          } else if constexpr (std::is_same_v<T, geom::GridSagTerm<double>>) {
            std::map<std::string, double> p{
                {"nx", static_cast<double>(t.nx())}, {"ny", static_cast<double>(t.ny())},
                {"spacing_x", t.spacing_x()},        {"spacing_y", t.spacing_y()},
                {"decenter_x", t.decenter_x()},      {"decenter_y", t.decenter_y()}};
            out.append(nb::make_tuple("grid_sag", p));
          } else {
            static_assert(kUnhandledShape<T>, "terms_data: unhandled term");
          }
        },
        term);
  }
  return out;
}

/// (aperture type, parameters in mm) or none.
std::optional<std::tuple<std::string, std::map<std::string, double>>> aperture_data(
    const std::optional<model::Aperture>& aperture) {
  if (!aperture) return std::nullopt;
  return std::visit(
      [](const auto& a) -> std::tuple<std::string, std::map<std::string, double>> {
        using A = std::decay_t<decltype(a)>;
        if constexpr (std::is_same_v<A, model::CircularAperture>) {
          return {"circular", {{"radius", a.radius.value}, {"inner_radius", a.inner_radius.value}}};
        } else if constexpr (std::is_same_v<A, model::RectangularAperture>) {
          return {"rectangular",
                  {{"half_width_x", a.half_width_x.value}, {"half_width_y", a.half_width_y.value}}};
        } else if constexpr (std::is_same_v<A, model::EllipticalAperture>) {
          return {"elliptical",
                  {{"semi_axis_x", a.semi_axis_x.value}, {"semi_axis_y", a.semi_axis_y.value}}};
        } else if constexpr (std::is_same_v<A, model::RegularPolygonAperture>) {
          return {"regular_polygon",
                  {{"sides", static_cast<double>(a.sides)},
                   {"circumradius", a.circumradius.value},
                   {"rotation_deg", a.rotation_deg}}};
        } else if constexpr (std::is_same_v<A, model::PolygonAperture>) {
          return {"polygon", {}};  // the vertices are in aperture_vertices (ADR 0050)
        } else {
          static_assert(kUnhandledAperture<A>, "aperture_data: unhandled aperture");
        }
      },
      *aperture);
}

/// The obscurations as a list of (type, parameters in mm and degree) in list order (ADR 0050,
/// point 3): "circular" (radius, decenter_x, decenter_y), "rectangular" (half_width_x,
/// half_width_y, decenter_x, decenter_y, rotation_deg), "elliptical" (semi_axis_x, semi_axis_y,
/// decenter_x, decenter_y, rotation_deg).
nb::list obscuration_data(const std::vector<model::Obscuration>& obscurations) {
  nb::list out;
  for (const model::Obscuration& o : obscurations) {
    std::visit(
        [&out](const auto& ob) {
          using O = std::decay_t<decltype(ob)>;
          std::map<std::string, double> p{{"decenter_x", ob.decenter_x.value},
                                          {"decenter_y", ob.decenter_y.value}};
          std::string type;
          if constexpr (std::is_same_v<O, model::CircularObscuration>) {
            type = "circular";
            p["radius"] = ob.radius.value;
          } else {
            p["rotation_deg"] = ob.rotation_deg;
            if constexpr (std::is_same_v<O, model::RectangularObscuration>) {
              type = "rectangular";
              p["half_width_x"] = ob.half_width_x.value;
              p["half_width_y"] = ob.half_width_y.value;
            } else {
              static_assert(std::is_same_v<O, model::EllipticalObscuration>,
                            "obscuration_data: unhandled obscuration");
              type = "elliptical";
              p["semi_axis_x"] = ob.semi_axis_x.value;
              p["semi_axis_y"] = ob.semi_axis_y.value;
            }
          }
          out.append(nb::make_tuple(type, p));
        },
        o);
  }
  return out;
}

/// The vertices of a polygon aperture as an (n, 2) array in mm, none for the other apertures
/// (ADR 0050, point 2).
std::optional<Out<double>> vertex_data(const std::optional<model::Aperture>& aperture) {
  if (!aperture) return std::nullopt;
  const std::vector<std::array<double, 2>> v = model::aperture_vertices(*aperture);
  if (v.empty()) return std::nullopt;
  std::vector<double> flat;
  flat.reserve(2 * v.size());
  for (const std::array<double, 2>& p : v) {
    flat.push_back(p[0]);
    flat.push_back(p[1]);
  }
  return owned(std::move(flat), {v.size(), 2});
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
              "only these have outlines.")
      .def_ro("active", &compile::CompiledElement::active,
              "False if the element or an ancestor has active false (ADR 0049): it keeps its "
              "place and geometry, but no path visits it.");

  m.def(
      "layout_node_frames",
      [](const CompiledSystem& s) {
        const auto matrix = [](const math::Isometry3& t) {
          std::vector<double> v(16, 0.0);
          for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j)
              v[(4 * static_cast<std::size_t>(i)) + static_cast<std::size_t>(j)] =
                  t.rotation()(i, j);
            v[(4 * static_cast<std::size_t>(i)) + 3] = t.translation()(i);
          }
          v[15] = 1.0;
          return owned(std::move(v), {4, 4});
        };
        nb::list out;
        for (const compile::NodeFrame& f : s.node_frames()) {
          out.append(nb::make_tuple(f.location, matrix(f.reference), matrix(f.to_global)));
        }
        return out;
      },
      "system"_a,
      "Frames of every node in pre-order as (pointer, reference, to_global), homogeneous 4 x 4 "
      "matrices (raytatouille.layout.node_frames).");

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
            for (int j = 0; j < 3; ++j)
              rotation[(3 * static_cast<std::size_t>(i)) + static_cast<std::size_t>(j)] = r(i, j);
          }
          const auto max_radius =
              std::visit([](const auto& sh) { return sh.max_radius(); }, c.shape);
          out.append(nb::make_tuple(
              c.id.str(), c.element, kind_name(c.element_kind), owned(std::move(rotation), {3, 3}),
              owned(std::vector<double>{t.x(), t.y(), t.z()}, {3}), shape_data(c.shape), max_radius,
              aperture_data(c.aperture), c.medium_front, c.medium_back, c.active,
              vertex_data(c.aperture), obscuration_data(c.obscurations), terms_data(c.shape)));
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
