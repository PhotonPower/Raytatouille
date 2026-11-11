// Irradiance map on a surface from traced rays (#379, ADR 0048): rtt/analysis/irradiance.hpp.
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/variant.h>
#include <nanobind/stl/vector.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "bindings.hpp"
#include "rtt/analysis/irradiance.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/math/isometry.hpp"
#include "rtt/model/ids.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/ray_paths.hpp"

namespace nb = nanobind;
using namespace nb::literals;

namespace rtt::py {
namespace {

using analysis::IrradianceExtent;
using analysis::IrradianceLosses;
using analysis::IrradianceMap;
using analysis::IrradianceOptions;

template <typename T, std::size_t D>
using Array = nb::ndarray<nb::numpy, const T, nb::ndim<D>, nb::c_contig>;

/// Read-only NumPy array of the given shape that owns a copy of `values` (row-major).
template <typename T, std::size_t D>
Array<T, D> owned(std::vector<T> values, const std::array<std::size_t, D>& shape) {
  auto data = std::make_unique<std::vector<T>>(std::move(values));
  const nb::capsule owner(data.get(), [](void* q) noexcept {
    const std::unique_ptr<std::vector<T>> freed(static_cast<std::vector<T>*>(q));
  });
  const std::vector<T>* by_capsule = data.release();  // freed by the capsule
  return {by_capsule->data(), D, shape.data(), owner};
}

/// Homogeneous 4 x 4 matrix of a pose, row-major: [R t; 0 0 0 1].
Array<double, 2> homogeneous(const math::Isometry3& pose) {
  const math::Mat3& r = pose.rotation();
  const math::Vec3 t = pose.apply_point(math::Vec3::Zero());
  std::vector<double> m(16, 0.0);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) m[static_cast<std::size_t>(4 * i + j)] = r(i, j);
    m[static_cast<std::size_t>(4 * i + 3)] = t[i];
  }
  m[15] = 1.0;
  return owned<double, 2>(std::move(m), {4, 4});
}

IrradianceOptions options(std::size_t nx,
                          std::size_t ny,
                          const std::optional<std::tuple<double, double, double, double>>& extent,
                          std::optional<std::size_t> event) {
  IrradianceOptions o;
  o.nx = nx;
  o.ny = ny;
  if (extent) {
    const auto& [x_min, x_max, y_min, y_max] = *extent;
    o.extent = IrradianceExtent{x_min, x_max, y_min, y_max};
  }
  o.event = event;
  return o;
}

}  // namespace

void bind_irradiance(nb::module_& m) {
  nb::class_<IrradianceLosses>(m, "IrradianceLosses",
                               "Rays of an irradiance map that did not count (#379).")
      .def_ro("count", &IrradianceLosses::count, "Number of rays (per counted event).")
      .def_ro("weight", &IrradianceLosses::weight, "Their summed incident weight.");

  nb::class_<IrradianceMap>(
      m, "IrradianceMap",
      "Irradiance on a surface from traced rays (#379, ADR 0048): the incident weight of the "
      "hits per cell of the local x-y plane, divided by the projected cell area.")
      .def_ro("nx", &IrradianceMap::nx, "Cells in local x.")
      .def_ro("ny", &IrradianceMap::ny, "Cells in local y.")
      .def_prop_ro(
          "values", [](const IrradianceMap& r) { return owned<double, 2>(r.values, {r.ny, r.nx}); },
          nb::rv_policy::reference,
          "(ny, nx) incident weight per mm^2; row iy, column ix (read-only copy).")
      .def_prop_ro(
          "x_edges", [](const IrradianceMap& r) { return owned<double, 1>(r.x_edges, {r.nx + 1}); },
          nb::rv_policy::reference, "nx + 1 cell edges in local x, mm (read-only copy).")
      .def_prop_ro(
          "y_edges", [](const IrradianceMap& r) { return owned<double, 1>(r.y_edges, {r.ny + 1}); },
          nb::rv_policy::reference, "ny + 1 cell edges in local y, mm (read-only copy).")
      .def_ro("surface", &IrradianceMap::surface, "Surface index in the compiled system.")
      .def_prop_ro(
          "surface_id", [](const IrradianceMap& r) { return r.surface_id.str(); },
          "Stable surface id, e.g. \"IMG\".")
      .def_ro("pointer", &IrradianceMap::pointer, "JSON pointer of the surface in the file.")
      .def_prop_ro(
          "to_global", [](const IrradianceMap& r) { return homogeneous(r.to_global); },
          nb::rv_policy::reference,
          "Local surface coordinates -> global, homogeneous 4 x 4 (read-only copy).")
      .def_ro("events", &IrradianceMap::events, "Path events that were counted.")
      .def_ro("rays_counted", &IrradianceMap::rays_counted,
              "Hits inside the aperture and the grid, over all counted events.")
      .def_ro("weight_counted", &IrradianceMap::weight_counted,
              "Their incident weight; values.sum() * dx * dy up to rounding.")
      .def_ro("outside_grid", &IrradianceMap::outside_grid,
              "Hits inside the aperture but outside the grid.")
      .def_ro("vignetted", &IrradianceMap::vignetted, "Hits outside the aperture.")
      .def_ro("not_reached", &IrradianceMap::not_reached,
              "Rays lost before the surface, or Missed / NoConvergence at it.");

  m.def(
      "irradiance_map_paths",
      [](const compile::CompiledSystem& system, const PathArg& path, const trace::RayPaths& paths,
         const std::string& surface, std::size_t nx, std::size_t ny,
         const std::optional<std::tuple<double, double, double, double>>& extent,
         std::optional<std::size_t> event, std::optional<std::uint16_t> wavelength) {
        IrradianceOptions o = options(nx, ny, extent, event);
        o.wavelength = wavelength;
        return analysis::irradiance_map(system, path_id(system, path).index, paths,
                                        model::SurfaceId{surface}, o);
      },
      "system"_a, "path"_a, "paths"_a, "surface"_a, "nx"_a, "ny"_a, "extent"_a.none(),
      "event"_a.none(), "wavelength"_a.none(),
      "Irradiance map on `surface` from recorded paths (raytatouille.analysis.irradiance_map).");

  m.def(
      "irradiance_map_rays",
      [](const compile::CompiledSystem& system, const PathArg& path, const trace::RayBatch& rays,
         std::size_t nx, std::size_t ny,
         const std::optional<std::tuple<double, double, double, double>>& extent) {
        return analysis::irradiance_map(system, path_id(system, path).index, rays,
                                        options(nx, ny, extent, std::nullopt));
      },
      "system"_a, "path"_a, "rays"_a, "nx"_a, "ny"_a, "extent"_a.none(),
      "Irradiance map on the last surface of `path` from the batch after the trace "
      "(raytatouille.analysis.irradiance_map).");
}

}  // namespace rtt::py
