#include "rtt/analysis/irradiance.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include "rtt/math/units.hpp"
#include "rtt/model/path.hpp"
#include "rtt/model/surface.hpp"

namespace rtt::analysis {

namespace {

[[noreturn]] void fail(const std::string& what) {
  throw std::invalid_argument("irradiance_map: " + what);
}

/// The bounding rectangle of an aperture in local coordinates (ADR 0048, point 3).
IrradianceExtent bounding(const model::Aperture& aperture) {
  return std::visit(
      [](const auto& a) -> IrradianceExtent {
        using A = std::decay_t<decltype(a)>;
        if constexpr (std::is_same_v<A, model::CircularAperture>) {
          return {-a.radius, a.radius, -a.radius, a.radius};
        } else if constexpr (std::is_same_v<A, model::RectangularAperture>) {
          return {-a.half_width_x, a.half_width_x, -a.half_width_y, a.half_width_y};
        } else {
          return {-a.semi_axis_x, a.semi_axis_x, -a.semi_axis_y, a.semi_axis_y};
        }
      },
      aperture);
}

/// Grid edges x_i = lo + i (hi - lo) / n, the last one exactly hi.
std::vector<double> edges(double lo, double hi, std::size_t n) {
  std::vector<double> e(n + 1);
  const double step = (hi - lo) / static_cast<double>(n);
  for (std::size_t i = 0; i < n; ++i) e[i] = lo + static_cast<double>(i) * step;
  e[n] = hi;
  return e;
}

/// Cell index of v for the edges e (n + 1 values, e[n] = hi), or nullopt outside [e[0], e[n]]:
/// the cell [e[i], e[i + 1]), the last one including e[n]. The estimate from the formula may be
/// off by one at an edge through rounding, so it is checked against the edges themselves.
std::optional<std::size_t> cell(double v, const std::vector<double>& e) {
  const std::size_t n = e.size() - 1;
  if (!(v >= e.front() && v <= e.back())) return std::nullopt;  // also NaN
  const double t = (v - e.front()) / (e.back() - e.front()) * static_cast<double>(n);
  std::size_t i = t < static_cast<double>(n) ? static_cast<std::size_t>(t) : n - 1;  // t >= 0
  while (i > 0 && v < e[i]) --i;
  while (i + 1 < n && v >= e[i + 1]) ++i;
  return i;
}

/// The empty map on surface `index`: grid, edges and the surface fields; the counts are added by
/// the accumulator below.
class Accumulator {
 public:
  Accumulator(const compile::CompiledSystem& system,
              std::uint32_t index,
              const IrradianceOptions& options)
      : surface_(system.surfaces()[index]) {
    if (options.nx == 0 || options.ny == 0) fail("nx and ny must be >= 1");
    IrradianceExtent e{};
    if (options.extent) {
      e = *options.extent;
    } else if (surface_.aperture) {
      e = bounding(*surface_.aperture);
    } else {
      fail("surface " + surface_.id.str() + " has no aperture; give an extent");
    }
    // Widths and the cell area must be finite and positive too: a finite extent of +-1e308 has an
    // infinite width (ADR 0048, point 3).
    const double width = e.x_max - e.x_min;
    const double height = e.y_max - e.y_min;
    const double area =
        width / static_cast<double>(options.nx) * (height / static_cast<double>(options.ny));
    if (!std::isfinite(e.x_min) || !std::isfinite(e.x_max) || !std::isfinite(e.y_min) ||
        !std::isfinite(e.y_max) || !(e.x_min < e.x_max) || !(e.y_min < e.y_max) ||
        !std::isfinite(width) || !std::isfinite(height) || !std::isfinite(area) || !(area > 0.0)) {
      fail(
          "the extent must be finite with x_min < x_max and y_min < y_max, and give finite "
          "cells of positive area");
    }
    area_ = area;
    map_.nx = options.nx;
    map_.ny = options.ny;
    map_.x_edges = edges(e.x_min, e.x_max, options.nx);
    map_.y_edges = edges(e.y_min, e.y_max, options.ny);
    map_.values.assign(options.nx * options.ny, 0.0);
    map_.surface = index;
    map_.surface_id = surface_.id;
    map_.pointer = surface_.location;
    map_.to_global = surface_.to_global;
  }

  /// A hit inside the aperture at the global point `p` with incident weight w.
  void hit(const math::Vec3& p, double w) {
    const math::Vec3 local = surface_.to_local.apply_point(p);
    const std::optional<std::size_t> ix = cell(local.x(), map_.x_edges);
    const std::optional<std::size_t> iy = cell(local.y(), map_.y_edges);
    if (!ix || !iy) {
      add(map_.outside_grid, w);
      return;
    }
    map_.values[*iy * map_.nx + *ix] += w;  // the sum; divided by the cell area in finish()
    ++map_.rays_counted;
    map_.weight_counted += w;
  }

  void vignetted(double w) { add(map_.vignetted, w); }
  void not_reached(double w) { add(map_.not_reached, w); }

  IrradianceMap finish(std::vector<std::size_t> events) {
    // Projected cell area dx dy (ADR 0048, point 4).
    for (double& v : map_.values) v /= area_;
    map_.events = std::move(events);
    return std::move(map_);
  }

 private:
  static void add(IrradianceLosses& l, double w) {
    ++l.count;
    l.weight += w;
  }

  const compile::CompiledSurface& surface_;
  double area_ = 0.0;
  IrradianceMap map_;
};

/// True for a status whose state carries the hit point inside the aperture (ray_paths.hpp):
/// the ray reached the surface and its incident power counts.
bool hit_inside(trace::RayStatus s) {
  switch (s) {
    case trace::RayStatus::Alive:
    case trace::RayStatus::Tir:
    case trace::RayStatus::Absorbed:
    case trace::RayStatus::EventImpossible:
    case trace::RayStatus::Evanescent:
      return true;
    case trace::RayStatus::Missed:
    case trace::RayStatus::Vignetted:
    case trace::RayStatus::NoConvergence:
      return false;
  }
  return false;
}

}  // namespace

IrradianceMap irradiance_map(const compile::CompiledSystem& system,
                             std::size_t path,
                             const trace::RayPaths& paths,
                             const model::SurfaceId& surface,
                             const IrradianceOptions& options) {
  if (path >= system.paths().size()) fail("no path " + std::to_string(path));
  const std::vector<compile::CompiledEvent>& path_events = system.paths()[path].events;
  const std::size_t n_events = paths.event_surfaces.size();
  bool same_path = paths.slots == n_events + 1 && n_events == path_events.size() &&
                   paths.count.size() == paths.ray_count() &&
                   paths.weight.size() == paths.ray_count() * paths.slots &&
                   paths.status.size() == paths.ray_count() * paths.slots &&
                   paths.position.size() == 3 * paths.ray_count() * paths.slots;
  for (std::size_t k = 0; same_path && k < n_events; ++k) {
    same_path = paths.event_surfaces[k] == path_events[k].surface;
  }
  if (!same_path) fail("the paths were not traced on path " + system.paths()[path].name);
  const std::optional<std::uint32_t> found = system.find_surface(surface);
  if (!found) fail("unknown surface " + surface.str());
  const std::uint32_t index = *found;
  std::vector<std::size_t> events;
  if (options.event) {
    const std::size_t e = *options.event;
    if (e >= n_events || paths.event_surfaces[e] != index) {
      fail("event " + std::to_string(e) + " of the path is not at surface " + surface.str());
    }
    events.push_back(e);
  } else {
    for (std::size_t k = 0; k < n_events; ++k) {
      if (paths.event_surfaces[k] == index) events.push_back(k);
    }
    if (events.empty()) fail("surface " + surface.str() + " is not an event of the path");
  }

  // Volume absorption on the way to each counted event (ADR 0048, point 1): kappa of the medium
  // before it at the wavelength of the rays; the trace applies it before the event (move_to_hit,
  // Byrnes, arXiv:1603.02720v5, Eqs. (1), (2): power ~ exp(-4 pi kappa z / lambda_vac)).
  std::vector<double> kappa(n_events, 0.0);
  double lambda_mm = 0.0;
  for (const std::size_t k : events) {
    const std::vector<math::Complex>& n = system.media()[path_events[k].medium_before].index;
    const bool absorbs =
        std::any_of(n.begin(), n.end(), [](const math::Complex& c) { return c.imag() != 0.0; });
    if (!absorbs) continue;
    if (!options.wavelength || *options.wavelength >= n.size()) {
      fail("the medium before event " + std::to_string(k) +
           " absorbs; give the wavelength index "
           "of the rays for the absorption on the last path");
    }
    kappa[k] = n[*options.wavelength].imag();
    lambda_mm = math::um_to_mm(system.wavelengths_um()[*options.wavelength]);
  }

  Accumulator acc(system, index, options);
  const std::size_t S = paths.slots;
  // Rays in batch order, events in path order per ray (ADR 0048, point 6).
  for (std::size_t r = 0; r < paths.ray_count(); ++r) {
    const std::size_t count = paths.count[r];
    for (const std::size_t k : events) {
      if (count < k + 2) {
        // Lost before event k (or not alive at the start): the weight before the event where it
        // was lost (slot count - 2), or its start weight.
        acc.not_reached(paths.weight[r * S + (count >= 2 ? count - 2 : 0)]);
        continue;
      }
      double w = paths.weight[r * S + k];  // the state after event k - 1
      if (kappa[k] != 0.0) {
        const double t = (paths.position_at(r, k + 1) - paths.position_at(r, k)).norm();
        w *= std::exp(-4.0 * std::numbers::pi * kappa[k] * t / lambda_mm);
      }
      const trace::RayStatus status = paths.status[r * S + k + 1];
      if (hit_inside(status)) {
        acc.hit(paths.position_at(r, k + 1), w);
      } else if (status == trace::RayStatus::Vignetted) {
        acc.vignetted(w);
      } else {
        acc.not_reached(w);
      }
    }
  }
  return acc.finish(std::move(events));
}

IrradianceMap irradiance_map(const compile::CompiledSystem& system,
                             std::size_t path,
                             const trace::RayBatch& rays,
                             const IrradianceOptions& options) {
  if (path >= system.paths().size()) fail("no path " + std::to_string(path));
  const std::vector<compile::CompiledEvent>& events = system.paths()[path].events;
  if (events.empty()) fail("path " + system.paths()[path].name + " has no event");
  const compile::CompiledEvent& last = events.back();
  const compile::CompiledSurface& s = system.surfaces()[last.surface];
  // The batch holds the weight after the last event; it is the incident one only if that event
  // cannot change it (ADR 0048, point 2).
  const bool lossless_interaction =
      std::holds_alternative<model::Fresnel>(s.interaction) ||
      std::holds_alternative<model::IdealAntiReflection>(s.interaction);
  if (last.kind != model::EventKind::Transmit || last.order != 0 || !lossless_interaction ||
      last.medium_before != last.medium_after || s.ideal_lens || s.diffraction_efficiency) {
    fail("the last event of path " + system.paths()[path].name + " at surface " + s.id.str() +
         " may change the weight; use the RayPaths overload for the incident weight");
  }

  Accumulator acc(system, last.surface, options);
  for (std::size_t r = 0; r < rays.size(); ++r) {
    const double w = rays.weight()[r];
    const trace::RayStatus status = rays.status()[r];
    if (status == trace::RayStatus::Alive) {
      acc.hit({rays.pos_x()[r], rays.pos_y()[r], rays.pos_z()[r]}, w);
    } else if (status == trace::RayStatus::Vignetted && rays.last_surface()[r] == last.surface) {
      acc.vignetted(w);
    } else {
      acc.not_reached(w);
    }
  }
  return acc.finish({events.size() - 1});
}

}  // namespace rtt::analysis
