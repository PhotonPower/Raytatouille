#include "rtt/optim/generators.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#include "rtt/analysis/opd.hpp"
#include "rtt/math/types.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/sequential.hpp"
#include "rtt/trace/sources.hpp"

namespace rtt::optim {

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

/// What both generator types have.
struct Spec {
  const std::string* path = nullptr;
  const std::optional<std::vector<std::uint16_t>>* fields = nullptr;
  const std::optional<std::vector<std::uint16_t>>* wavelengths = nullptr;
  int rings = 0;
  int arms = 0;
  double weight = 0.0;
};

// get_if, no std::visit returning a reference (GCC 13 -Wdangling-reference, see merit.cpp).
Spec spec_of(const model::Generator& g) {
  if (const auto* s = std::get_if<model::SpotGenerator>(&g)) {
    return {&s->path, &s->fields, &s->wavelengths, s->rings, s->arms, s->weight};
  }
  const auto& w = std::get<model::WavefrontGenerator>(g);
  return {&w.path, &w.fields, &w.wavelengths, w.rings, w.arms, w.weight};
}

trace::GaussPupil gauss_of(const Spec& spec) {
  if (spec.rings < 1 || spec.arms < 1) {
    throw std::invalid_argument("optimize: a generator needs rings >= 1 and arms >= 1");
  }
  return {spec.rings, spec.arms};
}

/// The chosen indices: the list of the generator, or 0 .. count - 1.
std::vector<std::uint16_t> chosen(const std::optional<std::vector<std::uint16_t>>& list,
                                  std::size_t count,
                                  const char* what) {
  std::vector<std::uint16_t> out;
  if (list) {
    for (const std::uint16_t i : *list) {
      if (i >= count) {
        throw std::invalid_argument(std::string("optimize: generator ") + what + " index " +
                                    std::to_string(i) + " out of range");
      }
    }
    out = *list;
  } else {
    for (std::size_t i = 0; i < count; ++i) out.push_back(static_cast<std::uint16_t>(i));
  }
  return out;
}

/// Model weights of the chosen entries, normalised to sum 1 over the choice (addendum #168).
std::vector<double> normalised(const std::vector<double>& weights, const char* what) {
  double total = 0.0;
  for (const double w : weights) total += w;
  if (!(total > 0.0)) {
    throw std::invalid_argument(std::string("optimize: the weights of the chosen ") + what +
                                " of a generator sum to 0");
  }
  std::vector<double> out;
  out.reserve(weights.size());
  for (const double w : weights) out.push_back(w / total);
  return out;
}

compile::PathId path_of(const compile::CompiledSystem& cs, const std::string& name) {
  const std::optional<compile::PathId> id = cs.find_path(name);
  if (!id) throw std::invalid_argument("optimize: unknown path '" + name + "' of a generator");
  return *id;
}

/// Traces the rays of one field and wavelength (real aiming) through the whole path.
trace::RayBatch traced(const compile::CompiledSystem& cs,
                       compile::PathId path,
                       std::uint16_t field,
                       std::uint16_t wavelength,
                       const trace::PupilSampling& sampling) {
  trace::RayBatch rays =
      trace::make_rays(cs, path, std::span<const std::uint16_t>(&field, 1), wavelength, sampling);
  static_cast<void>(trace::SequentialTracer().trace(cs, path, rays));
  return rays;
}

/// A point on the image surface in its local coordinates, mm; `arrived` false if there is none.
/// Plain values instead of std::optional<Vec3>: GCC at -O2 reports the Eigen payload of
/// an optional as maybe uninitialized.
struct Point {
  bool arrived = false;
  double x = 0.0;
  double y = 0.0;
};

/// Hit point of ray i in the local coordinates of the image surface; not arrived if the ray did
/// not end Alive there (as analysis::spot: the image surface is that of the last event).
Point local_hit(const compile::CompiledSystem& cs,
                const trace::RayBatch& rays,
                std::size_t i,
                std::uint32_t image) {
  if (rays.status()[i] != trace::RayStatus::Alive || rays.last_surface()[i] != image) return {};
  const math::Vec3 p = cs.surfaces()[image].to_local.apply_point(
      math::Vec3(rays.pos_x()[i], rays.pos_y()[i], rays.pos_z()[i]));
  return {true, p.x(), p.y()};
}

/// Choice and weights of one generator in one compiled system.
struct Expansion {
  compile::PathId path;
  std::vector<std::uint16_t> fields;
  std::vector<std::uint16_t> wavelengths;
  std::vector<double> wf;  ///< W_f, normalised
  std::vector<double> wl;  ///< W_l, normalised
  trace::GaussPupil gauss;
  std::vector<double> q;  ///< quadrature weights, pupil_points order
  double weight = 0.0;
};

Expansion expand(const model::Generator& g, const compile::CompiledSystem& cs) {
  const Spec spec = spec_of(g);
  Expansion e;
  e.path = path_of(cs, *spec.path);
  e.fields = chosen(*spec.fields, cs.fields().points.size(), "field");
  e.wavelengths = chosen(*spec.wavelengths, cs.wavelengths_um().size(), "wavelength");
  std::vector<double> w;
  w.reserve(e.fields.size());
  for (const std::uint16_t f : e.fields) w.push_back(cs.fields().points[f].weight);
  e.wf = normalised(w, "fields");
  w.clear();
  w.reserve(e.wavelengths.size());
  for (const std::uint16_t l : e.wavelengths) w.push_back(cs.wavelength_weights()[l]);
  e.wl = normalised(w, "wavelengths");
  e.gauss = gauss_of(spec);
  e.q = trace::gauss_pupil_weights(e.gauss);
  e.weight = spec.weight;
  return e;
}

/// rms_spot (ADR 0030, point 4; addendum #168). Two passes per field: hits of all chosen
/// wavelengths, then the reference and the residuals.
void spot_residuals(const model::SpotGenerator& g,
                    const compile::CompiledSystem& cs,
                    const Expansion& e,
                    std::span<double> out,
                    GeneratorStats& stats) {
  const std::uint32_t image = cs.path(e.path).events.back().surface;
  const std::size_t n = e.q.size();
  std::size_t o = 0;
  for (std::size_t f = 0; f < e.fields.size(); ++f) {
    // Hits in the order wavelengths, then pupil points.
    std::vector<Point> hits;
    hits.reserve(e.wavelengths.size() * n);
    for (const std::uint16_t l : e.wavelengths) {
      const trace::RayBatch rays = traced(cs, e.path, e.fields[f], l, e.gauss);
      for (std::size_t k = 0; k < n; ++k) hits.push_back(local_hit(cs, rays, k, image));
    }
    // Reference: centroid weighted with W_l q_k over the arrived rays, or the chief ray of the
    // reference wavelength (as analysis::spot).
    Point reference;
    std::string why;
    if (g.reference == model::SpotReference::Centroid) {
      double sw = 0.0;
      double cx = 0.0;
      double cy = 0.0;
      for (std::size_t l = 0; l < e.wavelengths.size(); ++l) {
        for (std::size_t k = 0; k < n; ++k) {
          const Point& h = hits[l * n + k];
          if (!h.arrived) continue;
          const double c = e.wl[l] * e.q[k];
          sw += c;
          cx += c * h.x;
          cy += c * h.y;
        }
      }
      if (sw > 0.0) {
        reference = {true, cx / sw, cy / sw};
      } else {
        why = "no ray of field " + std::to_string(e.fields[f]) + " reaches the image surface";
      }
    } else {
      const trace::RayBatch chief = traced(cs, e.path, e.fields[f], cs.reference_wavelength(),
                                           trace::SinglePupilPoint{0.0, 0.0});
      reference = local_hit(cs, chief, 0, image);
      if (!reference.arrived) {
        why = "the chief ray of field " + std::to_string(e.fields[f]) +
              " does not reach the image surface";
      }
      const bool any =
          std::any_of(hits.begin(), hits.end(), [](const Point& h) { return h.arrived; });
      if (reference.arrived && !any) {
        reference = {};
        why = "no ray of field " + std::to_string(e.fields[f]) + " reaches the image surface";
      }
    }
    for (std::size_t l = 0; l < e.wavelengths.size(); ++l) {
      for (std::size_t k = 0; k < n; ++k) {
        const Point& h = hits[l * n + k];
        ++stats.rays_launched;
        if (!reference.arrived) {
          if (!h.arrived) ++stats.rays_lost;
          out[o++] = kNaN;
          out[o++] = kNaN;
          continue;
        }
        if (!h.arrived) {
          ++stats.rays_lost;
          out[o++] = 0.0;
          out[o++] = 0.0;
          continue;
        }
        const double c = e.wf[f] * e.wl[l] * e.q[k];
        const double dx = h.x - reference.x;
        const double dy = h.y - reference.y;
        const double s = std::sqrt(e.weight * c);
        out[o++] = s * dx;
        out[o++] = s * dy;
        stats.mean_square += c * (dx * dx + dy * dy);
      }
    }
    if (!reference.arrived) {
      stats.mean_square = kNaN;
      if (stats.undefined.empty()) stats.undefined = why;
    }
  }
}

/// rms_wavefront (ADR 0030, point 4; addendum #168): W from analysis::opd_points, piston W_mean
/// per field and wavelength.
void wavefront_residuals(const compile::CompiledSystem& cs,
                         const Expansion& e,
                         std::span<double> out,
                         GeneratorStats& stats) {
  const std::size_t n = e.q.size();
  std::size_t o = 0;
  for (std::size_t f = 0; f < e.fields.size(); ++f) {
    for (std::size_t l = 0; l < e.wavelengths.size(); ++l) {
      const analysis::OpdPupilPoints p =
          analysis::opd_points(cs, e.path, e.fields[f], e.wavelengths[l], e.gauss);
      double sq = 0.0;
      double sw = 0.0;
      for (std::size_t k = 0; k < n; ++k) {
        if (p.points[k].status != trace::RayStatus::Alive) continue;
        sq += e.q[k];
        sw += e.q[k] * p.points[k].w;
      }
      const bool defined = sq > 0.0;
      const double mean = defined ? sw / sq : kNaN;
      for (std::size_t k = 0; k < n; ++k) {
        ++stats.rays_launched;
        const bool alive = p.points[k].status == trace::RayStatus::Alive;
        if (!alive) ++stats.rays_lost;
        if (!defined) {
          out[o++] = kNaN;
          continue;
        }
        if (!alive) {
          out[o++] = 0.0;
          continue;
        }
        const double c = e.wf[f] * e.wl[l] * e.q[k];
        const double d = p.points[k].w - mean;
        out[o++] = std::sqrt(e.weight * c) * d;
        stats.mean_square += c * d * d;
      }
      if (!defined) {
        stats.mean_square = kNaN;
        if (stats.undefined.empty()) {
          stats.undefined = "no ray of field " + std::to_string(e.fields[f]) + " at wavelength " +
                            std::to_string(e.wavelengths[l]) +
                            " reaches the image surface and the reference sphere";
        }
      }
    }
  }
}

}  // namespace

std::size_t generator_size(const model::Generator& generator,
                           const compile::CompiledSystem& system) {
  const Spec spec = spec_of(generator);
  const trace::GaussPupil gauss = gauss_of(spec);
  const std::size_t fields = *spec.fields ? (*spec.fields)->size() : system.fields().points.size();
  const std::size_t wavelengths =
      *spec.wavelengths ? (*spec.wavelengths)->size() : system.wavelengths_um().size();
  const std::size_t per_ray = std::holds_alternative<model::SpotGenerator>(generator) ? 2 : 1;
  return fields * wavelengths * static_cast<std::size_t>(gauss.rings) *
         static_cast<std::size_t>(gauss.arms) * per_ray;
}

void generator_residuals(const model::Generator& generator,
                         const compile::CompiledSystem& system,
                         std::span<double> out,
                         GeneratorStats& stats) {
  if (out.size() != generator_size(generator, system)) {
    throw std::invalid_argument("optimize: wrong number of generator residuals");
  }
  const Expansion e = expand(generator, system);
  stats = GeneratorStats{};
  if (const auto* spot = std::get_if<model::SpotGenerator>(&generator)) {
    spot_residuals(*spot, system, e, out, stats);
  } else {
    wavefront_residuals(system, e, out, stats);
  }
}

}  // namespace rtt::optim
