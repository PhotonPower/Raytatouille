#include "rtt/nonseq/source.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <numbers>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <variant>

#include "rtt/math/philox.hpp"

namespace rtt::nonseq {

namespace {

constexpr std::uint64_t kMaxRays = std::uint64_t{1} << 40;

/// Smallest n with n^k >= value, for k = 2 or 4 and 1 <= value <= 2^40 (exact integer search
/// around the floating-point root).
std::uint64_t ceil_root(std::uint64_t value, int k) noexcept {
  const auto power = [k](std::uint64_t n) { return k == 2 ? n * n : n * n * n * n; };
  auto n = static_cast<std::uint64_t>(std::pow(static_cast<double>(value), 1.0 / k));
  n = std::max<std::uint64_t>(n, 1);
  while (n > 1 && power(n - 1) >= value) --n;
  while (power(n) < value) ++n;
  return n;
}

bool is_grid_4d(const Source& s) noexcept {
  return s.sampling == Sampling::Grid && std::holds_alternative<LambertEmitter>(s.emitter);
}

/// Point on the shape for (xi_1, xi_2): disc PBRT 4ed, §A.5.1; rectangle linear.
math::Vec3 on_shape(const EmitterShape& shape, double x1, double x2) noexcept {
  if (const auto* d = std::get_if<Disc>(&shape)) {
    const double r = d->radius * std::sqrt(x1);
    const double phi = 2.0 * std::numbers::pi * x2;
    return {r * std::cos(phi), r * std::sin(phi), 0.0};
  }
  const auto& rect = std::get<Rectangle>(shape);
  return {(2.0 * x1 - 1.0) * rect.half_x, (2.0 * x2 - 1.0) * rect.half_y, 0.0};
}

void check_shape(const EmitterShape& shape) {
  const auto positive = [](double v) { return std::isfinite(v) && v > 0.0; };
  const bool ok = std::visit(
      [&](const auto& s) {
        using S = std::remove_cvref_t<decltype(s)>;
        if constexpr (std::is_same_v<S, Disc>) {
          return positive(s.radius);
        } else {
          return positive(s.half_x) && positive(s.half_y);
        }
      },
      shape);
  if (!ok) throw std::invalid_argument("nonseq: source shape dimensions must be finite and > 0");
}

}  // namespace

std::uint64_t source_ray_count(const Source& source) noexcept {
  if (source.sampling == Sampling::Random || source.rays == 0) return source.rays;
  const std::uint64_t n = ceil_root(source.rays, is_grid_4d(source) ? 4 : 2);
  return is_grid_4d(source) ? n * n * n * n : n * n;
}

SourceRay source_ray(const Source& source,
                     std::uint32_t source_number,
                     std::uint64_t index,
                     std::uint64_t seed) noexcept {
  // (xi_1, xi_2, xi_3, xi_4) of ADR 0057, Nachtrag PR 4, f.
  std::array<double, 4> xi{};
  if (source.sampling == Sampling::Random) {
    const math::PhiloxKey key{static_cast<std::uint32_t>(seed),
                              static_cast<std::uint32_t>(seed >> 32)};
    const std::uint32_t draws = std::holds_alternative<LambertEmitter>(source.emitter) ? 2 : 1;
    for (std::uint32_t j = 0; j < draws; ++j) {
      const math::PhiloxCounter counter{source_number << 16, static_cast<std::uint32_t>(index),
                                        static_cast<std::uint32_t>(index >> 32), j};
      const auto w = math::philox4x32_10(counter, key);
      xi[2 * j] = math::uniform53(w[0], w[1]);
      xi[2 * j + 1] = math::uniform53(w[2], w[3]);
    }
  } else {
    const bool four = is_grid_4d(source);
    const std::uint64_t n = ceil_root(source.rays, four ? 4 : 2);
    const double dn = static_cast<double>(n);
    const auto cell = [&](std::uint64_t c, std::size_t at) {
      xi[at] = (static_cast<double>(c % n) + 0.5) / dn;
      xi[at + 1] = (static_cast<double>(c / n) + 0.5) / dn;
    };
    if (four) {
      cell(index / (n * n), 0);  // position cell p
      cell(index % (n * n), 2);  // direction cell q
    } else {
      cell(index, 0);
    }
  }

  math::Vec3 pos = math::Vec3::Zero();
  math::Vec3 dir = math::Vec3::UnitZ();
  std::visit(
      [&](const auto& e) {
        using E = std::remove_cvref_t<decltype(e)>;
        if constexpr (std::is_same_v<E, PointEmitter>) {
          // PBRT 4ed, §A.5.4: cos theta = (1 - xi) + xi cos theta_max, with 1 - cos theta_max
          // = 2 sin^2(theta_max / 2) against cancellation (ADR 0057, point 6).
          const double s = std::sin(0.5 * e.max_angle);
          const double u = 2.0 * xi[0] * s * s;  // 1 - cos theta
          const double sin_t = std::sqrt(std::max(0.0, u * (2.0 - u)));
          const double phi = 2.0 * std::numbers::pi * xi[1];
          dir = {sin_t * std::cos(phi), sin_t * std::sin(phi), 1.0 - u};
        } else if constexpr (std::is_same_v<E, CollimatedEmitter>) {
          pos = on_shape(e.shape, xi[0], xi[1]);
        } else {
          pos = on_shape(e.shape, xi[0], xi[1]);
          // Malley (PBRT 4ed, §A.5.3): uniform point on the disc of radius sin theta_max
          // (§A.5.1), lifted to the hemisphere.
          const double rho = std::sin(e.max_angle) * std::sqrt(xi[2]);
          const double phi = 2.0 * std::numbers::pi * xi[3];
          dir = {rho * std::cos(phi), rho * std::sin(phi),
                 std::sqrt(std::max(0.0, 1.0 - rho * rho))};
        }
      },
      source.emitter);
  return {source.pose.apply_point(pos), source.pose.apply_vector(dir)};
}

void check_source(const Source& source, std::size_t wavelengths) {
  std::visit(
      [](const auto& e) {
        using E = std::remove_cvref_t<decltype(e)>;
        if constexpr (std::is_same_v<E, PointEmitter>) {
          if (!(e.max_angle > 0.0 && e.max_angle <= std::numbers::pi))
            throw std::invalid_argument("nonseq: point source max_angle must be in (0, pi]");
        } else if constexpr (std::is_same_v<E, CollimatedEmitter>) {
          check_shape(e.shape);
        } else {
          check_shape(e.shape);
          if (!(e.max_angle > 0.0 && e.max_angle <= std::numbers::pi / 2.0))
            throw std::invalid_argument("nonseq: Lambert source max_angle must be in (0, pi/2]");
        }
      },
      source.emitter);
  if (!(std::isfinite(source.power) && source.power > 0.0))
    throw std::invalid_argument("nonseq: source power must be finite and > 0");
  if (source.rays == 0 || source.rays > kMaxRays)
    throw std::invalid_argument("nonseq: source rays must be in [1, 2^40]");
  if (!source.pose.rotation().allFinite() || !source.pose.translation().allFinite())
    throw std::invalid_argument("nonseq: source pose not finite");
  for (std::size_t i = 0; i < source.wavelengths.size(); ++i) {
    if (source.wavelengths[i] >= wavelengths)
      throw std::invalid_argument("nonseq: source wavelength out of range");
    for (std::size_t j = 0; j < i; ++j) {
      if (source.wavelengths[j] == source.wavelengths[i])
        throw std::invalid_argument("nonseq: source wavelength listed twice");
    }
  }
}

}  // namespace rtt::nonseq
