// Obscurations in the tracer (#320, ADR 0050, point 3 and point 4, PR 4): a spider over a
// uniform disc, the sign of the rotation, a rotation by 90 deg, the decenter, the rim tolerance,
// the equality of a centred circle with the inner radius and the levels of the virtual rays.

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <optional>
#include <random>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/geom/conic.hpp"
#include "rtt/model/surface.hpp"
#include "rtt/trace/apply_event.hpp"
#include "rtt/trace/ray_batch.hpp"

using rtt::compile::CompiledSurface;
using rtt::math::Vec3;
using rtt::model::CircularAperture;
using rtt::model::CircularObscuration;
using rtt::model::EllipticalObscuration;
using rtt::model::Obscuration;
using rtt::model::Param;
using rtt::model::RectangularObscuration;
using rtt::trace::ObscurationTest;
using rtt::trace::RayState;
using rtt::trace::RayStatus;

namespace {

constexpr double kPi = std::numbers::pi;
constexpr double kTol = rtt::trace::kApertureTolerance;

/// A plane surface with an optional aperture and the given obscurations.
CompiledSurface surface(std::optional<rtt::model::Aperture> aperture,
                        std::vector<Obscuration> obscurations) {
  CompiledSurface s;
  s.shape = rtt::geom::Plane<double>{};
  s.aperture = std::move(aperture);
  s.obscurations = std::move(obscurations);
  return s;
}

/// Status after one sequential step of a ray along +z through the local point (x, y), starting
/// in front of the plane (z = -1, #151).
RayStatus status_at(const CompiledSurface& s,
                    double x,
                    double y,
                    ObscurationTest test = ObscurationTest::All) {
  RayState in;
  in.pos = Vec3(x, y, -1.0);
  in.obscuration_test = test;
  const RayState out =
      rtt::trace::sequential_step(in, s, 3, rtt::model::EventKind::Transmit, 1.0, 1.0);
  return out.status;
}

/// A uniform double in [0, 1) from 53 random bits (ADR 0050, point 4).
double unit(std::mt19937_64& g) {
  return static_cast<double>(g() >> 11U) * 0x1.0p-53;
}

/// Fraction of N points, uniform in the disc of radius `r_disc` (r = R sqrt(u), hergeleitet in
/// #320), that pass the surface.
double transmitted(const CompiledSurface& s, double r_disc, std::size_t n) {
  std::mt19937_64 g(20261119U);
  std::size_t pass = 0;
  for (std::size_t i = 0; i < n; ++i) {
    const double r = r_disc * std::sqrt(unit(g));
    const double phi = 2.0 * kPi * unit(g);
    if (status_at(s, r * std::cos(phi), r * std::sin(phi)) == RayStatus::Alive) ++pass;
  }
  return static_cast<double>(pass) / static_cast<double>(n);
}

}  // namespace

TEST_CASE("spider: a centred bar over a uniform disc blocks the strip area (ADR 0050)",
          "[obscuration]") {
  // Rectangle h_x = R, h_y = h over the disc of radius R blocks the strip |y| < h with
  // A = 2 (h sqrt(R^2 - h^2) + R^2 arcsin(h / R)) (hergeleitet in #320: the integral of
  // 2 sqrt(R^2 - y^2) over [-h, h]). Transmitted p = 1 - A / (pi R^2); binomial, 5 sigma.
  constexpr double kR = 12.0;
  constexpr double kH = 1.5;
  constexpr std::size_t kN = 1000000;
  const double area = 2.0 * (kH * std::sqrt(kR * kR - kH * kH) + kR * kR * std::asin(kH / kR));
  const double p = 1.0 - area / (kPi * kR * kR);
  const double bound = 5.0 * std::sqrt(p * (1.0 - p) / static_cast<double>(kN));
  CHECK(p > 0.83);  // hand value: A / (pi R^2) = 0.1587...
  CHECK(p < 0.85);
  const CompiledSurface spider =
      surface(CircularAperture{Param(kR), Param(0.0)},
              {RectangularObscuration{Param(kR), Param(kH), Param(0.0), Param(0.0), 0.0}});
  const double f = transmitted(spider, kR, kN);
  INFO("fraction " << f << ", expected " << p);
  CHECK(std::abs(f - p) <= bound);
  // Content guard: without the bar everything passes; a vertical bar (90 deg) blocks the same.
  CHECK(transmitted(surface(CircularAperture{Param(kR), Param(0.0)}, {}), kR, 10000) == 1.0);
  const CompiledSurface vertical = surface(
      std::nullopt, {RectangularObscuration{Param(kR), Param(kH), Param(0.0), Param(0.0), 90.0}});
  CHECK(std::abs(transmitted(vertical, kR, kN) - p) <= bound);
}

TEST_CASE("obscuration rotation: positive angles turn +x towards +y (ADR 0050)", "[obscuration]") {
  // Rectangle h_x = 2, h_y = 0.5 at theta = 30 deg. 1.5 (cos 30, sin 30) has u = 1.5, v = 0:
  // blocked. 1.5 (cos 30, -sin 30) has u = 1.5 cos 60 = 0.75, v = -1.5 sin 60 = -1.30 with
  // |v| > 0.5: passes. With the opposite sign it would be the other way round.
  const CompiledSurface s = surface(
      std::nullopt, {RectangularObscuration{Param(2.0), Param(0.5), Param(0.0), Param(0.0), 30.0}});
  const double c = std::cos(kPi / 6.0);
  const double sn = std::sin(kPi / 6.0);
  CHECK(status_at(s, 1.5 * c, 1.5 * sn) == RayStatus::Vignetted);
  CHECK(status_at(s, 1.5 * c, -1.5 * sn) == RayStatus::Alive);
  CHECK(std::abs(1.5 * std::sin(kPi / 3.0)) > 0.5);  // the hand value |v| = 1.30
  // The same for the ellipse with semi-axes 2 and 0.5.
  const CompiledSurface e = surface(
      std::nullopt, {EllipticalObscuration{Param(2.0), Param(0.5), Param(0.0), Param(0.0), 30.0}});
  CHECK(status_at(e, 1.5 * c, 1.5 * sn) == RayStatus::Vignetted);
  CHECK(status_at(e, 1.5 * c, -1.5 * sn) == RayStatus::Alive);
}

TEST_CASE("a rotation by 90 deg equals the swapped half widths (ADR 0050)", "[obscuration]") {
  // 201 x 201 grid over [-3, 3]^2, without a band of 1e-6 mm around the edges (cos 90 deg is
  // 6e-17, not 0, which moves the edges by about 1e-16 mm).
  for (const bool ellipse : {false, true}) {
    INFO((ellipse ? "ellipse" : "rectangle"));
    const auto make = [ellipse](double a, double b, double theta) -> Obscuration {
      if (ellipse) return EllipticalObscuration{Param(a), Param(b), Param(0.0), Param(0.0), theta};
      return RectangularObscuration{Param(a), Param(b), Param(0.0), Param(0.0), theta};
    };
    const CompiledSurface turned = surface(std::nullopt, {make(2.0, 0.7, 90.0)});
    const CompiledSurface swapped = surface(std::nullopt, {make(0.7, 2.0, 0.0)});
    std::size_t blocked = 0;
    std::size_t compared = 0;
    for (int i = 0; i <= 200; ++i) {
      for (int j = 0; j <= 200; ++j) {
        const double x = -3.0 + 0.03 * i;
        const double y = -3.0 + 0.03 * j;
        const double near_edge =
            ellipse ? std::abs(std::hypot(x / 0.7, y / 2.0) - 1.0) * 0.7
                    : std::min(std::abs(std::abs(x) - 0.7), std::abs(std::abs(y) - 2.0));
        if (near_edge < 1e-6) continue;
        const RayStatus a = status_at(turned, x, y);
        CHECK(a == status_at(swapped, x, y));
        if (a == RayStatus::Vignetted) ++blocked;
        ++compared;
      }
    }
    CHECK(blocked > 1000);  // content guard: the shapes block a part of the grid
    CHECK(blocked + 1000 < compared);
  }
}

TEST_CASE("a decentred circle blocks around its centre only (ADR 0050)", "[obscuration]") {
  // r = 1 < d = 3 at (d, 0): blocks (d, 0), passes (-d, 0) and the axis.
  const CompiledSurface s = surface(CircularAperture{Param(10.0), Param(0.0)},
                                    {CircularObscuration{Param(1.0), Param(3.0), Param(0.0)}});
  CHECK(status_at(s, 3.0, 0.0) == RayStatus::Vignetted);
  CHECK(status_at(s, -3.0, 0.0) == RayStatus::Alive);
  CHECK(status_at(s, 0.0, 0.0) == RayStatus::Alive);
  // Along y the same.
  const CompiledSurface t =
      surface(std::nullopt, {CircularObscuration{Param(1.0), Param(0.0), Param(-3.0)}});
  CHECK(status_at(t, 0.0, -3.0) == RayStatus::Vignetted);
  CHECK(status_at(t, 0.0, 3.0) == RayStatus::Alive);
  // A decentred circle is not a centred one: IgnoreCentredCircles still tests it.
  CHECK(status_at(s, 3.0, 0.0, ObscurationTest::IgnoreCentredCircles) == RayStatus::Vignetted);
  CHECK(status_at(s, 3.0, 0.0, ObscurationTest::IgnoreAll) == RayStatus::Alive);
}

TEST_CASE("obscuration edges belong to the aperture (ADR 0050)", "[obscuration]") {
  // Blocked strictly inside the shape shrunk by tol; a dimension - tol <= 0 blocks nothing.
  const CompiledSurface c =
      surface(std::nullopt, {CircularObscuration{Param(2.0), Param(0.0), Param(0.0)}});
  CHECK(status_at(c, 2.0, 0.0) == RayStatus::Alive);
  CHECK(status_at(c, 2.0 - 0.5 * kTol, 0.0) == RayStatus::Alive);
  CHECK(status_at(c, 2.0 - 2.0 * kTol, 0.0) == RayStatus::Vignetted);
  const CompiledSurface r = surface(
      std::nullopt, {RectangularObscuration{Param(2.0), Param(1.0), Param(0.0), Param(0.0), 0.0}});
  CHECK(status_at(r, 0.0, 1.0) == RayStatus::Alive);
  CHECK(status_at(r, 0.0, 1.0 - 2.0 * kTol) == RayStatus::Vignetted);
  CHECK(status_at(r, 2.0 - 0.5 * kTol, 0.0) == RayStatus::Alive);
  const CompiledSurface e = surface(
      std::nullopt, {EllipticalObscuration{Param(2.0), Param(1.0), Param(0.0), Param(0.0), 0.0}});
  CHECK(status_at(e, 0.0, 1.0) == RayStatus::Alive);
  CHECK(status_at(e, 0.0, 1.0 - 2.0 * kTol) == RayStatus::Vignetted);
  CHECK(status_at(e, 1.9, 0.5) == RayStatus::Alive);  // (1.9/2)^2 + (0.5/1)^2 > 1
  // Shapes thinner than the tolerance block nothing.
  const CompiledSurface thin =
      surface(std::nullopt,
              {CircularObscuration{Param(0.5 * kTol), Param(0.0), Param(0.0)},
               RectangularObscuration{Param(2.0), Param(0.5 * kTol), Param(0.0), Param(0.0), 0.0},
               EllipticalObscuration{Param(0.5 * kTol), Param(2.0), Param(0.0), Param(0.0), 0.0}});
  CHECK(status_at(thin, 0.0, 0.0) == RayStatus::Alive);
}

TEST_CASE("a centred circular obscuration decides exactly as the inner radius (ADR 0050)",
          "[obscuration]") {
  // Same arithmetic: x - 0 and y - 0 are exact, so every point gives the same status, also
  // within a few ulp of the inner rim (ADR 0050, point 3).
  constexpr double kRi = 3.0;
  const CompiledSurface ring = surface(CircularAperture{Param(10.0), Param(kRi)}, {});
  const CompiledSurface disc = surface(CircularAperture{Param(10.0), Param(0.0)},
                                       {CircularObscuration{Param(kRi), Param(0.0), Param(0.0)}});
  std::size_t blocked = 0;
  for (int k = 0; k < 720; ++k) {
    const double phi = 2.0 * kPi * k / 720.0;
    for (const double r : {0.0, 1.0, kRi - kTol, kRi - 1.5 * kTol, kRi - 0.5 * kTol, kRi,
                           std::nextafter(kRi - kTol, 0.0), std::nextafter(kRi - kTol, 10.0),
                           kRi + kTol, 5.0, 10.0, 10.0 + 2.0 * kTol}) {
      const double x = r * std::cos(phi);
      const double y = r * std::sin(phi);
      for (const ObscurationTest test :
           {ObscurationTest::All, ObscurationTest::IgnoreCentredCircles,
            ObscurationTest::IgnoreAll}) {
        const RayStatus a = status_at(ring, x, y, test);
        CHECK(a == status_at(disc, x, y, test));
        if (test == ObscurationTest::All && a == RayStatus::Vignetted) ++blocked;
      }
    }
  }
  CHECK(blocked > 2000);  // content guard: the obscuration and the outer rim block rays
}

TEST_CASE("levels of the obscuration test (ADR 0050, point 3)", "[obscuration]") {
  // A centred circle and a spider bar: All tests both, IgnoreCentredCircles only the bar,
  // IgnoreAll none. The outer rim counts at every level.
  const CompiledSurface s =
      surface(CircularAperture{Param(10.0), Param(1.0)},
              {CircularObscuration{Param(2.0), Param(0.0), Param(0.0)},
               RectangularObscuration{Param(10.0), Param(0.25), Param(0.0), Param(0.0), 0.0}});
  const auto at = [&](double x, double y, ObscurationTest t) { return status_at(s, x, y, t); };
  // (0, 1.5): in the centred circle, off the bar.
  CHECK(at(0.0, 1.5, ObscurationTest::All) == RayStatus::Vignetted);
  CHECK(at(0.0, 1.5, ObscurationTest::IgnoreCentredCircles) == RayStatus::Alive);
  CHECK(at(0.0, 1.5, ObscurationTest::IgnoreAll) == RayStatus::Alive);
  // (5, 0): on the bar.
  CHECK(at(5.0, 0.0, ObscurationTest::All) == RayStatus::Vignetted);
  CHECK(at(5.0, 0.0, ObscurationTest::IgnoreCentredCircles) == RayStatus::Vignetted);
  CHECK(at(5.0, 0.0, ObscurationTest::IgnoreAll) == RayStatus::Alive);
  // (0, 0.5): inside the inner radius, a centred circle as well.
  CHECK(at(0.0, 0.5, ObscurationTest::IgnoreCentredCircles) == RayStatus::Alive);
  // (0, 5): free; (0, 11): outside the outer rim at every level.
  for (const ObscurationTest t :
       {ObscurationTest::All, ObscurationTest::IgnoreCentredCircles, ObscurationTest::IgnoreAll}) {
    CHECK(at(0.0, 5.0, t) == RayStatus::Alive);
    CHECK(at(0.0, 11.0, t) == RayStatus::Vignetted);
  }
  // The compatible flag ignore_obscurations = true is IgnoreAll whatever the level says.
  RayState in;
  in.pos = Vec3(5.0, 0.0, -1.0);
  in.ignore_obscurations = true;
  in.obscuration_test = ObscurationTest::All;
  CHECK(rtt::trace::effective_obscuration_test(in) == ObscurationTest::IgnoreAll);
  CHECK(rtt::trace::sequential_step(in, s, 3, rtt::model::EventKind::Transmit, 1.0, 1.0).status ==
        RayStatus::Alive);
  // in_obscuration() directly, with its default level All.
  CHECK(rtt::trace::in_obscuration(s, 5.0, 0.0));
  CHECK_FALSE(rtt::trace::in_obscuration(s, 0.0, 5.0));
  CHECK_FALSE(rtt::trace::in_obscuration(s, 5.0, 0.0, ObscurationTest::IgnoreAll));
}

TEST_CASE("RayBatch: obscuration level and its compatible bool form (ADR 0050)", "[obscuration]") {
  rtt::trace::RayBatch rays(2);
  CHECK(rays.obscuration_test() == ObscurationTest::All);
  rays.set_obscuration_test(ObscurationTest::IgnoreCentredCircles);
  CHECK(rays.obscuration_test() == ObscurationTest::IgnoreCentredCircles);
  CHECK(rays.ignores_obscurations());  // != All
  rays.set_ignore_obscurations(true);
  CHECK(rays.obscuration_test() == ObscurationTest::IgnoreAll);
  rays.resize(5);
  CHECK(rays.obscuration_test() == ObscurationTest::IgnoreAll);  // kept by resize()
  rays.set_ignore_obscurations(false);
  CHECK(rays.obscuration_test() == ObscurationTest::All);
  CHECK_FALSE(rays.ignores_obscurations());
}

TEST_CASE("a surface without aperture may carry obscurations (ADR 0050)", "[obscuration]") {
  // A window with a spider: unlimited except for the bar.
  const CompiledSurface s =
      surface(std::nullopt,
              {RectangularObscuration{Param(50.0), Param(0.5), Param(0.0), Param(0.0), 45.0}});
  CHECK(status_at(s, 100.0, 0.0) == RayStatus::Alive);
  CHECK(status_at(s, 10.0, 10.0) == RayStatus::Vignetted);  // on the bar at 45 deg
  CHECK(status_at(s, 10.0, -10.0) == RayStatus::Alive);
  // The vignetted ray stands at its hit point on the surface.
  RayState in;
  in.pos = Vec3(10.0, 10.0, -1.0);
  const RayState out =
      rtt::trace::sequential_step(in, s, 7, rtt::model::EventKind::Transmit, 1.0, 1.0);
  CHECK(out.last_surface == 7);
  CHECK(out.pos.z() == 0.0);
}
