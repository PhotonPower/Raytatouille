#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <numbers>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/geom/asphere.hpp"
#include "rtt/geom/conic.hpp"
#include "rtt/math/isometry.hpp"
#include "rtt/trace/apply_event.hpp"

using rtt::compile::CompiledSurface;
using rtt::math::Isometry3;
using rtt::math::Vec3;
using rtt::model::EventKind;
using rtt::trace::RayState;
using rtt::trace::RayStatus;

// Laws of refraction and reflection: B. de Greve, Reflections and Refractions in Ray Tracing
// (2006), Eqs. (13), (22)-(25), Sec. 6; see also M. Born, E. Wolf, Principles of Optics, 7th ed.,
// Sec. 3.2.2. Refraction: n2 sin(theta2) = n1 sin(theta1), the refracted direction lies in the
// plane of incidence; critical angle arcsin(n2 / n1) (Eq. (25)). Reflection: r = d - 2 (d . n) n.

namespace {
constexpr double kAngleTol = 1e-14;  // rad, issue #6

bool near(const Vec3& a, const Vec3& b, double tol) {
  return (a - b).cwiseAbs().maxCoeff() <= tol;
}

/// Plane surface at the global origin (local = global).
CompiledSurface plane_surface() {
  CompiledSurface s;
  s.shape = rtt::geom::Plane<double>{};
  return s;
}

/// Ray in the y-z plane with angle `theta` to +z, aimed at the global origin from z = -10 mm.
RayState ray_at_angle(double theta) {
  RayState r;
  r.dir = Vec3(0.0, std::sin(theta), std::cos(theta));
  r.pos = -10.0 / r.dir.z() * r.dir;
  return r;
}

/// Angle between a unit direction and the +z axis (normal of the plane), rad.
double angle_to_axis(const Vec3& d) {
  return std::atan2(std::hypot(d.x(), d.y()), std::abs(d.z()));
}
}  // namespace

TEST_CASE("refraction at a plane follows Snell's law to 1e-14 rad", "[apply_event]") {
  const CompiledSurface plane = plane_surface();
  for (const auto& [n1, n2] : {std::pair{1.0, 1.5}, std::pair{1.5, 1.0}, std::pair{1.333, 1.7}}) {
    for (const double theta_deg : {0.0, 5.0, 20.0, 35.0, 40.0}) {
      const double theta = theta_deg * std::numbers::pi / 180.0;
      if (n1 * std::sin(theta) >= n2) {
        continue;
      }
      const RayState out =
          rtt::trace::sequential_step(ray_at_angle(theta), plane, 3, EventKind::Refract, n1, n2);
      INFO("n1 = " << n1 << ", n2 = " << n2 << ", theta = " << theta_deg << " deg");
      REQUIRE(out.status == RayStatus::Alive);
      REQUIRE(out.last_surface == 3);
      const double expected = std::asin(n1 * std::sin(theta) / n2);
      REQUIRE(std::abs(angle_to_axis(out.dir) - expected) <= kAngleTol);
      REQUIRE(std::abs(out.dir.norm() - 1.0) <= 1e-15);
      REQUIRE(out.dir.x() == 0.0);  // stays in the plane of incidence
      REQUIRE(out.dir.y() >= 0.0);  // on the other side of the normal
      REQUIRE(out.dir.z() > 0.0);   // transmitted, not reflected
      REQUIRE(near(out.pos, Vec3::Zero(), 1e-12));
    }
  }
}

TEST_CASE("refraction is independent of the orientation of the surface normal", "[apply_event]") {
  // Ray travelling -z hits the plane from the +z side; the surface normal is +z.
  const CompiledSurface plane = plane_surface();
  RayState in;
  const double theta = 0.3;
  in.dir = Vec3(std::sin(theta), 0.0, -std::cos(theta));
  in.pos = -5.0 / std::cos(theta) * in.dir;  // starts at z = +5 mm
  const RayState out = rtt::trace::sequential_step(in, plane, 0, EventKind::Refract, 1.0, 1.5);
  REQUIRE(out.status == RayStatus::Alive);
  REQUIRE(out.dir.z() < 0.0);
  REQUIRE(std::abs(angle_to_axis(out.dir) - std::asin(std::sin(theta) / 1.5)) <= kAngleTol);
}

TEST_CASE("rays travelling -z: refraction and TIR from Snell's law by hand", "[apply_event]") {
  // The surface normal from rtt-geom is +z, so for d_z < 0 it points into the incident medium
  // already; for d_z > 0 it has to be flipped (de Greve, Sec. 6). Both cases must agree.
  const CompiledSurface plane = plane_surface();
  // d = (0, sin a, -cos a) from glass n1 = 1.5 into n2 = 1.2, a = 0.5 rad. By hand:
  // sin a' = 1.5 sin(0.5) / 1.2, refracted d' = (0, sin a', -cos a').
  const double a = 0.5;
  RayState in;
  in.dir = Vec3(0.0, std::sin(a), -std::cos(a));
  in.pos = -4.0 / std::cos(a) * in.dir;  // starts at z = +4 mm
  const RayState out = rtt::trace::sequential_step(in, plane, 0, EventKind::Refract, 1.5, 1.2);
  REQUIRE(out.status == RayStatus::Alive);
  const double sin_t = 1.5 * std::sin(a) / 1.2;
  REQUIRE(near(out.dir, Vec3(0.0, sin_t, -std::sqrt(1.0 - sin_t * sin_t)), 1e-15));
  REQUIRE(std::abs(angle_to_axis(out.dir) - std::asin(sin_t)) <= kAngleTol);
  // Total internal reflection in -z direction: critical angle arcsin(1.2 / 1.5).
  const double critical = std::asin(1.2 / 1.5);
  for (const double delta : {-1e-9, 1e-9}) {
    RayState r;
    r.dir = Vec3(0.0, std::sin(critical + delta), -std::cos(critical + delta));
    r.pos = -4.0 / std::cos(critical + delta) * r.dir;  // starts at z = +4 mm
    const RayState t = rtt::trace::sequential_step(r, plane, 1, EventKind::Refract, 1.5, 1.2);
    REQUIRE(t.status == (delta < 0.0 ? RayStatus::Alive : RayStatus::Tir));
  }
}

TEST_CASE("total internal reflection just above the critical angle", "[apply_event]") {
  // Critical angle arcsin(n2 / n1) from glass n1 = 1.5 into air n2 = 1.
  const CompiledSurface plane = plane_surface();
  const double critical = std::asin(1.0 / 1.5);
  const RayState below = rtt::trace::sequential_step(ray_at_angle(critical - 1e-9), plane, 1,
                                                     EventKind::Refract, 1.5, 1.0);
  REQUIRE(below.status == RayStatus::Alive);
  const RayState above = rtt::trace::sequential_step(ray_at_angle(critical + 1e-9), plane, 1,
                                                     EventKind::Refract, 1.5, 1.0);
  REQUIRE(above.status == RayStatus::Tir);
  // The ray stops at the hit point and remembers the surface.
  REQUIRE(near(above.pos, Vec3::Zero(), 1e-12));
  REQUIRE(above.last_surface == 1);
}

TEST_CASE("reflection mirrors the direction about the normal", "[apply_event]") {
  const CompiledSurface plane = plane_surface();
  const double theta = 0.4;
  const RayState out =
      rtt::trace::sequential_step(ray_at_angle(theta), plane, 0, EventKind::Reflect, 1.0, 1.0);
  REQUIRE(out.status == RayStatus::Alive);
  REQUIRE(near(out.dir, Vec3(0.0, std::sin(theta), -std::cos(theta)), 1e-15));
}

TEST_CASE("transmit keeps the direction", "[apply_event]") {
  const CompiledSurface plane = plane_surface();
  const RayState in = ray_at_angle(0.2);
  const RayState out = rtt::trace::sequential_step(in, plane, 0, EventKind::Transmit, 1.0, 1.0);
  REQUIRE(out.status == RayStatus::Alive);
  REQUIRE(out.dir == in.dir);
}

TEST_CASE("optical path length accumulates n_before times the geometric path", "[apply_event]") {
  const CompiledSurface plane = plane_surface();
  RayState in = ray_at_angle(0.25);
  in.opl = 2.0;
  const double distance = (Vec3::Zero() - in.pos).norm();
  const RayState out = rtt::trace::sequential_step(in, plane, 0, EventKind::Refract, 1.333, 1.5);
  REQUIRE(std::abs(out.opl - (2.0 + 1.333 * distance)) <= 1e-12);
}

TEST_CASE("surface pose: intersection in local coordinates, result in global", "[apply_event]") {
  // Plane tilted by 45 deg about x at z = 20 mm; axial ray reflects into +y.
  CompiledSurface tilted = plane_surface();
  tilted.to_global = Isometry3::from_pose(Vec3(0.0, 0.0, 20.0), Vec3(45.0, 0.0, 0.0), Vec3::Zero());
  tilted.to_local = tilted.to_global.inverse();
  RayState in;
  in.pos = Vec3(0.0, 0.0, -5.0);
  const RayState out = rtt::trace::sequential_step(in, tilted, 7, EventKind::Reflect, 1.0, 1.0);
  REQUIRE(out.status == RayStatus::Alive);
  REQUIRE(near(out.pos, Vec3(0.0, 0.0, 20.0), 1e-12));
  // Rx(45 deg) turns the local normal +z into n = (0, -sin 45, cos 45). For d = (0, 0, 1):
  // d . n = cos 45 and r = d - 2 (d . n) n = (0, 2 sin 45 cos 45, 1 - 2 cos^2 45) = (0, 1, 0).
  REQUIRE(near(out.dir, Vec3(0.0, 1.0, 0.0), 1e-15));
  REQUIRE(std::abs(out.opl - 25.0) <= 1e-12);
}

TEST_CASE("sphere and even asphere surfaces are intersected", "[apply_event]") {
  // Axial ray onto the vertex of a sphere and of an even asphere: normal +z, no deviation.
  for (const rtt::compile::CompiledShape& shape :
       {rtt::compile::CompiledShape{rtt::geom::Conic<double>(0.02, 0.0)},
        rtt::compile::CompiledShape{rtt::geom::EvenAsphere<double>(0.02, -1.0, {1e-5})}}) {
    CompiledSurface s;
    s.shape = shape;
    RayState in;
    in.pos = Vec3(0.0, 0.0, -3.0);
    const RayState out = rtt::trace::sequential_step(in, s, 0, EventKind::Refract, 1.0, 1.5);
    REQUIRE(out.status == RayStatus::Alive);
    REQUIRE(near(out.pos, Vec3::Zero(), 1e-12));
    REQUIRE(near(out.dir, Vec3::UnitZ(), 1e-15));
  }
  // Off-axis on the asphere: the hit lies on the surface.
  CompiledSurface s;
  const rtt::geom::EvenAsphere<double> asphere(0.02, -1.0, {1e-5});
  s.shape = asphere;
  RayState in;
  in.pos = Vec3(3.0, -4.0, -3.0);
  const RayState out = rtt::trace::sequential_step(in, s, 0, EventKind::Refract, 1.0, 1.5);
  REQUIRE(out.status == RayStatus::Alive);
  REQUIRE(std::abs(out.pos.z() - asphere.sag(3.0, -4.0)) <= 1e-12);
}

TEST_CASE("intersection failures become Missed and keep the ray unchanged", "[apply_event]") {
  CompiledSurface sphere;
  sphere.shape = rtt::geom::Conic<double>(0.02, 0.0);  // R = 50
  RayState in;
  in.pos = Vec3(0.0, 60.0, -5.0);
  in.opl = 1.0;
  in.last_surface = 2;
  const RayState out = rtt::trace::sequential_step(in, sphere, 4, EventKind::Refract, 1.0, 1.5);
  REQUIRE(out.status == RayStatus::Missed);
  REQUIRE(out.pos == in.pos);
  REQUIRE(out.dir == in.dir);
  REQUIRE(out.opl == in.opl);
  REQUIRE(out.last_surface == 2);
}

TEST_CASE("apertures vignette in local coordinates", "[apply_event]") {
  // Axial rays at (x, y) onto a plane at the origin.
  const auto status_at = [](const rtt::model::Aperture& aperture, double x, double y) {
    CompiledSurface s = plane_surface();
    s.aperture = aperture;
    RayState in;
    in.pos = Vec3(x, y, -1.0);
    const RayState out = rtt::trace::sequential_step(in, s, 5, EventKind::Transmit, 1.0, 1.0);
    if (out.status == RayStatus::Vignetted) {
      // Stops at the hit point on the vignetting surface.
      REQUIRE(out.pos == Vec3(x, y, 0.0));
      REQUIRE(out.last_surface == 5);
    }
    return out.status;
  };
  const rtt::model::CircularAperture ring{10.0, 2.0};
  REQUIRE(status_at(ring, 0.0, 5.0) == RayStatus::Alive);
  REQUIRE(status_at(ring, 6.0, 8.0) == RayStatus::Alive);  // r = 10, on the edge
  REQUIRE(status_at(ring, 0.0, 10.5) == RayStatus::Vignetted);
  REQUIRE(status_at(ring, 1.0, 1.0) == RayStatus::Vignetted);  // inside the obscuration
  const rtt::model::RectangularAperture rect{4.0, 2.0};
  REQUIRE(status_at(rect, 3.9, -1.9) == RayStatus::Alive);
  REQUIRE(status_at(rect, 0.0, 2.5) == RayStatus::Vignetted);
  REQUIRE(status_at(rect, 4.5, 0.0) == RayStatus::Vignetted);
  const rtt::model::EllipticalAperture ellipse{4.0, 2.0};
  REQUIRE(status_at(ellipse, 3.9, 0.0) == RayStatus::Alive);
  REQUIRE(status_at(ellipse, 3.0, 1.5) == RayStatus::Vignetted);  // (3/4)^2 + (1.5/2)^2 > 1
}

TEST_CASE("events of later milestones are EventImpossible", "[apply_event]") {
  const CompiledSurface plane = plane_surface();
  for (const EventKind kind :
       {EventKind::Diffract, EventKind::Ordinary, EventKind::Extraordinary}) {
    const RayState out = rtt::trace::sequential_step(ray_at_angle(0.1), plane, 2, kind, 1.0, 1.5);
    REQUIRE(out.status == RayStatus::EventImpossible);
    REQUIRE(near(out.pos, Vec3::Zero(), 1e-12));
    REQUIRE(out.last_surface == 2);
  }
}

TEST_CASE("absorber stops the ray with status Absorbed", "[apply_event]") {
  CompiledSurface s = plane_surface();
  s.interaction = rtt::model::Absorber{};
  const RayState out =
      rtt::trace::sequential_step(ray_at_angle(0.1), s, 6, EventKind::Transmit, 1.0, 1.0);
  REQUIRE(out.status == RayStatus::Absorbed);
  REQUIRE(near(out.pos, Vec3::Zero(), 1e-12));
  REQUIRE(out.last_surface == 6);
}

TEST_CASE("rays that are not alive are left untouched", "[apply_event]") {
  const CompiledSurface plane = plane_surface();
  RayState in = ray_at_angle(0.1);
  in.status = RayStatus::Vignetted;
  const RayState out = rtt::trace::sequential_step(in, plane, 2, EventKind::Refract, 1.0, 1.5);
  REQUIRE(out.status == RayStatus::Vignetted);
  REQUIRE(out.pos == in.pos);
  REQUIRE(out.dir == in.dir);
}

TEST_CASE("Newton failure on an even asphere becomes NoConvergence", "[apply_event]") {
  // Spherical base R = 50 bounds the asphere at r = 50 mm; the strong negative A4 term moves the
  // surface far below the base sphere, so the Newton step of an oblique ray leaves the domain.
  CompiledSurface s;
  s.shape = rtt::geom::EvenAsphere<double>(0.02, 0.0, {-1e-4});
  RayState in;
  in.pos = Vec3(40.0, 0.0, -5.0);
  in.dir = Vec3(0.3, 0.0, 1.0).normalized();
  in.opl = 1.5;
  const RayState out = rtt::trace::sequential_step(in, s, 3, EventKind::Refract, 1.0, 1.5);
  REQUIRE(out.status == RayStatus::NoConvergence);
  REQUIRE(out.pos == in.pos);
  REQUIRE(out.dir == in.dir);
  REQUIRE(out.opl == in.opl);
  REQUIRE(out.last_surface == rtt::trace::kNoSurface);
}

TEST_CASE("apply_event does not check the aperture (reuse in non-sequential tracing)",
          "[apply_event]") {
  // The aperture decision belongs to the caller: sequential tracing vignettes, non-sequential
  // tracing (M9) will treat the hit as "surface not there".
  CompiledSurface s = plane_surface();
  s.aperture = rtt::model::CircularAperture{1.0, 0.0};
  RayState in;
  in.pos = Vec3(0.0, 5.0, -1.0);
  const auto hit = rtt::trace::intersect_surface(in, s);
  REQUIRE(hit.status == rtt::geom::HitStatus::Hit);
  REQUIRE_FALSE(rtt::trace::inside_aperture(s, hit));
  const RayState out = rtt::trace::apply_event(in, s, hit, 0, EventKind::Transmit, 1.0, 1.0);
  REQUIRE(out.status == RayStatus::Alive);
  REQUIRE(near(out.pos, Vec3(0.0, 5.0, 0.0), 1e-15));
  REQUIRE(rtt::trace::sequential_step(in, s, 0, EventKind::Transmit, 1.0, 1.0).status ==
          RayStatus::Vignetted);
}
