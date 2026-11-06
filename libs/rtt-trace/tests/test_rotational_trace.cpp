// Real rays through the odd asphere and the axicon (#314, ADR 0043, acceptance 6 and 8): an
// axicon refracts an axis-parallel ray by the law of refraction at a cone (de Greve, Reflections
// and Refractions in Ray Tracing, docs/quellen.md), and a ray through an odd asphere with A3
// agrees with an independent computation (bisection on the written-out sag, normal from its
// derivative by hand, vector form of Snell's law).

#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <string>
#include <variant>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/ray_paths.hpp"
#include "rtt/trace/sequential.hpp"

using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::math::Vec3;
using rtt::model::Element;
using rtt::model::Param;
using rtt::model::System;
using rtt::trace::RayBatch;
using rtt::trace::RayPaths;
using rtt::trace::RayStatus;
using rtt::trace::SequentialTracer;

namespace {

const double kPi = std::acos(-1.0);

/// m1/singlet_const in vacuum (stop at z = 0, L1.S1 at z = 5, L1.S2 at z = 9, image at
/// z = 106.363) with the bases `front` and `back` and the index `material`.
System singlet(const rtt::model::BaseShape& front,
               const rtt::model::BaseShape& back,
               const std::string& material) {
  System s = rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m1/singlet_const.rtt.json");
  s.environment.medium = "VACUUM";
  Element& lens = std::get<Element>(s.root.children[1].value);
  lens.material = material;
  lens.surfaces[0].shape.base = front;
  lens.surfaces[1].shape.base = back;
  return s;
}

/// Traces rays from `starts` along `directions` and returns their paths; every ray arrives.
RayPaths traced(const System& s, const std::vector<Vec3>& starts, const std::vector<Vec3>& dirs) {
  const rtt::material::MaterialLibrary lib;
  const CompiledSystem cs = rtt::compile::compile(s, lib);
  RayBatch rays(starts.size());
  for (std::size_t k = 0; k < rays.size(); ++k) {
    rays.pos_x()[k] = starts[k].x();
    rays.pos_y()[k] = starts[k].y();
    rays.pos_z()[k] = starts[k].z();
    rays.dir_x()[k] = dirs[k].x();
    rays.dir_y()[k] = dirs[k].y();
    rays.dir_z()[k] = dirs[k].z();
  }
  RayPaths paths;
  const auto stats = SequentialTracer().trace(cs, PathId{0}, rays, paths);
  REQUIRE(stats.count(RayStatus::Alive) == rays.size());
  REQUIRE(paths.slots == 5);  // start, stop, L1.S1, L1.S2, image
  return paths;
}

/// Vector form of Snell's law: d unit, n unit normal with n . d > 0, from index n1 into n2.
Vec3 refract(const Vec3& d, const Vec3& n, double n1, double n2) {
  const double mu = n1 / n2;
  const double ci = n.dot(d);
  const double ct = std::sqrt(1.0 - mu * mu * (1.0 - ci * ci));
  return mu * d + (ct - mu * ci) * n;
}

}  // namespace

TEST_CASE("an axicon refracts an axis-parallel ray by the cone angle", "[trace][rotational]") {
  // ADR 0043, acceptance 8: plano-convex, n = 1.5, the cone at the back with alpha = -5 deg (the
  // mantle falls towards -z: a convex cone). The ray in the glass parallel to the axis at
  // h = 2 mm meets the cone at z = -h tan 5 deg from the apex; the normal there is tilted by 5 deg,
  // so the ray leaves at beta = asin(n sin 5 deg) - 5 deg towards the axis and crosses it at
  // z = -h tan 5 deg + h / tan(beta) from the apex.
  const double n = 1.5;
  const double h = 2.0;
  const double alpha = 5.0 * kPi / 180.0;
  const System s = singlet(rtt::model::Plane{}, rtt::model::Axicon{Param(-5.0)}, "CONST:1.5");
  const RayPaths paths =
      traced(s, {Vec3(0.0, h, -1.0), Vec3(0.0, 0.0, -1.0)}, {Vec3::UnitZ(), Vec3::UnitZ()});
  const double apex = 9.0;  // global z of L1.S2
  const double z_hit = -h * std::tan(alpha);
  const double beta = std::asin(n * std::sin(alpha)) - alpha;
  const double z_axis = z_hit + h / std::tan(beta);
  // The closed values against the numbers of ADR 0043, acceptance 8 (content guard).
  REQUIRE(std::abs(z_hit + 0.174977327) <= 1e-9);
  REQUIRE(std::abs(beta * 180.0 / kPi - 2.511987114718) <= 1e-12);
  REQUIRE(std::abs(z_axis - 45.413683474) <= 1e-9);

  const Vec3 p = paths.position_at(0, 3);
  const Vec3 d = paths.direction_at(0, 3);
  CHECK(std::abs(p.y() - h) <= 1e-15);
  CHECK(std::abs(p.z() - apex - z_hit) <= 1e-13);
  CHECK(d.x() == 0.0);
  CHECK(d.y() < 0.0);  // towards the axis
  // Angle to the axis within 1e-13 rad.
  CHECK(std::abs(std::atan2(-d.y(), d.z()) - beta) <= 1e-13);
  // Axis crossing within 1e-9 mm: h / sin^2(beta) ~ 1.04e3 times the error of the direction.
  const double crossing = p.z() + (-p.y() / d.y()) * d.z() - apex;
  CHECK(std::abs(crossing - z_axis) <= 1e-9);

  // The axial ray passes the apex undeviated.
  CHECK(paths.position_at(1, 3) == Vec3(0.0, 0.0, apex));
  CHECK(paths.direction_at(1, 3) == Vec3::UnitZ());
}

TEST_CASE("a real ray through an odd asphere with A3 against an independent computation",
          "[trace][rotational]") {
  // ADR 0043, acceptance 6: the front L1.S1 at z = 5 is an odd asphere with A3 (order 2, no
  // Seidel cap); a ray from (0, 3, -1) along (0, 0.05, 1) in vacuum enters CONST:1.5168.
  const double c = 0.02;
  const double k = -0.5;
  const std::vector<double> a = {0.0, 0.0, 2e-4, 1e-6};
  const System s = singlet(
      rtt::model::OddAsphere{
          Param(c), Param(k), {Param(a[0]), Param(a[1]), Param(a[2]), Param(a[3])}},
      rtt::model::Plane{}, "CONST:1.5168");
  const Vec3 o(0.0, 3.0, -1.0);
  const Vec3 d = Vec3(0.0, 0.05, 1.0).normalized();
  const RayPaths paths = traced(s, {o}, {d});

  // Written-out sag and its radial derivative (OpticStudio "Odd Asphere"; Forbes 2011,
  // Eq. (2.1)): z = c r^2 / (1 + phi) + sum A_n r^n, dz/dr = c r / phi + sum n A_n r^(n - 1).
  const auto sag = [&](double r) {
    double z = c * r * r / (1.0 + std::sqrt(1.0 - (1.0 + k) * c * c * r * r));
    for (std::size_t i = 0; i < a.size(); ++i) z += a[i] * std::pow(r, static_cast<double>(i + 1));
    return z;
  };
  const auto slope = [&](double r) {
    double g = c * r / std::sqrt(1.0 - (1.0 + k) * c * c * r * r);
    for (std::size_t i = 0; i < a.size(); ++i) {
      g += static_cast<double>(i + 1) * a[i] * std::pow(r, static_cast<double>(i));
    }
    return g;
  };
  // F(t) = z_local(t) - sag(r(t)) in the local frame of L1.S1 (vertex at z = 5); bisection on
  // [0, 20], where F changes sign once (F(0) < 0 in front, F(20) > 0 behind the surface).
  const Vec3 vertex(0.0, 0.0, 5.0);
  const auto f = [&](double t) {
    const Vec3 q = o + t * d - vertex;
    return q.z() - sag(std::hypot(q.x(), q.y()));
  };
  double lo = 0.0;
  double hi = 20.0;
  REQUIRE(f(lo) < 0.0);
  REQUIRE(f(hi) > 0.0);
  for (int i = 0; i < 200; ++i) {
    const double mid = lo + (hi - lo) / 2.0;
    if (mid == lo || mid == hi) break;
    (f(mid) < 0.0 ? lo : hi) = mid;
  }
  const Vec3 hit = o + lo * d;
  const Vec3 q = hit - vertex;
  const double r = std::hypot(q.x(), q.y());
  const double g = slope(r);
  const Vec3 normal = Vec3(-g * q.x() / r, -g * q.y() / r, 1.0).normalized();
  const Vec3 out = refract(d, normal, 1.0, 1.5168);

  const Vec3 p = paths.position_at(0, 2);
  CHECK((p - hit).cwiseAbs().maxCoeff() <= 1e-12);
  CHECK((paths.direction_at(0, 2) - out).cwiseAbs().maxCoeff() <= 1e-12);
  // Content guard: A3 moves the hit by more than 1e-4 mm from the base conic's.
  const double conic_only = c * r * r / (1.0 + std::sqrt(1.0 - (1.0 + k) * c * c * r * r));
  CHECK(std::abs(q.z() - conic_only) > 1e-3);
}
